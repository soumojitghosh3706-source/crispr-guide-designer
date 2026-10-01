"""
g6pd_retry_flagged.py  --  Case study 5 (G6PD): independent re-verification of
"suspiciously Clean" guides.

WHY: in g6pd_blast_offtargets.py a guide with raw_hits > 1 but genomic_hits == 0
could mean (a) the 2nd hit really was G6PD itself / a transcript record (fine), or
(b) a fetch error during verification was swallowed and a real off-target was
dropped (false Clean). The progress CSV cannot tell (a) from (b). This script
re-BLASTs those guides and prints EVERY hit with the exact reason it was kept or
dropped. Nothing is swallowed: a failed step makes the guide UNVERIFIED, never Clean.

SELF-detection here is done OFFLINE: the aligned subject sequence of each hit is
compared with g6pd_gene.fasta (and its reverse complement). No network needed,
so no fetch error can hide a self-hit decision. Only the PAM check needs efetch.

USAGE (run from the SAME folder that holds the main script's CSVs):
    python g6pd_retry_flagged.py          # auto-selects flagged guides (raw_hits>1, genomic_hits==0)
    python g6pd_retry_flagged.py --all    # re-check all guides in the verdict CSV

It NEVER edits your existing CSVs. It writes:
    g6pd_retry_hits.csv      one row per BLAST hit, with the reason
    g6pd_retry_verdict.csv   old vs new verdict per guide
"""
import csv
import os
import re
import sys
import time

GENE_NAME = "G6PD"
GENE_FASTA = "g6pd_gene.fasta"
PROGRESS_CSV = "g6pd_blast_progress.csv"
VERDICT_CSV = "g6pd_offtarget_verdict.csv"
OUT_HITS = "g6pd_retry_hits.csv"
OUT_VERDICT = "g6pd_retry_verdict.csv"

EMAIL = "your_email@example.com"   # <-- use the SAME email as in g6pd_blast_offtargets.py

MAX_MISMATCH = 4
GUIDE_LEN = 20
GENOMIC_PREFIXES = ("NC_", "NG_", "NT_", "NW_", "AC_")   # same allow-list as main script
RETRIES = 3
PAUSE_BETWEEN_SEARCHES = 3
PAUSE_BETWEEN_FETCHES = 1


# ---------------------------------------------------------------- pure helpers
def revcomp(s):
    return s.upper().translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def read_fasta_seq(path):
    with open(path) as fh:
        return "".join(l.strip() for l in fh if not l.startswith(">")).upper()


def is_genomic_accession(acc):
    return acc.startswith(GENOMIC_PREFIXES)


def is_self_sequence(subject_seq, gene_seq):
    """True if the aligned subject DNA occurs verbatim in the gene (either strand)."""
    s = subject_seq.replace("-", "").upper()
    if len(s) < 12:
        return False
    return s in gene_seq or revcomp(s) in gene_seq


def pam_fetch_params(hit_start, hit_end):
    """Same convention as the main script (hit_end is the reference in BOTH orientations).
    plus strand  (start <= end): PAM = end+1..end+3, strand 1
    minus strand (start >  end): PAM = end-3..end-1, strand 2"""
    if hit_start <= hit_end:
        return hit_end + 1, hit_end + 3, 1
    return hit_end - 3, hit_end - 1, 2


def classify_hit(acc, sbjct_seq, align_len, identities, gaps, gene_seq):
    """Returns (category, mismatches, reason). Categories:
    SELF | EXCL_NONGENOMIC | PARTIAL | TOO_MANY_MM | CANDIDATE"""
    mm = align_len - identities if align_len else None
    if is_self_sequence(sbjct_seq, gene_seq):
        return "SELF", mm, f"aligned subject DNA occurs verbatim in {GENE_NAME} gene sequence"
    if not is_genomic_accession(acc):
        if align_len == GUIDE_LEN and gaps == 0 and mm <= MAX_MISMATCH:
            return ("EXCL_NONGENOMIC", mm,
                    "NOT self, <=4mm, but accession is not RefSeq genomic -> dropped by design; REVIEW MANUALLY")
        return "EXCL_NONGENOMIC", mm, "non-genomic accession (transcript/clone), weak or partial"
    if align_len != GUIDE_LEN or gaps != 0:
        return "PARTIAL", mm, "not a full-length ungapped 20-nt alignment (known limitation, not evaluated)"
    if mm > MAX_MISMATCH:
        return "TOO_MANY_MM", mm, f"{mm} mismatches > {MAX_MISMATCH}"
    return "CANDIDATE", mm, "non-self genomic hit within mismatch limit -> PAM check"


def collapse_loci(cands):
    """cands: list of dicts with acc, lo, hi, mm. Merge overlapping intervals per accession, keep lowest mm."""
    cands = sorted(cands, key=lambda c: (c["acc"], c["lo"]))
    out = []
    for c in cands:
        if out and out[-1]["acc"] == c["acc"] and c["lo"] <= out[-1]["hi"]:
            m = out[-1]
            m["hi"] = max(m["hi"], c["hi"])
            if c["mm"] < m["mm"]:
                m["mm"], m["start"], m["end"] = c["mm"], c["start"], c["end"]
        else:
            out.append(dict(c))
    return out


# --------------------------------------------------------------------- network
def blast_guide(guide):
    from Bio.Blast import NCBIWWW, NCBIXML
    last = None
    for attempt in range(1, RETRIES + 1):
        try:
            h = NCBIWWW.qblast("blastn", "nt", guide, entrez_query="Homo sapiens[Organism]",
                               word_size=7, expect=1000, hitlist_size=50)
            return NCBIXML.read(h)
        except Exception as e:                      # noqa
            last = e
            print(f"    BLAST attempt {attempt}/{RETRIES} failed: {e}")
            time.sleep(10 * attempt)
    raise RuntimeError(f"BLAST failed after {RETRIES} attempts: {last}")


def fetch_pam(acc, start, end, strand):
    from Bio import Entrez
    Entrez.email = EMAIL
    last = None
    for attempt in range(1, RETRIES + 1):
        try:
            h = Entrez.efetch(db="nuccore", id=acc, rettype="fasta", retmode="text",
                              seq_start=start, seq_stop=end, strand=strand)
            txt = h.read()
            h.close()
            seq = "".join(l.strip() for l in txt.splitlines() if not l.startswith(">")).upper()
            if len(seq) != 3:
                raise ValueError(f"expected 3 nt PAM, got {seq!r}")
            time.sleep(PAUSE_BETWEEN_FETCHES)
            return seq
        except Exception as e:                      # noqa
            last = e
            print(f"    PAM fetch attempt {attempt}/{RETRIES} failed: {e}")
            time.sleep(5 * attempt)
    raise RuntimeError(f"PAM fetch failed: {last}")


def versioned_accession(aln):
    m = re.search(r"([A-Z]{1,2}_?\d+\.\d+)", getattr(aln, "hit_id", "") or "")
    return m.group(1) if m else aln.accession


# ------------------------------------------------------------------------ main
def select_guides(run_all):
    with open(VERDICT_CSV) as fh:
        verdict = {r["guide"]: r["verdict"] for r in csv.DictReader(fh)}
    prog = {}
    with open(PROGRESS_CSV) as fh:
        for r in csv.DictReader(fh):
            prog[r["guide"]] = r                   # last row per guide wins (errors get retried)
    if run_all:
        return list(verdict), verdict
    flagged = [g for g in verdict
               if prog.get(g, {}).get("status") == "done"
               and int(prog[g]["raw_hits"]) > 1 and int(prog[g]["genomic_hits"]) == 0]
    return flagged, verdict


def process_guide(guide, gene_seq):
    rec = blast_guide(guide)
    rows, cands = [], []
    for aln in rec.alignments:
        acc = versioned_accession(aln)
        for hsp in aln.hsps:
            cat, mm, reason = classify_hit(acc, hsp.sbjct, hsp.align_length,
                                           hsp.identities, hsp.gaps, gene_seq)
            rows.append(dict(guide=guide, accession=acc, title=aln.title[:90],
                             sbjct_start=hsp.sbjct_start, sbjct_end=hsp.sbjct_end,
                             align_len=hsp.align_length, mismatches=mm,
                             category=cat, reason=reason, pam="", note=""))
            if cat == "CANDIDATE":
                cands.append(dict(acc=acc, lo=min(hsp.sbjct_start, hsp.sbjct_end),
                                  hi=max(hsp.sbjct_start, hsp.sbjct_end), mm=mm,
                                  start=hsp.sbjct_start, end=hsp.sbjct_end))
    loci = collapse_loci(cands)
    unverified = False
    any_pam = False
    for L in loci:
        s, e, strand = pam_fetch_params(L["start"], L["end"])
        try:
            pam = fetch_pam(L["acc"], s, e, strand)
            L["pam"] = pam
            L["has_pam"] = pam[1:] == "GG"
        except Exception as ex:                     # LOUD, never swallowed
            unverified = True
            L["pam"] = "?"
            L["has_pam"] = False
            print(f"    !! could not verify PAM for {L['acc']}: {ex}")
        any_pam |= L["has_pam"]
        for r in rows:
            if r["category"] == "CANDIDATE" and r["accession"] == L["acc"] \
                    and L["lo"] <= max(r["sbjct_start"], r["sbjct_end"]) \
                    and min(r["sbjct_start"], r["sbjct_end"]) <= L["hi"]:
                r["pam"] = L["pam"]
    if unverified:
        verdict = "UNVERIFIED"
    elif any_pam:
        verdict = "AT RISK"
    elif loci:
        verdict = "Low risk"
    else:
        verdict = "Clean"
    excl_close = sum(1 for r in rows if r["category"] == "EXCL_NONGENOMIC"
                     and "REVIEW" in r["reason"])
    counts = {}
    for r in rows:
        counts[r["category"]] = counts.get(r["category"], 0) + 1
    return verdict, len(rows), counts, loci, excl_close, rows


def main():
    print("cwd:", os.getcwd())
    print("files here:", sorted(f for f in os.listdir(".") if f.lower().startswith("g6pd")))
    for f in (GENE_FASTA, PROGRESS_CSV, VERDICT_CSV):
        if not os.path.exists(f):
            sys.exit(f"MISSING {f} in this folder -> you are in the wrong working folder.")
    gene_seq = read_fasta_seq(GENE_FASTA)
    print(f"{GENE_NAME} gene length: {len(gene_seq)} bp")
    guides, old = select_guides("--all" in sys.argv)
    print(f"Re-checking {len(guides)} guide(s)\n")

    hit_fields = ["guide", "accession", "title", "sbjct_start", "sbjct_end", "align_len",
                  "mismatches", "category", "reason", "pam", "note"]
    ver_fields = ["guide", "old_verdict", "new_verdict", "changed", "hsps_total",
                  "n_self", "n_nonself_kept_loci", "n_excluded_close_review", "category_counts"]
    with open(OUT_HITS, "w", newline="") as fh_h, open(OUT_VERDICT, "w", newline="") as fh_v:
        wh = csv.DictWriter(fh_h, fieldnames=hit_fields)
        wv = csv.DictWriter(fh_v, fieldnames=ver_fields)
        wh.writeheader()
        wv.writeheader()
        for i, g in enumerate(guides, 1):
            print(f"[{i}/{len(guides)}] {g}  (old verdict: {old[g]})")
            try:
                verdict, n, counts, loci, excl, rows = process_guide(g, gene_seq)
            except Exception as ex:
                print(f"    !! FAILED, guide left UNVERIFIED: {ex}")
                wv.writerow(dict(guide=g, old_verdict=old[g], new_verdict="UNVERIFIED",
                                 changed="?", category_counts=str(ex)))
                fh_v.flush()
                continue
            for r in rows:
                wh.writerow(r)
            changed = "YES" if verdict != old[g] else "no"
            print(f"    HSPs={n}  {counts}")
            print(f"    NEW verdict: {verdict}  (changed: {changed})  "
                  f"non-self loci: {len(loci)}  excluded-but-close (review): {excl}\n")
            wv.writerow(dict(guide=g, old_verdict=old[g], new_verdict=verdict, changed=changed,
                             hsps_total=n, n_self=counts.get("SELF", 0),
                             n_nonself_kept_loci=len(loci), n_excluded_close_review=excl,
                             category_counts=str(counts)))
            fh_h.flush()
            fh_v.flush()
            time.sleep(PAUSE_BETWEEN_SEARCHES)
    print(f"Done. Wrote {os.path.abspath(OUT_HITS)} and {os.path.abspath(OUT_VERDICT)}")
    print("Send me g6pd_retry_verdict.csv (and g6pd_retry_hits.csv).")


if __name__ == "__main__":
    main()
    