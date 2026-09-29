"""
kras_blast_offtargets.py  --  Case study 4, Step 4 of CRISPR-GuideDesigner

Off-target screen for the human KRAS shortlist, using NCBI BLAST online.
Built from the start with everything learned from PIP4K2C's Phase 2:

  1. GENOMIC-ONLY search: only genomic accessions are kept (NC_/NG_/NT_/
     NW_/AC_). Transcript/mRNA/cDNA records (NM_/XM_/XR_/NR_/BC_/AK_) are
     excluded up front, so the gene's own alternate transcript records
     never show up as false "off-targets" in the first place.
  2. SEQUENCE-VERIFIED self-hit exclusion: every remaining hit's own DNA is
     fetched and compared directly against KRAS's gene sequence (both
     strands). If the hit's sequence is genuinely part of the gene itself
     -- regardless of what accession or coordinates it was filed under --
     it is marked self and excluded. This is what PIP4K2C needed a second
     patch script (fix_self_hits.py) to add after the fact; here it is
     built in from the first run.
  3. Overlapping hits from the same accession are collapsed into one
     unique locus, so a single real off-target site isn't double-counted.
  4. Each surviving locus is checked for an adjacent NGG PAM (Cas9 can't
     cut without one).
  5. Each guide gets one verdict: Clean / Low risk / AT RISK.

This step is SLOW (each BLAST search can take 30s-3min) and needs a steady
internet connection. It is resumable: if it stops partway through, run it
again and it will skip guides already marked "done".

Files needed in the same folder: kras_gene.fasta, kras_shortlist.csv
Requires: biopython, pandas
"""

import os
import time

import pandas as pd
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

Entrez.email = "your_email@example.com"   # <-- CHANGE THIS

# ---- Settings you can change -------------------------------------------
DATABASE = "nt"
ORGANISM_FILTER = "Homo sapiens[Organism]"
WORD_SIZE = 7
EXPECT = 1000
HITLIST_SIZE = 50
MAX_MISMATCH = 4
MAX_HITS_CHECKED = 8   # how many closest-matching hits per guide to fully verify
PAUSE_BETWEEN_SEARCHES = 3
PAUSE_BETWEEN_FETCHES = 1

# Only these RefSeq/genomic prefixes are considered at all. Transcript and
# protein prefixes (NM_, XM_, XR_, NR_, NP_, XP_, BC_, AK_...) are excluded
# before any self-hit or PAM logic even runs.
GENOMIC_PREFIXES = ("NC_", "NG_", "NT_", "NW_", "AC_")

GENE_FASTA = "kras_gene.fasta"
SHORTLIST_FILE = "kras_shortlist.csv"
HITS_FILE = "kras_blast_hits_genomic.csv"
PROGRESS_FILE = "kras_blast_progress.csv"
LOCUS_FILE = "kras_unique_loci.csv"
VERDICT_FILE = "kras_offtarget_verdict.csv"
# -------------------------------------------------------------------------

HIT_COLUMNS = ["guide", "accession", "hit_title", "hit_start", "hit_end",
               "hit_strand", "mismatches", "pam", "has_pam"]
PROGRESS_COLUMNS = ["guide", "status", "raw_hits", "genomic_hits", "note"]


def is_genomic_accession(accession):
    return str(accession).strip().upper().startswith(GENOMIC_PREFIXES)


def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def fetch_sequence(accession, start, end):
    """Download a stretch of DNA in its own reported orientation."""
    lo, hi = sorted((int(start), int(end)))
    strand = 1 if start < end else 2
    handle = Entrez.efetch(
        db="nuccore", id=accession, rettype="fasta", retmode="text",
        seq_start=lo, seq_stop=hi, strand=strand,
    )
    try:
        return str(SeqIO.read(handle, "fasta").seq).upper()
    finally:
        handle.close()


def pam_fetch_params(hit_start, hit_end):
    """
    Where to look for the PAM next to a BLAST hit, and which strand to read
    it on. Returns (fetch_start, fetch_end, entrez_strand), 1-based genomic
    coordinates. entrez_strand: 1 = plus, 2 = minus (Entrez returns the
    reverse complement for us). Uses hit_end as the reference in both
    orientations, since BLAST's own convention already encodes the guide's
    3' side there for both strands (verified against Phase 2's PIP4K2C run).
    """
    if hit_start < hit_end:                # plus-strand hit
        return hit_end + 1, hit_end + 3, 1
    else:                                   # minus-strand hit
        return hit_end - 3, hit_end - 1, 2


def fetch_pam(accession, hit_start, hit_end):
    """Download the 3 bases next to a hit and return them (e.g. 'TGG')."""
    fs, fe, strand = pam_fetch_params(hit_start, hit_end)
    if fs < 1:
        return None
    return fetch_sequence_stranded(accession, fs, fe, strand)


def fetch_sequence_stranded(accession, fs, fe, strand):
    """Like fetch_sequence, but lets the caller force the Entrez strand (1/2)."""
    handle = Entrez.efetch(
        db="nuccore", id=accession, rettype="fasta", retmode="text",
        seq_start=fs, seq_stop=fe, strand=strand,
    )
    try:
        return str(SeqIO.read(handle, "fasta").seq).upper()
    finally:
        handle.close()


def is_true_self_match(hit_seq, gene_seq, gene_seq_rc):
    """True if a hit's own sequence genuinely is part of KRAS's own gene."""
    return bool(hit_seq) and (hit_seq in gene_seq or hit_seq in gene_seq_rc)


def hsp_mismatches(hsp, guide_len):
    if hsp.align_length != guide_len or hsp.gaps:
        return None
    return hsp.align_length - hsp.identities


def blast_one_guide(guide):
    result_handle = NCBIWWW.qblast(
        "blastn", DATABASE, guide,
        entrez_query=ORGANISM_FILTER,
        word_size=WORD_SIZE, expect=EXPECT, hitlist_size=HITLIST_SIZE,
    )
    record = NCBIXML.read(result_handle)
    result_handle.close()
    return record


def find_genomic_hits(guide, record, gene_seq, gene_seq_rc):
    """
    Walk a BLAST record: keep only genomic-accession hits within
    MAX_MISMATCH, verify each one's own sequence against the gene to drop
    real self-matches, and fetch the adjacent PAM for what remains.
    Returns (raw_hit_count, list_of_hit_dicts).
    """
    candidates = []
    raw_hits = 0
    for alignment in record.alignments:
        accession = alignment.accession
        if not is_genomic_accession(accession):
            continue
        title = alignment.hit_def
        for hsp in alignment.hsps:
            raw_hits += 1
            mm = hsp_mismatches(hsp, len(guide))
            if mm is None or mm > MAX_MISMATCH:
                continue
            candidates.append({
                "accession": accession, "hit_title": title,
                "hit_start": hsp.sbjct_start, "hit_end": hsp.sbjct_end,
                "hit_strand": "+" if hsp.sbjct_start < hsp.sbjct_end else "-",
                "mismatches": mm,
            })

    candidates.sort(key=lambda h: h["mismatches"])
    kept = []
    for hit in candidates[:MAX_HITS_CHECKED]:
        try:
            hit_seq = fetch_sequence(hit["accession"], hit["hit_start"], hit["hit_end"])
        except Exception as e:
            print(f"    (could not fetch {hit['accession']}: {e})")
            continue
        time.sleep(PAUSE_BETWEEN_FETCHES)
        if is_true_self_match(hit_seq, gene_seq, gene_seq_rc):
            continue   # confirmed self-match -- not an off-target

        try:
            pam = fetch_pam(hit["accession"], hit["hit_start"], hit["hit_end"])
        except Exception:
            pam = None
        time.sleep(PAUSE_BETWEEN_FETCHES)
        has_pam = bool(pam) and pam[1:] == "GG"
        kept.append({"guide": guide, **hit, "pam": pam, "has_pam": has_pam})
    return raw_hits, kept


def intervals_overlap(s1, e1, s2, e2):
    return max(s1, s2) <= min(e1, e2)


def collapse_overlapping_loci(hits):
    if hits.empty:
        return pd.DataFrame(columns=[
            "guide", "locus_id", "accession", "hit_title", "locus_start",
            "locus_end", "hit_strand", "best_mismatches", "pam", "has_pam",
            "duplicate_record_count",
        ])
    loci = []
    for guide, g_hits in hits.groupby("guide"):
        for accession, a_hits in g_hits.groupby("accession"):
            a_hits = a_hits.assign(
                lo=a_hits[["hit_start", "hit_end"]].min(axis=1),
                hi=a_hits[["hit_start", "hit_end"]].max(axis=1),
            ).sort_values("lo")
            current = None
            for row in a_hits.itertuples():
                if current is None:
                    current = {
                        "guide": guide, "accession": accession, "hit_title": row.hit_title,
                        "locus_start": row.lo, "locus_end": row.hi, "hit_strand": row.hit_strand,
                        "best_mismatches": row.mismatches, "pam": row.pam, "has_pam": row.has_pam,
                        "duplicate_record_count": 1,
                    }
                    continue
                if intervals_overlap(current["locus_start"], current["locus_end"], row.lo, row.hi):
                    current["locus_start"] = min(current["locus_start"], row.lo)
                    current["locus_end"] = max(current["locus_end"], row.hi)
                    current["duplicate_record_count"] += 1
                    if row.mismatches < current["best_mismatches"]:
                        current.update(best_mismatches=row.mismatches, pam=row.pam,
                                      has_pam=row.has_pam, hit_strand=row.hit_strand)
                else:
                    loci.append(current)
                    current = {
                        "guide": guide, "accession": accession, "hit_title": row.hit_title,
                        "locus_start": row.lo, "locus_end": row.hi, "hit_strand": row.hit_strand,
                        "best_mismatches": row.mismatches, "pam": row.pam, "has_pam": row.has_pam,
                        "duplicate_record_count": 1,
                    }
            if current is not None:
                loci.append(current)
    df = pd.DataFrame(loci)
    df["locus_id"] = df["accession"] + ":" + df["locus_start"].astype(str) + "-" + df["locus_end"].astype(str)
    return df


def classify(shortlist, unique_loci):
    rows = []
    for guide in shortlist["guide"]:
        g_loci = unique_loci[unique_loci["guide"] == guide]
        pam_count = int(g_loci["has_pam"].sum())
        locus_count = len(g_loci)
        if pam_count > 0:
            verdict, details = "AT RISK", f"{pam_count} unique genomic locus/loci with canonical NGG PAM."
        elif locus_count > 0:
            verdict, details = "Low risk", f"{locus_count} unique genomic locus/loci, none with canonical NGG PAM."
        else:
            verdict, details = "Clean", "No qualifying non-self genomic locus detected."
        rows.append({"guide": guide, "verdict": verdict, "details": details,
                    "genomic_hits": locus_count, "unique_loci": locus_count,
                    "pam_positive_loci": pam_count})
    return pd.DataFrame(rows)


def load_progress():
    if os.path.exists(PROGRESS_FILE):
        return pd.read_csv(PROGRESS_FILE)
    return pd.DataFrame(columns=PROGRESS_COLUMNS)


def load_hits():
    if os.path.exists(HITS_FILE):
        return pd.read_csv(HITS_FILE)
    return pd.DataFrame(columns=HIT_COLUMNS)


if __name__ == "__main__":
    print(f"Working folder: {os.getcwd()}")

    gene_seq = str(SeqIO.read(GENE_FASTA, "fasta").seq).upper()
    gene_seq_rc = reverse_complement(gene_seq)
    shortlist = pd.read_csv(SHORTLIST_FILE)
    progress = load_progress()
    hits = load_hits()
    done_guides = set(progress.loc[progress["status"] == "done", "guide"])

    print(f"{len(shortlist)} guides in shortlist, {len(done_guides)} already done.\n")

    for i, row in enumerate(shortlist.itertuples(), start=1):
        guide = row.guide
        if guide in done_guides:
            print(f"[{i}/{len(shortlist)}] {guide}  already done, skipping")
            continue

        print(f"[{i}/{len(shortlist)}] {guide}  submitting to NCBI BLAST...")
        try:
            record = blast_one_guide(guide)
            raw_hits, genomic_hits = find_genomic_hits(guide, record, gene_seq, gene_seq_rc)
            status, note = "done", ""
            print(f"    {raw_hits} raw hits -> {len(genomic_hits)} genomic, "
                  f"non-self off-target candidate(s) "
                  f"({sum(h['has_pam'] for h in genomic_hits)} with a PAM)")
        except Exception as e:
            raw_hits, genomic_hits = 0, []
            status, note = "error", str(e)[:200]
            print(f"    ERROR: {note}")

        if genomic_hits:
            hits = pd.concat([hits, pd.DataFrame(genomic_hits)], ignore_index=True)
            hits.to_csv(HITS_FILE, index=False)

        progress = pd.concat([progress, pd.DataFrame([{
            "guide": guide, "status": status, "raw_hits": raw_hits,
            "genomic_hits": len(genomic_hits), "note": note,
        }])], ignore_index=True)
        progress.to_csv(PROGRESS_FILE, index=False)

        time.sleep(PAUSE_BETWEEN_SEARCHES)

    unique_loci = collapse_overlapping_loci(hits)
    verdict = classify(shortlist, unique_loci)

    n_done = int((progress["status"] == "done").sum())
    n_error = int((progress["status"] == "error").sum())
    print(f"\nFinished: {n_done} guides screened, {n_error} errors.")
    if n_error:
        print("Re-run this script to retry the guides that errored.")
    print("\nVerdict counts:", verdict["verdict"].value_counts().to_dict())
    print("\nGuides that are not simply 'Clean':\n")
    print(verdict[verdict["verdict"] != "Clean"].to_string(index=False))

    unique_loci.to_csv(LOCUS_FILE, index=False)
    verdict.to_csv(VERDICT_FILE, index=False)
    print("\nSaved these files:")
    for name in (HITS_FILE, PROGRESS_FILE, LOCUS_FILE, VERDICT_FILE):
        print("  " + os.path.abspath(name))
        