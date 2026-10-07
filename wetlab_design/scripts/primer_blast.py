"""
primer_blast.py  --  Phase 4, Step 3b of CRISPR-GuideDesigner

Checks the 20 designed primers (primer_order.csv + amplicons.csv) for
specificity against the human genome (GRCh38 chromosomes only).

How it works:
  1. Sends the primers to NCBI BLAST in small batches (short-query settings:
     blastn, word size 7, E-value 1000), human GRCh38 only. Each batch is
     retried if NCBI errors (e.g. HTTP 502) and its raw reply is saved to
     primer_blast_raw_<n>.xml -- if the run stops, just run it again and it
     resumes from the batches already saved.
  2. Recognises each primer's INTENDED site from the coordinates in
     amplicons.csv (sequence-verified: the intended hit must be found).
  3. For every other hit it counts mismatches over the WHOLE primer (unaligned
     ends and gaps count as mismatches) and checks whether the primer's
     3' end (last 3 bases) matches -- the 3' end is what polymerase extends.
  4. Verdict per primer:
       AT RISK   a non-intended site with <= 3 mismatches AND a matching 3' end
       Low risk  non-intended sites with <= 5 mismatches, none of the above
       Clean     no non-intended site within 5 mismatches
       CHECK     the intended site was not recovered (settings problem)
  5. Looks for unintended PCR products: any two possible binding sites
     (<= 4 mismatches, matching 3' end) on opposite strands of the same
     chromosome, facing each other, <= 3000 bp apart, other than the
     intended pair.

Needs internet. Requires: biopython. BLAST on NCBI can take several minutes.
Saves: primer_blast_report.csv, amplicon_check.csv, primer_blast_raw_<n>.xml
"""

import csv
import os
import time

# ---- Settings ----------------------------------------------------------
PRIMER_FILE = "primer_order.csv"
AMPLICON_FILE = "amplicons.csv"
BATCH_SIZE = 5                             # primers per BLAST request
RETRIES, WAIT = 3, 45                      # tries per database setting, base wait (s)
RISK_MM, LOW_MM, RISK_3P = 3, 5, 3         # mismatch limits, 3' bases that must match
PRODUCT_MM, MAX_PRODUCT = 4, 3000          # unintended-product check
SLOP = 5                                   # bp tolerance when recognising the intended site
CHROM_ACC = {                              # GRCh38 chromosome accessions
    "chr1": "NC_000001.11", "chr2": "NC_000002.12", "chr3": "NC_000003.12",
    "chr4": "NC_000004.12", "chr5": "NC_000005.10", "chr6": "NC_000006.12",
    "chr7": "NC_000007.14", "chr8": "NC_000008.11", "chr9": "NC_000009.12",
    "chr10": "NC_000010.11", "chr11": "NC_000011.10", "chr12": "NC_000012.12",
    "chr13": "NC_000013.11", "chr14": "NC_000014.9", "chr15": "NC_000015.10",
    "chr16": "NC_000016.10", "chr17": "NC_000017.11", "chr18": "NC_000018.10",
    "chr19": "NC_000019.10", "chr20": "NC_000020.11", "chr21": "NC_000021.9",
    "chr22": "NC_000022.11", "chrX": "NC_000023.11", "chrY": "NC_000024.10",
}
# -------------------------------------------------------------------------


def three_prime_match(hsp, n, k=RISK_3P):
    """True if the alignment reaches the primer's 3' end and its last k bases match."""
    if hsp.query_end != n:
        return False
    counted = 0
    for q, m in zip(reversed(hsp.query), reversed(hsp.match)):
        if q == "-":
            continue
        if m != "|":
            return False
        counted += 1
        if counted == k:
            return True
    return False


def hsp_site(aln_id, hsp):
    """One binding site: chromosome, strand the primer matches, 5' end, 3' end (1-based)."""
    lead = hsp.query_start - 1                    # primer bases left unaligned at the 5' end
    if hsp.sbjct_start <= hsp.sbjct_end:          # primer = plus-strand sequence
        return dict(acc=aln_id, strand="+", five=hsp.sbjct_start - lead,
                    lo=hsp.sbjct_start, hi=hsp.sbjct_end)
    return dict(acc=aln_id, strand="-", five=hsp.sbjct_start + lead,
                lo=hsp.sbjct_end, hi=hsp.sbjct_start)


def analyse_primer(record, name, seq, target_acc, t_lo, t_hi, t_strand):
    """Classify every HSP of one primer's BLAST record."""
    n = len(seq)
    target_found, others, seen = False, [], set()
    for aln in record.alignments:
        for hsp in aln.hsps:
            site = hsp_site(aln.hit_id, hsp)
            mm = n - hsp.identities
            site.update(mm=mm, p3=three_prime_match(hsp, n))
            is_target = (target_acc in aln.hit_id and site["strand"] == t_strand
                         and site["lo"] <= t_hi + SLOP and site["hi"] >= t_lo - SLOP)
            if is_target:
                if mm == 0:
                    target_found = True
                site["target"] = True
            else:
                site["target"] = False
                key = (aln.hit_id, site["strand"], site["five"])
                if key in seen:
                    continue
                seen.add(key)
                others.append(site)
    risky = [s for s in others if s["mm"] <= RISK_MM and s["p3"]]
    low = [s for s in others if s["mm"] <= LOW_MM]
    if not target_found:
        verdict = "CHECK (intended site not recovered)"
    elif risky:
        verdict = "AT RISK"
    elif low:
        verdict = "Low risk"
    else:
        verdict = "Clean"
    worst = min(low, key=lambda s: (not s["p3"], s["mm"])) if low else None
    return dict(name=name, seq=seq, target_found=target_found, verdict=verdict,
                n_close_sites=len(low), n_risky_sites=len(risky), worst=worst,
                sites=others)


def sites_for_products(res, target_site):
    """Sites that could prime: the intended one plus close, 3'-matched off-targets."""
    out = [dict(target_site, target=True)]
    out += [s for s in res["sites"] if s["mm"] <= PRODUCT_MM and s["p3"]]
    return out


def unintended_products(sites):
    """Facing plus/minus sites on one chromosome <= MAX_PRODUCT apart, not the intended pair."""
    found = []
    plus = [s for s in sites if s["strand"] == "+"]
    minus = [s for s in sites if s["strand"] == "-"]
    for p in plus:
        for m in minus:
            if p["acc"] != m["acc"]:
                continue
            size = m["five"] - p["five"] + 1
            if 0 < size <= MAX_PRODUCT and not (p["target"] and m["target"]):
                found.append((p["acc"], p["five"], m["five"], size))
    return found


def blast_settings():
    """Database / filter combinations to try, in order."""
    acc_query = " OR ".join(f"{a}[Accession]" for a in CHROM_ACC.values())
    return [("refseq_genomic", acc_query),
            ("refseq_representative_genomes", "Homo sapiens[Organism]")]


def run_batch(n, batch):
    """BLAST one batch of (name, seq); retry, then fall back to the next setting."""
    from Bio.Blast import NCBIWWW
    fasta = "".join(f">{nm}\n{sq}\n" for nm, sq in batch)
    path = f"primer_blast_raw_{n}.xml"
    for db, query in blast_settings():
        for attempt in range(1, RETRIES + 1):
            try:
                print(f"  batch {n}: {db}, attempt {attempt} ...")
                handle = NCBIWWW.qblast(
                    "blastn", db, fasta, entrez_query=query,
                    word_size=7, expect=1000, nucl_reward=1, nucl_penalty=-3,
                    gapcosts="5 2", hitlist_size=100, megablast=False)
                xml = handle.read()
                handle.close()
                if "<BlastOutput" not in xml:
                    raise ValueError("reply was not BLAST XML")
                with open(path, "w") as f:
                    f.write(xml)
                return
            except Exception as e:
                print(f"    failed: {type(e).__name__}: {e}")
                if attempt < RETRIES:
                    print(f"    waiting {WAIT * attempt}s before retrying ...")
                    time.sleep(WAIT * attempt)
    raise SystemExit(f"NCBI BLAST kept failing on batch {n}. Batches already saved are kept; "
                     f"run again later (NCBI is busiest during US daytime).")


def run_all_batches(primers):
    files = []
    for n, i in enumerate(range(0, len(primers), BATCH_SIZE), start=1):
        path = f"primer_blast_raw_{n}.xml"
        if os.path.exists(path):
            print(f"  batch {n}: already saved, skipping")
        else:
            run_batch(n, primers[i:i + BATCH_SIZE])
        files.append(path)
    return files


def main():
    from Bio.Blast import NCBIXML
    primers = [(r["primer_name"], r["sequence_5to3"].strip().upper())
               for r in csv.DictReader(open(PRIMER_FILE))]
    amps = [a for a in csv.DictReader(open(AMPLICON_FILE)) if a.get("fwd_seq")]

    # intended site of every primer, from the amplicon coordinates
    intended = {}
    for a in amps:
        chrom_acc = CHROM_ACC[a["chrom"]]
        s, e = int(a["amplicon_chr_start"]), int(a["amplicon_chr_end"])
        f, r = a["fwd_seq"], a["rev_seq"]
        intended[a["fwd_name"]] = dict(acc=chrom_acc, lo=s, hi=s + len(f) - 1, strand="+", five=s)
        intended[a["rev_name"]] = dict(acc=chrom_acc, lo=e - len(r) + 1, hi=e, strand="-", five=e)

    print("Running BLAST in batches (each can take a few minutes) ...")
    files = run_all_batches(primers)
    records = []
    for path in files:
        with open(path) as fh:
            records += list(NCBIXML.parse(fh))
    if len(records) != len(primers):
        raise SystemExit(f"Expected {len(primers)} BLAST results, got {len(records)}. "
                         f"Delete the primer_blast_raw_*.xml files and run again.")

    results = {}
    for (name, seq), rec in zip(primers, records):
        t = intended[name]
        results[name] = analyse_primer(rec, name, seq, t["acc"], t["lo"], t["hi"], t["strand"])

    rows = []
    for name, seq in primers:
        r = results[name]
        w = r["worst"]
        rows.append(dict(
            primer_name=name, sequence_5to3=seq, verdict=r["verdict"],
            intended_site_found=r["target_found"],
            close_offtarget_sites=r["n_close_sites"], risky_offtarget_sites=r["n_risky_sites"],
            worst_mismatches=w["mm"] if w else "", worst_3prime_match=w["p3"] if w else "",
            worst_location=f"{w['acc']}:{w['lo']}-{w['hi']} ({w['strand']})" if w else ""))
    with open("primer_blast_report.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)

    pair_rows, done = [], set()
    for a in amps:
        key = (a["fwd_name"], a["rev_name"])
        if key in done:
            continue
        done.add(key)
        sites = (sites_for_products(results[a["fwd_name"]], intended[a["fwd_name"]])
                 + sites_for_products(results[a["rev_name"]], intended[a["rev_name"]]))
        extra = unintended_products(sites)
        pair_rows.append(dict(
            fwd=a["fwd_name"], rev=a["rev_name"], intended_product_bp=a["product_bp"],
            fwd_verdict=results[a["fwd_name"]]["verdict"],
            rev_verdict=results[a["rev_name"]]["verdict"],
            unintended_products=len(extra),
            details="; ".join(f"{acc}:{p5}-{m5} ({sz} bp)" for acc, p5, m5, sz in extra) or "none"))
    with open("amplicon_check.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(pair_rows[0].keys()))
        wr.writeheader()
        wr.writerows(pair_rows)

    print("\nPrimer verdicts:")
    for r in rows:
        print(f"  {r['primer_name']:<16} {r['verdict']}")
    print("\nAmplicon checks:")
    for p in pair_rows:
        print(f"  {p['fwd']} / {p['rev']}: unintended products = {p['unintended_products']}")
    print("\nSaved primer_blast_report.csv and amplicon_check.csv")


if __name__ == "__main__":
    main()
