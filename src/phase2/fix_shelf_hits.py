"""
fix_self_hits.py  --  Phase 2, Step 4 patch for CRISPR-GuideDesigner

WHY THIS EXISTS
----------------
blast_offtargets_v_2.py's is_self_hit() only recognised PIP4K2C's own
CHROMOSOME record (NC_000012) as "self". It missed NG_125907, an NCBI
"candidate regulatory element" record that -- we confirmed by direct
sequence comparison -- contains PIP4K2C's own exon 1 verbatim, just under
a different accession and a mislabelled coordinate range. That let 3
guides get wrongly flagged "AT RISK".

WHAT THIS SCRIPT DOES DIFFERENTLY
----------------------------------
Instead of matching accessions or coordinates (which is what let
NG_125907 slip through), this checks the SEQUENCE ITSELF:
for every hit blast_offtargets_v_2.py already found, it fetches that
hit's actual DNA (a small, fast NCBI call -- no BLAST re-run needed) and
asks: does this sequence appear anywhere in PIP4K2C's own gene, on
either strand? If yes, it's a self-match, regardless of what accession
or coordinates it was reported under.

This reuses the raw hits blast_offtargets_v_2.py already collected
(pip4k2c_blast_hits_genomic.csv) -- it does NOT resubmit anything to
BLAST, so it's fast and needs no queue wait.

Files needed in the same folder: pip4k2c_gene.fasta,
pip4k2c_blast_hits_genomic.csv, pip4k2c_shortlist.csv
Requires: biopython, pandas
"""

import os
import time

import pandas as pd
from Bio import Entrez, SeqIO

Entrez.email = "your_email@example.com"   # <-- CHANGE THIS

# ---- Settings -----------------------------------------------------------
PAUSE_BETWEEN_FETCHES = 1   # seconds, courtesy pause between Entrez calls

GENE_FASTA = "pip4k2c_gene.fasta"
RAW_HITS_FILE = "pip4k2c_blast_hits_genomic.csv"
SHORTLIST_FILE = "pip4k2c_shortlist.csv"

FIXED_HITS_FILE = "pip4k2c_blast_hits_genomic_fixed.csv"
FIXED_LOCUS_FILE = "pip4k2c_unique_loci_fixed.csv"
FIXED_VERDICT_FILE = "pip4k2c_offtarget_verdict_fixed.csv"
# ---------------------------------------------------------------------------


def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def fetch_hit_sequence(accession, hit_start, hit_end):
    """Download the exact stretch of DNA a hit covers, in its own reported orientation."""
    lo, hi = sorted((int(hit_start), int(hit_end)))
    strand = 1 if hit_start < hit_end else 2   # 2 = minus; Entrez returns the revcomp for us
    handle = Entrez.efetch(
        db="nuccore", id=accession, rettype="fasta", retmode="text",
        seq_start=lo, seq_stop=hi, strand=strand,
    )
    try:
        return str(SeqIO.read(handle, "fasta").seq).upper()
    finally:
        handle.close()


def is_true_self_match(hit_seq, gene_seq, gene_seq_rc):
    """
    True if a hit's own sequence is actually part of PIP4K2C's own gene,
    on either strand -- regardless of what accession it came from.
    """
    return hit_seq in gene_seq or hit_seq in gene_seq_rc


def rebuild_self_hit_flags(raw_hits, gene_seq):
    """
    Re-check every raw hit's real sequence against the gene. Returns the
    hits with a new `is_true_self` column, fetching only what's needed
    (skips hits already at >0 mismatches from a DIFFERENT accession that
    obviously can't be a verbatim self-match... actually we check all of
    them, since accession/mismatch count told us nothing reliable before).
    """
    gene_seq = gene_seq.upper()
    gene_seq_rc = reverse_complement(gene_seq)

    flags, errors = [], []
    for i, row in enumerate(raw_hits.itertuples(), start=1):
        try:
            hit_seq = fetch_hit_sequence(row.accession, row.hit_start, row.hit_end)
            is_self = is_true_self_match(hit_seq, gene_seq, gene_seq_rc)
            err = ""
        except Exception as e:
            hit_seq, is_self, err = "", False, str(e)[:200]
            print(f"    [{i}/{len(raw_hits)}] could not re-check "
                  f"{row.accession}:{row.hit_start}-{row.hit_end} ({err})")
        flags.append(is_self)
        errors.append(err)
        time.sleep(PAUSE_BETWEEN_FETCHES)

    out = raw_hits.copy()
    out["is_true_self"] = flags
    out["recheck_error"] = errors
    return out


def intervals_overlap(s1, e1, s2, e2):
    return max(s1, s2) <= min(e1, e2)


def collapse_overlapping_loci(hits):
    """Same collapsing rule as before: merge overlapping hits from the same accession."""
    if hits.empty:
        return pd.DataFrame(columns=[
            "guide", "locus_id", "accession", "hit_title", "locus_start",
            "locus_end", "hit_strand", "best_mismatches", "pam", "has_pam",
            "duplicate_record_count",
        ])

    loci = []
    for guide, g_hits in hits.groupby("guide"):
        for accession, a_hits in g_hits.groupby("accession"):
            a_hits = a_hits.sort_values(
                by=[a_hits[["hit_start", "hit_end"]].min(axis=1).name] if False else "hit_start"
            )
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
    """Same verdict rule as before: AT RISK > Low risk > Clean, per guide."""
    rows = []
    for guide in shortlist["guide"]:
        g_loci = unique_loci[unique_loci["guide"] == guide]
        pam_count = int(g_loci["has_pam"].sum())
        locus_count = len(g_loci)
        if pam_count > 0:
            verdict = "AT RISK"
            details = f"{pam_count} unique genomic locus/loci with canonical NGG PAM."
        elif locus_count > 0:
            verdict = "Low risk"
            details = f"{locus_count} unique genomic sequence-match locus/loci, but none with canonical NGG PAM."
        else:
            verdict = "Clean"
            details = "No qualifying non-self genomic locus detected."
        rows.append({"guide": guide, "verdict": verdict, "details": details,
                    "genomic_hits": locus_count, "unique_loci": locus_count,
                    "pam_positive_loci": pam_count})
    return pd.DataFrame(rows)


if __name__ == "__main__":
    print(f"Working folder: {os.getcwd()}")

    gene_seq = str(SeqIO.read(GENE_FASTA, "fasta").seq)
    raw_hits = pd.read_csv(RAW_HITS_FILE)
    shortlist = pd.read_csv(SHORTLIST_FILE)

    print(f"Re-checking {len(raw_hits)} previously found genomic hits "
         f"against PIP4K2C's own sequence (no BLAST re-run)...\n")
    rechecked = rebuild_self_hit_flags(raw_hits, gene_seq)

    n_newly_self = int(rechecked["is_true_self"].sum())
    print(f"\nHits that are actually the gene's own sequence, under a "
         f"different accession: {n_newly_self}")
    if n_newly_self:
        print(rechecked.loc[rechecked["is_true_self"],
                            ["guide", "accession", "hit_title"]].to_string(index=False))

    genuine_hits = rechecked[~rechecked["is_true_self"]].drop(
        columns=["is_true_self", "recheck_error"])
    unique_loci = collapse_overlapping_loci(genuine_hits)
    verdict = classify(shortlist, unique_loci)

    print("\nCorrected verdict counts:", verdict["verdict"].value_counts().to_dict())
    print("\nGuides that are not simply 'Clean':\n")
    print(verdict[verdict["verdict"] != "Clean"].to_string(index=False))

    rechecked.to_csv(FIXED_HITS_FILE, index=False)
    unique_loci.to_csv(FIXED_LOCUS_FILE, index=False)
    verdict.to_csv(FIXED_VERDICT_FILE, index=False)
    print("\nSaved these files:")
    for name in (FIXED_HITS_FILE, FIXED_LOCUS_FILE, FIXED_VERDICT_FILE):
        print("  " + os.path.abspath(name))
        