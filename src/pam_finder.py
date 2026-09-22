"""
pam_finder.py  --  Step 2 of CRISPR-GuideDesigner (Phase 1)

What it does:
  1. Reads lacZ.fasta (created by fetch.py).
  2. Scans BOTH strands for NGG PAM sites (SpCas9).
  3. For every PAM, takes the 20 nt immediately 5' of it as the guide.
  4. Saves all guides to guides_all.csv.

Coordinates are positions inside the gene (1-based), where position 1 is
the A of the ATG start codon. The '-' strand is the template strand.

Requires: biopython, pandas
"""

import pandas as pd
from Bio import SeqIO

GUIDE_LEN = 20
PAM = "NGG"          # SpCas9. The guide sits directly 5' of this PAM.

# Which bases each PAM letter allows (N = any base).
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T",
         "N": "ACGT", "R": "AG", "V": "ACG"}


def reverse_complement(seq):
    """Return the reverse complement of a plain A/C/G/T string."""
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def pam_matches(candidate, pattern=PAM):
    """True if the 3-letter candidate fits the PAM pattern (e.g. NGG)."""
    return all(base in IUPAC[p] for base, p in zip(candidate, pattern))


def scan(seq, pam=PAM, guide_len=GUIDE_LEN):
    """
    Scan ONE strand (a plain string). For every PAM that has at least
    guide_len bases in front of it, return (guide, pam_seq, pam_index),
    where pam_index is the 0-based position where the PAM starts.
    """
    hits = []
    for i in range(guide_len, len(seq) - len(pam) + 1):
        pam_seq = seq[i:i + len(pam)]
        if pam_matches(pam_seq, pam):
            guide = seq[i - guide_len:i]
            if set(guide) <= set("ACGT"):      # skip guides containing N etc.
                hits.append((guide, pam_seq, i))
    return hits


def find_guides(gene_seq, pam=PAM, guide_len=GUIDE_LEN):
    """Find guides on both strands. Returns a pandas DataFrame."""
    gene_seq = gene_seq.upper()
    L = len(gene_seq)
    rows = []

    # --- Forward strand (the gene as written) ---
    for guide, pam_seq, i in scan(gene_seq, pam, guide_len):
        start = i - guide_len + 1            # 1-based
        end = i
        rows.append({
            "guide": guide, "pam": pam_seq, "strand": "+",
            "start": start, "end": end,
            "cut_pos": end - 3,              # Cas9 cuts 3 bp upstream of the PAM
        })

    # --- Reverse strand: scan the reverse complement, then convert back ---
    rc = reverse_complement(gene_seq)
    for guide, pam_seq, j in scan(rc, pam, guide_len):
        start = L - j + 1                    # position in ORIGINAL gene coordinates
        end = L - j + guide_len
        rows.append({
            "guide": guide, "pam": pam_seq, "strand": "-",
            "start": start, "end": end,
            "cut_pos": start + 2,
        })

    df = pd.DataFrame(rows).sort_values("start").reset_index(drop=True)
    return df


if __name__ == "__main__":
    record = SeqIO.read("lacZ.fasta", "fasta")
    gene_seq = str(record.seq)

    df = find_guides(gene_seq)

    n_plus = int((df["strand"] == "+").sum())
    n_minus = int((df["strand"] == "-").sum())

    print(f"Gene length: {len(gene_seq)} bp")
    print(f"Guides found: {len(df)}  (+ strand: {n_plus}, - strand: {n_minus})")
    print(f"Rough expectation for ~50% GC: about {len(gene_seq) // 8}")
    print()
    print(df.head(10).to_string(index=False))

    df.to_csv("guides_all.csv", index=False)
    print("\nSaved guides_all.csv")
    
