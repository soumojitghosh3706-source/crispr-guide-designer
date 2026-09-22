"""
scoring.py  --  Step 3 of CRISPR-GuideDesigner (Phase 1)

What it does:
  1. Reads guides_all.csv (from pam_finder.py) and lacZ.fasta.
  2. Calculates the GC content of every guide.
  3. Applies filters:  GC between 40-60%, and no TTTT or GGGG runs.
  4. Gives every guide a simple 0-100 heuristic score.
  5. Picks a shortlist of the best guides whose cut sites are spread out.
  6. Saves guides_scored.csv (all guides) and guides_shortlist.csv (best ones).

NOTE: the score is a transparent rule-of-thumb, NOT a validated efficiency
model. In Phase 3 we compare our picks against CRISPOR.

Requires: biopython, pandas
"""

import pandas as pd
from Bio import SeqIO

try:
    from Bio.SeqUtils import gc_fraction          # Biopython 1.80 or newer
except ImportError:                                # older Biopython versions
    def gc_fraction(seq):
        seq = str(seq).upper()
        return (seq.count("G") + seq.count("C")) / len(seq)

# ---- Settings you can change -------------------------------------------
GC_MIN, GC_MAX = 0.40, 0.60       # allowed GC content of the 20-nt guide
BAD_MOTIFS = ["TTTT", "GGGG"]     # TTTT stops U6-promoter transcription; GGGG folds badly
SHORTLIST_SIZE = 20               # how many guides to send to the off-target step
MIN_SPACING = 20                  # min distance (bp) between shortlisted cut sites
# -------------------------------------------------------------------------


def gc_score(gc):
    """1.0 at 50% GC, falling to 0.5 at 40% or 60%."""
    return max(0.0, 1 - abs(gc - 0.5) / 0.2)


def position_score(cut_pos, gene_len):
    """
    Cuts in the first half of the gene are preferred, because an early
    disruption is more likely to knock out the whole gene.
    Cuts in the very first 5% get a reduced score.
    """
    frac = cut_pos / gene_len
    if frac < 0.05:
        return 0.5
    if frac <= 0.50:
        return 1.0
    return max(0.0, 1 - (frac - 0.5) / 0.5)


def score_guides(df, gene_len):
    """Add GC, filter flags, and a 0-100 score to the guide table."""
    df = df.copy()
    df["gc"] = df["guide"].apply(gc_fraction).round(3)
    df["gc_ok"] = df["gc"].between(GC_MIN, GC_MAX)
    df["motif_ok"] = ~df["guide"].apply(lambda g: any(m in g for m in BAD_MOTIFS))
    df["pass_filters"] = df["gc_ok"] & df["motif_ok"]
    df["pos_frac"] = (df["cut_pos"] / gene_len).round(3)
    df["score"] = df.apply(
        lambda r: round(
            100 * (0.6 * gc_score(r["gc"]) + 0.4 * position_score(r["cut_pos"], gene_len)),
            1,
        ),
        axis=1,
    )
    return df


def pick_shortlist(passed, n=SHORTLIST_SIZE, min_spacing=MIN_SPACING):
    """
    Go down the ranked list and keep a guide only if its cut site is at least
    min_spacing bp away from every guide already chosen. This avoids a
    shortlist made of near-identical, overlapping guides.
    """
    chosen = []
    for row in passed.itertuples():
        if all(abs(row.cut_pos - c.cut_pos) >= min_spacing for c in chosen):
            chosen.append(row)
        if len(chosen) == n:
            break
    return pd.DataFrame(chosen).drop(columns="Index", errors="ignore")


if __name__ == "__main__":
    gene_len = len(SeqIO.read("lacZ.fasta", "fasta").seq)
    df = pd.read_csv("guides_all.csv")

    scored = score_guides(df, gene_len)

    print(f"Total guides:             {len(scored)}")
    print(f"Fail GC (outside 40-60%): {int((~scored['gc_ok']).sum())}")
    print(f"Fail TTTT/GGGG motif:     {int((~scored['motif_ok']).sum())}")
    print(f"Pass all filters:         {int(scored['pass_filters'].sum())}")
    print("(a guide can fail both, so the fail counts overlap)")

    passed = scored[scored["pass_filters"]].sort_values(
        ["score", "start"], ascending=[False, True]
    )
    shortlist = pick_shortlist(passed).reset_index(drop=True)

    print(f"\nShortlist: {len(shortlist)} guides. Top 10:\n")
    cols = ["guide", "pam", "strand", "start", "end", "gc", "pos_frac", "score"]
    print(shortlist[cols].head(10).to_string(index=False))

    scored.to_csv("guides_scored.csv", index=False)
    shortlist.to_csv("guides_shortlist.csv", index=False)
    print("\nSaved guides_scored.csv and guides_shortlist.csv")
    
