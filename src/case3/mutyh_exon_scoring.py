"""
mutyh_exon_scoring.py  --  Case study 3, Step 3 of CRISPR-GuideDesigner

Same logic as exon_scoring.py for PIP4K2C, retargeted to MUTYH.

What it does:
  1. Reads mutyh_guides_coding.csv (from mutyh_exon_guides.py) and mutyh_exons.csv.
  2. Applies filters. A guide must pass ALL of them:
       - GC content 40-60%
       - no TTTT or GGGG run
       - cuts an exon shared by all curated isoforms
       - is not in an NMD-escape zone
       - cut is at least 5 bases from an exon edge
  3. Scores the survivors 0-100: 40% GC closeness, 40% position in the
     coding sequence, 20% distance from the exon edge.
  4. Picks a shortlist of 20 guides with at most 4 per exon and cut sites
     at least 20 bp apart, so the shortlist is spread across the gene.
  5. Saves mutyh_guides_scored.csv and mutyh_shortlist.csv.

Files needed in the same folder: scoring.py (Phase 1), mutyh_guides_coding.csv,
mutyh_exons.csv
Requires: pandas, biopython
"""

import os

import pandas as pd

from scoring import gc_fraction, gc_score, position_score, GC_MIN, GC_MAX, BAD_MOTIFS

# ---- Settings you can change -------------------------------------------
REQUIRE_SHARED = True
EXCLUDE_NMD_RISK = True
MIN_EDGE_DIST = 5
EDGE_FULL_SCORE_DIST = 15
SHORTLIST_SIZE = 20
MAX_PER_EXON = 4
MIN_SPACING = 20
# -------------------------------------------------------------------------


def score_coding_guides(df, cds_len):
    df = df.copy()
    df["gc"] = df["guide"].apply(gc_fraction).round(3)
    df["gc_ok"] = df["gc"].between(GC_MIN, GC_MAX)
    df["motif_ok"] = ~df["guide"].apply(lambda g: any(m in g for m in BAD_MOTIFS))
    df["shared_ok"] = df["shared_by_curated"] if REQUIRE_SHARED else True
    df["nmd_ok"] = ~df["nmd_risk"] if EXCLUDE_NMD_RISK else True
    df["edge_ok"] = df["dist_to_exon_edge"] >= MIN_EDGE_DIST
    flags = ["gc_ok", "motif_ok", "shared_ok", "nmd_ok", "edge_ok"]
    df["pass_filters"] = df[flags].all(axis=1)

    def row_score(r):
        edge = min(r["dist_to_exon_edge"] / EDGE_FULL_SCORE_DIST, 1.0)
        return round(100 * (0.4 * gc_score(r["gc"])
                            + 0.4 * position_score(r["cds_pos"], cds_len)
                            + 0.2 * edge), 1)

    df["score"] = df.apply(row_score, axis=1)
    return df


def pick_shortlist(passed, n=SHORTLIST_SIZE, max_per_exon=MAX_PER_EXON,
                   min_spacing=MIN_SPACING):
    chosen, per_exon = [], {}
    for row in passed.itertuples():
        if per_exon.get(row.coding_exon, 0) >= max_per_exon:
            continue
        if any(abs(row.cut_pos - c.cut_pos) < min_spacing for c in chosen):
            continue
        chosen.append(row)
        per_exon[row.coding_exon] = per_exon.get(row.coding_exon, 0) + 1
        if len(chosen) == n:
            break
    return pd.DataFrame(chosen).drop(columns="Index", errors="ignore")


if __name__ == "__main__":
    coding = pd.read_csv("mutyh_guides_coding.csv")
    exons = pd.read_csv("mutyh_exons.csv")
    cds_len = int(exons["cds_end"].max())

    scored = score_coding_guides(coding, cds_len)

    print(f"Guides cutting coding exons:   {len(scored)}")
    print(f"Fail GC (outside 40-60%):      {int((~scored['gc_ok']).sum())}")
    print(f"Fail TTTT/GGGG motif:          {int((~scored['motif_ok']).sum())}")
    print(f"Fail: exon not shared:         {int((~scored['shared_ok']).sum())}")
    print(f"Fail: NMD-escape zone:         {int((~scored['nmd_ok']).sum())}")
    print(f"Fail: too close to exon edge:  {int((~scored['edge_ok']).sum())}")
    print("(a guide can fail several rules, so these counts overlap)")
    print(f"Pass all filters:              {int(scored['pass_filters'].sum())}")

    passed = scored[scored["pass_filters"]].sort_values(
        ["score", "cds_pos"], ascending=[False, True]
    )
    shortlist = pick_shortlist(passed).reset_index(drop=True)

    print(f"\nShortlist: {len(shortlist)} guides. Top 10:\n")
    cols = ["guide", "pam", "strand", "coding_exon", "cds_pos", "gc",
            "dist_to_exon_edge", "score"]
    print(shortlist[cols].head(10).to_string(index=False))
    print("\nShortlist guides per coding exon:",
          shortlist["coding_exon"].value_counts().sort_index().to_dict())

    scored.to_csv("mutyh_guides_scored.csv", index=False)
    shortlist.to_csv("mutyh_shortlist.csv", index=False)
    print("\nSaved these files:")
    for name in ("mutyh_guides_scored.csv", "mutyh_shortlist.csv"):
        print("  " + os.path.abspath(name))
        