"""
phase3_compare_crispor.py  --  Phase 3 of CRISPR-GuideDesigner

Compares this project's own results for PIP4K2C's 5 top "Clean" guides
against CRISPOR (crispor.tefor.net, hg38), to check whether an independent,
published tool agrees with this pipeline's picks.

CRISPOR's web tool can't be automated from a tablet (no public batch API for
this use case), so its results below were read off its results page by hand
and are hardcoded here. This script's job is to line them up against our
own numbers, compute a plain-language verdict per guide, and write the
comparison table used in the README.

Requires: pandas
"""

import pandas as pd

# ---- Our own results (from pip4k2c_guides_scored.csv / pip4k2c_offtarget_verdict_fixed.csv) ----
OURS = [
    {"guide": "CAATGCCAAATCGATCACGG", "exon": 3, "our_score": 100.0, "our_verdict": "Clean"},
    {"guide": "CTTCGTTGTCCACACTGACT", "exon": 5, "our_score": 100.0, "our_verdict": "Clean"},
    {"guide": "ATGCTTCTTCTTGGTCTTGG", "exon": 1, "our_score": 90.0,  "our_verdict": "Clean"},
    {"guide": "TCAGATCAATGAGCTCAGCC", "exon": 2, "our_score": 97.3,  "our_verdict": "Clean"},
    {"guide": "GAGCTGGCCTTAAAGTCATC", "exon": 2, "our_score": 100.0, "our_verdict": "Clean"},
]

# ---- CRISPOR's results (hg38), read from its results table by hand ----
# off_targets = off-target site counts at 0/1/2/3/4 mismatches (CRISPOR's own column)
CRISPOR = [
    {"guide": "CAATGCCAAATCGATCACGG", "mit_specificity": 95, "cfd_specificity": 97,
     "off_targets": "0-0-0-2-27", "flag": ""},
    {"guide": "CTTCGTTGTCCACACTGACT", "mit_specificity": 79, "cfd_specificity": 92,
     "off_targets": "0-0-0-15-77", "flag": ""},
    {"guide": "ATGCTTCTTCTTGGTCTTGG", "mit_specificity": 68, "cfd_specificity": 71,
     "off_targets": "0-0-0-19-239", "flag": ""},
    {"guide": "TCAGATCAATGAGCTCAGCC", "mit_specificity": 67, "cfd_specificity": 84,
     "off_targets": "0-0-2-18-168", "flag": "Inefficient"},
    {"guide": "GAGCTGGCCTTAAAGTCATC", "mit_specificity": 44, "cfd_specificity": 90,
     "off_targets": "0-1-1-14-79", "flag": ""},
]

# Below this MIT specificity score, or any hit at 1 or 2 mismatches, a guide
# is treated as a genuine cross-tool disagreement worth flagging in the README.
MIT_SPECIFICITY_MIN = 50


def close_mismatch_hits(off_targets):
    """True if CRISPOR found any hit at 1 or 2 mismatches (the closest, riskiest kind)."""
    counts = [int(x) for x in off_targets.split("-")]
    return sum(counts[1:3]) > 0   # positions 1 and 2 in the 0-1-2-3-4 breakdown


def verdict(row):
    if close_mismatch_hits(row["off_targets"]) or row["mit_specificity"] < MIT_SPECIFICITY_MIN:
        return "Disagreement -- CRISPOR finds closer off-targets than our screen did"
    if row["flag"]:
        return f"Agreement on off-targets; CRISPOR flags: {row['flag']}"
    return "Agreement"


if __name__ == "__main__":
    ours = pd.DataFrame(OURS)
    crispor = pd.DataFrame(CRISPOR)
    table = ours.merge(crispor, on="guide")
    table["comparison"] = table.apply(verdict, axis=1)
    table = table.sort_values("mit_specificity", ascending=False).reset_index(drop=True)

    agree = table[table["comparison"] == "Agreement"]
    print(f"{len(agree)} of {len(table)} guides validated cleanly against CRISPOR.\n")
    print(table[["guide", "exon", "our_score", "mit_specificity", "cfd_specificity",
                "off_targets", "comparison"]].to_string(index=False))

    lines = [
        "## Phase 3: validation against CRISPOR",
        "",
        f"The top 5 \"Clean\" PIP4K2C guides were run through CRISPOR "
        f"(crispor.tefor.net, hg38) to check this pipeline's results against "
        f"an independent, published tool.",
        "",
        "| Guide | Exon | Our score | CRISPOR MIT spec. | CRISPOR CFD spec. | "
        "CRISPOR off-targets (0-1-2-3-4mm) | Result |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in table.itertuples():
        lines.append(f"| `{r.guide}` | {r.exon} | {r.our_score:.1f} | "
                     f"{r.mit_specificity} | {r.cfd_specificity} | {r.off_targets} | "
                     f"{r.comparison} |")
    lines += [
        "",
        f"**{len(agree)} of 5 guides validated cleanly.** The other "
        f"{len(table) - len(agree)} were flagged by CRISPOR with off-target "
        f"sites at 1-2 mismatches that this project's genomic BLAST screen "
        f"missed -- most likely because CRISPOR uses a purpose-built, "
        f"exhaustive genome index, while this project's screen relies on "
        f"NCBI's general-purpose BLAST search, which is less sensitive at "
        f"finding very-close near-matches. This is now recorded as a known "
        f"limitation of the off-target screen (see Limitations).",
        "",
        "**Final recommended guides** (validated on both methods): "
        + ", ".join(f"`{g}`" for g in agree["guide"]) + ".",
    ]
    with open("phase3_summary.md", "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    table.to_csv("phase3_comparison.csv", index=False)
    print("\nSaved phase3_comparison.csv and phase3_summary.md")
  
