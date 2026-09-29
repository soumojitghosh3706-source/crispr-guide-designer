"""
kras_compare_crispor.py  --  Phase 3 of CRISPR-GuideDesigner (KRAS)

Compares this project's own results for all 9 shortlisted KRAS guides
against CRISPOR (crispor.tefor.net, hg38), read from its results table by
hand across three small sequence windows (one per exon cluster: gene
positions 5400-5800, 23400-23700, 25050-25400).

Requires: pandas
"""

import pandas as pd

# ---- Our own results (from kras_shortlist.csv / kras_offtarget_verdict.csv) ----
OURS = [
    {"guide": "GTAGTTGGAGCTGGTGGCGT", "exon": 1, "our_score": 80.0,  "our_verdict": "Clean"},
    {"guide": "CTGAATTAGCTGTATCGTCA", "exon": 1, "our_score": 80.0,  "our_verdict": "Clean"},
    {"guide": "TGATGGAGAAACCTGTCTCT", "exon": 2, "our_score": 90.0,  "our_verdict": "Clean"},
    {"guide": "TCTCGACACAGCAGGTCAAG", "exon": 2, "our_score": 90.0,  "our_verdict": "Clean"},
    {"guide": "CAATGAGGGACCAGTACATG", "exon": 2, "our_score": 100.0, "our_verdict": "Clean"},
    {"guide": "TCCCTTCTCAGGATTCCTAC", "exon": 2, "our_score": 88.0,  "our_verdict": "Clean"},
    {"guide": "GGACTCTGAAGATGTACCTA", "exon": 3, "our_score": 90.0,  "our_verdict": "AT RISK"},
    {"guide": "AACATCAGCAAAGACAAGAC", "exon": 3, "our_score": 54.6,  "our_verdict": "Clean"},
    {"guide": "AGTAGACACAAAACAGGCTC", "exon": 3, "our_score": 84.6,  "our_verdict": "Clean"},
]

# ---- CRISPOR's results (hg38), read from its results table by hand ----
# nras_hit: True if a paralog (NRAS) exon appears explicitly among the
# top genome-browser-link annotations shown for that guide's off-targets.
CRISPOR = [
    {"guide": "GTAGTTGGAGCTGGTGGCGT", "mit_specificity": 54, "cfd_specificity": 75, "off_targets": "1-0-2-21-277", "flag": "",           "nras_hit": False},
    {"guide": "CTGAATTAGCTGTATCGTCA", "mit_specificity": 75, "cfd_specificity": 93, "off_targets": "1-0-0-14-62",  "flag": "",           "nras_hit": False},
    {"guide": "TGATGGAGAAACCTGTCTCT", "mit_specificity": 38, "cfd_specificity": 79, "off_targets": "1-0-6-37-211", "flag": "Inefficient", "nras_hit": False},
    {"guide": "TCTCGACACAGCAGGTCAAG", "mit_specificity": 80, "cfd_specificity": 90, "off_targets": "0-0-1-12-89",  "flag": "",           "nras_hit": False},
    {"guide": "CAATGAGGGACCAGTACATG", "mit_specificity": 79, "cfd_specificity": 91, "off_targets": "0-0-1-13-116","flag": "",           "nras_hit": True},
    {"guide": "TCCCTTCTCAGGATTCCTAC", "mit_specificity": 76, "cfd_specificity": 84, "off_targets": "0-0-3-17-151", "flag": "",           "nras_hit": False},
    {"guide": "GGACTCTGAAGATGTACCTA", "mit_specificity": 67, "cfd_specificity": 81, "off_targets": "0-1-3-28-189", "flag": "",           "nras_hit": True},
    {"guide": "AACATCAGCAAAGACAAGAC", "mit_specificity": 53, "cfd_specificity": 71, "off_targets": "0-0-1-40-353", "flag": "",           "nras_hit": False},
    {"guide": "AGTAGACACAAAACAGGCTC", "mit_specificity": 53, "cfd_specificity": 81, "off_targets": "0-0-7-103-213","flag": "",           "nras_hit": False},
]

MIT_SPECIFICITY_MIN = 50


def close_mismatch_hits(off_targets):
    """True if CRISPOR found any hit at 1 or 2 mismatches (the closest, riskiest kind)."""
    counts = [int(x) for x in off_targets.split("-")]
    return sum(counts[1:3]) > 0


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
    confirmed = table[(table["our_verdict"] == "AT RISK")]
    nras = table[table["nras_hit"]]

    print(f"{len(agree)} of {len(table)} guides validated cleanly against CRISPOR.")
    print(f"{len(nras)} guide(s) show a direct NRAS exon hit in CRISPOR's own results "
          f"(confirming the paralog hypothesis independently of our BLAST screen's "
          f"KRASP1-pseudogene finding).\n")
    print(table[["guide", "exon", "our_verdict", "mit_specificity", "cfd_specificity",
                "off_targets", "nras_hit", "comparison"]].to_string(index=False))

    lines = [
        "## Phase 3: validation against CRISPOR (KRAS)",
        "",
        "All 9 shortlisted KRAS guides were run through CRISPOR (crispor.tefor.net, "
        "hg38) across three small windows covering the gene's three exon clusters.",
        "",
        "| Guide | Exon | Our verdict | CRISPOR MIT spec. | CRISPOR CFD spec. | "
        "CRISPOR off-targets (0-1-2-3-4mm) | NRAS exon hit? | Result |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in table.itertuples():
        lines.append(f"| `{r.guide}` | {r.exon} | {r.our_verdict} | {r.mit_specificity} | "
                     f"{r.cfd_specificity} | {r.off_targets} | "
                     f"{'Yes' if r.nras_hit else ''} | {r.comparison} |")
    lines += [
        "",
        f"**Only {len(agree)} of {len(table)} guides ({100*len(agree)//len(table)}%) "
        f"validated cleanly** -- consistent with the same systematic BLAST-sensitivity "
        f"gap found on MUTYH (3/20) and PIP4K2C (3/5): our `word_size=7` BLAST search "
        f"misses off-targets whose 1-2 mismatches are spread too thinly for a 7-bp "
        f"perfect seed to exist.",
        "",
        "**The paralog hypothesis is independently confirmed by CRISPOR itself:** "
        "2 of the 9 guides (`CAATGAGGGACCAGTACATG`, and `GGACTCTGAAGATGTACCTA` -- the "
        "same guide our own BLAST screen already flagged AT RISK for a KRASP1 "
        "pseudogene hit) show an off-target site landing directly inside an **NRAS** "
        "exon in CRISPOR's own results. This means KRAS's off-target risk comes from "
        "*both* directions predicted: a near-identical pseudogene copy of itself "
        "(KRASP1, found by our own screen) AND its true functional paralog (NRAS, "
        "found only by CRISPOR's more sensitive search, missed by ours).",
        "",
        "**Final recommended KRAS guide** (the only one validated on both methods): "
        + ", ".join(f"`{g}`" for g in agree["guide"]) + ".",
    ]
    with open("kras_phase3_summary.md", "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    table.to_csv("kras_phase3_comparison.csv", index=False)
    print("\nSaved kras_phase3_comparison.csv and kras_phase3_summary.md")
    