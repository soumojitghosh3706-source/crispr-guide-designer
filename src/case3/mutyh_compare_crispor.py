"""
mutyh_compare_crispor.py  --  Phase 3 of CRISPR-GuideDesigner (MUTYH)

Compares this project's own results for all 20 shortlisted MUTYH guides
against CRISPOR (crispor.tefor.net, hg38), read from its results table by
hand (no batch API suited to a tablet workflow) across two sequence
windows covering the shortlist (gene positions 6900-8400 and 8300-9400).

Same structure as phase3_compare_crispor.py (PIP4K2C), reused directly.

Requires: pandas
"""

import pandas as pd

# ---- Our own results (from mutyh_guides_scored.csv / mutyh_offtarget_verdict.csv) ----
OURS = [
    {"guide": "AAGATGAGATGGACCTGGAC", "exon": 4, "our_score": 100.0, "our_verdict": "Clean"},
    {"guide": "TTGGTTGAAATCTCCTGGCC", "exon": 10, "our_score": 100.0, "our_verdict": "Clean"},
    {"guide": "CACAGGAGGTGAATCAACTC", "exon": 7, "our_score": 96.0, "our_verdict": "Clean"},
    {"guide": "ATCAGGTGGTAGAGGAGCTA", "exon": 8, "our_score": 96.0, "our_verdict": "Clean"},
    {"guide": "AGGAACAGCTCTTAGCCTCA", "exon": 11, "our_score": 93.6, "our_verdict": "Clean"},
    {"guide": "TGGGCTACTATTCTCGTGGC", "exon": 7, "our_score": 90.0, "our_verdict": "Clean"},
    {"guide": "GCCAAAGGCGATAGAGGCAA", "exon": 8, "our_score": 90.0, "our_verdict": "Clean"},
    {"guide": "TGGTGGATGGCAACGTAGCA", "exon": 9, "our_score": 90.0, "our_verdict": "Clean"},
    {"guide": "TGGGATCAGCACCAATGGCT", "exon": 9, "our_score": 90.0, "our_verdict": "Clean"},
    {"guide": "CCAGAGCTGCTGGGAAACAA", "exon": 9, "our_score": 90.0, "our_verdict": "Clean"},
    {"guide": "AGATTTCAACCAAGCAGCCA", "exon": 10, "our_score": 89.1, "our_verdict": "Clean"},
    {"guide": "AGTGGTCAACTTCCCCAGAA", "exon": 12, "our_score": 87.7, "our_verdict": "Clean"},
    {"guide": "TCACCACACTCCTCCACGTC", "exon": 11, "our_score": 69.6, "our_verdict": "Clean"},
    {"guide": "CCTCTGCACCAGCAGAATTT", "exon": 12, "our_score": 83.8, "our_verdict": "Clean"},
    {"guide": "TGCAGGGTCTCTGCTGTACG", "exon": 8, "our_score": 80.0, "our_verdict": "Clean"},
    {"guide": "AGGCTGTTCCAGAACACAGG", "exon": 12, "our_score": 75.6, "our_verdict": "Clean"},
    {"guide": "TCAAGTATATGGGCTGGCCT", "exon": 14, "our_score": 73.3, "our_verdict": "Clean"},
    {"guide": "AGGTGGCACTGTCCAGTGTT", "exon": 12, "our_score": 71.9, "our_verdict": "Clean"},
    {"guide": "TTGCGCTGAAGCTGCTCTGA", "exon": 13, "our_score": 70.1, "our_verdict": "Clean"},
    {"guide": "GCTTGATGTGAGAGAAGGTG", "exon": 14, "our_score": 70.0, "our_verdict": "Clean"},
]

# ---- CRISPOR's results (hg38), read from its results table by hand ----
CRISPOR = [
    {"guide": "AAGATGAGATGGACCTGGAC", "mit_specificity": 64, "cfd_specificity": 83, "off_targets": "0-1-4-24-208", "flag": ""},
    {"guide": "TTGGTTGAAATCTCCTGGCC", "mit_specificity": 63, "cfd_specificity": 81, "off_targets": "0-0-3-15-147", "flag": "Inefficient"},
    {"guide": "CACAGGAGGTGAATCAACTC", "mit_specificity": 73, "cfd_specificity": 92, "off_targets": "0-0-1-18-187", "flag": ""},
    {"guide": "ATCAGGTGGTAGAGGAGCTA", "mit_specificity": 66, "cfd_specificity": 80, "off_targets": "0-0-2-24-159", "flag": ""},
    {"guide": "AGGAACAGCTCTTAGCCTCA", "mit_specificity": 68, "cfd_specificity": 81, "off_targets": "0-0-6-32-175", "flag": ""},
    {"guide": "TGGGCTACTATTCTCGTGGC", "mit_specificity": 90, "cfd_specificity": 97, "off_targets": "0-0-0-5-63",   "flag": ""},
    {"guide": "GCCAAAGGCGATAGAGGCAA", "mit_specificity": 76, "cfd_specificity": 81, "off_targets": "0-0-2-10-130", "flag": ""},
    {"guide": "TGGTGGATGGCAACGTAGCA", "mit_specificity": 88, "cfd_specificity": 92, "off_targets": "0-0-0-8-89",   "flag": ""},
    {"guide": "TGGGATCAGCACCAATGGCT", "mit_specificity": 80, "cfd_specificity": 88, "off_targets": "0-0-2-14-129", "flag": ""},
    {"guide": "CCAGAGCTGCTGGGAAACAA", "mit_specificity": 55, "cfd_specificity": 74, "off_targets": "0-0-6-48-307", "flag": ""},
    {"guide": "AGATTTCAACCAAGCAGCCA", "mit_specificity": 64, "cfd_specificity": 85, "off_targets": "0-0-3-36-161", "flag": ""},
    {"guide": "AGTGGTCAACTTCCCCAGAA", "mit_specificity": 68, "cfd_specificity": 83, "off_targets": "0-0-2-17-160", "flag": ""},
    {"guide": "TCACCACACTCCTCCACGTC", "mit_specificity": 82, "cfd_specificity": 93, "off_targets": "0-0-1-6-106",  "flag": ""},
    {"guide": "CCTCTGCACCAGCAGAATTT", "mit_specificity": 72, "cfd_specificity": 89, "off_targets": "0-0-5-23-173", "flag": "Inefficient"},
    {"guide": "TGCAGGGTCTCTGCTGTACG", "mit_specificity": 87, "cfd_specificity": 92, "off_targets": "0-0-0-8-99",   "flag": ""},
    {"guide": "AGGCTGTTCCAGAACACAGG", "mit_specificity": 61, "cfd_specificity": 81, "off_targets": "0-0-7-24-205", "flag": ""},
    {"guide": "TCAAGTATATGGGCTGGCCT", "mit_specificity": 71, "cfd_specificity": 87, "off_targets": "0-0-1-20-124", "flag": "Inefficient"},
    {"guide": "AGGTGGCACTGTCCAGTGTT", "mit_specificity": 75, "cfd_specificity": 90, "off_targets": "0-0-4-12-138", "flag": "Inefficient"},
    {"guide": "TTGCGCTGAAGCTGCTCTGA", "mit_specificity": 68, "cfd_specificity": 84, "off_targets": "0-0-1-16-116", "flag": ""},
    {"guide": "GCTTGATGTGAGAGAAGGTG", "mit_specificity": 61, "cfd_specificity": 76, "off_targets": "0-0-1-25-242", "flag": ""},
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
    print(f"{len(agree)} of {len(table)} guides validated cleanly against CRISPOR.\n")
    print(table[["guide", "exon", "mit_specificity", "cfd_specificity",
                "off_targets", "comparison"]].to_string(index=False))

    lines = [
        "## Phase 3: validation against CRISPOR (MUTYH)",
        "",
        "All 20 shortlisted MUTYH guides -- our own pipeline had called all 20 "
        "\"Clean\" -- were run through CRISPOR (crispor.tefor.net, hg38) as a "
        "full cross-check.",
        "",
        "| Guide | Exon | CRISPOR MIT spec. | CRISPOR CFD spec. | "
        "CRISPOR off-targets (0-1-2-3-4mm) | Result |",
        "|---|---|---|---|---|---|",
    ]
    for r in table.itertuples():
        lines.append(f"| `{r.guide}` | {r.exon} | {r.mit_specificity} | "
                     f"{r.cfd_specificity} | {r.off_targets} | {r.comparison} |")
    lines += [
        "",
        f"**Only {len(agree)} of {len(table)} guides ({100*len(agree)//len(table)}%) "
        f"validated cleanly.** The other {len(table) - len(agree)} were flagged by "
        f"CRISPOR with off-target sites at 1-2 mismatches that this project's own "
        f"genomic BLAST screen missed entirely -- a systematic pattern, not an "
        f"isolated case (PIP4K2C showed the same gap on 2 of 5 checked guides). "
        f"The likely cause: our BLAST search uses `word_size=7`, so a guide with "
        f"2 mismatches spread evenly through its 20 nt often contains no single "
        f"7-bp perfectly-matching seed, and BLAST never registers it as a "
        f"candidate at all -- while CRISPOR's purpose-built aligner has no such "
        f"blind spot. This is recorded as the project's central methodological "
        f"limitation (see Limitations).",
        "",
        "**Final recommended MUTYH guides** (validated on both methods): "
        + ", ".join(f"`{g}`" for g in agree["guide"]) + ".",
    ]
    with open("mutyh_phase3_summary.md", "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    table.to_csv("mutyh_phase3_comparison.csv", index=False)
    print("\nSaved mutyh_phase3_comparison.csv and mutyh_phase3_summary.md")
    