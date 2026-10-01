"""
g6pd_compare_crispor.py  --  Phase 3 of CRISPR-GuideDesigner (G6PD)

Compares this project's own results for all 20 shortlisted G6PD guides
against CRISPOR (crispor.tefor.net, hg38), read from its results table by
hand across five sequence windows (gene positions 1350-1650, 11300-11750,
12100-12500, 13000-13600, 13850-14550).

Requires: pandas
"""

import pandas as pd

# ---- Our own results (from g6pd_shortlist.csv / g6pd_offtarget_verdict.csv) ----
OURS = [
    {"guide": "GTGTGTATCCGACTGATGGA", "exon": 1, "our_score": 100.0, "our_verdict": "Low risk"},
    {"guide": "GGTAGATCTTCTTCTTGGCC", "exon": 2, "our_score": 94.7, "our_verdict": "Clean"},
    {"guide": "GACACACTTACCAGATGGTG", "exon": 2, "our_score": 88.0, "our_verdict": "Clean"},
    {"guide": "TTTCGGGCAGAAGGCCATCC", "exon": 3, "our_score": 80.0, "our_verdict": "Clean"},
    {"guide": "CCCGAAAACACCTTCATCGT", "exon": 3, "our_score": 100.0, "our_verdict": "Clean"},
    {"guide": "AGAAGGGCTCACTCTGTTTG", "exon": 3, "our_score": 100.0, "our_verdict": "Clean"},
    {"guide": "AGAGGAGAAGCTCAAGCTGG", "exon": 4, "our_score": 90.0, "our_verdict": "Clean"},
    {"guide": "CTTTGCCCGCAACTCCTATG", "exon": 4, "our_score": 90.0, "our_verdict": "Clean"},
    {"guide": "GGGCATTCATGTGGCTGTTG", "exon": 4, "our_score": 90.0, "our_verdict": "Clean"},
    {"guide": "GCAGGACTCGTGAATGTTCT", "exon": 4, "our_score": 100.0, "our_verdict": "Clean"},
    {"guide": "ATCATCGTGGAGAAGCCCTT", "exon": 5, "our_score": 100.0, "our_verdict": "Clean"},
    {"guide": "AGGAGATGTGGTTGGACAGC", "exon": 5, "our_score": 90.0, "our_verdict": "Clean"},
    {"guide": "GATCTGGTCCTCACGGAACA", "exon": 5, "our_score": 90.0, "our_verdict": "Clean"},
    {"guide": "TACCGCATCGACCACTACCT", "exon": 5, "our_score": 90.0, "our_verdict": "Clean"},
    {"guide": "ATGTTGTCCCGGTTCCAGAT", "exon": 6, "our_score": 100.0, "our_verdict": "Clean"},
    {"guide": "CTTGAAGGTGAGGATAACGC", "exon": 6, "our_score": 100.0, "our_verdict": "Clean"},
    {"guide": "CTCACCTTCAAGGAGCCCTT", "exon": 6, "our_score": 90.0, "our_verdict": "Clean"},
    {"guide": "GACACAGCATCTGCAGTAGG", "exon": 7, "our_score": 89.1, "our_verdict": "Clean"},
    {"guide": "ATCACGGACGTCATCTGAGT", "exon": 7, "our_score": 96.5, "our_verdict": "Clean"},
    {"guide": "GAAATGCATCTCAGAGGTGC", "exon": 8, "our_score": 93.9, "our_verdict": "Clean"},
]

# ---- CRISPOR's results (hg38), read from its results table by hand ----
# variant: True if CRISPOR's "+ Variants" display shows a known sequence
# variant under the guide (G6PD is highly polymorphic, so this matters:
# a variant can stop the guide from cutting in carriers).
CRISPOR = [
    {"guide": "GTGTGTATCCGACTGATGGA", "mit_specificity": 85, "cfd_specificity": 92, "off_targets": "0-0-0-5-72", "flag": "", "variant": False},
    {"guide": "GGTAGATCTTCTTCTTGGCC", "mit_specificity": 68, "cfd_specificity": 84, "off_targets": "0-0-1-18-123", "flag": "Inefficient", "variant": False},
    {"guide": "GACACACTTACCAGATGGTG", "mit_specificity": 40, "cfd_specificity": 77, "off_targets": "0-0-3-23-229", "flag": "", "variant": False},
    {"guide": "TTTCGGGCAGAAGGCCATCC", "mit_specificity": 83, "cfd_specificity": 93, "off_targets": "0-0-1-7-72", "flag": "", "variant": False},
    {"guide": "CCCGAAAACACCTTCATCGT", "mit_specificity": 94, "cfd_specificity": 98, "off_targets": "0-0-0-3-34", "flag": "", "variant": True},
    {"guide": "AGAAGGGCTCACTCTGTTTG", "mit_specificity": 66, "cfd_specificity": 81, "off_targets": "0-0-1-27-203", "flag": "Inefficient", "variant": False},
    {"guide": "AGAGGAGAAGCTCAAGCTGG", "mit_specificity": 42, "cfd_specificity": 67, "off_targets": "0-1-5-90-405", "flag": "", "variant": False},
    {"guide": "CTTTGCCCGCAACTCCTATG", "mit_specificity": 84, "cfd_specificity": 93, "off_targets": "0-0-1-6-79", "flag": "", "variant": True},
    {"guide": "GGGCATTCATGTGGCTGTTG", "mit_specificity": 60, "cfd_specificity": 81, "off_targets": "0-0-4-26-186", "flag": "", "variant": False},
    {"guide": "GCAGGACTCGTGAATGTTCT", "mit_specificity": 78, "cfd_specificity": 87, "off_targets": "0-0-1-13-99", "flag": "Inefficient", "variant": False},
    {"guide": "ATCATCGTGGAGAAGCCCTT", "mit_specificity": 81, "cfd_specificity": 88, "off_targets": "0-0-4-8-125", "flag": "Inefficient", "variant": True},
    {"guide": "AGGAGATGTGGTTGGACAGC", "mit_specificity": 61, "cfd_specificity": 75, "off_targets": "0-0-3-25-287", "flag": "", "variant": True},
    {"guide": "GATCTGGTCCTCACGGAACA", "mit_specificity": 92, "cfd_specificity": 94, "off_targets": "0-0-0-7-61", "flag": "", "variant": False},
    {"guide": "TACCGCATCGACCACTACCT", "mit_specificity": 89, "cfd_specificity": 98, "off_targets": "0-0-0-3-24", "flag": "", "variant": False},
    {"guide": "ATGTTGTCCCGGTTCCAGAT", "mit_specificity": 86, "cfd_specificity": 94, "off_targets": "0-0-1-9-40", "flag": "", "variant": False},
    {"guide": "CTTGAAGGTGAGGATAACGC", "mit_specificity": 89, "cfd_specificity": 94, "off_targets": "0-0-2-5-110", "flag": "", "variant": False},
    {"guide": "CTCACCTTCAAGGAGCCCTT", "mit_specificity": 70, "cfd_specificity": 82, "off_targets": "0-0-3-19-156", "flag": "Inefficient", "variant": False},
    {"guide": "GACACAGCATCTGCAGTAGG", "mit_specificity": 48, "cfd_specificity": 75, "off_targets": "0-0-7-28-411", "flag": "", "variant": False},
    {"guide": "ATCACGGACGTCATCTGAGT", "mit_specificity": 94, "cfd_specificity": 96, "off_targets": "0-0-0-4-26", "flag": "", "variant": False},
    {"guide": "GAAATGCATCTCAGAGGTGC", "mit_specificity": 63, "cfd_specificity": 77, "off_targets": "0-0-2-29-214", "flag": "", "variant": False},
]

MIT_SPECIFICITY_MIN = 50


def close_mismatch_hits(off_targets):
    """True if CRISPOR found any hit at 1 or 2 mismatches (the closest, riskiest kind).
    Index 0 (the 0-mismatch slot) is deliberately ignored."""
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
    table = ours.merge(crispor, on="guide", validate="one_to_one")
    assert len(table) == len(OURS) == len(CRISPOR), "guide lists do not line up"
    table["comparison"] = table.apply(verdict, axis=1)
    table = table.sort_values("mit_specificity", ascending=False).reset_index(drop=True)

    agree = table[table["comparison"].str.startswith("Agreement")]
    clean_agree = table[table["comparison"] == "Agreement"]
    variant = table[table["variant"]]
    recommended = agree[~agree["variant"]]

    print(f"{len(agree)} of {len(table)} guides validated cleanly against CRISPOR.")
    print(f"{len(variant)} guide(s) overlap a known variant in CRISPOR's display.\n")
    print(table[["guide", "exon", "our_verdict", "mit_specificity", "cfd_specificity",
                 "off_targets", "flag", "variant", "comparison"]].to_string(index=False))

    lines = [
        "## Phase 3: validation against CRISPOR (G6PD)",
        "",
        "All 20 shortlisted G6PD guides were run through CRISPOR (crispor.tefor.net, "
        "hg38) across five small windows covering the gene's coding exons "
        "(guide positions in exons 1-8).",
        "",
        "| Guide | Exon | Our verdict | CRISPOR MIT spec. | CRISPOR CFD spec. | "
        "CRISPOR off-targets (0-1-2-3-4mm) | CRISPOR flag | Known variant? | Result |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for r in table.itertuples():
        lines.append(f"| `{r.guide}` | {r.exon} | {r.our_verdict} | {r.mit_specificity} | "
                     f"{r.cfd_specificity} | {r.off_targets} | {r.flag} | "
                     f"{'Yes' if r.variant else ''} | {r.comparison} |")
    lines += [
        "",
        f"**Only {len(agree)} of {len(table)} guides ({100*len(agree)//len(table)}%) "
        f"validated cleanly** -- the same systematic BLAST-sensitivity gap found on "
        f"MUTYH (3/20), KRAS (1/9) and PIP4K2C (3/5, a small hand-picked sample): our "
        f"`word_size=7` BLAST search misses off-targets whose 1-2 mismatches are spread "
        f"too thinly for a 7-bp perfect seed to exist. Of the {len(table) - len(agree)} "
        f"disagreements, every one has 2-mismatch hits and "
        f"{int(table['off_targets'].map(lambda s: int(s.split('-')[1]) > 0).sum())} "
        f"guide also has a 1-mismatch hit.",
        "",
        "**This gap is not explained by paralogs.** G6PD has only one paralog (H6PD), "
        "yet most guides still show close CRISPOR hits, and the genome-browser "
        "annotations CRISPOR displays (top three per guide) are mostly "
        "intergenic or intronic loci, and none of them is H6PD. The G6PD result (5/20) is "
        "similar to MUTYH's (3/20), not clearly better than the paralog-rich KRAS "
        "(1/9), so with four genes this project cannot claim that paralog count "
        "predicts CRISPOR-detected off-target risk; what it does show consistently is "
        "that a `word_size=7` BLAST screen alone under-reports close off-targets.",
        "",
        "**Variant caution:** " + f"{len(variant)} guides (" +
        ", ".join(f"`{g}`" for g in variant["guide"]) +
        ") overlap a known sequence variant in CRISPOR's display. G6PD is highly "
        "polymorphic, so these guides may fail to cut in carriers; check the variant "
        "position and allele frequency before using them.",
        "",
        f"**Validated on both methods:** " +
        ", ".join(f"`{g}`" for g in agree["guide"]) + ".",
        "",
        "**Final recommended G6PD guides** (validated on both methods, no overlapping "
        "known variant): " + ", ".join(f"`{g}`" for g in recommended["guide"]) + ".",
    ]
    with open("g6pd_phase3_summary.md", "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    table.to_csv("g6pd_phase3_comparison.csv", index=False)
    print("\nSaved g6pd_phase3_comparison.csv and g6pd_phase3_summary.md")
    