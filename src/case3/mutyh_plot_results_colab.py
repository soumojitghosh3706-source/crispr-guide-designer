"""
mutyh_plot_results_colab.py  --  Case study 3, Step 5 of CRISPR-GuideDesigner (Colab)

Same logic as exon_plot_results_colab.py for PIP4K2C, retargeted to MUTYH.

Draws MUTYH's exon/intron structure with every candidate and shortlisted
guide placed along it, ranks the shortlist by off-target verdict then
score, and writes a results summary.

HOW FILES MOVE
--------------
INPUT:  a file-picker dialog opens (Colab's upload widget) and you choose
        the files straight from your tablet's storage.
OUTPUT: once the plots and summary are built, a second dialog
        automatically downloads every output file back to your tablet.

Required input files:
  mutyh_gene.fasta
  mutyh_exons.csv
  mutyh_guides_all.csv
  mutyh_guides_coding.csv
  mutyh_guides_scored.csv
  mutyh_shortlist.csv
  mutyh_offtarget_verdict.csv

Output files (downloaded to your tablet automatically at the end):
  mutyh_gene_map.png
  mutyh_ranking.png
  mutyh_summary.md
"""

# ============================================================
# 1. INSTALL / IMPORT
# ============================================================
# Run once per Colab session if needed:
#   !pip -q install biopython pandas matplotlib

import os

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
from Bio import SeqIO

from google.colab import files

GENE_NAME = "MUTYH"
TOP_LABELS = 9

VERDICT_COLOR = {
    "Clean": "#2a9d8f", "Low risk": "#e9a23b",
    "AT RISK": "#e63946", "Not screened": "#999999",
}
VERDICT_ORDER = {"Clean": 0, "Low risk": 1, "AT RISK": 2, "Not screened": 9}

REQUIRED_FILES = [
    "mutyh_gene.fasta", "mutyh_exons.csv", "mutyh_guides_all.csv",
    "mutyh_guides_coding.csv", "mutyh_guides_scored.csv",
    "mutyh_shortlist.csv", "mutyh_offtarget_verdict.csv",
]

OUTPUT_FILES = ["mutyh_gene_map.png", "mutyh_ranking.png", "mutyh_summary.md"]


# ============================================================
# 2. GET INPUT FILES FROM YOUR TABLET
# ============================================================

def ensure_input_files(required=REQUIRED_FILES):
    missing = [f for f in required if not os.path.exists(f)]
    while missing:
        print("Please choose these files from your tablet's storage:")
        for f in missing:
            print("  -", f)
        uploaded = files.upload()
        for name, data in uploaded.items():
            with open(name, "wb") as fh:
                fh.write(data)
        missing = [f for f in required if not os.path.exists(f)]
        if missing:
            print("\nStill missing:", ", ".join(missing), "-- try again.\n")

    print("\nAll input files present:")
    for f in required:
        print("  \u2713", f)


# ============================================================
# 3. RANKING AND PLOTS
# ============================================================

def rank_shortlist(shortlist, verdict):
    merged = shortlist.merge(verdict[["guide", "verdict", "details"]], on="guide", how="left")
    merged["verdict"] = merged["verdict"].fillna("Not screened")
    merged["verdict_rank"] = merged["verdict"].map(VERDICT_ORDER).fillna(9).astype(int)
    merged = merged.sort_values(
        ["verdict_rank", "score"], ascending=[True, False], kind="stable"
    ).reset_index(drop=True)
    merged["rank"] = merged.index + 1
    merged["mid"] = (merged["start"] + merged["end"]) / 2
    return merged


def plot_gene_map(all_coding_guides, ranked, exons, gene_len, out):
    fig, ax = plt.subplots(figsize=(12, 4.6))
    y_of = {"+": 0.8, "-": -0.8}

    ax.plot([0, gene_len], [0, 0], color="#4c78a8", linewidth=1.5, zorder=1)

    for row in exons.itertuples():
        color = "#4c78a8" if row.shared_by_curated_isoforms else "#b0b0b0"
        ax.add_patch(Rectangle((row.gene_start, -0.15), row.length, 0.3,
                               color=color, zorder=2))
        ax.text((row.gene_start + row.gene_end) / 2, 0, str(row.coding_exon),
                color="white", ha="center", va="center", fontsize=6.5, zorder=3)

    for strand, y in y_of.items():
        sub = all_coding_guides[all_coding_guides["strand"] == strand]
        mid = (sub["start"] + sub["end"]) / 2
        lo, hi = sorted((y * 0.3, y * 0.55))
        ax.vlines(mid, lo, hi, color="#cfcfcf", linewidth=0.6, zorder=1)

    for strand, y in y_of.items():
        sub = ranked[ranked["strand"] == strand]
        if len(sub) == 0:
            continue
        colors = sub["verdict"].map(VERDICT_COLOR).fillna("#999999")
        ax.scatter(sub["mid"], [y] * len(sub), c=colors, s=70,
                  edgecolor="black", linewidth=0.6, zorder=4)

    for row in ranked.head(TOP_LABELS).itertuples():
        sign = 1 if row.strand == "+" else -1
        offset = 10 + 10 * (row.rank % 2)
        ax.annotate(str(row.rank), (row.mid, y_of[row.strand]),
                   xytext=(0, offset * sign), textcoords="offset points",
                   ha="center", va="center", fontsize=8, fontweight="bold")

    ax.set_xlim(-200, gene_len + 200)
    ax.set_ylim(-1.5, 1.5)
    ax.set_yticks(list(y_of.values()))
    ax.set_yticklabels(["+ strand", "- strand"])
    ax.set_xlabel("Position in gene (bp; 1 = the first base of the gene)")
    ax.set_title(f"CRISPR guides across {GENE_NAME} (human): "
                f"{len(exons)} coding exons, {len(all_coding_guides)} candidates, "
                f"{len(ranked)} shortlisted")
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)

    present = set(ranked["verdict"])
    handles = [
        Patch(color="#4c78a8", label="coding exon (shared by all curated isoforms)"),
        Patch(color="#b0b0b0", label="coding exon (not shared)"),
        Patch(color="#cfcfcf", label=f"all {len(all_coding_guides)} coding-exon candidates"),
    ] + [Patch(color=c, label=f"shortlisted: {v}") for v, c in VERDICT_COLOR.items() if v in present]
    ax.legend(handles=handles, loc="lower left", fontsize=7.5, frameon=False, ncol=2,
             bbox_to_anchor=(0.0, -0.32))

    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    return fig


def plot_ranking(ranked, out):
    fig, ax = plt.subplots(figsize=(10, 8))
    labels = [f"#{r.rank}  exon {r.coding_exon}, {r.start}-{r.end} ({r.strand})"
             for r in ranked.itertuples()]
    colors = [VERDICT_COLOR.get(v, "#999999") for v in ranked["verdict"]]
    ypos = list(range(len(ranked)))

    ax.barh(ypos, ranked["score"], color=colors)
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    for i, r in enumerate(ranked.itertuples()):
        ax.text(r.score + 0.8, i, f"{r.score:.1f}", va="center", fontsize=8)

    ax.set_xlim(0, 110)
    ax.set_xlabel("Step 3 score (GC + position + exon-edge distance)")
    ax.set_title(f"Ranked guides for {GENE_NAME}: BLAST verdict, then score")
    ax.legend(handles=[Patch(color=c, label=v) for v, c in VERDICT_COLOR.items()
                      if v in set(ranked["verdict"])],
             loc="lower right", fontsize=8)

    fig.tight_layout()
    fig.savefig(out, dpi=200)
    return fig


def write_summary(all_guides_count, all_coding, scored, ranked, exons,
                  gene_len, cds_len, out):
    counts = ranked["verdict"].value_counts()
    lines = [
        f"# CRISPR-GuideDesigner results: {GENE_NAME} (human)",
        "",
        f"- Target gene: {GENE_NAME}, {gene_len:,} bp, {cds_len} nt coding sequence "
        f"in {len(exons)} coding exons",
        f"- Candidate guides across the whole gene (NGG PAM, both strands): {all_guides_count}",
        f"- Guides cutting a coding exon: {len(all_coding)}",
        f"- Passed the exon-aware filters: {int(scored['pass_filters'].sum())}",
        f"- Shortlisted and screened with genomic-only BLAST "
        f"(self-matches sequence-verified from the first run): {len(ranked)}",
        f"- Verdict: **{int(counts.get('Clean', 0))} Clean**, "
        f"{int(counts.get('Low risk', 0))} Low risk, "
        f"{int(counts.get('AT RISK', 0))} AT RISK",
        "",
        "## Top guides (ranked by BLAST verdict, then Step 3 score)",
        "",
        "| Rank | Guide (5'-3') | PAM | Strand | Exon | CDS pos | Score | Verdict |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in ranked.head(10).itertuples():
        lines.append(f"| {r.rank} | `{r.guide}` | {r.pam} | {r.strand} | "
                     f"{r.coding_exon} | {r.cds_pos} | {r.score:.1f} | {r.verdict} |")
    lines += [
        "",
        "## Limitations",
        "",
        "- The score is a simple rule of thumb, not a validated efficiency model.",
        "- Off-target self-matches are identified by comparing each hit's actual "
        "sequence against the gene directly, built into the search from the "
        "start for this gene (unlike PIP4K2C, which needed a second patch "
        "script after the fact).",
        "- Guides have not been tested in the lab.",
    ]
    with open(out, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


# ============================================================
# 4. SAVE RESULTS BACK TO YOUR TABLET
# ============================================================

def download_outputs(paths=OUTPUT_FILES):
    print("\nDownloading results to your tablet...")
    for p in paths:
        if os.path.exists(p):
            files.download(p)
        else:
            print(f"  (skipped -- {p} was not created)")


# ============================================================
# 5. RUN EVERYTHING
# ============================================================

if __name__ == "__main__":
    ensure_input_files()

    gene_len = len(SeqIO.read("mutyh_gene.fasta", "fasta").seq)
    exons = pd.read_csv("mutyh_exons.csv")
    cds_len = int(exons["cds_end"].max())
    all_guides = pd.read_csv("mutyh_guides_all.csv")
    all_coding = pd.read_csv("mutyh_guides_coding.csv")
    scored = pd.read_csv("mutyh_guides_scored.csv")
    shortlist = pd.read_csv("mutyh_shortlist.csv")
    verdict = pd.read_csv("mutyh_offtarget_verdict.csv")

    ranked = rank_shortlist(shortlist, verdict)

    plot_gene_map(all_coding, ranked, exons, gene_len, "mutyh_gene_map.png")
    plot_ranking(ranked, "mutyh_ranking.png")
    write_summary(len(all_guides), all_coding, scored, ranked, exons,
                 gene_len, cds_len, "mutyh_summary.md")

    print("\nVerdict counts in the shortlist:",
         ranked["verdict"].value_counts().to_dict())
    print("\nTop 10:\n")
    cols = ["rank", "guide", "strand", "coding_exon", "score", "verdict"]
    print(ranked[cols].head(10).to_string(index=False))

    download_outputs()
    