"""
plot_results.py  --  Final step of CRISPR-GuideDesigner (Phase 1)

What it does:
  1. Draws a map of all guide candidates along lacZ, with the shortlisted
     guides coloured by final score       -> guide_map.png
  2. Draws the ranked shortlist as a bar chart, coloured by whether the
     guide has any off-target hits        -> guide_ranking.png
  3. Writes a short results summary that you can paste into your README
                                          -> summary.md

Files needed in the same folder: lacZ.fasta, guides_all.csv,
guides_scored.csv, guides_final.csv
Requires: pandas, matplotlib, biopython
"""

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from Bio import SeqIO

GENE_NAME = "lacZ"
TOP_LABELS = 9          # how many top guides to number on the gene map
CLEAN_COLOR = "#2a9d8f"
HIT_COLOR = "#e9a23b"


def plot_gene_map(all_guides, final, gene_len, out="guide_map.png"):
    """Map of guide candidates along the gene (+ strand above, - strand below)."""
    fig, ax = plt.subplots(figsize=(11, 4.4))
    y_of = {"+": 0.8, "-": -0.8}

    # Region preferred by the Step 3 position score (5%-50% of the gene)
    ax.axvspan(0.05 * gene_len, 0.5 * gene_len, color="#f0f0f0", zorder=0)

    # The gene itself
    ax.add_patch(plt.Rectangle((0, -0.15), gene_len, 0.3, color="#4c78a8", zorder=2))
    ax.text(gene_len / 2, 0, f"{GENE_NAME} ({gene_len:,} bp)", color="white",
            ha="center", va="center", fontsize=10, zorder=3)

    # All candidate guides as faint ticks
    for strand, y in y_of.items():
        sub = all_guides[all_guides["strand"] == strand]
        mid = (sub["start"] + sub["end"]) / 2
        lo, hi = sorted((y * 0.3, y * 0.55))
        ax.vlines(mid, lo, hi, color="#b5b5b5", linewidth=0.6, zorder=1)

    # Shortlisted guides, coloured by final score
    vmin, vmax = final["final_score"].min(), final["final_score"].max()
    sc = None
    for strand, y in y_of.items():
        sub = final[final["strand"] == strand]
        if len(sub) == 0:
            continue
        sc = ax.scatter(sub["mid"], [y] * len(sub), c=sub["final_score"],
                        cmap="viridis", vmin=vmin, vmax=vmax, s=75,
                        edgecolor="black", linewidth=0.6, zorder=4)

    # Number the top guides (numbers match the ranking chart)
    for row in final.head(TOP_LABELS).itertuples():
        sign = 1 if row.strand == "+" else -1
        ax.annotate(str(row.rank), (row.mid, y_of[row.strand]),
                    xytext=(0, 10 * sign), textcoords="offset points",
                    ha="center", va="center", fontsize=8, fontweight="bold")

    ax.set_xlim(-30, gene_len + 30)
    ax.set_ylim(-1.5, 1.5)
    ax.set_yticks(list(y_of.values()))
    ax.set_yticklabels(["+ strand", "- strand"])
    ax.set_xlabel("Position in gene (bp; 1 = the A of ATG)")
    ax.set_title(f"CRISPR guides across {GENE_NAME} (E. coli K-12): "
                 f"{len(all_guides)} candidates, {len(final)} shortlisted")
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)

    if sc is not None:
        fig.colorbar(sc, ax=ax, label="final score", pad=0.01)
    handles = [
        Patch(color="#b5b5b5", label=f"all {len(all_guides)} candidates"),
        Patch(color="#f0f0f0", label="preferred cut region (5-50% of gene)"),
    ]
    ax.legend(handles=handles, loc="lower left", fontsize=8, frameon=False,
              bbox_to_anchor=(0.0, -0.02))

    fig.tight_layout()
    fig.savefig(out, dpi=200)
    return fig


def plot_ranking(final, out="guide_ranking.png"):
    """Ranked bar chart of the shortlist."""
    fig, ax = plt.subplots(figsize=(10, 8))
    labels = [f"#{r.rank}  {r.start}-{r.end} ({r.strand})" for r in final.itertuples()]
    colors = [CLEAN_COLOR if h == 0 else HIT_COLOR for h in final["total_hits"]]
    ypos = list(range(len(final)))

    ax.barh(ypos, final["final_score"], color=colors)
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    for i, r in enumerate(final.itertuples()):
        ax.text(r.final_score + 0.8, i, f"{r.final_score:.1f}", va="center", fontsize=8)

    ax.set_xlim(0, 110)
    ax.set_xlabel("Final score (Step 3 score x off-target specificity)")
    ax.set_title(f"Ranked guides for {GENE_NAME}")
    ax.legend(handles=[
        Patch(color=CLEAN_COLOR, label="no off-target hits (up to 4 mismatches)"),
        Patch(color=HIT_COLOR, label="has off-target hit(s)"),
    ], loc="lower right", fontsize=8)

    fig.tight_layout()
    fig.savefig(out, dpi=200)
    return fig


def write_summary(all_guides, scored, final, gene_len, out="summary.md"):
    """Short Markdown summary of the whole run, ready for a README."""
    n_clean = int((final["total_hits"] == 0).sum())
    lines = [
        f"# CRISPR-GuideDesigner results: {GENE_NAME} (E. coli K-12 MG1655)",
        "",
        f"- Target gene: {GENE_NAME}, {gene_len:,} bp",
        f"- Candidate guides (NGG PAM, both strands): {len(all_guides)}",
        f"- Passed filters (GC 40-60%, no TTTT/GGGG): {int(scored['pass_filters'].sum())}",
        f"- Shortlisted for off-target screening: {len(final)}",
        f"- Shortlisted guides with no off-target site up to 4 mismatches: "
        f"{n_clean} of {len(final)}",
        "",
        "## Top guides",
        "",
        "| Rank | Guide (5'-3') | PAM | Strand | Position | Final score | "
        "Off-target sites (0/1/2/3/4 mismatches) |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in final.head(10).itertuples():
        hits = "/".join(str(int(getattr(r, f"off_{k}mm"))) for k in range(5))
        lines.append(f"| {r.rank} | `{r.guide}` | {r.pam} | {r.strand} | "
                     f"{r.start}-{r.end} | {r.final_score:.1f} | {hits} |")
    lines += [
        "",
        "## Limitations",
        "",
        "- The score is a simple rule-of-thumb, not a validated efficiency model.",
        "- The off-target screen checks NGG sites only, up to 4 mismatches, "
        "with no insertions or deletions, in a single genome.",
        "- Guides have not been tested in the lab.",
    ]
    with open(out, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    gene_len = len(SeqIO.read("lacZ.fasta", "fasta").seq)
    all_guides = pd.read_csv("guides_all.csv")
    scored = pd.read_csv("guides_scored.csv")
    final = pd.read_csv("guides_final.csv")

    # Prepare helper columns
    final = final.sort_values("final_score", ascending=False).reset_index(drop=True)
    final["rank"] = final.index + 1
    final["mid"] = (final["start"] + final["end"]) / 2
    final["total_hits"] = final[[f"off_{k}mm" for k in range(5)]].sum(axis=1)

    plot_gene_map(all_guides, final, gene_len)
    plot_ranking(final)
    write_summary(all_guides, scored, final, gene_len)

    print("Saved guide_map.png, guide_ranking.png and summary.md")
    plt.show()
    
