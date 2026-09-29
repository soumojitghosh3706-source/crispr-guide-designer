"""
kras_exon_guides.py  --  Case study 4, Step 2 of CRISPR-GuideDesigner

Same logic as exon_guides.py / mutyh_exon_guides.py, retargeted to KRAS.

What it does:
  1. Reads kras_gene.fasta and kras_exons.csv (from kras_fetch.py).
  2. Finds every NGG guide on BOTH strands of the whole gene, reusing the
     Phase 1 PAM finder (pam_finder.py).
  3. Keeps the guides whose CUT SITE falls inside a coding exon.
  4. Annotates each kept guide with:
       coding_exon         which coding exon it cuts in
       cds_pos, cds_frac   where the cut is in the coding sequence
       dist_to_exon_edge   bases between the cut and the nearest exon edge
       shared_by_curated   True if that exon is in every curated isoform
       nmd_risk            True if a frameshift here may escape nonsense-
                           mediated decay (last exon, or last 50 nt of the
                           second-to-last exon)
  5. Self-check: the base just before every cut must be the same base in
     the gene and in the spliced coding sequence.
  6. Saves kras_guides_all.csv and kras_guides_coding.csv.

NOTE: KRAS's exon 5 is only 1 nt long (a known feature of its alternative
terminal-exon splicing). No guide can have its cut site land inside a
1-nt exon with any margin, so it will naturally attract zero or very few
guides once Step 3's edge-distance filter is applied -- this is expected,
not a bug.

Files needed in the same folder: pam_finder.py, kras_gene.fasta, kras_exons.csv
Requires: pandas, biopython
"""

import os

import pandas as pd
from Bio import SeqIO

from pam_finder import find_guides        # reuse the Phase 1 code

NMD_WINDOW = 50   # nt at the end of the second-to-last exon where a cut may escape NMD


def annotate_guides(guides, exons, nmd_window=NMD_WINDOW):
    n_exons = len(exons)
    cds_len = int(exons["cds_end"].max())
    ann = []
    for cut in guides["cut_pos"]:
        hit = exons[(exons["gene_start"] <= cut) & (cut + 1 <= exons["gene_end"])]
        if hit.empty:
            ann.append({"in_coding_exon": False})
            continue
        ex = hit.iloc[0]
        g_start, g_end = int(ex["gene_start"]), int(ex["gene_end"])
        exon_no = int(ex["coding_exon"])
        cds_pos = int(ex["cds_start"]) + (cut - g_start)
        ann.append({
            "in_coding_exon": True,
            "coding_exon": exon_no,
            "cds_pos": cds_pos,
            "cds_frac": round(cds_pos / cds_len, 3),
            "dist_to_exon_edge": min(cut - g_start + 1, g_end - cut),
            "shared_by_curated": bool(ex["shared_by_curated_isoforms"]),
            "nmd_risk": exon_no == n_exons
                        or (exon_no == n_exons - 1 and g_end - cut <= nmd_window),
        })
    out = pd.concat([guides.reset_index(drop=True), pd.DataFrame(ann)], axis=1)
    out["in_coding_exon"] = out["in_coding_exon"].fillna(False).astype(bool)
    for col in ("coding_exon", "cds_pos", "dist_to_exon_edge"):
        out[col] = out[col].astype("Int64")
    return out


if __name__ == "__main__":
    gene_seq = str(SeqIO.read("kras_gene.fasta", "fasta").seq).upper()
    exons = pd.read_csv("kras_exons.csv")
    if "shared_by_curated_isoforms" not in exons.columns:
        raise SystemExit("kras_exons.csv is missing 'shared_by_curated_isoforms'. "
                         "Run kras_fetch.py again first.")

    guides = find_guides(gene_seq)
    annotated = annotate_guides(guides, exons)
    coding = annotated[annotated["in_coding_exon"]].reset_index(drop=True)

    cds_len = int(exons["cds_end"].max())
    n_plus = int((guides["strand"] == "+").sum())
    n_minus = int((guides["strand"] == "-").sum())
    print(f"Gene length: {len(gene_seq):,} bp | coding sequence: {cds_len} nt "
          f"in {len(exons)} coding exons")
    print(f"Guides in the whole gene: {len(guides)}  (+ strand: {n_plus}, - strand: {n_minus})")
    print(f"Guides cutting inside a coding exon: {len(coding)}  "
          f"(rough expectation: about {cds_len // 8})")

    per_exon = coding.groupby("coding_exon").size()
    table = exons[["coding_exon", "length", "shared_by_curated_isoforms"]].copy()
    table["guides"] = table["coding_exon"].map(per_exon).fillna(0).astype(int)
    table["nmd_risk"] = table["coding_exon"].map(
        coding.groupby("coding_exon")["nmd_risk"].sum()).fillna(0).astype(int)
    print("\nGuides per coding exon:\n")
    print(table.to_string(index=False))

    print(f"\nCuts that may escape NMD: {int(coding['nmd_risk'].sum())}")
    print(f"Cuts in exons shared by all curated isoforms: "
          f"{int(coding['shared_by_curated'].sum())} of {len(coding)}")

    cds_seq = "".join(gene_seq[int(s) - 1:int(e)]
                      for s, e in zip(exons["gene_start"], exons["gene_end"]))
    agree = all(gene_seq[int(c) - 1] == cds_seq[int(p) - 1]
                for c, p in zip(coding["cut_pos"], coding["cds_pos"]))
    print(f"\nSelf-check: gene and coding-sequence positions agree for every guide: {agree}")

    annotated.to_csv("kras_guides_all.csv", index=False)
    coding.to_csv("kras_guides_coding.csv", index=False)
    print("\nSaved these files:")
    for name in ("kras_guides_all.csv", "kras_guides_coding.csv"):
        print("  " + os.path.abspath(name))
        