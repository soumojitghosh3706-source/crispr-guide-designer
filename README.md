# CRISPR-GuideDesigner

A Python/Biopython pipeline that designs and ranks SpCas9 guide RNAs for a target gene, screens them for off-target sites, and validates the top picks against a published tool (CRISPOR).

**Case studies:**
1. *lacZ* in *E. coli* K-12 MG1655 (`NC_000913.3`) — proves the core pipeline with an exhaustive, genome-wide off-target search
2. **PIP4K2C** (human, chr12) — extends the pipeline to a eukaryotic gene with exons/introns, adds exon-aware filtering, a genomic BLAST off-target screen, and validation against CRISPOR
3. **MUTYH** (human, chr1) — a clinically established DNA-repair gene, contrasted with PIP4K2C's understudied kinase; full-shortlist CRISPOR validation revealed the project's central methodological limitation (see Phase 3 below)

Planned next: **KRAS** (autosomal, paralog-rich) and **G6PD** (X-linked, paralog-poor), to test whether paralog count predicts off-target rate.

> Developed on an Android tablet using Pydroid 3 and Google Colab.

## Why this project

Choosing a CRISPR guide is a genetics problem (where can Cas9 cut, and where might it cut by mistake?) that has to be solved computationally, because a single gene has hundreds of candidate sites and a genome has millions. This project encodes the biological rules of Cas9 targeting into code, step by step, and checks its own conclusions against real gene structure and an independent published tool rather than assuming its heuristics are correct.

---

## Case study 1: *E. coli lacZ*

### Pipeline

| Step | Script | What it does | Output |
|---|---|---|---|
| 1 | `src/fetch.py` | Downloads the *E. coli* genome from NCBI and extracts *lacZ* | `ecoli_genome.gb`, `lacZ.fasta` |
| 2 | `src/pam_finder.py` | Finds every NGG PAM on both strands and extracts the 20-nt guide in front of each | `guides_all.csv` |
| 3 | `src/scoring.py` | Filters (GC 40-60%, no TTTT/GGGG), scores, and shortlists 20 spaced-out guides | `guides_scored.csv`, `guides_shortlist.csv` |
| 4 | `src/offtarget.py` | Compares each guide with every NGG site in the genome (up to 4 mismatches, seed region weighted double) | `guides_final.csv`, `offtarget_hits.csv` |
| 5 | `src/plot_results.py` | Plots the guide map and ranking, and writes a summary | `guide_map.png`, `guide_ranking.png`, `summary.md` |

### Results (*lacZ*, 3,075 bp)

- **415** candidate guides (NGG PAM, both strands) -> **242** pass GC/motif filters -> **20** shortlisted
- **11 of 20** shortlisted guides have no off-target site up to 4 mismatches, checked against the **entire** *E. coli* genome

![Guide map](results/guide_map.png)

![Guide ranking](results/guide_ranking.png)

| Rank | Guide (5'-3') | PAM | Strand | Position | Score | Off-targets (0/1/2/3/4mm) |
|---|---|---|---|---|---|---|
| 1 | `CGTAATGGGATAGGTCACGT` | TGG | - | 308-327 | 100.0 | 0/0/0/0/0 |
| 2 | `TATGCAGCAACGAGACGTCA` | CGG | - | 633-652 | 100.0 | 0/0/0/0/0 |
| 3 | `GCAAATAATATCGGTGGCCG` | TGG | - | 1484-1503 | 100.0 | 0/0/0/0/0 |

Positions are inside the gene (1 = the A of the ATG start codon). Full table in `results/summary.md`.

---

## Case study 2: human PIP4K2C

PIP4K2C is a poorly-characterized lipid kinase -- chosen deliberately as an understudied gene, in contrast to a heavily-studied one like TP53. Unlike *lacZ*, it has introns, so a candidate guide's position relative to exon structure now matters.

### Pipeline

| Step | Script | What it does | Output |
|---|---|---|---|
| 1 | `src/phase2/phase2_fetch.py` | Downloads the PIP4K2C region from NCBI and builds a coding-exon table across all isoforms | `pip4k2c_gene.fasta`, `pip4k2c_exons.csv` |
| 2 | `src/phase2/exon_guides.py` | Finds every guide across the whole gene, keeps guides whose cut lands in a coding exon | `pip4k2c_guides_all.csv`, `pip4k2c_guides_coding.csv` |
| 3 | `src/phase2/exon_scoring.py` | Filters (GC, motifs, exon shared by all curated isoforms, not an NMD-escape zone, >=5 bp from a splice edge), scores, shortlists 20 | `pip4k2c_guides_scored.csv`, `pip4k2c_shortlist.csv` |
| 4 | `src/phase2/blast_offtargets_v_2.py` | Screens each guide with **genomic-only** NCBI BLAST (transcript/mRNA records excluded), collapses duplicate records into unique loci | `pip4k2c_blast_hits_genomic.csv`, `pip4k2c_unique_loci.csv` |
| 4b | `src/phase2/fix_self_hits.py` | Re-checks every hit's actual sequence against the gene itself, to catch self-matches under an unrelated accession that coordinate/text matching missed | `pip4k2c_offtarget_verdict_fixed.csv` |
| 5 | `src/phase2/exon_plot_results_colab.py` | Plots the exon/intron map and ranking, colored by verdict, and writes a summary | `pip4k2c_gene_map.png`, `pip4k2c_ranking.png`, `pip4k2c_summary.md` |

### Results (PIP4K2C, 12,227 bp, 1,266 nt CDS across 10 coding exons)

- **1,732** candidate guides across the whole gene -> **208** cut a coding exon -> **93** pass the exon-aware filters -> **20** shortlisted
- Off-target verdict after correction: **9 Clean, 9 Low risk, 2 AT RISK** (of 20)
  - AT RISK guides hit **PIP5K1C** (a related kinase) and **KCNQ1** (an unrelated cardiac ion-channel gene) with an adjacent PAM present
  - Exon 4 excluded (not shared by all curated isoforms); exons 9-10 largely excluded (nonsense-mediated-decay escape risk)

![PIP4K2C gene map](results/phase2/pip4k2c_gene_map.png)

![PIP4K2C ranking](results/phase2/pip4k2c_ranking.png)

| Rank | Guide (5'-3') | PAM | Strand | Exon | Score | Verdict |
|---|---|---|---|---|---|---|
| 1 | `GAGCTGGCCTTAAAGTCATC` | TGG | - | 2 | 100.0 | Clean |
| 2 | `CAATGCCAAATCGATCACGG` | AGG | - | 3 | 100.0 | Clean |
| 3 | `CTTCGTTGTCCACACTGACT` | CGG | - | 5 | 100.0 | Clean |
| 4 | `TCAGATCAATGAGCTCAGCC` | AGG | + | 2 | 97.3 | Clean |
| 5 | `ATGCTTCTTCTTGGTCTTGG` | AGG | - | 1 | 90.0 | Clean |

Full 20-guide table in `results/phase2/pip4k2c_summary.md`.

---

## Phase 3: validation against CRISPOR

The 5 top "Clean" PIP4K2C guides above were run through [CRISPOR](https://crispor.tefor.net) (hg38) to check this pipeline's results against an independent, published tool. CRISPOR has no batch API suited to a tablet workflow, so its results were read from its results table by hand (`src/phase3/phase3_compare_crispor.py` documents the readings and builds the comparison).

| Guide | Exon | Our score | CRISPOR MIT spec. | CRISPOR CFD spec. | CRISPOR off-targets (0-1-2-3-4mm) | Result |
|---|---|---|---|---|---|---|
| `CAATGCCAAATCGATCACGG` | 3 | 100.0 | 95 | 97 | 0-0-0-2-27 | Agreement |
| `CTTCGTTGTCCACACTGACT` | 5 | 100.0 | 79 | 92 | 0-0-0-15-77 | Agreement |
| `ATGCTTCTTCTTGGTCTTGG` | 1 | 90.0 | 68 | 71 | 0-0-0-19-239 | Agreement |
| `TCAGATCAATGAGCTCAGCC` | 2 | 97.3 | 67 | 84 | 0-0-2-18-168 | Disagreement (CRISPOR also flags "Inefficient") |
| `GAGCTGGCCTTAAAGTCATC` | 2 | 100.0 | 44 | 90 | 0-1-1-14-79 | Disagreement |

**3 of 5 guides validated cleanly.** The other 2 were flagged by CRISPOR with off-target sites at 1-2 mismatches that this project's own BLAST screen missed -- most likely because CRISPOR uses a purpose-built, exhaustive genome index, while this project relies on NCBI's general-purpose BLAST search, which is less sensitive at finding very-close near-matches. This is recorded as a known limitation below rather than papered over.

**Final recommended PIP4K2C guides** (validated on both methods): `CAATGCCAAATCGATCACGG`, `CTTCGTTGTCCACACTGACT`, `ATGCTTCTTCTTGGTCTTGG`.

---

## Case study 3: human MUTYH

MUTYH (mutY DNA glycosylase) is a base-excision-repair enzyme; loss-of-function mutations cause MUTYH-associated polyposis, a hereditary colorectal cancer syndrome. Chosen deliberately to contrast with PIP4K2C: a clinically established gene rather than an understudied one, with 16 coding exons rather than 10, and (as it turned out) no close paralogs, unlike PIP4K2C's PIP4K2A/PIP4K2B family.

### Pipeline

Same five-script structure as PIP4K2C (`src/phase2/` scripts, retargeted), plus two MUTYH-specific additions: `mutyh_blast_offtargets.py` builds the genomic-only, sequence-verified self-hit check directly into the first run (no separate patch script needed, unlike PIP4K2C), and `mutyh_retry_flagged.py` is a small helper for forcing a re-check of specific guides.

| Step | Script | What it does | Output |
|---|---|---|---|
| 1 | `src/case3/mutyh_fetch.py` | Downloads the MUTYH region from NCBI and builds the coding-exon table | `mutyh_gene.fasta`, `mutyh_exons.csv` |
| 2 | `src/case3/mutyh_exon_guides.py` | Finds every guide across the whole gene, keeps guides whose cut lands in a coding exon | `mutyh_guides_all.csv`, `mutyh_guides_coding.csv` |
| 3 | `src/case3/mutyh_exon_scoring.py` | Filters and scores, shortlists 20 | `mutyh_guides_scored.csv`, `mutyh_shortlist.csv` |
| 4 | `src/case3/mutyh_blast_offtargets.py` | Genomic-only BLAST off-target screen with sequence-verified self-hit exclusion, built in from the start | `mutyh_unique_loci.csv`, `mutyh_offtarget_verdict.csv` |
| 5 | `src/case3/mutyh_plot_results_colab.py` | Plots the exon/intron map and ranking, writes a summary | `mutyh_gene_map.png`, `mutyh_ranking.png`, `mutyh_summary.md` |

### Results (MUTYH, 11,199 bp, 1,650 nt CDS across 16 coding exons)

- **1,799** candidate guides across the whole gene → **338** cut a coding exon → **89** pass the exon-aware filters → **20** shortlisted
- This project's own off-target screen found **0 off-target candidates for all 20** shortlisted guides -- see Phase 3 below for why this number alone is misleading
- Exons 1-6 excluded from the shortlist (not shared by all curated isoforms, matching MUTYH's known alternative-first-exon structure); exons 15-16 mostly excluded (NMD-escape risk)

![MUTYH gene map](results/case3/mutyh_gene_map.png)

![MUTYH ranking](results/case3/mutyh_ranking.png)

Full 20-guide table in `results/case3/mutyh_summary.md`.

---

## Phase 3: validation against CRISPOR (both genes)

For PIP4K2C, the top 5 guides were checked (see above). For MUTYH, since this project's own screen found zero off-target candidates across the board, **all 20 shortlisted guides** were checked instead, to properly stress-test that result.

**MUTYH: only 3 of 20 guides (15%) validated cleanly** -- a much larger gap than PIP4K2C's 3-of-5 (60%).

| Guide | Exon | Our score | CRISPOR MIT spec. | CRISPOR CFD spec. | CRISPOR off-targets (0-1-2-3-4mm) | Result |
|---|---|---|---|---|---|---|
| `TGGGCTACTATTCTCGTGGC` | 7 | 90.0 | 90 | 97 | 0-0-0-5-63 | ✅ Agreement |
| `TGGTGGATGGCAACGTAGCA` | 9 | 90.0 | 88 | 92 | 0-0-0-8-89 | ✅ Agreement |
| `TGCAGGGTCTCTGCTGTACG` | 8 | 80.0 | 87 | 92 | 0-0-0-8-99 | ✅ Agreement |
| ...17 more guides | | | | | | ⚠️ Disagreement (2 also flagged "Inefficient") |

Full 20-row table in `results/case3/mutyh_phase3_comparison.csv` and `results/case3/mutyh_phase3_summary.md`.

**Final recommended MUTYH guides** (validated on both methods): `TGGGCTACTATTCTCGTGGC`, `TGGTGGATGGCAACGTAGCA`, `TGCAGGGTCTCTGCTGTACG`.

### The central methodological finding

Confirmed across **two independent genes** now, not a one-off: this project's BLAST-based off-target screen (`word_size=7`) systematically misses off-target sites where 2 mismatches are spread across the guide's 20 bp, because such a hit often contains no 7-bp perfectly-matching stretch for BLAST to seed on. CRISPOR's purpose-built aligner has no such blind spot. This means **"0 off-target candidates" from this project's own screen should be read as "0 candidates BLAST's seed-based search could find," not as a guarantee of true specificity** -- a materially different, more honest claim, and the single most important limitation of this project.

---

## How genetics and bioinformatics are integrated

| Genetics | Bioinformatics |
|---|---|
| Cas9 needs an NGG PAM next to the target | String scan of both strands using reverse complement |
| Guide efficiency depends on GC content | GC filter and scoring with `Bio.SeqUtils` |
| Coding exons, especially early ones, give the best knockouts | GenBank feature parsing, exon-based filtering and ranking |
| Mismatches near the PAM (seed region) block cutting most strongly | Seed-weighted mismatch counting (Phase 1) |
| A frameshift near the last exon may escape nonsense-mediated decay | Cut-position rule excluding NMD-escape zones (Phase 2/3) |
| Cutting near a splice junction risks abnormal splicing | Minimum distance-from-exon-edge filter (Phase 2/3) |
| Off-target cuts cause unwanted mutations | Genome-wide comparison against all NGG sites (Phase 1), genomic-only BLAST with sequence-verified self-hit exclusion (Phase 2/3) |
| Genes with close paralogs are inherently harder to target uniquely | Cross-gene comparison of off-target rates (PIP4K2C/paralog-rich vs. MUTYH/paralog-poor) -- planned to extend with KRAS and G6PD |
| Published tools are the field's benchmark | Comparison against CRISPOR's own scores and off-target calls (Phase 3, both genes) |

## Run it

```bash
pip install -r requirements.txt
# open each fetch/blast script and set Entrez.email to your own email address

# Case study 1: E. coli lacZ
python src/fetch.py
python src/pam_finder.py
python src/scoring.py
python src/offtarget.py
python src/plot_results.py

# Case study 2: human PIP4K2C
python src/phase2/phase2_fetch.py
python src/phase2/exon_guides.py
python src/phase2/exon_scoring.py
python src/phase2/blast_offtargets_v_2.py
python src/phase2/fix_self_hits.py
python src/phase2/exon_plot_results_colab.py   # Colab: prompts for input files, downloads outputs
python src/phase3/phase3_compare_crispor.py

# Case study 3: human MUTYH
python src/case3/mutyh_fetch.py
python src/case3/mutyh_exon_guides.py
python src/case3/mutyh_exon_scoring.py
python src/case3/mutyh_blast_offtargets.py
python src/case3/mutyh_retry_flagged.py        # only if a guide needs re-verification
python src/case3/mutyh_plot_results_colab.py   # Colab: prompts for input files, downloads outputs
python src/case3/mutyh_compare_crispor.py
```

Run scripts from the repository root; inputs/outputs are read from and written to the current folder. Large genome files (`ecoli_genome.gb`, `pip4k2c_region.gb`, `mutyh_region.gb`) are git-ignored, except a saved copy under `data/` kept for reproducibility.

## Limitations

- **This project's BLAST-based off-target screen (`word_size=7`) systematically misses off-targets where 2 mismatches are spread across the guide, confirmed across both PIP4K2C and MUTYH via CRISPOR validation.** A "0 off-target candidates" result from this project's own screen means no 7-bp-seedable hit was found, not that the guide is truly free of off-targets. This is the project's central methodological finding, not a per-gene quirk.
- The scores are a simple rule of thumb, not a validated efficiency model -- confirmed directly in Phase 3, where CRISPOR flagged several top-scoring guides "Inefficient."
- Phase 1's off-target screen checks NGG sites only, up to 4 mismatches, no indels, in a single genome (this search is exhaustive, unlike Phases 2/3's BLAST-based approach, so the word_size limitation above does not apply to Phase 1).
- Self-match detection (Phases 2/3) compares hit sequence against the gene directly, which is more reliable than accession/coordinate matching alone, but still not a full genome alignment.
- Guides have not been tested in the lab.

## Roadmap

- [x] **Case study 1:** core pipeline built and validated on *E. coli lacZ*
- [x] **Case study 2:** exon-aware targeting on human PIP4K2C, genomic BLAST off-target screen, validated against CRISPOR
- [x] **Case study 3:** exon-aware targeting on human MUTYH, validated against CRISPOR -- revealed the project's central BLAST-sensitivity limitation
- [ ] **Case study 4:** KRAS (autosomal, chr12, paralog-rich cancer driver gene) -- planned, to test whether paralog count predicts off-target rate
- [ ] **Case study 5:** G6PD (X-linked, Xq28, paralog-poor metabolic enzyme) -- planned, for a 3-autosomal + 1-X-linked comparison

## Author

Soumojit Ghosh, B.Sc. Biotechnology (Honours with Research), St. Xavier's College, Burdwan.

## License

MIT

