# CRISPR-GuideDesigner

> 🚧 **Work in progress** — Phase 1 is complete and validated below. Phases 2 (human gene, exon-aware targeting) and 3 (CRISPOR validation) are in progress; see [Roadmap](#roadmap).

A Python/Biopython pipeline that designs and ranks SpCas9 guide RNAs for a target gene and screens them for off-target sites across the whole genome.

**Demo target:** *lacZ* in *E. coli* K-12 MG1655 (NCBI `NC_000913.3`)

> Developed on an Android tablet using Pydroid 3.

## Why this project

Choosing a CRISPR guide is a genetics problem (where can Cas9 cut, and where might it cut by mistake?) that has to be solved computationally, because a single gene has hundreds of candidate sites and a genome has millions. This project encodes the biological rules of Cas9 targeting into code, step by step.

## Pipeline

| Step | Script | What it does | Output |
|---|---|---|---|
| 1 | `src/fetch.py` | Downloads the *E. coli* genome from NCBI and extracts *lacZ* | `ecoli_genome.gb`, `lacZ.fasta` |
| 2 | `src/pam_finder.py` | Finds every NGG PAM on both strands and extracts the 20-nt guide in front of each | `guides_all.csv` |
| 3 | `src/scoring.py` | Filters (GC 40-60%, no TTTT/GGGG), scores, and shortlists 20 spaced-out guides | `guides_scored.csv`, `guides_shortlist.csv` |
| 4 | `src/offtarget.py` | Compares each guide with every NGG site in the genome (up to 4 mismatches, seed region weighted double) | `guides_final.csv`, `offtarget_hits.csv` |
| 5 | `src/plot_results.py` | Plots the guide map and ranking, and writes a summary | `guide_map.png`, `guide_ranking.png`, `summary.md` |

## Results (*lacZ*, 3,075 bp)

- **415** candidate guides (NGG PAM, both strands)
- **242** pass the GC and motif filters
- **20** shortlisted for off-target screening
- **11 of 20** shortlisted guides have no off-target site up to 4 mismatches; the rest have only 4-mismatch hits, and none has a hit with 1-3 mismatches

![Guide map](results/guide_map.png)

![Guide ranking](results/guide_ranking.png)

### Top guides

| Rank | Guide (5'-3') | PAM | Strand | Position | Final score | Off-target sites (0/1/2/3/4 mismatches) |
|---|---|---|---|---|---|---|
| 1 | `CGTAATGGGATAGGTCACGT` | TGG | - | 308-327 | 100.0 | 0/0/0/0/0 |
| 2 | `TATGCAGCAACGAGACGTCA` | CGG | - | 633-652 | 100.0 | 0/0/0/0/0 |
| 3 | `GCAAATAATATCGGTGGCCG` | TGG | - | 1484-1503 | 100.0 | 0/0/0/0/0 |
| 4 | `ATTCATTGGCACCATGCCGT` | GGG | - | 1256-1275 | 100.0 | 0/0/0/0/0 |
| 5 | `ATTATCCGAACCATCCGCTG` | TGG | + | 1175-1194 | 100.0 | 0/0/0/0/0 |
| 6 | `GGATTCTCCGTGGGAACAAA` | CGG | - | 340-359 | 100.0 | 0/0/0/0/0 |
| 7 | `ACCACCGCACGATAGAGATT` | CGG | - | 922-941 | 100.0 | 0/0/0/0/0 |
| 8 | `GCGCTCAGGTCAAATTCAGA` | CGG | - | 508-527 | 100.0 | 0/0/0/0/0 |
| 9 | `GATGAAAGCTGGCTACAGGA` | AGG | + | 391-410 | 100.0 | 0/0/0/0/0 |
| 10 | `CGTATTCGCAAAGGATCAGC` | GGG | - | 1599-1618 | 98.3 | 0/0/0/0/0 |

Positions are inside the gene (1 = the A of the ATG start codon).

## How genetics and bioinformatics are integrated

The genetics decides the rules; the bioinformatics applies them at scale.

| Genetics | Bioinformatics |
|---|---|
| Cas9 needs an NGG PAM next to the target | String scan of both strands using reverse complement |
| Guide efficiency depends on GC content | GC filter and scoring with `Bio.SeqUtils` |
| Early cuts in a gene are more likely to knock it out | Position-based scoring along the gene |
| Mismatches near the PAM (seed region) block cutting most strongly | Seed-weighted mismatch counting across the genome with numpy |
| Off-target cuts cause unwanted mutations | Genome-wide comparison of every guide against all NGG sites |

## Run it

```bash
pip install -r requirements.txt
# open src/fetch.py and set Entrez.email to your own email address
python src/fetch.py
python src/pam_finder.py
python src/scoring.py
python src/offtarget.py
python src/plot_results.py
```

Run every script from the same folder (the repository root). Inputs and outputs are read from and written to the current folder. The genome file (`ecoli_genome.gb`, about 12 MB) is git-ignored.

## Limitations

- The score is a simple rule of thumb, not a validated efficiency model.
- The off-target screen checks NGG sites only, up to 4 mismatches, with no insertions or deletions, in a single genome.
- Guides have not been tested in the lab.

## Roadmap

- [x] **Phase 1:** core pipeline built and validated on *E. coli lacZ* (this repo)
- [ ] **Phase 2:** exon-aware targeting on a human gene (PIP4K2C) — restrict guides to coding exons, check off-targets with NCBI BLAST — *in progress*
- [ ] **Phase 3:** compare the top guides against the published CRISPOR tool and report where they agree and differ

## Author

Soumojit Ghosh, B.Sc. Biotechnology (Honours with Research), St. Xavier's College, Burdwan.

## License

MIT

