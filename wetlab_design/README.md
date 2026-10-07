# CRISPR guide design: from computational picks to a wet-lab-ready design package

**Author:** Soumojit Ghosh
**Status: designed, not yet validated.** Everything here is computational. No guide, oligo or primer in this package has been tested in cells or in the lab.

## Overview
This package takes 12 validated SpCas9 guide RNAs across four human genes (*PIP4K2C*, *MUTYH*, *KRAS*, *G6PD*) and turns them into an experiment that a lab could run:
cloning oligos for a Cas9 plasmid, genotyping primers around every cut site, and a predicted readout for T7E1 and Sanger sequencing.

## What is in the package
```
design_package/
  README.md                      this file
  order_sheet.csv                42 oligos/primers to order (24 cloning oligos + 18 genotyping primers)
  data/
    final_guides.csv             the 12 guides that passed off-target screening
    final_guides_located.csv     strand, PAM and genomic cut position for each guide (GRCh38)
    oligo_order.csv              pX330 cloning oligos per guide
    primer_specificity.csv       Primer-BLAST verdict for each primer pair
    t7e1_expected_bands.csv      predicted T7E1 band sizes per guide
    sanger_read_plan.csv         cut-site distance from each primer, which read to use
  docs/sanger_indel_plan.md      Sanger sequencing and indel-analysis plan, controls
  figures/expected_t7e1_gel_schematic.png
  scripts/                       Python scripts used to generate the tables
```

## Workflow
| Step | What | Script | Output |
|---|---|---|---|
| 1 | Keep guides where my off-target verdict and CRISPOR agree | `make_final_guides.py` | `final_guides.csv` (12 guides) |
| 2 | Cloning oligos for pX330 (BbsI overhangs `CACC` / `AAAC`; add a 5' G when the guide does not start with G; check for internal BbsI sites and `TTTT`) | `make_oligos.py` | `oligo_order.csv` |
| 3 | Locate each guide on GRCh38, find the cut site, design a genotyping primer pair with an off-centre cut | `primer_design.py` | amplicon and primer tables |
| 4 | Check each pair for genome-wide specificity in NCBI Primer-BLAST (run manually in the web tool) | n/a | `primer_specificity.csv` |
| 5 | Predict T7E1 bands and plan the Sanger readout | `docs/` | `t7e1_expected_bands.csv`, `sanger_read_plan.csv` |

Primer design rules: product 600-800 bp (up to 1000 bp in a fallback mode); cut site off-centre so each T7E1 fragment is at least 200 bp and the two fragments differ by at least 100 bp; primers 20-24 nt (up to 27 in fallback), Tm about 60 C, GC 40-60%; low self-/cross-dimer scores; no primers inside repeats (see lessons below).

## The 12 guides
| Guide | Gene | Exon | Off-target verdict | MIT / CFD | Overlaps known variant | Uncut / cleaved bands (bp) |
|---|---|---|---|---|---|---|
| PIP4K2C_g1 | PIP4K2C | 3 | Clean | 95 / 97 | not checked | excluded (see below) |
| PIP4K2C_g2 | PIP4K2C | 5 | Clean | 79 / 92 | not checked | 665 / 401 + 264 |
| PIP4K2C_g3 | PIP4K2C | 1 | Clean | 68 / 71 | not checked | 722 / 418 + 304 |
| MUTYH_g1 | MUTYH | 7 | Clean | 90 / 97 | not checked | 704 / 412 + 292 |
| MUTYH_g2 | MUTYH | 9 | Clean | 88 / 92 | not checked | 679 / 396 + 283 |
| MUTYH_g3 | MUTYH | 8 | Clean | 87 / 92 | not checked | 704 / 430 + 274 |
| KRAS_g1 | KRAS | 1 | Clean | 75 / 93 | not checked | 708 / 447 + 261 |
| G6PD_g1 | G6PD | 7 | Clean | 94 / 96 | no | 690 / 450 + 240 |
| G6PD_g2 | G6PD | 3 | Clean | 94 / 98 | yes | 748 / 441 + 307 |
| G6PD_g3 | G6PD | 5 | Clean | 92 / 94 | no | 701 / 460 + 241 |
| G6PD_g4 | G6PD | 5 | Clean | 89 / 98 | no | 701 / 426 + 275 |
| G6PD_g5 | G6PD | 1 | Low risk | 85 / 92 | no | 637 / 428 + 209 |

- Off-target verdicts and MIT/CFD scores come from my off-target screen compared with CRISPOR. "Low risk" is weaker than "Clean".
- **G6PD_g2** overlaps a known variant; check the cell line's genotype before using it.
- **PIP4K2C_g1** has the best specificity score of its gene but sits in a repeat-rich region. Genotyping primers could not be made specific (every attempt matched many genomic sites), so it is kept for cloning (oligos designed) but excluded from the T7E1/Sanger plan.
- **MUTYH_g3** shares MUTYH_g1's primer pair, and **G6PD_g4** shares G6PD_g3's.

## Expected readout
![Expected T7E1 products](figures/expected_t7e1_gel_schematic.png)

An unedited control gives one band at the uncut size. An edited sample gives the uncut band plus two cleaved bands. T7E1 only cuts mismatched DNA, so it underestimates high editing rates. Sanger sequencing of the same PCR product, analysed with a tool such as TIDE, DECODR or ICE, gives the editing percentage and indel spectrum (see `docs/sanger_indel_plan.md`). The cut site sits 188-438 bp past both primers for every guide, so each amplicon can be read from either end.

## Primer specificity (summary)
Pairs checked in NCBI Primer-BLAST against human GRCh38 (RefSeq representative genomes): **Clean** (KRAS_g1, G6PD_g2, G6PD_g5) means only the intended product was reported. **Low risk** (the other pairs) means weak extra products were reported, each needing several mismatches per primer. Details are in `data/primer_specificity.csv`. For G6PD_g3/g4 only the first part of the product list was reviewed.

## Design iterations and lessons learned
- **Repeat-derived primers.** Primers copied from Alu repeat sequence matched hundreds of genomic sites in Primer-BLAST. My first designer had no repeat awareness. I added (a) an in-region uniqueness check, (b) an approximate Alu-consensus screen and (c) a check for near-copies elsewhere in the gene region. The Alu consensus in the script was written from memory, so treat it as approximate. Primer-BLAST remains the real specificity check.
- **Palindromic primers.** One reverse primer contained a 12-base palindrome, scored high for self-dimer formation and bound many weak genomic sites. I added a Primer3-style self- and cross-complementarity filter. Overall, 3 of the 10 pairs from the first design round failed Primer-BLAST (PIP4K2C_g1, G6PD_g3/g4, G6PD_g5); the two G6PD pairs were redesigned and passed, while PIP4K2C_g1 never produced a specific pair.
- **NCBI BLAST API timeouts.** `primer_blast.py` (batched BLAST through Biopython) never completed because NCBI's service timed out, so specificity was checked manually in the Primer-BLAST web tool. The script is kept for reference but is not part of the validated workflow.
- **Hand-entered rows.** The G6PD_g3/g4 and G6PD_g5 final primer pairs were entered from the redesign output and Primer-BLAST screenshots rather than generated in a single script run.

## Limitations
- No wet-lab validation of any kind.
- Cell line and delivery method are not chosen, so editing efficiency is unknown. G6PD is X-linked, so X-chromosome copy number in the chosen line affects how many alleles are edited.
- Primers were not checked against SNP databases; a variant under a primer's 3' end could block amplification of one allele.
- Primer-BLAST verdicts reflect one run per pair at default stringency.
- Sanger decomposition is unreliable below about 5% editing.

## How to reproduce
Requirements: Python 3 with Biopython and an internet connection (NCBI downloads).
1. Place the four `*_phase3_comparison.csv` files next to `make_final_guides.py` and run it.
2. Run `make_oligos.py` in the folder with `final_guides.csv`.
3. Set `Entrez.email` in `primer_design.py` and run it with `final_guides.csv`. To redesign only some guides, list them in `ONLY_GUIDES`.
4. Check each primer pair in NCBI Primer-BLAST (organism *Homo sapiens*, database RefSeq representative genomes).

## What a real experiment would involve
Order `order_sheet.csv`, anneal and phosphorylate each oligo pair, ligate into BbsI-digested pX330, confirm inserts by Sanger sequencing, transfect, extract genomic DNA, PCR with the genotyping primers, then run T7E1 and Sanger with TIDE/DECODR/ICE against an untreated control.

## Tools and references
Biopython; CRISPOR (guide scoring); NCBI Primer-BLAST; pX330 (Cong et al., Science 2013); TIDE (Brinkman et al., Nucleic Acids Res 2014).
