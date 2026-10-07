# Sanger sequencing and indel analysis plan

Status: **designed, not yet validated** (no experiment has been run).
Covers the 11 guides with Primer-BLAST-checked genotyping primers (PIP4K2C_g1 is excluded: its site is repeat-rich and no specific primers were found).

## Purpose
T7E1 shows *that* a guide cut and was repaired imperfectly. Sanger sequencing of the same PCR product, analysed with a decomposition tool, shows *how much* editing occurred and *which* indels were made.

## 1. Which primer to sequence from
Sanger traces are unreliable for roughly the first 30-50 bases after the primer, so the cut site should sit well downstream of the sequencing primer. Common guidance is about 150 bp or more; here "good" means the cut is 150-500 bp past the primer's 3' end.

All 11 guides are "good" from **both** primers, so each amplicon can be read in either direction. The table shows the distance, and the primer to use if only one read is affordable (the one placing the cut closest to ~300 bp).

| Guide | Cut is this far past F primer | ...past R primer | Single-read pick |
|---|---|---|---|
| PIP4K2C_g2 | 242 bp | 380 bp | F |
| PIP4K2C_g3 | 397 bp | 281 bp | R |
| MUTYH_g1 | 391 bp | 272 bp | R |
| MUTYH_g2 | 375 bp | 262 bp | R |
| MUTYH_g3 | 253 bp | 410 bp | F |
| KRAS_g1 | 238 bp | 426 bp | F |
| G6PD_g1 | 218 bp | 429 bp | F |
| G6PD_g2 | 284 bp | 419 bp | F |
| G6PD_g3 | 438 bp | 220 bp | R |
| G6PD_g4 | 404 bp | 254 bp | R |
| G6PD_g5 | 188 bp | 406 bp | R |

Full numbers are in `sanger_read_plan.csv`. Reading both directions is better when possible: the two traces confirm each other.

## 2. Sample preparation
1. Extract genomic DNA from edited cells and from a matching control.
2. PCR with the designed primer pair using a high-fidelity polymerase.
3. Confirm a **single** band of the expected size on a gel (sizes in `t7e1_expected_bands.csv`). Several bands give an unreadable trace.
4. Clean up the product (column or enzymatic) and send it for Sanger sequencing with the F primer, the R primer, or both.

## 3. Controls (needed for the analysis)
- **Wild-type control:** the same cells with no editing, PCR'd and sequenced with the same primer. Decomposition tools use it as the reference trace.
- **No-guide control:** Cas9 plasmid without a guide inserted, to confirm that editing depends on the guide.
- **No-template PCR control:** water instead of DNA, to detect contamination.
- **Positive control (optional):** a guide already known to edit in your cell line.

## 4. Analysing the trace
Tools that estimate editing from a mixed trace: **TIDE**, **DECODR**, and **ICE**. Check which are currently online. Each takes:
- the control `.ab1` file,
- the edited-sample `.ab1` file,
- the 20-nt guide sequence (no PAM; sequences are in `final_guides.csv`).

Output: an overall editing percentage, an indel spectrum, and a fit score (R-squared). Treat results with a poor fit as unreliable and repeat the sequencing.

## 5. What to expect in the trace
- **Control:** clean single peaks all the way through.
- **Edited sample:** clean peaks up to the cut site, then overlapping peaks (a mixture of sequences) from the cut site onward.
- The cut lies about 3 bp upstream of the PAM. Check that the mixed region starts there; if it starts elsewhere, the primer or guide assignment is wrong.

## 6. Caveats
- Editing percentage depends on the cell line and delivery method, which are not yet chosen.
- G6PD is X-linked, so the number of X chromosomes in the chosen cell line affects how many alleles are edited.
- Decomposition of Sanger traces is less sensitive than amplicon sequencing and is unreliable below roughly 5% editing.
- All primers and guides are computational designs. Nothing here has been tested in cells.
