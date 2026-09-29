# CRISPR-GuideDesigner results: KRAS (human)

- Target gene: KRAS, 45,684 bp, 684 nt coding sequence in 6 coding exons
- Candidate guides across the whole gene (NGG PAM, both strands): 3875
- Guides cutting a coding exon: 51
- Passed the exon-aware filters: 21
- Shortlisted and screened with genomic-only BLAST (self-matches sequence-verified from the first run): 9
- Verdict: **8 Clean**, 0 Low risk, 1 AT RISK

## Top guides (ranked by BLAST verdict, then Step 3 score)

| Rank | Guide (5'-3') | PAM | Strand | Exon | CDS pos | Score | Verdict |
|---|---|---|---|---|---|---|---|
| 1 | `CAATGAGGGACCAGTACATG` | AGG | + | 2 | 213 | 100.0 | Clean |
| 2 | `TGATGGAGAAACCTGTCTCT` | TGG | + | 2 | 154 | 90.0 | Clean |
| 3 | `TCTCGACACAGCAGGTCAAG` | AGG | + | 2 | 181 | 90.0 | Clean |
| 4 | `TCCCTTCTCAGGATTCCTAC` | AGG | + | 2 | 117 | 88.0 | Clean |
| 5 | `AGTAGACACAAAACAGGCTC` | AGG | + | 3 | 388 | 84.6 | Clean |
| 6 | `GTAGTTGGAGCTGGTGGCGT` | AGG | + | 1 | 38 | 80.0 | Clean |
| 7 | `CTGAATTAGCTGTATCGTCA` | AGG | - | 1 | 58 | 80.0 | Clean |
| 8 | `AACATCAGCAAAGACAAGAC` | AGG | + | 3 | 445 | 54.6 | Clean |
| 9 | `GGACTCTGAAGATGTACCTA` | TGG | + | 3 | 328 | 90.0 | AT RISK |

## Limitations

- The score is a simple rule of thumb, not a validated efficiency model.
- Off-target self-matches are identified by comparing each hit's actual sequence against the gene directly, built into the search from the start for this gene (unlike PIP4K2C, which needed a second patch script after the fact).
- Guides have not been tested in the lab.
