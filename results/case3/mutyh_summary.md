# CRISPR-GuideDesigner results: MUTYH (human)

- Target gene: MUTYH, 11,199 bp, 1650 nt coding sequence in 16 coding exons
- Candidate guides across the whole gene (NGG PAM, both strands): 1799
- Guides cutting a coding exon: 338
- Passed the exon-aware filters: 89
- Shortlisted and screened with genomic-only BLAST (self-matches sequence-verified from the first run): 20
- Verdict: **20 Clean**, 0 Low risk, 0 AT RISK

## Top guides (ranked by BLAST verdict, then Step 3 score)

| Rank | Guide (5'-3') | PAM | Strand | Exon | CDS pos | Score | Verdict |
|---|---|---|---|---|---|---|---|
| 1 | `AAGATGAGATGGACCTGGAC` | AGG | + | 4 | 369 | 100.0 | Clean |
| 2 | `TTGGTTGAAATCTCCTGGCC` | GGG | - | 10 | 823 | 100.0 | Clean |
| 3 | `CACAGGAGGTGAATCAACTC` | TGG | + | 7 | 516 | 96.0 | Clean |
| 4 | `ATCAGGTGGTAGAGGAGCTA` | GGG | + | 8 | 588 | 96.0 | Clean |
| 5 | `AGGAACAGCTCTTAGCCTCA` | GGG | + | 11 | 957 | 93.6 | Clean |
| 6 | `TGGGCTACTATTCTCGTGGC` | CGG | + | 7 | 546 | 90.0 | Clean |
| 7 | `GCCAAAGGCGATAGAGGCAA` | TGG | - | 8 | 670 | 90.0 | Clean |
| 8 | `TGGTGGATGGCAACGTAGCA` | CGG | + | 9 | 717 | 90.0 | Clean |
| 9 | `TGGGATCAGCACCAATGGCT` | CGG | - | 9 | 743 | 90.0 | Clean |
| 10 | `CCAGAGCTGCTGGGAAACAA` | GGG | - | 9 | 772 | 90.0 | Clean |

## Limitations

- The score is a simple rule of thumb, not a validated efficiency model.
- Off-target self-matches are identified by comparing each hit's actual sequence against the gene directly, built into the search from the start for this gene (unlike PIP4K2C, which needed a second patch script after the fact).
- Guides have not been tested in the lab.
