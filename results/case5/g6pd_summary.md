# CRISPR-GuideDesigner results: G6PD (human)

- Target gene: G6PD, 16,180 bp, 1548 nt coding sequence in 12 coding exons
- Candidate guides across the whole gene (NGG PAM, both strands): 3057
- Guides cutting a coding exon: 301
- Passed the exon-aware filters: 138
- Shortlisted and screened with genomic-only BLAST (self-matches sequence-verified from the first run): 20
- Verdict: **19 Clean**, 1 Low risk, 0 AT RISK

## Top guides (ranked by BLAST verdict, then Step 3 score)

| Rank | Guide (5'-3') | PAM | Strand | Exon | CDS pos | Score | Verdict |
|---|---|---|---|---|---|---|---|
| 1 | `CCCGAAAACACCTTCATCGT` | GGG | + | 3 | 200 | 100.0 | Clean |
| 2 | `AGAAGGGCTCACTCTGTTTG` | CGG | - | 3 | 245 | 100.0 | Clean |
| 3 | `GCAGGACTCGTGAATGTTCT` | TGG | - | 4 | 457 | 100.0 | Clean |
| 4 | `ATCATCGTGGAGAAGCCCTT` | CGG | + | 5 | 515 | 100.0 | Clean |
| 5 | `ATGTTGTCCCGGTTCCAGAT` | GGG | - | 6 | 672 | 100.0 | Clean |
| 6 | `CTTGAAGGTGAGGATAACGC` | AGG | - | 6 | 697 | 100.0 | Clean |
| 7 | `ATCACGGACGTCATCTGAGT` | TGG | - | 7 | 841 | 96.5 | Clean |
| 8 | `GGTAGATCTTCTTCTTGGCC` | AGG | - | 2 | 131 | 94.7 | Clean |
| 9 | `GAAATGCATCTCAGAGGTGC` | AGG | + | 8 | 892 | 93.9 | Clean |
| 10 | `AGAGGAGAAGCTCAAGCTGG` | AGG | + | 4 | 292 | 90.0 | Clean |

## Limitations

- The score is a simple rule of thumb, not a validated efficiency model.
- Off-target self-matches are identified by comparing each hit's actual sequence against the gene directly, built into the search from the start for this gene.
- G6PD is heavily sequenced: for the 7 guides re-checked, all 50 BLAST hits returned were G6PD's own records, so weaker off-target sites may lie beyond the returned hit list. Clean means clean by BLAST only; verdicts are cross-checked against CRISPOR in Phase 3.
- Guides have not been tested in the lab.
