# CRISPR-GuideDesigner results: PIP4K2C (human)

- Target gene: PIP4K2C, 12,227 bp, 1266 nt coding sequence in 10 coding exons
- Candidate guides across the whole gene (NGG PAM, both strands): 1732
- Guides cutting a coding exon: 208
- Passed the exon-aware filters: 93
- Shortlisted and screened with genomic BLAST (self-matches sequence-verified, see fix_self_hits.py): 20
- Corrected verdict: **9 Clean**, 9 Low risk, 2 AT RISK

## Top guides (ranked by corrected BLAST verdict, then Step 3 score)

| Rank | Guide (5'-3') | PAM | Strand | Exon | CDS pos | Score | Verdict |
|---|---|---|---|---|---|---|---|
| 1 | `GAGCTGGCCTTAAAGTCATC` | TGG | - | 2 | 222 | 100.0 | Clean |
| 2 | `CAATGCCAAATCGATCACGG` | AGG | - | 3 | 335 | 100.0 | Clean |
| 3 | `CTTCGTTGTCCACACTGACT` | CGG | - | 5 | 572 | 100.0 | Clean |
| 4 | `TCAGATCAATGAGCTCAGCC` | AGG | + | 2 | 187 | 97.3 | Clean |
| 5 | `ATGCTTCTTCTTGGTCTTGG` | AGG | - | 1 | 79 | 90.0 | Clean |
| 6 | `GCACAGGAAGTATGACCTCA` | AGG | + | 5 | 655 | 85.3 | Clean |
| 7 | `AGGATGAGTCAGAGGTGGAT` | GGG | + | 8 | 927 | 81.4 | Clean |
| 8 | `CTCCTTGAACTTGAAATGAC` | TGG | - | 3 | 289 | 80.0 | Clean |
| 9 | `GGACTCTCACCAAGTAATCT` | TGG | - | 3 | 362 | 79.3 | Clean |
| 10 | `GAAGCATTTCGTGCAGCAGA` | AGG | + | 1 | 106 | 100.0 | Low risk |

## Limitations

- The score is a simple rule of thumb, not a validated efficiency model.
- Off-target self-matches are identified by comparing each hit's actual sequence against the gene (fix_self_hits.py), which is more reliable than matching on accession or coordinates alone, but still not a substitute for a full genome alignment.
- Guides have not been tested in the lab.
