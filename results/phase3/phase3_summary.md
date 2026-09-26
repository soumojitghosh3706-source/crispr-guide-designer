## Phase 3: validation against CRISPOR

The top 5 "Clean" PIP4K2C guides were run through CRISPOR (crispor.tefor.net, hg38) to check this pipeline's results against an independent, published tool.

| Guide | Exon | Our score | CRISPOR MIT spec. | CRISPOR CFD spec. | CRISPOR off-targets (0-1-2-3-4mm) | Result |
|---|---|---|---|---|---|---|
| `CAATGCCAAATCGATCACGG` | 3 | 100.0 | 95 | 97 | 0-0-0-2-27 | Agreement |
| `CTTCGTTGTCCACACTGACT` | 5 | 100.0 | 79 | 92 | 0-0-0-15-77 | Agreement |
| `ATGCTTCTTCTTGGTCTTGG` | 1 | 90.0 | 68 | 71 | 0-0-0-19-239 | Agreement |
| `TCAGATCAATGAGCTCAGCC` | 2 | 97.3 | 67 | 84 | 0-0-2-18-168 | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `GAGCTGGCCTTAAAGTCATC` | 2 | 100.0 | 44 | 90 | 0-1-1-14-79 | Disagreement -- CRISPOR finds closer off-targets than our screen did |

**3 of 5 guides validated cleanly.** The other 2 were flagged by CRISPOR with off-target sites at 1-2 mismatches that this project's genomic BLAST screen missed -- most likely because CRISPOR uses a purpose-built, exhaustive genome index, while this project's screen relies on NCBI's general-purpose BLAST search, which is less sensitive at finding very-close near-matches. This is now recorded as a known limitation of the off-target screen (see Limitations).

**Final recommended guides** (validated on both methods): `CAATGCCAAATCGATCACGG`, `CTTCGTTGTCCACACTGACT`, `ATGCTTCTTCTTGGTCTTGG`.
