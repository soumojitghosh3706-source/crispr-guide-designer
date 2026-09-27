## Phase 3: validation against CRISPOR (MUTYH)

All 20 shortlisted MUTYH guides -- our own pipeline had called all 20 "Clean" -- were run through CRISPOR (crispor.tefor.net, hg38) as a full cross-check.

| Guide | Exon | CRISPOR MIT spec. | CRISPOR CFD spec. | CRISPOR off-targets (0-1-2-3-4mm) | Result |
|---|---|---|---|---|---|
| `TGGGCTACTATTCTCGTGGC` | 7 | 90 | 97 | 0-0-0-5-63 | Agreement |
| `TGGTGGATGGCAACGTAGCA` | 9 | 88 | 92 | 0-0-0-8-89 | Agreement |
| `TGCAGGGTCTCTGCTGTACG` | 8 | 87 | 92 | 0-0-0-8-99 | Agreement |
| `TCACCACACTCCTCCACGTC` | 11 | 82 | 93 | 0-0-1-6-106 | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `TGGGATCAGCACCAATGGCT` | 9 | 80 | 88 | 0-0-2-14-129 | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `GCCAAAGGCGATAGAGGCAA` | 8 | 76 | 81 | 0-0-2-10-130 | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `AGGTGGCACTGTCCAGTGTT` | 12 | 75 | 90 | 0-0-4-12-138 | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `CACAGGAGGTGAATCAACTC` | 7 | 73 | 92 | 0-0-1-18-187 | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `CCTCTGCACCAGCAGAATTT` | 12 | 72 | 89 | 0-0-5-23-173 | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `TCAAGTATATGGGCTGGCCT` | 14 | 71 | 87 | 0-0-1-20-124 | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `AGTGGTCAACTTCCCCAGAA` | 12 | 68 | 83 | 0-0-2-17-160 | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `TTGCGCTGAAGCTGCTCTGA` | 13 | 68 | 84 | 0-0-1-16-116 | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `AGGAACAGCTCTTAGCCTCA` | 11 | 68 | 81 | 0-0-6-32-175 | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `ATCAGGTGGTAGAGGAGCTA` | 8 | 66 | 80 | 0-0-2-24-159 | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `AGATTTCAACCAAGCAGCCA` | 10 | 64 | 85 | 0-0-3-36-161 | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `AAGATGAGATGGACCTGGAC` | 4 | 64 | 83 | 0-1-4-24-208 | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `TTGGTTGAAATCTCCTGGCC` | 10 | 63 | 81 | 0-0-3-15-147 | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `AGGCTGTTCCAGAACACAGG` | 12 | 61 | 81 | 0-0-7-24-205 | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `GCTTGATGTGAGAGAAGGTG` | 14 | 61 | 76 | 0-0-1-25-242 | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `CCAGAGCTGCTGGGAAACAA` | 9 | 55 | 74 | 0-0-6-48-307 | Disagreement -- CRISPOR finds closer off-targets than our screen did |

**Only 3 of 20 guides (15%) validated cleanly.** The other 17 were flagged by CRISPOR with off-target sites at 1-2 mismatches that this project's own genomic BLAST screen missed entirely -- a systematic pattern, not an isolated case (PIP4K2C showed the same gap on 2 of 5 checked guides). The likely cause: our BLAST search uses `word_size=7`, so a guide with 2 mismatches spread evenly through its 20 nt often contains no single 7-bp perfectly-matching seed, and BLAST never registers it as a candidate at all -- while CRISPOR's purpose-built aligner has no such blind spot. This is recorded as the project's central methodological limitation (see Limitations).

**Final recommended MUTYH guides** (validated on both methods): `TGGGCTACTATTCTCGTGGC`, `TGGTGGATGGCAACGTAGCA`, `TGCAGGGTCTCTGCTGTACG`.
