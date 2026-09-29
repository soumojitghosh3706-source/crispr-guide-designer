## Phase 3: validation against CRISPOR (KRAS)

All 9 shortlisted KRAS guides were run through CRISPOR (crispor.tefor.net, hg38) across three small windows covering the gene's three exon clusters.

| Guide | Exon | Our verdict | CRISPOR MIT spec. | CRISPOR CFD spec. | CRISPOR off-targets (0-1-2-3-4mm) | NRAS exon hit? | Result |
|---|---|---|---|---|---|---|---|
| `TCTCGACACAGCAGGTCAAG` | 2 | Clean | 80 | 90 | 0-0-1-12-89 |  | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `CAATGAGGGACCAGTACATG` | 2 | Clean | 79 | 91 | 0-0-1-13-116 | Yes | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `TCCCTTCTCAGGATTCCTAC` | 2 | Clean | 76 | 84 | 0-0-3-17-151 |  | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `CTGAATTAGCTGTATCGTCA` | 1 | Clean | 75 | 93 | 1-0-0-14-62 |  | Agreement |
| `GGACTCTGAAGATGTACCTA` | 3 | AT RISK | 67 | 81 | 0-1-3-28-189 | Yes | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `GTAGTTGGAGCTGGTGGCGT` | 1 | Clean | 54 | 75 | 1-0-2-21-277 |  | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `AACATCAGCAAAGACAAGAC` | 3 | Clean | 53 | 71 | 0-0-1-40-353 |  | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `AGTAGACACAAAACAGGCTC` | 3 | Clean | 53 | 81 | 0-0-7-103-213 |  | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `TGATGGAGAAACCTGTCTCT` | 2 | Clean | 38 | 79 | 1-0-6-37-211 |  | Disagreement -- CRISPOR finds closer off-targets than our screen did |

**Only 1 of 9 guides (11%) validated cleanly** -- consistent with the same systematic BLAST-sensitivity gap found on MUTYH (3/20) and PIP4K2C (3/5): our `word_size=7` BLAST search misses off-targets whose 1-2 mismatches are spread too thinly for a 7-bp perfect seed to exist.

**The paralog hypothesis is independently confirmed by CRISPOR itself:** 2 of the 9 guides (`CAATGAGGGACCAGTACATG`, and `GGACTCTGAAGATGTACCTA` -- the same guide our own BLAST screen already flagged AT RISK for a KRASP1 pseudogene hit) show an off-target site landing directly inside an **NRAS** exon in CRISPOR's own results. This means KRAS's off-target risk comes from *both* directions predicted: a near-identical pseudogene copy of itself (KRASP1, found by our own screen) AND its true functional paralog (NRAS, found only by CRISPOR's more sensitive search, missed by ours).

**Final recommended KRAS guide** (the only one validated on both methods): `CTGAATTAGCTGTATCGTCA`.
