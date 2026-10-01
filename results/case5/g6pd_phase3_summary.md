## Phase 3: validation against CRISPOR (G6PD)

All 20 shortlisted G6PD guides were run through CRISPOR (crispor.tefor.net, hg38) across five small windows covering the gene's coding exons (guide positions in exons 1-8).

| Guide | Exon | Our verdict | CRISPOR MIT spec. | CRISPOR CFD spec. | CRISPOR off-targets (0-1-2-3-4mm) | CRISPOR flag | Known variant? | Result |
|---|---|---|---|---|---|---|---|---|
| `ATCACGGACGTCATCTGAGT` | 7 | Clean | 94 | 96 | 0-0-0-4-26 |  |  | Agreement |
| `CCCGAAAACACCTTCATCGT` | 3 | Clean | 94 | 98 | 0-0-0-3-34 |  | Yes | Agreement |
| `GATCTGGTCCTCACGGAACA` | 5 | Clean | 92 | 94 | 0-0-0-7-61 |  |  | Agreement |
| `CTTGAAGGTGAGGATAACGC` | 6 | Clean | 89 | 94 | 0-0-2-5-110 |  |  | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `TACCGCATCGACCACTACCT` | 5 | Clean | 89 | 98 | 0-0-0-3-24 |  |  | Agreement |
| `ATGTTGTCCCGGTTCCAGAT` | 6 | Clean | 86 | 94 | 0-0-1-9-40 |  |  | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `GTGTGTATCCGACTGATGGA` | 1 | Low risk | 85 | 92 | 0-0-0-5-72 |  |  | Agreement |
| `CTTTGCCCGCAACTCCTATG` | 4 | Clean | 84 | 93 | 0-0-1-6-79 |  | Yes | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `TTTCGGGCAGAAGGCCATCC` | 3 | Clean | 83 | 93 | 0-0-1-7-72 |  |  | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `ATCATCGTGGAGAAGCCCTT` | 5 | Clean | 81 | 88 | 0-0-4-8-125 | Inefficient | Yes | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `GCAGGACTCGTGAATGTTCT` | 4 | Clean | 78 | 87 | 0-0-1-13-99 | Inefficient |  | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `CTCACCTTCAAGGAGCCCTT` | 6 | Clean | 70 | 82 | 0-0-3-19-156 | Inefficient |  | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `GGTAGATCTTCTTCTTGGCC` | 2 | Clean | 68 | 84 | 0-0-1-18-123 | Inefficient |  | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `AGAAGGGCTCACTCTGTTTG` | 3 | Clean | 66 | 81 | 0-0-1-27-203 | Inefficient |  | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `GAAATGCATCTCAGAGGTGC` | 8 | Clean | 63 | 77 | 0-0-2-29-214 |  |  | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `AGGAGATGTGGTTGGACAGC` | 5 | Clean | 61 | 75 | 0-0-3-25-287 |  | Yes | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `GGGCATTCATGTGGCTGTTG` | 4 | Clean | 60 | 81 | 0-0-4-26-186 |  |  | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `GACACAGCATCTGCAGTAGG` | 7 | Clean | 48 | 75 | 0-0-7-28-411 |  |  | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `AGAGGAGAAGCTCAAGCTGG` | 4 | Clean | 42 | 67 | 0-1-5-90-405 |  |  | Disagreement -- CRISPOR finds closer off-targets than our screen did |
| `GACACACTTACCAGATGGTG` | 2 | Clean | 40 | 77 | 0-0-3-23-229 |  |  | Disagreement -- CRISPOR finds closer off-targets than our screen did |

**Only 5 of 20 guides (25%) validated cleanly** -- the same systematic BLAST-sensitivity gap found on MUTYH (3/20), KRAS (1/9) and PIP4K2C (3/5, a small hand-picked sample): our `word_size=7` BLAST search misses off-targets whose 1-2 mismatches are spread too thinly for a 7-bp perfect seed to exist. Of the 15 disagreements, every one has 2-mismatch hits and 1 guide also has a 1-mismatch hit.

**This gap is not explained by paralogs.** G6PD has only one paralog (H6PD), yet most guides still show close CRISPOR hits, and the genome-browser annotations CRISPOR displays (top three per guide) are mostly intergenic or intronic loci, and none of them is H6PD. The G6PD result (5/20) is similar to MUTYH's (3/20), not clearly better than the paralog-rich KRAS (1/9), so with four genes this project cannot claim that paralog count predicts CRISPOR-detected off-target risk; what it does show consistently is that a `word_size=7` BLAST screen alone under-reports close off-targets.

**Variant caution:** 4 guides (`CCCGAAAACACCTTCATCGT`, `CTTTGCCCGCAACTCCTATG`, `ATCATCGTGGAGAAGCCCTT`, `AGGAGATGTGGTTGGACAGC`) overlap a known sequence variant in CRISPOR's display. G6PD is highly polymorphic, so these guides may fail to cut in carriers; check the variant position and allele frequency before using them.

**Validated on both methods:** `ATCACGGACGTCATCTGAGT`, `CCCGAAAACACCTTCATCGT`, `GATCTGGTCCTCACGGAACA`, `TACCGCATCGACCACTACCT`, `GTGTGTATCCGACTGATGGA`.

**Final recommended G6PD guides** (validated on both methods, no overlapping known variant): `ATCACGGACGTCATCTGAGT`, `GATCTGGTCCTCACGGAACA`, `TACCGCATCGACCACTACCT`, `GTGTGTATCCGACTGATGGA`.
