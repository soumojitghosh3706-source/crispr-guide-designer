"""
offtarget.py  --  Step 4 of CRISPR-GuideDesigner (Phase 1)

What it does:
  1. Loads the E. coli genome (ecoli_genome.gb, downloaded by fetch.py).
  2. Builds a table of EVERY possible NGG site in the genome (both strands),
     storing the 20 nt in front of each PAM.
  3. Compares each shortlisted guide against all sites at once (numpy).
  4. Reports sites with up to 4 mismatches. Mismatches in the 12 nt nearest
     the PAM (the "seed") count double, because they block cutting most.
  5. Self-check: each guide must have exactly ONE perfect match inside lacZ
     (its intended site). Anything else is reported as an off-target.
  6. Combines the Step 3 score with an off-target penalty -> guides_final.csv.

Files needed in the same folder: fetch.py, ecoli_genome.gb, guides_shortlist.csv
Requires: numpy, pandas, biopython
"""

import numpy as np
import pandas as pd

from fetch import load_genome, get_gene      # reuse the code from Step 1

# ---- Settings you can change -------------------------------------------
GUIDE_LEN = 20
MAX_MISMATCH = 4      # report sites with up to this many mismatches
SEED_LEN = 12         # PAM-proximal bases that matter most
SEED_WEIGHT = 2       # a seed mismatch counts this many times
RISK_FACTOR = 10      # bigger = off-targets are punished harder
# -------------------------------------------------------------------------

LETTERS = np.array(list("ACGTN"))
LUT = np.full(256, 4, dtype=np.uint8)         # anything that is not A/C/G/T -> 4
for _code, _base in enumerate("ACGT"):
    LUT[ord(_base)] = _code


def encode(seq):
    """DNA string -> numpy array of numbers: A=0, C=1, G=2, T=3, other=4."""
    return LUT[np.frombuffer(seq.upper().encode(), dtype=np.uint8)]


def reverse_complement_codes(codes):
    """Reverse complement of an encoded sequence (A<->T, C<->G)."""
    rev = codes[::-1].astype(np.int8)
    return np.where(rev < 4, 3 - rev, 4).astype(np.uint8)


def find_sites(codes, strand, genome_len):
    """
    Find every NGG PAM on ONE strand and grab the 20 nt in front of it.
    Returns the sites plus their coordinates on the ORIGINAL genome (1-based).
    """
    # PAM starts at i when positions i+1 and i+2 are both G (code 2)
    is_gg = (codes[1:-1] == 2) & (codes[2:] == 2)
    pam_pos = np.nonzero(is_gg)[0]
    pam_pos = pam_pos[pam_pos >= GUIDE_LEN]

    all_windows = np.lib.stride_tricks.sliding_window_view(codes, GUIDE_LEN)
    win_start = pam_pos - GUIDE_LEN              # 0-based start on this strand
    windows = all_windows[win_start]             # shape (n_sites, 20)
    pam_n = codes[pam_pos]                       # the "N" of NGG

    if strand == "+":
        start = win_start + 1
        end = win_start + GUIDE_LEN
    else:                                        # convert back to original coordinates
        start = genome_len - win_start - GUIDE_LEN + 1
        end = genome_len - win_start
    return windows, pam_n, start, end


def build_site_table(genome_seq):
    """Table of all NGG sites on both strands of the genome."""
    codes = encode(genome_seq)
    n = len(codes)
    w1, p1, s1, e1 = find_sites(codes, "+", n)
    w2, p2, s2, e2 = find_sites(reverse_complement_codes(codes), "-", n)
    return {
        "windows": np.vstack([w1, w2]),
        "pam_n": np.concatenate([p1, p2]),
        "start": np.concatenate([s1, s2]),
        "end": np.concatenate([e1, e2]),
        "strand": np.array(["+"] * len(w1) + ["-"] * len(w2)),
    }


def screen_guide(guide, sites, gene_start, gene_end, max_mm=MAX_MISMATCH):
    """
    Compare one guide against every site in the genome.
    Returns (number_of_on_target_sites, list_of_off_target_hits).
    """
    g = encode(guide)
    mism = sites["windows"] != g                 # True wherever a base differs
    n_mm = mism.sum(axis=1)                      # mismatches per site
    hit_idx = np.nonzero(n_mm <= max_mm)[0]

    on_target = 0
    hits = []
    for h in hit_idx:
        start, end = int(sites["start"][h]), int(sites["end"][h])
        k = int(n_mm[h])

        # A perfect match inside lacZ is the intended site, not an off-target
        if k == 0 and start >= gene_start and end <= gene_end:
            on_target += 1
            continue

        m = mism[h]
        seed = int(m[GUIDE_LEN - SEED_LEN:].sum())      # last 12 = next to PAM
        distal = k - seed
        weighted = SEED_WEIGHT * seed + distal

        letters = LETTERS[sites["windows"][h]]
        site_seq = "".join(c.lower() if bad else c for c, bad in zip(letters, m))

        hits.append({
            "guide": guide,
            "hit_start": start,
            "hit_end": end,
            "hit_strand": str(sites["strand"][h]),
            "hit_pam": str(LETTERS[sites["pam_n"][h]]) + "GG",
            "hit_seq": site_seq,                        # mismatches in lowercase
            "mismatches": k,
            "seed_mm": seed,
            "distal_mm": distal,
            "mm_positions": ",".join(str(i + 1) for i in np.nonzero(m)[0]),
            "weighted": weighted,
            "hit_weight": 0.5 ** weighted,              # 1.0 = as bad as a perfect match
        })
    return on_target, hits


if __name__ == "__main__":
    genome = load_genome()
    genome_seq = str(genome.seq)
    _, gene0, gene_end, _ = get_gene(genome)
    gene_start = gene0 + 1                               # 1-based, inclusive

    print("Building site table for both strands...")
    sites = build_site_table(genome_seq)
    print(f"Genome: {len(genome_seq):,} bp | candidate NGG sites: {len(sites['start']):,}")

    shortlist = pd.read_csv("guides_shortlist.csv")
    all_hits, summaries = [], []

    for row in shortlist.itertuples():
        on_target, hits = screen_guide(row.guide, sites, gene_start, gene_end)
        all_hits.extend(hits)

        counts = {k: sum(1 for h in hits if h["mismatches"] == k)
                  for k in range(MAX_MISMATCH + 1)}
        risk = sum(h["hit_weight"] for h in hits)
        specificity = 1 / (1 + RISK_FACTOR * risk)

        summaries.append({
            "guide": row.guide,
            "on_target_sites": on_target,
            **{f"off_{k}mm": counts[k] for k in range(MAX_MISMATCH + 1)},
            "ot_risk": round(risk, 3),
            "specificity": round(specificity, 3),
        })

        flag = "OK" if on_target == 1 else "CHECK!"
        print(f"{row.guide}  self-check {flag}  "
              f"off-targets (0..{MAX_MISMATCH} mm): "
              f"{[counts[k] for k in range(MAX_MISMATCH + 1)]}")

    summary_df = pd.DataFrame(summaries)
    cols = ["guide", "pam", "strand", "start", "end", "gc", "score"]
    final = shortlist[cols].merge(summary_df, on="guide")
    final["final_score"] = (final["score"] * final["specificity"]).round(1)
    final = final.sort_values("final_score", ascending=False).reset_index(drop=True)

    print("\nTop 10 guides after off-target screening:\n")
    show = ["guide", "strand", "start", "end", "score", "ot_risk", "final_score"]
    print(final[show].head(10).to_string(index=False))

    final.to_csv("guides_final.csv", index=False)
    pd.DataFrame(all_hits).to_csv("offtarget_hits.csv", index=False)
    print("\nSaved guides_final.csv and offtarget_hits.csv")
    
