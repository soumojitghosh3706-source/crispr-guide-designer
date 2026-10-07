"""
primer_design.py  --  Phase 4, Step 3 of CRISPR-GuideDesigner

For each of the final validated guides (final_guides.csv) this script:
  1. Downloads the same NCBI region your fetch scripts use.
  2. Locates the guide on the genome (either strand, must be followed by NGG)
     and works out its strand, PAM and Cas9 cut site (3 bp upstream of PAM).
  3. Designs a genotyping primer pair around the cut site:
       - product 600-800 bp
       - cut site OFF-centre, each T7E1 fragment >= 200 bp, fragments >= 100 bp
         different in size (so the two cleavage bands are clearly separate)
       - primers 20-24 nt, GC 40-60%, Tm 58-62 C, Tm difference <= 2 C
       - no runs of 4+ identical bases, no dinucleotide repeats,
         G/C at the 3' end, at most 3 G/C in the last 5 bases
       - no 3'-end self/cross-dimer match between the two primers
       - self- and cross-complementarity (primer-dimer / hairpin) scores low
         (Primer3-style score; any <= 8, 3' end <= 5, same for the F+R pair)
       - the primer's 3' 12-mer is UNIQUE within the downloaded region, and the
         primer is not similar to an Alu consensus (<=4 mismatches over 18 bases),
         and matches the region (<=3 mismatches, either strand) only at its own site
  4. Re-uses an existing amplicon when a later guide in the same gene already
     sits well inside it (fewer primers to order).
  5. Saves final_guides_located.csv, amplicons.csv and primer_order.csv.

Primer specificity (BLAST) is a separate next step -- this script only designs.

Needs internet every run. Requires: biopython
"""

import csv
import re
import time

from Bio import Entrez, SeqIO
from Bio.SeqUtils import MeltingTemp as mt

Entrez.email = "soumojitghosh3706@gmail.com"   # <-- CHANGE THIS

# ---- Settings ----------------------------------------------------------
GUIDE_FILE = "final_guides.csv"
GENES = {
    "PIP4K2C": dict(acc="NC_000012.12", start=57_585_000, stop=57_610_000, chrom="chr12"),
    "MUTYH":   dict(acc="NC_000001.11", start=45_325_000, stop=45_345_000, chrom="chr1"),
    "KRAS":    dict(acc="NC_000012.12", start=25_200_000, stop=25_255_000, chrom="chr12"),
    "G6PD":    dict(acc="NC_000023.11", start=154_520_000, stop=154_560_000, chrom="chrX"),
}
# To redesign only some guides (e.g. after a failed Primer-BLAST), list them here.
# Outputs then get a "_redesign" suffix so the first-round files are not overwritten.
ONLY_GUIDES = {"PIP4K2C_g1"}      # e.g. {"G6PD_g3", "G6PD_g4", "G6PD_g5", "PIP4K2C_g1"}
PROD_MIN, PROD_MAX, PROD_TARGET = 600, 800, 700
MIN_ARM, MIN_ARM_DIFF = 200, 100       # T7E1 fragment rules
TM_TARGET = 60.0
# Tier 1 = strict rules. If no pair is found (AT-rich introns are common),
# the script falls back to the next, more relaxed tier and says so in the output.
TIERS = [
    dict(name="strict",  lens=range(20, 25), tm=(58.0, 62.0), tm_diff=2.0, gc=(40.0, 60.0), self_any=8, self_end=5),
    dict(name="relaxed", lens=range(20, 28), tm=(56.0, 63.0), tm_diff=3.0, gc=(35.0, 65.0), self_any=9, self_end=6),
    # Last resort for repeat-rich regions: look further from the cut and allow a longer product.
    dict(name="wide", lens=range(20, 28), tm=(56.0, 63.0), tm_diff=3.0, gc=(35.0, 65.0), self_any=9, self_end=6,
         reach=700, prod=(600, 1000)),
]
# -------------------------------------------------------------------------

COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(COMP)[::-1]


def gc_percent(s):
    return 100.0 * sum(b in "GC" for b in s) / len(s)


def tm(s):
    return mt.Tm_NN(s, Na=50, Mg=1.5, dNTPs=0.2, saltcorr=7)


def primer_ok(s, cfg):
    """Single-primer filters. Returns Tm if it passes, otherwise None."""
    if set(s) - set("ACGT"):
        return None
    if not (cfg["gc"][0] <= gc_percent(s) <= cfg["gc"][1]):
        return None
    if re.search(r"(A{4}|C{4}|G{4}|T{4})", s):
        return None
    if re.search(r"(..)\1\1", s):                       # e.g. ATATAT
        return None
    if s[-1] not in "GC" or sum(b in "GC" for b in s[-5:]) > 3:
        return None
    t = tm(s)
    return t if cfg["tm"][0] <= t <= cfg["tm"][1] else None


PAIR = {"A": "T", "T": "A", "G": "C", "C": "G"}

# Approximate Alu consensus (written from memory -- best effort). Primer-BLAST showed that
# primers copied from Alu elements match hundreds of sites, so any primer sharing a
# 12-base stretch with this sequence (either strand) is rejected.
ALU_CONSENSUS = ("GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACCTGAGGTCAGGAG"
                 "TTCGAGACCAGCCTGGCCAACATGGTGAAACCCCGTCTCTACTAAAAATACAAAAATTAGCCGGGCGTGGTGGCGCG"
                 "CGCCTGTAATCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATCGCTTGAACCCGGGAGGCGGAGGTTGCAGTGAGCC"
                 "GAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAAAA")
_ALU_12MERS = set()
for _strand in (ALU_CONSENSUS, ALU_CONSENSUS.translate(str.maketrans("ACGT", "TGCA"))[::-1]):
    for _i in range(len(_strand) - 11):
        _ALU_12MERS.add(_strand[_i:_i + 12])


def alu_like(p):
    """True if the primer shares any 12-mer with the Alu consensus."""
    return any(p[i:i + 12] in _ALU_12MERS for i in range(len(p) - 11))


_ALU_18MERS = []
for _strand in (ALU_CONSENSUS, ALU_CONSENSUS.translate(str.maketrans("ACGT", "TGCA"))[::-1]):
    _ALU_18MERS += [_strand[_i:_i + 18] for _i in range(len(_strand) - 17)]


def alu_similar(p, max_mm=4, w=18):
    """Approximate Alu check: any 18-base window of the primer within max_mm mismatches
    of the consensus (Alu copies differ from it, so exact matching misses many).
    Measured on real Primer-BLAST results: Alu-derived primers score 0-4, good ones 6-8."""
    for i in range(len(p) - w + 1):
        win = p[i:i + w]
        for a in _ALU_18MERS:
            if sum(x != y for x, y in zip(win, a)) <= max_mm:
                return True
    return False


def build_index(seq, k=5):
    idx = {}
    for i in range(len(seq) - k + 1):
        idx.setdefault(seq[i:i + k], []).append(i)
    return idx


def similar_sites(seq, idx, p, max_mm=3):
    """Number of places in the region (both strands) where the whole primer matches with
    <= max_mm mismatches. A unique primer gives exactly 1 (its own site); repeat-derived
    primers give more. Uses 5-base seeds from 4 pieces of the primer, then verifies."""
    n = len(p)
    hits = set()
    for strand, q in (("+", p), ("-", revcomp(p))):
        for c in range(4):
            off = c * (n // 4)
            for pos in idx.get(q[off:off + 5], ()):
                start = pos - off
                if start < 0 or start + n > len(seq):
                    continue
                if sum(x != y for x, y in zip(q, seq[start:start + n])) <= max_mm:
                    hits.add((strand, start))
    return len(hits)


def complementarity(a, b):
    """Primer3-style ungapped score of primer a binding primer b (anti-parallel).
    Returns (any, end): best local score, and best score that includes a's 3' base.
    +1 per complementary pair, -1 per mismatch. a == b gives self-dimer/hairpin risk."""
    n, m = len(a), len(b)
    best_any = best_end = 0
    for k in range(n + m - 1):                      # pairs (i, k - i)
        lo, hi = max(0, k - (m - 1)), min(n - 1, k)
        sc = [1 if PAIR[a[i]] == b[k - i] else -1 for i in range(lo, hi + 1)]
        cur = 0
        for v in sc:
            cur = max(0, cur + v)
            best_any = max(best_any, cur)
        if hi == n - 1:                             # diagonal reaches a's 3' end
            tot, best = 0, -9
            for v in reversed(sc):
                tot += v
                best = max(best, tot)
            best_end = max(best_end, best)
    return best_any, best_end


def dimer_ok(a, b, cfg):
    """False if the primers (or either one alone) can pair too strongly."""
    for x, y in ((a, a), (b, b), (a, b), (b, a)):
        sa, se = complementarity(x, y)
        if sa > cfg["self_any"] or se > cfg["self_end"]:
            return False
    return True


def three_prime_dimer(a, b):
    """Kept for the old simple check: last 5 bases of one primer pair inside either primer."""
    for x, y in ((a, b), (a, a), (b, b)):
        if revcomp(x[-5:]) in y:
            return True
    return False


def locate_guide(seq, guide):
    """All places the guide sits next to an NGG PAM, on either strand.
    Returns (strand, protospacer_start_index, PAM, cut_index).
    cut_index = boundary: bases seq[:cut_index] are left of the cut."""
    hits = []
    for m in re.finditer("(?=%s)" % guide, seq):
        i = m.start()
        pam = seq[i + 20:i + 23]
        if len(pam) == 3 and pam[1:] == "GG":
            hits.append(("+", i, pam, i + 17))
    rc = revcomp(guide)
    for m in re.finditer("(?=%s)" % rc, seq):
        j = m.start()
        if j >= 3 and seq[j - 3:j - 1] == "CC":
            hits.append(("-", j, revcomp(seq[j - 3:j]), j + 3))
    return hits


def candidates(seq, lo, hi, forward, cfg):
    """Passing primers. Forward: start s. Reverse: end e (exclusive)."""
    out = []
    for pos in range(max(lo, 0), min(hi, len(seq))):
        for n in cfg["lens"]:
            if forward:
                if pos + n > len(seq):
                    continue
                p = seq[pos:pos + n]
            else:
                if pos - n < 0:
                    continue
                p = revcomp(seq[pos - n:pos])
            t = primer_ok(p, cfg)
            if t is None:
                continue
            sa, se = complementarity(p, p)
            if sa > cfg["self_any"] or se > cfg["self_end"]:
                continue
            if alu_like(p):                          # Alu-derived primers match hundreds of sites
                continue
            k = p[-12:]                              # 3' 12-mer must be unique in the region
            if seq.count(k) + seq.count(revcomp(k)) != 1:
                continue
            out.append((pos, p, t))
    return out


def design_pair(seq, cut, idx):
    """Try each tier in turn; return the best primer pair of the first tier that works."""
    repeat_cache = {}

    def repeat_free(p):
        if p not in repeat_cache:
            repeat_cache[p] = (not alu_similar(p)) and similar_sites(seq, idx, p) == 1
        return repeat_cache[p]

    notes = []
    for tier_no, cfg in enumerate(TIERS, start=1):
        reach = cfg.get("reach", 450)
        pmin, pmax = cfg.get("prod", (PROD_MIN, PROD_MAX))
        fwd = candidates(seq, cut - reach, cut - MIN_ARM + 1, True, cfg)
        rev = candidates(seq, cut + MIN_ARM, cut + reach + 1, False, cfg)
        pairs = []
        for s, pf, tf in fwd:
            for e, pr, tr in rev:
                prod, left, right = e - s, cut - s, e - cut
                if not (pmin <= prod <= pmax):
                    continue
                if min(left, right) < MIN_ARM or abs(left - right) < MIN_ARM_DIFF:
                    continue
                if abs(tf - tr) > cfg["tm_diff"]:
                    continue
                score = (abs(prod - PROD_TARGET) / 100 + abs(min(left, right) - 300) / 100
                         + abs(tf - TM_TARGET) + abs(tr - TM_TARGET) + abs(tf - tr)
                         + 0.2 * (abs(len(pf) - 22) + abs(len(pr) - 22)))
                pairs.append((score, s, e, pf, pr, tf, tr))
        pairs.sort(key=lambda x: x[0])
        dimer_fail = repeat_fail = 0
        for score, s, e, pf, pr, tf, tr in pairs:
            if three_prime_dimer(pf, pr) or not dimer_ok(pf, pr, cfg):
                dimer_fail += 1
                continue
            if not (repeat_free(pf) and repeat_free(pr)):
                repeat_fail += 1
                continue
            return dict(s=s, e=e, fwd=pf, rev=pr, tm_f=tf, tm_r=tr, tier=cfg["name"])
        notes.append(f"    tier {cfg['name']}: {len(fwd)} forward / {len(rev)} reverse primers passed single-primer "
                     f"rules; {len(pairs)} pairs fit size rules; {dimer_fail} lost to dimers, "
                     f"{repeat_fail} lost to repeat checks")
    print("\n".join(notes))
    return None


def design_gene(gene, seq, seq_start, chrom, guides):
    """guides: list of (guide_id, guide_seq). Returns (located_rows, amplicon_rows)."""
    located, amplicons, designed = [], [], []
    kmer_idx = build_index(seq)
    for gid, g in guides:
        if ONLY_GUIDES and gid not in ONLY_GUIDES:
            continue
        hits = locate_guide(seq, g)
        if len(hits) != 1:
            print(f"  {gid}: found {len(hits)} PAM-adjacent sites (need exactly 1) -- skipped")
            located.append(dict(guide_id=gid, gene=gene, guide_seq=g, strand="", pam="",
                                chrom=chrom, cut_chr_pos="", note=f"{len(hits)} sites found"))
            continue
        strand, idx, pam, cut = hits[0]
        located.append(dict(guide_id=gid, gene=gene, guide_seq=g, strand=strand, pam=pam,
                            chrom=chrom, cut_chr_pos=seq_start + cut - 1, note=""))
        # reuse an existing amplicon of this gene if the cut sits well inside it
        reuse = None
        for d in designed:
            left, right = cut - d["s"], d["e"] - cut
            if min(left, right) >= MIN_ARM and abs(left - right) >= MIN_ARM_DIFF:
                reuse = d
                break
        if reuse:
            d, shared = reuse, reuse["owner"]
        else:
            d = design_pair(seq, cut, kmer_idx)
            if d is None:
                print(f"  {gid}: no primer pair met all rules -- relax settings or check region")
                amplicons.append(dict(guide_id=gid, gene=gene, note="no pair found"))
                continue
            d["owner"] = gid
            d["names"] = (f"{gene}_{gid.split('_')[-1]}_F", f"{gene}_{gid.split('_')[-1]}_R")
            designed.append(d)
            shared = ""
        left, right = cut - d["s"], d["e"] - cut
        amplicons.append(dict(
            guide_id=gid, gene=gene, strand=strand, pam=pam, chrom=chrom,
            cut_chr_pos=seq_start + cut - 1,
            amplicon_chr_start=seq_start + d["s"], amplicon_chr_end=seq_start + d["e"] - 1,
            product_bp=d["e"] - d["s"], t7e1_left_bp=left, t7e1_right_bp=right,
            fwd_name=d["names"][0], rev_name=d["names"][1],
            fwd_seq=d["fwd"], rev_seq=d["rev"],
            tm_fwd=round(d["tm_f"], 1), tm_rev=round(d["tm_r"], 1),
            shares_primers_with=shared,
            note=("relaxed primer rules used" if d["tier"] == "relaxed" else "")))
        print(f"  {gid}: {strand} strand, PAM {pam}, amplicon {d['e'] - d['s']} bp, "
              f"fragments {left}+{right}, rules: {d['tier']}"
              + (f" (shares primers with {shared})" if shared else ""))
    return located, amplicons


def fetch_seq(acc, start, stop):
    h = Entrez.efetch(db="nuccore", id=acc, rettype="gb", retmode="text",
                      seq_start=start, seq_stop=stop, strand=1)
    rec = SeqIO.read(h, "genbank")
    h.close()
    return str(rec.seq).upper()


def write_csv(path, rows, cols):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    if Entrez.email == "your_email@example.com":
        print("NOTE: put your own email address in Entrez.email at the top of this file.")
    guides = list(csv.DictReader(open(GUIDE_FILE)))
    all_loc, all_amp = [], []
    for gene, cfg in GENES.items():
        gl = [(r["guide_id"], r["guide_seq"].strip().upper()) for r in guides if r["gene"] == gene
              and (not ONLY_GUIDES or r["guide_id"] in ONLY_GUIDES)]
        if not gl:
            continue
        print(f"\n{gene}: downloading {cfg['acc']}:{cfg['start']:,}-{cfg['stop']:,} ...")
        seq = fetch_seq(cfg["acc"], cfg["start"], cfg["stop"])
        time.sleep(0.5)
        loc, amp = design_gene(gene, seq, cfg["start"], cfg["chrom"], gl)
        all_loc += loc
        all_amp += amp

    sfx = "_redesign" if ONLY_GUIDES else ""
    write_csv(f"final_guides_located{sfx}.csv", all_loc,
              ["guide_id", "gene", "guide_seq", "strand", "pam", "chrom", "cut_chr_pos", "note"])
    amp_cols = ["guide_id", "gene", "strand", "pam", "chrom", "cut_chr_pos",
                "amplicon_chr_start", "amplicon_chr_end", "product_bp",
                "t7e1_left_bp", "t7e1_right_bp", "fwd_name", "rev_name",
                "fwd_seq", "rev_seq", "tm_fwd", "tm_rev", "shares_primers_with", "note"]
    write_csv(f"amplicons{sfx}.csv", all_amp, amp_cols)

    seen, primers = set(), []
    for a in all_amp:
        if a.get("shares_primers_with") or "fwd_seq" not in a:
            continue
        for nm, sq, t in ((a["fwd_name"], a["fwd_seq"], a["tm_fwd"]),
                          (a["rev_name"], a["rev_seq"], a["tm_rev"])):
            if sq not in seen:
                seen.add(sq)
                primers.append(dict(primer_name=nm, sequence_5to3=sq, length=len(sq),
                                    gc_percent=round(gc_percent(sq), 1), tm_c=t))
    write_csv(f"primer_order{sfx}.csv", primers,
              ["primer_name", "sequence_5to3", "length", "gc_percent", "tm_c"])
    print(f"\nSaved final_guides_located{sfx}.csv, amplicons{sfx}.csv, primer_order{sfx}.csv "
          f"({len(primers)} primers to order)")
          
