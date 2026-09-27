"""
mutyh_fetch.py  --  Case study 3, Step 1 of CRISPR-GuideDesigner

Target: human MUTYH (mutY DNA glycosylase, base-excision DNA repair;
mutations cause MUTYH-associated polyposis, a hereditary colorectal cancer
syndrome), chromosome 1, GRCh38.

Same logic as phase2_fetch.py for PIP4K2C, retargeted to MUTYH.

What it does:
  1. Downloads a window around MUTYH from NCBI straight into memory (no
     local GenBank file needed -- set SAVE_GENBANK = True to keep a copy).
  2. Finds the gene and all of its protein-coding isoforms (CDS features).
  3. Picks the longest isoform, unless PREFERRED_PROTEIN is set below.
  4. Builds a table of its coding exons: position inside the gene, on
     chromosome 1, and inside the coding sequence, plus whether each exon
     is shared by every isoform / every curated (NP_) isoform.
  5. Self-checks: the spliced coding sequence must match NCBI's own, start
     with ATG, be a multiple of 3, and translate to one complete protein.
  6. Saves mutyh_gene.fasta and mutyh_exons.csv.

Needs internet every run. Requires: biopython, pandas
"""

import io
import os

import pandas as pd
from Bio import Entrez, SeqIO
from Bio.Seq import Seq

Entrez.email = "your_email@example.com"   # <-- CHANGE THIS

# ---- Settings ----------------------------------------------------------
GENE_NAME = "MUTYH"
ACCESSION = "NC_000001.11"     # human chromosome 1 (GRCh38)
SEQ_START = 45_325_000         # window around the gene, with margin
SEQ_STOP = 45_345_000
PREFERRED_PROTEIN = None       # e.g. "NP_036354" to force a specific isoform
SAVE_GENBANK = False            # True = also keep the downloaded GenBank file
GENBANK_COPY = "mutyh_region.gb"
# -------------------------------------------------------------------------


def fetch_region():
    """Download the chromosome window from NCBI and return it as a Biopython record."""
    if Entrez.email == "your_email@example.com":
        print("NOTE: put your own email address in Entrez.email at the top of this file.")
    print("Downloading the region from NCBI (about 20 kb, a few seconds)...")
    handle = Entrez.efetch(
        db="nuccore", id=ACCESSION, rettype="gb", retmode="text",
        seq_start=SEQ_START, seq_stop=SEQ_STOP, strand=1,
    )
    data = handle.read()
    handle.close()
    if not data.startswith("LOCUS"):
        raise RuntimeError("NCBI did not return a GenBank record. Reply was:\n" + data[:300])
    if SAVE_GENBANK:
        with open(GENBANK_COPY, "w") as f:
            f.write(data)
        print(f"Also saved a copy: {os.path.abspath(GENBANK_COPY)}")
    return SeqIO.read(io.StringIO(data), "genbank")


# ---------- coordinate helpers (all intervals are 0-based, end exclusive) ---

def to_gene_coords(start, end, gene_start, gene_end, strand):
    if strand == 1:
        return start - gene_start, end - gene_start
    return gene_end - end, gene_end - start


def to_region_coords(start, end, gene_start, gene_end, strand):
    if strand == 1:
        return gene_start + start, gene_start + end
    return gene_end - end, gene_end - start


def overlaps(a, b):
    return a[0] < b[1] and b[0] < a[1]


def feature_parts(feature):
    return [(int(p.start), int(p.end)) for p in feature.location.parts]


def find_gene(record, name):
    for f in record.features:
        if f.type == "gene" and name in f.qualifiers.get("gene", []):
            return f
    found = sorted({q for f in record.features if f.type == "gene"
                    for q in f.qualifiers.get("gene", [])})
    raise ValueError(f"{name} not found. Genes in this window: {', '.join(found)}")


def collect_isoforms(record, gene_name, gene_start, gene_end, strand):
    isoforms = {}
    for f in record.features:
        if f.type != "CDS" or gene_name not in f.qualifiers.get("gene", []):
            continue
        parts = sorted(
            to_gene_coords(s, e, gene_start, gene_end, strand) for s, e in feature_parts(f)
        )
        entry = isoforms.setdefault(
            tuple(parts), {"parts": parts, "proteins": [], "feature": f}
        )
        entry["proteins"].append(f.qualifiers.get("protein_id", ["?"])[0])
    for entry in isoforms.values():
        entry["curated"] = any(p.startswith("NP_") for p in entry["proteins"])
    return list(isoforms.values())


def choose_canonical(isoforms, preferred=None):
    if preferred:
        for iso in isoforms:
            if any(p.startswith(preferred) for p in iso["proteins"]):
                return iso
    return max(isoforms, key=lambda iso: sum(e - s for s, e in iso["parts"]))


def build_exon_table(canonical, isoforms, gene_start, gene_end, strand, seq_start):
    rows = []
    cds_pos = 0
    curated = [iso for iso in isoforms if iso["curated"]] or isoforms
    for i, (s, e) in enumerate(canonical["parts"], start=1):
        length = e - s
        r_start, r_end = to_region_coords(s, e, gene_start, gene_end, strand)
        shared = all(
            any(overlaps((s, e), q) for q in iso["parts"]) for iso in isoforms
        )
        shared_cur = all(
            any(overlaps((s, e), q) for q in iso["parts"]) for iso in curated
        )
        rows.append({
            "coding_exon": i,
            "gene_start": s + 1,
            "gene_end": e,
            "length": length,
            "chr_start": seq_start + r_start,
            "chr_end": seq_start + r_end - 1,
            "cds_start": cds_pos + 1,
            "cds_end": cds_pos + length,
            "shared_by_all_isoforms": shared,
            "shared_by_curated_isoforms": shared_cur,
        })
        cds_pos += length
    return pd.DataFrame(rows)


if __name__ == "__main__":
    print(f"Working folder: {os.getcwd()}")
    record = fetch_region()
    print(f"Window: {ACCESSION}:{SEQ_START:,}-{SEQ_STOP:,} ({len(record.seq):,} bp)")

    gene = find_gene(record, GENE_NAME)
    g_start, g_end = int(gene.location.start), int(gene.location.end)
    strand = gene.location.strand
    gene_seq = str(gene.extract(record.seq)).upper()
    strand_symbol = "+" if strand == 1 else "-"

    print(f"\n{GENE_NAME}: chr1:{SEQ_START + g_start:,}-{SEQ_START + g_end - 1:,} "
          f"(strand {strand_symbol}), gene length {len(gene_seq):,} bp")

    isoforms = collect_isoforms(record, GENE_NAME, g_start, g_end, strand)
    if not isoforms:
        raise SystemExit("No CDS features found for this gene in the window.")
    print(f"\nDistinct protein-coding isoforms: {len(isoforms)}")
    for iso in isoforms:
        cds_len = sum(e - s for s, e in iso["parts"])
        kind = "curated" if iso["curated"] else "predicted"
        print(f"  {', '.join(iso['proteins'])} [{kind}]: CDS {cds_len} nt, "
              f"{len(iso['parts'])} coding exons")

    canonical = choose_canonical(isoforms, PREFERRED_PROTEIN)
    print(f"\nUsing isoform: {', '.join(canonical['proteins'])}")

    exons = build_exon_table(canonical, isoforms, g_start, g_end, strand, SEQ_START)
    print("\nCoding exons (gene-oriented coordinates):\n")
    print(exons.to_string(index=False))

    cds_seq = "".join(gene_seq[s:e] for s, e in canonical["parts"])
    ncbi_cds = str(canonical["feature"].extract(record.seq)).upper()
    protein = Seq(cds_seq).translate(to_stop=True)
    print("\nSelf-checks:")
    print(f"  Our spliced CDS matches NCBI's own: {cds_seq == ncbi_cds}")
    print(f"  CDS length: {len(cds_seq)} nt (multiple of 3: {len(cds_seq) % 3 == 0})")
    print(f"  Starts with ATG: {cds_seq[:3] == 'ATG'}")
    print(f"  Protein: {len(protein)} aa, complete ORF: "
          f"{len(protein) * 3 + 3 == len(cds_seq)}, starts {str(protein)[:10]}...")

    with open("mutyh_gene.fasta", "w") as f:
        f.write(f">{GENE_NAME} gene, oriented 5'->3', "
                f"chr1:{SEQ_START + g_start}-{SEQ_START + g_end - 1} strand {strand_symbol}\n")
        f.write(gene_seq + "\n")
    exons.to_csv("mutyh_exons.csv", index=False)
    print("\nSaved these files:")
    for name in ("mutyh_gene.fasta", "mutyh_exons.csv"):
        print("  " + os.path.abspath(name))
        
