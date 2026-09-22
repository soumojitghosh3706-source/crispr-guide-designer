"""
fetch.py  --  Step 1 of CRISPR-GuideDesigner (Phase 1)

What it does:
  1. Downloads the E. coli K-12 MG1655 genome from NCBI (only once).
  2. Saves it locally as a GenBank file.
  3. Finds the lacZ gene and extracts its sequence.
  4. Saves lacZ to lacZ.fasta and prints a quick sanity check.

Needs internet the first time. Requires: biopython
"""

import os
from Bio import Entrez, SeqIO

# NCBI asks for a real email so they can contact you if something goes wrong.
Entrez.email = "soumojitghosh3706@gmail.com"   # <-- CHANGE THIS

GENOME_ID = "NC_000913.3"        # E. coli K-12 MG1655 reference genome
GENOME_FILE = "ecoli_genome.gb"  # local copy, so we download only once
TARGET_GENE = "lacZ"


def download_genome(accession=GENOME_ID, out_file=GENOME_FILE):
    """Download the genome as a GenBank file, unless it already exists."""
    if os.path.exists(out_file):
        print(f"Using existing file: {out_file}")
        return
    print("Downloading genome from NCBI (about 10-15 MB, please wait)...")
    handle = Entrez.efetch(
        db="nucleotide", id=accession, rettype="gbwithparts", retmode="text"
    )
    data = handle.read()
    handle.close()
    with open(out_file, "w") as f:
        f.write(data)
    print(f"Saved to {out_file}")


def load_genome(gb_file=GENOME_FILE):
    """Read the GenBank file into a Biopython SeqRecord."""
    return SeqIO.read(gb_file, "genbank")


def get_gene(record, gene_name=TARGET_GENE):
    """
    Find a gene by name in the genome annotations.
    Returns (sequence, start, end, strand).
      - sequence is given in the gene's own (coding) direction
      - start/end are 0-based genome coordinates (start inclusive, end exclusive)
      - strand is +1 or -1
    """
    for feature in record.features:
        if feature.type == "gene" and gene_name in feature.qualifiers.get("gene", []):
            seq = feature.extract(record.seq)  # reverse-complements automatically if needed
            start = int(feature.location.start)
            end = int(feature.location.end)
            strand = feature.location.strand
            return seq, start, end, strand
    raise ValueError(f"Gene '{gene_name}' not found in the record")


if __name__ == "__main__":
    download_genome()

    genome = load_genome()
    print(f"Genome: {genome.id}, length {len(genome.seq):,} bp")

    gene_seq, start, end, strand = get_gene(genome)
    strand_symbol = "+" if strand == 1 else "-"

    print(f"{TARGET_GENE}: {len(gene_seq)} bp")
    print(f"Genome position: {start + 1}-{end} (1-based), strand {strand_symbol}")
    print(f"First 60 bp: {str(gene_seq[:60])}")
    print(f"Starts with ATG? {str(gene_seq[:3]) == 'ATG'}")

    with open("lacZ.fasta", "w") as f:
        f.write(f">{TARGET_GENE}\n{str(gene_seq)}\n")
    print("Saved lacZ.fasta")
    
