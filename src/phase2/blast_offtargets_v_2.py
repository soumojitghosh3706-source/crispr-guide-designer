"""
blast_offtarget.py  --  Phase 2, Step 4 of CRISPR-GuideDesigner

GENOMIC-ONLY off-target screen for the human PIP4K2C shortlist.

What it does, for EACH shortlisted guide:
  1. Submits a short BLAST search against NCBI nucleotide records.
  2. Keeps ONLY genomic accessions:
       NC_, NT_, NW_, NG_, AC_
     and discards transcript/mRNA/cDNA accessions such as:
       NM_, XM_, XR_, NR_, BC_, AK_
  3. Keeps full-length, ungapped hits with <= MAX_MISMATCH mismatches.
  4. Excludes the guide's own PIP4K2C genomic locus.
  5. Fetches the 3 bases next to EVERY retained genomic hit.
  6. Checks for a canonical SpCas9 NGG PAM.
  7. Collapses overlapping hits from the SAME genomic accession into
     a single locus.
  8. Saves both raw genomic hits and unique-locus results.

IMPORTANT:
  - This is a genomic BLAST pre-screen, not a replacement for CRISPOR.
  - Transcript records are excluded entirely.
  - Different NCBI accessions are not assumed to represent the same
    physical genomic locus unless their coordinate system is demonstrably
    the same.
"""

import os
import re
import time

import pandas as pd
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML


# =====================================================================
# NCBI SETTINGS
# =====================================================================

Entrez.email = "your_email@example.com"   # <-- CHANGE THIS


# =====================================================================
# BLAST SETTINGS
# =====================================================================

# Keep nt for broad NCBI coverage, then strictly post-filter to genomic
# accessions. If your NCBI setup supports refseq_genomic, that can be
# used instead for a more genomic-focused search.
DATABASE = "nt"

ORGANISM_FILTER = "Homo sapiens[Organism]"

WORD_SIZE = 7

EXPECT = 1000

HITLIST_SIZE = 100

MAX_MISMATCH = 4

PAUSE_BETWEEN_SEARCHES = 3


# =====================================================================
# GENOMIC ACCESSION FILTER
# =====================================================================

# Genomic accession prefixes to KEEP.
#
# NC_ = RefSeq chromosome
# NT_ = RefSeq genomic contig
# NW_ = RefSeq genomic scaffold
# NG_ = RefSeq genomic region / gene-oriented genomic record
# AC_ = GenBank genomic clone/contig
#
# Everything else is rejected by is_genomic_accession().

GENOMIC_PREFIXES = (
    "NC_",
    "NT_",
    "NW_",
    "NG_",
    "AC_",
)


# Explicit transcript / non-genomic prefixes.
# These are not strictly necessary because we already use an allow-list,
# but keeping them documented makes the logic obvious.

TRANSCRIPT_PREFIXES = (
    "NM_",
    "XM_",
    "XR_",
    "NR_",
    "BC_",
    "AK_",
)


# =====================================================================
# PIP4K2C GENOMIC LOCATION
# =====================================================================

GENE_ACCESSION_PREFIX = "NC_000012"

GENE_CHR_START = 57_591_192

GENE_CHR_END = 57_603_418


# =====================================================================
# FILES
# =====================================================================

SHORTLIST_FILE = "pip4k2c_shortlist.csv"

HITS_FILE = "pip4k2c_blast_hits_genomic.csv"

LOCUS_FILE = "pip4k2c_unique_loci.csv"

PROGRESS_FILE = "pip4k2c_blast_progress_genomic.csv"


# =====================================================================
# OUTPUT COLUMNS
# =====================================================================

HIT_COLUMNS = [
    "guide",
    "accession",
    "hit_title",
    "hit_start",
    "hit_end",
    "hit_strand",
    "mismatches",
    "pam",
    "has_pam",
    "locus_id",
    "duplicate_record_count",
]


LOCUS_COLUMNS = [
    "guide",
    "locus_id",
    "accession",
    "hit_title",
    "locus_start",
    "locus_end",
    "hit_strand",
    "best_mismatches",
    "pam",
    "has_pam",
    "duplicate_record_count",
]


PROGRESS_COLUMNS = [
    "guide",
    "status",
    "raw_blast_hits",
    "genomic_hits",
    "unique_loci",
    "pam_positive_loci",
    "note",
]


# =====================================================================
# ACCESSION FILTER
# =====================================================================

def is_genomic_accession(accession):
    """
    Return True ONLY for genomic accessions.

    Allowed:
        NC_, NT_, NW_, NG_, AC_

    Everything else is rejected.

    This is deliberately an allow-list rather than trying to identify
    transcripts by name.
    """

    accession = str(accession).strip().upper()

    if accession.startswith(TRANSCRIPT_PREFIXES):
        return False

    return accession.startswith(GENOMIC_PREFIXES)


# =====================================================================
# SELF-HIT DETECTION
# =====================================================================

def is_self_hit(accession, sbjct_start, sbjct_end):
    """
    True if a genomic BLAST hit falls inside PIP4K2C's known genomic
    location.

    Only NC_000012 is used for the current PIP4K2C reference location.
    """

    accession = str(accession).strip().upper()

    if not accession.startswith(GENE_ACCESSION_PREFIX):
        return False

    lo, hi = sorted((sbjct_start, sbjct_end))

    return (
        lo >= GENE_CHR_START - 50
        and hi <= GENE_CHR_END + 50
    )


# =====================================================================
# PAM LOCATION
# =====================================================================

def pam_fetch_params(sbjct_start, sbjct_end):
    """
    Determine where the 3-nt SpCas9 PAM lies relative to the BLAST hit.

    Returns:
        fetch_start
        fetch_end
        Entrez strand

    Coordinates are 1-based.

    For a plus-strand hit:
        guide + PAM
        PAM is immediately downstream.

    For a minus-strand hit:
        PAM lies on the opposite genomic strand.
        Entrez strand=2 returns the reverse complement.
    """

    if sbjct_start < sbjct_end:

        # Plus-strand hit
        return (
            sbjct_end + 1,
            sbjct_end + 3,
            1
        )

    else:

        # Minus-strand hit
        return (
            sbjct_end - 3,
            sbjct_end - 1,
            2
        )


# =====================================================================
# FETCH PAM
# =====================================================================

def fetch_pam(accession, sbjct_start, sbjct_end):
    """
    Fetch the three bases immediately adjacent to the guide.

    Returns something like:
        AGG
        TGG
        CGG
        TCT
        etc.
    """

    fs, fe, strand = pam_fetch_params(
        sbjct_start,
        sbjct_end
    )

    if fs < 1:
        return None

    handle = Entrez.efetch(
        db="nuccore",
        id=accession,
        rettype="fasta",
        retmode="text",
        seq_start=fs,
        seq_stop=fe,
        strand=strand,
    )

    try:
        seq = str(
            SeqIO.read(
                handle,
                "fasta"
            ).seq
        ).upper()

    finally:
        handle.close()

    return seq


# =====================================================================
# MISMATCH CALCULATION
# =====================================================================

def hsp_mismatches(hsp, guide_len):
    """
    Count mismatches for a full-length ungapped HSP.

    Partial or gapped alignments are excluded.
    """

    if hsp.align_length != guide_len:
        return None

    if hsp.gaps:
        return None

    return hsp.align_length - hsp.identities


# =====================================================================
# BLAST SEARCH
# =====================================================================

def blast_one_guide(guide):
    """
    Submit one guide to NCBI BLAST.
    """

    result_handle = NCBIWWW.qblast(
        "blastn",
        DATABASE,
        guide,
        entrez_query=ORGANISM_FILTER,
        word_size=WORD_SIZE,
        expect=EXPECT,
        hitlist_size=HITLIST_SIZE,
    )

    try:
        record = NCBIXML.read(result_handle)

    finally:
        result_handle.close()

    return record


# =====================================================================
# CREATE LOCUS KEY
# =====================================================================

def make_locus_id(
    accession,
    hit_start,
    hit_end,
):
    """
    Create a stable locus identifier.

    IMPORTANT:
    Coordinates from different accessions are NOT assumed to be
    interchangeable.

    Therefore the accession is part of the locus key.
    """

    lo, hi = sorted(
        (
            int(hit_start),
            int(hit_end),
        )
    )

    return f"{accession}:{lo}-{hi}"


# =====================================================================
# OVERLAP TEST
# =====================================================================

def intervals_overlap(
    start1,
    end1,
    start2,
    end2,
):
    """
    Return True when two genomic intervals overlap.
    """

    return (
        max(start1, start2)
        <=
        min(end1, end2)
    )


# =====================================================================
# COLLAPSE SAME-ACCESSION OVERLAPPING HITS
# =====================================================================

def collapse_overlapping_loci(hits):
    """
    Collapse overlapping hits belonging to the SAME genomic accession.

    Multiple NCBI records can sometimes represent overlapping portions
    of the same genomic sequence.

    We therefore merge overlapping intervals when their accession is
    identical.

    Different accessions remain separate because their coordinates may
    belong to different coordinate systems.
    """

    if not hits:
        return []

    # Group by accession first.
    by_accession = {}

    for hit in hits:

        accession = hit["accession"]

        by_accession.setdefault(
            accession,
            []
        ).append(hit)

    unique_loci = []

    for accession, accession_hits in by_accession.items():

        accession_hits.sort(
            key=lambda x: (
                min(
                    x["hit_start"],
                    x["hit_end"]
                ),
                max(
                    x["hit_start"],
                    x["hit_end"]
                ),
            )
        )

        current = None

        for hit in accession_hits:

            start = min(
                hit["hit_start"],
                hit["hit_end"]
            )

            end = max(
                hit["hit_start"],
                hit["hit_end"]
            )

            if current is None:

                current = {
                    "guide": hit["guide"],
                    "accession": accession,
                    "hit_title": hit["hit_title"],
                    "locus_start": start,
                    "locus_end": end,
                    "hit_strand": hit["hit_strand"],
                    "best_mismatches": hit["mismatches"],
                    "pam": hit["pam"],
                    "has_pam": hit["has_pam"],
                    "duplicate_record_count": 1,
                }

                continue

            # Same accession, overlapping coordinate.
            if intervals_overlap(
                current["locus_start"],
                current["locus_end"],
                start,
                end,
            ):

                current["locus_start"] = min(
                    current["locus_start"],
                    start
                )

                current["locus_end"] = max(
                    current["locus_end"],
                    end
                )

                current["duplicate_record_count"] += 1

                # Keep the best mismatch count.
                if hit["mismatches"] < current["best_mismatches"]:

                    current["best_mismatches"] = hit["mismatches"]
                    current["pam"] = hit["pam"]
                    current["has_pam"] = hit["has_pam"]
                    current["hit_strand"] = hit["hit_strand"]

            else:

                unique_loci.append(current)

                current = {
                    "guide": hit["guide"],
                    "accession": accession,
                    "hit_title": hit["hit_title"],
                    "locus_start": start,
                    "locus_end": end,
                    "hit_strand": hit["hit_strand"],
                    "best_mismatches": hit["mismatches"],
                    "pam": hit["pam"],
                    "has_pam": hit["has_pam"],
                    "duplicate_record_count": 1,
                }

        if current is not None:
            unique_loci.append(current)

    # Add stable locus IDs.
    for locus in unique_loci:

        locus["locus_id"] = make_locus_id(
            locus["accession"],
            locus["locus_start"],
            locus["locus_end"],
        )

    return unique_loci


# =====================================================================
# FIND GENOMIC OFF-TARGETS
# =====================================================================

def find_offtargets(
    guide,
    record,
):
    """
    Process a BLAST result.

    Returns:

        raw_blast_hits
        genomic_candidates
        unique_loci
    """

    raw_blast_hits = 0

    genomic_candidates = []

    for alignment in record.alignments:

        accession = str(
            alignment.accession
        ).strip().upper()

        title = alignment.hit_def

        for hsp in alignment.hsps:

            raw_blast_hits += 1

            # ---------------------------------------------------------
            # GENOMIC ACCESSION FILTER
            # ---------------------------------------------------------

            if not is_genomic_accession(
                accession
            ):
                continue

            # ---------------------------------------------------------
            # MISMATCH FILTER
            # ---------------------------------------------------------

            mm = hsp_mismatches(
                hsp,
                len(guide)
            )

            if mm is None:
                continue

            if mm > MAX_MISMATCH:
                continue

            # ---------------------------------------------------------
            # REMOVE OWN TARGET
            # ---------------------------------------------------------

            if (
                mm == 0
                and
                is_self_hit(
                    accession,
                    hsp.sbjct_start,
                    hsp.sbjct_end,
                )
            ):
                continue

            # ---------------------------------------------------------
            # STORE GENOMIC HIT
            # ---------------------------------------------------------

            genomic_candidates.append(
                {
                    "guide": guide,
                    "accession": accession,
                    "hit_title": title,
                    "hit_start": hsp.sbjct_start,
                    "hit_end": hsp.sbjct_end,
                    "hit_strand":
                        "+"
                        if hsp.sbjct_start < hsp.sbjct_end
                        else "-",
                    "mismatches": mm,
                }
            )

    # Best sequence matches first.
    genomic_candidates.sort(
        key=lambda x: (
            x["mismatches"],
            x["accession"],
            min(
                x["hit_start"],
                x["hit_end"]
            ),
        )
    )

    # -------------------------------------------------------------
    # FETCH PAM FOR EVERY GENOMIC CANDIDATE
    # -------------------------------------------------------------

    checked_hits = []

    for hit in genomic_candidates:

        try:

            pam = fetch_pam(
                hit["accession"],
                hit["hit_start"],
                hit["hit_end"],
            )

        except Exception as e:

            pam = None

            print(
                "    PAM fetch failed for "
                f"{hit['accession']} "
                f"{hit['hit_start']}: {e}"
            )

        has_pam = (
            bool(pam)
            and len(pam) == 3
            and pam[1:] == "GG"
        )

        checked_hits.append(
            {
                **hit,
                "pam": pam,
                "has_pam": has_pam,
            }
        )

    # -------------------------------------------------------------
    # COLLAPSE OVERLAPPING LOCI
    # -------------------------------------------------------------

    unique_loci = collapse_overlapping_loci(
        checked_hits
    )

    return (
        raw_blast_hits,
        checked_hits,
        unique_loci,
    )


# =====================================================================
# LOAD PROGRESS
# =====================================================================

def load_progress():

    if os.path.exists(
        PROGRESS_FILE
    ):

        return pd.read_csv(
            PROGRESS_FILE
        )

    return pd.DataFrame(
        columns=PROGRESS_COLUMNS
    )


# =====================================================================
# LOAD RAW HITS
# =====================================================================

def load_hits():

    if os.path.exists(
        HITS_FILE
    ):

        return pd.read_csv(
            HITS_FILE
        )

    return pd.DataFrame(
        columns=HIT_COLUMNS
    )


# =====================================================================
# LOAD UNIQUE LOCI
# =====================================================================

def load_loci():

    if os.path.exists(
        LOCUS_FILE
    ):

        return pd.read_csv(
            LOCUS_FILE
        )

    return pd.DataFrame(
        columns=LOCUS_COLUMNS
    )


# =====================================================================
# MAIN
# =====================================================================

if __name__ == "__main__":

    shortlist = pd.read_csv(
        SHORTLIST_FILE
    )

    progress = load_progress()

    hits = load_hits()

    loci = load_loci()

    done_guides = set(
        progress.loc[
            progress["status"] == "done",
            "guide"
        ]
    )

    print(
        f"Working folder: {os.getcwd()}"
    )

    print(
        f"{len(shortlist)} guides in shortlist, "
        f"{len(done_guides)} already done.\n"
    )

    # -------------------------------------------------------------
    # GUIDE LOOP
    # -------------------------------------------------------------

    for i, row in enumerate(
        shortlist.itertuples(),
        start=1,
    ):

        guide = row.guide

        if guide in done_guides:

            print(
                f"[{i}/{len(shortlist)}] "
                f"{guide} already done, skipping"
            )

            continue

        print(
            f"[{i}/{len(shortlist)}] "
            f"{guide} submitting to NCBI BLAST..."
        )

        try:

            # -----------------------------------------------------
            # BLAST
            # -----------------------------------------------------

            record = blast_one_guide(
                guide
            )

            (
                raw_blast_hits,
                genomic_hits,
                unique_loci,
            ) = find_offtargets(
                guide,
                record,
            )

            status = "done"

            note = ""

            # -----------------------------------------------------
            # COUNTS
            # -----------------------------------------------------

            pam_positive = sum(
                bool(x["has_pam"])
                for x in unique_loci
            )

            print(
                f"    {raw_blast_hits} raw BLAST HSPs"
       )

            print(
                f"    {len(genomic_hits)} genomic "
                f"candidate hits"
            )

            print(
                f"    {len(unique_loci)} unique genomic loci"
            )

            print(
                f"    {pam_positive} unique loci "
                f"with canonical NGG PAM"
            )

        except Exception as e:

            raw_blast_hits = 0

            genomic_hits = []

            unique_loci = []

            status = "error"

            note = str(e)[:300]

            print(
                f"    ERROR: {note}"
            )

        # ---------------------------------------------------------
        # SAVE RAW GENOMIC HITS
        # ---------------------------------------------------------

        if genomic_hits:

            new_hits = pd.DataFrame(
                genomic_hits
            )

            hits = pd.concat(
                [
                    hits,
                    new_hits
                ],
                ignore_index=True,
            )

            hits.to_csv(
                HITS_FILE,
                index=False,
            )

        # ---------------------------------------------------------
        # SAVE UNIQUE LOCI
        # ---------------------------------------------------------

        if unique_loci:

            new_loci = pd.DataFrame(
                unique_loci
            )

            loci = pd.concat(
                [
                    loci,
                    new_loci
                ],
                ignore_index=True,
            )

            loci.to_csv(
                LOCUS_FILE,
                index=False,
            )

        # ---------------------------------------------------------
        # SAVE PROGRESS
        # ---------------------------------------------------------

        progress = pd.concat(
            [
                progress,
                pd.DataFrame(
                    [
                        {
                            "guide": guide,
                            "status": status,
                            "raw_blast_hits":
                                raw_blast_hits,
                            "genomic_hits":
                                len(genomic_hits),
                            "unique_loci":
                                len(unique_loci),
                            "pam_positive_loci":
                                sum(
                                    bool(x["has_pam"])
                                    for x in unique_loci
                                ),
                            "note": note,
                        }
                    ]
                ),
            ],
            ignore_index=True,
        )

        progress.to_csv(
            PROGRESS_FILE,
            index=False,
        )

        time.sleep(
            PAUSE_BETWEEN_SEARCHES
        )

    # =================================================================
    # FINAL SUMMARY
    # =================================================================

    n_done = int(
        (
            progress["status"]
            == "done"
        ).sum()
    )

    n_error = int(
        (
            progress["status"]
            == "error"
        ).sum()
    )

    total_genomic_hits = int(
        progress["genomic_hits"].sum()
    )

    total_unique_loci = int(
        progress["unique_loci"].sum()
    )

    total_pam_positive = int(
        progress["pam_positive_loci"].sum()
    )

    print(
        "\n=========================================="
    )

    print(
        "GENOMIC OFF-TARGET SCREEN COMPLETE"
    )

    print(
        "=========================================="
    )

    print(
        f"Guides completed: {n_done}"
    )

    print(
        f"Errors: {n_error}"
    )

    print(
        f"Genomic candidate hits: "
        f"{total_genomic_hits}"
    )

    print(
        f"Unique genomic loci: "
        f"{total_unique_loci}"
    )

    print(
        f"Canonical NGG-positive loci: "
        f"{total_pam_positive}"
    )

    print(
        "\nSaved files:"
    )

    print(
        f"  {os.path.abspath(PROGRESS_FILE)}"
    )

    print(
        f"  {os.path.abspath(HITS_FILE)}"
    )

    print(
        f"  {os.path.abspath(LOCUS_FILE)}"
    )