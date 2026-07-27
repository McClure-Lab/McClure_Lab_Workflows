#!/usr/bin/env python3

import pysam
import argparse
import subprocess
import sys
import os
from difflib import SequenceMatcher


def parse_args():
    """
    Parse command-line arguments for RFB motif detection.

    The script takes a BAM file and a genomic region, then searches reads in
    that region for an approximate match to the RFB motif.

    Returns
    -------
    argparse.Namespace
        Parsed arguments containing the BAM path, chromosome, start coordinate,
        end coordinate, and optional output file path.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("bam")
    parser.add_argument("-c", "--chrom", required=True)
    parser.add_argument("-s", "--start", required=True, type=int)
    parser.add_argument("-e", "--end", required=True, type=int)
    parser.add_argument("-o", "--output", default=None)
    return parser.parse_args()


def get_project_paths():
    """
    Build paths to the workflow-managed BAM directory.

    Sorted BAM files and BAM indexes are stored next to source BAMs in data/bam
    so users have only one BAM location to check.

    Returns
    -------
    tuple[str, str]
        A tuple containing:
        - bam_dir: directory for BAM files
        - index_dir: directory for BAM index files
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    bam_dir = os.path.abspath(os.path.join(script_dir, "../../data/bam"))
    return bam_dir, bam_dir


def check_and_sort_bam(bam_path, bam_dir):
    """
    Ensure the input BAM is coordinate-sorted and stored with a standard name.

    Region-based read fetching requires a coordinate-sorted BAM. This function
    checks the BAM header to see whether the file is already sorted. If it is,
    the input path is used directly. If it is not, samtools sort is used to
    create a sorted copy in data/bam.

    Parameters
    ----------
    bam_path : str
        Path to the input BAM file.

    bam_dir : str
        Directory where the standardized sorted BAM should be stored.

    Returns
    -------
    str
        Path to the standardized sorted BAM file.
    """
    bam = pysam.AlignmentFile(bam_path, "rb", check_sq=False)
    header = bam.header.to_dict()
    bam.close()

    sort_order = header.get("HD", {}).get("SO", "unknown")

    os.makedirs(bam_dir, exist_ok=True)

    # If the BAM is already coordinate-sorted, use the user supplied file
    # directly instead of redirecting to a renamed cache file.
    if sort_order == "coordinate":
        print(f"[INFO] Using provided sorted BAM: {bam_path}", file=sys.stderr)
        return bam_path

    # Avoid creating names like sample.sorted.indexed.sorted.indexed.bam when
    # the input BAM filename already has the workflow suffix.
    original_base = os.path.splitext(os.path.basename(bam_path))[0]
    if original_base.endswith(".sorted.indexed"):
        original_base = original_base[:-15]
    elif original_base.endswith(".sorted"):
        original_base = original_base[:-7]

    standardized_path = os.path.join(
        bam_dir,
        f"{original_base}.sorted.indexed.bam"
    )

    # Reuse an existing standardized BAM so repeated workflow runs do not waste
    # time copying or sorting the same file again.
    if (
        os.path.exists(standardized_path)
        and os.path.getmtime(standardized_path) >= os.path.getmtime(bam_path)
    ):
        print(f"[INFO] Reusing standardized BAM: {standardized_path}", file=sys.stderr)
        return standardized_path

    subprocess.run(["samtools", "sort", "-o", standardized_path, bam_path], check=True)
    return standardized_path


def check_and_index_bam(bam_path, index_dir):
    """
    Ensure the sorted BAM has a current index file.

    pysam needs a BAM index to fetch reads from a specific genomic region. This
    function creates the index if it does not exist, or rebuilds it if the BAM
    file is newer than the existing index.

    Parameters
    ----------
    bam_path : str
        Path to the coordinate-sorted BAM file.

    index_dir : str
        Directory where the BAM index should be stored.

    Returns
    -------
    str
        Path to the BAM index file.
    """
    index_path = os.path.join(index_dir, f"{os.path.basename(bam_path)}.bai")
    os.makedirs(index_dir, exist_ok=True)

    # Rebuild if there is no index, or if the BAM was modified after the index.
    # This prevents pysam from using a stale index.
    rebuild = (
        not os.path.exists(index_path)
        or os.path.getmtime(index_path) < os.path.getmtime(bam_path)
    )

    if rebuild:
        print("[INFO] Creating/Rebuilding index...", file=sys.stderr)
        subprocess.run(["samtools", "index", "-o", index_path, bam_path], check=True)
    else:
        print(f"[INFO] Index found: {index_path}", file=sys.stderr)

    return index_path


def similar(a, b):
    """
    Calculate sequence similarity between two strings.

    SequenceMatcher returns a value between 0 and 1, where 1 means the sequences
    match exactly. This is used to allow approximate motif matches instead of
    requiring a perfect RFB motif match.

    Parameters
    ----------
    a : str
        First sequence.

    b : str
        Second sequence.

    Returns
    -------
    float
        Similarity ratio between the two sequences.
    """
    return SequenceMatcher(None, a, b).ratio()


def extract_rfb_reads(bam_path, chrom, start, end):
    """
    Find reads containing an approximate RFB motif match.

    This function scans reads overlapping the requested genomic region and
    searches each read sequence for the RFB motif:

        TTTACCAAGAAAGATGTAAG

    The match threshold allows up to about two mismatches. If a read contains at
    least one approximate motif match, one BED-like row is returned for that
    read. The output uses the read's reference start and end positions so the
    motif-positive read can be overlaid on the rain plot.

    Parameters
    ----------
    bam_path : str
        Path to the sorted and indexed BAM file.

    chrom : str
        Chromosome or contig to search.

    start : int
        Start coordinate of the region.

    end : int
        End coordinate of the region.

    Returns
    -------
    list[tuple]
        A list of motif-positive reads. Each tuple contains:
        chrom, read_start, read_end, read_id, base_label, score.
    """
    motif = "TTTACCAAGAAAGATGTAAG"

    # Allow approximate RFB matches instead of requiring a perfect sequence.
    # The threshold is set so the motif can differ by about two bases.
    threshold = (len(motif) - 2) / len(motif)

    _, index_dir = get_project_paths()
    index_path = os.path.join(index_dir, f"{os.path.basename(bam_path)}.bai")

    bam = pysam.AlignmentFile(bam_path, "rb", index_filename=index_path)

    results = []

    for read in bam.fetch(chrom, start, end):
        if read.is_unmapped:
            continue

        seq = read.query_sequence
        if not seq:
            continue

        # Map read/query positions back to reference coordinates so motif
        # matches can be written at their genomic location rather than at the
        # read's full alignment span.
        aligned_pairs = dict(read.get_aligned_pairs(matches_only=True))

        # Slide a motif-sized window across the read sequence and compare each
        # window to the expected RFB motif.
        for i in range(len(seq) - len(motif) + 1):
            motif_window = seq[i:i + len(motif)]
            if similar(motif_window, motif) >= threshold:
                motif_query_positions = range(i, i + len(motif))
                motif_ref_positions = [
                    aligned_pairs[qpos]
                    for qpos in motif_query_positions
                    if qpos in aligned_pairs
                ]

                if not motif_ref_positions:
                    continue

                motif_start = min(motif_ref_positions)
                motif_end = max(motif_ref_positions) + 1

                results.append((
                    chrom,
                    motif_start,
                    motif_end,
                    read.query_name,
                    "T",
                    1.0
                ))

                # Only report one row per read. Once a motif is found, there is
                # no need to keep scanning that same read.
                break

    bam.close()
    return results


def main():
    """
    Run the RFB motif extraction workflow.

    This function prepares a sorted and indexed BAM, searches reads in the
    requested region for the RFB motif, and writes motif-positive reads in a
    BED-like format.

    Output columns are:
        chrom, start, end, read_id, base_label, score
    """
    args = parse_args()

    sorted_dir, index_dir = get_project_paths()

    bam_path = check_and_sort_bam(args.bam, sorted_dir)
    check_and_index_bam(bam_path, index_dir)

    # Search against the standardized sorted/indexed BAM so the BAM and BAI
    # naming conventions stay aligned across workflow runs.
    results = extract_rfb_reads(bam_path, args.chrom, args.start, args.end)

    print(f"[INFO] {len(results)} reads containing RFB motif found.", file=sys.stderr)

    if len(results) == 0:
        print("[INFO] No RFB found — expected for non S-phase datasets.", file=sys.stderr)
        return

    out = open(args.output, "w") if args.output else sys.stdout

    for r in results:
        # Write a BED-like row so the rain plot script can use this file as an
        # overlay track.
        out.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]:.7f}\n")

    if args.output:
        out.close()


if __name__ == "__main__":
    main()
