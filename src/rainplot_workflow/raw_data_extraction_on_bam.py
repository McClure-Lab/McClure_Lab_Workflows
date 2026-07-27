#!/usr/bin/env python3

import pysam
import argparse
import subprocess
import sys
import os
import re


def parse_args():
    """
    Parse command-line arguments for extracting BrdU modification data.

    The script takes a BAM file, a genomic region, and an optional read ID filter.
    It outputs per-base BrdU modification calls in a BED-like format.

    Returns
    -------
    argparse.Namespace
        Parsed arguments containing the BAM path, chromosome, start coordinate,
        end coordinate, optional read ID, and optional output file path.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Extract BrdU modification data from a BAM file in a specific region "
            "or for an exact read-ID list."
        )
    )
    parser.add_argument("bam")
    parser.add_argument("-c", "--chrom", default=None)
    parser.add_argument("-s", "--start", default=None, type=int)
    parser.add_argument("-e", "--end", default=None, type=int)
    parser.add_argument("-r", "--read_id", default=None)
    parser.add_argument(
        "--read_ids_file",
        default=None,
        help=(
            "Optional TSV or plain-text file containing read IDs to include. "
            "TSV files default to the read_id column."
        ),
    )
    parser.add_argument(
        "--read_id_column",
        default="read_id",
        help="Column name to read from --read_ids_file when it has a header.",
    )
    parser.add_argument(
        "--read_ids_limit",
        type=int,
        default=None,
        help="Maximum number of read IDs to load from --read_ids_file.",
    )
    parser.add_argument("-o", "--output", default=None)
    return parser.parse_args()


def load_read_ids(read_ids_file, read_id_column="read_id", limit=None):
    """
    Load read IDs from a TSV/CSV/plain-text file while preserving file order.

    The rank files produced by the BrdU summary utilities contain a read_id
    column. For simpler one-ID-per-line files, the first non-empty field on each
    line is used.
    """
    read_ids = []
    seen = set()

    if limit is not None and limit < 1:
        raise ValueError("read ID limit must be at least 1 when provided.")

    def add_read_id(value):
        value = value.strip()
        if value and value not in seen:
            read_ids.append(value)
            seen.add(value)

    def split_fields(line, delimiter):
        stripped = line.rstrip("\n")
        if delimiter:
            return stripped.split(delimiter)
        return re.split(r"\s+", stripped.strip())

    with open(read_ids_file, "r", encoding="utf-8") as handle:
        first_line = handle.readline()
        if not first_line:
            return read_ids

        if "\t" in first_line:
            delimiter = "\t"
        elif "," in first_line:
            delimiter = ","
        else:
            delimiter = None

        first_fields = split_fields(first_line, delimiter)

        if read_id_column in first_fields:
            read_id_index = first_fields.index(read_id_column)
        else:
            read_id_index = 0
            value = first_fields[read_id_index].strip()
            add_read_id(value)

        for line in handle:
            if limit is not None and len(read_ids) >= limit:
                break
            if not line.strip():
                continue
            fields = split_fields(line, delimiter)
            if read_id_index >= len(fields):
                continue
            add_read_id(fields[read_id_index])

    return read_ids


def get_project_paths():
    """
    Build paths to the workflow-managed BAM directory.

    Sorted BAMs and BAM indexes are stored next to the source BAMs in data/bam
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
    Ensure the BAM file is coordinate-sorted and saved with a standardized name.

    Region-based BAM fetching requires a coordinate-sorted and indexed BAM. This
    function checks the BAM header to see whether the input is already sorted.
    If it is already coordinate-sorted, the input path is used directly. If it
    is not sorted, samtools sort is used to create a sorted copy in data/bam.

    The output name is standardized as:
        <original_basename>.sorted.indexed.bam

    This makes later indexing and reuse predictable.

    Parameters
    ----------
    bam_path : str
        Path to the input BAM file.

    bam_dir : str
        Directory where the sorted/standardized BAM should be stored.

    Returns
    -------
    str
        Path to the standardized coordinate-sorted BAM file.
    """
    bam = pysam.AlignmentFile(bam_path, "rb", check_sq=False)
    header = bam.header.to_dict()
    bam.close()

    sort_order = header.get("HD", {}).get("SO", "unknown")

    os.makedirs(bam_dir, exist_ok=True)

    # If the BAM header says it is already coordinate-sorted, use the user
    # supplied file directly instead of redirecting to a renamed cache file.
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

    standardized_filename = f"{original_base}.sorted.indexed.bam"
    standardized_path = os.path.join(bam_dir, standardized_filename)

    # Reuse an existing standardized BAM to save time on repeated workflow runs.
    if (
        os.path.exists(standardized_path)
        and os.path.getmtime(standardized_path) >= os.path.getmtime(bam_path)
    ):
        print(f"[INFO] Reusing standardized BAM: {standardized_path}", file=sys.stderr)
        return standardized_path

    print(f"[INFO] Sorting BAM → {standardized_path}", file=sys.stderr)
    subprocess.run(["samtools", "sort", "-o", standardized_path, bam_path], check=True)
    return standardized_path


def check_and_index_bam(bam_path, index_dir):
    """
    Ensure the sorted BAM has a current BAM index.

    pysam needs a BAM index to fetch reads from a specific genomic region. This
    function creates the index if it does not exist, or rebuilds it if the BAM
    file is newer than the index.

    The index is written next to the BAM file in data/bam.

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
    bam_basename = os.path.basename(bam_path)
    index_path = os.path.join(index_dir, f"{bam_basename}.bai")

    os.makedirs(index_dir, exist_ok=True)

    rebuild = False

    if not os.path.exists(index_path):
        print("[INFO] No index found. Creating index...", file=sys.stderr)
        rebuild = True
    else:
        # Rebuild the index if the BAM has been modified more recently than the
        # index. This prevents using a stale index with a newer BAM.
        if os.path.getmtime(index_path) < os.path.getmtime(bam_path):
            print("[INFO] Index older than BAM. Rebuilding index...", file=sys.stderr)
            rebuild = True
        else:
            print(f"[INFO] Index found: {index_path}", file=sys.stderr)
            return index_path

    subprocess.run(["samtools", "index", "-o", index_path, bam_path], check=True)
    print(f"[INFO] Indexing complete: {index_path}", file=sys.stderr)
    return index_path


def iter_candidate_reads(bam, chrom=None, start=None, end=None):
    """Yield reads either from a requested region or from the whole BAM."""
    if chrom is None:
        yield from bam.fetch(until_eof=True)
    else:
        yield from bam.fetch(chrom, start, end)


def extract_brdu(bam_path, chrom=None, start=None, end=None, read_id_filter=None, read_ids=None):
    """
    Extract BrdU modification calls from reads overlapping a genomic region.

    This function reads a sorted/indexed BAM file and searches for modified base
    calls with modification code "b", which represents BrdU in this workflow.
    It maps each modified query base back to its reference coordinate and keeps
    only BrdU calls that fall inside the requested region.

    Only thymine bases are kept because BrdU is interpreted as a thymidine analog
    in this pipeline.

    Parameters
    ----------
    bam_path : str
        Path to the sorted and indexed BAM file.

    chrom : str, optional
        Chromosome or contig name to fetch from the BAM.

    start : int, optional
        Start coordinate of the region, using 0-based BED-style coordinates.

    end : int, optional
        End coordinate of the region, using a half-open interval.

    read_id_filter : str, optional
        Optional substring used to restrict extraction to a specific read ID.

    read_ids : list[str], optional
        Optional exact read ID list used to restrict extraction to many reads.

    Returns
    -------
    list[tuple]
        A sorted list of BrdU calls. Each tuple contains:
        chrom, start, end, read_id, base, normalized_probability.
    """
    results = []
    read_id_set = set(read_ids) if read_ids else None
    read_id_order = {read_id: i for i, read_id in enumerate(read_ids or [])}

    _, index_dir = get_project_paths()
    index_path = os.path.join(index_dir, f"{os.path.basename(bam_path)}.bai")

    bam = pysam.AlignmentFile(bam_path, "rb", index_filename=index_path)

    for read in iter_candidate_reads(bam, chrom, start, end):
        if read.is_unmapped:
            continue

        # Allow partial matching for the single-ID mode so users can provide
        # either the full read ID or a unique substring from the read name.
        if read_id_filter and read_id_filter not in read.query_name:
            continue

        if read_id_set is not None and read.query_name not in read_id_set:
            continue

        if not read.modified_bases:
            continue

        for (_, _, mod_code), mod_list in read.modified_bases.items():
            # DNAscent/modBAM BrdU calls are represented with modification code
            # "b". Depending on pysam/version behavior, the code may appear as
            # the character "b" or its ASCII integer value.
            if mod_code != ord("b") and mod_code != "b":
                continue

            # modified_bases reports positions in read/query coordinates.
            # To plot them on the genome, convert query positions to reference
            # positions using the read alignment.
            aligned_pairs = dict(read.get_aligned_pairs(matches_only=True))

            for query_pos, raw_prob in mod_list:
                ref_pos = aligned_pairs.get(query_pos)

                if ref_pos is None:
                    continue

                if start is not None and end is not None and not (start <= ref_pos < end):
                    continue

                base = read.query_sequence[query_pos]

                # BrdU replaces thymidine, so this keeps the output focused on
                # biologically relevant T positions instead of all modified calls.
                if base != "T":
                    continue

                results.append((
                    read.reference_name,
                    ref_pos,
                    ref_pos + 1,
                    read.query_name,
                    base,
                    raw_prob / 255.0
                ))

    bam.close()

    # Sort by read ID and genomic position so each read's BrdU calls are grouped
    # together in a stable order for downstream plotting. When a ranked list is
    # supplied, preserve that ranking in the downstream plot numbering.
    if read_id_order:
        results.sort(key=lambda x: (read_id_order.get(x[3], len(read_id_order)), x[1]))
    else:
        results.sort(key=lambda x: (x[3], x[1]))
    return results


def main():
    """
    Run the BrdU extraction workflow.

    This function validates the input BAM, prepares a sorted and indexed BAM,
    extracts BrdU calls for the requested region, and writes the results either
    to a user-provided output file or to standard output.

    Output columns are:
        chrom, start, end, read_id, base, BrdU_probability
    """
    args = parse_args()

    if args.read_ids_file:
        if any(value is not None for value in (args.chrom, args.start, args.end)):
            if args.chrom is None or args.start is None or args.end is None:
                print(
                    "[ERROR] Provide --chrom, --start, and --end together, or omit all three "
                    "when using --read_ids_file across the whole BAM.",
                    file=sys.stderr,
                )
                sys.exit(1)
    else:
        if args.chrom is None or args.start is None or args.end is None:
            print(
                "[ERROR] --chrom, --start, and --end are required unless --read_ids_file is provided.",
                file=sys.stderr,
            )
            sys.exit(1)

    if args.start is not None and args.end is not None and args.start >= args.end:
        print("[ERROR] --start must be less than --end.", file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(args.bam):
        print(f"[ERROR] BAM not found: {args.bam}", file=sys.stderr)
        sys.exit(1)

    sorted_dir, index_dir = get_project_paths()

    bam_path = check_and_sort_bam(args.bam, sorted_dir)
    check_and_index_bam(bam_path, index_dir)

    read_ids = None
    if args.read_ids_file:
        if not os.path.exists(args.read_ids_file):
            print(f"[ERROR] Read IDs file not found: {args.read_ids_file}", file=sys.stderr)
            sys.exit(1)

        try:
            read_ids = load_read_ids(
                args.read_ids_file,
                args.read_id_column,
                args.read_ids_limit,
            )
        except ValueError as error:
            print(f"[ERROR] {error}", file=sys.stderr)
            sys.exit(1)
        if not read_ids:
            print(f"[ERROR] No read IDs found in {args.read_ids_file}", file=sys.stderr)
            sys.exit(1)
        print(f"[INFO] Loaded {len(read_ids)} read IDs from {args.read_ids_file}", file=sys.stderr)

    results = extract_brdu(
        bam_path,
        args.chrom,
        args.start,
        args.end,
        args.read_id,
        read_ids,
    )

    out = open(args.output, "w") if args.output else sys.stdout

    for r in results:
        # Write BED-like rows:
        # chrom, 0-based start, 1-base end, read ID, base, normalized BrdU score.
        out.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]:.7f}\n")

    if args.output:
        out.close()
        print(f"[INFO] Results written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
