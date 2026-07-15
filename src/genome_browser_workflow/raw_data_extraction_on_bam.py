#!/usr/bin/env python3
"""
Create strand-specific BrdU bedgraph files from a BAM using modkit pileup.

Outputs are written to data/bedgraph by default and contain six tab-delimited
columns:
    chrom, start, end, frac_mod, Nmod, Nvalid_cov
"""

import argparse
import subprocess
import sys
from pathlib import Path

import pysam


def parse_args():
    """
    Parse command-line arguments for BAM-to-bedGraph generation.

    Collects the input BAM, reference FASTA, output prefix, thread count, and
    optional bedGraph output directory.
    """
    parser = argparse.ArgumentParser(
        description="Extract positive and negative strand BrdU bedgraphs from a BAM."
    )
    parser.add_argument(
        "bam",
        help="BAM filename in data/bam or an explicit path to a BAM file.",
    )
    parser.add_argument(
        "--ref",
        required=True,
        help="Reference FASTA used by modkit pileup.",
    )
    parser.add_argument(
        "-o",
        "--output-prefix",
        default=None,
        help="Output prefix. Defaults to the input BAM basename without .bam.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=12,
        help="Threads passed to modkit pileup.",
    )
    parser.add_argument(
        "--mod-threshold",
        default="B:0.5",
        help="Modification threshold passed to modkit pileup.",
    )
    parser.add_argument(
        "--filter-threshold",
        default="0.0",
        help=(
            "Global filter threshold passed to modkit pileup. Defaults to 0.0 "
            "to prevent modkit from estimating its default filtering threshold."
        ),
    )
    parser.add_argument(
        "--bedgraph-dir",
        default=None,
        help="Output directory for generated bedgraph files.",
    )
    return parser.parse_args()


def repo_root():
    """
    Return the repository root directory.

    Assumes this script is located two directories below the project root.
    """
    return Path(__file__).resolve().parents[2]


def project_paths():
    """
    Return standard project directories used by the workflow.

    Provides paths for input BAMs, sorted BAMs, BAM indexes, and generated
    bedGraph files.
    """
    root = repo_root()
    return {
        "bam_dir": root / "data" / "bam",
        "sorted_bam_dir": root / "data" / "sorted_bam",
        "bedgraph_dir": root / "data" / "bedgraph",
    }


def resolve_existing_path(value, fallback_dir=None, description="file"):
    """
    Resolve an input path and confirm that it exists.

    Checks the provided path first, then optionally checks inside a fallback
    directory. Raises FileNotFoundError if neither path exists.
    """
    path = Path(value).expanduser()

    if path.exists():
        return path.resolve()

    if fallback_dir is not None:
        candidate = fallback_dir / value
        if candidate.exists():
            return candidate.resolve()

    raise FileNotFoundError(f"{description} not found: {value}")


def strip_bam_suffix(path):
    """
    Return a BAM filename stem without the .bam suffix.

    If the name does not end in .bam, returns the normal Path stem.
    """
    name = Path(path).name
    return name[:-4] if name.endswith(".bam") else Path(name).stem


def prepared_bam_name(input_bam):
    """
    Return the canonical prepared BAM name for this workflow.

    The prepared BAM is always named as sorted, indexed, and BrdU-detected so it
    is clear which file should be passed to modkit pileup.
    """
    base = strip_bam_suffix(input_bam)
    known_suffixes = (
        ".sorted.indexed.BrdU.detect",
        ".BrdU.detect.sorted.indexed",
        ".sorted.indexed",
        ".BrdU.detect",
        ".sorted",
    )

    changed = True
    while changed:
        changed = False
        for suffix in known_suffixes:
            if base.endswith(suffix):
                base = base[: -len(suffix)]
                changed = True
                break

    return f"{base}.sorted.indexed.BrdU.detect.bam"


def check_and_sort_bam(input_bam, sorted_bam_dir):
    """
    Create or reuse the workflow's prepared coordinate-sorted BAM.

    The BAM header's sort-order flag is not trusted here because an incorrectly
    marked BAM will still fail indexing. A prepared BAM is reused only when both
    the BAM and adjacent index are newer than the source BAM.
    """
    sorted_bam_dir.mkdir(parents=True, exist_ok=True)
    output_bam = sorted_bam_dir / prepared_bam_name(input_bam)
    output_index = Path(f"{output_bam}.bai")

    if (
        output_bam.exists()
        and output_index.exists()
        and output_bam.stat().st_mtime >= input_bam.stat().st_mtime
        and output_index.stat().st_mtime >= output_bam.stat().st_mtime
    ):
        print(f"[INFO] Reusing prepared BAM: {output_bam}", file=sys.stderr)
        return output_bam

    print(f"[INFO] Sorting BAM with samtools: {output_bam}", file=sys.stderr)
    subprocess.run(
        ["samtools", "sort", "-o", str(output_bam), str(input_bam)],
        check=True,
    )

    return output_bam


def check_and_index_bam(sorted_bam):
    """
    Ensure a sorted BAM has an up-to-date index.

    Creates or reuses the adjacent .bai index and returns its path.
    """
    adjacent_index = Path(f"{sorted_bam}.bai")

    needs_index = (
        not adjacent_index.exists()
        or adjacent_index.stat().st_mtime < sorted_bam.stat().st_mtime
    )

    if needs_index:
        print(f"[INFO] Indexing BAM with samtools: {sorted_bam}", file=sys.stderr)
        subprocess.run(["samtools", "index", str(sorted_bam)], check=True)
    else:
        print(f"[INFO] Reusing adjacent BAM index: {adjacent_index}", file=sys.stderr)

    return adjacent_index


def normalize_mod_code(mod_code):
    """
    Return a readable modified-base code from pysam's modified_bases key.
    """
    if isinstance(mod_code, int):
        try:
            return chr(mod_code)
        except ValueError:
            return str(mod_code)

    return str(mod_code)


def summarize_modified_base_tags(bam_path, sample_limit=50000):
    """
    Inspect a BAM sample for modified-base tags expected by modkit.

    modkit reads modified-base calls from MM/ML SAM tags. This workflow then
    keeps BrdU calls with modification code "b".
    """
    summary = {
        "reads_checked": 0,
        "reads_with_mm": 0,
        "reads_with_ml": 0,
        "reads_with_mod_calls": 0,
        "reads_with_brdu": 0,
        "mod_codes": {},
        "parse_errors": 0,
    }

    with pysam.AlignmentFile(str(bam_path), "rb", check_sq=False) as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            summary["reads_checked"] += 1

            if read.has_tag("MM"):
                summary["reads_with_mm"] += 1
            if read.has_tag("ML"):
                summary["reads_with_ml"] += 1

            try:
                modified_bases = read.modified_bases or {}
            except ValueError:
                summary["parse_errors"] += 1
                modified_bases = {}

            if modified_bases:
                summary["reads_with_mod_calls"] += 1

            read_has_brdu = False
            for (_, _, mod_code), mod_list in modified_bases.items():
                code = normalize_mod_code(mod_code)
                summary["mod_codes"][code] = summary["mod_codes"].get(code, 0) + len(
                    mod_list
                )
                if code == "b":
                    read_has_brdu = True

            if read_has_brdu:
                summary["reads_with_brdu"] += 1

            if summary["reads_checked"] >= sample_limit:
                break

    return summary


def validate_brdu_mod_tags(bam_path):
    """
    Fail early when the BAM does not appear to contain BrdU modBAM calls.
    """
    summary = summarize_modified_base_tags(bam_path)
    codes = ", ".join(
        f"{code}:{count}" for code, count in sorted(summary["mod_codes"].items())
    )
    codes = codes or "none"

    print(
        "[INFO] BAM modified-base preflight: "
        f"checked {summary['reads_checked']} primary mapped reads; "
        f"MM tags in {summary['reads_with_mm']}; "
        f"ML tags in {summary['reads_with_ml']}; "
        f"modified-base calls in {summary['reads_with_mod_calls']}; "
        f"BrdU code b in {summary['reads_with_brdu']}; "
        f"codes observed: {codes}",
        file=sys.stderr,
    )

    if summary["parse_errors"]:
        print(
            "[WARN] Some reads had modified-base tags that pysam could not parse: "
            f"{summary['parse_errors']}",
            file=sys.stderr,
        )

    if summary["reads_checked"] == 0:
        raise ValueError("No primary mapped reads were found in the BAM.")

    if summary["reads_with_mm"] == 0 or summary["reads_with_ml"] == 0:
        raise ValueError(
            "No MM/ML modified-base tags were found in the sampled reads. "
            "modkit pileup requires a modBAM with modified-base calls; a normal "
            "aligned BAM can be large and still produce 0 modkit rows."
        )

    if summary["reads_with_brdu"] == 0:
        raise ValueError(
            'No BrdU modification code "b" was found in the sampled reads. '
            "This genome browser workflow only plots BrdU calls encoded as "
            '"b"; run it on the DNAscent/Dorado BrdU-detected modBAM or check '
            "which modification codes are present."
        )


def run_modkit_pileup(
    sorted_bam,
    reference,
    output_prefix,
    bedgraph_dir,
    threads,
    mod_threshold,
    filter_threshold,
):
    """
    Run modkit pileup and create strand-specific bedGraphs.

    Generates a full bedMethyl file, writes modkit logs, splits BrdU calls into
    positive and negative strand bedGraphs, and returns the output paths.
    """
    bedgraph_dir.mkdir(parents=True, exist_ok=True)

    bedmethyl = bedgraph_dir / f"{output_prefix}.full.bedmethyl"
    positive = bedgraph_dir / f"{output_prefix}.positive.bedgraph"
    negative = bedgraph_dir / f"{output_prefix}.negative.bedgraph"
    log_path = bedgraph_dir / f"{output_prefix}.modkit.log"

    for generated_path in (bedmethyl, positive, negative, log_path):
        if generated_path.exists():
            generated_path.unlink()

    cmd = [
        "modkit",
        "pileup",
        str(sorted_bam),
        str(bedmethyl),
        "--ref",
        str(reference),
        "--threads",
        str(threads),
        "--mod-threshold",
        mod_threshold,
        "--filter-threshold",
        filter_threshold,
        "--only-tabs",
        "--log-filepath",
        str(log_path),
    ]

    print("[INFO] Running modkit pileup.", file=sys.stderr)
    subprocess.run(cmd, check=True)

    write_strand_bedgraphs(bedmethyl, positive, negative)

    return positive, negative, bedmethyl


def write_strand_bedgraphs(bedmethyl, positive_output, negative_output):
    """
    Split a bedMethyl file into positive and negative BrdU bedGraphs.

    Keeps BrdU modification rows on + or - strands, calculates frac_mod as
    Nmod / valid coverage, and writes chrom, start, end, frac_mod, Nmod, and
    valid coverage columns.
    """
    with open(bedmethyl, "r") as source, open(positive_output, "w") as pos, open(
        negative_output, "w"
    ) as neg:
        for line in source:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 12:
                continue

            mod_code = fields[3]
            strand = fields[5]

            if mod_code != "b" or strand not in {"+", "-"}:
                continue

            try:
                valid_cov = float(fields[9])
                n_mod = float(fields[11])
            except ValueError:
                continue

            frac_mod = n_mod / valid_cov if valid_cov > 0 else 0.0
            output_line = (
                f"{fields[0]}\t{fields[1]}\t{fields[2]}\t"
                f"{frac_mod:.8g}\t{n_mod:.0f}\t{valid_cov:.0f}\n"
            )

            if strand == "+":
                pos.write(output_line)
            else:
                neg.write(output_line)


def main():
    """
    Run the BAM-to-bedGraph workflow from command-line arguments.

    Resolves inputs, sorts and indexes the BAM, runs modkit pileup, writes
    strand-specific bedGraphs, and reports the generated output paths.
    """
    args = parse_args()
    paths = project_paths()

    try:
        input_bam = resolve_existing_path(args.bam, paths["bam_dir"], "BAM")
        reference = resolve_existing_path(args.ref, None, "reference FASTA")
    except FileNotFoundError as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        sys.exit(1)

    bedgraph_dir = (
        Path(args.bedgraph_dir).expanduser().resolve()
        if args.bedgraph_dir
        else paths["bedgraph_dir"]
    )

    output_prefix = args.output_prefix or strip_bam_suffix(input_bam)

    sorted_bam = check_and_sort_bam(input_bam, paths["sorted_bam_dir"])
    check_and_index_bam(sorted_bam)

    try:
        validate_brdu_mod_tags(sorted_bam)
    except ValueError as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        sys.exit(1)

    positive, negative, bedmethyl = run_modkit_pileup(
        sorted_bam=sorted_bam,
        reference=reference,
        output_prefix=output_prefix,
        bedgraph_dir=bedgraph_dir,
        threads=args.threads,
        mod_threshold=args.mod_threshold,
        filter_threshold=args.filter_threshold,
    )

    print(f"[INFO] Bedmethyl written to {bedmethyl}", file=sys.stderr)
    print(f"[INFO] Positive strand bedgraph written to {positive}", file=sys.stderr)
    print(f"[INFO] Negative strand bedgraph written to {negative}", file=sys.stderr)


if __name__ == "__main__":
    main()
