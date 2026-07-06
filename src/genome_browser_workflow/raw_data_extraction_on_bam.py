#!/usr/bin/env python3
"""
Create strand-specific BrdU bedgraph files from a BAM using modkit pileup.

Outputs are written to data/bedgraph by default and contain six tab-delimited
columns:
    chrom, start, end, frac_mod, Nmod, Nvalid_cov
"""

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pysam


def parse_args():
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
        "--bedgraph-dir",
        default=None,
        help="Output directory for generated bedgraph files.",
    )
    return parser.parse_args()


def repo_root():
    return Path(__file__).resolve().parents[2]


def project_paths():
    root = repo_root()
    return {
        "bam_dir": root / "data" / "bam",
        "sorted_bam_dir": root / "data" / "sorted_bam",
        "index_dir": root / "data" / "index_sorted_bam_bai",
        "bedgraph_dir": root / "data" / "bedgraph",
    }


def resolve_existing_path(value, fallback_dir=None, description="file"):
    path = Path(value).expanduser()

    if path.exists():
        return path.resolve()

    if fallback_dir is not None:
        candidate = fallback_dir / value
        if candidate.exists():
            return candidate.resolve()

    raise FileNotFoundError(f"{description} not found: {value}")


def strip_bam_suffix(path):
    name = Path(path).name
    return name[:-4] if name.endswith(".bam") else Path(name).stem


def standardized_bam_name(input_bam):
    base = strip_bam_suffix(input_bam)

    if base.endswith(".sorted.indexed"):
        base = base[: -len(".sorted.indexed")]
    elif base.endswith(".sorted"):
        base = base[: -len(".sorted")]

    return f"{base}.sorted.indexed.bam"


def check_and_sort_bam(input_bam, sorted_bam_dir):
    sorted_bam_dir.mkdir(parents=True, exist_ok=True)
    output_bam = sorted_bam_dir / standardized_bam_name(input_bam)

    if output_bam.exists() and output_bam.stat().st_mtime >= input_bam.stat().st_mtime:
        print(f"[INFO] Reusing sorted BAM: {output_bam}", file=sys.stderr)
        return output_bam

    with pysam.AlignmentFile(str(input_bam), "rb", check_sq=False) as bam:
        sort_order = bam.header.to_dict().get("HD", {}).get("SO", "unknown")

    if sort_order == "coordinate":
        print(f"[INFO] Copying coordinate-sorted BAM to {output_bam}", file=sys.stderr)
        shutil.copy2(input_bam, output_bam)
    else:
        print(f"[INFO] Sorting BAM with samtools: {output_bam}", file=sys.stderr)
        subprocess.run(
            ["samtools", "sort", "-o", str(output_bam), str(input_bam)],
            check=True,
        )

    return output_bam


def check_and_index_bam(sorted_bam, index_dir):
    index_dir.mkdir(parents=True, exist_ok=True)

    adjacent_index = Path(f"{sorted_bam}.bai")
    archived_index = index_dir / f"{sorted_bam.name}.bai"

    needs_index = (
        not adjacent_index.exists()
        or adjacent_index.stat().st_mtime < sorted_bam.stat().st_mtime
    )

    if needs_index:
        print(f"[INFO] Indexing BAM with samtools: {sorted_bam}", file=sys.stderr)
        subprocess.run(["samtools", "index", str(sorted_bam)], check=True)
    else:
        print(f"[INFO] Reusing adjacent BAM index: {adjacent_index}", file=sys.stderr)

    if (
        not archived_index.exists()
        or archived_index.stat().st_mtime < adjacent_index.stat().st_mtime
    ):
        shutil.copy2(adjacent_index, archived_index)
        print(f"[INFO] Archived BAM index: {archived_index}", file=sys.stderr)

    return adjacent_index


def run_modkit_pileup(sorted_bam, reference, output_prefix, bedgraph_dir, threads):
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
        "--only-tabs",
        "--log-filepath",
        str(log_path),
    ]

    print("[INFO] Running modkit pileup.", file=sys.stderr)
    subprocess.run(cmd, check=True)

    write_strand_bedgraphs(bedmethyl, positive, negative)

    return positive, negative, bedmethyl


def write_strand_bedgraphs(bedmethyl, positive_output, negative_output):
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
    check_and_index_bam(sorted_bam, paths["index_dir"])

    positive, negative, bedmethyl = run_modkit_pileup(
        sorted_bam=sorted_bam,
        reference=reference,
        output_prefix=output_prefix,
        bedgraph_dir=bedgraph_dir,
        threads=args.threads,
    )

    print(f"[INFO] Bedmethyl written to {bedmethyl}", file=sys.stderr)
    print(f"[INFO] Positive strand bedgraph written to {positive}", file=sys.stderr)
    print(f"[INFO] Negative strand bedgraph written to {negative}", file=sys.stderr)


if __name__ == "__main__":
    main()
