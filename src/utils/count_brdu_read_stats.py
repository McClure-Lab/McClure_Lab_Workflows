#!/usr/bin/env python3
"""
Count BrdU-positive reads and BrdU calls from BAM modified-base tags.
"""

import argparse
import math
import sys

import pysam


def parse_args():
    parser = argparse.ArgumentParser(
        description="Count BrdU-positive primary mapped reads in a BAM."
    )
    parser.add_argument("bam", help="Input BAM path.")
    parser.add_argument(
        "--mod-threshold",
        type=float,
        required=True,
        help="Minimum BrdU modification probability from 0 to 1.",
    )
    parser.add_argument(
        "--mod-code",
        default="b",
        help='Modified-base code to count. Defaults to "b" for BrdU.',
    )
    return parser.parse_args()


def normalize_mod_code(mod_code):
    if isinstance(mod_code, int):
        try:
            return chr(mod_code)
        except ValueError:
            return str(mod_code)

    return str(mod_code)


def count_brdu_stats(bam_path, mod_threshold, mod_code):
    """
    Return counts for primary mapped reads with passing BrdU calls.

    pysam reports modified-base probabilities as integer ML values from 0 to
    255. A call passes when ML / 255 is greater than or equal to the requested
    probability threshold.
    """
    min_ml = math.ceil(mod_threshold * 255)
    total_primary_mapped = 0
    reads_with_mm = 0
    reads_with_ml = 0
    reads_with_mod_calls = 0
    brdu_positive_reads = 0
    brdu_calls = 0
    parse_errors = 0

    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            total_primary_mapped += 1

            if read.has_tag("MM"):
                reads_with_mm += 1
            if read.has_tag("ML"):
                reads_with_ml += 1

            try:
                modified_bases = read.modified_bases or {}
            except ValueError:
                parse_errors += 1
                continue

            if modified_bases:
                reads_with_mod_calls += 1

            read_brdu_calls = 0
            for (_, _, observed_mod_code), mod_list in modified_bases.items():
                if normalize_mod_code(observed_mod_code) != mod_code:
                    continue

                for _, ml_value in mod_list:
                    if ml_value >= min_ml:
                        read_brdu_calls += 1

            if read_brdu_calls > 0:
                brdu_positive_reads += 1
                brdu_calls += read_brdu_calls

    return {
        "total_primary_mapped_from_pysam": total_primary_mapped,
        "reads_with_mm": reads_with_mm,
        "reads_with_ml": reads_with_ml,
        "reads_with_mod_calls": reads_with_mod_calls,
        "reads_with_brdu_calls": brdu_positive_reads,
        "brdu_calls": brdu_calls,
        "parse_errors": parse_errors,
        "min_ml": min_ml,
    }


def main():
    args = parse_args()

    if args.mod_threshold < 0 or args.mod_threshold > 1:
        print("[ERROR] --mod-threshold must be between 0 and 1.", file=sys.stderr)
        return 1

    stats = count_brdu_stats(args.bam, args.mod_threshold, args.mod_code)

    for key, value in stats.items():
        print(f"{key}\t{value}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
