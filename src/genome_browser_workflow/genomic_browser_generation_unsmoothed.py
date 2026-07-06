#!/usr/bin/env python3
"""
Generate unsmoothed per-chromosome BrdU genome browser plots.
"""

from genomic_browser_generation import add_common_args, generate_genome_browser

import argparse
from pathlib import Path


def repo_root():
    return Path(__file__).resolve().parents[2]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate unsmoothed genome browser plots from BrdU bedgraphs."
    )
    add_common_args(parser)
    parser.set_defaults(
        output_dir=str(
            repo_root() / "results" / "genome_browser_results" / "unsmoothed"
        )
    )
    return parser.parse_args()


def main():
    args = parse_args()
    generate_genome_browser(
        positive_bedgraph=args.positive_bedgraph,
        negative_bedgraph=args.negative_bedgraph,
        output_dir=args.output_dir,
        prefix=args.prefix,
        g4_bed=args.g4_bed,
        trna_bed=args.trna_bed,
        te_bed=args.te_bed,
        smoothed=False,
    )


if __name__ == "__main__":
    main()
