#!/usr/bin/env python3
"""
Rank reads by their average BrdU modification probability.

Purpose
-------
This script reads a tab-separated file produced by:

    modkit extract full

It identifies rows associated with a selected modification code, calculates the
average modification probability for every read, and writes a ranked TSV with
the highest average probability first.

The default modification code is:

    b

which is used here to represent BrdU calls.

Basic usage
-----------
Run the script with an input Modkit TSV:

    python rank_reads_by_average_brdu.py modkit_extract_full.tsv

Compressed input is also supported automatically by pandas:

    python rank_reads_by_average_brdu.py modkit_extract_full.tsv.gz

Specify a custom output filename:

    python rank_reads_by_average_brdu.py \
        modkit_extract_full.tsv \
        --output ranked_reads.tsv

Specify a different modification code:

    python rank_reads_by_average_brdu.py \
        modkit_extract_full.tsv \
        --mod-code b

Change the number of input rows processed at once:

    python rank_reads_by_average_brdu.py \
        modkit_extract_full.tsv \
        --chunksize 500000

Default output filename
-----------------------
When --output is not provided, the input prefix is followed by:

    .average_brdu_probability.sorted.tsv

Examples:

    sample.tsv
        becomes:
    sample.average_brdu_probability.sorted.tsv

    sample.tsv.gz
        becomes:
    sample.average_brdu_probability.sorted.tsv

Chunked processing
------------------
Modkit extract files can be extremely large. The script therefore reads the
input in chunks rather than loading the entire TSV into memory.

The default chunk size is:

    1,000,000 rows

Only four columns are loaded:

    read_id
    read_length
    mod_qual
    mod_code

This reduces memory use compared with loading all Modkit columns.

BrdU probability calculation
----------------------------
For every read, the script accumulates:

    probability_sum
        Sum of all mod_qual values for rows matching the requested mod_code.

    evaluated_positions
        Number of matching modification rows for the read.

    read_length
        Read length reported by Modkit.

The average probability is calculated as:

    average_brdu_probability =
        probability_sum / evaluated_positions

Important interpretation
------------------------
This script ranks reads by the average probability across all rows having the
selected modification code.

It does not:

    - Apply a minimum modification-probability threshold.
    - Count only positions classified as positive.
    - Rank reads by the total number of BrdU calls.
    - Rank reads primarily by read length.
    - Normalize the probability sum by read length.

A read with a high average probability but relatively few evaluated positions
can rank above a read with many evaluated positions but a lower average.

Ranking
-------
Reads are sorted by:

    1. Highest average_brdu_probability.
    2. Highest evaluated_brdu_positions as a tie-breaker.

A stable mergesort is used.

Output columns
--------------
The output TSV contains:

    rank
        One-based ranking after sorting.

    read_id
        Read identifier from the Modkit TSV.

    read_length
        Read length reported by Modkit.

    evaluated_brdu_positions
        Number of rows matching the selected modification code.

    brdu_probability_sum
        Sum of matching mod_qual values.

    average_brdu_probability
        Mean matching mod_qual value for the read.

Numeric output is written with six digits after the decimal point.
"""

import argparse
from pathlib import Path

import pandas as pd


# =============================================================================
# modkit extract full columns
# =============================================================================
#
# This list documents the expected column order from a full Modkit extraction.
#
# The current analysis reads the input using the header row and only loads:
#
#     read_id
#     read_length
#     mod_qual
#     mod_code
#
# Therefore, COLUMN_NAMES is retained as documentation of the full schema but
# is not passed directly to pandas.
#
COLUMN_NAMES = [
    "read_id",
    "forward_read_position",
    "ref_position",
    "chrom",
    "mod_strand",
    "ref_strand",
    "ref_mod_strand",
    "fw_soft_clipped_start",
    "fw_soft_clipped_end",
    "alignment_start",
    "alignment_end",
    "read_length",
    "mod_qual",
    "mod_code",
    "base_qual",
    "ref_kmer",
    "query_kmer",
    "canonical_base",
    "modified_primary_base",
    "inferred",
    "flag",
]


# =============================================================================
# Command-line argument parsing
# =============================================================================


def parse_arguments():
    """
    Define and parse command-line arguments.

    Positional arguments
    --------------------
    input_tsv
        TSV or compressed TSV created by ``modkit extract full``.

    Optional arguments
    ------------------
    -o, --output
        Explicit output TSV path.

    --mod-code
        Modification code retained for analysis. The default is ``b``.

    --chunksize
        Number of input rows read by pandas at one time.

    Returns
    -------
    argparse.Namespace
        Parsed command-line settings.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Calculate the average BrdU probability for each read from a "
            "modkit extract full TSV and sort reads from highest to lowest."
        )
    )

    parser.add_argument(
        "input_tsv",
        help="TSV file produced by modkit extract full.",
    )

    parser.add_argument(
        "-o",
        "--output",
        help=(
            "Output TSV filename. By default, the input prefix is followed "
            "by .average_brdu_probability.sorted.tsv"
        ),
    )

    parser.add_argument(
        "--mod-code",
        default="b",
        help="Modification code representing BrdU. Default: b",
    )

    parser.add_argument(
        "--chunksize",
        type=int,
        default=1_000_000,
        help="Number of TSV rows processed at once. Default: 1000000",
    )

    return parser.parse_args()


# =============================================================================
# Output-path construction
# =============================================================================


def determine_output_path(
    input_path,
    output_argument,
):
    """
    Determine the output TSV path.

    When an explicit --output value is supplied, that value is used directly.

    Otherwise, the function removes a recognized input extension and appends:

        .average_brdu_probability.sorted.tsv

    Recognized input suffixes
    -------------------------
    .tsv.gz
    .tsv

    Parameters
    ----------
    input_path
        Input path represented as a pathlib.Path.

    output_argument
        Optional output path supplied through --output.

    Returns
    -------
    pathlib.Path
        Final output path.
    """
    if output_argument:
        return Path(output_argument)

    input_name = input_path.name

    # Remove the complete compressed TSV suffix.
    if input_name.endswith(".tsv.gz"):
        prefix = input_name[:-7]

    # Remove a normal TSV suffix.
    elif input_name.endswith(".tsv"):
        prefix = input_name[:-4]

    # Fall back to pathlib's final-suffix removal.
    else:
        prefix = input_path.stem

    # By default, write the output beside the input file.
    return input_path.with_name(
        f"{prefix}.average_brdu_probability.sorted.tsv"
    )


# =============================================================================
# Main analysis
# =============================================================================


def main():
    """
    Read the Modkit TSV, aggregate statistics by read, rank reads, and write TSV.

    Processing stages
    -----------------
    1. Parse command-line arguments.
    2. Validate the input file.
    3. Determine the output filename.
    4. Read selected columns in chunks.
    5. Retain rows matching the requested modification code.
    6. Calculate per-chunk read summaries.
    7. Accumulate those summaries across all chunks.
    8. Calculate each read's average modification probability.
    9. Sort and rank the reads.
    10. Write the final TSV and print the top ten reads.
    """
    args = parse_arguments()

    input_path = Path(
        args.input_tsv
    )

    # Stop immediately when the requested input does not exist.
    if not input_path.is_file():
        raise FileNotFoundError(
            f"Input TSV was not found: {input_path}"
        )

    output_path = determine_output_path(
        input_path,
        args.output,
    )

    print(f"Input TSV: {input_path}")
    print(f"BrdU modification code: {args.mod_code}")
    print(f"Chunk size: {args.chunksize:,} rows")
    print("Processing TSV...")

    # -------------------------------------------------------------------------
    # Per-read cumulative statistics
    # -------------------------------------------------------------------------
    #
    # The dictionary is keyed by read ID.
    #
    # Value format:
    #
    # {
    #     read_id: {
    #         "probability_sum": float,
    #         "evaluated_positions": int,
    #         "read_length": int or float
    #     }
    # }
    #
    # Statistics are accumulated across chunks because rows for one read may
    # appear in more than one input chunk.
    #
    read_statistics = {}

    # Number of all input rows read, regardless of modification code.
    rows_processed = 0

    # Number of valid retained rows matching args.mod_code.
    brdu_rows_processed = 0

    # -------------------------------------------------------------------------
    # Columns loaded from the input
    # -------------------------------------------------------------------------
    #
    # Loading only required columns substantially reduces memory use.
    #
    columns_to_use = [
        "read_id",
        "read_length",
        "mod_qual",
        "mod_code",
    ]

    # -------------------------------------------------------------------------
    # Create the chunked pandas reader
    # -------------------------------------------------------------------------
    #
    # header=0:
    #     Use the first row as column names.
    #
    # usecols:
    #     Load only the columns required by this analysis.
    #
    # chunksize:
    #     Return an iterator of DataFrames rather than one complete DataFrame.
    #
    # dtype:
    #     Preserve read IDs and modification codes as strings.
    #
    # compression="infer":
    #     Automatically recognize compression from extensions such as .gz.
    #
    chunks = pd.read_csv(
        input_path,
        sep="\t",
        header=0,
        usecols=columns_to_use,
        chunksize=args.chunksize,
        dtype={
            "read_id": "string",
            "mod_code": "string",
        },
        compression="infer",
        low_memory=False,
    )

    # -------------------------------------------------------------------------
    # Process each input chunk
    # -------------------------------------------------------------------------

    for chunk_number, chunk in enumerate(
        chunks,
        start=1,
    ):
        # Count every row loaded from the TSV.
        rows_processed += len(chunk)

        # Keep only records for the requested modification code.
        #
        # With the default settings, only rows where mod_code == "b" remain.
        chunk = chunk.loc[
            chunk["mod_code"] == args.mod_code
        ].copy()

        # A chunk may contain no matching modification records.
        if chunk.empty:
            print(
                f"Chunk {chunk_number}: "
                f"{rows_processed:,} total rows processed; "
                "no BrdU rows in this chunk."
            )
            continue

        # Convert the modification-quality or probability values to numeric.
        #
        # Invalid values become NaN and are removed below.
        chunk["mod_qual"] = pd.to_numeric(
            chunk["mod_qual"],
            errors="coerce",
        )

        # Convert read length to numeric.
        #
        # A missing read length does not cause the row to be discarded because
        # read length is not required for the probability calculation.
        chunk["read_length"] = pd.to_numeric(
            chunk["read_length"],
            errors="coerce",
        )

        # Remove rows that cannot contribute to the analysis.
        #
        # read_id is needed for grouping.
        # mod_qual is needed for the probability calculation.
        chunk = chunk.dropna(
            subset=[
                "read_id",
                "mod_qual",
            ]
        )

        # Count valid matching rows retained after numeric conversion.
        brdu_rows_processed += len(chunk)

        # ---------------------------------------------------------------------
        # Summarize the current chunk by read ID
        # ---------------------------------------------------------------------
        #
        # probability_sum:
        #     Sum of mod_qual for this read in this chunk.
        #
        # evaluated_positions:
        #     Number of matching modification rows for this read in the chunk.
        #
        # read_length:
        #     First read length encountered for the read in the chunk.
        #
        chunk_summary = (
            chunk.groupby(
                "read_id",
                sort=False,
                observed=True,
            )
            .agg(
                probability_sum=(
                    "mod_qual",
                    "sum",
                ),
                evaluated_positions=(
                    "mod_qual",
                    "size",
                ),
                read_length=(
                    "read_length",
                    "first",
                ),
            )
        )

        # ---------------------------------------------------------------------
        # Merge the current chunk summary into the cumulative dictionary
        # ---------------------------------------------------------------------

        for read_id, row in chunk_summary.iterrows():
            # Initialize a read when it has not appeared in an earlier chunk.
            if read_id not in read_statistics:
                read_statistics[read_id] = {
                    "probability_sum": float(
                        row["probability_sum"]
                    ),
                    "evaluated_positions": int(
                        row["evaluated_positions"]
                    ),
                    "read_length": row["read_length"],
                }

            # Add chunk-level values when the read was observed previously.
            else:
                read_statistics[read_id][
                    "probability_sum"
                ] += float(
                    row["probability_sum"]
                )

                read_statistics[read_id][
                    "evaluated_positions"
                ] += int(
                    row["evaluated_positions"]
                )

                # Update read length only when the current value is valid.
                if pd.notna(row["read_length"]):
                    read_statistics[read_id][
                        "read_length"
                    ] = row["read_length"]

        # Print cumulative progress after every chunk.
        print(
            f"Chunk {chunk_number}: "
            f"{rows_processed:,} total rows processed; "
            f"{brdu_rows_processed:,} BrdU rows retained; "
            f"{len(read_statistics):,} reads observed."
        )

    # -------------------------------------------------------------------------
    # Require at least one matching modification row
    # -------------------------------------------------------------------------

    if not read_statistics:
        raise RuntimeError(
            f"No rows were found with mod_code '{args.mod_code}'. "
            "Check the mod_code column in the input TSV."
        )

    # -------------------------------------------------------------------------
    # Calculate final per-read values
    # -------------------------------------------------------------------------

    results = []

    for read_id, statistics in read_statistics.items():
        evaluated_positions = statistics[
            "evaluated_positions"
        ]

        # Calculate the arithmetic mean of mod_qual for the read.
        average_probability = (
            statistics["probability_sum"]
            / evaluated_positions
        )

        read_length = statistics[
            "read_length"
        ]

        # Preserve valid read lengths as integers.
        if pd.notna(read_length):
            read_length = int(
                read_length
            )
        else:
            read_length = pd.NA

        results.append(
            {
                "read_id": read_id,
                "read_length": read_length,
                "evaluated_brdu_positions": (
                    evaluated_positions
                ),
                "brdu_probability_sum": statistics[
                    "probability_sum"
                ],
                "average_brdu_probability": (
                    average_probability
                ),
            }
        )

    # Convert the list of per-read dictionaries into a DataFrame.
    ranked_reads = pd.DataFrame(
        results
    )

    # -------------------------------------------------------------------------
    # Rank reads
    # -------------------------------------------------------------------------
    #
    # Highest average BrdU probability appears first.
    #
    # evaluated_brdu_positions is only used as a tie-breaker. It does not
    # determine the primary ranking.
    #
    # mergesort is stable, which gives deterministic ordering when sort values
    # are otherwise equal.
    #
    ranked_reads = ranked_reads.sort_values(
        by=[
            "average_brdu_probability",
            "evaluated_brdu_positions",
        ],
        ascending=[
            False,
            False,
        ],
        kind="mergesort",
    ).reset_index(
        drop=True
    )

    # Insert a one-based rank as the first output column.
    ranked_reads.insert(
        0,
        "rank",
        range(
            1,
            len(ranked_reads) + 1,
        ),
    )

    # -------------------------------------------------------------------------
    # Write the final ranked TSV
    # -------------------------------------------------------------------------
    #
    # float_format applies six decimal places to floating-point output values.
    #
    ranked_reads.to_csv(
        output_path,
        sep="\t",
        index=False,
        float_format="%.6f",
    )

    # -------------------------------------------------------------------------
    # Print final analysis summary
    # -------------------------------------------------------------------------

    print()
    print("Analysis complete.")
    print(f"Total TSV rows processed: {rows_processed:,}")
    print(f"BrdU rows processed: {brdu_rows_processed:,}")
    print(f"Reads ranked: {len(ranked_reads):,}")
    print(f"Output TSV: {output_path}")

    # -------------------------------------------------------------------------
    # Print the first ten ranked reads
    # -------------------------------------------------------------------------

    print()
    print("Top 10 reads:")

    print(
        ranked_reads[
            [
                "rank",
                "read_id",
                "read_length",
                "evaluated_brdu_positions",
                "average_brdu_probability",
            ]
        ]
        .head(10)
        .to_string(
            index=False
        )
    )


# =============================================================================
# Script entry point
# =============================================================================
#
# Calling main only when this file is executed directly allows the functions to
# be imported by another Python module without automatically starting analysis.
#
if __name__ == "__main__":
    main()