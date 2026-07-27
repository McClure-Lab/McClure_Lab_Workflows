#!/bin/bash

# ==============================================================================
# TOP-20 POINTS-OF-INTEREST COMPARISON WORKFLOW
# ==============================================================================
#
# Purpose
# -------
# This script compares two CSV files containing genomic points of interest
# (POIs). It identifies POIs that occur on the same chromosome and whose
# representative genomic positions are within a configurable distance.
#
# The workflow is intended for top-POI result files stored under:
#
#   /beevol/home/pineirok/workflows/data/csv
#
# Interactive execution
# ---------------------
#
# Run the script from the workflow root with:
#
#   bash src/utils/compare_top_20_poi.sh
#
# During interactive execution, the script:
#
#   1. Lists all CSV files under data/csv.
#   2. Prompts for the first POI CSV.
#   3. Prompts for the second POI CSV.
#   4. Prevents the same file from being selected twice.
#   5. Constructs descriptive output filenames.
#   6. Submits the comparison to SLURM.
#
# Internal SLURM execution
# ------------------------
#
# The submitted job calls this same script using:
#
#   --run-comparison FIRST_CSV SECOND_CSV OUTPUT_PREFIX
#
# Users normally do not need to invoke --run-comparison manually.
#
# POI overlap definition
# ----------------------
#
# Two POIs are considered candidate overlaps when:
#
#   - They are located on the same chromosome.
#   - Their representative positions differ by no more than
#     POSITION_TOLERANCE_BP.
#
# The default tolerance is:
#
#   +/- 1,000 bp
#
# Representative genomic position
# --------------------------------
#
# When a CSV contains a position column, that column is used.
#
# When position is absent, the script calculates the midpoint:
#
#   position = (start + end) // 2
#
# One-to-one matching
# -------------------
#
# A single POI cannot be matched to multiple POIs from the other file.
#
# The workflow:
#
#   1. Generates all same-chromosome candidate pairs within the tolerance.
#   2. Sorts candidate pairs by:
#        - Smallest distance.
#        - First-file rank.
#        - Second-file rank.
#   3. Selects the closest available pair.
#   4. Marks both POIs as used.
#   5. Skips later candidates containing either used POI.
#
# This prevents one POI from being counted as several independent overlaps.
#
# Required CSV columns
# --------------------
#
# Each input CSV must contain:
#
#   chrom
#   start
#   end
#
# Optional useful columns:
#
#   position
#   rank
#   frac_mod
#   BrdU_pct
#
# If rank is absent, ranks are generated from the current row order beginning
# with rank 1.
#
# Invalid-row handling
# --------------------
#
# Rows are removed when:
#
#   - start, end, or position cannot be converted to numeric values.
#   - start is negative.
#   - end is less than or equal to start.
#   - position is negative.
#
# Output files
# ------------
#
# Given input prefixes FIRST and SECOND, the output prefix is:
#
#   FIRST__vs__SECOND
#
# The following CSV files are created:
#
#   FIRST__vs__SECOND_overlaps.csv
#       Contains the selected one-to-one overlapping POI pairs.
#
#   FIRST__vs__SECOND_first_unmatched.csv
#       Contains valid POIs from the first CSV that were not matched.
#
#   FIRST__vs__SECOND_second_unmatched.csv
#       Contains valid POIs from the second CSV that were not matched.
#
#   FIRST__vs__SECOND_summary.csv
#       Contains input counts, overlap counts, unmatched counts, overlap
#       percentages, and the selected tolerance.
#
# Combined SLURM log:
#
#   FIRST.SECOND.poi.log
#
# Python environment
# ------------------
#
# The workflow creates or reuses:
#
#   /beevol/home/pineirok/workflows/.venv
#
# NumPy and pandas are installed only when they cannot already be imported.
#
# ==============================================================================


# ==============================================================================
# SLURM RESOURCE REQUESTS
# ==============================================================================
#
# Resources:
#   - 1 task
#   - 4 CPU cores
#   - 8 GB RAM
#   - 1-hour runtime limit
#
#SBATCH --job-name=compare_top20_poi
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00


# ==============================================================================
# BASH SAFETY SETTINGS
# ==============================================================================
#
# -e:
#   Exit when an unhandled command fails.
#
# -u:
#   Treat unset variables as errors.
#
# -o pipefail:
#   Treat a pipeline as failed when any command in that pipeline fails.
#
set -euo pipefail


# ==============================================================================
# PROJECT SETTINGS
# ==============================================================================

# Root of the workflow project.
PROJECT_DIR="/beevol/home/pineirok/workflows"

# Directory containing input POI CSV files and generated comparison outputs.
CSV_DIR="${PROJECT_DIR}/data/csv"

# Shared Python virtual environment used by this utility.
VENV_DIR="${PROJECT_DIR}/.venv"

# POIs must occur on the same chromosome and their representative positions
# must be no more than this many bases apart.
POSITION_TOLERANCE_BP=1000

# Create the CSV directory when it does not already exist.
mkdir -p "${CSV_DIR}"


# ==============================================================================
# FUNCTION: list_csv_files
# ==============================================================================
#
# Purpose:
#   List CSV files located directly under CSV_DIR.
#
# Behavior:
#   - Does not search subdirectories.
#   - Matches .csv case-insensitively.
#   - Prints filenames rather than full paths.
#   - Sorts the filenames alphabetically.
#
list_csv_files() {
    find "${CSV_DIR}" \
        -maxdepth 1 \
        -type f \
        -iname "*.csv" \
        -printf "%f\n" |
        sort
}


# ==============================================================================
# FUNCTION: resolve_selected_file
# ==============================================================================
#
# Purpose:
#   Resolve a numeric menu selection, filename, or complete file path.
#
# Arguments:
#   $1
#       User selection.
#
#   $2
#       Name of the Bash array containing the available CSV filenames.
#
# Resolution order:
#   1. Numeric input is interpreted as a one-based menu selection.
#   2. An existing complete or relative path is accepted.
#   3. A filename under CSV_DIR is accepted.
#
# Output:
#   Resolved CSV path written to standard output.
#
# Return status:
#   0 when a file is resolved.
#   1 when the selection is invalid or the file is missing.
#
# Implementation note:
#   local -n creates a nameref to the array whose name was passed as $2.
#
resolve_selected_file() {
    local selection="$1"
    local -n files_ref="$2"

    # A numeric selection refers to the displayed one-based menu.
    if [[ "${selection}" =~ ^[0-9]+$ ]]; then
        local index=$((selection - 1))

        if (( index < 0 || index >= ${#files_ref[@]} )); then
            echo "[ERROR] Selection ${selection} is outside the available range." >&2
            return 1
        fi

        printf "%s\n" "${CSV_DIR}/${files_ref[index]}"
        return 0
    fi

    # Accept the selection directly when it points to an existing file.
    if [[ -f "${selection}" ]]; then
        realpath "${selection}"
        return 0
    fi

    # Otherwise, interpret the selection as a filename under CSV_DIR.
    if [[ -f "${CSV_DIR}/${selection}" ]]; then
        realpath "${CSV_DIR}/${selection}"
        return 0
    fi

    echo "[ERROR] CSV file not found: ${selection}" >&2
    return 1
}


# ==============================================================================
# FUNCTION: remove_csv_extension
# ==============================================================================
#
# Purpose:
#   Return a CSV filename without its directory or extension.
#
# Arguments:
#   $1
#       CSV path or filename.
#
# Supported extensions:
#   .csv
#   .CSV
#
remove_csv_extension() {
    local filename

    # Remove the directory component.
    filename="$(basename "$1")"

    # Remove a lowercase or uppercase CSV extension.
    filename="${filename%.csv}"
    filename="${filename%.CSV}"

    printf "%s\n" "${filename}"
}


# ==============================================================================
# FUNCTION: prepare_python_environment
# ==============================================================================
#
# Purpose:
#   Prepare the Python runtime used by the embedded comparison analysis.
#
# Steps:
#   1. Purge currently loaded modules when modules are available.
#   2. Attempt to load anaconda3, python, or Python.
#   3. Locate python3 or python.
#   4. Create VENV_DIR when its Python executable is missing.
#   5. Ensure pip is installed.
#   6. Install NumPy and pandas when needed.
#   7. Verify both packages can be imported.
#
# Global variables assigned:
#   SYSTEM_PYTHON
#       Python executable used to create the environment.
#
#   PYTHON_CMD
#       Python executable inside VENV_DIR.
#
prepare_python_environment() {
    echo "[INFO] Preparing Python environment..."

    # Use the cluster module system when available.
    if command -v module >/dev/null 2>&1; then
        # Clear currently loaded modules to reduce environment conflicts.
        module purge >/dev/null 2>&1 || true

        # Attempt several likely Python module names.
        if module load anaconda3 >/dev/null 2>&1; then
            echo "[INFO] Loaded module: anaconda3"
        elif module load python >/dev/null 2>&1; then
            echo "[INFO] Loaded module: python"
        elif module load Python >/dev/null 2>&1; then
            echo "[INFO] Loaded module: Python"
        else
            echo "[WARNING] Could not load a Python or Anaconda module."
            echo "[WARNING] Using the Python installation already available."
        fi
    fi

    # Prefer python3 and fall back to python.
    if command -v python3 >/dev/null 2>&1; then
        SYSTEM_PYTHON="$(command -v python3)"
    elif command -v python >/dev/null 2>&1; then
        SYSTEM_PYTHON="$(command -v python)"
    else
        echo "[ERROR] Python could not be found." >&2
        exit 1
    fi

    echo "[INFO] System Python: ${SYSTEM_PYTHON}"

    # Create the virtual environment only when it is missing or incomplete.
    if [[ ! -x "${VENV_DIR}/bin/python" ]]; then
        echo "[INFO] Creating virtual environment:"
        echo "       ${VENV_DIR}"

        rm -rf "${VENV_DIR}"
        "${SYSTEM_PYTHON}" -m venv "${VENV_DIR}"
    else
        echo "[INFO] Reusing virtual environment:"
        echo "       ${VENV_DIR}"
    fi

    PYTHON_CMD="${VENV_DIR}/bin/python"

    # Confirm the environment was created successfully.
    if [[ ! -x "${PYTHON_CMD}" ]]; then
        echo "[ERROR] Virtual-environment Python was not created correctly." >&2
        exit 1
    fi

    # Install pip when it is unavailable inside the environment.
    if ! "${PYTHON_CMD}" -m pip --version >/dev/null 2>&1; then
        echo "[INFO] Installing pip inside the virtual environment..."
        "${PYTHON_CMD}" -m ensurepip --upgrade
    fi

    # Install NumPy and pandas only when either package cannot be imported.
    if ! "${PYTHON_CMD}" -c "import numpy, pandas" >/dev/null 2>&1; then
        echo "[INFO] Installing NumPy and pandas..."
        "${PYTHON_CMD}" -m pip install numpy pandas
    else
        echo "[INFO] NumPy and pandas are already installed."
    fi

    # Perform a final dependency validation.
    if ! "${PYTHON_CMD}" -c "import numpy, pandas" >/dev/null 2>&1; then
        echo "[ERROR] NumPy and pandas could not be imported." >&2
        exit 1
    fi

    echo "[INFO] Python environment is ready:"
    echo "       ${PYTHON_CMD}"
}


# ==============================================================================
# FUNCTION: submit_interactive_job
# ==============================================================================
#
# Purpose:
#   Run the user-facing file-selection and SLURM-submission stage.
#
# This function:
#   - Requires at least two CSV files.
#   - Displays a numbered file-selection menu.
#   - Accepts a number, filename, or complete path.
#   - Prevents comparing a file with itself.
#   - Builds output and log names.
#   - Submits this script in --run-comparison mode.
#
submit_interactive_job() {
    # Confirm the CSV directory exists.
    if [[ ! -d "${CSV_DIR}" ]]; then
        echo "[ERROR] CSV directory does not exist:" >&2
        echo "        ${CSV_DIR}" >&2
        exit 1
    fi

    # Read all available CSV filenames into an indexed Bash array.
    mapfile -t CSV_FILES < <(list_csv_files)

    # At least two files are required for a comparison.
    if (( ${#CSV_FILES[@]} < 2 )); then
        echo "[ERROR] At least two CSV files are required in:" >&2
        echo "        ${CSV_DIR}" >&2
        exit 1
    fi

    echo
    echo "============================================================"
    echo "Available POI CSV files"
    echo "Directory: ${CSV_DIR}"
    echo "============================================================"

    # Display a one-based numbered menu.
    local i

    for i in "${!CSV_FILES[@]}"; do
        printf "  %3d) %s\n" \
            "$((i + 1))" \
            "${CSV_FILES[i]}"
    done

    # --------------------------------------------------------------------------
    # Select the first CSV
    # --------------------------------------------------------------------------

    echo
    read -r -p "Enter the number or filename for the FIRST POI CSV: " FIRST_SELECTION

    FIRST_CSV="$(
        resolve_selected_file \
            "${FIRST_SELECTION}" \
            CSV_FILES
    )"

    # --------------------------------------------------------------------------
    # Select the second CSV
    # --------------------------------------------------------------------------

    echo
    read -r -p "Enter the number or filename for the SECOND POI CSV: " SECOND_SELECTION

    SECOND_CSV="$(
        resolve_selected_file \
            "${SECOND_SELECTION}" \
            CSV_FILES
    )"

    # Prevent a file from being compared with itself.
    if [[ "${FIRST_CSV}" == "${SECOND_CSV}" ]]; then
        echo "[ERROR] You selected the same CSV file twice." >&2
        exit 1
    fi

    # --------------------------------------------------------------------------
    # Construct comparison output names
    # --------------------------------------------------------------------------

    FIRST_PREFIX="$(remove_csv_extension "${FIRST_CSV}")"
    SECOND_PREFIX="$(remove_csv_extension "${SECOND_CSV}")"

    # Example:
    #
    #   first_sample__vs__second_sample
    #
    OUTPUT_PREFIX="${FIRST_PREFIX}__vs__${SECOND_PREFIX}"

    # Standard output and standard error are written to the same combined log.
    LOG_FILE="${CSV_DIR}/${FIRST_PREFIX}.${SECOND_PREFIX}.poi.log"

    # --------------------------------------------------------------------------
    # Display the final submission configuration
    # --------------------------------------------------------------------------

    echo
    echo "============================================================"
    echo "POI comparison"
    echo "============================================================"
    echo "First CSV:       ${FIRST_CSV}"
    echo "Second CSV:      ${SECOND_CSV}"
    echo "Tolerance:       +/- ${POSITION_TOLERANCE_BP} bp"
    echo "Output prefix:   ${OUTPUT_PREFIX}"
    echo "Log file:        ${LOG_FILE}"
    echo "Resources:       4 CPUs, 8 GB RAM"
    echo "============================================================"
    echo

    # Confirm the script is running in an environment with SLURM access.
    if ! command -v sbatch >/dev/null 2>&1; then
        echo "[ERROR] The sbatch command is unavailable." >&2
        echo "[ERROR] Run this script from an HPC node with SLURM access." >&2
        exit 1
    fi

    # --------------------------------------------------------------------------
    # Submit this script in internal comparison mode
    # --------------------------------------------------------------------------
    #
    # --output and --error:
    #   Send standard output and standard error to the same combined log.
    #
    # --open-mode=truncate:
    #   Replace an existing log with the same name rather than appending.
    #
    # "$0":
    #   Submit this same Bash script.
    #
    SUBMISSION_OUTPUT="$(
        sbatch \
            --output="${LOG_FILE}" \
            --error="${LOG_FILE}" \
            --open-mode=truncate \
            "$0" \
            --run-comparison \
            "${FIRST_CSV}" \
            "${SECOND_CSV}" \
            "${OUTPUT_PREFIX}"
    )"

    echo "${SUBMISSION_OUTPUT}"
    echo
    echo "[INFO] Comparison job submitted."
    echo "[INFO] Log file:"
    echo "       ${LOG_FILE}"
}


# ==============================================================================
# FUNCTION: run_comparison
# ==============================================================================
#
# Purpose:
#   Execute the comparison inside the submitted SLURM job.
#
# Required arguments:
#   $1
#       First CSV path.
#
#   $2
#       Second CSV path.
#
#   $3
#       Output prefix.
#
# This function:
#   - Validates the internal argument count.
#   - Validates both input files.
#   - Prepares the Python environment.
#   - Exports configuration to the embedded Python analysis.
#   - Runs the comparison.
#
run_comparison() {
    # The internal mode requires exactly three arguments.
    if [[ $# -ne 3 ]]; then
        echo "[ERROR] Internal comparison mode requires:" >&2
        echo "        first CSV, second CSV, and output prefix" >&2
        exit 1
    fi

    FIRST_CSV="$1"
    SECOND_CSV="$2"
    OUTPUT_PREFIX="$3"

    # Validate the first input on the compute node.
    if [[ ! -f "${FIRST_CSV}" ]]; then
        echo "[ERROR] First CSV does not exist:" >&2
        echo "        ${FIRST_CSV}" >&2
        exit 1
    fi

    # Validate the second input on the compute node.
    if [[ ! -f "${SECOND_CSV}" ]]; then
        echo "[ERROR] Second CSV does not exist:" >&2
        echo "        ${SECOND_CSV}" >&2
        exit 1
    fi

    # Create or reuse the Python environment.
    prepare_python_environment

    # Export all values required by the quoted Python here-document.
    export FIRST_CSV
    export SECOND_CSV
    export CSV_DIR
    export OUTPUT_PREFIX
    export POSITION_TOLERANCE_BP

    echo
    echo "============================================================"
    echo "Starting POI comparison"
    echo "============================================================"
    echo "SLURM job ID:  ${SLURM_JOB_ID:-not available}"
    echo "First CSV:     ${FIRST_CSV}"
    echo "Second CSV:    ${SECOND_CSV}"
    echo "Tolerance:     +/- ${POSITION_TOLERANCE_BP} bp"
    echo "Output prefix: ${OUTPUT_PREFIX}"
    echo "============================================================"

    # --------------------------------------------------------------------------
    # Embedded Python analysis
    # --------------------------------------------------------------------------
    #
    # The delimiter is quoted as 'PYTHON'. Bash therefore does not expand
    # variables inside the Python code. Configuration is passed safely through
    # environment variables instead.
    #
    "${PYTHON_CMD}" <<'PYTHON'
"""
Compare two genomic point-of-interest CSV files.

The analysis:

1. Reads and validates both CSV files.
2. Determines a representative position for each POI.
3. Generates all candidate same-chromosome pairs within the tolerance.
4. Selects unique one-to-one matches beginning with the closest pair.
5. Writes overlap, unmatched, and summary CSV files.
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd


# =============================================================================
# CONFIGURATION FROM THE BASH ENVIRONMENT
# =============================================================================

first_path = Path(os.environ["FIRST_CSV"])
second_path = Path(os.environ["SECOND_CSV"])
output_dir = Path(os.environ["CSV_DIR"])
output_prefix = os.environ["OUTPUT_PREFIX"]
tolerance = int(os.environ["POSITION_TOLERANCE_BP"])


def read_poi_csv(path: Path, label: str) -> pd.DataFrame:
    """
    Read and validate a top-POI CSV file.

    Required columns
    ----------------
    chrom
        Chromosome identifier.

    start
        Zero-based interval start.

    end
        Interval end.

    Optional columns
    ----------------
    position
        Representative genomic position. When absent, the midpoint of start
        and end is calculated.

    rank
        POI rank. When absent, ranks are assigned from current row order.

    Parameters
    ----------
    path
        Path to the input CSV.

    label
        Human-readable label used in status and error messages.

    Returns
    -------
    pandas.DataFrame
        Validated POI rows with integer start, end, and position columns.
    """
    print(f"[INFO] Reading {label} CSV: {path}")

    # Load the CSV using pandas' standard comma-separated parser.
    df = pd.read_csv(path)

    # Verify the minimum required schema.
    required_columns = {"chrom", "start", "end"}
    missing_columns = required_columns.difference(df.columns)

    if missing_columns:
        raise ValueError(
            f"{label} CSV is missing required columns: "
            f"{', '.join(sorted(missing_columns))}"
        )

    # Work on a copy so the original object is not modified unexpectedly.
    df = df.copy()

    # Normalize chromosome values as strings.
    df["chrom"] = df["chrom"].astype(str)

    # Convert coordinates to numeric values.
    #
    # Invalid values become NaN and are removed below.
    df["start"] = pd.to_numeric(
        df["start"],
        errors="coerce",
    )

    df["end"] = pd.to_numeric(
        df["end"],
        errors="coerce",
    )

    # Prefer the explicit representative-position column.
    if "position" in df.columns:
        df["position"] = pd.to_numeric(
            df["position"],
            errors="coerce",
        )
    else:
        # Use integer midpoint when position is unavailable.
        df["position"] = (
            (df["start"] + df["end"]) // 2
        )

    original_count = len(df)

    # Remove rows missing required genomic values.
    df = df.dropna(
        subset=[
            "chrom",
            "start",
            "end",
            "position",
        ]
    ).copy()

    # Retain only valid genomic intervals and positions.
    df = df[
        (df["start"] >= 0)
        & (df["end"] > df["start"])
        & (df["position"] >= 0)
    ].copy()

    removed_count = original_count - len(df)

    if removed_count > 0:
        print(
            f"[WARNING] Removed {removed_count} invalid rows "
            f"from the {label} CSV."
        )

    # Store genomic coordinates as 64-bit integers.
    df["start"] = df["start"].astype(np.int64)
    df["end"] = df["end"].astype(np.int64)
    df["position"] = df["position"].astype(np.int64)

    # Generate ranks from current row order when no rank column exists.
    if "rank" not in df.columns:
        df.insert(
            0,
            "rank",
            range(1, len(df) + 1),
        )

    print(f"[INFO] Valid {label} POIs: {len(df)}")

    return df


def build_candidate_matches(
    first_df: pd.DataFrame,
    second_df: pd.DataFrame,
    max_distance: int,
) -> pd.DataFrame:
    """
    Find all same-chromosome POI pairs within the distance tolerance.

    This function does not enforce one-to-one matching. It generates every
    possible candidate pair satisfying the chromosome and distance rules.

    Parameters
    ----------
    first_df
        Validated POIs from the first file.

    second_df
        Validated POIs from the second file.

    max_distance
        Maximum allowed absolute position difference in base pairs.

    Returns
    -------
    pandas.DataFrame
        One row per candidate pair. An empty DataFrame is returned when no
        candidate pairs are found.
    """
    candidate_rows = []

    # Process one chromosome at a time.
    for chrom, first_chrom_df in first_df.groupby(
        "chrom",
        sort=False,
    ):
        # Retain only second-file POIs from the same chromosome.
        second_chrom_df = second_df[
            second_df["chrom"] == chrom
        ]

        if second_chrom_df.empty:
            continue

        # Compare each first-file POI with all second-file POIs from the same
        # chromosome.
        for first_index, first_row in first_chrom_df.iterrows():
            distances = (
                second_chrom_df["position"]
                - first_row["position"]
            ).abs()

            # Candidate rows are those within the configured tolerance.
            matching_rows = second_chrom_df[
                distances <= max_distance
            ]

            for second_index, second_row in matching_rows.iterrows():
                candidate_rows.append(
                    {
                        # Original DataFrame indices are retained so matched
                        # and unmatched rows can be identified later.
                        "first_index": first_index,
                        "second_index": second_index,

                        "chrom": chrom,

                        "first_rank": first_row["rank"],
                        "second_rank": second_row["rank"],

                        "first_start": first_row["start"],
                        "first_end": first_row["end"],
                        "first_position": first_row["position"],

                        "second_start": second_row["start"],
                        "second_end": second_row["end"],
                        "second_position": second_row["position"],

                        # Absolute genomic distance between representative
                        # positions.
                        "distance_bp": abs(
                            int(first_row["position"])
                            - int(second_row["position"])
                        ),

                        # Preserve optional modification statistics when they
                        # are present in the source CSV files.
                        "first_frac_mod": first_row.get(
                            "frac_mod",
                            np.nan,
                        ),
                        "second_frac_mod": second_row.get(
                            "frac_mod",
                            np.nan,
                        ),
                        "first_BrdU_pct": first_row.get(
                            "BrdU_pct",
                            np.nan,
                        ),
                        "second_BrdU_pct": second_row.get(
                            "BrdU_pct",
                            np.nan,
                        ),
                    }
                )

    if not candidate_rows:
        return pd.DataFrame()

    return pd.DataFrame(candidate_rows)


def choose_one_to_one_matches(
    candidates: pd.DataFrame,
) -> pd.DataFrame:
    """
    Select one-to-one POI matches beginning with the closest candidate pairs.

    A first-file POI and a second-file POI can each appear in no more than one
    selected match.

    Candidate sort order
    --------------------
    1. Smallest distance.
    2. Smallest first-file rank.
    3. Smallest second-file rank.

    A stable mergesort is used so tied rows retain deterministic ordering.

    Parameters
    ----------
    candidates
        Candidate same-chromosome pairs within the distance tolerance.

    Returns
    -------
    pandas.DataFrame
        Selected unique matches with a generated overlap_id column.
    """
    if candidates.empty:
        return pd.DataFrame()

    # Sort closest pairs first. Rank columns provide deterministic tie-breakers.
    candidates = candidates.sort_values(
        by=[
            "distance_bp",
            "first_rank",
            "second_rank",
        ],
        ascending=True,
        kind="mergesort",
    )

    # Track POIs that have already been assigned to a selected pair.
    used_first = set()
    used_second = set()

    selected_rows = []

    for row in candidates.itertuples(index=False):
        # Skip candidates whose first POI is already matched.
        if row.first_index in used_first:
            continue

        # Skip candidates whose second POI is already matched.
        if row.second_index in used_second:
            continue

        selected_rows.append(row._asdict())

        used_first.add(row.first_index)
        used_second.add(row.second_index)

    matches = pd.DataFrame(selected_rows)

    # Assign a human-readable sequential overlap identifier.
    if not matches.empty:
        matches.insert(
            0,
            "overlap_id",
            range(1, len(matches) + 1),
        )

    return matches


# =============================================================================
# READ AND VALIDATE INPUT TABLES
# =============================================================================

first_df = read_poi_csv(
    first_path,
    "first",
)

second_df = read_poi_csv(
    second_path,
    "second",
)

print("[INFO] Searching for overlapping POIs...")


# =============================================================================
# GENERATE ALL CANDIDATE OVERLAPS
# =============================================================================

candidate_matches = build_candidate_matches(
    first_df=first_df,
    second_df=second_df,
    max_distance=tolerance,
)

print(
    f"[INFO] Candidate pairs within +/- {tolerance:,} bp: "
    f"{len(candidate_matches)}"
)


# =============================================================================
# SELECT UNIQUE ONE-TO-ONE MATCHES
# =============================================================================

matches = choose_one_to_one_matches(candidate_matches)

if matches.empty:
    matched_first_indices = set()
    matched_second_indices = set()
else:
    matched_first_indices = set(matches["first_index"])
    matched_second_indices = set(matches["second_index"])


# =============================================================================
# IDENTIFY UNMATCHED POIs
# =============================================================================

first_unmatched = first_df.loc[
    ~first_df.index.isin(matched_first_indices)
].copy()

second_unmatched = second_df.loc[
    ~second_df.index.isin(matched_second_indices)
].copy()

# Add the source filename to make unmatched tables self-describing.
first_unmatched.insert(
    0,
    "source_file",
    first_path.name,
)

second_unmatched.insert(
    0,
    "source_file",
    second_path.name,
)


# =============================================================================
# CONSTRUCT OUTPUT PATHS
# =============================================================================

overlap_output = (
    output_dir
    / f"{output_prefix}_overlaps.csv"
)

first_unmatched_output = (
    output_dir
    / f"{output_prefix}_first_unmatched.csv"
)

second_unmatched_output = (
    output_dir
    / f"{output_prefix}_second_unmatched.csv"
)

summary_output = (
    output_dir
    / f"{output_prefix}_summary.csv"
)


# =============================================================================
# STANDARDIZE OVERLAP OUTPUT COLUMNS
# =============================================================================

match_columns = [
    "overlap_id",
    "chrom",
    "first_rank",
    "second_rank",
    "first_start",
    "first_end",
    "first_position",
    "second_start",
    "second_end",
    "second_position",
    "distance_bp",
    "first_frac_mod",
    "second_frac_mod",
    "first_BrdU_pct",
    "second_BrdU_pct",
]

# Even when there are no overlaps, create an empty CSV with the expected header.
if matches.empty:
    matches = pd.DataFrame(columns=match_columns)
else:
    matches = matches[match_columns]


# =============================================================================
# WRITE OVERLAP AND UNMATCHED TABLES
# =============================================================================

matches.to_csv(
    overlap_output,
    index=False,
)

first_unmatched.to_csv(
    first_unmatched_output,
    index=False,
)

second_unmatched.to_csv(
    second_unmatched_output,
    index=False,
)


# =============================================================================
# CALCULATE OVERLAP PERCENTAGES
# =============================================================================
#
# The same overlap count is divided independently by each file's valid POI
# total. The two percentages may differ when the files contain different
# numbers of valid POIs.
#

first_overlap_pct = (
    100.0 * len(matches) / len(first_df)
    if len(first_df) > 0
    else 0.0
)

second_overlap_pct = (
    100.0 * len(matches) / len(second_df)
    if len(second_df) > 0
    else 0.0
)


# =============================================================================
# CREATE AND WRITE THE SUMMARY TABLE
# =============================================================================

summary = pd.DataFrame(
    [
        {
            "first_file": first_path.name,
            "second_file": second_path.name,
            "tolerance_bp": tolerance,
            "first_total_pois": len(first_df),
            "second_total_pois": len(second_df),
            "overlapping_pairs": len(matches),
            "first_unmatched_pois": len(first_unmatched),
            "second_unmatched_pois": len(second_unmatched),
            "first_overlap_pct": first_overlap_pct,
            "second_overlap_pct": second_overlap_pct,
        }
    ]
)

summary.to_csv(
    summary_output,
    index=False,
)


# =============================================================================
# PRINT A HUMAN-READABLE SUMMARY TO THE SLURM LOG
# =============================================================================

print()
print("============================================================")
print("POI overlap summary")
print("============================================================")
print(f"First file:          {first_path.name}")
print(f"Second file:         {second_path.name}")
print(f"Tolerance:           +/- {tolerance:,} bp")
print(f"First POIs:          {len(first_df)}")
print(f"Second POIs:         {len(second_df)}")
print(f"Overlapping pairs:   {len(matches)}")
print(f"First overlap rate:  {first_overlap_pct:.2f}%")
print(f"Second overlap rate: {second_overlap_pct:.2f}%")
print()

if matches.empty:
    print("No overlapping POIs were found.")
else:
    # Print a concise subset of overlap columns in the log.
    print(
        matches[
            [
                "overlap_id",
                "chrom",
                "first_position",
                "second_position",
                "distance_bp",
                "first_BrdU_pct",
                "second_BrdU_pct",
            ]
        ].to_string(index=False)
    )

print()
print("[SUCCESS] Output files:")
print(f"  Overlaps:          {overlap_output}")
print(f"  First unmatched:   {first_unmatched_output}")
print(f"  Second unmatched:  {second_unmatched_output}")
print(f"  Summary:           {summary_output}")
PYTHON

    echo
    echo "============================================================"
    echo "Comparison complete"
    echo "============================================================"
}


# ==============================================================================
# MAIN SCRIPT ROUTING
# ==============================================================================
#
# --run-comparison:
#   Execute the internal SLURM comparison stage.
#
# Any other invocation:
#   Enter the interactive file-selection and submission stage.
#

if [[ "${1:-}" == "--run-comparison" ]]; then
    shift
    run_comparison "$@"
else
    submit_interactive_job
fi