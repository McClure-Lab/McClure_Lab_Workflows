#!/bin/bash

# ==============================================================================
# TOP BrdU POINTS-OF-INTEREST SELECTION WORKFLOW
# ==============================================================================
#
# Purpose
# -------
# This SLURM workflow combines positive- and negative-strand BrdU BEDGRAPH
# files, calculates a strand-combined fraction modified at each genomic
# coordinate, and selects up to 20 high-ranking points of interest (POIs).
#
# The workflow is designed for BEDGRAPH files stored under:
#
#   /beevol/home/pineirok/workflows/data/bedgraph
#
# Run interactively with:
#
#   bash src/utils/find_top_20_poi.sh
#
# Interactive stage
# -----------------
# The launcher:
#
#   1. Lists all BEDGRAPH files directly under data/bedgraph.
#   2. Prompts for the positive-strand BEDGRAPH.
#   3. Prompts for the negative-strand BEDGRAPH.
#   4. Prevents the same file from being selected twice.
#   5. Determines a shared sample prefix from the filenames.
#   6. Constructs the output CSV name.
#   7. Submits this same script to SLURM with --run-analysis.
#
# Analysis stage
# --------------
# The submitted job:
#
#   1. Loads or locates Python.
#   2. Creates or reuses the project virtual environment.
#   3. Installs NumPy and pandas when missing.
#   4. Reads and validates both BEDGRAPH files.
#   5. Normalizes W303 GenBank chromosome accessions.
#   6. Collapses duplicate coordinates within each strand.
#   7. Outer-merges positive and negative coordinates.
#   8. Sums Nmod and Nvalid_cov across strands.
#   9. Recalculates combined frac_mod and BrdU percentage.
#  10. Filters candidate positions by the configured fraction threshold.
#  11. Applies a same-chromosome exclusion window around selected POIs.
#  12. Writes the ranked POI table to data/csv.
#
# BEDGRAPH formats
# ----------------
# Preferred six-column format:
#
#   chrom  start  end  frac_mod  Nmod  Nvalid_cov
#
# A five-column format is also accepted:
#
#   chrom  start  end  frac_mod  Nmod
#
# For five-column inputs, coverage is estimated as:
#
#   Nvalid_cov = Nmod / frac_mod
#
# Six-column input is preferred because it preserves the original coverage
# rather than reconstructing it from rounded values.
#
# Strand-combined calculation
# ---------------------------
# Positive and negative rows are matched by:
#
#   chrom, start, end
#
# Missing strand counts are filled with zero. The combined values are:
#
#   combined_Nmod =
#       positive_Nmod + negative_Nmod
#
#   combined_Nvalid_cov =
#       positive_Nvalid_cov + negative_Nvalid_cov
#
#   frac_mod =
#       combined_Nmod / combined_Nvalid_cov
#
#   BrdU_pct =
#       frac_mod * 100
#
# POI selection rules
# -------------------
# Candidates must satisfy:
#
#   0 < frac_mod <= 0.75
#
# The default maximum therefore excludes positions above 75% modified.
#
# Candidates are ranked by:
#
#   1. Higher frac_mod.
#   2. Higher combined valid coverage.
#   3. Higher combined modified count.
#   4. Chromosome.
#   5. Start coordinate.
#
# After selecting one POI, another candidate on the same chromosome is skipped
# when its midpoint lies within +/- 10,000 bp of any previously selected POI.
#
# Output
# ------
# The output is written as:
#
#   data/csv/<shared_prefix>_top20_poi.csv
#
# The CSV includes ranks, coordinates, combined and strand-specific counts,
# BrdU percentage, and the exclusion window around each selected midpoint.
#
# ==============================================================================


# ==============================================================================
# SLURM RESOURCE REQUESTS
# ==============================================================================
#
# Requested resources:
#
#   Job name:       top20_brdU_poi
#   Tasks:          1
#   CPU cores:      32
#   Memory:         32 GB
#   Runtime limit:  2 hours
#
# Standard output and standard error are written under data/csv.
#
#SBATCH --job-name=top20_brdU_poi
#SBATCH --output=/beevol/home/pineirok/workflows/data/csv/top20_brdU_poi_%j.out
#SBATCH --error=/beevol/home/pineirok/workflows/data/csv/top20_brdU_poi_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=02:00:00


# ==============================================================================
# BASH SAFETY SETTINGS
# ==============================================================================
#
# -e
#   Exit when an unhandled command fails.
#
# -u
#   Treat references to unset variables as errors.
#
# -o pipefail
#   Treat a pipeline as failed when any command in that pipeline fails.
#
set -euo pipefail


# =============================================================================
# Project settings
# =============================================================================

# Root directory of the workflow project.
PROJECT_DIR="/beevol/home/pineirok/workflows"

# Directory containing positive- and negative-strand BEDGRAPH files.
BEDGRAPH_DIR="${PROJECT_DIR}/data/bedgraph"

# Directory receiving the ranked POI CSV and SLURM logs.
CSV_DIR="${PROJECT_DIR}/data/csv"

# Shared project Python virtual environment.
VENV_DIR="${PROJECT_DIR}/.venv"

# Maximum number of POIs to select.
TOP_N=20

# Maximum allowed fraction modified.
#
# 0.75 corresponds to 75%.
THRESHOLD_FRAC=0.75

# After selecting a POI, skip any other candidate within 10 kb
# to the left or right on the same chromosome.
EXCLUSION_WINDOW_BP=10000

# Create the output directory when necessary.
mkdir -p "${CSV_DIR}"


# =============================================================================
# Helper functions
# =============================================================================


# ==============================================================================
# FUNCTION: load_required_modules
# ==============================================================================
#
# Purpose:
#   Locate Python, create or reuse the project virtual environment, and ensure
#   NumPy and pandas are installed.
#
# Module-loading order:
#
#   1. anaconda3
#   2. python
#   3. Python
#
# When environment modules are unavailable, the function uses Python already
# available through PATH.
#
# Global variables assigned:
#
#   SYSTEM_PYTHON
#       Python executable used to create the virtual environment.
#
#   PYTHON_CMD
#       Python executable inside VENV_DIR.
#
load_required_modules() {
    echo "[INFO] Loading Python module and preparing virtual environment..."

    # Load an available Python or Anaconda module when modules are supported.
    if command -v module >/dev/null 2>&1; then
        module purge >/dev/null 2>&1 || true

        if module load anaconda3 >/dev/null 2>&1; then
            echo "[INFO] Loaded module: anaconda3"
        elif module load python >/dev/null 2>&1; then
            echo "[INFO] Loaded module: python"
        elif module load Python >/dev/null 2>&1; then
            echo "[INFO] Loaded module: Python"
        else
            echo "[WARNING] Could not load an Anaconda or Python module."
            echo "[WARNING] Using the Python installation already available."
        fi
    else
        echo "[WARNING] The module command is unavailable."
        echo "[WARNING] Using the Python installation already available."
    fi

    # Find a usable system Python executable.
    if command -v python3 >/dev/null 2>&1; then
        SYSTEM_PYTHON="$(command -v python3)"
    elif command -v python >/dev/null 2>&1; then
        SYSTEM_PYTHON="$(command -v python)"
    else
        echo "[ERROR] Python was not found." >&2
        exit 1
    fi

    echo "[INFO] System Python: ${SYSTEM_PYTHON}"

    # Confirm that this Python installation supports virtual environments.
    if ! "${SYSTEM_PYTHON}" -m venv --help >/dev/null 2>&1; then
        echo "[ERROR] The selected Python installation does not support:" >&2
        echo "        python -m venv" >&2
        echo "[ERROR] Load a Python module that includes the venv module." >&2
        exit 1
    fi

    # Create the virtual environment only when it does not already exist.
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

    # Confirm that the environment Python was created successfully.
    if [[ ! -x "${PYTHON_CMD}" ]]; then
        echo "[ERROR] Virtual-environment Python was not created correctly:" >&2
        echo "        ${PYTHON_CMD}" >&2
        exit 1
    fi

    echo "[INFO] Virtual-environment Python: ${PYTHON_CMD}"

    # Ensure pip exists inside the virtual environment.
    if ! "${PYTHON_CMD}" -m pip --version >/dev/null 2>&1; then
        echo "[INFO] pip was not found in the virtual environment."
        echo "[INFO] Attempting to install pip with ensurepip..."

        "${PYTHON_CMD}" -m ensurepip --upgrade
    fi

    # Upgrade packaging tools when possible.
    #
    # Failure to upgrade is treated as a warning because the currently
    # installed versions may still be able to install the required packages.
    echo "[INFO] Checking pip, setuptools, and wheel..."

    if ! "${PYTHON_CMD}" -m pip install --upgrade pip setuptools wheel; then
        echo "[WARNING] Could not upgrade pip, setuptools, or wheel."
        echo "[WARNING] Continuing with the versions currently installed."
    fi

    # Determine which required libraries are missing.
    MISSING_PACKAGES=()

    if ! "${PYTHON_CMD}" -c "import numpy" >/dev/null 2>&1; then
        MISSING_PACKAGES+=("numpy")
    fi

    if ! "${PYTHON_CMD}" -c "import pandas" >/dev/null 2>&1; then
        MISSING_PACKAGES+=("pandas")
    fi

    # Install only missing packages.
    if (( ${#MISSING_PACKAGES[@]} > 0 )); then
        echo "[INFO] Installing missing Python packages:"
        printf "       %s\n" "${MISSING_PACKAGES[@]}"

        "${PYTHON_CMD}" -m pip install "${MISSING_PACKAGES[@]}"
    else
        echo "[INFO] NumPy and pandas are already installed."
    fi

    # Final import validation.
    if ! "${PYTHON_CMD}" -c "import numpy, pandas" >/dev/null 2>&1; then
        echo "[ERROR] NumPy and pandas could not be imported from:" >&2
        echo "        ${VENV_DIR}" >&2
        exit 1
    fi

    # Report the installed library versions.
    #
    # The quoted here-document prevents Bash from expanding any Python text.
    "${PYTHON_CMD}" - <<'PYTHON_VERSION_CHECK'
import numpy
import pandas

print(f"[INFO] NumPy version: {numpy.__version__}")
print(f"[INFO] pandas version: {pandas.__version__}")
PYTHON_VERSION_CHECK
}


# ==============================================================================
# FUNCTION: list_bedgraph_files
# ==============================================================================
#
# Purpose:
#   List BEDGRAPH files directly under BEDGRAPH_DIR.
#
# Behavior:
#
#   - Does not search subdirectories.
#   - Matches the .bedgraph extension case-insensitively.
#   - Prints filenames without directory components.
#   - Sorts the results alphabetically.
#
list_bedgraph_files() {
    find "${BEDGRAPH_DIR}" \
        -maxdepth 1 \
        -type f \
        -iname "*.bedgraph" \
        -printf "%f\n" |
        sort
}


# ==============================================================================
# FUNCTION: resolve_selected_file
# ==============================================================================
#
# Purpose:
#   Resolve a menu number, filename, or complete BEDGRAPH path.
#
# Arguments:
#
#   $1
#       User selection.
#
#   $2
#       Name of the Bash array containing available filenames.
#
# Resolution order:
#
#   1. Numeric input is treated as a one-based menu selection.
#   2. An existing supplied path is accepted.
#   3. A filename under BEDGRAPH_DIR is accepted.
#
# Output:
#   Resolved file path.
#
# Return status:
#
#   0 when the file is resolved.
#   1 when the selection is invalid.
#
resolve_selected_file() {
    local selection="$1"

    # Create a nameref to the array whose name was passed as the second
    # argument.
    local -n files_ref="$2"

    # Resolve a numeric menu selection.
    if [[ "${selection}" =~ ^[0-9]+$ ]]; then
        local index=$((selection - 1))

        if (( index < 0 || index >= ${#files_ref[@]} )); then
            echo "[ERROR] Selection ${selection} is outside the available range." >&2
            return 1
        fi

        printf "%s\n" "${BEDGRAPH_DIR}/${files_ref[index]}"
        return 0
    fi

    # Accept a complete or relative path.
    if [[ -f "${selection}" ]]; then
        realpath "${selection}"
        return 0
    fi

    # Otherwise, interpret the selection as a filename under BEDGRAPH_DIR.
    if [[ -f "${BEDGRAPH_DIR}/${selection}" ]]; then
        realpath "${BEDGRAPH_DIR}/${selection}"
        return 0
    fi

    echo "[ERROR] BEDGRAPH file not found: ${selection}" >&2
    return 1
}


# ==============================================================================
# FUNCTION: remove_bedgraph_extension
# ==============================================================================
#
# Purpose:
#   Remove the directory and a common BEDGRAPH filename extension.
#
# Supported extensions:
#
#   .bedgraph
#   .bedGraph
#   .BEDGRAPH
#
remove_bedgraph_extension() {
    local filename

    filename="$(basename "$1")"

    filename="${filename%.bedgraph}"
    filename="${filename%.bedGraph}"
    filename="${filename%.BEDGRAPH}"

    printf "%s\n" "${filename}"
}


# ==============================================================================
# FUNCTION: remove_strand_suffix
# ==============================================================================
#
# Purpose:
#   Remove common positive- or negative-strand suffixes from a filename prefix.
#
# Examples handled include:
#
#   sample_positive
#   sample_negative
#   sample_pos
#   sample_neg
#   sample.positive
#   sample-negative
#   sample_positive_strand
#   sample_negative_strand
#   sample_strand_positive
#   sample_strand_negative
#
# Arguments:
#
#   $1
#       Filename without its BEDGRAPH extension.
#
# Output:
#   Sample prefix with the strand suffix removed.
#
remove_strand_suffix() {
    local filename="$1"

    filename="$(
        printf "%s\n" "${filename}" |
            sed -E '
                s/([._-])(positive|negative|pos|neg)([._-]strand)?$//I;
                s/([._-])strand([._-])(positive|negative|pos|neg)$//I;
                s/[._-]+$//
            '
    )"

    printf "%s\n" "${filename}"
}


# ==============================================================================
# FUNCTION: determine_shared_prefix
# ==============================================================================
#
# Purpose:
#   Determine the shared sample prefix for the selected positive and negative
#   BEDGRAPH files.
#
# Preferred behavior:
#
#   - Remove each BEDGRAPH extension.
#   - Remove a recognized strand suffix.
#   - Require the cleaned prefixes to match exactly.
#
# Fallback behavior:
#
#   - Calculate the longest character-by-character shared prefix.
#   - Remove trailing positive/negative strand words and separators.
#   - Emit a warning explaining that fallback naming was used.
#
# Arguments:
#
#   $1
#       Positive-strand BEDGRAPH path.
#
#   $2
#       Negative-strand BEDGRAPH path.
#
# Output:
#   Shared output prefix.
#
determine_shared_prefix() {
    local positive_file="$1"
    local negative_file="$2"

    local positive_base
    local negative_base
    local positive_prefix
    local negative_prefix

    positive_base="$(remove_bedgraph_extension "${positive_file}")"
    negative_base="$(remove_bedgraph_extension "${negative_file}")"

    positive_prefix="$(remove_strand_suffix "${positive_base}")"
    negative_prefix="$(remove_strand_suffix "${negative_base}")"

    # Use the cleaned prefix when both filenames reduce to the same value.
    if [[ -n "${positive_prefix}" &&
          -n "${negative_prefix}" &&
          "${positive_prefix}" == "${negative_prefix}" ]]; then

        printf "%s\n" "${positive_prefix}"
        return 0
    fi

    # Fallback: calculate the longest shared filename prefix.
    local common=""
    local max_length=${#positive_base}

    if (( ${#negative_base} < max_length )); then
        max_length=${#negative_base}
    fi

    local i

    for ((i = 0; i < max_length; i++)); do
        if [[ "${positive_base:i:1}" == "${negative_base:i:1}" ]]; then
            common+="${positive_base:i:1}"
        else
            break
        fi
    done

    # Clean strand words and trailing separators from the fallback prefix.
    common="$(
        printf "%s\n" "${common}" |
            sed -E '
                s/(positive|negative|pos|neg)$//I;
                s/[._-]+$//
            '
    )"

    if [[ -z "${common}" ]]; then
        echo "[ERROR] Could not determine a shared filename prefix." >&2
        echo "[ERROR] Positive file: ${positive_base}" >&2
        echo "[ERROR] Negative file: ${negative_base}" >&2
        return 1
    fi

    echo "[WARNING] Exact strand suffix matching was not found." >&2
    echo "[WARNING] Using the longest shared prefix: ${common}" >&2

    printf "%s\n" "${common}"
}


# ==============================================================================
# FUNCTION: submit_interactive_job
# ==============================================================================
#
# Purpose:
#   Run the user-facing BEDGRAPH selection and SLURM submission stage.
#
# This function:
#
#   - Requires at least two BEDGRAPH files.
#   - Displays a numbered file menu.
#   - Prompts separately for positive and negative files.
#   - Prevents selecting the same file twice.
#   - Determines the shared output prefix.
#   - Displays the active POI-selection settings.
#   - Submits this script using --run-analysis.
#
submit_interactive_job() {
    # Validate the BEDGRAPH input directory.
    if [[ ! -d "${BEDGRAPH_DIR}" ]]; then
        echo "[ERROR] BEDGRAPH directory does not exist:" >&2
        echo "        ${BEDGRAPH_DIR}" >&2
        exit 1
    fi

    # Read available BEDGRAPH filenames into an indexed Bash array.
    mapfile -t BEDGRAPH_FILES < <(list_bedgraph_files)

    # Require at least two candidate files.
    if (( ${#BEDGRAPH_FILES[@]} < 2 )); then
        echo "[ERROR] At least two BEDGRAPH files are required in:" >&2
        echo "        ${BEDGRAPH_DIR}" >&2
        exit 1
    fi

    # Display the available files.
    echo
    echo "============================================================"
    echo "Available BEDGRAPH files"
    echo "Directory: ${BEDGRAPH_DIR}"
    echo "============================================================"

    local i

    for i in "${!BEDGRAPH_FILES[@]}"; do
        printf "  %3d) %s\n" \
            "$((i + 1))" \
            "${BEDGRAPH_FILES[i]}"
    done

    # Select the positive-strand BEDGRAPH.
    echo
    read -r -p "Enter the number or filename for the POSITIVE BEDGRAPH: " POS_SELECTION

    POSITIVE_BEDGRAPH="$(
        resolve_selected_file \
            "${POS_SELECTION}" \
            BEDGRAPH_FILES
    )"

    # Select the negative-strand BEDGRAPH.
    echo
    read -r -p "Enter the number or filename for the NEGATIVE BEDGRAPH: " NEG_SELECTION

    NEGATIVE_BEDGRAPH="$(
        resolve_selected_file \
            "${NEG_SELECTION}" \
            BEDGRAPH_FILES
    )"

    # Prevent the same input from being used for both strands.
    if [[ "${POSITIVE_BEDGRAPH}" == "${NEGATIVE_BEDGRAPH}" ]]; then
        echo "[ERROR] The positive and negative selections are the same file." >&2
        exit 1
    fi

    # Determine the sample prefix used for the output CSV.
    SHARED_PREFIX="$(
        determine_shared_prefix \
            "${POSITIVE_BEDGRAPH}" \
            "${NEGATIVE_BEDGRAPH}"
    )"

    OUTPUT_CSV="${CSV_DIR}/${SHARED_PREFIX}_top20_poi.csv"

    # Display the final submission configuration.
    echo
    echo "============================================================"
    echo "Selected files"
    echo "============================================================"
    echo "Positive BEDGRAPH: ${POSITIVE_BEDGRAPH}"
    echo "Negative BEDGRAPH: ${NEGATIVE_BEDGRAPH}"
    echo "Shared prefix:     ${SHARED_PREFIX}"
    echo "Output CSV:        ${OUTPUT_CSV}"
    echo "Virtual env:       ${VENV_DIR}"
    echo "Threshold:         <= 75%"
    echo "POI exclusion:     +/- 10,000 bp"
    echo "Number of POIs:    ${TOP_N}"
    echo "Resources:         32 CPUs, 32 GB RAM"
    echo "============================================================"
    echo

    # Confirm SLURM submission is available.
    if ! command -v sbatch >/dev/null 2>&1; then
        echo "[ERROR] The sbatch command is unavailable." >&2
        echo "[ERROR] Run this script from an HPC login node with SLURM access." >&2
        exit 1
    fi

    # Submit this same script in internal analysis mode.
    SUBMISSION_OUTPUT="$(
        sbatch \
            "$0" \
            --run-analysis \
            "${POSITIVE_BEDGRAPH}" \
            "${NEGATIVE_BEDGRAPH}" \
            "${OUTPUT_CSV}"
    )"

    echo "${SUBMISSION_OUTPUT}"
    echo
    echo "[INFO] The analysis was submitted to SLURM."
    echo "[INFO] Final CSV will be written to:"
    echo "       ${OUTPUT_CSV}"
}


# ==============================================================================
# FUNCTION: run_analysis
# ==============================================================================
#
# Purpose:
#   Execute BEDGRAPH combination and POI selection inside the submitted SLURM
#   job.
#
# Required arguments:
#
#   $1
#       Positive-strand BEDGRAPH path.
#
#   $2
#       Negative-strand BEDGRAPH path.
#
#   $3
#       Final output CSV path.
#
# The numerical analysis is implemented as embedded Python. Bash configuration
# values are passed through exported environment variables.
#
run_analysis() {
    # Require exactly three internal analysis arguments.
    if [[ $# -ne 3 ]]; then
        echo "[ERROR] Internal analysis mode requires three arguments:" >&2
        echo "        positive BEDGRAPH, negative BEDGRAPH, output CSV" >&2
        exit 1
    fi

    POSITIVE_BEDGRAPH="$1"
    NEGATIVE_BEDGRAPH="$2"
    OUTPUT_CSV="$3"

    # Validate the positive-strand input.
    if [[ ! -f "${POSITIVE_BEDGRAPH}" ]]; then
        echo "[ERROR] Positive BEDGRAPH does not exist:" >&2
        echo "        ${POSITIVE_BEDGRAPH}" >&2
        exit 1
    fi

    # Validate the negative-strand input.
    if [[ ! -f "${NEGATIVE_BEDGRAPH}" ]]; then
        echo "[ERROR] Negative BEDGRAPH does not exist:" >&2
        echo "        ${NEGATIVE_BEDGRAPH}" >&2
        exit 1
    fi

    # Create the output directory when necessary.
    mkdir -p "$(dirname "${OUTPUT_CSV}")"

    # Prepare the Python environment.
    load_required_modules

    # Display the analysis configuration.
    echo
    echo "============================================================"
    echo "Top 20 combined BrdU POI analysis"
    echo "============================================================"
    echo "SLURM job ID:       ${SLURM_JOB_ID:-not available}"
    echo "Positive BEDGRAPH:  ${POSITIVE_BEDGRAPH}"
    echo "Negative BEDGRAPH:  ${NEGATIVE_BEDGRAPH}"
    echo "Output CSV:         ${OUTPUT_CSV}"
    echo "Virtual env:        ${VENV_DIR}"
    echo "Maximum BrdU:       75%"
    echo "Exclusion window:   +/- ${EXCLUSION_WINDOW_BP} bp"
    echo "Requested POIs:     ${TOP_N}"
    echo "CPUs:               ${SLURM_CPUS_PER_TASK:-32}"
    echo "============================================================"

    # Export configuration for the quoted Python here-document.
    export POSITIVE_BEDGRAPH
    export NEGATIVE_BEDGRAPH
    export OUTPUT_CSV
    export THRESHOLD_FRAC
    export EXCLUSION_WINDOW_BP
    export TOP_N

    # --------------------------------------------------------------------------
    # Embedded Python analysis
    # --------------------------------------------------------------------------
    #
    # The quoted PYTHON delimiter prevents Bash variable expansion inside the
    # Python code. Values are read through os.environ instead.
    #
    "${PYTHON_CMD}" <<'PYTHON'
"""
Combine positive- and negative-strand BrdU BEDGRAPH files and select top POIs.

The analysis validates BEDGRAPH coordinates and counts, normalizes W303
chromosome accessions, combines strand-specific counts, recalculates fraction
modified, filters candidate positions, and applies a same-chromosome exclusion
window before writing the ranked output CSV.
"""

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd


# =============================================================================
# Configuration supplied by the Bash wrapper
# =============================================================================

positive_path = Path(
    os.environ["POSITIVE_BEDGRAPH"]
)

negative_path = Path(
    os.environ["NEGATIVE_BEDGRAPH"]
)

output_path = Path(
    os.environ["OUTPUT_CSV"]
)

threshold_frac = float(
    os.environ["THRESHOLD_FRAC"]
)

exclusion_window = int(
    os.environ["EXCLUSION_WINDOW_BP"]
)

top_n = int(
    os.environ["TOP_N"]
)


# =============================================================================
# W303 chromosome-name normalization
# =============================================================================
#
# Convert the W303 GenBank accessions into concise chromosome labels used in
# the POI output.
#
GENBANK_TO_CHR = {
    "CM007964.1": "1",
    "CM007965.1": "2",
    "CM007966.1": "3",
    "CM007967.1": "4",
    "CM007968.1": "5",
    "CM007969.1": "6",
    "CM007970.1": "7",
    "CM007971.1": "8",
    "CM007972.1": "9",
    "CM007973.1": "10",
    "CM007974.1": "11",
    "CM007975.1": "12",
    "CM007976.1": "13",
    "CM007977.1": "14",
    "CM007978.1": "15",
    "CM007979.1": "16",
    "CM007980.1": "p2-micron",
    "CM007981.1": "MT",
}


def read_bedgraph(
    path: Path,
    source_name: str,
) -> pd.DataFrame:
    """
    Read, validate, normalize, and collapse one strand BEDGRAPH.

    Preferred six-column format
    ---------------------------
    chrom
        Chromosome or reference-sequence identifier.

    start
        Interval start.

    end
        Interval end.

    frac_mod
        Fraction modified for the strand.

    Nmod
        Number of modified observations.

    Nvalid_cov
        Valid modification-aware coverage.

    Five-column compatibility
    -------------------------
    A five-column file lacking Nvalid_cov is accepted. Coverage is estimated as:

        Nvalid_cov = Nmod / frac_mod

    The estimate is undefined when frac_mod is zero, so those rows receive NaN
    coverage and are removed during validation.

    Duplicate coordinates
    ---------------------
    Rows sharing chrom, start, and end within one strand file are combined by
    summing Nmod and Nvalid_cov.

    Parameters
    ----------
    path
        BEDGRAPH input path.

    source_name
        Descriptive strand label used in messages.

    Returns
    -------
    pandas.DataFrame
        Unique valid coordinates with Nmod and Nvalid_cov columns.
    """
    print(
        f"[INFO] Reading {source_name} BEDGRAPH: "
        f"{path}"
    )

    try:
        df = pd.read_csv(
            path,
            sep="\t",
            comment="#",
            header=None,
            dtype={
                0: str,
            },
            low_memory=False,
        )
    except Exception as exc:
        raise RuntimeError(
            f"Could not read the {source_name} BEDGRAPH: "
            f"{path}"
        ) from exc

    if df.empty:
        raise ValueError(
            f"The {source_name} BEDGRAPH contains no data: "
            f"{path}"
        )

    column_count = df.shape[1]

    print(
        f"[INFO] {source_name.capitalize()} BEDGRAPH columns detected: "
        f"{column_count}"
    )

    if column_count < 5:
        raise ValueError(
            f"The {source_name} BEDGRAPH has only {column_count} columns. "
            "At least five columns are required: "
            "chrom, start, end, frac_mod, Nmod."
        )

    # Use the first six columns when original valid coverage is present.
    if column_count >= 6:
        df = df.iloc[:, :6].copy()

        df.columns = [
            "chrom",
            "start",
            "end",
            "frac_mod",
            "Nmod",
            "Nvalid_cov",
        ]
    else:
        # Support a five-column file by reconstructing valid coverage.
        df = df.iloc[:, :5].copy()

        df.columns = [
            "chrom",
            "start",
            "end",
            "frac_mod",
            "Nmod",
        ]

        df["frac_mod"] = pd.to_numeric(
            df["frac_mod"],
            errors="coerce",
        )

        df["Nmod"] = pd.to_numeric(
            df["Nmod"],
            errors="coerce",
        )

        df["Nvalid_cov"] = np.where(
            df["frac_mod"] > 0,
            df["Nmod"] / df["frac_mod"],
            np.nan,
        )

        print(
            f"[WARNING] The {source_name} BEDGRAPH contains five columns.",
            file=sys.stderr,
        )

        print(
            "[WARNING] Nvalid_cov is being estimated as Nmod / frac_mod.",
            file=sys.stderr,
        )

        print(
            "[WARNING] A six-column BEDGRAPH is preferred because it "
            "contains the original Nvalid_cov.",
            file=sys.stderr,
        )

    # Convert all coordinate and count fields to numeric values.
    #
    # Unparseable values become NaN and are removed below.
    for column in [
        "start",
        "end",
        "frac_mod",
        "Nmod",
        "Nvalid_cov",
    ]:
        df[column] = pd.to_numeric(
            df[column],
            errors="coerce",
        )

    original_rows = len(df)

    # Replace positive and negative infinity with NaN.
    df = df.replace(
        [
            np.inf,
            -np.inf,
        ],
        np.nan,
    )

    # Remove rows lacking essential values.
    df = df.dropna(
        subset=[
            "chrom",
            "start",
            "end",
            "Nmod",
            "Nvalid_cov",
        ]
    ).copy()

    # Retain only valid genomic coordinates and nonnegative counts.
    df = df[
        (df["start"] >= 0)
        & (df["end"] > df["start"])
        & (df["Nmod"] >= 0)
        & (df["Nvalid_cov"] >= 0)
    ].copy()

    removed_rows = (
        original_rows
        - len(df)
    )

    if removed_rows > 0:
        print(
            f"[WARNING] Removed {removed_rows:,} invalid rows from the "
            f"{source_name} BEDGRAPH.",
            file=sys.stderr,
        )

    # Normalize chromosome identifiers and numeric data types.
    df["chrom"] = (
        df["chrom"]
        .astype(str)
        .replace(GENBANK_TO_CHR)
    )

    df["start"] = df["start"].astype(
        np.int64
    )

    df["end"] = df["end"].astype(
        np.int64
    )

    df["Nmod"] = df["Nmod"].astype(
        float
    )

    df["Nvalid_cov"] = df[
        "Nvalid_cov"
    ].astype(float)

    # Combine duplicate coordinates inside the same strand file.
    #
    # Aggregating before the positive/negative merge ensures each strand has
    # no more than one row per genomic coordinate.
    df = (
        df.groupby(
            [
                "chrom",
                "start",
                "end",
            ],
            as_index=False,
            sort=False,
        )
        .agg(
            Nmod=(
                "Nmod",
                "sum",
            ),
            Nvalid_cov=(
                "Nvalid_cov",
                "sum",
            ),
        )
    )

    print(
        f"[INFO] Valid unique {source_name} positions: "
        f"{len(df):,}"
    )

    return df


def combine_bedgraphs(
    positive_df: pd.DataFrame,
    negative_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Combine positive- and negative-strand counts at matching coordinates.

    An outer merge is used so a coordinate present on only one strand is
    retained. Missing strand counts are filled with zero.

    Combined fraction modified
    --------------------------
    The strand-combined fraction is calculated as:

        combined_Nmod / combined_Nvalid_cov

    Parameters
    ----------
    positive_df
        Validated unique positive-strand coordinates.

    negative_df
        Validated unique negative-strand coordinates.

    Returns
    -------
    pandas.DataFrame
        Combined coordinates with strand-specific and total counts.
    """
    # Rename the strand-specific count columns before the merge.
    positive_df = positive_df.rename(
        columns={
            "Nmod": "positive_Nmod",
            "Nvalid_cov": "positive_Nvalid_cov",
        }
    )

    negative_df = negative_df.rename(
        columns={
            "Nmod": "negative_Nmod",
            "Nvalid_cov": "negative_Nvalid_cov",
        }
    )

    # Outer-merge the two strand tables by exact genomic coordinate.
    #
    # one_to_one validation is possible because duplicate coordinates were
    # already collapsed by read_bedgraph().
    combined = positive_df.merge(
        negative_df,
        on=[
            "chrom",
            "start",
            "end",
        ],
        how="outer",
        validate="one_to_one",
    )

    count_columns = [
        "positive_Nmod",
        "positive_Nvalid_cov",
        "negative_Nmod",
        "negative_Nvalid_cov",
    ]

    # A missing strand contributes zero counts at that coordinate.
    combined[count_columns] = (
        combined[count_columns]
        .fillna(0.0)
    )

    # Sum modified counts across strands.
    combined["combined_Nmod"] = (
        combined["positive_Nmod"]
        + combined["negative_Nmod"]
    )

    # Sum valid coverage across strands.
    combined["combined_Nvalid_cov"] = (
        combined["positive_Nvalid_cov"]
        + combined["negative_Nvalid_cov"]
    )

    # Recalculate fraction modified from the combined counts.
    combined["frac_mod"] = np.where(
        combined["combined_Nvalid_cov"] > 0,
        (
            combined["combined_Nmod"]
            / combined["combined_Nvalid_cov"]
        ),
        np.nan,
    )

    # Convert the fraction to percentage units for reporting.
    combined["BrdU_pct"] = (
        combined["frac_mod"]
        * 100.0
    )

    combined = combined.replace(
        [
            np.inf,
            -np.inf,
        ],
        np.nan,
    )

    # Remove coordinates where a valid fraction could not be calculated.
    combined = combined.dropna(
        subset=[
            "chrom",
            "start",
            "end",
            "frac_mod",
        ]
    ).copy()

    print(
        "[INFO] Combined unique genomic positions: "
        f"{len(combined):,}"
    )

    return combined


def select_top_pois(
    combined_df: pd.DataFrame,
    maximum_fraction: float,
    exclusion_bp: int,
    number_to_select: int,
) -> pd.DataFrame:
    """
    Select the highest-ranking spatially separated BrdU POIs.

    Candidate requirement
    ---------------------
    Candidates must satisfy:

        0 < frac_mod <= maximum_fraction

    Ranking
    -------
    Candidates are sorted by:

    1. Higher frac_mod.
    2. Higher combined_Nvalid_cov.
    3. Higher combined_Nmod.
    4. Chromosome.
    5. Start coordinate.

    Exclusion rule
    --------------
    Once a candidate is selected, another candidate on the same chromosome is
    skipped when its midpoint lies within ``exclusion_bp`` bases of any
    previously selected midpoint.

    Parameters
    ----------
    combined_df
        Combined positive- and negative-strand positions.

    maximum_fraction
        Maximum allowed fraction modified.

    exclusion_bp
        Same-chromosome exclusion distance in base pairs.

    number_to_select
        Maximum number of POIs to retain.

    Returns
    -------
    pandas.DataFrame
        Ranked selected POIs, or an empty DataFrame when none qualify.
    """
    # Apply the positive-fraction and maximum-threshold requirements.
    candidates = combined_df[
        (combined_df["frac_mod"] > 0)
        & (
            combined_df["frac_mod"]
            <= maximum_fraction
        )
    ].copy()

    print(
        "[INFO] Candidates passing 0 < frac_mod <= "
        f"{maximum_fraction:.2f}: "
        f"{len(candidates):,}"
    )

    if candidates.empty:
        return pd.DataFrame()

    # Use the interval midpoint as the representative POI position.
    candidates["position"] = (
        (
            candidates["start"]
            + candidates["end"]
        )
        // 2
    ).astype(np.int64)

    # Highest frac_mod is selected first.
    #
    # Coverage and modified counts break ties, followed by genomic position.
    candidates = candidates.sort_values(
        by=[
            "frac_mod",
            "combined_Nvalid_cov",
            "combined_Nmod",
            "chrom",
            "start",
        ],
        ascending=[
            False,
            False,
            False,
            True,
            True,
        ],
        kind="mergesort",
    )

    selected_rows = []

    # Store previously selected midpoint positions separately for each
    # chromosome.
    selected_positions: dict[
        str,
        list[int],
    ] = {}

    for row in candidates.itertuples(
        index=False
    ):
        chrom = str(row.chrom)
        position = int(row.position)

        chromosome_positions = (
            selected_positions.get(
                chrom,
                [],
            )
        )

        # Determine whether the candidate falls inside an existing exclusion
        # window on the same chromosome.
        too_close = any(
            abs(
                position
                - previously_selected
            )
            <= exclusion_bp
            for previously_selected
            in chromosome_positions
        )

        if too_close:
            continue

        selected_rows.append(
            row._asdict()
        )

        selected_positions.setdefault(
            chrom,
            [],
        ).append(position)

        # Stop after selecting the requested maximum number of POIs.
        if len(selected_rows) >= number_to_select:
            break

    if not selected_rows:
        return pd.DataFrame()

    selected = pd.DataFrame(
        selected_rows
    )

    # Assign one-based rank values in selection order.
    selected.insert(
        0,
        "rank",
        range(
            1,
            len(selected) + 1,
        ),
    )

    # Report the exclusion interval surrounding each POI.
    selected["window_start"] = (
        selected["position"]
        - exclusion_bp
    ).clip(lower=0)

    selected["window_end"] = (
        selected["position"]
        + exclusion_bp
    )

    return selected


# =============================================================================
# Read and validate both strand BEDGRAPH files
# =============================================================================

positive_df = read_bedgraph(
    positive_path,
    "positive",
)

negative_df = read_bedgraph(
    negative_path,
    "negative",
)


# =============================================================================
# Combine positive- and negative-strand counts
# =============================================================================

print(
    "[INFO] Combining positive and negative BEDGRAPH files..."
)

combined_df = combine_bedgraphs(
    positive_df,
    negative_df,
)


# =============================================================================
# Select the ranked POIs
# =============================================================================

print(
    f"[INFO] Selecting up to {top_n} POIs with a "
    f"+/- {exclusion_window:,} bp exclusion window..."
)

top_pois = select_top_pois(
    combined_df=combined_df,
    maximum_fraction=threshold_frac,
    exclusion_bp=exclusion_window,
    number_to_select=top_n,
)


# =============================================================================
# Standardize output columns
# =============================================================================

output_columns = [
    "rank",
    "chrom",
    "start",
    "end",
    "position",
    "frac_mod",
    "BrdU_pct",
    "combined_Nmod",
    "combined_Nvalid_cov",
    "positive_Nmod",
    "positive_Nvalid_cov",
    "negative_Nmod",
    "negative_Nvalid_cov",
    "window_start",
    "window_end",
]

if top_pois.empty:
    print(
        "[WARNING] No POIs passed the selection requirements.",
        file=sys.stderr,
    )

    # Create an empty output with the expected header.
    top_pois = pd.DataFrame(
        columns=output_columns
    )
else:
    top_pois = top_pois[
        output_columns
    ]


# =============================================================================
# Write the output CSV
# =============================================================================

output_path.parent.mkdir(
    parents=True,
    exist_ok=True,
)

top_pois.to_csv(
    output_path,
    index=False,
)


# =============================================================================
# Print a concise POI summary to the SLURM log
# =============================================================================

print()
print("============================================================")
print("Selected POIs")
print("============================================================")

if top_pois.empty:
    print("No POIs were selected.")
else:
    display_columns = [
        "rank",
        "chrom",
        "start",
        "end",
        "frac_mod",
        "BrdU_pct",
        "combined_Nmod",
        "combined_Nvalid_cov",
    ]

    print(
        top_pois[
            display_columns
        ].to_string(
            index=False,
        )
    )

print()
print(
    f"[SUCCESS] CSV saved to: "
    f"{output_path}"
)
PYTHON

    # Report successful completion after the embedded Python analysis exits.
    echo
    echo "============================================================"
    echo "Analysis complete"
    echo "============================================================"
    echo "Output CSV:"
    echo "${OUTPUT_CSV}"
}


# =============================================================================
# Main script behavior
# =============================================================================
#
# --run-analysis
#   Execute the internal SLURM analysis stage.
#
# Any other invocation
#   Enter the interactive BEDGRAPH-selection and job-submission stage.
#

if [[ "${1:-}" == "--run-analysis" ]]; then
    shift
    run_analysis "$@"
else
    submit_interactive_job
fi