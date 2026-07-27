]#!/bin/bash

# ==============================================================================
# BrdU READ-PERCENTAGE WORKFLOW
# ==============================================================================
#
# Purpose
# -------
# This SLURM workflow calculates how many mapped primary reads in a BAM file
# contain at least one passing BrdU modified-base call.
#
# The workflow also reports:
#
#   - Total mapped primary reads.
#   - Reads containing MM tags.
#   - Reads containing ML tags.
#   - Reads containing any parsed modified-base calls.
#   - Reads containing at least one passing BrdU call.
#   - Total number of passing BrdU calls.
#   - Percentage of reads containing BrdU.
#   - Number of BrdU calls per 100 total mapped primary reads.
#   - Modified-base parsing errors.
#
# This Bash wrapper uses:
#
#   samtools
#       Validates the BAM and counts mapped primary alignments.
#
#   count_brdu_read_stats.py
#       Uses pysam to inspect MM/ML modified-base information and count BrdU
#       calls that pass the selected modification-probability threshold.
#
# Execution modes
# ---------------
#
# Submission mode:
#
#   bash src/utils/calculate_brdu_read_percentage.sh
#
#   The script lists BAM files, prompts for a BAM and threshold, validates the
#   choices, and submits itself to SLURM.
#
# Compute-job mode:
#
#   --run-job
#
#   This internal mode is normally entered automatically by the submitted
#   SLURM job. It performs validation, counting, calculations, and reporting.
#
# Read-count definition
# ---------------------
#
# The denominator is the number of mapped primary alignments returned by:
#
#   samtools view -c -F 2308
#
# Flag value 2308 excludes:
#
#   0x4    unmapped alignments
#   0x100  secondary alignments
#   0x800  supplementary alignments
#
# Therefore, secondary and supplementary representations of a read are not
# counted as separate primary reads.
#
# BrdU definition
# ---------------
#
# The Python helper interprets BrdU as:
#
#   Canonical base: T
#   Modification code: b
#
# A read is counted as BrdU-positive when it contains at least one BrdU call
# whose ML-derived probability passes the selected threshold.
#
# Important metric distinction
# ----------------------------
#
# pct_reads_with_brdu_calls:
#
#   reads_with_brdu_calls / total_reads * 100
#
# This is the percentage of mapped primary reads containing one or more passing
# BrdU calls.
#
# pct_brdu_calls_per_total_reads:
#
#   brdu_calls / total_reads * 100
#
# This is not a percentage of unique reads. A single read may contribute more
# than one BrdU call, so this value may exceed pct_reads_with_brdu_calls and can
# theoretically exceed 100.
#
# Outputs
# -------
#
# All outputs are written under:
#
#   workflow_root/logs/read_pct
#
# Human-readable report:
#
#   <bam_prefix>.BrdU_read_percentage.threshold_<threshold>.<timestamp>.log
#
# Machine-readable TSV:
#
#   <bam_prefix>.BrdU_read_percentage.threshold_<threshold>.<timestamp>.tsv
#
# Combined workflow log:
#
#   <bam_prefix>.threshold_<threshold>.<job_id>.workflow.log
#
# SLURM standard-output and standard-error files are created during submission.
#
# ==============================================================================


# ==============================================================================
# SLURM RESOURCE REQUESTS
# ==============================================================================
#
# Requested resources:
#   - 4 CPU cores
#   - 8 GB RAM
#   - 2-hour runtime limit
#   - normal partition
#
#SBATCH --job-name=brdu_read_pct
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --partition=normal


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
# RESOLVE THE WORKFLOW ROOT
# ==============================================================================
#
# Resolution order:
#
#   1. A path supplied through --workflow-root.
#   2. BRDU_READ_PCT_WORKFLOW_ROOT from the environment.
#   3. Two directories above this script.
#
# The third method assumes this script is stored at:
#
#   workflow_root/src/utils/calculate_brdu_read_percentage.sh
#

WORKFLOW_ROOT_ARG=""

# Search all arguments for an explicitly supplied --workflow-root value.
for ((arg_i = 1; arg_i <= $#; arg_i++)); do
    if [[ "${!arg_i}" == "--workflow-root" ]]; then
        next_arg_i=$((arg_i + 1))
        WORKFLOW_ROOT_ARG="${!next_arg_i:-}"
        break
    fi
done

# Resolve the final workflow-root path.
if [[ -n "$WORKFLOW_ROOT_ARG" ]]; then
    WORKFLOW_ROOT="$(cd "$WORKFLOW_ROOT_ARG" && pwd)"
elif [[ -n "${BRDU_READ_PCT_WORKFLOW_ROOT:-}" ]]; then
    WORKFLOW_ROOT="$BRDU_READ_PCT_WORKFLOW_ROOT"
else
    WORKFLOW_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
fi


# ==============================================================================
# WORKFLOW PATHS AND DEFAULTS
# ==============================================================================
#
# All input, helper-script, environment, and output paths are derived from the
# resolved workflow root.
#

# Directory containing BAM files available to the workflow.
BAM_DIR="$WORKFLOW_ROOT/data/bam"

# Absolute path to this Bash wrapper.
#
# The script submits itself back to SLURM using --run-job.
SCRIPT_PATH="$WORKFLOW_ROOT/src/utils/calculate_brdu_read_percentage.sh"

# Python helper that parses MM/ML tags and counts BrdU calls.
PYTHON_STATS_SCRIPT="$WORKFLOW_ROOT/src/utils/count_brdu_read_stats.py"

# Directory containing all logs and summary TSV files.
LOG_DIR="$WORKFLOW_ROOT/logs/read_pct"

# Workflow-specific Python virtual environment.
VENV_DIR="$WORKFLOW_ROOT/src/utils/.brdu_read_pct_env"

# Default BrdU modification-probability threshold.
DEFAULT_MOD_THRESHOLD="0.60"


# ==============================================================================
# FUNCTION: usage
# ==============================================================================
#
# Purpose:
#   Print command-line usage information.
#
# Optional positional arguments:
#   $1
#       BAM filename or path.
#
#   $2
#       BrdU modification threshold.
#
usage() {
    echo "Usage: bash src/utils/calculate_brdu_read_percentage.sh [BAM] [mod_threshold]"
    echo
    echo "BAM is resolved by filename under data/bam."
    echo "If arguments are omitted, the script prompts for inputs before submitting a SLURM job."
    echo "mod_threshold defaults to $DEFAULT_MOD_THRESHOLD."
}


# ==============================================================================
# FUNCTION: print_error
# ==============================================================================
#
# Purpose:
#   Write a consistently formatted error message to standard error.
#
# Arguments:
#   $*
#       Error-message text.
#
print_error() {
    echo "[ERROR] $*" >&2
}


# ==============================================================================
# FUNCTION: print_info
# ==============================================================================
#
# Purpose:
#   Write a consistently formatted informational message to standard output.
#
# Arguments:
#   $*
#       Informational message text.
#
print_info() {
    echo "[INFO] $*"
}


# ==============================================================================
# FUNCTION: absolute_existing_file
# ==============================================================================
#
# Purpose:
#   Convert an existing file path to an absolute path.
#
# Arguments:
#   $1
#       Existing file path.
#
# Output:
#   Absolute path written to standard output.
#
absolute_existing_file() {
    local path="$1"
    local dir
    local file

    # Resolve the containing directory to an absolute path.
    dir="$(cd "$(dirname "$path")" && pwd)"

    # Preserve the original filename.
    file="$(basename "$path")"

    # Print the normalized path.
    printf '%s/%s\n' "$dir" "$file"
}


# ==============================================================================
# FUNCTION: resolve_bam_file
# ==============================================================================
#
# Purpose:
#   Locate the requested BAM file.
#
# Search order:
#   1. Use the supplied path directly when it exists.
#   2. Search workflow_root/data/bam using the supplied basename.
#
# Arguments:
#   $1
#       BAM filename or path.
#
# Output:
#   Absolute BAM path when found.
#
# Return status:
#   0 when the BAM is found.
#   1 when the BAM cannot be found.
#
resolve_bam_file() {
    local bam_input="$1"
    local bam_path

    # Accept a complete path when it points to an existing file.
    if [[ -f "$bam_input" ]]; then
        absolute_existing_file "$bam_input"
        return 0
    fi

    # Otherwise, search by basename under the workflow BAM directory.
    bam_path="$BAM_DIR/$(basename "$bam_input")"

    if [[ -f "$bam_path" ]]; then
        absolute_existing_file "$bam_path"
        return 0
    fi

    return 1
}


# ==============================================================================
# FUNCTION: require_existing_dir
# ==============================================================================
#
# Purpose:
#   Verify that a required directory exists.
#
# Arguments:
#   $1
#       Human-readable directory name used in error messages.
#
#   $2
#       Directory path.
#
# Behavior:
#   The workflow exits when the directory is unavailable.
#
require_existing_dir() {
    local name="$1"
    local path="$2"

    if [[ ! -d "$path" ]]; then
        print_error "$name does not exist: $path"
        print_error "Confirm the project path is mounted on the compute node."
        exit 1
    fi
}


# ==============================================================================
# FUNCTION: normalize_mod_threshold
# ==============================================================================
#
# Purpose:
#   Validate and normalize the BrdU probability threshold.
#
# Accepted examples:
#   0
#   0.5
#   .5
#   0.60
#   1
#   1.0
#   B:0.6
#   b:0.6
#
# Processing:
#   - Removes an optional B: or b: prefix.
#   - Adds a leading zero to values such as .6.
#   - Requires a numeric value from 0 through 1.
#   - Prints a normalized numeric value without the B: prefix.
#
# Return status:
#   0 for a valid threshold.
#   1 for an invalid threshold.
#
normalize_mod_threshold() {
    local threshold_input="$1"

    # Remove an uppercase B: prefix when present.
    local threshold_value="${threshold_input#B:}"

    # Remove a lowercase b: prefix when present.
    threshold_value="${threshold_value#b:}"

    # Convert values such as ".6" into "0.6".
    if [[ "$threshold_value" =~ ^\.[0-9]+$ ]]; then
        threshold_value="0$threshold_value"
    fi

    # Require a valid representation between 0 and 1.
    if ! [[ "$threshold_value" =~ ^(0(\.[0-9]+)?|1(\.0+)?)$ ]]; then
        return 1
    fi

    # Perform a final numeric range check and normalize formatting.
    awk -v value="$threshold_value" 'BEGIN {
        if (value < 0 || value > 1) {
            exit 1
        }

        printf "%.12g\n", value + 0
    }'
}


# ==============================================================================
# FUNCTION: threshold_label
# ==============================================================================
#
# Purpose:
#   Convert a numeric threshold into a filename-safe label.
#
# Examples:
#   0.5  -> 0p5
#   0.60 -> 0p6
#   0.75 -> 0p75
#
# Arguments:
#   $1
#       Normalized numeric threshold.
#
# Output:
#   Filename-safe threshold label.
#
threshold_label() {
    local threshold="$1"

    awk -v value="$threshold" 'BEGIN {
        formatted = sprintf("%.12g", value + 0)
        gsub(/\./, "p", formatted)
        printf "%s\n", formatted
    }'
}


# ==============================================================================
# FUNCTION: load_required_modules
# ==============================================================================
#
# Purpose:
#   Load Python and Samtools through the cluster module system.
#
# Behavior:
#   - Sources the environment-modules initialization script when available.
#   - Attempts to load Python 3.13.7.
#   - Falls back to the generic Python module.
#   - Loads Samtools.
#   - Verifies that python3 and samtools are executable.
#
# If the module command is unavailable, the function continues and expects
# Python and Samtools to already be available through PATH.
#
load_required_modules() {
    # Initialize the module command when the cluster provides this script.
    if [[ -f /etc/profile.d/modules.sh ]]; then
        # shellcheck source=/dev/null
        source /etc/profile.d/modules.sh
    fi

    if command -v module >/dev/null 2>&1; then
        module load python/3.13.7 || module load python
        module load samtools
    else
        echo "[WARN] Environment modules are not available in this shell."
        echo "[WARN] Continuing with python and samtools from PATH."
    fi

    # Confirm Python is available after module loading.
    if ! command -v python3 >/dev/null 2>&1; then
        print_error "python3 was not found."
        exit 1
    fi

    # Confirm Samtools is available after module loading.
    if ! command -v samtools >/dev/null 2>&1; then
        print_error "samtools was not found."
        exit 1
    fi
}


# ==============================================================================
# FUNCTION: prepare_python_environment
# ==============================================================================
#
# Purpose:
#   Create or reuse a workflow-specific Python environment.
#
# Environment location:
#
#   workflow_root/src/utils/.brdu_read_pct_env
#
# Behavior:
#   - Creates the environment only when its Python executable is missing.
#   - Uses --system-site-packages so cluster-provided packages remain visible.
#   - Activates the environment.
#   - Installs pysam 0.23.3 only when pysam cannot already be imported.
#   - Verifies that pysam can be imported after installation.
#
prepare_python_environment() {
    if [[ ! -x "$VENV_DIR/bin/python" ]]; then
        print_info "Creating Python virtual environment: $VENV_DIR"

        python3 -m venv \
            --system-site-packages \
            --clear \
            "$VENV_DIR"
    fi

    # Activate the workflow-specific environment.
    # shellcheck source=/dev/null
    source "$VENV_DIR/bin/activate"

    # Install the pinned pysam version only when pysam is unavailable.
    if ! python -c "import pysam" >/dev/null 2>&1; then
        print_info "Installing pysam in utility environment."

        python -m pip install \
            "pysam==0.23.3" \
            --quiet
    fi

    # Confirm installation succeeded.
    if ! python -c "import pysam" >/dev/null 2>&1; then
        print_error "pysam is not available in $VENV_DIR."
        exit 1
    fi
}


# ==============================================================================
# FUNCTION: cleanup
# ==============================================================================
#
# Purpose:
#   Remove the temporary TSV generated by count_brdu_read_stats.py.
#
# Notes:
#   This function is registered with an EXIT trap. Cleanup is therefore
#   attempted during successful completion and most error exits.
#
cleanup() {
    if [[ -n "${TEMP_STATS_TSV:-}" && -f "$TEMP_STATS_TSV" ]]; then
        rm -f "$TEMP_STATS_TSV"
    fi
}

# Register cleanup for script exit.
trap cleanup EXIT


# ==============================================================================
# FUNCTION: run_job
# ==============================================================================
#
# Purpose:
#   Perform the compute stage inside the SLURM allocation.
#
# Internal named arguments:
#   --workflow-root PATH
#   --bam PATH
#   --mod-threshold VALUE
#
# Major stages:
#   1. Parse and validate internal job arguments.
#   2. Validate required directories, BAM, and Python helper.
#   3. Construct timestamped output paths.
#   4. Load software and prepare the Python environment.
#   5. Validate the BAM with samtools quickcheck.
#   6. Count mapped primary reads with Samtools.
#   7. Count MM, ML, modified-base, and BrdU statistics with pysam.
#   8. Calculate the reported percentages.
#   9. Write a machine-readable TSV.
#  10. Write a human-readable report.
#
run_job() {
    # --------------------------------------------------------------------------
    # Job arguments and defaults
    # --------------------------------------------------------------------------

    local bam_path=""

    # Default modification threshold used when no explicit value is passed.
    local mod_threshold="$DEFAULT_MOD_THRESHOLD"
    local raw_mod_threshold

    # BAM naming variables.
    local bam_name
    local bam_prefix

    # Filename-safe threshold label.
    local safe_threshold

    # Output naming and paths.
    local timestamp
    local log_file
    local summary_tsv
    local workflow_log

    # Software-version information.
    local samtools_version
    local python_version
    local pysam_version

    # Primary analysis results.
    local total_reads
    local brdu_reads
    local brdu_calls
    local pct_reads_with_brdu
    local pct_brdu_calls_per_read

    # Additional pysam statistics.
    local total_reads_from_pysam
    local reads_with_mm
    local reads_with_ml
    local reads_with_mod_calls
    local parse_errors
    local min_ml

    # Remove the leading --run-job argument.
    shift

    # --------------------------------------------------------------------------
    # Parse named arguments passed by submit_workflow()
    # --------------------------------------------------------------------------

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --workflow-root)
                shift 2
                ;;
            --bam)
                bam_path="$2"
                shift 2
                ;;
            --mod-threshold)
                mod_threshold="$2"
                shift 2
                ;;
            *)
                print_error "Unknown job argument: $1"
                exit 1
                ;;
        esac
    done

    # --------------------------------------------------------------------------
    # Validate the modification threshold
    # --------------------------------------------------------------------------

    raw_mod_threshold="$mod_threshold"

    if ! mod_threshold="$(normalize_mod_threshold "$raw_mod_threshold")"; then
        print_error "Invalid mod threshold in job: $raw_mod_threshold"
        print_error "Expected a number from 0 to 1, for example 0.5 or 0.6."
        exit 1
    fi

    # --------------------------------------------------------------------------
    # Validate workflow paths and required inputs
    # --------------------------------------------------------------------------

    require_existing_dir "Workflow root" "$WORKFLOW_ROOT"
    require_existing_dir "BAM directory" "$BAM_DIR"

    # Create the output directory when it does not already exist.
    mkdir -p "$LOG_DIR"

    # The analysis cannot run without the Python tag-parsing helper.
    if [[ ! -f "$PYTHON_STATS_SCRIPT" ]]; then
        print_error "Python helper script not found: $PYTHON_STATS_SCRIPT"
        exit 1
    fi

    # Validate the BAM again on the compute node.
    if [[ ! -f "$bam_path" ]]; then
        print_error "BAM file not found on compute node: $bam_path"
        exit 1
    fi

    # --------------------------------------------------------------------------
    # Construct timestamped and threshold-specific output paths
    # --------------------------------------------------------------------------

    bam_name="$(basename "$bam_path")"
    bam_prefix="${bam_name%.bam}"

    # Replace the decimal point in the threshold for safe filenames.
    safe_threshold="$(threshold_label "$mod_threshold")"

    # Timestamp format:
    #
    #   YYYYMMDD_HHMMSS
    #
    timestamp="$(date +"%Y%m%d_%H%M%S")"

    # Human-readable result report.
    log_file="$LOG_DIR/${bam_prefix}.BrdU_read_percentage.threshold_${safe_threshold}.${timestamp}.log"

    # Machine-readable one-row result table.
    summary_tsv="$LOG_DIR/${bam_prefix}.BrdU_read_percentage.threshold_${safe_threshold}.${timestamp}.tsv"

    # Combined runtime log containing messages from the full compute stage.
    workflow_log="$LOG_DIR/${bam_prefix}.threshold_${safe_threshold}.${SLURM_JOB_ID:-manual}.workflow.log"

    # Temporary key/value TSV produced by the Python helper.
    TEMP_STATS_TSV="$(mktemp "$LOG_DIR/.${bam_prefix}.brdu_stats.XXXXXX.tsv")"

    # --------------------------------------------------------------------------
    # Configure combined workflow logging
    # --------------------------------------------------------------------------
    #
    # All subsequent standard output and standard error are:
    #
    #   - Preserved in the active SLURM streams.
    #   - Appended to the workflow log through tee.
    #
    exec > >(tee -a "$workflow_log") 2>&1

    print_info "BrdU read-percentage analysis started: $(date)"
    print_info "Workflow log: $workflow_log"

    # --------------------------------------------------------------------------
    # Load software and prepare the isolated Python environment
    # --------------------------------------------------------------------------

    load_required_modules
    prepare_python_environment

    # Record exact software versions for reproducibility.
    samtools_version="$(samtools --version | head -n 1)"
    python_version="$(python --version 2>&1)"
    pysam_version="$(python -c "import pysam; print(pysam.__version__)")"

    print_info "Using $samtools_version"
    print_info "Using $python_version"
    print_info "Using pysam $pysam_version"
    print_info "Checking BAM file."

    # --------------------------------------------------------------------------
    # Validate the BAM structure
    # --------------------------------------------------------------------------
    #
    # samtools quickcheck checks the BAM header and EOF structure. It is a fast
    # integrity check and does not read every alignment record.
    #
    if ! samtools quickcheck -v "$bam_path"; then
        print_error "The BAM failed samtools quickcheck validation."
        exit 1
    fi

    # --------------------------------------------------------------------------
    # Count mapped primary alignments
    # --------------------------------------------------------------------------
    #
    # samtools options:
    #
    #   -@
    #       Number of additional Samtools threads.
    #
    #   -c
    #       Return only the number of alignments.
    #
    #   -F 2308
    #       Exclude unmapped, secondary, and supplementary alignments.
    #
    # This total is used as the denominator for the reported metrics.
    #
    print_info "Counting mapped primary reads."

    total_reads="$(
        samtools view \
            -@ "${SLURM_CPUS_PER_TASK:-4}" \
            -c \
            -F 2308 \
            "$bam_path"
    )"

    # Require a nonnegative integer result.
    if [[ ! "$total_reads" =~ ^[0-9]+$ ]]; then
        print_error "samtools did not return a valid read count."
        exit 1
    fi

    # Prevent division by zero and reject BAMs without primary mapped reads.
    if (( total_reads == 0 )); then
        print_error "No mapped primary reads were found in the selected BAM."
        exit 1
    fi

    # --------------------------------------------------------------------------
    # Count modified-base and BrdU statistics with pysam
    # --------------------------------------------------------------------------

    print_info "Mapped primary reads: $total_reads"
    print_info "Counting BrdU-positive reads and total passing BrdU calls."

    # The Python helper writes a key/value TSV to standard output.
    python "$PYTHON_STATS_SCRIPT" \
        "$bam_path" \
        --mod-threshold "$mod_threshold" \
        > "$TEMP_STATS_TSV"

    # Require nonempty helper output.
    if [[ ! -s "$TEMP_STATS_TSV" ]]; then
        print_error "The Python BrdU counter did not produce output."
        exit 1
    fi

    # --------------------------------------------------------------------------
    # Parse statistics returned by the Python helper
    # --------------------------------------------------------------------------
    #
    # Expected keys:
    #
    #   total_primary_mapped_from_pysam
    #       Primary mapped alignments independently counted by pysam.
    #
    #   reads_with_mm
    #       Reads containing an MM modified-base tag.
    #
    #   reads_with_ml
    #       Reads containing an ML probability tag.
    #
    #   reads_with_mod_calls
    #       Reads for which at least one modified-base call was parsed.
    #
    #   reads_with_brdu_calls
    #       Reads containing at least one passing T+b BrdU call.
    #
    #   brdu_calls
    #       Total number of passing BrdU calls across all selected reads.
    #
    #   parse_errors
    #       Number of reads or records whose modification tags could not be
    #       parsed successfully.
    #
    #   min_ml
    #       Minimum integer ML value corresponding to the probability threshold.
    #

    total_reads_from_pysam="$(
        awk -F'\t' \
            '$1 == "total_primary_mapped_from_pysam" {print $2}' \
            "$TEMP_STATS_TSV"
    )"

    reads_with_mm="$(
        awk -F'\t' \
            '$1 == "reads_with_mm" {print $2}' \
            "$TEMP_STATS_TSV"
    )"

    reads_with_ml="$(
        awk -F'\t' \
            '$1 == "reads_with_ml" {print $2}' \
            "$TEMP_STATS_TSV"
    )"

    reads_with_mod_calls="$(
        awk -F'\t' \
            '$1 == "reads_with_mod_calls" {print $2}' \
            "$TEMP_STATS_TSV"
    )"

    brdu_reads="$(
        awk -F'\t' \
            '$1 == "reads_with_brdu_calls" {print $2}' \
            "$TEMP_STATS_TSV"
    )"

    brdu_calls="$(
        awk -F'\t' \
            '$1 == "brdu_calls" {print $2}' \
            "$TEMP_STATS_TSV"
    )"

    parse_errors="$(
        awk -F'\t' \
            '$1 == "parse_errors" {print $2}' \
            "$TEMP_STATS_TSV"
    )"

    min_ml="$(
        awk -F'\t' \
            '$1 == "min_ml" {print $2}' \
            "$TEMP_STATS_TSV"
    )"

    # Require valid integer BrdU counts before calculations proceed.
    if [[ ! "$brdu_reads" =~ ^[0-9]+$ || ! "$brdu_calls" =~ ^[0-9]+$ ]]; then
        print_error "Could not determine BrdU read/call counts."
        exit 1
    fi

    # Samtools and pysam should normally return the same number of primary
    # mapped reads. A difference is reported but does not automatically stop
    # the workflow.
    if [[ "$total_reads_from_pysam" != "$total_reads" ]]; then
        echo "[WARN] samtools and pysam primary mapped read counts differ."
        echo "[WARN] samtools: $total_reads; pysam: $total_reads_from_pysam"
    fi

    # --------------------------------------------------------------------------
    # Calculate summary metrics
    # --------------------------------------------------------------------------
    #
    # Metric 1: percentage of unique reads containing BrdU
    #
    #   reads_with_brdu_calls / total_reads * 100
    #
    pct_reads_with_brdu="$(
        awk \
            -v brdu_reads="$brdu_reads" \
            -v total_reads="$total_reads" \
            'BEGIN {
                printf "%.6f", (brdu_reads / total_reads) * 100
            }'
    )"

    # Metric 2: number of BrdU calls per 100 total reads
    #
    #   brdu_calls / total_reads * 100
    #
    # This metric does not represent unique reads because one read can contain
    # multiple passing BrdU calls.
    pct_brdu_calls_per_read="$(
        awk \
            -v brdu_calls="$brdu_calls" \
            -v total_reads="$total_reads" \
            'BEGIN {
                printf "%.6f", (brdu_calls / total_reads) * 100
            }'
    )"

    # --------------------------------------------------------------------------
    # Write the machine-readable one-row summary TSV
    # --------------------------------------------------------------------------

    {
        echo "sample_bam	total_reads	reads_with_mm	reads_with_ml	reads_with_mod_calls	reads_with_brdu_calls	brdu_calls	pct_reads_with_brdu_calls	pct_brdu_calls_per_total_reads	mod_threshold	min_ml	parse_errors"

        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$bam_name" \
            "$total_reads" \
            "$reads_with_mm" \
            "$reads_with_ml" \
            "$reads_with_mod_calls" \
            "$brdu_reads" \
            "$brdu_calls" \
            "$pct_reads_with_brdu" \
            "$pct_brdu_calls_per_read" \
            "$mod_threshold" \
            "$min_ml" \
            "$parse_errors"
    } > "$summary_tsv"

    # --------------------------------------------------------------------------
    # Write the human-readable analysis report
    # --------------------------------------------------------------------------
    #
    # tee writes the report to the timestamped result log and also sends it to
    # the combined workflow/SLURM output stream.
    #

    {
        echo "============================================================"
        echo "BrdU read-percentage analysis"
        echo "============================================================"
        echo
        echo "Analysis date:                         $(date +"%Y-%m-%d %H:%M:%S")"
        echo "Workflow root:                         $WORKFLOW_ROOT"
        echo "Input BAM:                             $bam_path"
        echo
        echo "Samtools version:                      $samtools_version"
        echo "Python version:                        $python_version"
        echo "pysam version:                         $pysam_version"
        echo
        echo "BrdU canonical base:                   T"
        echo "BrdU modification code:                b"
        echo "BrdU modification threshold:           $mod_threshold"
        echo "Minimum ML value counted:              $min_ml"
        echo
        echo "Read-count definition:"
        echo "  Primary mapped reads only"
        echo "  Excluded flags: unmapped, secondary, and supplementary"
        echo
        echo "Total reads:                           $total_reads"
        echo "Reads with MM tags:                    $reads_with_mm"
        echo "Reads with ML tags:                    $reads_with_ml"
        echo "Reads with any modified-base calls:    $reads_with_mod_calls"
        echo "Reads with at least one BrdU call:     $brdu_reads"
        echo "Total passing BrdU calls:              $brdu_calls"
        echo "Modified-base parse errors:            $parse_errors"
        echo
        echo "Percent reads with BrdU calls:"
        echo "  ($brdu_reads / $total_reads) x 100 = $pct_reads_with_brdu%"
        echo
        echo "Percent BrdU calls per total reads:"
        echo "  ($brdu_calls / $total_reads) x 100 = $pct_brdu_calls_per_read%"
        echo
        echo "Summary TSV:                           $summary_tsv"
        echo "============================================================"
    } | tee "$log_file"

    # --------------------------------------------------------------------------
    # Finish the compute job
    # --------------------------------------------------------------------------

    print_info "Analysis completed successfully."
    print_info "Log written to: $log_file"
    print_info "Summary TSV written to: $summary_tsv"

    # Exit the workflow-specific Python environment.
    deactivate

    print_info "BrdU read-percentage analysis finished: $(date)"
}


# ==============================================================================
# FUNCTION: submit_workflow
# ==============================================================================
#
# Purpose:
#   Provide the user-facing submission flow.
#
# Positional arguments:
#   $1
#       Optional BAM filename or path.
#
#   $2
#       Optional BrdU modification threshold.
#
# Behavior:
#   - Lists BAM files under data/bam.
#   - Accepts either a menu number or filename.
#   - Prompts for and validates the BrdU threshold.
#   - Confirms sbatch is available.
#   - Submits this script back to SLURM in --run-job mode.
#
submit_workflow() {
    # --------------------------------------------------------------------------
    # Optional positional arguments
    # --------------------------------------------------------------------------

    local bam_input="${1:-}"
    local mod_threshold_input="${2:-}"

    # --------------------------------------------------------------------------
    # Submission variables
    # --------------------------------------------------------------------------

    local bam_path
    local bam_name
    local bam_prefix
    local mod_threshold
    local safe_threshold
    local job_id
    local log_job_id
    local bam_files
    local selection

    # --------------------------------------------------------------------------
    # Validate required workflow directories
    # --------------------------------------------------------------------------

    require_existing_dir "Workflow root" "$WORKFLOW_ROOT"
    require_existing_dir "BAM directory" "$BAM_DIR"

    mkdir -p "$LOG_DIR"

    # --------------------------------------------------------------------------
    # Select the BAM file
    # --------------------------------------------------------------------------

    if [[ -z "$bam_input" ]]; then
        # Read all BAM filenames directly under data/bam into an indexed array.
        mapfile -t bam_files < <(
            find "$BAM_DIR" \
                -maxdepth 1 \
                -type f \
                -name "*.bam" \
                -printf "%f\n" |
                sort
        )

        # Stop when no candidate BAM files are available.
        if [[ ${#bam_files[@]} -eq 0 ]]; then
            print_error "No BAM files were found in: $BAM_DIR"
            exit 1
        fi

        # Display a numbered BAM-selection menu.
        echo "[INFO] Available BAM files in $BAM_DIR:"

        for index in "${!bam_files[@]}"; do
            printf "  %3d) %s\n" \
                "$((index + 1))" \
                "${bam_files[index]}"
        done

        echo

        # Accept either a menu number or a BAM filename.
        while true; do
            read -r -p "Select a BAM file by number or enter its filename: " selection

            if [[ -z "$selection" ]]; then
                echo "Please provide a selection."
                continue
            fi

            # Numeric input is interpreted as a menu position.
            if [[ "$selection" =~ ^[0-9]+$ ]]; then
                if (( selection >= 1 && selection <= ${#bam_files[@]} )); then
                    bam_input="${bam_files[selection - 1]}"
                    break
                fi

                echo "Selection must be between 1 and ${#bam_files[@]}."
                continue
            fi

            # Non-numeric input is interpreted as a filename or path.
            bam_input="$selection"
            break
        done
    fi

    # Resolve the selected BAM to an absolute path.
    if ! bam_path="$(resolve_bam_file "$bam_input")"; then
        print_error "BAM file not found under $BAM_DIR: $(basename "$bam_input")"
        exit 1
    fi

    # --------------------------------------------------------------------------
    # Select and validate the BrdU modification threshold
    # --------------------------------------------------------------------------

    if [[ -z "$mod_threshold_input" ]]; then
        read -r -p "Enter the BrdU modification threshold between 0 and 1 [$DEFAULT_MOD_THRESHOLD]: " mod_threshold_input

        # Pressing Enter accepts the configured default.
        mod_threshold_input="${mod_threshold_input:-$DEFAULT_MOD_THRESHOLD}"
    fi

    if ! mod_threshold="$(normalize_mod_threshold "$mod_threshold_input")"; then
        print_error "Invalid mod threshold: $mod_threshold_input"
        print_error "Expected a number from 0 to 1, for example 0.5 or 0.6."
        exit 1
    fi

    # Derive output naming components.
    bam_name="$(basename "$bam_path")"
    bam_prefix="${bam_name%.bam}"
    safe_threshold="$(threshold_label "$mod_threshold")"

    # --------------------------------------------------------------------------
    # Confirm SLURM submission is available
    # --------------------------------------------------------------------------

    if ! command -v sbatch >/dev/null 2>&1; then
        print_error "The sbatch command is unavailable."
        print_error "Run this script on the cluster login node or submit the --run-job mode manually."
        exit 1
    fi

    # --------------------------------------------------------------------------
    # Display the final submission configuration
    # --------------------------------------------------------------------------

    echo
    echo "[INFO] Workflow root: $WORKFLOW_ROOT"
    echo "[INFO] BAM: $bam_path"
    echo "[INFO] BrdU mod threshold: $mod_threshold"
    echo
    echo "Submitting BrdU read-percentage analysis to SLURM..."

    # --------------------------------------------------------------------------
    # Submit this script back to SLURM in compute-job mode
    # --------------------------------------------------------------------------
    #
    # sbatch options:
    #
    #   --parsable
    #       Return a machine-readable job ID.
    #
    #   --job-name
    #       Assign a sample-specific job name.
    #
    #   --chdir
    #       Run the compute job from the workflow root.
    #
    #   --output
    #       Write SLURM standard output to the read_pct log directory.
    #
    #   --error
    #       Write SLURM standard error to a separate file.
    #
    #   --export
    #       Export the current environment and explicitly supply the workflow
    #       root to the compute node.
    #
    job_id="$(sbatch --parsable \
        --job-name="brdu_read_pct_${bam_prefix}" \
        --chdir="$WORKFLOW_ROOT" \
        --output="$LOG_DIR/${bam_prefix}.threshold_${safe_threshold}.%j.slurm.log" \
        --error="$LOG_DIR/${bam_prefix}.threshold_${safe_threshold}.%j.slurm.err" \
        --export=ALL,BRDU_READ_PCT_WORKFLOW_ROOT="$WORKFLOW_ROOT" \
        "$SCRIPT_PATH" \
        --run-job \
        --workflow-root "$WORKFLOW_ROOT" \
        --bam "$bam_path" \
        --mod-threshold "$mod_threshold")"

    # Some SLURM configurations append cluster information after a semicolon.
    # Retain only the leading job-ID component for output path reporting.
    log_job_id="${job_id%%;*}"

    # --------------------------------------------------------------------------
    # Report the submitted job and expected output paths
    # --------------------------------------------------------------------------

    echo "[INFO] Submitted SLURM job: $job_id"

    echo "[INFO] Expected SLURM log:  $LOG_DIR/${bam_prefix}.threshold_${safe_threshold}.${log_job_id}.slurm.log"

    echo "[INFO] Expected SLURM err:  $LOG_DIR/${bam_prefix}.threshold_${safe_threshold}.${log_job_id}.slurm.err"

    echo "[INFO] Expected workflow log pattern:"

    echo "[INFO]   $LOG_DIR/${bam_prefix}.threshold_${safe_threshold}.${log_job_id}.workflow.log"

    echo "[INFO] Expected output log/TSV pattern:"

    echo "[INFO]   $LOG_DIR/${bam_prefix}.BrdU_read_percentage.threshold_${safe_threshold}.<timestamp>.log"

    echo "[INFO]   $LOG_DIR/${bam_prefix}.BrdU_read_percentage.threshold_${safe_threshold}.<timestamp>.tsv"
}


# ==============================================================================
# MAIN ROUTING
# ==============================================================================
#
# -h or --help:
#   Print usage information and exit.
#
# --run-job:
#   Execute the compute stage. This option is normally supplied internally by
#   submit_workflow().
#
# Any other invocation:
#   Enter the interactive submission stage.
#

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

if [[ "${1:-}" == "--run-job" ]]; then
    run_job "$@"
else
    submit_workflow "$@"
fi