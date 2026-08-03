#!/bin/bash

# ==============================================================================
# SLURM RESOURCE REQUESTS
# ==============================================================================
#
# These directives are read by SLURM when the script is submitted with sbatch.
#
# Resources requested:
#   - 12 CPU cores
#   - 16 GB of RAM
#   - Maximum runtime of 8 hours
#   - The normal SLURM partition
#
# The default job name may be overridden later by the sbatch command inside
# submit_workflow().
#
#SBATCH --job-name=genome_browser_brdu
#SBATCH --cpus-per-task=12
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --partition=normal


# ==============================================================================
# BASH SAFETY SETTINGS
# ==============================================================================
#
# -e
#   Exit immediately when a command returns a nonzero exit status.
#
# -u
#   Treat references to unset variables as errors.
#
# -o pipefail
#   If any command in a pipeline fails, the entire pipeline is considered failed.
#
# These settings help prevent the workflow from continuing after an unexpected
# error or from silently using an undefined variable.
#
set -euo pipefail


# ==============================================================================
# WORKFLOW OVERVIEW
# ==============================================================================
#
# This script has two operating modes:
#
#   1. Submission mode
#      The user runs the script directly from the workflow root. The script:
#
#        - Lists available BAM files.
#        - Prompts for the BAM file.
#        - Prompts for the reference FASTA.
#        - Prompts for the Modkit modification threshold.
#        - Prompts for the cell-cycle phase.
#        - Prompts for smoothed or unsmoothed plotting.
#        - Submits a SLURM job that runs this same script in job mode.
#
#   2. Job mode
#      The submitted SLURM job runs this script with --run-job. The script:
#
#        - Loads Python, Samtools, and Modkit.
#        - Creates or reuses a Python virtual environment.
#        - Installs the Python requirements.
#        - Extracts positive- and negative-strand BrdU bedgraph files.
#        - Validates that the bedgraph files were created.
#        - Runs the selected genome-browser plotting script.
#
# Expected project structure:
#
#   workflow_root/
#   ├── data/
#   │   ├── bam/
#   │   ├── bedgraph/
#   │   ├── sorted_bam/
#   │   ├── bed/
#   │   └── ncbi/
#   ├── logs/
#   │   └── genome_browser_workflow/
#   ├── results/
#   │   └── genome_browser_results/
#   └── src/
#       └── genome_browser_workflow/
#           ├── genome_browser_workflow_script.sh
#           ├── raw_data_extraction_on_bam.py
#           ├── genomic_browser_generation.py
#           ├── genomic_browser_generation_unsmoothed.py
#           └── requirements.txt
#
# Typical interactive command:
#
#   bash src/genome_browser_workflow/genome_browser_workflow_script.sh
#
# Optional positional arguments:
#
#   bash src/genome_browser_workflow/genome_browser_workflow_script.sh \
#       sample.bam \
#       output_prefix \
#       0.5
#
# Help:
#
#   bash src/genome_browser_workflow/genome_browser_workflow_script.sh --help
#
# Important:
#   Users should normally not call --run-job manually. That option is used
#   internally when the script submits itself to SLURM.
# ==============================================================================


# ==============================================================================
# DETERMINE THE WORKFLOW ROOT
# ==============================================================================
#
# The workflow root is resolved using the following priority:
#
#   1. A path provided through:
#
#        --workflow-root /path/to/workflow
#
#   2. The GENOME_BROWSER_WORKFLOW_ROOT environment variable.
#
#   3. The directory two levels above this script.
#
# The third method assumes this script is stored at:
#
#   workflow_root/src/genome_browser_workflow/
#
# Scanning all command-line arguments here is necessary because --workflow-root
# may not be the first argument when the script is running in job mode.
#

# Stores a workflow-root path provided directly on the command line.
WORKFLOW_ROOT_ARG=""

# Examine each command-line argument to find --workflow-root.
for ((arg_i = 1; arg_i <= $#; arg_i++)); do
    if [[ "${!arg_i}" == "--workflow-root" ]]; then
        # The value immediately after --workflow-root is the requested path.
        next_arg_i=$((arg_i + 1))
        WORKFLOW_ROOT_ARG="${!next_arg_i:-}"
        break
    fi
done

# Resolve the final workflow-root directory.
if [[ -n "$WORKFLOW_ROOT_ARG" ]]; then
    # Convert the explicitly supplied directory to an absolute path.
    WORKFLOW_ROOT="$(cd "$WORKFLOW_ROOT_ARG" && pwd)"
elif [[ -n "${GENOME_BROWSER_WORKFLOW_ROOT:-}" ]]; then
    # Use the exported environment variable when present.
    WORKFLOW_ROOT="$GENOME_BROWSER_WORKFLOW_ROOT"
else
    # Infer the root from the script's expected location:
    # workflow_root/src/genome_browser_workflow/script.sh
    WORKFLOW_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
fi


# ==============================================================================
# WORKFLOW PATH CONFIGURATION
# ==============================================================================
#
# All workflow directories and default input paths are defined here so the rest
# of the script can use consistent absolute paths.
#

# Directory containing the BAM files presented to the user.
BAM_DIR="$WORKFLOW_ROOT/data/bam"

# Directory where positive- and negative-strand bedgraph files are written.
BEDGRAPH_DIR="$WORKFLOW_ROOT/data/bedgraph"

# Directory reserved for sorted BAM outputs or intermediate sorted BAM files.
#
# This script creates the directory, although sorting behavior may be performed
# by one of the associated Python scripts.
SORTED_BAM_DIR="$WORKFLOW_ROOT/data/sorted_bam"

# Directory containing the Bash wrapper and Python workflow scripts.
SRC_DIR="$WORKFLOW_ROOT/src/genome_browser_workflow"

# Absolute path to this workflow wrapper.
#
# This path is passed back to sbatch so the script can invoke itself in job mode.
SCRIPT_PATH="$SRC_DIR/genome_browser_workflow_script.sh"

# Parent directory for generated genome-browser images.
RESULTS_DIR="$WORKFLOW_ROOT/results/genome_browser_results"

# Directory containing SLURM logs and the combined workflow log.
LOG_DIR="$WORKFLOW_ROOT/logs/genome_browser_workflow"

# Default W303 reference FASTA used by Modkit pileup.
DEFAULT_REF="$WORKFLOW_ROOT/data/ncbi/W303/ncbi_dataset/GCA_002163515.1_ASM216351v1_genomic.fna"

# Optional genomic feature files used by the plotting scripts.
#
# Each feature file is only supplied to the Python plotting script when the
# corresponding file exists.
DEFAULT_G4_BED="$WORKFLOW_ROOT/data/bed/W303_g4_motifs.bed"
DEFAULT_TE_BED="$WORKFLOW_ROOT/data/bed/w303_te_and_ltrs.bed"
DEFAULT_TRNA_BED="$WORKFLOW_ROOT/data/bed/trna_coordinates.bed"


# ==============================================================================
# FUNCTION: usage
# ==============================================================================
#
# Purpose:
#   Print the basic command-line syntax for the workflow.
#
# Arguments:
#   None.
#
# Output:
#   Human-readable usage information written to standard output.
#
# Return status:
#   0
#
usage() {
    echo "Usage: bash src/genome_browser_workflow/genome_browser_workflow_script.sh [BAM] [output_prefix] [mod_threshold]"
    echo
    echo "BAM is resolved by filename under data/bam."
    echo "If arguments are omitted, the workflow prompts for inputs before submitting a SLURM job."
}

warn_unexpected_input() {
    local received="$1"
    local expected="$2"

    echo "[WARN] Unexpected input: ${received:-<blank>}"
    echo "[WARN] Expected input: $expected"
}


# ==============================================================================
# FUNCTION: absolute_existing_file
# ==============================================================================
#
# Purpose:
#   Convert a file path into an absolute path.
#
# Arguments:
#   $1
#       Path to an existing file.
#
# Output:
#   The absolute file path is written to standard output.
#
# Return status:
#   0 when dirname can enter the containing directory.
#   Nonzero when the containing directory cannot be accessed.
#
# Notes:
#   This function does not independently test whether the file exists. Calling
#   functions verify file existence before using it.
#
absolute_existing_file() {
    local path="$1"
    local dir
    local file

    # Resolve the containing directory to an absolute path.
    dir="$(cd "$(dirname "$path")" && pwd)"

    # Preserve the original filename.
    file="$(basename "$path")"

    # Print the reconstructed absolute path.
    printf '%s/%s\n' "$dir" "$file"
}


# ==============================================================================
# FUNCTION: resolve_bam_file
# ==============================================================================
#
# Purpose:
#   Resolve a BAM filename under the workflow's data/bam directory.
#
# Arguments:
#   $1
#       BAM filename or path supplied by the user.
#
# Behavior:
#   Only the basename of the supplied value is used. This intentionally limits
#   BAM selection to files stored directly under:
#
#       workflow_root/data/bam
#
# Example:
#   Input:
#       /some/other/path/sample.bam
#
#   Resolved candidate:
#       workflow_root/data/bam/sample.bam
#
# Output:
#   Absolute BAM path when the file exists.
#
# Return status:
#   0 when the BAM is found.
#   1 when the BAM is not found.
#
resolve_bam_file() {
    local bam_input="$1"
    local bam_filename
    local bam_path

    # Remove any directory components supplied by the user.
    bam_filename="$(basename "$bam_input")"

    # Construct the expected BAM location.
    bam_path="$BAM_DIR/$bam_filename"

    if [[ -f "$bam_path" ]]; then
        absolute_existing_file "$bam_path"
        return 0
    fi

    return 1
}


# ==============================================================================
# FUNCTION: resolve_reference_file
# ==============================================================================
#
# Purpose:
#   Locate and return an absolute path for the requested reference FASTA.
#
# Arguments:
#   $1
#       Reference filename or path supplied by the user.
#
# Search order:
#   1. Use the supplied value directly when it is an existing file.
#   2. Check the path relative to the workflow root.
#   3. Check the path relative to workflow_root/data.
#   4. Recursively search workflow_root/data for a matching filename.
#
# Output:
#   Absolute reference-file path when found.
#
# Return status:
#   0 when the reference is found.
#   1 when no matching reference is found.
#
# Notes:
#   The recursive search returns the first matching filename found by find.
#   Therefore, filenames should ideally be unique under workflow_root/data.
#
resolve_reference_file() {
    local reference_input="$1"
    local reference_path

    # First, test the supplied path exactly as entered.
    if [[ -f "$reference_input" ]]; then
        absolute_existing_file "$reference_input"
        return 0
    fi

    # Next, interpret the value as a path relative to the workflow root.
    reference_path="$WORKFLOW_ROOT/$reference_input"
    if [[ -f "$reference_path" ]]; then
        absolute_existing_file "$reference_path"
        return 0
    fi

    # Next, interpret the value as a path relative to workflow_root/data.
    reference_path="$WORKFLOW_ROOT/data/$reference_input"
    if [[ -f "$reference_path" ]]; then
        absolute_existing_file "$reference_path"
        return 0
    fi

    # Finally, search recursively under the data directory by filename.
    reference_path="$(find "$WORKFLOW_ROOT/data" -type f -name "$reference_input" -print -quit)"
    if [[ -n "$reference_path" && -f "$reference_path" ]]; then
        absolute_existing_file "$reference_path"
        return 0
    fi

    return 1
}


# ==============================================================================
# FUNCTION: normalize_phase_label
# ==============================================================================
#
# Purpose:
#   Convert several accepted cell-cycle phase inputs into the exact label used
#   in genome-browser plot titles.
#
# Arguments:
#   $1
#       User-provided phase value.
#
# Accepted M-phase values:
#   m
#   mitosis
#   mphase
#   m-phase
#
# Accepted S-phase values:
#   s
#   sphase
#   s-phase
#
# Output:
#   "Mitosis" or "S Phase"
#
# Return status:
#   0 for a recognized phase.
#   1 for an unrecognized phase.
#
# Notes:
#   Matching is case-insensitive because ${phase_input,,} converts the input
#   to lowercase before the case statement is evaluated.
#
normalize_phase_label() {
    local phase_input="$1"

    case "${phase_input,,}" in
        m|mitosis|mphase|m-phase)
            printf '%s\n' "Mitosis"
            ;;
        s|sphase|s-phase)
            printf '%s\n' "S Phase"
            ;;
        *)
            return 1
            ;;
    esac
}


# ==============================================================================
# FUNCTION: normalize_plot_mode
# ==============================================================================
#
# Purpose:
#   Convert several accepted plot-mode inputs into one of the two standardized
#   mode names used by the workflow.
#
# Arguments:
#   $1
#       User-provided plotting mode.
#
# Accepted smoothed values:
#   smoothed
#   smooth
#   s
#   1
#
# Accepted unsmoothed values:
#   unsmoothed
#   unsmooth
#   raw
#   u
#   2
#
# Output:
#   "smoothed" or "unsmoothed"
#
# Return status:
#   0 for a recognized plotting mode.
#   1 for an unrecognized plotting mode.
#
normalize_plot_mode() {
    local mode_input="$1"

    case "${mode_input,,}" in
        smoothed|smooth|s|1)
            printf '%s\n' "smoothed"
            ;;
        unsmoothed|unsmooth|raw|u|2)
            printf '%s\n' "unsmoothed"
            ;;
        *)
            return 1
            ;;
    esac
}


# ==============================================================================
# FUNCTION: normalize_mod_threshold
# ==============================================================================
#
# Purpose:
#   Validate and standardize the BrdU modification threshold passed to Modkit.
#
# Arguments:
#   $1
#       Threshold supplied by the user.
#
# Accepted examples:
#   0
#   0.5
#   .5
#   0.6
#   1
#   1.0
#   B:0.5
#   b:0.5
#
# Validation rules:
#   - The numeric value must be between 0 and 1, inclusive.
#   - Values greater than 1 are rejected.
#   - Negative values are rejected.
#   - Non-numeric text is rejected.
#
# Output:
#   A normalized Modkit threshold in the following format:
#
#       B:<value>
#
# Example:
#   Input:
#       .6
#
#   Output:
#       B:0.6
#
# Return status:
#   0 for a valid threshold.
#   1 for an invalid threshold.
#
# Implementation details:
#   - A leading B: or b: is removed before validation.
#   - A value beginning with a decimal point receives a leading zero.
#   - awk performs a final numeric range test and normalized formatting.
#
normalize_mod_threshold() {
    local threshold_input="$1"

    # Remove an uppercase Modkit base-prefix when present.
    local threshold_value="${threshold_input#B:}"

    # Remove a lowercase Modkit base-prefix when present.
    threshold_value="${threshold_value#b:}"

    # Convert values such as ".5" into "0.5".
    if [[ "$threshold_value" =~ ^\.[0-9]+$ ]]; then
        threshold_value="0$threshold_value"
    fi

    # Require a valid numeric representation from 0 through 1.
    if ! [[ "$threshold_value" =~ ^(0(\.[0-9]+)?|1(\.0+)?)$ ]]; then
        return 1
    fi

    # Perform a final numeric check and print the standardized Modkit syntax.
    awk -v value="$threshold_value" 'BEGIN {
        if (value < 0 || value > 1) {
            exit 1
        }
        printf "B:%.12g\n", value + 0
    }'
}


# ==============================================================================
# FUNCTION: mod_threshold_suffix
# ==============================================================================
#
# Purpose:
#   Generate a filename suffix that identifies a nondefault modification
#   threshold.
#
# Arguments:
#   $1
#       A normalized Modkit threshold, such as B:0.5 or B:0.6.
#
# Output behavior:
#   Threshold 0.5:
#       No suffix is printed.
#
#   Threshold 0.6:
#       _06
#
#   Threshold 0.25:
#       _025
#
# Rationale:
#   The default threshold is 0.5, so existing output naming is preserved when
#   the default is used. Nondefault thresholds are included in output names so
#   files produced with different thresholds do not overwrite one another.
#
# Return status:
#   0
#
# Notes:
#   The decimal point is removed from the formatted threshold value.
#
mod_threshold_suffix() {
    local mod_threshold="$1"

    # Remove the Modkit B: prefix before formatting the filename suffix.
    local threshold_value="${mod_threshold#B:}"

    awk -v value="$threshold_value" 'BEGIN {
        # The default threshold does not receive an added suffix.
        if (value + 0 == 0.5) {
            exit 0
        }

        # Format the number and remove the decimal point.
        formatted = sprintf("%.12g", value + 0)
        gsub(/\./, "", formatted)

        # Prefix the threshold identifier with an underscore.
        printf "_%s\n", formatted
    }'
}


# ==============================================================================
# FUNCTION: require_existing_dir
# ==============================================================================
#
# Purpose:
#   Verify that a required directory exists before the compute job proceeds.
#
# Arguments:
#   $1
#       Human-readable directory name used in the error message.
#
#   $2
#       Directory path to validate.
#
# Behavior:
#   If the directory does not exist, the function prints an error and exits the
#   entire workflow.
#
# Return status:
#   Does not return when validation fails.
#   Returns 0 implicitly when the directory exists.
#
require_existing_dir() {
    local name="$1"
    local path="$2"

    if [[ ! -d "$path" ]]; then
        echo "[ERROR] $name does not exist: $path"
        echo "[ERROR] Confirm the project path is mounted on the compute node."
        exit 1
    fi
}


# ==============================================================================
# FUNCTION: load_genome_browser_modules
# ==============================================================================
#
# Purpose:
#   Load the software modules required by the genome-browser workflow.
#
# Modules:
#   - Python 3.13.7, with a fallback to the default Python module.
#   - Samtools.
#   - Modkit.
#
# Behavior when the module command is unavailable:
#   The function prints warnings and continues. In that situation, python,
#   samtools, and modkit must already be available through PATH.
#
# Notes:
#   The expression:
#
#       module load python/3.13.7 || module load python
#
#   first attempts to load the specific Python version and falls back to the
#   cluster's generic Python module if that version cannot be loaded.
#
load_genome_browser_modules() {
    if command -v module >/dev/null 2>&1; then
        module load python/3.13.7 || module load python
        module load samtools
        module load modkit
    else
        echo "[WARN] Environment modules are not available in this shell."
        echo "[WARN] Continuing with samtools, modkit, and python from PATH."
    fi
}


# ==============================================================================
# FUNCTION: run_job
# ==============================================================================
#
# Purpose:
#   Execute the computational portion of the genome-browser workflow inside the
#   submitted SLURM job.
#
# This function:
#   1. Parses the named arguments passed by submit_workflow().
#   2. Normalizes the Modkit threshold.
#   3. Configures workflow logging.
#   4. Validates required input paths.
#   5. Selects the smoothed or unsmoothed Python plotting script.
#   6. Loads required software modules.
#   7. Creates or reuses a Python virtual environment.
#   8. Installs Python requirements.
#   9. Runs BrdU bedgraph extraction.
#  10. Validates the positive and negative bedgraph outputs.
#  11. Adds optional genomic feature tracks when available.
#  12. Generates the genome-browser plots.
#
# Expected internal arguments:
#   --workflow-root PATH
#   --bam PATH
#   --ref PATH
#   --output-prefix PREFIX
#   --mod-threshold VALUE
#   --plot-mode smoothed|unsmoothed
#   --phase-label LABEL
#
# Notes:
#   The first argument is expected to be --run-job. The initial shift removes
#   that mode-selection argument before named argument parsing begins.
#
run_job() {
    # --------------------------------------------------------------------------
    # Job argument defaults and local variables
    # --------------------------------------------------------------------------

    local bam_path=""
    local reference=""
    local output_prefix=""

    # The default threshold is represented in Modkit syntax.
    local mod_threshold="B:0.5"

    local raw_mod_threshold
    local threshold_suffix
    local plot_mode=""
    local phase_label=""
    local generation_script
    local mode_results_dir
    local venv_dir
    local positive_output
    local negative_output
    local workflow_log

    # Remove the leading --run-job argument.
    shift

    # --------------------------------------------------------------------------
    # Parse named arguments passed to the compute job
    # --------------------------------------------------------------------------
    #
    # Each recognized option consumes itself and its following value with:
    #
    #   shift 2
    #
    # --workflow-root was already used near the top of the script to establish
    # WORKFLOW_ROOT, so its value does not need to be assigned again here.
    #
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --workflow-root)
                shift 2
                ;;
            --bam)
                bam_path="$2"
                shift 2
                ;;
            --ref)
                reference="$2"
                shift 2
                ;;
            --output-prefix)
                output_prefix="$2"
                shift 2
                ;;
            --mod-threshold)
                mod_threshold="$2"
                shift 2
                ;;
            --plot-mode)
                plot_mode="$2"
                shift 2
                ;;
            --phase-label)
                phase_label="$2"
                shift 2
                ;;
            *)
                echo "[ERROR] Unknown job argument: $1"
                exit 1
                ;;
        esac
    done

    # --------------------------------------------------------------------------
    # Normalize and validate the modification threshold
    # --------------------------------------------------------------------------

    # Preserve the original value for a useful error message.
    raw_mod_threshold="$mod_threshold"

    if ! mod_threshold="$(normalize_mod_threshold "$raw_mod_threshold")"; then
        echo "[ERROR] Invalid mod threshold in job: $raw_mod_threshold"
        echo "[ERROR] Expected a number from 0 to 1, for example 0.5 or 0.6."
        exit 1
    fi

    # Generate an output suffix for nondefault thresholds.
    threshold_suffix="$(mod_threshold_suffix "$mod_threshold")"

    # Add the suffix only when one is required and the prefix does not already
    # end with that threshold identifier.
    if [[ -n "$threshold_suffix" && "$output_prefix" != *"$threshold_suffix" ]]; then
        output_prefix="${output_prefix}${threshold_suffix}"
    fi

    # --------------------------------------------------------------------------
    # Configure the combined workflow log
    # --------------------------------------------------------------------------
    #
    # Standard output and standard error from this point forward are:
    #
    #   - Displayed in the active SLURM log stream.
    #   - Appended to a separate workflow log through tee.
    #
    # When SLURM_JOB_ID is unavailable, such as during manual execution, the
    # identifier "manual" is used.
    #
    mkdir -p "$LOG_DIR"
    workflow_log="$LOG_DIR/${output_prefix}.${SLURM_JOB_ID:-manual}.workflow.log"

    exec > >(tee -a "$workflow_log") 2>&1

    echo "[INFO] Genome browser workflow started: $(date)"
    echo "[INFO] Workflow log: $workflow_log"

    # --------------------------------------------------------------------------
    # Validate required directories
    # --------------------------------------------------------------------------

    require_existing_dir "Workflow root" "$WORKFLOW_ROOT"
    require_existing_dir "Source directory" "$SRC_DIR"
    require_existing_dir "BAM directory" "$BAM_DIR"

    # Create all workflow output directories when they do not already exist.
    mkdir -p "$BEDGRAPH_DIR" "$SORTED_BAM_DIR" "$RESULTS_DIR" "$LOG_DIR"

    # --------------------------------------------------------------------------
    # Validate required input files on the compute node
    # --------------------------------------------------------------------------
    #
    # This second validation is important because a path may have existed in
    # the submission shell but may not be mounted or visible on the compute node.
    #
    if [[ ! -f "$bam_path" ]]; then
        echo "[ERROR] BAM file not found on compute node: $bam_path"
        exit 1
    fi

    if [[ ! -f "$reference" ]]; then
        echo "[ERROR] Reference FASTA not found on compute node: $reference"
        exit 1
    fi

    # --------------------------------------------------------------------------
    # Select the plotting script and results directory
    # --------------------------------------------------------------------------

    case "$plot_mode" in
        smoothed)
            generation_script="$SRC_DIR/genomic_browser_generation.py"
            mode_results_dir="$RESULTS_DIR/smoothed"
            ;;
        unsmoothed)
            generation_script="$SRC_DIR/genomic_browser_generation_unsmoothed.py"
            mode_results_dir="$RESULTS_DIR/unsmoothed"
            ;;
        *)
            echo "[ERROR] Invalid plot mode in job: $plot_mode"
            exit 1
            ;;
    esac

    mkdir -p "$mode_results_dir"

    # --------------------------------------------------------------------------
    # Load cluster software
    # --------------------------------------------------------------------------

    load_genome_browser_modules

    # --------------------------------------------------------------------------
    # Create or reuse the Python virtual environment
    # --------------------------------------------------------------------------
    #
    # The environment is stored inside the source directory:
    #
    #   src/genome_browser_workflow/.genome_browser_env
    #
    # A new environment is created only when its Python executable is missing.
    #
    venv_dir="$SRC_DIR/.genome_browser_env"

    if [[ ! -x "$venv_dir/bin/python" ]]; then
        echo "[INFO] Creating virtual environment: $venv_dir"
        python3 -m venv --clear "$venv_dir"
    fi

    echo "[INFO] Activating virtual environment."
    source "$venv_dir/bin/activate"

    # Install or verify all packages defined in requirements.txt.
    #
    # Running pip install on each workflow execution allows missing packages or
    # updated requirement versions to be installed without manually rebuilding
    # the environment.
    echo "[INFO] Installing Python requirements."
    python -m pip install -r "$SRC_DIR/requirements.txt" --quiet

    # --------------------------------------------------------------------------
    # Define expected strand-specific bedgraph outputs
    # --------------------------------------------------------------------------

    positive_output="$BEDGRAPH_DIR/$output_prefix.positive.bedgraph"
    negative_output="$BEDGRAPH_DIR/$output_prefix.negative.bedgraph"

    # --------------------------------------------------------------------------
    # Extract raw BrdU data from the BAM file
    # --------------------------------------------------------------------------
    #
    # raw_data_extraction_on_bam.py is responsible for running the underlying
    # BrdU/Modkit extraction and producing separate bedgraph files for the
    # positive and negative strands.
    #
    # Arguments:
    #   Positional BAM path
    #       Input modBAM file.
    #
    #   --ref
    #       Reference FASTA used for Modkit pileup.
    #
    #   --output-prefix
    #       Prefix shared by bedgraph and log outputs.
    #
    #   --threads
    #       Number of threads allocated by SLURM. Falls back to 12 when the
    #       SLURM_CPUS_PER_TASK environment variable is unavailable.
    #
    #   --mod-threshold
    #       Normalized BrdU probability threshold, such as B:0.5.
    #
    #   --bedgraph-dir
    #       Destination directory for bedgraph outputs.
    #
    echo "[INFO] Running raw BrdU bedgraph extraction."

    python "$SRC_DIR/raw_data_extraction_on_bam.py" \
        "$bam_path" \
        --ref "$reference" \
        --output-prefix "$output_prefix" \
        --threads "${SLURM_CPUS_PER_TASK:-12}" \
        --mod-threshold "$mod_threshold" \
        --bedgraph-dir "$BEDGRAPH_DIR"

    # --------------------------------------------------------------------------
    # Validate the bedgraph outputs
    # --------------------------------------------------------------------------
    #
    # -s requires the file to exist and have a size greater than zero.
    #
    # The workflow stops before plotting if either strand-specific file is
    # absent or empty.
    #
    if [[ ! -s "$positive_output" ]]; then
        echo "[ERROR] Positive bedgraph file is empty or missing: $positive_output"
        deactivate
        exit 1
    fi

    if [[ ! -s "$negative_output" ]]; then
        echo "[ERROR] Negative bedgraph file is empty or missing: $negative_output"
        deactivate
        exit 1
    fi

    echo "[INFO] Positive strand bedgraph: $positive_output"
    echo "[INFO] Negative strand bedgraph: $negative_output"
    echo "[INFO] Modkit log: $BEDGRAPH_DIR/$output_prefix.modkit.log"

    # --------------------------------------------------------------------------
    # Build the plotting-script argument array
    # --------------------------------------------------------------------------
    #
    # An array is used so each path remains a separate argument, including paths
    # that might contain spaces.
    #
    generation_args=(
        --positive-bedgraph "$positive_output"
        --negative-bedgraph "$negative_output"
        --output-dir "$mode_results_dir"
        --prefix "$output_prefix"
        --phase-label "$phase_label"
    )

    # --------------------------------------------------------------------------
    # Add optional genomic feature tracks
    # --------------------------------------------------------------------------
    #
    # Missing optional BED files do not stop the workflow. Their corresponding
    # command-line arguments are simply omitted.
    #

    # Add the G-quadruplex motif BED file when available.
    if [[ -f "$DEFAULT_G4_BED" ]]; then
        generation_args+=(--g4-bed "$DEFAULT_G4_BED")
    fi

    # Add the transposable element and LTR BED file when available.
    if [[ -f "$DEFAULT_TE_BED" ]]; then
        generation_args+=(--te-bed "$DEFAULT_TE_BED")
    fi

    # Add the tRNA-coordinate BED file when available.
    if [[ -f "$DEFAULT_TRNA_BED" ]]; then
        generation_args+=(--trna-bed "$DEFAULT_TRNA_BED")
    fi

    # --------------------------------------------------------------------------
    # Generate the genome-browser plots
    # --------------------------------------------------------------------------

    echo "[INFO] Generating $plot_mode genome browser plots."
    python "$generation_script" "${generation_args[@]}"

    # Exit the workflow-specific Python environment after plotting is complete.
    deactivate

    echo "[INFO] Genome browser workflow complete."
    echo "[INFO] Results directory: $mode_results_dir"
    echo "[INFO] Genome browser workflow finished: $(date)"
}


# ==============================================================================
# FUNCTION: submit_workflow
# ==============================================================================
#
# Purpose:
#   Collect and validate user input, then submit the computational workflow as a
#   SLURM job.
#
# Positional arguments:
#   $1
#       Optional BAM filename.
#
#   $2
#       Optional output prefix.
#
#   $3
#       Optional Modkit modification threshold.
#
# Interactive prompts:
#   - BAM filename, when not supplied as $1.
#   - Reference FASTA.
#   - Modification threshold, when not supplied as $3.
#   - Cell-cycle phase.
#   - Smoothed or unsmoothed plotting mode.
#
# Output:
#   Prints the selected configuration, submitted job ID, and expected output
#   locations.
#
# Notes:
#   This function submits the same script back to SLURM using --run-job. This
#   separates lightweight user interaction from the compute-intensive work.
#
submit_workflow() {
    # --------------------------------------------------------------------------
    # Read optional positional arguments
    # --------------------------------------------------------------------------

    local bam_input="${1:-}"
    local output_prefix="${2:-}"
    local mod_threshold_input="${3:-}"

    # --------------------------------------------------------------------------
    # Declare remaining submission variables
    # --------------------------------------------------------------------------

    local ref_input=""
    local bam_path=""
    local reference=""
    local bam_name
    local phase_input
    local phase_label
    local plot_mode_input
    local plot_mode
    local mod_threshold
    local threshold_suffix
    local effective_output_prefix
    local mode_results_dir
    local job_id
    local log_job_id

    # Create the expected input/output directory structure when directories are
    # missing. Creating BAM_DIR allows the user to see the expected location
    # even when no BAM files have been copied into it yet.
    mkdir -p "$BAM_DIR" "$BEDGRAPH_DIR" "$SORTED_BAM_DIR" "$RESULTS_DIR" "$LOG_DIR"

    # --------------------------------------------------------------------------
    # Select the input BAM file
    # --------------------------------------------------------------------------

    while true; do
        if [[ -z "$bam_input" ]]; then
            echo "[INFO] Available BAM files in $BAM_DIR:"

            # List BAM files directly inside data/bam, sorted alphabetically.
            find "$BAM_DIR" -maxdepth 1 -type f -name "*.bam" -printf "  %f\n" | sort

            echo
            read -r -p "Enter the BAM filename from data/bam: " bam_input
        fi

        if [[ -z "$bam_input" ]]; then
            warn_unexpected_input "$bam_input" "a BAM filename under $BAM_DIR, for example sample.bam"
            continue
        fi

        if bam_path="$(resolve_bam_file "$bam_input")"; then
            break
        fi

        warn_unexpected_input "$bam_input" "an existing BAM filename under $BAM_DIR, for example sample.bam"
        bam_input=""
    done

    # --------------------------------------------------------------------------
    # Determine the base output prefix
    # --------------------------------------------------------------------------
    #
    # When the user does not provide an output prefix, the .bam extension is
    # removed from the selected BAM filename.
    #
    # Example:
    #   sample.sorted.bam
    #
    # Default output prefix:
    #   sample.sorted
    #
    if [[ -z "$output_prefix" ]]; then
        bam_name="$(basename "$bam_path")"
        output_prefix="${bam_name%.bam}"
    fi

    # --------------------------------------------------------------------------
    # Select and validate the reference FASTA
    # --------------------------------------------------------------------------
    #
    # Pressing Enter accepts DEFAULT_REF.
    #
    while true; do
        read -r -p "Reference FASTA for modkit pileup [$DEFAULT_REF]: " ref_input
        ref_input="${ref_input:-$DEFAULT_REF}"

        if reference="$(resolve_reference_file "$ref_input")"; then
            break
        fi

        warn_unexpected_input "$ref_input" "an existing FASTA path, for example $DEFAULT_REF"
    done

    # --------------------------------------------------------------------------
    # Select and validate the Modkit threshold
    # --------------------------------------------------------------------------
    #
    # The default threshold is 0.5. This can also be supplied as the third
    # positional argument to avoid the threshold prompt.
    #
    while true; do
        if [[ -z "$mod_threshold_input" ]]; then
            read -r -p "Modkit mod threshold [0.5]: " mod_threshold_input
            mod_threshold_input="${mod_threshold_input:-0.5}"
        fi

        if mod_threshold="$(normalize_mod_threshold "$mod_threshold_input")"; then
            break
        fi

        warn_unexpected_input "$mod_threshold_input" "a number from 0 to 1, for example 0.5 or 0.6"
        mod_threshold_input=""
    done

    # Add a filename suffix for nondefault thresholds.
    threshold_suffix="$(mod_threshold_suffix "$mod_threshold")"
    effective_output_prefix="${output_prefix}${threshold_suffix}"

    # --------------------------------------------------------------------------
    # Select the cell-cycle phase label
    # --------------------------------------------------------------------------
    #
    # This value controls the label shown in generated genome-browser titles.
    #
    echo
    echo "Choose cell-cycle phase for genome browser plot titles:"
    echo "  M) Mitosis"
    echo "  S) S Phase"
    while true; do
        read -r -p "Enter M or S: " phase_input

        if phase_label="$(normalize_phase_label "$phase_input")"; then
            break
        fi

        warn_unexpected_input "$phase_input" "M or S, for example M"
    done

    # --------------------------------------------------------------------------
    # Select the plot-generation mode
    # --------------------------------------------------------------------------
    #
    # The default is smoothed when the user presses Enter.
    #
    echo
    echo "Choose genome browser output mode:"
    echo "  1) smoothed"
    echo "  2) unsmoothed"
    while true; do
        read -r -p "Enter smoothed or unsmoothed [smoothed]: " plot_mode_input
        plot_mode_input="${plot_mode_input:-smoothed}"

        if plot_mode="$(normalize_plot_mode "$plot_mode_input")"; then
            break
        fi

        warn_unexpected_input "$plot_mode_input" "smoothed or unsmoothed, for example smoothed"
    done

    # Determine where the selected plot type will be written.
    mode_results_dir="$RESULTS_DIR/$plot_mode"
    mkdir -p "$mode_results_dir"

    # Load modules in the submission shell so the sbatch command and expected
    # software environment are available.
    load_genome_browser_modules

    # --------------------------------------------------------------------------
    # Display the final configuration before submission
    # --------------------------------------------------------------------------

    echo
    echo "[INFO] Workflow root: $WORKFLOW_ROOT"
    echo "[INFO] BAM: $bam_path"
    echo "[INFO] Reference FASTA: $reference"
    echo "[INFO] Output prefix: $effective_output_prefix"
    echo "[INFO] Mod threshold: $mod_threshold"
    echo "[INFO] Phase label: $phase_label"
    echo "[INFO] Plot mode: $plot_mode"
    echo
    echo "Submitting genome browser workflow to SLURM..."

    # --------------------------------------------------------------------------
    # Submit this script back to SLURM in compute-job mode
    # --------------------------------------------------------------------------
    #
    # sbatch options:
    #
    #   --parsable
    #       Return the job ID in a machine-readable format.
    #
    #   --job-name
    #       Assign a job name containing the effective output prefix.
    #
    #   --chdir
    #       Run the job from the workflow root.
    #
    #   --output
    #       Save standard output to a SLURM log containing the job ID.
    #
    #   --error
    #       Save standard error to a separate SLURM error log.
    #
    #   --export
    #       Export the current environment and explicitly provide the resolved
    #       workflow root to the compute job.
    #
    # Script arguments:
    #
    #   --run-job
    #       Select the run_job() execution path.
    #
    #   Remaining named arguments
    #       Pass all validated paths and selections to the compute job.
    #
    job_id="$(sbatch --parsable \
        --job-name="genome_browser_${effective_output_prefix}" \
        --chdir="$WORKFLOW_ROOT" \
        --output="$LOG_DIR/${effective_output_prefix}.%j.slurm.log" \
        --error="$LOG_DIR/${effective_output_prefix}.%j.slurm.err" \
        --export=ALL,GENOME_BROWSER_WORKFLOW_ROOT="$WORKFLOW_ROOT" \
        "$SCRIPT_PATH" \
        --run-job \
        --workflow-root "$WORKFLOW_ROOT" \
        --bam "$bam_path" \
        --ref "$reference" \
        --output-prefix "$effective_output_prefix" \
        --mod-threshold "$mod_threshold" \
        --plot-mode "$plot_mode" \
        --phase-label "$phase_label")"

    # Some SLURM configurations may append cluster information after a
    # semicolon. Keep only the numeric job-ID portion for constructing paths.
    log_job_id="${job_id%%;*}"

    # --------------------------------------------------------------------------
    # Report the submitted job and expected outputs
    # --------------------------------------------------------------------------

    echo "[INFO] Submitted SLURM job: $job_id"
    echo "[INFO] Expected positive bedgraph: $BEDGRAPH_DIR/$effective_output_prefix.positive.bedgraph"
    echo "[INFO] Expected negative bedgraph: $BEDGRAPH_DIR/$effective_output_prefix.negative.bedgraph"
    echo "[INFO] Expected results dir:       $mode_results_dir"
    echo "[INFO] SLURM log:                  $LOG_DIR/${effective_output_prefix}.${log_job_id}.slurm.log"
    echo "[INFO] SLURM err:                  $LOG_DIR/${effective_output_prefix}.${log_job_id}.slurm.err"
    echo "[INFO] Workflow log:               $LOG_DIR/${effective_output_prefix}.${log_job_id}.workflow.log"
    echo "[INFO] Modkit log:                 $BEDGRAPH_DIR/$effective_output_prefix.modkit.log"
}


# ==============================================================================
# MAIN SCRIPT ROUTING
# ==============================================================================
#
# The first command-line argument determines which execution path is used:
#
#   -h or --help
#       Print usage information and exit.
#
#   --run-job
#       Execute the computational job. This is normally supplied internally by
#       submit_workflow() when the SLURM job is created.
#
#   Anything else, including no arguments
#       Enter submission mode, collect user input, and submit the SLURM job.
#

# Print help without running or submitting the workflow.
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

# Select compute-job mode or user-facing submission mode.
if [[ "${1:-}" == "--run-job" ]]; then
    run_job "$@"
else
    submit_workflow "$@"
fi
