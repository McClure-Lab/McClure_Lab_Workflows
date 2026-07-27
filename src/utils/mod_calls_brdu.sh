#!/bin/bash

# ==============================================================================
# BrdU-POSITIVE READ EXTRACTION WORKFLOW
# ==============================================================================
#
# Purpose
# -------
# This workflow identifies reads containing at least one passing BrdU modified-
# base call and creates a new BAM containing the complete primary mapped
# alignments for those reads.
#
# BrdU is identified using:
#
#   Canonical base:      T
#   Modification code:   b
#
# A read is selected when Modkit reports at least one BrdU call whose:
#
#   call_code == "b"
#   fail      == "false"
#
# The Modkit --mod-threshold value controls the minimum modification probability
# required for a call to pass.
#
# Main workflow stages
# --------------------
#
# Interactive submission stage:
#
#   1. Determine the workflow root.
#   2. List BAM files under data/bam.
#   3. Prompt for an input BAM.
#   4. Prompt for a BrdU modification threshold.
#   5. Construct descriptive output filenames.
#   6. Prevent accidental overwriting of existing outputs.
#   7. Submit this script to SLURM in --run-analysis mode.
#
# Submitted analysis stage:
#
#   1. Load Modkit and Samtools.
#   2. Validate required programs and input files.
#   3. Create an input BAM index when one is missing.
#   4. Count mapped primary reads in the input BAM.
#   5. Run modkit extract calls.
#   6. Parse the Modkit output header dynamically.
#   7. Collect unique read IDs containing passing BrdU calls.
#   8. Extract complete primary mapped alignments for those read IDs.
#   9. Coordinate-sort the selected alignments.
#  10. Index and validate the final BAM.
#  11. Compare the output read count with the selected read-ID count.
#
# Run command
# -----------
#
# Run from the workflow root:
#
#   bash src/utils/mod_calls_brdu.sh
#
# Optional positional arguments:
#
#   bash src/utils/mod_calls_brdu.sh \
#       sample.bam \
#       0.50
#
# Help:
#
#   bash src/utils/mod_calls_brdu.sh --help
#
# Expected project structure
# --------------------------
#
# workflow_root/
# ├── data/
# │   └── bam/
# ├── logs/
# │   └── mod_calls_brdu/
# └── src/
#     └── utils/
#         └── mod_calls_brdu.sh
#
# Output files
# ------------
#
# The output prefix is based on the input sample name and threshold:
#
#   <sample>.Brdu_positive.threshold_<threshold>
#
# Example at threshold 0.50:
#
#   sample.Brdu_positive.threshold_0p5.bam
#   sample.Brdu_positive.threshold_0p5.bam.bai
#   sample.Brdu_positive.threshold_0p5.read_ids.txt
#
# The threshold is made filename-safe by replacing the decimal point with p.
#
# Primary mapped read definition
# ------------------------------
#
# Samtools flag filter:
#
#   -F 2308
#
# excludes:
#
#   0x4    unmapped alignments
#   0x100  secondary alignments
#   0x800  supplementary alignments
#
# Therefore, this workflow extracts primary mapped alignments only.
#
# Temporary files
# ---------------
#
# Temporary files are created under:
#
#   $SLURM_TMPDIR
#
# when that variable is available. Otherwise, /tmp is used.
#
# The temporary directory is removed automatically when the analysis exits.
#
# ==============================================================================


# ==============================================================================
# SLURM RESOURCE REQUESTS
# ==============================================================================
#
# Requested resources:
#
#   Job name:       BrdU_mod_calls
#   Tasks:          1
#   CPU cores:      8
#   Memory:         128 GB
#   Runtime limit:  12 hours
#
#SBATCH --job-name=BrdU_mod_calls
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00


# ==============================================================================
# BASH SAFETY SETTINGS
# ==============================================================================
#
# -e
#   Exit when an unhandled command fails.
#
# -u
#   Treat unset variables as errors.
#
# -o pipefail
#   Treat a pipeline as failed when any command in that pipeline fails.
#
set -euo pipefail


# ==============================================================================
# STATUS-MESSAGE FUNCTIONS
# ==============================================================================

# Print a formatted informational message.
print_info() {
    printf '[INFO] %s\n' "$1"
}

# Print a formatted warning to standard error.
print_warning() {
    printf '[WARNING] %s\n' "$1" >&2
}

# Print a formatted error to standard error.
print_error() {
    printf '[ERROR] %s\n' "$1" >&2
}


# ==============================================================================
# FUNCTION: determine_workflow_root
# ==============================================================================
#
# Purpose:
#   Determine the workflow project root when no explicit path or environment
#   variable is supplied.
#
# Resolution strategy:
#
#   1. If the current directory contains both data/bam and src/utils, use the
#      current directory.
#
#   2. Otherwise, resolve the directory containing this script and move two
#      directory levels upward.
#
# Output:
#   Workflow-root path written to standard output.
#
determine_workflow_root() {
    local script_directory

    # Prefer the current working directory when it looks like the workflow root.
    if [[ -d "$PWD/data/bam" && -d "$PWD/src/utils" ]]; then
        printf '%s\n' "$PWD"
        return 0
    fi

    # Resolve the directory containing the current script.
    script_directory="$(
        cd "$(dirname "${BASH_SOURCE[0]}")" &&
        pwd
    )"

    # This script is expected under workflow_root/src/utils.
    cd "${script_directory}/../.." && pwd
}


# ==============================================================================
# RESOLVE THE WORKFLOW ROOT
# ==============================================================================
#
# Resolution order:
#
#   1. --workflow-root PATH
#   2. MOD_CALLS_BRDU_WORKFLOW_ROOT environment variable
#   3. determine_workflow_root()
#

WORKFLOW_ROOT_ARG=""

# Search all command-line arguments for --workflow-root.
for ((arg_i = 1; arg_i <= $#; arg_i++)); do
    if [[ "${!arg_i}" == "--workflow-root" ]]; then
        next_arg_i=$((arg_i + 1))
        WORKFLOW_ROOT_ARG="${!next_arg_i:-}"
        break
    fi
done

if [[ -n "$WORKFLOW_ROOT_ARG" ]]; then
    WORKFLOW_ROOT="$(cd "$WORKFLOW_ROOT_ARG" && pwd)"
elif [[ -n "${MOD_CALLS_BRDU_WORKFLOW_ROOT:-}" ]]; then
    WORKFLOW_ROOT="$MOD_CALLS_BRDU_WORKFLOW_ROOT"
else
    WORKFLOW_ROOT="$(determine_workflow_root)"
fi


# ==============================================================================
# WORKFLOW PATHS AND RESOURCE SETTINGS
# ==============================================================================

# Directory containing input and output BAM files.
BAM_DIRECTORY="${WORKFLOW_ROOT}/data/bam"

# Directory containing SLURM output and error logs.
LOG_DIRECTORY="${WORKFLOW_ROOT}/logs/mod_calls_brdu"

# Absolute expected location of this script.
#
# The interactive stage submits this path back to SLURM.
SCRIPT_PATH="${WORKFLOW_ROOT}/src/utils/mod_calls_brdu.sh"

# Memory supplied explicitly during sbatch submission.
JOB_MEMORY="128G"

# Create the log directory when it does not already exist.
mkdir -p "$LOG_DIRECTORY"


# ==============================================================================
# FUNCTION: strip_bam_suffixes
# ==============================================================================
#
# Purpose:
#   Generate a clean sample prefix from a BAM filename.
#
# The function removes:
#
#   - The final .bam extension.
#   - Common BrdU-detection suffixes.
#   - Common sorted/indexed suffixes.
#
# Examples:
#
#   sample.sorted.indexed.BrdU.detect.bam
#       becomes:
#   sample
#
#   sample.BrdU.detect.bam
#       becomes:
#   sample
#
# Arguments:
#   $1
#       BAM filename.
#
# Output:
#   Clean sample prefix.
#
strip_bam_suffixes() {
    local filename="$1"

    # Remove the final BAM extension.
    filename="${filename%.bam}"

    # Remove combined sorted/indexed/BrdU-detect suffixes.
    filename="${filename%.sorted.indexed.BrdU.detect}"
    filename="${filename%.sorted.indexed.Brdu.detect}"
    filename="${filename%.sorted.indexed.brdu.detect}"

    # Remove indexed/BrdU-detect suffixes.
    filename="${filename%.indexed.BrdU.detect}"
    filename="${filename%.indexed.Brdu.detect}"
    filename="${filename%.indexed.brdu.detect}"

    # Remove sorted/BrdU-detect suffixes.
    filename="${filename%.sorted.BrdU.detect}"
    filename="${filename%.sorted.Brdu.detect}"
    filename="${filename%.sorted.brdu.detect}"

    # Remove BrdU-detect suffixes.
    filename="${filename%.BrdU.detect}"
    filename="${filename%.Brdu.detect}"
    filename="${filename%.brdu.detect}"

    # Remove remaining sorted or indexed suffixes.
    filename="${filename%.sorted.indexed}"
    filename="${filename%.indexed}"
    filename="${filename%.sorted}"

    printf '%s\n' "$filename"
}


# ==============================================================================
# FUNCTION: normalize_threshold_for_filename
# ==============================================================================
#
# Purpose:
#   Convert a numeric threshold into a compact filename-safe label.
#
# Processing:
#
#   - Remove trailing zeros.
#   - Remove a trailing decimal point.
#   - Replace the decimal point with p.
#
# Examples:
#
#   0.50  -> 0p5
#   0.60  -> 0p6
#   0.625 -> 0p625
#
normalize_threshold_for_filename() {
    local threshold="$1"

    printf '%s\n' "$threshold" |
        sed -E 's/0+$//; s/\.$//; s/\./p/g'
}


# ==============================================================================
# FUNCTION: validate_threshold
# ==============================================================================
#
# Purpose:
#   Validate a numeric modification threshold between 0 and 1 inclusive.
#
# Accepted examples:
#
#   0
#   0.5
#   0.50
#   0.75
#   1
#   1.0
#
# Return status:
#
#   0 for a valid threshold.
#   1 for an invalid threshold.
#
validate_threshold() {
    local threshold="$1"

    awk -v value="$threshold" 'BEGIN {
        if (value ~ /^[0-9]+([.][0-9]+)?$/ && value >= 0 && value <= 1) {
            exit 0
        }

        exit 1
    }'
}


# ==============================================================================
# FUNCTION: check_required_command
# ==============================================================================
#
# Purpose:
#   Verify that a required executable is available through PATH.
#
# Arguments:
#   $1
#       Command name.
#
# Behavior:
#   The workflow exits when the command cannot be found.
#
check_required_command() {
    local command_name="$1"

    if ! command -v "$command_name" >/dev/null 2>&1; then
        print_error "Required command is unavailable: ${command_name}"
        exit 1
    fi
}


# ==============================================================================
# FUNCTION: load_required_modules
# ==============================================================================
#
# Purpose:
#   Load Modkit and Samtools through the cluster module system.
#
# Behavior:
#
#   - Sources /etc/profile.d/modules.sh when available.
#   - Loads the Modkit module.
#   - Loads the Samtools module.
#   - Falls back to commands already available through PATH when the module
#     command is unavailable.
#
load_required_modules() {
    if [[ -f /etc/profile.d/modules.sh ]]; then
        # shellcheck source=/dev/null
        source /etc/profile.d/modules.sh
    fi

    if command -v module >/dev/null 2>&1; then
        module load modkit
        module load samtools
    else
        print_warning "Environment modules are not available in this shell."
        print_warning "Continuing with modkit and samtools from PATH."
    fi
}


# ==============================================================================
# FUNCTION: usage
# ==============================================================================
#
# Purpose:
#   Print command-line usage information.
#
# Optional positional arguments:
#
#   $1
#       BAM filename.
#
#   $2
#       BrdU modification threshold.
#
usage() {
    echo "Usage: bash src/utils/mod_calls_brdu.sh [BAM] [mod_threshold]"
    echo
    echo "BAM is resolved by filename under data/bam."
    echo "If arguments are omitted, the script prompts for inputs before submitting a SLURM job."
    echo "mod_threshold defaults to 0.50."
}


# ==============================================================================
# FUNCTION: submit_workflow
# ==============================================================================
#
# Purpose:
#   Run the interactive launcher and submit the analysis to SLURM.
#
# Positional arguments:
#
#   $1
#       Optional BAM filename.
#
#   $2
#       Optional BrdU modification threshold.
#
# Major steps:
#
#   1. Validate the workflow directories and script path.
#   2. Confirm sbatch is available.
#   3. List BAM files when one was not supplied.
#   4. Prompt for a numeric BAM selection.
#   5. Validate the modification threshold.
#   6. Construct output filenames.
#   7. Prevent overwriting existing outputs.
#   8. Submit the analysis to SLURM.
#
submit_workflow() {
    # --------------------------------------------------------------------------
    # Optional positional arguments
    # --------------------------------------------------------------------------

    local bam_input="${1:-}"
    local mod_threshold="${2:-}"

    # --------------------------------------------------------------------------
    # Interactive-stage variables
    # --------------------------------------------------------------------------

    local bam_files=()
    local bam_selection
    local selected_bam
    local sample_prefix
    local threshold_tag
    local output_prefix
    local output_bam
    local output_index
    local output_read_ids
    local existing_output_found
    local slurm_log_prefix
    local submission_output
    local job_id

    # --------------------------------------------------------------------------
    # Print the launcher heading
    # --------------------------------------------------------------------------

    printf '\n'
    printf '%s\n' '============================================================'
    printf '%s\n' 'Extract reads with passing BrdU modification calls'
    printf '%s\n' '============================================================'
    printf '\n'

    print_info "Workflow root: ${WORKFLOW_ROOT}"
    print_info "BAM directory: ${BAM_DIRECTORY}"

    # --------------------------------------------------------------------------
    # Validate required paths and SLURM access
    # --------------------------------------------------------------------------

    if [[ ! -d "$BAM_DIRECTORY" ]]; then
        print_error "BAM directory does not exist:"
        print_error "$BAM_DIRECTORY"
        exit 1
    fi

    if [[ ! -f "$SCRIPT_PATH" ]]; then
        print_error "Could not find this script:"
        print_error "$SCRIPT_PATH"
        exit 1
    fi

    if ! command -v sbatch >/dev/null 2>&1; then
        print_error "The sbatch command is unavailable."
        exit 1
    fi

    # --------------------------------------------------------------------------
    # Select the input BAM
    # --------------------------------------------------------------------------

    if [[ -z "$bam_input" ]]; then
        # Read available BAM filenames into an indexed Bash array.
        mapfile -t bam_files < <(
            find "$BAM_DIRECTORY" \
                -maxdepth 1 \
                -type f \
                -name '*.bam' \
                -printf '%f\n' |
            sort
        )

        if (( ${#bam_files[@]} == 0 )); then
            print_error "No BAM files were found under:"
            print_error "$BAM_DIRECTORY"
            exit 1
        fi

        printf '\n'
        printf '%s\n' 'Available BAM files:'
        printf '\n'

        # Display a one-based numeric selection menu.
        for index in "${!bam_files[@]}"; do
            printf '  %3d) %s\n' \
                "$((index + 1))" \
                "${bam_files[$index]}"
        done

        printf '\n'

        # Repeat until the user chooses a valid BAM number.
        while true; do
            read -r -p "Select a BAM file by number: " bam_selection

            if [[ ! "$bam_selection" =~ ^[0-9]+$ ]]; then
                print_warning "Enter a numeric selection."
                continue
            fi

            if (( bam_selection < 1 || bam_selection > ${#bam_files[@]} )); then
                print_warning "Selection must be between 1 and ${#bam_files[@]}."
                continue
            fi

            bam_input="${bam_files[$((bam_selection - 1))]}"
            break
        done
    fi

    # Use only the filename component.
    selected_bam="$(basename "$bam_input")"

    # BAM inputs are resolved under BAM_DIRECTORY.
    if [[ ! -f "${BAM_DIRECTORY}/${selected_bam}" ]]; then
        print_error "Input BAM does not exist:"
        print_error "${BAM_DIRECTORY}/${selected_bam}"
        exit 1
    fi

    # --------------------------------------------------------------------------
    # Select and validate the BrdU modification threshold
    # --------------------------------------------------------------------------

    if [[ -z "$mod_threshold" ]]; then
        printf '\n'

        while true; do
            read -r -p "Enter the BrdU modification threshold [0.50]: " mod_threshold

            # Pressing Enter accepts the default.
            mod_threshold="${mod_threshold:-0.50}"

            if validate_threshold "$mod_threshold"; then
                break
            fi

            print_warning "The modification threshold must be a number between 0 and 1."
        done
    elif ! validate_threshold "$mod_threshold"; then
        print_error "Invalid modification threshold: ${mod_threshold}"
        exit 1
    fi

    # --------------------------------------------------------------------------
    # Construct output paths
    # --------------------------------------------------------------------------

    sample_prefix="$(strip_bam_suffixes "$selected_bam")"

    threshold_tag="$(normalize_threshold_for_filename "$mod_threshold")"

    output_prefix="${BAM_DIRECTORY}/${sample_prefix}.Brdu_positive.threshold_${threshold_tag}"

    output_bam="${output_prefix}.bam"

    output_index="${output_bam}.bai"

    output_read_ids="${output_prefix}.read_ids.txt"

    # --------------------------------------------------------------------------
    # Display the submission summary
    # --------------------------------------------------------------------------

    printf '\n'
    printf '%s\n' '------------------------------------------------------------'
    printf '%s\n' 'Submission summary'
    printf '%s\n' '------------------------------------------------------------'
    printf 'Input BAM:       %s\n' "$selected_bam"
    printf 'BrdU code:       b\n'
    printf 'Mod threshold:   %s\n' "$mod_threshold"
    printf 'Output BAM:      %s\n' "$output_bam"
    printf 'Read-ID file:    %s\n' "$output_read_ids"
    printf '%s\n' '------------------------------------------------------------'
    printf '\n'

    # --------------------------------------------------------------------------
    # Prevent accidental output replacement
    # --------------------------------------------------------------------------
    #
    # Unlike the merge utility, this script does not offer an overwrite prompt.
    # Existing outputs must be renamed or removed before rerunning.
    #

    existing_output_found=0

    for output_file in "$output_bam" "$output_index" "$output_read_ids"; do
        if [[ -e "$output_file" ]]; then
            print_error "Output already exists: ${output_file}"
            existing_output_found=1
        fi
    done

    if (( existing_output_found == 1 )); then
        printf '\n' >&2
        print_error "Rename or remove the existing output before rerunning."
        exit 1
    fi

    # --------------------------------------------------------------------------
    # Construct the SLURM log prefix
    # --------------------------------------------------------------------------

    slurm_log_prefix="${LOG_DIRECTORY}/${sample_prefix}.Brdu_positive.threshold_${threshold_tag}"

    # --------------------------------------------------------------------------
    # Submit the analysis job
    # --------------------------------------------------------------------------
    #
    # sbatch options:
    #
    #   --job-name
    #       Adds the sample prefix to the job name.
    #
    #   --chdir
    #       Runs the job from the workflow root.
    #
    #   --mem
    #       Explicitly requests the configured memory.
    #
    #   --output / --error
    #       Write sample- and threshold-specific SLURM logs.
    #
    #   --export
    #       Export the current environment and workflow-root variable.
    #
    # Script arguments:
    #
    #   --run-analysis
    #       Select the analysis execution path.
    #
    #   --workflow-root
    #       Preserve the resolved project root on the compute node.
    #
    #   selected_bam
    #       Input BAM filename.
    #
    #   mod_threshold
    #       Selected BrdU modification threshold.
    #

    submission_output="$(
        sbatch \
            --job-name="BrdU_mod_calls_${sample_prefix}" \
            --chdir="$WORKFLOW_ROOT" \
            --mem="$JOB_MEMORY" \
            --output="${slurm_log_prefix}.%j.log" \
            --error="${slurm_log_prefix}.%j.err" \
            --export=ALL,MOD_CALLS_BRDU_WORKFLOW_ROOT="$WORKFLOW_ROOT" \
            "$SCRIPT_PATH" \
            --run-analysis \
            --workflow-root "$WORKFLOW_ROOT" \
            "$selected_bam" \
            "$mod_threshold"
    )"

    printf '%s\n' "$submission_output"

    # Extract the final whitespace-delimited field as the submitted job ID.
    job_id="$(
        printf '%s\n' "$submission_output" |
        awk '{print $NF}'
    )"

    # --------------------------------------------------------------------------
    # Report expected output locations
    # --------------------------------------------------------------------------

    printf '\n'

    print_info "Submitted job ID: ${job_id}"

    print_info "Expected output BAM:"
    printf '  %s\n' "$output_bam"

    print_info "Expected read-ID file:"
    printf '  %s\n' "$output_read_ids"

    print_info "SLURM logs:"
    printf '  %s.%s.log\n' "$slurm_log_prefix" "$job_id"
    printf '  %s.%s.err\n' "$slurm_log_prefix" "$job_id"
}


# ==============================================================================
# FUNCTION: run_analysis
# ==============================================================================
#
# Purpose:
#   Execute BrdU-positive read extraction inside the submitted SLURM job.
#
# Arguments:
#
#   $1
#       Input BAM filename.
#
#   $2
#       BrdU modification threshold.
#
# Important behavior:
#
#   - This function refuses to run outside a SLURM job.
#   - It uses modkit extract calls rather than modkit extract full.
#   - It parses the Modkit column names instead of relying on fixed positions.
#   - It selects unique read IDs containing at least one passing b call.
#   - It extracts complete primary mapped BAM records for those read IDs.
#
run_analysis() {
    # --------------------------------------------------------------------------
    # Function arguments
    # --------------------------------------------------------------------------

    local selected_bam="$1"
    local mod_threshold="$2"

    # --------------------------------------------------------------------------
    # Analysis variables
    # --------------------------------------------------------------------------

    local input_bam
    local sample_prefix
    local threshold_tag
    local output_prefix
    local output_bam
    local output_index
    local read_ids
    local threads
    local temp_base_directory
    local temp_directory
    local temp_read_ids
    local temp_unsorted_bam
    local input_primary_reads
    local brdu_read_id_count
    local output_primary_reads
    local percent_selected
    local count_status

    # --------------------------------------------------------------------------
    # Require SLURM execution
    # --------------------------------------------------------------------------

    if [[ -z "${SLURM_JOB_ID:-}" ]]; then
        print_error "The --run-analysis mode must run through a submitted SLURM job."
        exit 1
    fi

    # --------------------------------------------------------------------------
    # Validate threshold and input BAM
    # --------------------------------------------------------------------------

    if ! validate_threshold "$mod_threshold"; then
        print_error "Invalid modification threshold: ${mod_threshold}"
        exit 1
    fi

    input_bam="${BAM_DIRECTORY}/${selected_bam}"

    if [[ ! -f "$input_bam" ]]; then
        print_error "Input BAM does not exist:"
        print_error "$input_bam"
        exit 1
    fi

    # --------------------------------------------------------------------------
    # Reconstruct output paths
    # --------------------------------------------------------------------------
    #
    # These calculations must match submit_workflow() so the reported and
    # generated output paths are identical.
    #

    sample_prefix="$(strip_bam_suffixes "$selected_bam")"

    threshold_tag="$(normalize_threshold_for_filename "$mod_threshold")"

    output_prefix="${BAM_DIRECTORY}/${sample_prefix}.Brdu_positive.threshold_${threshold_tag}"

    output_bam="${output_prefix}.bam"

    output_index="${output_bam}.bai"

    read_ids="${output_prefix}.read_ids.txt"

    # --------------------------------------------------------------------------
    # Load and validate required software
    # --------------------------------------------------------------------------

    load_required_modules

    check_required_command modkit
    check_required_command samtools
    check_required_command awk
    check_required_command sort
    check_required_command wc
    check_required_command mktemp

    # --------------------------------------------------------------------------
    # Configure threads and temporary workspace
    # --------------------------------------------------------------------------

    # Use the CPU count assigned by SLURM, with an 8-thread fallback.
    threads="${SLURM_CPUS_PER_TASK:-8}"

    # Prefer node-local SLURM temporary storage.
    temp_base_directory="${SLURM_TMPDIR:-/tmp}"

    # Create a unique temporary directory for this job.
    temp_directory="$(
        mktemp \
            -d \
            "${temp_base_directory}/mod_calls_brdu.${SLURM_JOB_ID}.XXXXXX"
    )"

    # Temporary unique read-ID file.
    temp_read_ids="${temp_directory}/Brdu_positive.read_ids.txt"

    # Temporary unsorted BAM containing selected primary records.
    temp_unsorted_bam="${temp_directory}/Brdu_positive.unsorted.bam"

    # --------------------------------------------------------------------------
    # FUNCTION: cleanup
    # --------------------------------------------------------------------------
    #
    # Remove the complete temporary workspace when the job exits.
    #

    cleanup() {
        rm -rf "$temp_directory"
    }

    trap cleanup EXIT

    # --------------------------------------------------------------------------
    # Print the analysis configuration
    # --------------------------------------------------------------------------

    printf '\n'
    printf '%s\n' '============================================================'
    printf '%s\n' 'BrdU-positive read extraction'
    printf '%s\n' '============================================================'
    printf '\n'
    printf 'Analysis date:                 %s\n' "$(date '+%Y-%m-%d %H:%M:%S')"
    printf 'SLURM job ID:                  %s\n' "$SLURM_JOB_ID"
    printf 'Workflow root:                 %s\n' "$WORKFLOW_ROOT"
    printf 'Input BAM:                     %s\n' "$input_bam"
    printf 'BrdU canonical base:           T\n'
    printf 'BrdU modification code:        b\n'
    printf 'BrdU modification threshold:   %s\n' "$mod_threshold"
    printf 'Threads:                       %s\n' "$threads"
    printf 'Output BAM:                    %s\n' "$output_bam"
    printf 'Read-ID file:                  %s\n' "$read_ids"
    printf '\n'
    printf 'Requested memory:              %s\n' "$JOB_MEMORY"
    printf 'Modkit version:                %s\n' "$(modkit --version 2>&1 | head -n 1)"
    printf 'Samtools version:              %s\n' "$(samtools --version | head -n 1)"
    printf '\n'

    # --------------------------------------------------------------------------
    # Recheck output existence on the compute node
    # --------------------------------------------------------------------------

    for output_file in "$output_bam" "$output_index" "$read_ids"; do
        if [[ -e "$output_file" ]]; then
            print_error "Output already exists: ${output_file}"
            exit 1
        fi
    done

    # --------------------------------------------------------------------------
    # Create the input BAM index when missing
    # --------------------------------------------------------------------------
    #
    # Both common index naming patterns are accepted:
    #
    #   input.bam.bai
    #   input.bai
    #

    if [[ ! -f "${input_bam}.bai" && ! -f "${input_bam%.bam}.bai" ]]; then
        print_info "Input BAM index was not found."
        print_info "Creating input BAM index."

        samtools index \
            -@ "$threads" \
            "$input_bam"
    fi

    # --------------------------------------------------------------------------
    # Count mapped primary reads in the input BAM
    # --------------------------------------------------------------------------

    print_info "Counting primary mapped reads in the input BAM."

    input_primary_reads="$(
        samtools view \
            -@ "$threads" \
            -c \
            -F 2308 \
            "$input_bam"
    )"

    print_info "Input primary mapped reads: ${input_primary_reads}"

    # --------------------------------------------------------------------------
    # Find unique reads containing passing BrdU calls
    # --------------------------------------------------------------------------
    #
    # Modkit command:
    #
    #   modkit extract calls INPUT -
    #
    # writes call-level information to standard output.
    #
    # Options:
    #
    #   --mapped-only
    #       Include calls from mapped reads only.
    #
    #   --threads
    #       Use the SLURM CPU allocation.
    #
    #   --mod-threshold "b:<threshold>"
    #       Apply the selected threshold specifically to modification code b.
    #
    # AWK processing:
    #
    #   1. Ignore comment lines beginning with #.
    #   2. Read the first non-comment line as the header.
    #   3. Locate read_id, call_code, and fail by column name.
    #   4. Reject output when required columns are missing.
    #   5. Print read IDs where call_code is b and fail is false.
    #
    # sort --unique:
    #
    #   A read may have many passing BrdU calls. Sorting and deduplicating
    #   ensures each selected read ID appears exactly once.
    #

    print_info "Finding reads with at least one passing BrdU call using modkit extract calls."

    modkit extract calls \
        "$input_bam" \
        - \
        --mapped-only \
        --threads "$threads" \
        --mod-threshold "b:${mod_threshold}" \
    | awk -F'\t' '
        BEGIN {
            read_id_column = 0
            call_code_column = 0
            fail_column = 0
            header_found = 0
        }

        # Ignore Modkit metadata/comment lines.
        /^#/ {
            next
        }

        # Parse the first non-comment line as the table header.
        header_found == 0 {
            for (column = 1; column <= NF; column++) {
                header_value = $column
                gsub(/\r/, "", header_value)

                if (header_value == "read_id") {
                    read_id_column = column
                } else if (header_value == "call_code") {
                    call_code_column = column
                } else if (header_value == "fail") {
                    fail_column = column
                }
            }

            # Stop if the expected schema is unavailable.
            if (
                read_id_column == 0 ||
                call_code_column == 0 ||
                fail_column == 0
            ) {
                print \
                    "[ERROR] Could not find read_id, call_code, and fail columns." \
                    > "/dev/stderr"

                exit 2
            }

            header_found = 1
            next
        }

        # Process call-level records.
        {
            read_id_value = $read_id_column
            call_code_value = $call_code_column
            fail_value = $fail_column

            # Remove possible carriage-return characters.
            gsub(/\r/, "", read_id_value)
            gsub(/\r/, "", call_code_value)
            gsub(/\r/, "", fail_value)

            # Retain passing BrdU calls only.
            if (
                call_code_value == "b" &&
                tolower(fail_value) == "false"
            ) {
                print read_id_value
            }
        }

        END {
            # Report completely missing or empty Modkit tabular output.
            if (header_found == 0) {
                print \
                    "[ERROR] No modkit extract calls header was found." \
                    > "/dev/stderr"

                exit 3
            }
        }
    ' \
    | LC_ALL=C sort \
        --unique \
        --temporary-directory="$temp_directory" \
    > "$temp_read_ids"

    # --------------------------------------------------------------------------
    # Validate and preserve the selected read-ID list
    # --------------------------------------------------------------------------

    if [[ ! -s "$temp_read_ids" ]]; then
        print_error "No reads had a passing BrdU call at threshold ${mod_threshold}."
        exit 1
    fi

    brdu_read_id_count="$(wc -l < "$temp_read_ids")"

    print_info "Unique reads with passing BrdU calls: ${brdu_read_id_count}"

    # Move the completed read-ID list from temporary storage to its final path.
    mv "$temp_read_ids" "$read_ids"

    # --------------------------------------------------------------------------
    # Extract complete primary mapped records for selected reads
    # --------------------------------------------------------------------------
    #
    # samtools view options:
    #
    #   -F 2308
    #       Exclude unmapped, secondary, and supplementary alignments.
    #
    #   -N
    #       Retain alignments whose read names appear in the read-ID file.
    #
    #   -b
    #       Write BAM output.
    #
    # The output contains complete BAM alignment records, not individual Modkit
    # call rows.
    #

    print_info "Extracting complete primary mapped BrdU-positive reads."

    samtools view \
        -@ "$threads" \
        -F 2308 \
        -N "$read_ids" \
        -b \
        -o "$temp_unsorted_bam" \
        "$input_bam"

    # --------------------------------------------------------------------------
    # Coordinate-sort the selected BAM
    # --------------------------------------------------------------------------

    print_info "Coordinate-sorting the output BAM."

    samtools sort \
        -@ "$threads" \
        -o "$output_bam" \
        "$temp_unsorted_bam"

    # --------------------------------------------------------------------------
    # Index the final BAM
    # --------------------------------------------------------------------------

    print_info "Indexing the output BAM."

    samtools index \
        -@ "$threads" \
        "$output_bam"

    # --------------------------------------------------------------------------
    # Validate the final BAM
    # --------------------------------------------------------------------------

    print_info "Validating the output BAM."

    samtools quickcheck -v "$output_bam"

    # --------------------------------------------------------------------------
    # Count primary mapped reads in the output BAM
    # --------------------------------------------------------------------------

    output_primary_reads="$(
        samtools view \
            -@ "$threads" \
            -c \
            -F 2308 \
            "$output_bam"
    )"

    # --------------------------------------------------------------------------
    # Calculate the percentage of input reads selected
    # --------------------------------------------------------------------------

    percent_selected="$(
        awk \
            -v selected="$output_primary_reads" \
            -v total="$input_primary_reads" \
            'BEGIN {
                if (total > 0) {
                    printf "%.6f", (selected / total) * 100
                } else {
                    printf "0.000000"
                }
            }'
    )"

    # --------------------------------------------------------------------------
    # Compare selected read IDs with output primary records
    # --------------------------------------------------------------------------
    #
    # These counts should normally match because one primary mapped alignment is
    # expected for each unique selected read ID.
    #
    # A mismatch may indicate duplicate primary records for one or more read
    # names in the input BAM.
    #

    if [[ "$output_primary_reads" -eq "$brdu_read_id_count" ]]; then
        count_status="PASS"
    else
        count_status="WARNING"
    fi

    # --------------------------------------------------------------------------
    # Print the final analysis summary
    # --------------------------------------------------------------------------

    printf '\n'
    printf '%s\n' '============================================================'
    printf '%s\n' 'BrdU-positive read extraction complete'
    printf '%s\n' '============================================================'
    printf '\n'
    printf 'Input primary mapped reads:    %s\n' "$input_primary_reads"
    printf 'BrdU-positive read IDs:        %s\n' "$brdu_read_id_count"
    printf 'Output primary mapped reads:   %s\n' "$output_primary_reads"
    printf 'Percent of reads selected:     %s%%\n' "$percent_selected"
    printf 'Read-count agreement:          %s\n' "$count_status"
    printf '\n'
    printf 'Output BAM:                    %s\n' "$output_bam"
    printf 'Output BAM index:              %s\n' "$output_index"
    printf 'Selected read IDs:             %s\n' "$read_ids"
    printf '\n'

    if [[ "$count_status" == "PASS" ]]; then
        print_info "Output BAM count matches the selected read-ID count."
    else
        print_warning "Output BAM count does not match the selected read-ID count."
        print_warning "Check whether the input BAM contains duplicate primary records."
    fi

    printf '\n'
    printf '%s\n' '============================================================'
}


# ==============================================================================
# MAIN SCRIPT ROUTING
# ==============================================================================
#
# -h or --help
#   Print usage information.
#
# --run-analysis
#   Execute the submitted SLURM analysis stage.
#
# Any other invocation
#   Enter interactive submission mode.
#

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

if [[ "${1:-}" == "--run-analysis" ]]; then
    # Remove the execution-mode argument.
    shift

    # Skip the optional --workflow-root argument because the root was already
    # resolved near the beginning of the script.
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --workflow-root)
                shift 2
                ;;
            *)
                break
                ;;
        esac
    done

    # The analysis requires exactly:
    #
    #   BAM filename
    #   Modification threshold
    #
    if (( $# != 2 )); then
        print_error "Expected: --run-analysis [--workflow-root PATH] <bam_filename> <mod_threshold>"
        exit 1
    fi

    run_analysis "$1" "$2"
else
    submit_workflow "$@"
fi