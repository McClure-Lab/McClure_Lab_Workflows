#!/bin/bash
#SBATCH --job-name=dnascent_workflow
#SBATCH --partition=gpu
#SBATCH --account=gpu_rbi
#SBATCH --gpus=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=36G
#SBATCH --time=12:00:00

# Stop immediately if a command fails.
# Treat unset variables as errors.
# Make a pipeline fail when any command within the pipeline fails.
set -euo pipefail

# Search the original command-line arguments for an optional
# --workflow-root argument before constructing the workflow paths.
WORKFLOW_ROOT_ARG=""
for ((arg_i = 1; arg_i <= $#; arg_i++)); do
    if [[ "${!arg_i}" == "--workflow-root" ]]; then
        next_arg_i=$((arg_i + 1))
        WORKFLOW_ROOT_ARG="${!next_arg_i:-}"
        break
    fi
done

# Determine the root directory of the workflow.
#
# Priority:
# 1. A path supplied through --workflow-root.
# 2. The DNASCENT_WORKFLOW_ROOT environment variable.
# 3. Two directories above the location of this script.
if [[ -n "$WORKFLOW_ROOT_ARG" ]]; then
    WORKFLOW_ROOT="$(cd "$WORKFLOW_ROOT_ARG" && pwd)"
elif [[ -n "${DNASCENT_WORKFLOW_ROOT:-}" ]]; then
    WORKFLOW_ROOT="$DNASCENT_WORKFLOW_ROOT"
else
    WORKFLOW_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
fi

# Define the directories and files used by the workflow.
SRC_DIR="$WORKFLOW_ROOT/src/nanopore_sequence_workflow"
SCRIPT_PATH="$SRC_DIR/dnascent_workflow.sh"
LOG_DIR="$WORKFLOW_ROOT/logs/dnascent_workflow"
BAM_DIR="$WORKFLOW_ROOT/data/bam"
FASTQ_DIR="$WORKFLOW_ROOT/data/fastq"
POD5_DIR="$WORKFLOW_ROOT/data/pod5"

# Use the DNASCENT_IMAGE environment variable when it is set.
# Otherwise, use the cluster's default DNAscent Singularity image.
DEFAULT_DNASCENT_IMAGE="${DNASCENT_IMAGE:-/cluster/singularity_images/DNAscent.sif}"

# Print instructions explaining how to run the workflow and where
# the workflow expects the input files to be located.
usage() {
    echo "Usage: bash src/nanopore_sequence_workflow/dnascent_workflow.sh [BAM] [BAM_INDEX] [REFERENCE_FASTA] [DNASCENT_INDEX] [POD5]"
    echo
    echo "BAM and BAM index may be filenames in data/bam or explicit paths."
    echo "Reference FASTA may be a filename in data/fastq or an explicit path."
    echo "DNAscent index and POD5 may be filenames in data/pod5 or explicit paths."
    echo
    echo "If arguments are omitted, the workflow prompts for them before submitting a SLURM job."
}

# Convert an existing file path into an absolute path.
#
# The directory and filename are handled separately so the resulting
# path points to the exact existing file.
absolute_existing_file() {
    local path="$1"
    local dir
    local file

    dir="$(cd "$(dirname "$path")" && pwd)"
    file="$(basename "$path")"
    printf '%s/%s\n' "$dir" "$file"
}

# Resolve an input filename or path.
#
# The function checks for the file in this order:
# 1. The input exactly as entered.
# 2. The expected default directory.
# 3. Anywhere under WORKFLOW_ROOT/data when fallback searching is enabled.
resolve_file_in_dir() {
    local input="$1"
    local default_dir="$2"
    local fallback_data_search="${3:-no}"
    local path

    # Check whether the input is already a valid file path.
    if [[ -f "$input" ]]; then
        absolute_existing_file "$input"
        return 0
    fi

    # Check whether the input is a filename inside the expected directory.
    path="$default_dir/$input"
    if [[ -f "$path" ]]; then
        absolute_existing_file "$path"
        return 0
    fi

    # Optionally search all directories under WORKFLOW_ROOT/data.
    # This is enabled for reference FASTA files.
    if [[ "$fallback_data_search" == "yes" ]]; then
        path="$(find "$WORKFLOW_ROOT/data" -type f -name "$input" -print -quit)"
        if [[ -n "$path" && -f "$path" ]]; then
            absolute_existing_file "$path"
            return 0
        fi
    fi

    # Return a failure status when the file cannot be found.
    return 1
}

# Resolve a BAM file from either an explicit path or data/bam.
resolve_bam_file() {
    resolve_file_in_dir "$1" "$BAM_DIR" "no"
}

# Resolve a BAM index from either an explicit path or data/bam.
resolve_bam_index_file() {
    local index_input="$1"
    resolve_file_in_dir "$index_input" "$BAM_DIR" "no"
}

# Resolve a reference FASTA.
#
# Reference files are first checked directly and in data/fastq.
# If they are not found there, the complete data directory is searched.
resolve_reference_file() {
    resolve_file_in_dir "$1" "$FASTQ_DIR" "yes"
}

# Resolve a POD5 file from either an explicit path or data/pod5.
resolve_pod5_file() {
    resolve_file_in_dir "$1" "$POD5_DIR" "no"
}

# Resolve a DNAscent index from either an explicit path or data/pod5.
resolve_dnascent_index_file() {
    resolve_file_in_dir "$1" "$POD5_DIR" "no"
}

# Create a sample prefix from the BAM filename.
#
# For example:
# sample.sorted.indexed.bam becomes sample
derive_bam_prefix() {
    local bam_basename="$1"
    local prefix

    prefix="${bam_basename%.bam}"
    prefix="${prefix%.sorted.indexed}"
    printf '%s\n' "$prefix"
}

# Prompt the user for a value.
#
# Pressing Enter without typing a value returns the supplied default.
prompt_with_default() {
    local prompt="$1"
    local default="$2"
    local value

    read -r -p "$prompt [$default]: " value
    printf '%s\n' "${value:-$default}"
}

# Confirm that a value is a positive integer greater than zero.
validate_positive_int() {
    local name="$1"
    local value="$2"

    if ! [[ "$value" =~ ^[1-9][0-9]*$ ]]; then
        echo "[ERROR] $name must be a positive integer. Received: $value"
        exit 1
    fi
}

# Confirm that a value is a non-negative integer.
# Zero is allowed by this validation function.
validate_nonnegative_int() {
    local name="$1"
    local value="$2"

    if ! [[ "$value" =~ ^[0-9]+$ ]]; then
        echo "[ERROR] $name must be a non-negative integer. Received: $value"
        exit 1
    fi
}

# Confirm that a required directory exists.
#
# The compute-node job intentionally does not create these directories.
# This prevents output from being written to an unexpected location when
# the project directory is not mounted correctly.
require_existing_dir() {
    local name="$1"
    local path="$2"

    if [[ ! -d "$path" ]]; then
        echo "[ERROR] $name does not exist: $path"
        echo "[ERROR] The SLURM job will not create this directory on the compute node."
        echo "[ERROR] Confirm the project path is mounted on the compute node and the workflow was submitted from the repository."
        exit 1
    fi
}

# Load the programs needed by the DNAscent workflow when the cluster's
# environment-module system is available.
load_dnascent_modules() {
    if command -v module >/dev/null 2>&1; then
        echo "[INFO] Loading DNAscent workflow modules: DNAscent, samtools, singularity"

        # Loading DNAscent is allowed to fail because DNAscent is ultimately
        # run from the supplied container image.
        module load DNAscent || true

        # Samtools is required outside the container for BAM validation.
        module load samtools

        # Try Singularity first and Apptainer second.
        # A runtime already available in PATH may also be used.
        module load singularity || module load apptainer || true
    else
        echo "[WARN] Environment modules are not available in this shell."
        echo "[WARN] Continuing with DNAscent container runtime and samtools from PATH."
    fi
}

# Determine which supported container runtime is available.
#
# Singularity is preferred when both programs are installed.
container_runtime() {
    if command -v singularity >/dev/null 2>&1; then
        printf '%s\n' "singularity"
    elif command -v apptainer >/dev/null 2>&1; then
        printf '%s\n' "apptainer"
    else
        echo "[ERROR] Neither singularity nor apptainer was found in PATH." >&2
        exit 1
    fi
}

# Run the compute-node portion of the workflow.
#
# This function is called after the script submits itself to SLURM
# using the --run-job argument.
run_job() {
    # Initialize all job parameters.
    local bam_file=""
    local bam_index=""
    local reference=""
    local dnascent_index=""
    local pod5_file=""
    local dnascent_image=""
    local threads=""
    local gpu=""
    local min_qscore=""
    local min_read_length=""
    local output_dir=""

    # Remove the initial --run-job argument.
    shift

    # Parse the named arguments passed by submit_workflow.
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --workflow-root) shift 2 ;;
            --bam) bam_file="$2"; shift 2 ;;
            --bam-index) bam_index="$2"; shift 2 ;;
            --reference) reference="$2"; shift 2 ;;
            --dnascent-index) dnascent_index="$2"; shift 2 ;;
            --pod5) pod5_file="$2"; shift 2 ;;
            --dnascent-image) dnascent_image="$2"; shift 2 ;;
            --threads) threads="$2"; shift 2 ;;
            --gpu) gpu="$2"; shift 2 ;;
            --min-qscore) min_qscore="$2"; shift 2 ;;
            --min-read-length) min_read_length="$2"; shift 2 ;;
            --output-dir) output_dir="$2"; shift 2 ;;
            *)
                echo "[ERROR] Unknown job argument: $1"
                exit 1
                ;;
        esac
    done

    # Confirm that the workflow, log, and output directories are visible
    # from the SLURM compute node.
    require_existing_dir "Workflow root" "$WORKFLOW_ROOT"
    require_existing_dir "Log directory" "$LOG_DIR"
    require_existing_dir "Output directory" "$output_dir"

    # Verify that every required input file exists before running DNAscent.
    for required_file in \
        "$bam_file" \
        "$bam_index" \
        "$reference" \
        "$dnascent_index" \
        "$pod5_file" \
        "$dnascent_image"
    do
        if [[ ! -f "$required_file" ]]; then
            echo "[ERROR] Required file not found: $required_file" >&2
            exit 1
        fi
    done

    # Prepare the required software environment.
    load_dnascent_modules

    local runtime
    local bam_basename
    local bam_prefix
    local output_bam
    local detect_log

    # Select Singularity or Apptainer.
    runtime="$(container_runtime)"

    # Derive the sample prefix from the input BAM filename.
    bam_basename="$(basename "$bam_file")"
    bam_prefix="$(derive_bam_prefix "$bam_basename")"

    # Construct the output BAM and DNAscent log filenames.
    output_bam="$output_dir/${bam_prefix}.sorted.indexed.BrdU.detect.bam"
    detect_log="$LOG_DIR/${bam_prefix}.${SLURM_JOB_ID:-manual}.dnascent.log"

    # Write a complete summary of the run configuration to the DNAscent log.
    {
        echo "========================================="
        echo "  DNASCENT BRDU DETECTION"
        echo "========================================="
        echo "Started:          $(date)"
        echo "Workflow root:    $WORKFLOW_ROOT"
        echo "BAM:              $bam_file"
        echo "BAM index:        $bam_index"
        echo "Reference:        $reference"
        echo "DNAscent index:   $dnascent_index"
        echo "POD5:             $pod5_file"
        echo "Container image:  $dnascent_image"
        echo "Output BAM:       $output_bam"
        echo "Threads:          $threads"
        echo "GPU:              $gpu"
        echo "Min q-score:      $min_qscore"
        echo "Min read length:  $min_read_length"
        echo "SLURM job ID:     ${SLURM_JOB_ID:-not_set}"
        echo "CUDA devices:     ${CUDA_VISIBLE_DEVICES:-not_set}"
        echo "========================================="
    } > "$detect_log"

    # Check the structural integrity of the input BAM before running DNAscent.
    echo "[INFO] Checking input BAM." | tee -a "$detect_log"
    samtools quickcheck -v "$bam_file" >> "$detect_log" 2>&1

    # Record information about the assigned GPU when nvidia-smi is available.
    # Failure to run nvidia-smi does not stop the workflow.
    if command -v nvidia-smi >/dev/null 2>&1; then
        {
            echo ""
            echo "========================================="
            echo "  GPU STATUS"
            echo "========================================="
            nvidia-smi
        } >> "$detect_log" 2>&1 || true
    fi

    # Store the complete Singularity/Apptainer and DNAscent command in an array.
    #
    # run:
    #   Starts the container's default run action.
    #
    # --nv:
    #   Makes NVIDIA GPU libraries and devices available inside the container.
    #
    # --bind:
    #   Makes the workflow directory available at the same path inside the container.
    #
    # detect:
    #   Runs DNAscent's BrdU detection command.
    #
    # -t:
    #   Number of CPU threads.
    #
    # --GPU:
    #   GPU device DNAscent should use.
    #
    # -q:
    #   Minimum quality threshold.
    #
    # -l:
    #   Minimum read length.
    #
    # -b:
    #   Input aligned BAM.
    #
    # -r:
    #   Reference FASTA.
    #
    # -i:
    #   DNAscent index connecting read identifiers to POD5 signal.
    #
    # -o:
    #   Output BAM containing DNAscent BrdU calls.
    local dnascent_args=(
        run
        --nv
        --bind "$WORKFLOW_ROOT:$WORKFLOW_ROOT"
        "$dnascent_image"
        detect
        -t "$threads"
        --GPU "$gpu"
        -q "$min_qscore"
        -l "$min_read_length"
        -b "$bam_file"
        -r "$reference"
        -i "$dnascent_index"
        -o "$output_bam"
    )

    # Add a shell-escaped representation of the command to the log.
    # This makes it easier to reproduce or troubleshoot the exact command.
    {
        echo ""
        echo "[INFO] Running DNAscent detect."
        printf '[INFO] Command: %s' "$runtime"
        printf ' %q' "${dnascent_args[@]}"
        echo
    } >> "$detect_log"

    # Run DNAscent inside the selected GPU-enabled container.
    # Standard output and errors are both appended to the DNAscent log.
    "$runtime" "${dnascent_args[@]}" >> "$detect_log" 2>&1

    # Confirm that DNAscent produced a nonempty output BAM.
    if [[ -s "$output_bam" ]]; then
        {
            echo ""
            echo "[INFO] DNAscent detect completed: $(date)"
            ls -lh "$output_bam"
        } | tee -a "$detect_log"
    else
        echo "[ERROR] DNAscent did not produce a nonempty output BAM: $output_bam" | tee -a "$detect_log" >&2
        exit 1
    fi
}

# Run the interactive login-node portion of the workflow.
#
# This function:
# 1. Collects input files and parameters.
# 2. Resolves all input paths.
# 3. Validates the supplied values.
# 4. Submits the actual DNAscent job to SLURM.
submit_workflow() {
    # The five primary inputs may be supplied as positional arguments.
    # Missing values are collected interactively later in the function.
    local bam_input="${1:-}"
    local bam_index_input="${2:-}"
    local reference_input="${3:-}"
    local dnascent_index_input="${4:-}"
    local pod5_input="${5:-}"

    # Variables that will store resolved absolute file paths.
    local bam_file=""
    local bam_index=""
    local reference=""
    local dnascent_index=""
    local pod5_file=""

    # Additional workflow and SLURM parameters.
    local bam_basename
    local bam_prefix
    local default_bam_index
    local dnascent_image
    local threads
    local gpu
    local min_qscore
    local min_read_length
    local slurm_mem
    local slurm_time
    local output_dir
    local job_id
    local log_job_id

    # Create the workflow's standard directories on the submission node.
    mkdir -p "$LOG_DIR" "$BAM_DIR" "$FASTQ_DIR" "$POD5_DIR"

    # When the BAM was not provided as an argument, display available
    # BAM files and prompt the user to select one.
    if [[ -z "$bam_input" ]]; then
        echo "[INFO] Available BAM files in $BAM_DIR:"
        find "$BAM_DIR" -maxdepth 1 -type f -name "*.bam" -printf "  %f\n" | sort
        echo
        read -r -p "Enter BAM filename from data/bam or an explicit path: " bam_input
    fi

    # Stop if the user did not provide a BAM.
    if [[ -z "$bam_input" ]]; then
        echo "[ERROR] BAM file is required."
        exit 1
    fi

    # Resolve the BAM to an absolute existing path.
    if ! bam_file="$(resolve_bam_file "$bam_input")"; then
        echo "[ERROR] BAM file not found: $bam_input"
        echo "[ERROR] Expected a valid path or a filename in $BAM_DIR"
        exit 1
    fi

    # Use BAM.bai as the default index path.
    #
    # For example:
    # sample.bam becomes sample.bam.bai
    default_bam_index="${bam_file}.bai"

    # Prompt for the BAM index when one was not supplied as an argument.
    if [[ -z "$bam_index_input" ]]; then
        bam_index_input="$(prompt_with_default "Enter BAM index filename from data/bam or an explicit path" "$default_bam_index")"
    fi

    # Display available FASTA files and prompt for the reference when
    # one was not supplied as an argument.
    if [[ -z "$reference_input" ]]; then
        echo "[INFO] Available reference FASTA files in $FASTQ_DIR:"
        find "$FASTQ_DIR" -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" \) -printf "  %f\n" | sort
        echo
        read -r -p "Enter reference FASTA filename from data/fastq, data/, or an explicit path: " reference_input
    fi

    # Display files in data/pod5 and prompt for the DNAscent index when
    # one was not supplied as an argument.
    if [[ -z "$dnascent_index_input" ]]; then
        echo "[INFO] Available DNAscent index files in $POD5_DIR:"
        find "$POD5_DIR" -maxdepth 1 -type f -printf "  %f\n" | sort
        echo
        read -r -p "Enter DNAscent index filename from data/pod5 or an explicit path: " dnascent_index_input
    fi

    # Display available POD5 files and prompt for the raw signal file
    # when one was not supplied as an argument.
    if [[ -z "$pod5_input" ]]; then
        echo "[INFO] Available POD5 files in $POD5_DIR:"
        find "$POD5_DIR" -maxdepth 1 -type f -name "*.pod5" -printf "  %f\n" | sort
        echo
        read -r -p "Enter POD5 filename from data/pod5 or an explicit path: " pod5_input
    fi

    # Confirm that all required input values were provided.
    if [[ -z "$bam_index_input" || -z "$reference_input" || -z "$dnascent_index_input" || -z "$pod5_input" ]]; then
        echo "[ERROR] BAM index, reference FASTA, DNAscent index, and POD5 are required."
        exit 1
    fi

    # Resolve the BAM index to an absolute existing path.
    if ! bam_index="$(resolve_bam_index_file "$bam_index_input")"; then
        echo "[ERROR] BAM index not found: $bam_index_input"
        echo "[ERROR] Expected a valid path or a filename in $BAM_DIR"
        exit 1
    fi

    # Resolve the reference FASTA to an absolute existing path.
    if ! reference="$(resolve_reference_file "$reference_input")"; then
        echo "[ERROR] Reference FASTA not found: $reference_input"
        echo "[ERROR] Expected a valid path, a filename in $FASTQ_DIR, or an exact filename under $WORKFLOW_ROOT/data"
        exit 1
    fi

    # Resolve the DNAscent index to an absolute existing path.
    if ! dnascent_index="$(resolve_dnascent_index_file "$dnascent_index_input")"; then
        echo "[ERROR] DNAscent index not found: $dnascent_index_input"
        echo "[ERROR] Expected a valid path or a filename in $POD5_DIR"
        exit 1
    fi

    # Resolve the POD5 file to an absolute existing path.
    if ! pod5_file="$(resolve_pod5_file "$pod5_input")"; then
        echo "[ERROR] POD5 file not found: $pod5_input"
        echo "[ERROR] Expected a valid path or a filename in $POD5_DIR"
        exit 1
    fi

    # Derive the output prefix from the BAM name.
    bam_basename="$(basename "$bam_file")"
    bam_prefix="$(derive_bam_prefix "$bam_basename")"

    # Write the DNAscent output BAM to data/bam.
    output_dir="$BAM_DIR"

    echo
    echo "DNAscent detect parameters"

    # Collect the DNAscent container and runtime parameters.
    # The values in brackets are used when the user presses Enter.
    dnascent_image="$(prompt_with_default "Enter DNAscent container image path" "$DEFAULT_DNASCENT_IMAGE")"
    threads="$(prompt_with_default "Enter DNAscent threads / SLURM cpus-per-task" "12")"
    gpu="$(prompt_with_default "Enter DNAscent --GPU value" "0")"
    min_qscore="$(prompt_with_default "Enter DNAscent -q minimum alignment/read quality" "20")"
    min_read_length="$(prompt_with_default "Enter DNAscent -l minimum read length" "1000")"
    slurm_mem="$(prompt_with_default "Enter SLURM memory" "36G")"
    slurm_time="$(prompt_with_default "Enter SLURM time" "12:00:00")"

    # Validate the parameters that must be numeric.
    validate_positive_int "--threads" "$threads"
    validate_nonnegative_int "--GPU" "$gpu"
    validate_nonnegative_int "-q" "$min_qscore"
    validate_positive_int "-l" "$min_read_length"

    # Confirm that the selected DNAscent container exists before
    # submitting the job.
    if [[ ! -f "$dnascent_image" ]]; then
        echo "[ERROR] DNAscent container image not found: $dnascent_image"
        exit 1
    fi

    # Display all resolved inputs and the expected output path.
    echo
    echo "[INFO] Workflow root: $WORKFLOW_ROOT"
    echo "[INFO] BAM: $bam_file"
    echo "[INFO] BAM index: $bam_index"
    echo "[INFO] Reference: $reference"
    echo "[INFO] DNAscent index: $dnascent_index"
    echo "[INFO] POD5: $pod5_file"
    echo "[INFO] Expected output BAM: $output_dir/${bam_prefix}.sorted.indexed.BrdU.detect.bam"
    echo
    echo "Submitting DNAscent workflow to SLURM..."

    # Submit this same script to SLURM in --run-job mode.
    #
    # --parsable:
    #   Returns a machine-readable job ID.
    #
    # --job-name:
    #   Creates a job name based on the BAM prefix.
    #
    # --chdir:
    #   Runs the job from the workflow root.
    #
    # --output and --error:
    #   Store SLURM output and errors in the DNAscent log directory.
    #
    # --export:
    #   Preserves the environment and explicitly passes the workflow root.
    #
    # --cpus-per-task, --mem, and --time:
    #   Override the matching SBATCH header values with the user's selections.
    #
    # Everything after SCRIPT_PATH is passed to this script as command-line
    # arguments when the compute-node job begins.
    job_id="$(sbatch --parsable \
        --job-name="dnascent_${bam_prefix}" \
        --chdir="$WORKFLOW_ROOT" \
        --output="$LOG_DIR/${bam_prefix}.%j.slurm.log" \
        --error="$LOG_DIR/${bam_prefix}.%j.slurm.err" \
        --export=ALL,DNASCENT_WORKFLOW_ROOT="$WORKFLOW_ROOT" \
        --cpus-per-task="$threads" \
        --mem="$slurm_mem" \
        --time="$slurm_time" \
        "$SCRIPT_PATH" \
        --run-job \
        --workflow-root "$WORKFLOW_ROOT" \
        --bam "$bam_file" \
        --bam-index "$bam_index" \
        --reference "$reference" \
        --dnascent-index "$dnascent_index" \
        --pod5 "$pod5_file" \
        --dnascent-image "$dnascent_image" \
        --threads "$threads" \
        --gpu "$gpu" \
        --min-qscore "$min_qscore" \
        --min-read-length "$min_read_length" \
        --output-dir "$output_dir")"

    # A parsable SLURM job ID may contain a cluster suffix separated
    # from the numeric ID by a semicolon. Keep only the local job ID
    # when constructing the expected log filenames.
    log_job_id="${job_id%%;*}"

    # Display the submitted job ID and all expected output/log locations.
    echo "[INFO] Submitted SLURM job: $job_id"
    echo "[INFO] Expected BAM:  $output_dir/${bam_prefix}.sorted.indexed.BrdU.detect.bam"
    echo "[INFO] SLURM log:     $LOG_DIR/${bam_prefix}.${log_job_id}.slurm.log"
    echo "[INFO] SLURM err:     $LOG_DIR/${bam_prefix}.${log_job_id}.slurm.err"
    echo "[INFO] DNAscent log:  $LOG_DIR/${bam_prefix}.${log_job_id}.dnascent.log"
}

# Main script entry point.
#
# -h or --help:
#   Print the usage instructions.
#
# --run-job:
#   Run the compute-node portion after SLURM starts the job.
#
# Any other invocation:
#   Collect inputs and submit the workflow to SLURM.
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

if [[ "${1:-}" == "--run-job" ]]; then
    run_job "$@"
else
    submit_workflow "$@"
fi