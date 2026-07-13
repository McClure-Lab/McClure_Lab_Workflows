#!/bin/bash
#SBATCH --job-name=nanopore_workflow
#SBATCH --partition=gpu
#SBATCH --account=gpu_rbi
#SBATCH --gpus=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00

# Stop the script immediately if:
# - a command exits with an error,
# - an unset variable is referenced,
# - any command within a pipeline fails.
set -euo pipefail

# Search the command-line arguments for an optional --workflow-root value
# before the workflow constructs its standard directory paths.
WORKFLOW_ROOT_ARG=""
for ((arg_i = 1; arg_i <= $#; arg_i++)); do
    if [[ "${!arg_i}" == "--workflow-root" ]]; then
        next_arg_i=$((arg_i + 1))
        WORKFLOW_ROOT_ARG="${!next_arg_i:-}"
        break
    fi
done

# Determine the root directory of the nanopore workflow.
#
# Priority:
# 1. A path supplied with --workflow-root.
# 2. The NANOPORE_WORKFLOW_ROOT environment variable.
# 3. Two directories above this script's location.
if [[ -n "$WORKFLOW_ROOT_ARG" ]]; then
    WORKFLOW_ROOT="$(cd "$WORKFLOW_ROOT_ARG" && pwd)"
elif [[ -n "${NANOPORE_WORKFLOW_ROOT:-}" ]]; then
    WORKFLOW_ROOT="$NANOPORE_WORKFLOW_ROOT"
else
    WORKFLOW_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
fi

# Define the workflow's script, log, input, and output directories.
SRC_DIR="$WORKFLOW_ROOT/src/nanopore_sequence_workflow"
SCRIPT_PATH="$SRC_DIR/nanopore_workflow.sh"
LOG_DIR="$WORKFLOW_ROOT/logs/nanopore_sequence_workflow"
POD5_DIR="$WORKFLOW_ROOT/data/pod5"
DEFAULT_FASTQ_DIR="$WORKFLOW_ROOT/data/fastq"
DEFAULT_BAM_DIR="$WORKFLOW_ROOT/data/bam"

# Display the expected positional arguments and default workflow behavior.
usage() {
    echo "Usage: bash src/nanopore_sequence_workflow/nanopore_workflow.sh [POD5] [REFERENCE_FASTA] [FASTQ_OUTPUT_DIR] [BAM_OUTPUT_DIR]"
    echo
    echo "POD5 may be a filename in data/pod5 or an explicit path."
    echo "FASTQ and BAM outputs default to data/fastq and data/bam."
    echo
    echo "If arguments are omitted, the workflow prompts for them before submitting a SLURM job."
}

# Convert an output-directory value into a complete path.
#
# Behavior:
# - Absolute paths are returned unchanged.
# - Paths beginning with data/ are placed under WORKFLOW_ROOT.
# - Other relative values are placed under WORKFLOW_ROOT/data.
resolve_data_output_dir() {
    local path="$1"

    if [[ "$path" = /* ]]; then
        printf '%s\n' "$path"
    elif [[ "$path" == data/* ]]; then
        printf '%s/%s\n' "$WORKFLOW_ROOT" "$path"
    else
        printf '%s/data/%s\n' "$WORKFLOW_ROOT" "$path"
    fi
}

# Convert an existing file path into an absolute path.
absolute_existing_file() {
    local path="$1"
    local dir
    local file

    dir="$(cd "$(dirname "$path")" && pwd)"
    file="$(basename "$path")"
    printf '%s/%s\n' "$dir" "$file"
}

# Resolve the POD5 input.
#
# The function first checks the input exactly as entered.
# If it is not found, it checks data/pod5.
resolve_pod5_file() {
    local pod5_input="$1"
    local pod5_path

    if [[ -f "$pod5_input" ]]; then
        absolute_existing_file "$pod5_input"
        return 0
    fi

    pod5_path="$POD5_DIR/$pod5_input"
    if [[ -f "$pod5_path" ]]; then
        absolute_existing_file "$pod5_path"
        return 0
    fi

    return 1
}

# Resolve the reference FASTA path.
#
# The function checks:
# 1. The exact input path.
# 2. data/fastq.
# 3. The top level of data/.
# 4. Any matching filename below data/.
resolve_reference_file() {
    local reference_input="$1"
    local reference_path

    if [[ -f "$reference_input" ]]; then
        absolute_existing_file "$reference_input"
        return 0
    fi

    reference_path="$DEFAULT_FASTQ_DIR/$reference_input"
    if [[ -f "$reference_path" ]]; then
        absolute_existing_file "$reference_path"
        return 0
    fi

    reference_path="$WORKFLOW_ROOT/data/$reference_input"
    if [[ -f "$reference_path" ]]; then
        absolute_existing_file "$reference_path"
        return 0
    fi

    reference_path="$(find "$WORKFLOW_ROOT/data" -type f -name "$reference_input" -print -quit)"
    if [[ -n "$reference_path" && -f "$reference_path" ]]; then
        absolute_existing_file "$reference_path"
        return 0
    fi

    return 1
}

# Prompt the user for a value and display a default.
#
# Pressing Enter without typing a value returns the default.
prompt_with_default() {
    local prompt="$1"
    local default="$2"
    local value

    read -r -p "$prompt [$default]: " value
    printf '%s\n' "${value:-$default}"
}

# Prompt the user for an optional value without supplying a default.
prompt_optional() {
    local prompt="$1"
    local value

    read -r -p "$prompt: " value
    printf '%s\n' "$value"
}

# Validate options that must be exactly yes or no.
validate_yes_no() {
    local name="$1"
    local value="$2"

    case "$value" in
        yes|no) ;;
        *)
            echo "[ERROR] $name must be yes or no. Received: $value"
            exit 1
            ;;
    esac
}

# Confirm that an output directory is located under WORKFLOW_ROOT/data.
#
# This prevents the workflow from writing FASTQ or BAM output outside
# the project's standard data directory.
validate_data_output_dir() {
    local name="$1"
    local path="$2"

    case "$path" in
        "$WORKFLOW_ROOT/data"/*) ;;
        *)
            echo "[ERROR] $name must be under $WORKFLOW_ROOT/data"
            echo "[ERROR] Resolved path was: $path"
            exit 1
            ;;
    esac
}

# Confirm that a directory exists on the compute node.
#
# The SLURM job intentionally does not create these directories,
# helping identify project-mounting or path problems.
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

# Load the software modules required for the nanopore workflow.
#
# If environment modules are unavailable, the workflow continues and
# expects dorado, minimap2, and samtools to already be in PATH.
load_nanopore_modules() {
    if command -v module >/dev/null 2>&1; then
        echo "[INFO] Loading nanopore workflow modules: dorado, minimap2, samtools"
        module load dorado
        module load minimap2
        module load samtools
    else
        echo "[WARN] Environment modules are not available in this shell."
        echo "[WARN] Continuing with dorado, minimap2, and samtools from PATH."
    fi
}

# Run the compute-node portion of the nanopore workflow.
#
# This function is called after the script submits itself to SLURM
# using the --run-job argument.
run_job() {
    # Initialize all job inputs and processing parameters.
    local pod5=""
    local reference=""
    local fastq_output_dir=""
    local bam_output_dir=""
    local dorado_model=""
    local device=""
    local min_qscore=""
    local kit_name=""
    local barcode_mode=""
    local modified_bases=""
    local preset=""
    local threads=""
    local secondary=""
    local sort_bam=""
    local index_bam=""
    local min_mapq=""
    local min_read_length=""
    local primary_only=""

    # Remove the initial --run-job argument.
    shift

    # Parse the named parameters passed by submit_workflow.
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --pod5) pod5="$2"; shift 2 ;;
            --workflow-root) shift 2 ;;
            --reference) reference="$2"; shift 2 ;;
            --fastq-output-dir) fastq_output_dir="$2"; shift 2 ;;
            --bam-output-dir) bam_output_dir="$2"; shift 2 ;;
            --dorado-model) dorado_model="$2"; shift 2 ;;
            --device) device="$2"; shift 2 ;;
            --min-qscore) min_qscore="$2"; shift 2 ;;
            --kit-name) kit_name="$2"; shift 2 ;;
            --barcode-mode) barcode_mode="$2"; shift 2 ;;
            --modified-bases) modified_bases="$2"; shift 2 ;;
            --preset) preset="$2"; shift 2 ;;
            --threads) threads="$2"; shift 2 ;;
            --secondary) secondary="$2"; shift 2 ;;
            --sort) sort_bam="$2"; shift 2 ;;
            --index) index_bam="$2"; shift 2 ;;
            --min-mapq) min_mapq="$2"; shift 2 ;;
            --min-read-length) min_read_length="$2"; shift 2 ;;
            --primary-only) primary_only="$2"; shift 2 ;;
            *)
                echo "[ERROR] Unknown job argument: $1"
                exit 1
                ;;
        esac
    done

    # Confirm that the FASTQ and BAM output paths are under data/.
    validate_data_output_dir "FASTQ output directory" "$fastq_output_dir"
    validate_data_output_dir "BAM output directory" "$bam_output_dir"

    # Confirm that all required directories are available on the compute node.
    require_existing_dir "Workflow root" "$WORKFLOW_ROOT"
    require_existing_dir "Log directory" "$LOG_DIR"
    require_existing_dir "FASTQ output directory" "$fastq_output_dir"
    require_existing_dir "BAM output directory" "$bam_output_dir"

    # Prepare Dorado, Minimap2, and samtools.
    load_nanopore_modules

    local pod5_basename
    local pod5_prefix
    local dorado_log
    local minimap2_log
    local fastq_file
    local bam_file

    # Derive the sample name from the POD5 filename.
    #
    # Example:
    # sample.pod5 becomes sample
    pod5_basename="$(basename "$pod5")"
    pod5_prefix="${pod5_basename%.pod5}"

    # Construct separate log files for Dorado and Minimap2.
    #
    # The SLURM job ID is included so repeated runs do not overwrite
    # one another's logs.
    dorado_log="$LOG_DIR/${pod5_prefix}.${SLURM_JOB_ID:-manual}.dorado.log"
    minimap2_log="$LOG_DIR/${pod5_prefix}.${SLURM_JOB_ID:-manual}.minimap2.log"

    # Print the main workflow status and log locations.
    echo "[INFO] Nanopore workflow started: $(date)"
    echo "[INFO] POD5 prefix: $pod5_prefix"
    echo "[INFO] Dorado log: $dorado_log"
    echo "[INFO] Minimap2 log: $minimap2_log"

    # Run the Dorado basecalling script.
    #
    # The dorado_basecall.sh script prints its generated FASTQ path
    # to standard output. Command substitution captures that path in
    # the fastq_file variable.
    fastq_file="$("$SRC_DIR/dorado_basecall.sh" \
        --pod5 "$pod5" \
        --output-dir "$fastq_output_dir" \
        --log-file "$dorado_log" \
        --dorado-model "$dorado_model" \
        --device "$device" \
        --min-qscore "$min_qscore" \
        --kit-name "$kit_name" \
        --barcode-mode "$barcode_mode" \
        --modified-bases "$modified_bases")"

    # Report the FASTQ produced by Dorado.
    echo "[INFO] Dorado complete: $fastq_file"

    # Run the Minimap2 alignment script using the FASTQ generated above.
    #
    # The minimap2_alignment.sh script prints its final sorted and
    # indexed BAM path to standard output. Command substitution captures
    # that path in the bam_file variable.
    bam_file="$("$SRC_DIR/minimap2_alignment.sh" \
        --fastq "$fastq_file" \
        --reference "$reference" \
        --output-dir "$bam_output_dir" \
        --log-file "$minimap2_log" \
        --preset "$preset" \
        --threads "$threads" \
        --secondary "$secondary" \
        --sort "$sort_bam" \
        --index "$index_bam" \
        --min-mapq "$min_mapq" \
        --min-read-length "$min_read_length" \
        --primary-only "$primary_only")"

    # Report the final BAM path and workflow completion time.
    echo "[INFO] Minimap2 complete: $bam_file"
    echo "[INFO] Nanopore workflow complete: $(date)"
}

# Run the interactive submission-node portion of the workflow.
#
# This function:
# 1. Collects the input files and settings.
# 2. Resolves and validates paths.
# 3. Creates output directories.
# 4. Submits the complete workflow to SLURM.
submit_workflow() {
    # Positional arguments may optionally provide the POD5, reference,
    # FASTQ output directory, and BAM output directory.
    local pod5_input="${1:-}"
    local reference_input="${2:-}"
    local fastq_output_input="${3:-}"
    local bam_output_input="${4:-}"

    # Variables that will hold resolved absolute paths.
    local pod5=""
    local reference=""
    local fastq_output_dir
    local bam_output_dir

    # Variables used to construct sample-specific names.
    local pod5_basename
    local pod5_prefix

    # Create the workflow's standard directories on the submission node.
    mkdir -p "$LOG_DIR" "$POD5_DIR" "$DEFAULT_FASTQ_DIR" "$DEFAULT_BAM_DIR"

    # When no POD5 was provided as a positional argument, display the
    # available POD5 files and prompt the user to select one.
    if [[ -z "$pod5_input" ]]; then
        echo "[INFO] Available POD5 files in $POD5_DIR:"
        find "$POD5_DIR" -maxdepth 1 -type f -name "*.pod5" -printf "  %f\n" | sort
        echo
        read -r -p "Enter POD5 filename from data/pod5 or an explicit path: " pod5_input
    fi

    # Prompt for the reference FASTA when it was not supplied.
    if [[ -z "$reference_input" ]]; then
        read -r -p "Enter reference FASTA filename from data/fastq, data/, or an explicit path: " reference_input
    fi

    # Prompt for the FASTQ output directory when it was not supplied.
    if [[ -z "$fastq_output_input" ]]; then
        fastq_output_input="$(prompt_with_default "Enter FASTQ output directory under data/" "$DEFAULT_FASTQ_DIR")"
    fi

    # Prompt for the BAM output directory when it was not supplied.
    if [[ -z "$bam_output_input" ]]; then
        bam_output_input="$(prompt_with_default "Enter BAM output directory under data/" "$DEFAULT_BAM_DIR")"
    fi

    # Confirm that both primary input files were provided.
    if [[ -z "$pod5_input" || -z "$reference_input" ]]; then
        echo "[ERROR] POD5 and reference FASTA are required."
        exit 1
    fi

    # Resolve the POD5 input to an absolute existing path.
    if ! pod5="$(resolve_pod5_file "$pod5_input")"; then
        echo "[ERROR] POD5 file not found: $pod5_input"
        echo "[ERROR] Expected a valid path or a filename in $POD5_DIR"
        exit 1
    fi

    # Resolve the reference FASTA to an absolute existing path.
    if ! reference="$(resolve_reference_file "$reference_input")"; then
        echo "[ERROR] Reference FASTA not found: $reference_input"
        echo "[ERROR] Expected a valid path, a filename in $DEFAULT_FASTQ_DIR, or an exact filename under $WORKFLOW_ROOT/data"
        exit 1
    fi

    # Convert the supplied output-directory values into complete paths.
    fastq_output_dir="$(resolve_data_output_dir "$fastq_output_input")"
    bam_output_dir="$(resolve_data_output_dir "$bam_output_input")"

    # Confirm that both output paths remain under WORKFLOW_ROOT/data.
    validate_data_output_dir "FASTQ output directory" "$fastq_output_dir"
    validate_data_output_dir "BAM output directory" "$bam_output_dir"

    # Derive the sample prefix from the POD5 filename.
    pod5_basename="$(basename "$pod5")"
    pod5_prefix="${pod5_basename%.pod5}"

    echo
    echo "Dorado/basecalling parameters"

    # Collect the Dorado model.
    #
    # The default "sup" requests the super-accuracy model.
    dorado_model="$(prompt_with_default "Enter --dorado-model (sup, hac, fast, or explicit model path)" "sup")"

    # Collect the Dorado compute device.
    #
    # cuda:0 requests the first visible GPU.
    device="$(prompt_with_default "Enter --device" "cuda:0")"

    # Collect the minimum Dorado quality threshold.
    min_qscore="$(prompt_with_default "Enter --min-qscore" "6")"

    # Inform the user that FASTQ output is required because Minimap2
    # receives FASTQ as its input.
    echo "[INFO] --emit-fastq is enabled because the alignment step requires FASTQ input."

    # Determine whether barcode demultiplexing should be performed.
    barcode_mode="$(prompt_with_default "Enter --barcode-mode (none or demux)" "none")"

    # Only request a barcode-kit name when demultiplexing is selected.
    kit_name=""
    if [[ "$barcode_mode" == "demux" ]]; then
        kit_name="$(prompt_optional "Enter --kit-name")"
    fi

    # Optionally collect a Dorado modified-base model or code.
    modified_bases="$(prompt_optional "Enter --modified-bases (optional; press Enter to skip)")"

    # Validate the barcode mode.
    case "$barcode_mode" in
        none|demux) ;;
        *)
            echo "[ERROR] --barcode-mode must be none or demux."
            exit 1
            ;;
    esac

    # Require a sequencing-kit name for demultiplexing.
    if [[ "$barcode_mode" == "demux" && -z "$kit_name" ]]; then
        echo "[ERROR] --kit-name is required when --barcode-mode demux is selected."
        exit 1
    fi

    echo
    echo "Minimap2/alignment parameters"

    # Collect the Minimap2 preset.
    #
    # map-ont is the standard preset for Oxford Nanopore reads.
    preset="$(prompt_with_default "Enter --preset" "map-ont")"

    # Collect the CPU thread count used by alignment and BAM sorting.
    threads="$(prompt_with_default "Enter --threads" "8")"

    # Determine whether secondary alignments should be retained.
    secondary="$(prompt_with_default "Enter --secondary (yes or no)" "yes")"

    # Sorting and indexing are required by the alignment QC steps.
    sort_bam="$(prompt_with_default "Enter --sort (yes required for QC)" "yes")"
    index_bam="$(prompt_with_default "Enter --index (yes required for QC)" "yes")"

    # Collect the mapping-quality and read-length thresholds used by QC.
    min_mapq="$(prompt_with_default "Enter --min-mapq" "20")"
    min_read_length="$(prompt_with_default "Enter --min-read-length" "1000")"

    # Determine whether the output should be filtered to primary
    # alignments only during BAM creation.
    primary_only="$(prompt_with_default "Enter --primary-only (yes or no)" "no")"

    # Validate all yes/no alignment parameters.
    validate_yes_no "--secondary" "$secondary"
    validate_yes_no "--sort" "$sort_bam"
    validate_yes_no "--index" "$index_bam"
    validate_yes_no "--primary-only" "$primary_only"

    # Enforce sorting and indexing because the later QC checks rely
    # on a coordinate-sorted, indexed BAM.
    if [[ "$sort_bam" != "yes" || "$index_bam" != "yes" ]]; then
        echo "[ERROR] This nanopore workflow requires --sort yes and --index yes."
        echo "[ERROR] The requested coverage, MAPQ/read-length, and chromosome 12/rDNA checks run on the sorted indexed BAM."
        exit 1
    fi

    # Create the selected FASTQ and BAM output directories.
    mkdir -p "$fastq_output_dir" "$bam_output_dir"

    # Load the required software on the submission node.
    #
    # This checks that the software environment is available before
    # the workflow is submitted.
    load_nanopore_modules

    # Display the resolved workflow and output locations.
    echo
    echo "[INFO] Workflow root: $WORKFLOW_ROOT"
    echo "[INFO] FASTQ output directory: $fastq_output_dir"
    echo "[INFO] BAM output directory: $bam_output_dir"
    echo
    echo "Submitting nanopore workflow to SLURM..."

    # Submit this same script to SLURM in --run-job mode.
    #
    # --parsable:
    #   Returns a machine-readable SLURM job ID.
    #
    # --job-name:
    #   Creates a sample-specific job name.
    #
    # --chdir:
    #   Runs the job from the workflow root.
    #
    # --output and --error:
    #   Store SLURM standard output and error logs.
    #
    # --export:
    #   Preserves the current environment and explicitly passes the
    #   nanopore workflow root.
    #
    # --cpus-per-task:
    #   Overrides the SBATCH header using the selected thread count.
    #
    # All arguments after SCRIPT_PATH are passed to the compute-node
    # execution of this script.
    job_id="$(sbatch --parsable \
        --job-name="nanopore_${pod5_prefix}" \
        --chdir="$WORKFLOW_ROOT" \
        --output="$LOG_DIR/${pod5_prefix}.%j.slurm.log" \
        --error="$LOG_DIR/${pod5_prefix}.%j.slurm.err" \
        --export=ALL,NANOPORE_WORKFLOW_ROOT="$WORKFLOW_ROOT" \
        --cpus-per-task="$threads" \
        "$SCRIPT_PATH" \
        --run-job \
        --workflow-root "$WORKFLOW_ROOT" \
        --pod5 "$pod5" \
        --reference "$reference" \
        --fastq-output-dir "$fastq_output_dir" \
        --bam-output-dir "$bam_output_dir" \
        --dorado-model "$dorado_model" \
        --device "$device" \
        --min-qscore "$min_qscore" \
        --kit-name "$kit_name" \
        --barcode-mode "$barcode_mode" \
        --modified-bases "$modified_bases" \
        --preset "$preset" \
        --threads "$threads" \
        --secondary "$secondary" \
        --sort "$sort_bam" \
        --index "$index_bam" \
        --min-mapq "$min_mapq" \
        --min-read-length "$min_read_length" \
        --primary-only "$primary_only")"

    # A parsable SLURM job ID may include a cluster name after a semicolon.
    # Keep only the local job ID for constructing the expected log paths.
    log_job_id="${job_id%%;*}"

    # Display the submitted job and all expected output/log paths.
    echo "[INFO] Submitted SLURM job: $job_id"
    echo "[INFO] Expected FASTQ: $fastq_output_dir/${pod5_prefix}.fastq"
    echo "[INFO] Expected BAM:   $bam_output_dir/${pod5_prefix}.sorted.indexed.bam"
    echo "[INFO] SLURM log:      $LOG_DIR/${pod5_prefix}.${log_job_id}.slurm.log"
    echo "[INFO] SLURM err:      $LOG_DIR/${pod5_prefix}.${log_job_id}.slurm.err"
    echo "[INFO] Dorado log:     $LOG_DIR/${pod5_prefix}.${log_job_id}.dorado.log"
    echo "[INFO] Minimap2 log:   $LOG_DIR/${pod5_prefix}.${log_job_id}.minimap2.log"
}

# Main script entry point.
#
# -h or --help:
#   Display usage instructions.
#
# --run-job:
#   Execute the basecalling and alignment steps on the SLURM compute node.
#
# Any other invocation:
#   Collect parameters interactively and submit the workflow to SLURM.
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

if [[ "${1:-}" == "--run-job" ]]; then
    run_job "$@"
else
    submit_workflow "$@"
fi