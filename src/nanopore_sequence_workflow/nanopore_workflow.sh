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

# Search the command-line arguments for optional workflow path values
# before the workflow constructs its standard directory paths.
WORKFLOW_ROOT_ARG=""
WORKFLOW_SRC_DIR_ARG=""
for ((arg_i = 1; arg_i <= $#; arg_i++)); do
    if [[ "${!arg_i}" == "--workflow-root" ]]; then
        next_arg_i=$((arg_i + 1))
        WORKFLOW_ROOT_ARG="${!next_arg_i:-}"
    elif [[ "${!arg_i}" == "--workflow-src-dir" ]]; then
        next_arg_i=$((arg_i + 1))
        WORKFLOW_SRC_DIR_ARG="${!next_arg_i:-}"
    fi
done

# Determine the directory containing this workflow script.
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

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
    WORKFLOW_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
fi

# Define the workflow's script, log, input, and output directories.
#
# In a SLURM job, $0 points to a copied script under /var/spool/slurmd.
# The submission step passes --workflow-src-dir so helper scripts are
# still resolved from the real repository directory.
if [[ -n "$WORKFLOW_SRC_DIR_ARG" ]]; then
    SRC_DIR="$(cd "$WORKFLOW_SRC_DIR_ARG" && pwd)"
elif [[ -n "${NANOPORE_WORKFLOW_SRC_DIR:-}" ]]; then
    SRC_DIR="$NANOPORE_WORKFLOW_SRC_DIR"
else
    SRC_DIR="$SCRIPT_DIR"
fi
SCRIPT_PATH="$SRC_DIR/nanopore_workflow.sh"
LOG_DIR="$WORKFLOW_ROOT/logs/nanopore_sequence_workflow"
POD5_DIR="$WORKFLOW_ROOT/data/pod5"
DEFAULT_FASTQ_DIR="$WORKFLOW_ROOT/data/fastq"
DEFAULT_BAM_DIR="$WORKFLOW_ROOT/data/bam"

# Display the expected positional arguments and default workflow behavior.
usage() {
    echo "Usage: bash src/nanopore_sequence_workflow/nanopore_workflow.sh [POD5] [REFERENCE_FASTA] [FASTQ_OUTPUT_DIR] [BAM_OUTPUT_DIR]"
    echo
    echo "POD5 may be a single .pod5 file or a directory of .pod5 files."
    echo "Relative POD5 inputs are resolved under data/pod5."
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

# Convert an existing file or directory path into an absolute path.
absolute_existing_path() {
    local path="$1"
    local dir
    local file

    dir="$(cd "$(dirname "$path")" && pwd)"
    file="$(basename "$path")"
    printf '%s/%s\n' "$dir" "$file"
}

absolute_existing_file() {
    absolute_existing_path "$1"
}

# Resolve the POD5 input.
#
# The function first checks the input exactly as entered.
# If it is not found, it checks data/pod5.
resolve_pod5_path() {
    local pod5_input="$1"
    local pod5_path

    if [[ -f "$pod5_input" || -d "$pod5_input" ]]; then
        absolute_existing_path "$pod5_input"
        return 0
    fi

    pod5_path="$POD5_DIR/$pod5_input"
    if [[ -f "$pod5_path" || -d "$pod5_path" ]]; then
        absolute_existing_path "$pod5_path"
        return 0
    fi

    return 1
}

pod5_input_has_reads() {
    local pod5_path="$1"

    if [[ -f "$pod5_path" ]]; then
        [[ "$pod5_path" == *.pod5 ]]
        return
    fi

    find "$pod5_path" -maxdepth 1 -type f -name "*.pod5" -print -quit | grep -q .
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

warn_unexpected_input() {
    local received="$1"
    local expected="$2"

    echo "[WARN] Unexpected input: ${received:-<blank>}"
    echo "[WARN] Expected input: $expected"
}

is_yes_no() {
    case "$1" in
        yes|no) return 0 ;;
        *) return 1 ;;
    esac
}

is_optional_positive_integer() {
    [[ -z "$1" || "$1" =~ ^[1-9][0-9]*$ ]]
}

is_positive_integer() {
    [[ "$1" =~ ^[1-9][0-9]*$ ]]
}

is_nonnegative_integer() {
    [[ "$1" =~ ^[0-9]+$ ]]
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

# Validate optional positive-integer options.
#
# A blank value is allowed when the option should be omitted.
validate_optional_positive_integer() {
    local name="$1"
    local value="$2"

    if [[ -n "$value" && ! "$value" =~ ^[1-9][0-9]*$ ]]; then
        echo "[ERROR] $name must be a positive integer or blank. Received: $value"
        exit 1
    fi
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
    local output_format="${1:-fastq}"

    if command -v module >/dev/null 2>&1; then
        if [[ "$output_format" == "fastq" ]]; then
            echo "[INFO] Loading nanopore workflow modules: dorado, minimap2, samtools"
        else
            echo "[INFO] Loading nanopore workflow modules: dorado, samtools"
        fi
        module load dorado
        if [[ "$output_format" == "fastq" ]]; then
            module load minimap2
        fi
        module load samtools
    else
        echo "[WARN] Environment modules are not available in this shell."
        if [[ "$output_format" == "fastq" ]]; then
            echo "[WARN] Continuing with dorado, minimap2, and samtools from PATH."
        else
            echo "[WARN] Continuing with dorado and samtools from PATH."
        fi
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
    local max_reads=""
    local output_format=""
    local emit_moves=""
    local kit_name=""
    local barcode_mode=""
    local barcode_output_dir=""
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
            --workflow-src-dir) shift 2 ;;
            --reference) reference="$2"; shift 2 ;;
            --fastq-output-dir) fastq_output_dir="$2"; shift 2 ;;
            --bam-output-dir) bam_output_dir="$2"; shift 2 ;;
            --dorado-model) dorado_model="$2"; shift 2 ;;
            --device) device="$2"; shift 2 ;;
            --min-qscore) min_qscore="$2"; shift 2 ;;
            --max-reads) max_reads="$2"; shift 2 ;;
            --output-format) output_format="$2"; shift 2 ;;
            --emit-moves) emit_moves="$2"; shift 2 ;;
            --kit-name) kit_name="$2"; shift 2 ;;
            --barcode-mode) barcode_mode="$2"; shift 2 ;;
            --barcode-output-dir) barcode_output_dir="$2"; shift 2 ;;
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

    output_format="${output_format:-fastq}"
    emit_moves="${emit_moves:-yes}"
    barcode_output_dir="${barcode_output_dir:-$DEFAULT_BAM_DIR}"

    case "$output_format" in
        fastq|bam) ;;
        *)
            echo "[ERROR] Dorado output must be fastq or bam."
            exit 1
            ;;
    esac

    validate_yes_no "--emit-moves" "$emit_moves"
    validate_optional_positive_integer "--max-reads" "$max_reads"

    # Confirm that the selected output paths are under data/.
    if [[ "$output_format" == "fastq" ]]; then
        validate_data_output_dir "FASTQ output directory" "$fastq_output_dir"
    fi
    validate_data_output_dir "BAM output directory" "$bam_output_dir"
    if [[ "$barcode_mode" == "demux" ]]; then
        validate_data_output_dir "Barcode output directory" "$barcode_output_dir"
    fi

    # Confirm that all required directories are available on the compute node.
    require_existing_dir "Workflow root" "$WORKFLOW_ROOT"
    require_existing_dir "Log directory" "$LOG_DIR"
    if [[ "$output_format" == "fastq" ]]; then
        require_existing_dir "FASTQ output directory" "$fastq_output_dir"
    fi
    require_existing_dir "BAM output directory" "$bam_output_dir"
    if [[ "$barcode_mode" == "demux" ]]; then
        require_existing_dir "Barcode output directory" "$barcode_output_dir"
    fi

    # Prepare the software required for the selected workflow branch.
    load_nanopore_modules "$output_format"

    local pod5_basename
    local pod5_prefix
    local dorado_log
    local alignment_log
    local dorado_output
    local fastq_file
    local bam_file
    local dorado_output_dir
    local dorado_mm2_opts
    local barcode_dir
    local barcode_file
    local barcode_label
    local barcode_bam
    local barcode_outputs=()
    local barcode_files=()

    # Derive the sample name from the POD5 file or directory name.
    #
    # Example:
    # sample.pod5 becomes sample
    pod5_basename="$(basename "$pod5")"
    pod5_prefix="${pod5_basename%.pod5}"

    # Construct separate log files for Dorado and alignment/QC.
    #
    # The SLURM job ID is included so repeated runs do not overwrite
    # one another's logs.
    dorado_log="$LOG_DIR/${pod5_prefix}.${SLURM_JOB_ID:-manual}.dorado.log"
    if [[ "$barcode_mode" == "demux" ]]; then
        alignment_log="$LOG_DIR/${pod5_prefix}.${SLURM_JOB_ID:-manual}.barcode_alignment.log"
    elif [[ "$output_format" == "bam" ]]; then
        alignment_log="$LOG_DIR/${pod5_prefix}.${SLURM_JOB_ID:-manual}.dorado_alignment.log"
    else
        alignment_log="$LOG_DIR/${pod5_prefix}.${SLURM_JOB_ID:-manual}.minimap2.log"
    fi

    # Print the main workflow status and log locations.
    echo "[INFO] Nanopore workflow started: $(date)"
    echo "[INFO] POD5 prefix: $pod5_prefix"
    echo "[INFO] Dorado log: $dorado_log"
    echo "[INFO] Alignment/QC log: $alignment_log"
    echo "[INFO] Dorado output format: $output_format"

    dorado_output_dir="$fastq_output_dir"
    if [[ "$output_format" == "bam" ]]; then
        dorado_output_dir="$bam_output_dir"
    fi

    # Dorado exposes supported minimap2 alignment settings through
    # --mm2-opts. The workflow's sort, index, MAPQ/read-length, and
    # primary-only settings are handled by samtools after Dorado writes BAM.
    dorado_mm2_opts=""
    if [[ "$output_format" == "bam" ]]; then
        dorado_mm2_opts="-x $preset"
        if [[ "$secondary" == "yes" ]]; then
            dorado_mm2_opts="$dorado_mm2_opts --secondary"
        fi
    fi

    # Run the Dorado basecalling script.
    #
    # The dorado_basecall.sh script prints its generated output path
    # to standard output. Command substitution captures that path in
    # the dorado_output variable.
    dorado_output="$("$SRC_DIR/dorado_basecall.sh" \
        --pod5 "$pod5" \
        --output-dir "$dorado_output_dir" \
        --log-file "$dorado_log" \
        --dorado-model "$dorado_model" \
        --device "$device" \
        --min-qscore "$min_qscore" \
        --max-reads "$max_reads" \
        --output-format "$output_format" \
        --reference "$reference" \
        --emit-moves "$emit_moves" \
        --mm2-opts "$dorado_mm2_opts" \
        --kit-name "$kit_name" \
        --barcode-mode "$barcode_mode" \
        --barcode-output-dir "$barcode_output_dir" \
        --modified-bases "$modified_bases")"

    echo "[INFO] Dorado complete: $dorado_output"

    if [[ "$barcode_mode" == "demux" ]]; then
        barcode_dir="$barcode_output_dir/${pod5_prefix}_${SLURM_JOB_ID:-manual}"
        echo "[INFO] Aligning demuxed barcode files in: $barcode_dir"

        if [[ "$output_format" == "fastq" ]]; then
            for barcode_file in "$barcode_dir"/*.fastq "$barcode_dir"/*.fq; do
                [[ -f "$barcode_file" ]] || continue
                barcode_files+=("$barcode_file")
            done
        else
            for barcode_file in "$barcode_dir"/*.bam; do
                [[ -f "$barcode_file" ]] || continue
                barcode_files+=("$barcode_file")
            done
        fi

        if [[ "${#barcode_files[@]}" -eq 0 ]]; then
            echo "[ERROR] No demuxed barcode files found in: $barcode_dir"
            exit 1
        fi

        {
            echo "========================================="
            echo "  BARCODE ALIGNMENT/QC RUN"
            echo "========================================="
            echo "Started:       $(date)"
            echo "Barcode dir:   $barcode_dir"
            echo "Output format: $output_format"
            echo "Reference:     $reference"
            echo "========================================="
        } > "$alignment_log"

        for barcode_file in "${barcode_files[@]}"; do
            barcode_label="$(basename "$barcode_file")"
            barcode_label="${barcode_label%.fastq}"
            barcode_label="${barcode_label%.fq}"
            barcode_label="${barcode_label%.bam}"

            echo "[INFO] Processing barcode file: $barcode_file"

            if [[ "$output_format" == "fastq" ]]; then
                barcode_bam="$("$SRC_DIR/minimap2_alignment.sh" \
                    --fastq "$barcode_file" \
                    --reference "$reference" \
                    --output-dir "$barcode_dir" \
                    --log-file "$alignment_log" \
                    --preset "$preset" \
                    --threads "$threads" \
                    --secondary "$secondary" \
                    --sort "$sort_bam" \
                    --index "$index_bam" \
                    --min-mapq "$min_mapq" \
                    --min-read-length "$min_read_length" \
                    --primary-only "$primary_only" \
                    --append-log yes \
                    --log-label "$barcode_label")"
            else
                barcode_bam="$("$SRC_DIR/minimap2_alignment.sh" \
                    --input-bam "$barcode_file" \
                    --output-dir "$barcode_dir" \
                    --log-file "$alignment_log" \
                    --preset "$preset" \
                    --threads "$threads" \
                    --secondary "$secondary" \
                    --sort "$sort_bam" \
                    --index "$index_bam" \
                    --min-mapq "$min_mapq" \
                    --min-read-length "$min_read_length" \
                    --primary-only "$primary_only" \
                    --append-log yes \
                    --log-label "$barcode_label")"
            fi

            barcode_outputs+=("$barcode_bam")
            echo "[INFO] Barcode alignment/QC complete: $barcode_bam"
        done

        {
            echo ""
            echo "========================================="
            echo "  BARCODE ALIGNMENT/QC COMPLETE"
            echo "========================================="
            echo "Completed: $(date)"
            echo "Final BAMs:"
            printf '%s\n' "${barcode_outputs[@]}"
            echo "========================================="
        } >> "$alignment_log"
    elif [[ "$output_format" == "fastq" ]]; then
        fastq_file="$dorado_output"

        # Run the Minimap2 alignment script using the FASTQ generated above.
        #
        # The minimap2_alignment.sh script prints its final sorted and
        # indexed BAM path to standard output. Command substitution captures
        # that path in the bam_file variable.
        bam_file="$("$SRC_DIR/minimap2_alignment.sh" \
            --fastq "$fastq_file" \
            --reference "$reference" \
            --output-dir "$bam_output_dir" \
            --log-file "$alignment_log" \
            --preset "$preset" \
            --threads "$threads" \
            --secondary "$secondary" \
            --sort "$sort_bam" \
            --index "$index_bam" \
            --min-mapq "$min_mapq" \
            --min-read-length "$min_read_length" \
            --primary-only "$primary_only")"

        echo "[INFO] Minimap2 complete: $bam_file"
    else
        # Sort, index, optionally filter, and QC the Dorado-aligned BAM.
        bam_file="$("$SRC_DIR/minimap2_alignment.sh" \
            --input-bam "$dorado_output" \
            --output-dir "$bam_output_dir" \
            --log-file "$alignment_log" \
            --preset "$preset" \
            --threads "$threads" \
            --secondary "$secondary" \
            --sort "$sort_bam" \
            --index "$index_bam" \
            --min-mapq "$min_mapq" \
            --min-read-length "$min_read_length" \
            --primary-only "$primary_only")"

        echo "[INFO] Dorado alignment QC complete: $bam_file"
    fi

    # Report workflow completion time.
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
    local output_format
    local emit_moves
    local run_demux
    local max_reads

    # Variables used to construct sample-specific names.
    local pod5_basename
    local pod5_prefix

    # Create the workflow's standard directories on the submission node.
    mkdir -p "$LOG_DIR" "$POD5_DIR" "$DEFAULT_FASTQ_DIR" "$DEFAULT_BAM_DIR"

    # When no POD5 was provided as a positional argument, display the
    # available POD5 files/directories and prompt the user to select one.
    while true; do
        if [[ -z "$pod5_input" ]]; then
            echo "[INFO] Available POD5 files and directories in $POD5_DIR:"
            find "$POD5_DIR" -maxdepth 1 -type f -name "*.pod5" -printf "  %f\n" | sort
            find "$POD5_DIR" -mindepth 1 -maxdepth 1 -type d -printf "  %f/\n" | sort
            echo
            read -r -p "Enter POD5 file/directory from data/pod5 or an explicit path: " pod5_input
        fi

        if [[ -z "$pod5_input" ]]; then
            warn_unexpected_input "$pod5_input" "a .pod5 file, a directory containing .pod5 files, or a name in $POD5_DIR"
            continue
        fi

        if pod5="$(resolve_pod5_path "$pod5_input")" && pod5_input_has_reads "$pod5"; then
            break
        fi

        warn_unexpected_input "$pod5_input" "a valid .pod5 file/directory or a name in $POD5_DIR"
        pod5_input=""
    done

    # Prompt for the reference FASTA when it was not supplied.
    while true; do
        if [[ -z "$reference_input" ]]; then
            read -r -p "Enter reference FASTA filename from data/fastq, data/, or an explicit path: " reference_input
        fi

        if [[ -z "$reference_input" ]]; then
            warn_unexpected_input "$reference_input" "a FASTA filename or path, for example reference.fna"
            continue
        fi

        if reference="$(resolve_reference_file "$reference_input")"; then
            break
        fi

        warn_unexpected_input "$reference_input" "a valid path, a filename in $DEFAULT_FASTQ_DIR, or an exact filename under $WORKFLOW_ROOT/data"
        reference_input=""
    done

    echo
    echo "Dorado output selection"
    while true; do
        output_format="$(prompt_with_default "Should Dorado produce fastq or bam?" "bam")"
        case "$output_format" in
            fastq|bam) break ;;
            *) warn_unexpected_input "$output_format" "fastq or bam, for example bam" ;;
        esac
    done

    # Prompt for the FASTQ output directory only when FASTQ is selected.
    if [[ "$output_format" == "fastq" && -z "$fastq_output_input" ]]; then
        fastq_output_input="$(prompt_with_default "Enter FASTQ output directory under data/" "$DEFAULT_FASTQ_DIR")"
    elif [[ -z "$fastq_output_input" ]]; then
        fastq_output_input="$DEFAULT_FASTQ_DIR"
    fi

    # Prompt for the BAM output directory when it was not supplied.
    if [[ -z "$bam_output_input" ]]; then
        bam_output_input="$(prompt_with_default "Enter BAM output directory under data/" "$DEFAULT_BAM_DIR")"
    fi

    # Convert the supplied output-directory values into complete paths.
    while true; do
        fastq_output_dir="$(resolve_data_output_dir "$fastq_output_input")"
        if [[ "$output_format" != "fastq" || "$fastq_output_dir" == "$WORKFLOW_ROOT/data"/* ]]; then
            break
        fi
        warn_unexpected_input "$fastq_output_input" "a FASTQ output directory under $WORKFLOW_ROOT/data"
        fastq_output_input="$(prompt_with_default "Enter FASTQ output directory under data/" "$DEFAULT_FASTQ_DIR")"
    done

    while true; do
        bam_output_dir="$(resolve_data_output_dir "$bam_output_input")"
        if [[ "$bam_output_dir" == "$WORKFLOW_ROOT/data"/* ]]; then
            break
        fi
        warn_unexpected_input "$bam_output_input" "a BAM output directory under $WORKFLOW_ROOT/data"
        bam_output_input="$(prompt_with_default "Enter BAM output directory under data/" "$DEFAULT_BAM_DIR")"
    done

    # Derive the sample prefix from the POD5 file or directory name.
    pod5_basename="$(basename "$pod5")"
    pod5_prefix="${pod5_basename%.pod5}"

    echo
    echo "Dorado demultiplexing"

    # Ask about Dorado demux before collecting the rest of the
    # basecalling parameters so the barcode kit can be supplied to
    # basecalling for barcode classification.
    while true; do
        run_demux="$(prompt_with_default "Run dorado demux after basecalling? (yes or no)" "no")"
        if is_yes_no "$run_demux"; then
            break
        fi
        warn_unexpected_input "$run_demux" "yes or no"
    done

    barcode_mode="none"
    kit_name=""
    if [[ "$run_demux" == "yes" ]]; then
        barcode_mode="demux"
        while [[ -z "$kit_name" ]]; do
            kit_name="$(prompt_optional "Enter barcode-kit / --kit-name")"
            if [[ -z "$kit_name" ]]; then
                warn_unexpected_input "$kit_name" "a barcode kit name, for example SQK-NBD114-24"
            fi
        done
    fi

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
    while true; do
        min_qscore="$(prompt_with_default "Enter --min-qscore" "6")"
        if is_nonnegative_integer "$min_qscore"; then
            break
        fi
        warn_unexpected_input "$min_qscore" "a non-negative integer, for example 6"
    done

    # Limit basecalling to the first N reads when supplied.
    #
    # Leaving this blank omits --max-reads so Dorado basecalls all reads.
    while true; do
        max_reads="$(prompt_optional "Enter --max-reads (optional; press Enter to basecall all reads)")"
        if is_optional_positive_integer "$max_reads"; then
            break
        fi
        warn_unexpected_input "$max_reads" "a positive integer, for example 1000, or blank"
    done

    # Move tables are retained by BAM output and are useful for
    # downstream signal-aware tools. FASTQ cannot store this metadata.
    while true; do
        emit_moves="$(prompt_with_default "Enter --emit-moves (yes or no; BAM output only)" "yes")"
        if is_yes_no "$emit_moves"; then
            break
        fi
        warn_unexpected_input "$emit_moves" "yes or no"
    done

    if [[ "$output_format" == "fastq" ]]; then
        echo "[INFO] --emit-fastq is enabled because the FASTQ path uses Minimap2 alignment."
        echo "[INFO] --emit-moves is ignored for FASTQ output because FASTQ cannot store move-table tags."
    else
        echo "[INFO] Dorado will basecall and align with --reference, producing BAM output."
    fi

    # Optionally collect a Dorado modified-base model or code.
    modified_bases="$(prompt_optional "Enter --modified-bases (optional; press Enter to skip)")"

    echo
    echo "Alignment and QC parameters"

    # Collect the minimap2/Dorado alignment preset.
    #
    # map-ont is the standard preset for Oxford Nanopore reads.
    preset="$(prompt_with_default "Enter --preset" "map-ont")"

    # Collect the CPU thread count used by alignment and BAM sorting.
    while true; do
        threads="$(prompt_with_default "Enter --threads" "8")"
        if is_positive_integer "$threads"; then
            break
        fi
        warn_unexpected_input "$threads" "a positive integer, for example 8"
    done

    # Determine whether secondary alignments should be retained.
    while true; do
        secondary="$(prompt_with_default "Enter --secondary (yes or no)" "no")"
        if is_yes_no "$secondary"; then
            break
        fi
        warn_unexpected_input "$secondary" "yes or no"
    done

    # Sorting and indexing are required by the alignment QC steps.
    while true; do
        sort_bam="$(prompt_with_default "Enter --sort (yes required for QC)" "yes")"
        if [[ "$sort_bam" == "yes" ]]; then
            break
        fi
        warn_unexpected_input "$sort_bam" "yes because this workflow requires sorted BAM output for QC"
    done
    while true; do
        index_bam="$(prompt_with_default "Enter --index (yes required for QC)" "yes")"
        if [[ "$index_bam" == "yes" ]]; then
            break
        fi
        warn_unexpected_input "$index_bam" "yes because this workflow requires indexed BAM output for QC"
    done

    # Collect the mapping-quality and read-length thresholds used by QC.
    while true; do
        min_mapq="$(prompt_with_default "Enter --min-mapq" "20")"
        if is_nonnegative_integer "$min_mapq"; then
            break
        fi
        warn_unexpected_input "$min_mapq" "a non-negative integer, for example 20"
    done
    while true; do
        min_read_length="$(prompt_with_default "Enter --min-read-length" "1000")"
        if is_positive_integer "$min_read_length"; then
            break
        fi
        warn_unexpected_input "$min_read_length" "a positive integer, for example 1000"
    done

    # Determine whether the output should be filtered to primary
    # alignments only during BAM creation.
    while true; do
        primary_only="$(prompt_with_default "Enter --primary-only (yes or no)" "yes")"
        if is_yes_no "$primary_only"; then
            break
        fi
        warn_unexpected_input "$primary_only" "yes or no"
    done

    # Create the selected output directories.
    if [[ "$output_format" == "fastq" ]]; then
        mkdir -p "$fastq_output_dir"
    fi
    mkdir -p "$bam_output_dir"

    # Load the required software on the submission node.
    #
    # This checks that the software environment is available before
    # the workflow is submitted.
    load_nanopore_modules "$output_format"

    # Display the resolved workflow and output locations.
    echo
    echo "[INFO] Workflow root: $WORKFLOW_ROOT"
    echo "[INFO] Dorado output format: $output_format"
    if [[ "$output_format" == "fastq" ]]; then
        echo "[INFO] FASTQ output directory: $fastq_output_dir"
    fi
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
    #   nanopore workflow root and source-script directory.
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
        --export=ALL,NANOPORE_WORKFLOW_ROOT="$WORKFLOW_ROOT",NANOPORE_WORKFLOW_SRC_DIR="$SRC_DIR" \
        --cpus-per-task="$threads" \
        "$SCRIPT_PATH" \
        --run-job \
        --workflow-root "$WORKFLOW_ROOT" \
        --workflow-src-dir "$SRC_DIR" \
        --pod5 "$pod5" \
        --reference "$reference" \
        --fastq-output-dir "$fastq_output_dir" \
        --bam-output-dir "$bam_output_dir" \
        --dorado-model "$dorado_model" \
        --device "$device" \
        --min-qscore "$min_qscore" \
        --max-reads "$max_reads" \
        --output-format "$output_format" \
        --emit-moves "$emit_moves" \
        --kit-name "$kit_name" \
        --barcode-mode "$barcode_mode" \
        --barcode-output-dir "$DEFAULT_BAM_DIR" \
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
    if [[ "$output_format" == "fastq" ]]; then
        echo "[INFO] Expected FASTQ: $fastq_output_dir/${pod5_prefix}_${log_job_id}.fastq"
        echo "[INFO] Expected raw BAM: $bam_output_dir/${pod5_prefix}_${log_job_id}.bam"
    else
        echo "[INFO] Expected Dorado aligned BAM: $bam_output_dir/${pod5_prefix}_${log_job_id}.bam"
        echo "[INFO] Expected raw BAM: $bam_output_dir/${pod5_prefix}_${log_job_id}.dorado.raw.bam"
    fi
    if [[ "$barcode_mode" == "demux" ]]; then
        echo "[INFO] Expected barcode output directory: $DEFAULT_BAM_DIR/${pod5_prefix}_${log_job_id}"
        echo "[INFO] Expected barcode alignment/QC log: $LOG_DIR/${pod5_prefix}.${log_job_id}.barcode_alignment.log"
    else
        echo "[INFO] Expected sorted indexed BAM: $bam_output_dir/${pod5_prefix}.sorted.indexed_${log_job_id}.bam"
    fi
    echo "[INFO] SLURM log:      $LOG_DIR/${pod5_prefix}.${log_job_id}.slurm.log"
    echo "[INFO] SLURM err:      $LOG_DIR/${pod5_prefix}.${log_job_id}.slurm.err"
    echo "[INFO] Dorado log:     $LOG_DIR/${pod5_prefix}.${log_job_id}.dorado.log"
    if [[ "$barcode_mode" == "demux" ]]; then
        echo "[INFO] Alignment/QC log: $LOG_DIR/${pod5_prefix}.${log_job_id}.barcode_alignment.log"
    elif [[ "$output_format" == "fastq" ]]; then
        echo "[INFO] Alignment/QC log: $LOG_DIR/${pod5_prefix}.${log_job_id}.minimap2.log"
    else
        echo "[INFO] Alignment/QC log: $LOG_DIR/${pod5_prefix}.${log_job_id}.dorado_alignment.log"
    fi
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
