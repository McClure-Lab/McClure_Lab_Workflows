#!/bin/bash
#SBATCH --job-name=nanopore_workflow
#SBATCH --partition=gpu
#SBATCH --account=gpu_rbi
#SBATCH --gpus=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
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
POD5_DIR="$WORKFLOW_ROOT/data/pod5"

# Reference FASTA files are still searched for under a fixed location,
# independent of wherever a particular POD5 input happens to live.
# (FASTQ/BAM/log *output* locations, by contrast, are derived from the
# POD5 input's own location -- see derive_pod5_data_root() below. LOG_DIR
# itself is set later, in submit_workflow/run_job, once the POD5 input has
# been resolved and that root is known.)
REFERENCE_SEARCH_DIR="$WORKFLOW_ROOT/data/fastq"

# Display the expected positional arguments and default workflow behavior.
usage() {
    echo "Usage: bash src/nanopore_sequence_workflow/nanopore_workflow.sh [POD5] [REFERENCE_FASTA] [FASTQ_OUTPUT_DIR] [BAM_OUTPUT_DIR]"
    echo
    echo "POD5 may be a single .pod5 file, a directory of .pod5 files, or a .zip archive of .pod5 files."
    echo "A directory may also contain a .zip archive, or nested subfolders of .pod5 files, instead of"
    echo "flat .pod5 files directly inside it -- any of these are found and normalized automatically."
    echo "Relative POD5 inputs are resolved under data/pod5."
    echo
    echo "FASTQ, BAM, and log outputs default to folders created one directory above the POD5 input, e.g.:"
    echo "  <one directory above POD5>/data/fastq"
    echo "  <one directory above POD5>/data/bam_aligned"
    echo "  <one directory above POD5>/data/bam_basecalled   (only used with --output-format both)"
    echo "  <one directory above POD5>/logs/nanopore_sequence_workflow"
    echo
    echo "Dorado output format may be fastq, bam (aligned only), or both (basecalled BAM + derived FASTQ + aligned BAM)."
    echo "Note: barcode demultiplexing always uses --output-format both internally, so every barcode gets both a"
    echo "demuxed FASTQ and a demuxed aligned BAM automatically."
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
    local root="$2"

    if [[ "$path" = /* ]]; then
        printf '%s\n' "$path"
    elif [[ "$path" == data/* ]]; then
        printf '%s/%s\n' "$root" "$path"
    else
        printf '%s/data/%s\n' "$root" "$path"
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

# Extract a .zip archive of POD5 files into a sibling directory (named
# after the archive, minus .zip) and return that directory's path.
#
# An existing extraction is reused rather than re-extracted, so running
# the workflow again against the same archive does not redo the work.
extract_pod5_zip() {
    local zip_path="$1"
    local extract_dir="${zip_path%.zip}"

    if [[ ! -d "$extract_dir" ]]; then
        if ! command -v unzip >/dev/null 2>&1; then
            echo "[ERROR] 'unzip' is required to extract $zip_path but was not found in PATH." >&2
            exit 1
        fi
        echo "[INFO] Extracting POD5 archive: $zip_path -> $extract_dir" >&2
        mkdir -p "$extract_dir"
        unzip -q -o "$zip_path" -d "$extract_dir"
    else
        echo "[INFO] Reusing already-extracted POD5 directory: $extract_dir" >&2
    fi

    printf '%s\n' "$extract_dir"
}

# Normalize a POD5 directory so Dorado can be pointed straight at it.
#
# Handles, in any combination, at any depth under the given directory:
# - .pod5 files already sitting directly in it (the simple, common case),
# - .pod5 files nested in subfolders (e.g. a zip that unpacked into its
#   own named folder instead of unpacking flat),
# - .zip archives found anywhere under it (each is extracted in place via
#   extract_pod5_zip).
#
# If everything is already flat (only top-level .pod5 files, no nested
# .pod5, no .zip anywhere below), the directory is returned unchanged.
# Otherwise every .pod5 file discovered (original or freshly extracted)
# is symlinked into a flat staging directory ("<dir>_pod5_staged"), and
# that staging directory is returned instead, so the caller always ends
# up with something Dorado can read directly.
normalize_pod5_directory() {
    local dir="$1"
    local top_level_pod5_count
    local nested_pod5_count
    local zip_count
    local zip_file
    local pod5_file
    local staged_dir

    top_level_pod5_count=$(find "$dir" -maxdepth 1 -type f -name "*.pod5" | wc -l)
    nested_pod5_count=$(find "$dir" -mindepth 2 -type f -name "*.pod5" | wc -l)
    zip_count=$(find "$dir" -type f -name "*.zip" | wc -l)

    if [[ "$top_level_pod5_count" -gt 0 && "$nested_pod5_count" -eq 0 && "$zip_count" -eq 0 ]]; then
        printf '%s\n' "$dir"
        return 0
    fi

    # Extract every .zip found anywhere under this directory (each next to
    # itself; extract_pod5_zip reuses an existing extraction if present).
    while IFS= read -r zip_file; do
        extract_pod5_zip "$zip_file" >/dev/null
    done < <(find "$dir" -type f -name "*.zip")

    # Collect every .pod5 file now present anywhere under this directory,
    # originals and anything just extracted from a zip alike.
    local pod5_files=()
    while IFS= read -r pod5_file; do
        pod5_files+=("$pod5_file")
    done < <(find "$dir" -type f -name "*.pod5")

    if [[ "${#pod5_files[@]}" -eq 0 ]]; then
        return 1
    fi

    staged_dir="${dir%/}_pod5_staged"
    mkdir -p "$staged_dir"
    for pod5_file in "${pod5_files[@]}"; do
        ln -sf "$pod5_file" "$staged_dir/$(basename "$pod5_file")"
    done

    printf '%s\n' "$staged_dir"
    return 0
}

# Resolve the POD5 input to an existing file or directory path, WITHOUT
# extracting any .zip archive or otherwise touching its contents.
#
# This is the only POD5 resolution step meant to run on a login/submission
# node: it just confirms *something* exists at the given path (checked as
# entered, then under data/pod5), which is cheap regardless of how large
# the eventual POD5 data turns out to be. Actual extraction/normalization
# (extract_pod5_zip / normalize_pod5_directory, below) is deliberately
# deferred to run_job, so it only ever happens on a compute node, after
# the SLURM job has been submitted -- not while a user is sitting at an
# interactive prompt on the login node.
resolve_pod5_input_path() {
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

# A cheap, non-extracting plausibility check: does this resolved path
# look like it could contain POD5 data? Only inspects filenames (find
# never opens/decompresses anything), so it stays safe to run on a login
# node even against a very large POD5 directory tree. A single file must
# itself be .pod5 or .zip; a directory must contain a .pod5 or .zip
# *somewhere* under it (any depth), matching what normalize_pod5_directory
# is later able to find and extract on the compute node.
pod5_input_looks_valid() {
    local pod5_path="$1"

    if [[ -f "$pod5_path" ]]; then
        [[ "$pod5_path" == *.pod5 || "$pod5_path" == *.zip ]]
        return
    fi

    find "$pod5_path" -type f \( -name "*.pod5" -o -name "*.zip" \) -print -quit | grep -q .
}

# Fully resolve the POD5 input to something Dorado can be pointed at
# directly: extracts a .zip archive given directly (extract_pod5_zip),
# and/or normalizes a directory that contains nested .pod5 files or .zip
# archives anywhere under it (normalize_pod5_directory). This does real
# extraction I/O and is only ever called from run_job, i.e. on the
# compute node, never during the interactive submission prompts.
resolve_pod5_path() {
    local pod5_input="$1"
    local resolved

    resolved="$(resolve_pod5_input_path "$pod5_input")" || return 1

    if [[ -f "$resolved" && "$resolved" == *.zip ]]; then
        resolved="$(extract_pod5_zip "$resolved")"
    fi

    if [[ -d "$resolved" ]]; then
        resolved="$(normalize_pod5_directory "$resolved")" || return 1
    fi

    printf '%s\n' "$resolved"
    return 0
}

# Determine the "data root" used for default FASTQ/BAM output locations:
# the directory one level above wherever the resolved POD5 input lives.
#
# Example:
# /projects/run1/pod5/sample.pod5  ->  pod5 directory: /projects/run1/pod5
#                                   ->  data root:      /projects/run1
#
# Works the same whether the resolved POD5 input is a file or a directory
# (including a directory produced by extract_pod5_zip above).
derive_pod5_data_root() {
    local resolved_pod5="$1"
    local pod5_container_dir="$resolved_pod5"

    if [[ -f "$resolved_pod5" ]]; then
        pod5_container_dir="$(dirname "$resolved_pod5")"
    fi

    dirname "$pod5_container_dir"
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

    reference_path="$REFERENCE_SEARCH_DIR/$reference_input"
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
    local allowed_root="$3"

    case "$path" in
        "$allowed_root"/*) ;;
        *)
            echo "[ERROR] $name must be under $allowed_root"
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
    local needs_minimap2="no"

    # fastq aligns via minimap2; both derives a FASTQ from the basecalled
    # BAM and also aligns it via minimap2. Plain bam mode has Dorado align
    # internally and does not need minimap2 loaded separately.
    if [[ "$output_format" == "fastq" || "$output_format" == "both" ]]; then
        needs_minimap2="yes"
    fi

    if command -v module >/dev/null 2>&1; then
        if [[ "$needs_minimap2" == "yes" ]]; then
            echo "[INFO] Loading nanopore workflow modules: dorado, minimap2, samtools"
        else
            echo "[INFO] Loading nanopore workflow modules: dorado, samtools"
        fi
        module load dorado
        if [[ "$needs_minimap2" == "yes" ]]; then
            module load minimap2
        fi
        module load samtools
    else
        echo "[WARN] Environment modules are not available in this shell."
        if [[ "$needs_minimap2" == "yes" ]]; then
            echo "[WARN] Continuing with dorado, minimap2, and samtools from PATH."
        else
            echo "[WARN] Continuing with dorado and samtools from PATH."
        fi
    fi
}

# Print the barcode kit names Dorado currently recognizes, parsed
# best-effort from 'dorado demux --help'. This deliberately does not fall
# back to a hardcoded kit list -- ONT's supported kits change over time,
# and a stale hardcoded list would be worse than no list. If dorado isn't
# on PATH yet, or its help output can't be parsed, this just warns and
# lets the user type a kit name manually.
list_available_barcode_kits() {
    if ! command -v dorado >/dev/null 2>&1; then
        echo "[WARN] 'dorado' was not found in PATH; cannot list barcode kits automatically."
        echo "[WARN] Load the dorado module first, or check the kit name against Oxford Nanopore's documentation."
        return
    fi

    local kits
    kits="$(dorado demux --help 2>&1 | grep -oE '[A-Z]+[A-Z0-9]*-[A-Z0-9-]*[0-9]+' | sort -u || true)"

    if [[ -z "$kits" ]]; then
        echo "[WARN] Could not automatically determine barcode kit names from 'dorado demux --help'."
        echo "[WARN] Run 'dorado demux --help' manually, or check Oxford Nanopore's documentation for valid --kit-name values."
        return
    fi

    echo "[INFO] Available barcode kits:"
    printf '  %s\n' $kits
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
    local basecalled_bam_output_dir=""
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
            --basecalled-bam-output-dir) basecalled_bam_output_dir="$2"; shift 2 ;;
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

    if [[ -z "$pod5" ]]; then
        echo "[ERROR] --pod5 is required."
        exit 1
    fi

    local pod5_data_root
    pod5_data_root="$(derive_pod5_data_root "$pod5")"

    # Logs live one directory above the POD5 input too, alongside the
    # data/ output directories.
    LOG_DIR="$pod5_data_root/logs/nanopore_sequence_workflow"

    output_format="${output_format:-fastq}"
    emit_moves="${emit_moves:-yes}"
    if [[ "$output_format" == "both" ]]; then
        # Demux against a "both" run splits the basecalled (unaligned) BAM,
        # so barcode output belongs alongside it, not the aligned BAM dir.
        barcode_output_dir="${barcode_output_dir:-$basecalled_bam_output_dir}"
    else
        barcode_output_dir="${barcode_output_dir:-$pod5_data_root/data/bam_aligned}"
    fi

    case "$output_format" in
        fastq|bam|both) ;;
        *)
            echo "[ERROR] Dorado output must be fastq, bam, or both."
            exit 1
            ;;
    esac

    validate_yes_no "--emit-moves" "$emit_moves"
    validate_optional_positive_integer "--max-reads" "$max_reads"

    # Confirm that the selected output paths are under data/.
    if [[ "$output_format" == "fastq" || "$output_format" == "both" ]]; then
        validate_data_output_dir "FASTQ output directory" "$fastq_output_dir" "$pod5_data_root/data"
    fi
    validate_data_output_dir "BAM output directory" "$bam_output_dir" "$pod5_data_root/data"
    if [[ "$output_format" == "both" ]]; then
        validate_data_output_dir "Basecalled BAM output directory" "$basecalled_bam_output_dir" "$pod5_data_root/data"
    fi
    if [[ "$barcode_mode" == "demux" ]]; then
        validate_data_output_dir "Barcode output directory" "$barcode_output_dir" "$pod5_data_root/data"
    fi

    # Confirm that all required directories are available on the compute node.
    require_existing_dir "Workflow root" "$WORKFLOW_ROOT"
    require_existing_dir "Log directory" "$LOG_DIR"
    if [[ "$output_format" == "fastq" || "$output_format" == "both" ]]; then
        require_existing_dir "FASTQ output directory" "$fastq_output_dir"
    fi
    require_existing_dir "BAM output directory" "$bam_output_dir"
    if [[ "$output_format" == "both" ]]; then
        require_existing_dir "Basecalled BAM output directory" "$basecalled_bam_output_dir"
    fi
    if [[ "$barcode_mode" == "demux" ]]; then
        require_existing_dir "Barcode output directory" "$barcode_output_dir"
    fi

    # Prepare the software required for the selected workflow branch.
    load_nanopore_modules "$output_format"

    local pod5_basename
    local pod5_prefix
    local actual_pod5
    local dorado_log
    local alignment_log
    local dorado_output
    local dorado_call_format
    local fastq_file
    local bam_file
    local basecalled_bam_file
    local derived_fastq
    local dorado_output_dir
    local dorado_mm2_opts
    local barcode_dir
    local barcode_fastq_dir
    local barcode_aligned_dir
    local barcode_derived_fastq
    local barcode_file
    local barcode_label
    local barcode_bam
    local barcode_outputs=()
    local barcode_files=()

    # Derive the sample name from the POD5 file or directory name. Uses
    # $pod5 as originally given (before any extraction/normalization
    # below), so this matches the name submit_workflow already used to
    # build the SLURM job name and output paths.
    #
    # Examples:
    # sample.pod5 becomes sample
    # sample.zip  becomes sample
    pod5_basename="$(basename "$pod5")"
    pod5_prefix="${pod5_basename%.pod5}"
    pod5_prefix="${pod5_prefix%.zip}"

    # Construct separate log files for Dorado and alignment/QC.
    #
    # The SLURM job ID is included so repeated runs do not overwrite
    # one another's logs.
    dorado_log="$LOG_DIR/${pod5_prefix}.${SLURM_JOB_ID:-manual}.dorado.log"
    if [[ "$barcode_mode" == "demux" ]]; then
        alignment_log="$LOG_DIR/${pod5_prefix}.${SLURM_JOB_ID:-manual}.barcode_alignment.log"
    elif [[ "$output_format" == "bam" ]]; then
        alignment_log="$LOG_DIR/${pod5_prefix}.${SLURM_JOB_ID:-manual}.dorado_alignment.log"
    elif [[ "$output_format" == "both" ]]; then
        alignment_log="$LOG_DIR/${pod5_prefix}.${SLURM_JOB_ID:-manual}.both_alignment.log"
    else
        alignment_log="$LOG_DIR/${pod5_prefix}.${SLURM_JOB_ID:-manual}.minimap2.log"
    fi

    # Print the main workflow status and log locations.
    echo "[INFO] Nanopore workflow started: $(date)"
    echo "[INFO] POD5 prefix: $pod5_prefix"
    echo "[INFO] Dorado log: $dorado_log"
    echo "[INFO] Alignment/QC log: $alignment_log"
    echo "[INFO] Dorado output format: $output_format"

    # Pick where Dorado itself should write, and which format to ask it
    # for. "both" basecalls once into an unaligned BAM (bam_basecalled);
    # the FASTQ and aligned BAM are produced afterward, below.
    dorado_call_format="$output_format"
    dorado_output_dir="$fastq_output_dir"
    if [[ "$output_format" == "bam" ]]; then
        dorado_output_dir="$bam_output_dir"
    elif [[ "$output_format" == "both" ]]; then
        dorado_output_dir="$basecalled_bam_output_dir"
        dorado_call_format="bam_basecalled"
    fi

    # Dorado exposes supported minimap2 alignment settings through
    # --mm2-opts. This only applies to Dorado's own basecall+align path
    # (plain --output-format bam); "both" aligns separately below, and
    # the workflow's sort, index, MAPQ/read-length, and primary-only
    # settings are handled by samtools after alignment either way.
    dorado_mm2_opts=""
    if [[ "$output_format" == "bam" ]]; then
        dorado_mm2_opts="-x $preset"
        if [[ "$secondary" == "yes" ]]; then
            dorado_mm2_opts="$dorado_mm2_opts --secondary"
        fi
    fi

    # Fully resolve the POD5 input now, on the compute node. This is
    # where any .zip archive actually gets extracted (extract_pod5_zip)
    # and any nested/zipped directory gets normalized into something
    # Dorado can read directly (normalize_pod5_directory) -- deliberately
    # deferred until here so that I/O never runs on the login node while
    # a user is sitting at the interactive submit_workflow prompts.
    if ! actual_pod5="$(resolve_pod5_path "$pod5")"; then
        echo "[ERROR] Could not resolve POD5 input on the compute node: $pod5"
        exit 1
    fi
    echo "[INFO] Resolved POD5 input: $actual_pod5"

    # Run the Dorado basecalling script.
    #
    # The dorado_basecall.sh script prints its generated output path
    # to standard output. Command substitution captures that path in
    # the dorado_output variable.
    dorado_output="$("$SRC_DIR/dorado_basecall.sh" \
        --pod5 "$actual_pod5" \
        --output-dir "$dorado_output_dir" \
        --log-file "$dorado_log" \
        --dorado-model "$dorado_model" \
        --device "$device" \
        --min-qscore "$min_qscore" \
        --max-reads "$max_reads" \
        --output-format "$dorado_call_format" \
        --reference "$reference" \
        --emit-moves "$emit_moves" \
        --mm2-opts "$dorado_mm2_opts" \
        --kit-name "$kit_name" \
        --barcode-mode "$barcode_mode" \
        --barcode-output-dir "$barcode_output_dir" \
        --modified-bases "$modified_bases")"

    echo "[INFO] Dorado complete: $dorado_output"

    if [[ "$barcode_mode" == "demux" && "$output_format" == "both" ]]; then
        # Demultiplexing a "both" run: Dorado already basecalled once into
        # a barcode-tagged unaligned BAM and demuxed it by barcode
        # (dorado_basecall.sh, above). Each barcode's unaligned BAM is now
        # turned into a derived FASTQ (samtools, no re-basecalling) and an
        # aligned BAM, exactly like the non-demux "both" path -- this is
        # what gives you both a demuxed FASTQ and a demuxed aligned BAM
        # per barcode automatically.
        barcode_dir="$barcode_output_dir/${pod5_prefix}_${SLURM_JOB_ID:-manual}"
        barcode_fastq_dir="$fastq_output_dir/${pod5_prefix}_${SLURM_JOB_ID:-manual}"
        barcode_aligned_dir="$bam_output_dir/${pod5_prefix}_${SLURM_JOB_ID:-manual}"
        echo "[INFO] Processing demuxed basecalled BAMs in: $barcode_dir"
        mkdir -p "$barcode_fastq_dir" "$barcode_aligned_dir"

        for barcode_file in "$barcode_dir"/*.bam; do
            [[ -f "$barcode_file" ]] || continue
            barcode_files+=("$barcode_file")
        done

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
            echo "Output format: both"
            echo "Reference:     $reference"
            echo "========================================="
        } > "$alignment_log"

        for barcode_file in "${barcode_files[@]}"; do
            barcode_label="$(basename "$barcode_file")"
            barcode_label="${barcode_label%.bam}"

            echo "[INFO] Processing barcode file: $barcode_file"

            # Sort and index this barcode's basecalled (unaligned) BAM
            # in place before deriving a FASTQ from it or aligning it.
            echo "[INFO] Sorting and indexing basecalled BAM for $barcode_label: $barcode_file" >> "$alignment_log"
            samtools sort -@ "$threads" -o "${barcode_file}.sorted.tmp" "$barcode_file" >> "$alignment_log" 2>&1
            mv "${barcode_file}.sorted.tmp" "$barcode_file"
            samtools index "$barcode_file" >> "$alignment_log" 2>&1

            barcode_derived_fastq="$barcode_fastq_dir/${barcode_label}.fastq"
            echo "[INFO] Deriving FASTQ for $barcode_label: $barcode_derived_fastq" >> "$alignment_log"
            samtools fastq "$barcode_file" > "$barcode_derived_fastq" 2>> "$alignment_log"

            if [[ ! -s "$barcode_derived_fastq" ]]; then
                echo "[ERROR] Derived FASTQ is empty or was not created: $barcode_derived_fastq" | tee -a "$alignment_log"
                exit 1
            fi

            # Align the sorted, indexed basecalled BAM directly via
            # 'dorado aligner' (--basecalled-bam) rather than the derived
            # FASTQ, so move tables / modification tags survive alignment.
            barcode_bam="$("$SRC_DIR/minimap2_alignment.sh" \
                --basecalled-bam "$barcode_file" \
                --reference "$reference" \
                --output-dir "$barcode_aligned_dir" \
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

            barcode_outputs+=("$barcode_bam")
            echo "[INFO] Barcode alignment/QC complete: $barcode_bam"
        done

        {
            echo ""
            echo "========================================="
            echo "  BARCODE ALIGNMENT/QC COMPLETE"
            echo "========================================="
            echo "Completed: $(date)"
            echo "Basecalled (unaligned) BAMs: $barcode_dir"
            echo "Derived FASTQs:              $barcode_fastq_dir"
            echo "Final aligned BAMs:"
            printf '%s\n' "${barcode_outputs[@]}"
            echo "========================================="
        } >> "$alignment_log"
    elif [[ "$barcode_mode" == "demux" ]]; then
        # Legacy single-format demux (fastq-only or bam-only), kept for
        # direct invocation of this script outside the interactive
        # submit_workflow path, which now always uses --output-format both
        # for demux so every barcode gets both a FASTQ and an aligned BAM.
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
    elif [[ "$output_format" == "both" ]]; then
        # "both": Dorado already basecalled once into an unaligned BAM
        # (dorado_output, above). Sort and index that basecalled BAM in
        # place, derive a FASTQ from it with samtools (no re-basecalling,
        # kept purely as a FASTQ deliverable), then align the sorted
        # basecalled BAM itself via 'dorado aligner' (not the FASTQ) so
        # its tags -- move tables, base-modification calls -- survive
        # alignment instead of being discarded by a FASTQ round-trip.
        basecalled_bam_file="$dorado_output"
        derived_fastq="$fastq_output_dir/${pod5_prefix}_${SLURM_JOB_ID:-manual}.fastq"

        echo "[INFO] Basecalled (unaligned) BAM: $basecalled_bam_file"
        echo "[INFO] Sorting and indexing basecalled BAM: $basecalled_bam_file" >> "$alignment_log"

        samtools sort -@ "$threads" -o "${basecalled_bam_file}.sorted.tmp" "$basecalled_bam_file" >> "$alignment_log" 2>&1
        mv "${basecalled_bam_file}.sorted.tmp" "$basecalled_bam_file"
        samtools index "$basecalled_bam_file" >> "$alignment_log" 2>&1

        echo "[INFO] Deriving FASTQ from basecalled BAM: $derived_fastq" >> "$alignment_log"
        samtools fastq "$basecalled_bam_file" > "$derived_fastq" 2>> "$alignment_log"

        if [[ ! -s "$derived_fastq" ]]; then
            echo "[ERROR] Derived FASTQ is empty or was not created: $derived_fastq" | tee -a "$alignment_log"
            exit 1
        fi

        echo "[INFO] Derived FASTQ: $derived_fastq"

        bam_file="$("$SRC_DIR/minimap2_alignment.sh" \
            --basecalled-bam "$basecalled_bam_file" \
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

        echo "[INFO] Dorado aligner complete: $bam_file"
        echo "[INFO] Basecalled BAM (sorted, indexed): $basecalled_bam_file"
        echo "[INFO] FASTQ:                            $derived_fastq"
        echo "[INFO] Aligned BAM:                      $bam_file"
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
    local basecalled_bam_output_dir=""
    local basecalled_bam_output_input=""
    local output_format
    local emit_moves
    local run_demux
    local max_reads

    # Variables used to construct sample-specific names.
    local pod5_basename
    local pod5_prefix
    local pod5_data_root
    local default_fastq_output_dir
    local default_basecalled_bam_dir
    local default_aligned_bam_dir

    # Create the workflow's standard directories on the submission node.
    # (FASTQ/BAM output directories, and the log directory, are created
    # later, once the POD5 input is resolved and the pod5-relative
    # defaults can be computed.)
    mkdir -p "$POD5_DIR" "$REFERENCE_SEARCH_DIR"

    # When no POD5 was provided as a positional argument, display the
    # available POD5 files/directories/archives and prompt the user to
    # select one.
    while true; do
        if [[ -z "$pod5_input" ]]; then
            echo "[INFO] Available POD5 files, directories, and .zip archives in default path $POD5_DIR:"
            find "$POD5_DIR" -maxdepth 1 -type f \( -name "*.pod5" -o -name "*.zip" \) -printf "  %f\n" | sort
            find "$POD5_DIR" -mindepth 1 -maxdepth 1 -type d -printf "  %f/\n" | sort
            echo
            read -r -p "Enter a POD5 file/directory/.zip from data/pod5, or an explicit path plus file/directory: " pod5_input
        fi

        if [[ -z "$pod5_input" ]]; then
            warn_unexpected_input "$pod5_input" "a .pod5 file, directory, or .zip archive in $POD5_DIR, or an explicit path like /path/to/sample.pod5"
            continue
        fi

        if pod5="$(resolve_pod5_input_path "$pod5_input")" && pod5_input_looks_valid "$pod5"; then
            break
        fi

        warn_unexpected_input "$pod5_input" "a valid .pod5 file/directory/.zip in $POD5_DIR, or an explicit path like /path/to/sample.pod5"
        pod5_input=""
    done

    # Default FASTQ/BAM output locations live one directory above wherever
    # the resolved POD5 input actually is, e.g.:
    #   /projects/run1/pod5/  ->  /projects/run1/data/fastq
    #                             /projects/run1/data/bam_aligned
    #                             /projects/run1/data/bam_basecalled
    pod5_data_root="$(derive_pod5_data_root "$pod5")"
    default_fastq_output_dir="$pod5_data_root/data/fastq"
    default_basecalled_bam_dir="$pod5_data_root/data/bam_basecalled"
    default_aligned_bam_dir="$pod5_data_root/data/bam_aligned"

    # Logs live one directory above the POD5 input too, alongside the
    # data/ output directories, rather than in a fixed repo-relative spot.
    LOG_DIR="$pod5_data_root/logs/nanopore_sequence_workflow"
    mkdir -p "$LOG_DIR"

    # Prompt for the reference FASTA when it was not supplied.
    while true; do
        if [[ -z "$reference_input" ]]; then
            echo "[INFO] Available reference FASTA files in default path $REFERENCE_SEARCH_DIR:"
            find "$REFERENCE_SEARCH_DIR" \
                -maxdepth 1 \
                -type f \
                \( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" \) \
                -printf "  %f\n" |
                sort
            echo
            read -r -p "Enter a reference FASTA filename from data/fastq, or an explicit path plus filename: " reference_input
        fi

        if [[ -z "$reference_input" ]]; then
            warn_unexpected_input "$reference_input" "a FASTA filename in $REFERENCE_SEARCH_DIR, for example reference.fna, or an explicit path like /path/to/reference.fna"
            continue
        fi

        if reference="$(resolve_reference_file "$reference_input")"; then
            break
        fi

        warn_unexpected_input "$reference_input" "a valid FASTA filename in $REFERENCE_SEARCH_DIR, or an explicit path like /path/to/reference.fna"
        reference_input=""
    done

    echo
    echo "Dorado demultiplexing"

    # Ask about demux before the Dorado output-format question below,
    # because demultiplexing always produces both a demuxed FASTQ and a
    # demuxed (aligned) BAM per barcode -- the format question is skipped
    # and forced to "both" whenever demux is on.
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
        list_available_barcode_kits
        while [[ -z "$kit_name" ]]; do
            kit_name="$(prompt_optional "Enter barcode-kit / --kit-name")"
            if [[ -z "$kit_name" ]]; then
                warn_unexpected_input "$kit_name" "a barcode kit name, for example SQK-NBD114-24"
            fi
        done
    fi

    echo
    echo "Dorado output selection"
    if [[ "$barcode_mode" == "demux" ]]; then
        output_format="both"
        echo "[INFO] Demultiplexing always produces both a demuxed FASTQ and a demuxed (aligned) BAM per barcode."
        echo "[INFO] Dorado output format is set to 'both' automatically."
    else
        while true; do
            output_format="$(prompt_with_default "Should Dorado produce fastq, bam, or both?" "bam")"
            case "$output_format" in
                fastq|bam|both) break ;;
                *) warn_unexpected_input "$output_format" "fastq, bam, or both, for example bam" ;;
            esac
        done
    fi

    # Prompt for the FASTQ output directory when FASTQ will actually be
    # produced: the fastq-only path, or "both" (which derives a FASTQ
    # from the basecalled BAM -- always the case when demuxing).
    if [[ ( "$output_format" == "fastq" || "$output_format" == "both" ) && -z "$fastq_output_input" ]]; then
        fastq_output_input="$(prompt_with_default "Enter FASTQ output directory under data/" "$default_fastq_output_dir")"
    elif [[ -z "$fastq_output_input" ]]; then
        fastq_output_input="$default_fastq_output_dir"
    fi

    # Prompt for the aligned BAM output directory when it was not supplied.
    if [[ -z "$bam_output_input" ]]; then
        bam_output_input="$(prompt_with_default "Enter aligned BAM output directory under data/" "$default_aligned_bam_dir")"
    fi

    # Prompt for a separate basecalled (unaligned) BAM output directory
    # only when producing both (always the case when demuxing).
    if [[ "$output_format" == "both" ]]; then
        basecalled_bam_output_input="$(prompt_with_default "Enter basecalled (unaligned) BAM output directory under data/" "$default_basecalled_bam_dir")"
    fi

    # Convert the supplied output-directory values into complete paths.
    while true; do
        fastq_output_dir="$(resolve_data_output_dir "$fastq_output_input" "$pod5_data_root")"
        if [[ ( "$output_format" != "fastq" && "$output_format" != "both" ) || "$fastq_output_dir" == "$pod5_data_root/data"/* ]]; then
            break
        fi
        warn_unexpected_input "$fastq_output_input" "a FASTQ output directory under $pod5_data_root/data"
        fastq_output_input="$(prompt_with_default "Enter FASTQ output directory under data/" "$default_fastq_output_dir")"
    done

    while true; do
        bam_output_dir="$(resolve_data_output_dir "$bam_output_input" "$pod5_data_root")"
        if [[ "$bam_output_dir" == "$pod5_data_root/data"/* ]]; then
            break
        fi
        warn_unexpected_input "$bam_output_input" "a BAM output directory under $pod5_data_root/data"
        bam_output_input="$(prompt_with_default "Enter aligned BAM output directory under data/" "$default_aligned_bam_dir")"
    done

    if [[ "$output_format" == "both" ]]; then
        while true; do
            basecalled_bam_output_dir="$(resolve_data_output_dir "$basecalled_bam_output_input" "$pod5_data_root")"
            if [[ "$basecalled_bam_output_dir" == "$pod5_data_root/data"/* ]]; then
                break
            fi
            warn_unexpected_input "$basecalled_bam_output_input" "a basecalled BAM output directory under $pod5_data_root/data"
            basecalled_bam_output_input="$(prompt_with_default "Enter basecalled (unaligned) BAM output directory under data/" "$default_basecalled_bam_dir")"
        done
    fi

    # Derive the sample prefix from the POD5 file or directory name.
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

    case "$output_format" in
        fastq)
            echo "[INFO] --emit-fastq is enabled because the FASTQ path uses Minimap2 alignment."
            echo "[INFO] --emit-moves is ignored for FASTQ output because FASTQ cannot store move-table tags."
            ;;
        bam)
            echo "[INFO] Dorado will basecall and align with --reference, producing an aligned BAM."
            ;;
        both)
            echo "[INFO] Dorado will basecall once into an unaligned BAM (basecalled BAM directory)."
            echo "[INFO] A FASTQ is then derived from that BAM (no re-basecalling) and aligned separately to produce the aligned BAM."
            ;;
    esac

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
    if [[ "$output_format" == "fastq" || "$output_format" == "both" ]]; then
        mkdir -p "$fastq_output_dir"
    fi
    mkdir -p "$bam_output_dir"
    if [[ "$output_format" == "both" ]]; then
        mkdir -p "$basecalled_bam_output_dir"
    fi

    # Load the required software on the submission node.
    #
    # This checks that the software environment is available before
    # the workflow is submitted.
    load_nanopore_modules "$output_format"

    # Display the resolved workflow and output locations.
    echo
    echo "[INFO] Workflow root: $WORKFLOW_ROOT"
    echo "[INFO] POD5 data root (used for default output paths): $pod5_data_root"
    echo "[INFO] Dorado output format: $output_format"
    if [[ "$output_format" == "fastq" || "$output_format" == "both" ]]; then
        echo "[INFO] FASTQ output directory: $fastq_output_dir"
    fi
    echo "[INFO] Aligned BAM output directory: $bam_output_dir"
    if [[ "$output_format" == "both" ]]; then
        echo "[INFO] Basecalled BAM output directory: $basecalled_bam_output_dir"
    fi
    echo
    echo "Submitting nanopore workflow to SLURM..."

    # Demux against a "both" run splits the basecalled (unaligned) BAM, so
    # barcode output belongs alongside it; otherwise it belongs alongside
    # whatever aligned BAM directory you actually chose (not a hardcoded
    # default -- this is the fix for the earlier "always uses the default
    # aligned-BAM directory" issue).
    local barcode_output_dir_to_pass="$bam_output_dir"
    if [[ "$output_format" == "both" ]]; then
        barcode_output_dir_to_pass="$basecalled_bam_output_dir"
    fi

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
        --basecalled-bam-output-dir "$basecalled_bam_output_dir" \
        --dorado-model "$dorado_model" \
        --device "$device" \
        --min-qscore "$min_qscore" \
        --max-reads "$max_reads" \
        --output-format "$output_format" \
        --emit-moves "$emit_moves" \
        --kit-name "$kit_name" \
        --barcode-mode "$barcode_mode" \
        --barcode-output-dir "$barcode_output_dir_to_pass" \
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
    case "$output_format" in
        fastq)
            echo "[INFO] Expected FASTQ: $fastq_output_dir/${pod5_prefix}_${log_job_id}.fastq"
            echo "[INFO] Expected raw BAM: $bam_output_dir/${pod5_prefix}_${log_job_id}.bam"
            ;;
        bam)
            echo "[INFO] Expected Dorado aligned BAM: $bam_output_dir/${pod5_prefix}_${log_job_id}.bam"
            echo "[INFO] Expected raw BAM: $bam_output_dir/${pod5_prefix}_${log_job_id}.dorado.raw.bam"
            ;;
        both)
            echo "[INFO] Expected basecalled (unaligned) BAM: $basecalled_bam_output_dir/${pod5_prefix}_${log_job_id}.bam"
            echo "[INFO] Expected derived FASTQ: $fastq_output_dir/${pod5_prefix}_${log_job_id}.fastq"
            echo "[INFO] Expected raw aligned BAM: $bam_output_dir/${pod5_prefix}_${log_job_id}.bam"
            ;;
    esac
    if [[ "$barcode_mode" == "demux" ]]; then
        echo "[INFO] Expected demuxed basecalled BAMs: $barcode_output_dir_to_pass/${pod5_prefix}_${log_job_id}"
        if [[ "$output_format" == "both" ]]; then
            echo "[INFO] Expected demuxed FASTQs:         $fastq_output_dir/${pod5_prefix}_${log_job_id}"
            echo "[INFO] Expected demuxed aligned BAMs:   $bam_output_dir/${pod5_prefix}_${log_job_id}"
        fi
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
    elif [[ "$output_format" == "both" ]]; then
        echo "[INFO] Alignment/QC log: $LOG_DIR/${pod5_prefix}.${log_job_id}.both_alignment.log"
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
