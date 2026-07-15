#!/bin/bash
#SBATCH --job-name=genome_browser_brdu
#SBATCH --cpus-per-task=12
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --partition=normal

set -euo pipefail

WORKFLOW_ROOT_ARG=""
for ((arg_i = 1; arg_i <= $#; arg_i++)); do
    if [[ "${!arg_i}" == "--workflow-root" ]]; then
        next_arg_i=$((arg_i + 1))
        WORKFLOW_ROOT_ARG="${!next_arg_i:-}"
        break
    fi
done

if [[ -n "$WORKFLOW_ROOT_ARG" ]]; then
    WORKFLOW_ROOT="$(cd "$WORKFLOW_ROOT_ARG" && pwd)"
elif [[ -n "${GENOME_BROWSER_WORKFLOW_ROOT:-}" ]]; then
    WORKFLOW_ROOT="$GENOME_BROWSER_WORKFLOW_ROOT"
else
    WORKFLOW_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
fi

BAM_DIR="$WORKFLOW_ROOT/data/bam"
BEDGRAPH_DIR="$WORKFLOW_ROOT/data/bedgraph"
SORTED_BAM_DIR="$WORKFLOW_ROOT/data/sorted_bam"
SRC_DIR="$WORKFLOW_ROOT/src/genome_browser_workflow"
SCRIPT_PATH="$SRC_DIR/genome_browser_workflow_script.sh"
RESULTS_DIR="$WORKFLOW_ROOT/results/genome_browser_results"
LOG_DIR="$WORKFLOW_ROOT/logs/genome_browser_workflow"
DEFAULT_REF="$WORKFLOW_ROOT/data/ncbi/W303/ncbi_dataset/GCA_002163515.1_ASM216351v1_genomic.fna"
DEFAULT_G4_BED="$WORKFLOW_ROOT/data/bed/W303_g4_motifs.bed"
DEFAULT_TE_BED="$WORKFLOW_ROOT/data/bed/w303_te_and_ltrs.bed"
DEFAULT_TRNA_BED="$WORKFLOW_ROOT/data/bed/trna_coordinates.bed"

usage() {
    echo "Usage: bash src/genome_browser_workflow/genome_browser_workflow_script.sh [BAM] [output_prefix]"
    echo
    echo "BAM is resolved by filename under data/bam."
    echo "If arguments are omitted, the workflow prompts for inputs before submitting a SLURM job."
}

absolute_existing_file() {
    local path="$1"
    local dir
    local file

    dir="$(cd "$(dirname "$path")" && pwd)"
    file="$(basename "$path")"
    printf '%s/%s\n' "$dir" "$file"
}

resolve_bam_file() {
    local bam_input="$1"
    local bam_filename
    local bam_path

    bam_filename="$(basename "$bam_input")"
    bam_path="$BAM_DIR/$bam_filename"

    if [[ -f "$bam_path" ]]; then
        absolute_existing_file "$bam_path"
        return 0
    fi

    return 1
}

resolve_reference_file() {
    local reference_input="$1"
    local reference_path

    if [[ -f "$reference_input" ]]; then
        absolute_existing_file "$reference_input"
        return 0
    fi

    reference_path="$WORKFLOW_ROOT/$reference_input"
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

require_existing_dir() {
    local name="$1"
    local path="$2"

    if [[ ! -d "$path" ]]; then
        echo "[ERROR] $name does not exist: $path"
        echo "[ERROR] Confirm the project path is mounted on the compute node."
        exit 1
    fi
}

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

run_job() {
    local bam_path=""
    local reference=""
    local output_prefix=""
    local plot_mode=""
    local phase_label=""
    local generation_script
    local mode_results_dir
    local venv_dir
    local positive_output
    local negative_output
    local workflow_log

    shift

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --workflow-root) shift 2 ;;
            --bam) bam_path="$2"; shift 2 ;;
            --ref) reference="$2"; shift 2 ;;
            --output-prefix) output_prefix="$2"; shift 2 ;;
            --plot-mode) plot_mode="$2"; shift 2 ;;
            --phase-label) phase_label="$2"; shift 2 ;;
            *)
                echo "[ERROR] Unknown job argument: $1"
                exit 1
                ;;
        esac
    done

    mkdir -p "$LOG_DIR"
    workflow_log="$LOG_DIR/${output_prefix}.${SLURM_JOB_ID:-manual}.workflow.log"
    exec > >(tee -a "$workflow_log") 2>&1

    echo "[INFO] Genome browser workflow started: $(date)"
    echo "[INFO] Workflow log: $workflow_log"

    require_existing_dir "Workflow root" "$WORKFLOW_ROOT"
    require_existing_dir "Source directory" "$SRC_DIR"
    require_existing_dir "BAM directory" "$BAM_DIR"

    mkdir -p "$BEDGRAPH_DIR" "$SORTED_BAM_DIR" "$RESULTS_DIR" "$LOG_DIR"

    if [[ ! -f "$bam_path" ]]; then
        echo "[ERROR] BAM file not found on compute node: $bam_path"
        exit 1
    fi

    if [[ ! -f "$reference" ]]; then
        echo "[ERROR] Reference FASTA not found on compute node: $reference"
        exit 1
    fi

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

    load_genome_browser_modules

    venv_dir="$SRC_DIR/.genome_browser_env"
    if [[ ! -x "$venv_dir/bin/python" ]]; then
        echo "[INFO] Creating virtual environment: $venv_dir"
        python3 -m venv --clear "$venv_dir"
    fi

    echo "[INFO] Activating virtual environment."
    source "$venv_dir/bin/activate"

    echo "[INFO] Installing Python requirements."
    python -m pip install -r "$SRC_DIR/requirements.txt" --quiet

    positive_output="$BEDGRAPH_DIR/$output_prefix.positive.bedgraph"
    negative_output="$BEDGRAPH_DIR/$output_prefix.negative.bedgraph"

    echo "[INFO] Running raw BrdU bedgraph extraction."
    python "$SRC_DIR/raw_data_extraction_on_bam.py" \
        "$bam_path" \
        --ref "$reference" \
        --output-prefix "$output_prefix" \
        --threads "${SLURM_CPUS_PER_TASK:-12}" \
        --bedgraph-dir "$BEDGRAPH_DIR"

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

    generation_args=(
        --positive-bedgraph "$positive_output"
        --negative-bedgraph "$negative_output"
        --output-dir "$mode_results_dir"
        --prefix "$output_prefix"
        --phase-label "$phase_label"
    )

    if [[ -f "$DEFAULT_G4_BED" ]]; then
        generation_args+=(--g4-bed "$DEFAULT_G4_BED")
    fi

    if [[ -f "$DEFAULT_TE_BED" ]]; then
        generation_args+=(--te-bed "$DEFAULT_TE_BED")
    fi

    if [[ -f "$DEFAULT_TRNA_BED" ]]; then
        generation_args+=(--trna-bed "$DEFAULT_TRNA_BED")
    fi

    echo "[INFO] Generating $plot_mode genome browser plots."
    python "$generation_script" "${generation_args[@]}"

    deactivate

    echo "[INFO] Genome browser workflow complete."
    echo "[INFO] Results directory: $mode_results_dir"
    echo "[INFO] Genome browser workflow finished: $(date)"
}

submit_workflow() {
    local bam_input="${1:-}"
    local output_prefix="${2:-}"
    local ref_input=""
    local bam_path=""
    local reference=""
    local bam_name
    local phase_input
    local phase_label
    local plot_mode_input
    local plot_mode
    local mode_results_dir
    local job_id
    local log_job_id

    mkdir -p "$BAM_DIR" "$BEDGRAPH_DIR" "$SORTED_BAM_DIR" "$RESULTS_DIR" "$LOG_DIR"

    if [[ -z "$bam_input" ]]; then
        echo "[INFO] Available BAM files in $BAM_DIR:"
        find "$BAM_DIR" -maxdepth 1 -type f -name "*.bam" -printf "  %f\n" | sort
        echo
        read -r -p "Enter the BAM filename from data/bam: " bam_input
    fi

    if [[ -z "$bam_input" ]]; then
        echo "[ERROR] No BAM file was provided."
        exit 1
    fi

    if ! bam_path="$(resolve_bam_file "$bam_input")"; then
        echo "[ERROR] BAM file not found under $BAM_DIR: $(basename "$bam_input")"
        exit 1
    fi

    if [[ -z "$output_prefix" ]]; then
        bam_name="$(basename "$bam_path")"
        output_prefix="${bam_name%.bam}"
    fi

    read -r -p "Reference FASTA for modkit pileup [$DEFAULT_REF]: " ref_input
    ref_input="${ref_input:-$DEFAULT_REF}"

    if ! reference="$(resolve_reference_file "$ref_input")"; then
        echo "[ERROR] Reference FASTA not found: $ref_input"
        exit 1
    fi

    echo
    echo "Choose cell-cycle phase for genome browser plot titles:"
    echo "  M) Mitosis"
    echo "  S) S Phase"
    read -r -p "Enter M or S: " phase_input

    if ! phase_label="$(normalize_phase_label "$phase_input")"; then
        echo "[ERROR] Invalid phase: $phase_input"
        echo "[ERROR] Expected M or S."
        exit 1
    fi

    echo
    echo "Choose genome browser output mode:"
    echo "  1) smoothed"
    echo "  2) unsmoothed"
    read -r -p "Enter smoothed or unsmoothed [smoothed]: " plot_mode_input
    plot_mode_input="${plot_mode_input:-smoothed}"

    if ! plot_mode="$(normalize_plot_mode "$plot_mode_input")"; then
        echo "[ERROR] Invalid mode: $plot_mode_input"
        echo "[ERROR] Expected smoothed or unsmoothed."
        exit 1
    fi

    mode_results_dir="$RESULTS_DIR/$plot_mode"
    mkdir -p "$mode_results_dir"

    load_genome_browser_modules

    echo
    echo "[INFO] Workflow root: $WORKFLOW_ROOT"
    echo "[INFO] BAM: $bam_path"
    echo "[INFO] Reference FASTA: $reference"
    echo "[INFO] Output prefix: $output_prefix"
    echo "[INFO] Phase label: $phase_label"
    echo "[INFO] Plot mode: $plot_mode"
    echo
    echo "Submitting genome browser workflow to SLURM..."

    job_id="$(sbatch --parsable \
        --job-name="genome_browser_${output_prefix}" \
        --chdir="$WORKFLOW_ROOT" \
        --output="$LOG_DIR/${output_prefix}.%j.slurm.log" \
        --error="$LOG_DIR/${output_prefix}.%j.slurm.err" \
        --export=ALL,GENOME_BROWSER_WORKFLOW_ROOT="$WORKFLOW_ROOT" \
        "$SCRIPT_PATH" \
        --run-job \
        --workflow-root "$WORKFLOW_ROOT" \
        --bam "$bam_path" \
        --ref "$reference" \
        --output-prefix "$output_prefix" \
        --plot-mode "$plot_mode" \
        --phase-label "$phase_label")"

    log_job_id="${job_id%%;*}"

    echo "[INFO] Submitted SLURM job: $job_id"
    echo "[INFO] Expected positive bedgraph: $BEDGRAPH_DIR/$output_prefix.positive.bedgraph"
    echo "[INFO] Expected negative bedgraph: $BEDGRAPH_DIR/$output_prefix.negative.bedgraph"
    echo "[INFO] Expected results dir:       $mode_results_dir"
    echo "[INFO] SLURM log:                  $LOG_DIR/${output_prefix}.${log_job_id}.slurm.log"
    echo "[INFO] SLURM err:                  $LOG_DIR/${output_prefix}.${log_job_id}.slurm.err"
    echo "[INFO] Workflow log:               $LOG_DIR/${output_prefix}.${log_job_id}.workflow.log"
    echo "[INFO] Modkit log:                 $BEDGRAPH_DIR/$output_prefix.modkit.log"
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

if [[ "${1:-}" == "--run-job" ]]; then
    run_job "$@"
else
    submit_workflow "$@"
fi
