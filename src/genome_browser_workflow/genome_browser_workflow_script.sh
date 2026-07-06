#!/bin/bash
#SBATCH --job-name=genome_browser_brdu
#SBATCH --output=genome_browser_brdu_%j.out
#SBATCH --error=genome_browser_brdu_%j.err
#SBATCH --cpus-per-task=12
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --partition=normal

set -euo pipefail

if command -v module >/dev/null 2>&1; then
    module load python/3.13.7 || module load python
    module load samtools
    module load modkit
else
    echo "[WARN] Environment modules are not available in this shell."
    echo "[WARN] Continuing with samtools, modkit, and python from PATH."
fi

WORKFLOW_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"

BAM_DIR="$WORKFLOW_ROOT/data/bam"
BEDGRAPH_DIR="$WORKFLOW_ROOT/data/bedgraph"
SORTED_BAM_DIR="$WORKFLOW_ROOT/data/sorted_bam"
INDEX_DIR="$WORKFLOW_ROOT/data/index_sorted_bam_bai"
SRC_DIR="$WORKFLOW_ROOT/src/genome_browser_workflow"
RESULTS_DIR="$WORKFLOW_ROOT/results/genome_browser_results"
DEFAULT_REF="$WORKFLOW_ROOT/data/ncbi/W303/ncbi_dataset/GCA_002163515.1_ASM216351v1_genomic.fna"
DEFAULT_G4_BED="$WORKFLOW_ROOT/data/bed/W303_g4_motifs.bed"
DEFAULT_TE_BED="$WORKFLOW_ROOT/data/bed/w303_te_and_ltrs.bed"
DEFAULT_TRNA_BED="$WORKFLOW_ROOT/data/bed/trna_coordinates.bed"

mkdir -p \
    "$BAM_DIR" \
    "$BEDGRAPH_DIR" \
    "$SORTED_BAM_DIR" \
    "$INDEX_DIR" \
    "$RESULTS_DIR"

usage() {
    echo "Usage: bash src/genome_browser_workflow/genome_browser_workflow_script.sh [BAM] [output_prefix]"
    echo
    echo "BAM may be a filename in data/bam or an explicit path."
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

BAM_INPUT="${1:-}"
OUTPUT_PREFIX="${2:-}"

if [[ -z "$BAM_INPUT" ]]; then
    echo "[INFO] Available BAM files in $BAM_DIR:"
    find "$BAM_DIR" -maxdepth 1 -type f -name "*.bam" -printf "  %f\n" | sort
    echo
    read -r -p "Enter the BAM filename from data/bam: " BAM_INPUT
fi

if [[ -z "$BAM_INPUT" ]]; then
    echo "[ERROR] No BAM file was provided."
    exit 1
fi

if [[ -f "$BAM_INPUT" ]]; then
    BAM_PATH="$(cd "$(dirname "$BAM_INPUT")" && pwd)/$(basename "$BAM_INPUT")"
else
    BAM_PATH="$BAM_DIR/$BAM_INPUT"
fi

if [[ ! -f "$BAM_PATH" ]]; then
    echo "[ERROR] BAM file not found: $BAM_PATH"
    exit 1
fi

if [[ -z "$OUTPUT_PREFIX" ]]; then
    BAM_NAME="$(basename "$BAM_PATH")"
    OUTPUT_PREFIX="${BAM_NAME%.bam}"
fi

read -r -p "Reference FASTA for modkit pileup [$DEFAULT_REF]: " REF_INPUT
REF_PATH="${REF_INPUT:-$DEFAULT_REF}"

if [[ ! -f "$REF_PATH" ]]; then
    echo "[ERROR] Reference FASTA not found: $REF_PATH"
    exit 1
fi

echo
echo "Choose genome browser output mode:"
echo "  1) smoothed"
echo "  2) unsmoothed"
read -r -p "Enter smoothed or unsmoothed [smoothed]: " PLOT_MODE
PLOT_MODE="${PLOT_MODE:-smoothed}"

case "${PLOT_MODE,,}" in
    smoothed|smooth|s|1)
        PLOT_MODE="smoothed"
        GENERATION_SCRIPT="$SRC_DIR/genomic_browser_generation.py"
        MODE_RESULTS_DIR="$RESULTS_DIR/smoothed"
        ;;
    unsmoothed|unsmooth|raw|u|2)
        PLOT_MODE="unsmoothed"
        GENERATION_SCRIPT="$SRC_DIR/genomic_browser_generation_unsmoothed.py"
        MODE_RESULTS_DIR="$RESULTS_DIR/unsmoothed"
        ;;
    *)
        echo "[ERROR] Invalid mode: $PLOT_MODE"
        echo "[ERROR] Expected smoothed or unsmoothed."
        exit 1
        ;;
esac

VENV_DIR="$SRC_DIR/.genome_browser_env"

if [[ ! -x "$VENV_DIR/bin/python" ]]; then
    echo "[INFO] Creating virtual environment: $VENV_DIR"
    python3 -m venv --clear "$VENV_DIR"
fi

echo "[INFO] Activating virtual environment."
source "$VENV_DIR/bin/activate"

echo "[INFO] Installing Python requirements."
python -m pip install -r "$SRC_DIR/requirements.txt" --quiet

echo "[INFO] Running raw BrdU bedgraph extraction."
python "$SRC_DIR/raw_data_extraction_on_bam.py" \
    "$BAM_PATH" \
    --ref "$REF_PATH" \
    --output-prefix "$OUTPUT_PREFIX" \
    --threads "${SLURM_CPUS_PER_TASK:-12}" \
    --bedgraph-dir "$BEDGRAPH_DIR"

POSITIVE_OUTPUT="$BEDGRAPH_DIR/$OUTPUT_PREFIX.positive.bedgraph"
NEGATIVE_OUTPUT="$BEDGRAPH_DIR/$OUTPUT_PREFIX.negative.bedgraph"

if [[ ! -s "$POSITIVE_OUTPUT" ]]; then
    echo "[ERROR] Positive bedgraph file is empty or missing: $POSITIVE_OUTPUT"
    deactivate
    exit 1
fi

if [[ ! -s "$NEGATIVE_OUTPUT" ]]; then
    echo "[ERROR] Negative bedgraph file is empty or missing: $NEGATIVE_OUTPUT"
    deactivate
    exit 1
fi

echo "[INFO] Positive strand bedgraph: $POSITIVE_OUTPUT"
echo "[INFO] Negative strand bedgraph: $NEGATIVE_OUTPUT"

mkdir -p "$MODE_RESULTS_DIR"

GENERATION_ARGS=(
    --positive-bedgraph "$POSITIVE_OUTPUT"
    --negative-bedgraph "$NEGATIVE_OUTPUT"
    --output-dir "$MODE_RESULTS_DIR"
    --prefix "$OUTPUT_PREFIX"
)

if [[ -f "$DEFAULT_G4_BED" ]]; then
    GENERATION_ARGS+=(--g4-bed "$DEFAULT_G4_BED")
fi

if [[ -f "$DEFAULT_TE_BED" ]]; then
    GENERATION_ARGS+=(--te-bed "$DEFAULT_TE_BED")
fi

if [[ -f "$DEFAULT_TRNA_BED" ]]; then
    GENERATION_ARGS+=(--trna-bed "$DEFAULT_TRNA_BED")
fi

echo "[INFO] Generating $PLOT_MODE genome browser plots."
python "$GENERATION_SCRIPT" "${GENERATION_ARGS[@]}"

deactivate

echo "[INFO] Genome browser workflow complete."
echo "[INFO] Results directory: $MODE_RESULTS_DIR"
