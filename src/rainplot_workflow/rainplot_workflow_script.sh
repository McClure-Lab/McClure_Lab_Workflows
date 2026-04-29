#!/bin/bash
#SBATCH --job-name=extract_brdu
#SBATCH --output=extract_brdu_%j.out
#SBATCH --error=extract_brdu_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --partition=normal

module load python/3.13.7
module load samtools
module load R/4.5.2

# --- Directory structure ---
# workflows/
# ├── data/
# │   ├── bam/
# │   ├── bed/
# │   ├── sorted_bam/
# │   ├── index_sorted_bam_bai/
# │   ├── liftover_chains/
# │   └── ncbi/
# ├── src/
# │   ├── rainplot_workflow/
# │   └── utils/
# ├── tools/
# └── results/
#     └── rainplot_results/

WORKFLOW_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"

BAM_DIR="$WORKFLOW_ROOT/data/bam"
BED_DIR="$WORKFLOW_ROOT/data/bed"
SORTED_BAM_DIR="$WORKFLOW_ROOT/data/sorted_bam"
INDEX_DIR="$WORKFLOW_ROOT/data/index_sorted_bam_bai"
CHAIN_DIR="$WORKFLOW_ROOT/data/liftover_chains"
NCBI_DIR="$WORKFLOW_ROOT/data/ncbi"
SRC_DIR="$WORKFLOW_ROOT/src/rainplot_workflow"
UTILS_DIR="$WORKFLOW_ROOT/src/utils"
TOOLS_DIR="$WORKFLOW_ROOT/tools"
RESULTS_DIR="$WORKFLOW_ROOT/results/rainplot_results"
GFF_FILE="$NCBI_DIR/sacCer3/genomic.gff"

LIFTOVER_BIN="$TOOLS_DIR/liftOver"
LIFTOVER_URL="https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver"

mkdir -p \
  "$BAM_DIR" \
  "$BED_DIR" \
  "$SORTED_BAM_DIR" \
  "$INDEX_DIR" \
  "$CHAIN_DIR" \
  "$NCBI_DIR" \
  "$TOOLS_DIR" \
  "$RESULTS_DIR"

# ------------------------------------------------------------------
# Expected positional arguments:
#   1 = BAM file
#   2 = chromosome
#   3 = start
#   4 = end
#   5 = read_id
#   6 = output .bed filename
#
# Example:
# bash src/rainplot_workflow/rainplot_workflow_script.sh \
#   Mitosis.sorted.bam CM007964.1 150000 200000 "" chr1_150000_to_200000bases.bed
# ------------------------------------------------------------------

BAM_FILE="$1"
CHROM="$2"
START="$3"
END="$4"
READ_ID="$5"
OUTPUT_NAME="${6:-${CHROM}_${START}_${END}_brdu.bed}"

if [ -z "$BAM_FILE" ] || [ -z "$CHROM" ] || [ -z "$START" ] || [ -z "$END" ]; then
  echo "[ERROR] Missing required arguments."
  echo "Usage: bash src/rainplot_workflow/rainplot_workflow_script.sh BAM CHROM START END [READ_ID] [OUTPUT_BED]"
  exit 1
fi

BAM="$BAM_DIR/$BAM_FILE"
OUTPUT="$BED_DIR/$OUTPUT_NAME"
RFB_OUTPUT="$BED_DIR/rfb_bases${START}_to_${END}.bed"

OUTPUT_BASENAME="$(basename "$OUTPUT_NAME" .bed)"
LIFTOVER_OUTPUT="$BED_DIR/liftover_${OUTPUT_BASENAME}.bed"
LIFTOVER_UNMAPPED_OUTPUT="$BED_DIR/liftover_${OUTPUT_BASENAME}_unmapped.bed"
RAINPLOT_MANIFEST="$RESULTS_DIR/${OUTPUT_BASENAME}_rainplots_manifest.txt"
GENOMIC_FEATURE_PNG="$RESULTS_DIR/genomic_feature_${OUTPUT_BASENAME}.png"

BED_TO_USE="$OUTPUT"
USE_RFB_OVERLAY="yes"
GENERATE_GENOMIC_FEATURES="yes"

if [ ! -f "$BAM" ]; then
  echo "[ERROR] BAM file not found: $BAM"
  exit 1
fi

# --- Auto-detect phase label from BAM filename ---
BAM_LOWER="$(echo "$BAM_FILE" | tr '[:upper:]' '[:lower:]')"
if [[ "$BAM_LOWER" == *mitosis* ]]; then
  PHASE_LABEL="Mitosis"
elif [[ "$BAM_LOWER" == *s_phase* ]] || [[ "$BAM_LOWER" == *sphase* ]] || [[ "$BAM_LOWER" == *s-phase* ]]; then
  PHASE_LABEL="S Phase"
else
  PHASE_LABEL="Unknown"
fi

echo "[INFO] Phase label detected: $PHASE_LABEL"

# --- Virtual environment setup ---
VENV_DIR="$SRC_DIR/.rainplot_env"
RENV_DIR="$SRC_DIR/.renv"

if [ ! -d "$VENV_DIR" ]; then
  echo "[INFO] Virtual environment not found. Creating $VENV_DIR..."
  python3 -m venv "$VENV_DIR"
  if [ $? -ne 0 ]; then
    echo "[ERROR] Failed to create virtual environment. Exiting."
    exit 1
  fi
  echo "[INFO] Virtual environment created."
fi

echo "[INFO] Activating virtual environment..."
source "$VENV_DIR/bin/activate"
if [ $? -ne 0 ]; then
  echo "[ERROR] Failed to activate virtual environment. Exiting."
  exit 1
fi

echo "[INFO] Installing required libraries from requirements.txt..."
pip install -r "$SRC_DIR/requirements.txt" --quiet
if [ $? -ne 0 ]; then
  echo "[ERROR] pip install failed. Exiting."
  deactivate
  exit 1
fi
echo "[INFO] Libraries installed successfully."

# --- Local R library setup ---
if [ ! -d "$RENV_DIR" ]; then
  echo "[INFO] renv project directory not found. Creating $RENV_DIR..."
  mkdir -p "$RENV_DIR"
  if [ $? -ne 0 ]; then
    echo "[ERROR] Failed to create renv project directory. Exiting."
    deactivate
    exit 1
  fi
  echo "[INFO] renv project directory created."
fi

echo "[INFO] Ensuring renv is initialized and restored for $SRC_DIR..."
Rscript "$SRC_DIR/ensure_r_environment.R" "$SRC_DIR"
if [ $? -ne 0 ]; then
  echo "[ERROR] renv setup failed. Exiting."
  deactivate
  exit 1
fi
echo "[INFO] renv environment ready."

# --- BrdU extraction ---
CMD="python $SRC_DIR/raw_data_extraction_on_bam.py $BAM -c $CHROM -s $START -e $END -o $OUTPUT"

if [ -n "$READ_ID" ]; then
  CMD="$CMD -r $READ_ID"
fi

echo "[INFO] Running BrdU extraction..."
eval "$CMD"

if [ ! -s "$OUTPUT" ]; then
  echo "[ERROR] Output file $OUTPUT is empty or was not created. Exiting."
  deactivate
  exit 1
fi

echo "[INFO] Extraction complete: $OUTPUT"

# --- Ask user if they want LiftOver ---
while true; do
  read -r -p "Would you like to do a liftover? [y/n]: " DO_LIFTOVER
  case "$DO_LIFTOVER" in
    [Yy] )
      if [ ! -x "$LIFTOVER_BIN" ]; then
        echo "[INFO] liftOver not found at $LIFTOVER_BIN"
        echo "[INFO] Downloading UCSC liftOver into $TOOLS_DIR ..."
        wget -O "$LIFTOVER_BIN" "$LIFTOVER_URL"
        if [ $? -ne 0 ]; then
          echo "[ERROR] Failed to download liftOver from $LIFTOVER_URL"
          deactivate
          exit 1
        fi
        chmod +x "$LIFTOVER_BIN"
        if [ $? -ne 0 ]; then
          echo "[ERROR] Failed to make liftOver executable."
          deactivate
          exit 1
        fi
        echo "[INFO] liftOver installed successfully at: $LIFTOVER_BIN"
      else
        echo "[INFO] Using existing liftOver binary: $LIFTOVER_BIN"
      fi

      while true; do
        read -r -p "Enter the full path to the chain file: " CHAIN_PATH
        if [ -z "$CHAIN_PATH" ]; then
          echo "[WARN] No chain file path entered. Please try again."
          continue
        fi
        if [ ! -f "$CHAIN_PATH" ]; then
          echo "[WARN] Chain file not found: $CHAIN_PATH"
          continue
        fi
        break
      done

      TEMP_LIFTOVER_INPUT="$(mktemp "$BED_DIR/.liftover_input_${OUTPUT_BASENAME}_XXXXXX.bed")"

      echo "[INFO] Converting BED chromosome names from GenBank to UCSC for LiftOver..."
      python "$UTILS_DIR/convert_genbank_bed_to_ucsc.py" "$OUTPUT" "$TEMP_LIFTOVER_INPUT"

      if [ ! -s "$TEMP_LIFTOVER_INPUT" ]; then
        echo "[ERROR] Temporary UCSC-converted BED file is empty or was not created."
        rm -f "$TEMP_LIFTOVER_INPUT"
        deactivate
        exit 1
      fi

      echo "[INFO] Running LiftOver..."
      python "$UTILS_DIR/liftover.py" \
        "$TEMP_LIFTOVER_INPUT" \
        --chain "$CHAIN_PATH" \
        --mapped "$LIFTOVER_OUTPUT" \
        --unmapped "$LIFTOVER_UNMAPPED_OUTPUT" \
        --liftover_bin "$LIFTOVER_BIN"

      rm -f "$TEMP_LIFTOVER_INPUT"

      if [ ! -s "$LIFTOVER_OUTPUT" ]; then
        echo "[ERROR] LiftOver output file is empty or was not created: $LIFTOVER_OUTPUT"
        deactivate
        exit 1
      fi

      BED_TO_USE="$LIFTOVER_OUTPUT"
      USE_RFB_OVERLAY="no"
      GENERATE_GENOMIC_FEATURES="no"

      echo "[INFO] LiftOver complete: $LIFTOVER_OUTPUT"
      echo "[INFO] Unmapped LiftOver intervals written to: $LIFTOVER_UNMAPPED_OUTPUT"
      echo "[INFO] The workflow will use the lifted BED file for plotting."
      echo "[INFO] RFB overlay will be skipped because the RFB BED file was not lifted and would be in a different coordinate system."
      echo "[INFO] Genomic feature annotation will also be skipped because the local GFF coordinates do not match the lifted coordinate system."
      break
      ;;
    [Nn] )
      echo "[INFO] LiftOver skipped. The workflow will use the original BED file for plotting."
      break
      ;;
    * )
      echo "[WARN] Invalid response. Please enter 'y' or 'n'."
      ;;
  esac
done

# --- RFB extraction ---
echo "[INFO] Running RFB motif extraction..."
python "$SRC_DIR/rfb_seq_matcher.py" \
  "$BAM" \
  -c "$CHROM" \
  -s "$START" \
  -e "$END" \
  -o "$RFB_OUTPUT"

if [ ! -s "$RFB_OUTPUT" ]; then
  echo "[WARN] RFB output file is empty or was not created: $RFB_OUTPUT"
else
  echo "[INFO] RFB extraction complete: $RFB_OUTPUT"
fi

# --- Generate rain plots ---
echo "[INFO] Generating rain plots..."

if [ "$USE_RFB_OVERLAY" = "yes" ]; then
  python "$SRC_DIR/rainplot_generation.py" "$BED_TO_USE" \
    -o "$RESULTS_DIR" \
    --region_start "$START" \
    --region_end "$END" \
    --rfb_dir "$BED_DIR" \
    --phase "$PHASE_LABEL" \
    --output_manifest "$RAINPLOT_MANIFEST"
else
  python "$SRC_DIR/rainplot_generation.py" "$BED_TO_USE" \
    -o "$RESULTS_DIR" \
    --region_start "$START" \
    --region_end "$END" \
    --phase "$PHASE_LABEL" \
    --output_manifest "$RAINPLOT_MANIFEST"
fi

if [ $? -ne 0 ]; then
  echo "[ERROR] Rain plot generation failed."
  deactivate
  exit 1
fi

echo "[INFO] Rain plots complete. Saved to: $RESULTS_DIR"

if [ ! -s "$RAINPLOT_MANIFEST" ]; then
  echo "[ERROR] Rain plot manifest is empty or was not created: $RAINPLOT_MANIFEST"
  deactivate
  exit 1
fi

if [ "$GENERATE_GENOMIC_FEATURES" = "yes" ]; then
  if [ ! -f "$GFF_FILE" ]; then
    echo "[ERROR] Genomic annotation file not found: $GFF_FILE"
    deactivate
    exit 1
  fi

  echo "[INFO] Generating genomic feature plot with Gviz..."
  Rscript "$SRC_DIR/genomic_feature_plot.R" \
    "$GFF_FILE" \
    "$CHROM" \
    "$START" \
    "$END" \
    "$GENOMIC_FEATURE_PNG"

  if [ $? -ne 0 ] || [ ! -s "$GENOMIC_FEATURE_PNG" ]; then
    echo "[ERROR] Genomic feature plot generation failed."
    deactivate
    exit 1
  fi

  echo "[INFO] Genomic feature plot complete: $GENOMIC_FEATURE_PNG"

  echo "[INFO] Combining rain plots with genomic feature plot..."
  python "$SRC_DIR/combine_rainplot_images.py" \
    --manifest "$RAINPLOT_MANIFEST" \
    --annotation_png "$GENOMIC_FEATURE_PNG" \
    --output_dir "$RESULTS_DIR" \
    --delete_inputs

  if [ $? -ne 0 ]; then
    echo "[ERROR] Combined plot generation failed."
    deactivate
    exit 1
  fi

  echo "[INFO] Combined rain plots saved to: $RESULTS_DIR"
else
  echo "[INFO] Combined genomic feature plots were skipped for this run."
fi

rm -f "$RAINPLOT_MANIFEST"

deactivate
echo "[INFO] Virtual environment deactivated."
