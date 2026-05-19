#!/bin/bash
#SBATCH --job-name=extract_brdu
#SBATCH --output=extract_brdu_%j.out
#SBATCH --error=extract_brdu_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --partition=normal

# Exit on errors, undefined variables, or failed commands inside pipes.
# This makes the workflow fail early instead of continuing with missing files.
set -euo pipefail

module load python/3.13.7
module load modkit
module load samtools
module load R/4.5.2

# Resolve the workflow root based on this script location so the workflow can be
# launched from different directories without breaking relative paths.
WORKFLOW_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"

BAM_DIR="$WORKFLOW_ROOT/data/bam"
BED_DIR="$WORKFLOW_ROOT/data/bed"
SORTED_BAM_DIR="$WORKFLOW_ROOT/data/sorted_bam"
INDEX_DIR="$WORKFLOW_ROOT/data/index_sorted_bam_bai"
CHAIN_DIR="$WORKFLOW_ROOT/data/liftover_chains"
NCBI_ROOT="$WORKFLOW_ROOT/data/ncbi"
SRC_DIR="$WORKFLOW_ROOT/src/rainplot_workflow"
UTILS_DIR="$WORKFLOW_ROOT/src/utils"
TOOLS_DIR="$WORKFLOW_ROOT/tools"
RESULTS_DIR="$WORKFLOW_ROOT/results/rainplot_results"

LIFTOVER_BIN="$TOOLS_DIR/liftOver"
LIFTOVER_URL="https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver"

# Create required workflow directories up front so later steps can write outputs
# without failing because a directory is missing.
mkdir -p \
  "$BAM_DIR" \
  "$BED_DIR" \
  "$SORTED_BAM_DIR" \
  "$INDEX_DIR" \
  "$CHAIN_DIR" \
  "$NCBI_ROOT" \
  "$TOOLS_DIR" \
  "$RESULTS_DIR"

BAM_FILE="${1:-}"
CHROM="${2:-}"
START="${3:-}"
END="${4:-}"
READ_ID="${5:-}"

if [ -z "$BAM_FILE" ] || [ -z "$CHROM" ] || [ -z "$START" ] || [ -z "$END" ]; then
  echo "[ERROR] Missing required arguments."
  echo "Usage: bash src/rainplot_workflow/rainplot_workflow_script.sh BAM CHROM START END [READ_ID] [OUTPUT_BED]"
  exit 1
fi

PHASE_LABEL=""
RFB_PLOT_MODE="none"
FILTER_RFB_READS="no"
REQUEST_RFB_OVERLAY="yes"
ENABLE_RFB_WORKFLOW="no"

while true; do
  read -r -p "Which phase is this run for? [M/S]: " PHASE_CHOICE

  case "$PHASE_CHOICE" in
    [Mm] )
      PHASE_LABEL="Mitosis"
      RFB_PLOT_MODE="none"
      FILTER_RFB_READS="no"
      REQUEST_RFB_OVERLAY="no"
      ENABLE_RFB_WORKFLOW="no"
      break
      ;;
    [Ss] )
      PHASE_LABEL="S Phase"
      ENABLE_RFB_WORKFLOW="yes"

      while true; do
        read -r -p "For S Phase, choose rain plot mode: [a] RFB coords only, [b] without RFB, [c] with and without RFB coords: " S_PHASE_PLOT_CHOICE

        case "$S_PHASE_PLOT_CHOICE" in
          [Aa] )
            RFB_PLOT_MODE="rfb_only"
            FILTER_RFB_READS="yes"
            REQUEST_RFB_OVERLAY="yes"
            break
            ;;
          [Bb] )
            RFB_PLOT_MODE="without_rfb"
            FILTER_RFB_READS="no"
            REQUEST_RFB_OVERLAY="no"
            break
            ;;
          [Cc] )
            RFB_PLOT_MODE="mixed"
            FILTER_RFB_READS="no"
            REQUEST_RFB_OVERLAY="yes"
            break
            ;;
          * )
            echo "[WARN] Invalid response. Please enter a, b, or c."
            ;;
        esac
      done

      break
      ;;
    * )
      echo "[WARN] Invalid response. Please enter M or S."
      ;;
  esac
done

echo "[INFO] Phase selected: $PHASE_LABEL"
echo "[INFO] RFB plot mode selected: $RFB_PLOT_MODE"
echo "[INFO] RFB workflow enabled: $ENABLE_RFB_WORKFLOW"

VENV_DIR="$SRC_DIR/.rainplot_env"
R_ENV_DIR="$SRC_DIR/.r_library"

# Use a workflow-local Python environment so package versions are controlled by
# this project instead of whatever happens to be installed on the HPC system.
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

echo "[INFO] Installing required Python libraries from requirements.txt..."
pip install -r "$SRC_DIR/requirements.txt" --quiet

if [ $? -ne 0 ]; then
  echo "[ERROR] pip install failed. Exiting."
  deactivate
  exit 1
fi

echo "[INFO] Python libraries installed successfully."

# Use a local R library for the same reason as the Python venv: reproducible
# package installs without needing write access to system-wide R libraries.
if [ ! -d "$R_ENV_DIR" ]; then
  echo "[INFO] Local R library not found. Creating $R_ENV_DIR..."
  mkdir -p "$R_ENV_DIR"

  if [ $? -ne 0 ]; then
    echo "[ERROR] Failed to create local R library directory. Exiting."
    deactivate
    exit 1
  fi

  echo "[INFO] Local R library directory created."
fi

export R_LIBS_USER="$R_ENV_DIR"
export WORKFLOW_ROOT="$WORKFLOW_ROOT"

chromosome_label_from_id() {
  local chrom_id="$1"

  # Convert common yeast chromosome identifiers into short labels used in output
  # filenames. This keeps output names readable even when the input uses GenBank
  # or RefSeq accession IDs.
  case "$chrom_id" in
    CM007964.1|NC_001133.9|1|chrI|chri|chr1) printf '%s\n' "chr1" ;;
    CM007965.1|NC_001134.8|2|chrII|chrii|chr2) printf '%s\n' "chr2" ;;
    CM007966.1|NC_001135.5|3|chrIII|chriii|chr3) printf '%s\n' "chr3" ;;
    CM007967.1|NC_001136.10|4|chrIV|chriv|chr4) printf '%s\n' "chr4" ;;
    CM007968.1|NC_001137.3|5|chrV|chrv|chr5) printf '%s\n' "chr5" ;;
    CM007969.1|NC_001138.5|6|chrVI|chrvi|chr6) printf '%s\n' "chr6" ;;
    CM007970.1|NC_001139.9|7|chrVII|chrvii|chr7) printf '%s\n' "chr7" ;;
    CM007971.1|NC_001140.6|8|chrVIII|chrviii|chr8) printf '%s\n' "chr8" ;;
    CM007972.1|NC_001141.2|9|chrIX|chrix|chr9) printf '%s\n' "chr9" ;;
    CM007973.1|NC_001142.9|10|chrX|chrx|chr10) printf '%s\n' "chr10" ;;
    CM007974.1|NC_001143.9|11|chrXI|chrxi|chr11) printf '%s\n' "chr11" ;;
    CM007975.1|NC_001144.5|12|chrXII|chrxii|chr12) printf '%s\n' "chr12" ;;
    CM007976.1|NC_001145.3|13|chrXIII|chrxiii|chr13) printf '%s\n' "chr13" ;;
    CM007977.1|NC_001146.8|14|chrXIV|chrxiv|chr14) printf '%s\n' "chr14" ;;
    CM007978.1|NC_001147.6|15|chrXV|chrxv|chr15) printf '%s\n' "chr15" ;;
    CM007979.1|NC_001148.4|16|chrXVI|chrxvi|chr16) printf '%s\n' "chr16" ;;
    *) printf '%s\n' "$chrom_id" ;;
  esac
}

resolve_liftover_chain() {
  local target_strain="$1"
  local expected_chain="$CHAIN_DIR/W303TosacCer${target_strain#sacCer}.over.chain.gz"

  if [ -f "$expected_chain" ]; then
    printf '%s\n' "$expected_chain"
    return 0
  fi

  return 1
}

resolve_annotation_file() {
  local target_strain="$1"

  # Pick the annotation source that matches the coordinate system being plotted.
  # W303 currently falls back to sacCer3 annotations unless liftover changes the
  # target strain.
  case "$target_strain" in
    W303) printf '%s\n' "$WORKFLOW_ROOT/data/ncbi/sacCer3/genomic.gff" ;;
    sacCer1) printf '%s\n' "$WORKFLOW_ROOT/data/ncbi/sacCer1/sacCer1_features.bed" ;;
    sacCer2) printf '%s\n' "$WORKFLOW_ROOT/data/ncbi/sacCer2/sacCer2_features.bed" ;;
    sacCer3) printf '%s\n' "$WORKFLOW_ROOT/data/ncbi/sacCer3/genomic.gff" ;;
    *) return 1 ;;
  esac
}

set_named_paths() {
  local basename_root="$1"

  # Keep all downstream filenames tied to one basename so the BED, rain plots,
  # manifest, and annotation PNG are easy to trace back to the same run.
  OUTPUT_NAME="${basename_root}.bed"
  OUTPUT="$BED_DIR/$OUTPUT_NAME"
  OUTPUT_BASENAME="$basename_root"
  LIFTOVER_OUTPUT="$BED_DIR/liftover_${OUTPUT_BASENAME}.bed"
  LIFTOVER_UNMAPPED_OUTPUT="$BED_DIR/liftover_${OUTPUT_BASENAME}_unmapped.bed"
  RAINPLOT_MANIFEST="$RESULTS_DIR/${OUTPUT_BASENAME}_rainplots_manifest.txt"
  GENOMIC_FEATURE_PNG="$RESULTS_DIR/genomic_feature_${OUTPUT_BASENAME}.png"
}

DEFAULT_CHROM_LABEL="$(chromosome_label_from_id "$CHROM")"
DEFAULT_OUTPUT_BASENAME="${DEFAULT_CHROM_LABEL}_${START}_to_${END}bases"
USER_OUTPUT_NAME="${6:-}"

if [ -n "$USER_OUTPUT_NAME" ]; then
  OUTPUT_BASENAME="$(basename "$USER_OUTPUT_NAME" .bed)"
else
  OUTPUT_BASENAME="$DEFAULT_OUTPUT_BASENAME"
fi

ANNOTATION_FILE="$(resolve_annotation_file "W303")"
ANNOTATION_LABEL="W303"
G4_BED_FILE="$WORKFLOW_ROOT/data/bed/g4.motifs.bed"
TARGET_STRAIN="W303"

BAM="$BAM_DIR/$BAM_FILE"
set_named_paths "$OUTPUT_BASENAME"
RFB_OUTPUT="$BED_DIR/rfb_bases${START}_to_${END}.bed"
BED_TO_USE="$OUTPUT"
USE_RFB_OVERLAY="yes"
GENERATE_GENOMIC_FEATURES="yes"
REGION_METADATA_FILE="$RESULTS_DIR/rainplot_regions.tsv"

if [ ! -f "$BAM" ]; then
  echo "[ERROR] BAM file not found: $BAM"
  exit 1
fi

echo "[INFO] Ensuring local R library is ready at $R_LIBS_USER..."
Rscript "$SRC_DIR/ensure_r_environment.R" "$R_LIBS_USER"

if [ $? -ne 0 ]; then
  echo "[ERROR] Local R library setup failed. Exiting."
  deactivate
  exit 1
fi

echo "[INFO] Local R environment ready."

# BrdU extraction stays in the original BAM coordinate system.
# This matters because W303 BAMs and sacCer annotations may not share coordinates
# unless a liftover step is explicitly performed later.
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

echo "[INFO] BrdU extraction complete: $OUTPUT"
echo "[INFO] BrdU data will remain in the BAM coordinate system."
echo "[INFO] For W303 BAMs, this means BrdU data remains in W303 coordinates."

while true; do
  read -r -p "Would you like to do a UCSC liftOver on the BrdU BED? [y/n]: " DO_LIFTOVER

  case "$DO_LIFTOVER" in
    [Yy] )
      # Download liftOver only when needed so the workflow does not require the
      # binary to be preinstalled.
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
        read -r -p "Which target yeast strain do you want to lift over to? [sacCer1/sacCer2/sacCer3]: " TARGET_STRAIN

        case "$TARGET_STRAIN" in
          sacCer1|sacCer2|sacCer3)
            CHAIN_PATH="$(resolve_liftover_chain "$TARGET_STRAIN" || true)"
            ANNOTATION_FILE="$(resolve_annotation_file "$TARGET_STRAIN" || true)"
            ANNOTATION_LABEL="$TARGET_STRAIN"

            if [ -z "$CHAIN_PATH" ]; then
              echo "[WARN] No chain file found for target strain $TARGET_STRAIN in $CHAIN_DIR"
              echo "[WARN] Expected W303TosacCer${TARGET_STRAIN#sacCer}.over.chain.gz"
              continue
            fi

            if [ -z "$ANNOTATION_FILE" ] || [ ! -s "$ANNOTATION_FILE" ]; then
              echo "[WARN] Annotation file not found for target strain $TARGET_STRAIN"
              continue
            fi

            echo "[INFO] Using chain file: $CHAIN_PATH"
            echo "[INFO] Using annotation file: $ANNOTATION_FILE"
            break
            ;;
          * )
            echo "[WARN] Invalid strain. Please enter sacCer1, sacCer2, or sacCer3."
            ;;
        esac
      done

      LIFTOVER_BASENAME="W303_to_${TARGET_STRAIN}_${OUTPUT_BASENAME}"
      LIFTOVER_OUTPUT="$BED_DIR/${LIFTOVER_BASENAME}.bed"
      LIFTOVER_UNMAPPED_OUTPUT="$BED_DIR/${LIFTOVER_BASENAME}_unmapped.bed"
      RAINPLOT_MANIFEST="$RESULTS_DIR/${LIFTOVER_BASENAME}_rainplots_manifest.txt"
      GENOMIC_FEATURE_PNG="$RESULTS_DIR/genomic_feature_${LIFTOVER_BASENAME}.png"

      echo "[INFO] Running chain-based liftOver on BrdU BED using W303 GenBank chromosome names..."
      python "$UTILS_DIR/liftover_brdu_bed.py" \
        "$OUTPUT" \
        --chain "$CHAIN_PATH" \
        --mapped "$LIFTOVER_OUTPUT" \
        --unmapped "$LIFTOVER_UNMAPPED_OUTPUT" \
        --liftover_bin "$LIFTOVER_BIN"

      if [ ! -s "$LIFTOVER_OUTPUT" ]; then
        echo "[ERROR] LiftOver output file is empty or was not created: $LIFTOVER_OUTPUT"
        deactivate
        exit 1
      fi

      BED_TO_USE="$LIFTOVER_OUTPUT"

      # RFB coordinates are extracted from the original BAM. After liftover, the
      # BrdU BED is in a different coordinate system, so overlaying the original
      # RFB BED would put the motif in the wrong location.
      USE_RFB_OVERLAY="no"

      GENERATE_GENOMIC_FEATURES="yes"

      # G4 motifs are disabled after chain liftover because the default G4 BED
      # may not match the lifted target strain coordinate system.
      G4_BED_FILE=""

      echo "[INFO] Chain-based LiftOver complete: $LIFTOVER_OUTPUT"
      echo "[INFO] Unmapped LiftOver intervals written to: $LIFTOVER_UNMAPPED_OUTPUT"
      echo "[INFO] The workflow will use the lifted BrdU BED file for plotting."
      echo "[INFO] RFB overlay is being skipped because the RFB BED remains in the original coordinate system."
      echo "[INFO] Genomic feature annotation will use the selected $TARGET_STRAIN annotation file."
      break
      ;;

    [Nn] )
      echo "[INFO] UCSC liftOver skipped."
      echo "[INFO] The workflow will continue in the original BAM/W303 coordinate system."
      break
      ;;

    * )
      echo "[WARN] Invalid response. Please enter 'y' or 'n'."
      ;;
  esac
done

if [ "$ENABLE_RFB_WORKFLOW" = "yes" ]; then
  echo "[INFO] Running RFB motif extraction for S phase..."
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
else
  echo "[INFO] Skipping RFB motif extraction because this run is not S phase."
fi

echo "[INFO] Generating rain plots..."
rm -f "$RAINPLOT_MANIFEST" "$REGION_METADATA_FILE"

RAINPLOT_FILENAME_PREFIX="${LIFTOVER_BASENAME:-$OUTPUT_BASENAME}"
SHOW_RFB_OVERLAY="no"

if [ "$ENABLE_RFB_WORKFLOW" = "yes" ] && [ "$REQUEST_RFB_OVERLAY" = "yes" ] && [ "$USE_RFB_OVERLAY" = "yes" ]; then
  SHOW_RFB_OVERLAY="yes"
elif [ "$ENABLE_RFB_WORKFLOW" = "yes" ] && [ "$REQUEST_RFB_OVERLAY" = "yes" ] && [ "$USE_RFB_OVERLAY" != "yes" ]; then
  echo "[WARN] RFB overlay was requested, but it is unavailable after liftOver because the RFB BED remains in the original coordinate system."
fi

run_rainplot_generation() {
  local manifest_path="$1"
  local filename_prefix="$2"
  local enable_rfb_overlay="$3"
  local filter_rfb_reads="$4"

  local cmd=(
    python "$SRC_DIR/rainplot_generation.py" "$BED_TO_USE"
    -o "$RESULTS_DIR"
    --region_start "$START"
    --region_end "$END"
    --phase "$PHASE_LABEL"
    --filename_prefix "$filename_prefix"
    --output_manifest "$manifest_path"
  )

  if [ "$ENABLE_RFB_WORKFLOW" = "yes" ] && { [ "$enable_rfb_overlay" = "yes" ] || [ "$filter_rfb_reads" = "yes" ]; }; then
    cmd+=(--rfb_dir "$BED_DIR")
  fi

  if [ "$filter_rfb_reads" = "yes" ]; then
    cmd+=(--filter_reads_with_rfb)
  fi

  if [ "$enable_rfb_overlay" = "yes" ]; then
    cmd+=(--show_rfb_overlay)
  fi

  "${cmd[@]}"
}

run_rainplot_generation "$RAINPLOT_MANIFEST" "$RAINPLOT_FILENAME_PREFIX" "$SHOW_RFB_OVERLAY" "$FILTER_RFB_READS"

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
  if [ ! -s "$ANNOTATION_FILE" ]; then
    echo "[ERROR] Annotation file not found:"
    echo "$ANNOTATION_FILE"
    deactivate
    exit 1
  fi

  echo "[INFO] Generating genomic feature plot from selected annotation source..."
  echo "[INFO] Using annotation file: $ANNOTATION_FILE"

  # Pass G4_BED_FILE only for this Rscript call so the R script can decide
  # whether to draw G4 motifs without changing the global shell environment.
  G4_BED_FILE="$G4_BED_FILE" Rscript "$SRC_DIR/genomic_feature_plot.R" \
    "$CHROM" \
    "$START" \
    "$END" \
    "$GENOMIC_FEATURE_PNG" \
    "$ANNOTATION_FILE" \
    "$ANNOTATION_LABEL"

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

# The manifest is only an intermediate file used by the image-combining step.
rm -f "$RAINPLOT_MANIFEST"

deactivate
echo "[INFO] Virtual environment deactivated."
echo "[INFO] Workflow complete."
