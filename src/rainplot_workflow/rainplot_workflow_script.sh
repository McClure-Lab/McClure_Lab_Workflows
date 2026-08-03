#!/bin/bash

# ==============================================================================
# RAINPLOT BrdU WORKFLOW
# ==============================================================================
#
# Purpose
# -------
# This script is the SLURM wrapper for the BrdU rainplot workflow. It supports
# both interactive job submission and non-interactive compute-job execution.
#
# The workflow can:
#
#   1. Select a BAM file from workflow_root/data/bam.
#   2. Extract BrdU calls for:
#        - A chromosome interval.
#        - One specific read.
#        - A ranked or predefined set of reads supplied in a TSV/list file.
#   3. Calculate BrdU signal in 100-thymidine windows using:
#        - Binary mode: threshold each probability before averaging.
#        - Mean mode: average the raw probabilities directly.
#   4. Limit the number of reads plotted.
#   5. Generate S-phase RFB-specific rainplots.
#   6. Optionally lift W303 BED coordinates to sacCer1, sacCer2, or sacCer3.
#   7. Generate genomic-feature annotation plots.
#   8. Combine each rainplot with its genomic-feature panel.
#
# Execution modes
# ---------------
#
# Submission mode:
#
#   bash src/rainplot_workflow/rainplot_workflow_script.sh
#
#   The script prompts for input values, validates them, and submits itself to
#   SLURM using the internal --run-job option.
#
# Job mode:
#
#   --run-job
#
#   This mode is normally entered automatically by the submitted SLURM job.
#   It performs extraction, optional liftOver, RFB processing, rainplot
#   generation, annotation generation, and image combination.
#
# Expected project structure
# --------------------------
#
# workflow_root/
# ├── data/
# │   ├── bam/
# │   ├── bed/
# │   ├── liftover_chains/
# │   └── ncbi/
# ├── logs/
# │   └── rainplot_workflow/
# ├── results/
# │   └── rainplot_results/
# ├── src/
# │   ├── rainplot_workflow/
# │   │   ├── rainplot_workflow_script.sh
# │   │   ├── raw_data_extraction_on_bam.py
# │   │   ├── rainplot_generation.py
# │   │   ├── rfb_seq_matcher.py
# │   │   ├── genomic_feature_plot.R
# │   │   ├── combine_rainplot_images.py
# │   │   ├── ensure_r_environment.R
# │   │   └── requirements.txt
# │   └── utils/
# │       └── liftover_brdu_bed.py
# └── tools/
#
# Important coordinate-system behavior
# ------------------------------------
#
# Without liftOver:
#   BrdU data remains in the coordinate system used by the BAM file. For the
#   W303 BAM files used by this workflow, this means W303 coordinates.
#
# With liftOver:
#   BrdU BED intervals are converted from W303 coordinates to a selected
#   sacCer assembly. The RFB BED remains in its original coordinate system,
#   so RFB overlays are disabled after liftOver to avoid mixing coordinate
#   systems.
#
# Output naming
# -------------
#
# Output names include the selected BrdU T-window mode:
#
#   *.binary.bed
#   *.mean.bed
#
# The generated image names and manifest names use the same output basename.
#
# ==============================================================================


# ==============================================================================
# SLURM RESOURCE REQUESTS
# ==============================================================================
#
# Requested resources:
#   - 4 CPU cores
#   - 8 GB RAM
#   - 2-hour runtime limit
#   - normal partition
#
#SBATCH --job-name=rainplot_brdu
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --partition=normal


# ==============================================================================
# BASH SAFETY SETTINGS
# ==============================================================================
#
# -e
#   Stop the script when a command fails.
#
# -u
#   Treat references to unset variables as errors.
#
# -o pipefail
#   Return a failure when any command in a pipeline fails.
#
set -euo pipefail


# ==============================================================================
# RESOLVE THE WORKFLOW ROOT
# ==============================================================================
#
# The workflow root is selected in this order:
#
#   1. --workflow-root supplied anywhere on the command line.
#   2. RAINPLOT_WORKFLOW_ROOT exported in the environment.
#   3. Two directories above this script.
#
# The third option assumes this file is located at:
#
#   workflow_root/src/rainplot_workflow/rainplot_workflow_script.sh
#

# Store an explicitly supplied workflow-root argument.
WORKFLOW_ROOT_ARG=""

# Search all command-line arguments for --workflow-root.
for ((arg_i = 1; arg_i <= $#; arg_i++)); do
  if [[ "${!arg_i}" == "--workflow-root" ]]; then
    next_arg_i=$((arg_i + 1))
    WORKFLOW_ROOT_ARG="${!next_arg_i:-}"
    break
  fi
done

# Resolve the final workflow-root path.
if [[ -n "$WORKFLOW_ROOT_ARG" ]]; then
  WORKFLOW_ROOT="$(cd "$WORKFLOW_ROOT_ARG" && pwd)"
elif [[ -n "${RAINPLOT_WORKFLOW_ROOT:-}" ]]; then
  WORKFLOW_ROOT="$RAINPLOT_WORKFLOW_ROOT"
else
  WORKFLOW_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
fi


# ==============================================================================
# WORKFLOW PATH CONFIGURATION
# ==============================================================================
#
# Centralized paths used for inputs, intermediate files, software, logs, and
# final results.
#

# Directory containing input BAM files.
BAM_DIR="$WORKFLOW_ROOT/data/bam"

# Directory containing generated BED files, RFB files, and other BED resources.
BED_DIR="$WORKFLOW_ROOT/data/bed"

# Directory containing W303-to-sacCer liftOver chain files.
CHAIN_DIR="$WORKFLOW_ROOT/data/liftover_chains"

# Root directory containing NCBI or UCSC genome annotations.
NCBI_ROOT="$WORKFLOW_ROOT/data/ncbi"

# Directory containing the rainplot workflow scripts.
SRC_DIR="$WORKFLOW_ROOT/src/rainplot_workflow"

# Absolute path to this wrapper script.
#
# The script submits itself back to SLURM with --run-job.
SCRIPT_PATH="$SRC_DIR/rainplot_workflow_script.sh"

# Directory containing shared workflow utility scripts.
UTILS_DIR="$WORKFLOW_ROOT/src/utils"

# Directory containing locally installed command-line tools.
TOOLS_DIR="$WORKFLOW_ROOT/tools"

# Directory containing final rainplot images.
RESULTS_DIR="$WORKFLOW_ROOT/results/rainplot_results"

# Directory containing SLURM and workflow logs.
LOG_DIR="$WORKFLOW_ROOT/logs/rainplot_workflow"

# UCSC liftOver executable and official download location.
#
# The executable is downloaded automatically only when liftOver is requested
# and the local executable is missing.
LIFTOVER_BIN="$TOOLS_DIR/liftOver"
LIFTOVER_URL="https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver"


# ==============================================================================
# FUNCTION: usage
# ==============================================================================
#
# Purpose:
#   Print the supported command-line syntax.
#
# Positional arguments:
#   $1
#       BAM filename.
#
#   $2
#       Chromosome ID.
#
#   $3
#       Start coordinate.
#
#   $4
#       End coordinate.
#
#   $5
#       Optional single read ID.
#
#   $6
#       Optional output BED filename.
#
#   $7
#       Optional TSV or list file containing read IDs.
#
usage() {
  echo "Usage: bash src/rainplot_workflow/rainplot_workflow_script.sh [BAM] [CHROM] [START] [END] [READ_ID] [OUTPUT_BED] [READ_IDS_FILE]"
  echo
  echo "BAM is resolved by filename under data/bam."
  echo "If arguments are omitted, the workflow prompts for inputs before submitting a SLURM job."
  echo "READ_ID may be empty to include all reads in the requested region."
  echo "READ_IDS_FILE may be a TSV with a read_id column to plot a ranked read set in one job."
  echo "BrdU T-window mode is prompted interactively and defaults to binary."
  echo "Rainplot read count is prompted interactively; blank plots all reads in the region."
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
#   Convert an existing file path into a normalized absolute path.
#
# Arguments:
#   $1
#       Path to an existing file.
#
# Output:
#   Absolute path written to standard output.
#
absolute_existing_file() {
  local path="$1"
  local dir
  local file

  # Convert the containing directory into an absolute path.
  dir="$(cd "$(dirname "$path")" && pwd)"

  # Preserve the original filename.
  file="$(basename "$path")"

  # Print the complete absolute file path.
  printf '%s/%s\n' "$dir" "$file"
}


# ==============================================================================
# FUNCTION: resolve_bam_file
# ==============================================================================
#
# Purpose:
#   Resolve a BAM by filename under workflow_root/data/bam.
#
# Arguments:
#   $1
#       BAM filename or path supplied by the user.
#
# Behavior:
#   Only the basename is used. Supplying a BAM path outside data/bam does not
#   bypass the workflow's BAM-directory restriction.
#
# Output:
#   Absolute BAM path when the file exists.
#
# Return status:
#   0 when the BAM file is found.
#   1 when the BAM file is not found.
#
resolve_bam_file() {
  local bam_input="$1"
  local bam_filename
  local bam_path

  # Remove directory components from the supplied BAM value.
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
# FUNCTION: require_existing_dir
# ==============================================================================
#
# Purpose:
#   Verify that a required directory is available on the compute node.
#
# Arguments:
#   $1
#       Human-readable directory name.
#
#   $2
#       Directory path.
#
# Behavior:
#   The workflow exits when the directory is missing.
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
# FUNCTION: load_rainplot_modules
# ==============================================================================
#
# Purpose:
#   Load the software required by the rainplot workflow.
#
# Modules:
#   - Python 3.13.7, with a fallback to the default Python module.
#   - Modkit.
#   - Samtools.
#   - R 4.5.2, with a fallback to the default R module.
#
# Behavior when modules are unavailable:
#   The workflow continues and expects python, modkit, samtools, and Rscript to
#   already be available through PATH.
#
load_rainplot_modules() {
  if command -v module >/dev/null 2>&1; then
    module load python/3.13.7 || module load python
    module load modkit
    module load samtools
    module load R/4.5.2 || module load R
  else
    echo "[WARN] Environment modules are not available in this shell."
    echo "[WARN] Continuing with python, samtools, modkit, and Rscript from PATH."
  fi
}


# ==============================================================================
# FUNCTION: chromosome_label_from_id
# ==============================================================================
#
# Purpose:
#   Convert accepted W303, sacCer, numeric, and Roman-numeral chromosome
#   identifiers into standardized chr1 through chr16 labels.
#
# Supported identifier types:
#   - W303 GenBank accessions.
#   - sacCer RefSeq accessions.
#   - Numeric chromosome values.
#   - Roman-numeral chromosome labels.
#   - Existing chr1 through chr16 labels.
#
# Arguments:
#   $1
#       Chromosome identifier.
#
# Output:
#   Standardized chromosome label.
#
# Notes:
#   Unrecognized identifiers are returned unchanged.
#
chromosome_label_from_id() {
  local chrom_id="$1"

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


# ==============================================================================
# FUNCTION: normalize_phase_label
# ==============================================================================
#
# Purpose:
#   Normalize accepted M-phase or S-phase values to the labels used in plots.
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
#   0 for a valid phase.
#   1 for an invalid phase.
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
# FUNCTION: normalize_brdu_window_mode
# ==============================================================================
#
# Purpose:
#   Select how BrdU probabilities are summarized within each 100-thymidine
#   window.
#
# Binary mode:
#   Each probability is converted to 0 or 1 using the selected probability
#   threshold. The binary calls are then averaged within each window.
#
# Mean mode:
#   Raw BrdU probabilities are averaged directly within each window.
#
# Accepted binary values:
#   blank
#   binary
#   b
#   1
#
# Accepted mean values:
#   mean
#   m
#   2
#
# Output:
#   "binary" or "mean"
#
normalize_brdu_window_mode() {
  local mode_input="$1"

  case "${mode_input,,}" in
    ""|binary|b|1)
      printf '%s\n' "binary"
      ;;
    mean|m|2)
      printf '%s\n' "mean"
      ;;
    *)
      return 1
      ;;
  esac
}


# ==============================================================================
# FUNCTION: normalize_probability_threshold
# ==============================================================================
#
# Purpose:
#   Validate a BrdU probability threshold from 0 through 1.
#
# Accepted examples:
#   0
#   0.5
#   .5
#   0.7
#   1
#   1.0
#
# Output:
#   Normalized numeric threshold.
#
# Return status:
#   0 for a valid threshold.
#   1 for an invalid threshold.
#
normalize_probability_threshold() {
  local threshold_input="$1"
  local threshold_value="$threshold_input"

  # Convert values such as ".5" to "0.5".
  if [[ "$threshold_value" =~ ^\.[0-9]+$ ]]; then
    threshold_value="0$threshold_value"
  fi

  # Require a valid numeric value from 0 through 1.
  if ! [[ "$threshold_value" =~ ^(0(\.[0-9]+)?|1(\.0+)?)$ ]]; then
    return 1
  fi

  # Perform a final numeric range check and normalize the formatting.
  awk -v value="$threshold_value" 'BEGIN {
    if (value < 0 || value > 1) {
      exit 1
    }
    printf "%.12g\n", value + 0
  }'
}


# ==============================================================================
# FUNCTION: normalize_max_reads
# ==============================================================================
#
# Purpose:
#   Validate the optional maximum number of reads to plot.
#
# Blank value:
#   Plot every selected read.
#
# Positive integer:
#   Plot no more than the specified number of reads.
#
# Output:
#   Blank text or the validated positive integer.
#
# Return status:
#   0 for a valid value.
#   1 for an invalid value.
#
normalize_max_reads() {
  local max_reads_input="$1"

  if [[ -z "$max_reads_input" ]]; then
    printf '%s\n' ""
    return 0
  fi

  if ! [[ "$max_reads_input" =~ ^[1-9][0-9]*$ ]]; then
    return 1
  fi

  printf '%s\n' "$max_reads_input"
}


# ==============================================================================
# FUNCTION: resolve_liftover_chain
# ==============================================================================
#
# Purpose:
#   Resolve a W303-to-sacCer liftOver chain file.
#
# Arguments:
#   $1
#       Target strain: sacCer1, sacCer2, or sacCer3.
#
# Expected filename:
#
#   W303TosacCer<version>.over.chain.gz
#
# Examples:
#
#   W303TosacCer1.over.chain.gz
#   W303TosacCer2.over.chain.gz
#   W303TosacCer3.over.chain.gz
#
# Output:
#   Chain-file path when found.
#
# Return status:
#   0 when found.
#   1 when missing.
#
resolve_liftover_chain() {
  local target_strain="$1"
  local expected_chain="$CHAIN_DIR/W303TosacCer${target_strain#sacCer}.over.chain.gz"

  if [[ -f "$expected_chain" ]]; then
    printf '%s\n' "$expected_chain"
    return 0
  fi

  return 1
}


# ==============================================================================
# FUNCTION: resolve_annotation_file
# ==============================================================================
#
# Purpose:
#   Return the genomic annotation file associated with W303 or the selected
#   sacCer target assembly.
#
# Arguments:
#   $1
#       W303, sacCer1, sacCer2, or sacCer3.
#
# Output:
#   Expected annotation-file path.
#
# Notes:
#   File existence is validated later by the calling workflow.
#
resolve_annotation_file() {
  local target_strain="$1"

  case "$target_strain" in
    W303)
      printf '%s\n' "$WORKFLOW_ROOT/data/ncbi/sacCer3/genomic.gff"
      ;;
    sacCer1)
      printf '%s\n' "$WORKFLOW_ROOT/data/ncbi/sacCer1/sacCer1_features.bed"
      ;;
    sacCer2)
      printf '%s\n' "$WORKFLOW_ROOT/data/ncbi/sacCer2/sacCer2_features.bed"
      ;;
    sacCer3)
      printf '%s\n' "$WORKFLOW_ROOT/data/ncbi/sacCer3/genomic.gff"
      ;;
    *)
      return 1
      ;;
  esac
}


# ==============================================================================
# FUNCTION: set_named_paths
# ==============================================================================
#
# Purpose:
#   Derive all major intermediate and output paths from one output basename.
#
# Arguments:
#   $1
#       Output basename without a .bed extension.
#
# Global variables assigned:
#   OUTPUT_NAME
#   OUTPUT
#   OUTPUT_BASENAME
#   LIFTOVER_OUTPUT
#   LIFTOVER_UNMAPPED_OUTPUT
#   RAINPLOT_MANIFEST
#   GENOMIC_FEATURE_PNG
#
set_named_paths() {
  local basename_root="$1"

  OUTPUT_NAME="${basename_root}.bed"
  OUTPUT="$BED_DIR/$OUTPUT_NAME"
  OUTPUT_BASENAME="$basename_root"
  LIFTOVER_OUTPUT="$BED_DIR/liftover_${OUTPUT_BASENAME}.bed"
  LIFTOVER_UNMAPPED_OUTPUT="$BED_DIR/liftover_${OUTPUT_BASENAME}_unmapped.bed"
  RAINPLOT_MANIFEST="$RESULTS_DIR/${OUTPUT_BASENAME}_rainplots_manifest.txt"
  GENOMIC_FEATURE_PNG="$RESULTS_DIR/genomic_feature_${OUTPUT_BASENAME}.png"
}


# ==============================================================================
# FUNCTION: default_output_basename
# ==============================================================================
#
# Purpose:
#   Build the default basename for a chromosome-interval run.
#
# Arguments:
#   $1
#       BAM path.
#
#   $2
#       Chromosome identifier.
#
#   $3
#       Start coordinate.
#
#   $4
#       End coordinate.
#
# Output format:
#
#   <bam_prefix>.<chromosome>.<start>.<end>based
#
default_output_basename() {
  local bam_path="$1"
  local chrom="$2"
  local start="$3"
  local end="$4"
  local bam_name
  local bam_prefix
  local chrom_label

  bam_name="$(basename "$bam_path")"
  bam_prefix="${bam_name%.bam}"
  chrom_label="$(chromosome_label_from_id "$chrom")"

  printf '%s.%s.%s.%sbased\n' \
    "$bam_prefix" \
    "$chrom_label" \
    "$start" \
    "$end"
}


# ==============================================================================
# FUNCTION: default_read_ids_output_basename
# ==============================================================================
#
# Purpose:
#   Build the default basename for a read-ID-file run.
#
# Arguments:
#   $1
#       BAM path.
#
#   $2
#       Optional maximum read count.
#
# Output when max_reads is set:
#
#   <bam_prefix>.top<max_reads>
#
# Output when max_reads is blank:
#
#   <bam_prefix>.read_ids
#
default_read_ids_output_basename() {
  local bam_path="$1"
  local max_reads="$2"
  local bam_name
  local bam_prefix

  bam_name="$(basename "$bam_path")"
  bam_prefix="${bam_name%.bam}"

  if [[ -n "$max_reads" ]]; then
    printf '%s.top%s\n' "$bam_prefix" "$max_reads"
  else
    printf '%s.read_ids\n' "$bam_prefix"
  fi
}


# ==============================================================================
# FUNCTION: append_brdu_mode_to_basename
# ==============================================================================
#
# Purpose:
#   Add the selected BrdU window mode to the output basename.
#
# Arguments:
#   $1
#       Base output name.
#
#   $2
#       binary or mean.
#
# Behavior:
#   If the basename already ends in .binary or .mean, it is returned unchanged.
#   Otherwise, the selected mode is appended.
#
append_brdu_mode_to_basename() {
  local basename_root="$1"
  local brdu_mode="$2"

  case "$basename_root" in
    *.binary|*.mean)
      printf '%s\n' "$basename_root"
      ;;
    *)
      printf '%s.%s\n' "$basename_root" "$brdu_mode"
      ;;
  esac
}


# ==============================================================================
# FUNCTION: run_job
# ==============================================================================
#
# Purpose:
#   Execute the compute-intensive portion of the workflow inside the SLURM
#   allocation.
#
# Major stages:
#   1. Parse named job arguments.
#   2. Validate BrdU mode, probability threshold, and read count.
#   3. Validate BAM and read-list inputs on the compute node.
#   4. Prepare isolated Python and R environments.
#   5. Extract BrdU BED data.
#   6. Optionally lift W303 coordinates into a sacCer assembly.
#   7. Optionally extract and use RFB motifs for S-phase runs.
#   8. Generate rainplots.
#   9. Generate genomic-feature panels.
#  10. Combine rainplots and annotation panels.
#
run_job() {
  # ---------------------------------------------------------------------------
  # Job arguments and defaults
  # ---------------------------------------------------------------------------

  local bam_path=""
  local chrom=""
  local start=""
  local end=""
  local read_id=""
  local read_ids_file=""
  local output_basename=""
  local phase_label=""

  # Default BrdU window calculation mode.
  local brdu_window_mode="binary"
  local raw_brdu_window_mode

  # Default probability threshold used in binary mode.
  local binary_threshold="0.5"
  local raw_binary_threshold

  # Blank means all selected reads.
  local max_reads=""
  local raw_max_reads

  # RFB-related settings.
  local rfb_plot_mode=""
  local filter_rfb_reads=""
  local request_rfb_overlay=""
  local enable_rfb_workflow=""

  # liftOver-related settings.
  local do_liftover=""
  local target_strain="W303"

  # Annotation defaults.
  local annotation_file
  local annotation_label="W303"
  local g4_bed_file="$WORKFLOW_ROOT/data/bed/g4.motifs.bed"

  # Runtime paths and command arrays.
  local cmd
  local r_env_dir
  local venv_dir
  local workflow_log

  # Tracks whether selection is controlled by a read IDs file.
  local uses_read_ids_file="no"

  # Remove the leading --run-job argument.
  shift

  # ---------------------------------------------------------------------------
  # Parse named arguments supplied by submit_workflow()
  # ---------------------------------------------------------------------------
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --workflow-root)
        shift 2
        ;;
      --bam)
        bam_path="$2"
        shift 2
        ;;
      --chrom)
        chrom="$2"
        shift 2
        ;;
      --start)
        start="$2"
        shift 2
        ;;
      --end)
        end="$2"
        shift 2
        ;;
      --read-id)
        read_id="$2"
        shift 2
        ;;
      --read-ids-file)
        read_ids_file="$2"
        shift 2
        ;;
      --output-basename)
        output_basename="$2"
        shift 2
        ;;
      --phase-label)
        phase_label="$2"
        shift 2
        ;;
      --brdu-window-mode)
        brdu_window_mode="$2"
        shift 2
        ;;
      --binary-threshold)
        binary_threshold="$2"
        shift 2
        ;;
      --max-reads)
        max_reads="$2"
        shift 2
        ;;
      --rfb-plot-mode)
        rfb_plot_mode="$2"
        shift 2
        ;;
      --filter-rfb-reads)
        filter_rfb_reads="$2"
        shift 2
        ;;
      --request-rfb-overlay)
        request_rfb_overlay="$2"
        shift 2
        ;;
      --enable-rfb-workflow)
        enable_rfb_workflow="$2"
        shift 2
        ;;
      --do-liftover)
        do_liftover="$2"
        shift 2
        ;;
      --target-strain)
        target_strain="$2"
        shift 2
        ;;
      *)
        echo "[ERROR] Unknown job argument: $1"
        exit 1
        ;;
    esac
  done

  # ---------------------------------------------------------------------------
  # Validate the BrdU calculation settings
  # ---------------------------------------------------------------------------

  raw_brdu_window_mode="$brdu_window_mode"

  if ! brdu_window_mode="$(normalize_brdu_window_mode "$raw_brdu_window_mode")"; then
    echo "[ERROR] Invalid BrdU T-window mode in job: $raw_brdu_window_mode"
    echo "[ERROR] Expected binary or mean."
    exit 1
  fi

  raw_binary_threshold="$binary_threshold"

  if ! binary_threshold="$(normalize_probability_threshold "$raw_binary_threshold")"; then
    echo "[ERROR] Invalid binary BrdU probability threshold in job: $raw_binary_threshold"
    echo "[ERROR] Expected a number from 0 to 1, for example 0.5 or 0.7."
    exit 1
  fi

  # Ensure the output name identifies the selected BrdU window mode.
  output_basename="$(
    append_brdu_mode_to_basename \
      "$output_basename" \
      "$brdu_window_mode"
  )"

  raw_max_reads="$max_reads"

  if ! max_reads="$(normalize_max_reads "$raw_max_reads")"; then
    echo "[ERROR] Invalid max reads in job: $raw_max_reads"
    echo "[ERROR] Expected a positive integer, or blank for all reads."
    exit 1
  fi

  # A single-read filter and a read IDs file are mutually exclusive.
  if [[ -n "$read_id" && -n "$read_ids_file" ]]; then
    echo "[ERROR] Use either --read-id or --read-ids-file, not both."
    exit 1
  fi

  if [[ -n "$read_ids_file" ]]; then
    uses_read_ids_file="yes"
  fi

  # ---------------------------------------------------------------------------
  # Configure combined workflow logging
  # ---------------------------------------------------------------------------
  #
  # From this point onward, standard output and standard error are:
  #
  #   - Written to the SLURM output/error streams.
  #   - Appended to the workflow-specific log through tee.
  #
  mkdir -p "$LOG_DIR"

  workflow_log="$LOG_DIR/${output_basename}.${SLURM_JOB_ID:-manual}.workflow.log"

  exec > >(tee -a "$workflow_log") 2>&1

  echo "[INFO] Rainplot workflow started: $(date)"
  echo "[INFO] Workflow log: $workflow_log"

  # ---------------------------------------------------------------------------
  # Validate required directories
  # ---------------------------------------------------------------------------

  require_existing_dir "Workflow root" "$WORKFLOW_ROOT"
  require_existing_dir "Source directory" "$SRC_DIR"
  require_existing_dir "BAM directory" "$BAM_DIR"

  # Create output and tool directories when missing.
  mkdir -p \
    "$BAM_DIR" \
    "$BED_DIR" \
    "$CHAIN_DIR" \
    "$NCBI_ROOT" \
    "$TOOLS_DIR" \
    "$RESULTS_DIR" \
    "$LOG_DIR"

  # ---------------------------------------------------------------------------
  # Validate required input files on the compute node
  # ---------------------------------------------------------------------------

  if [[ ! -f "$bam_path" ]]; then
    echo "[ERROR] BAM file not found on compute node: $bam_path"
    exit 1
  fi

  if [[ -n "$read_ids_file" && ! -f "$read_ids_file" ]]; then
    echo "[ERROR] Read IDs file not found on compute node: $read_ids_file"
    exit 1
  fi

  # ---------------------------------------------------------------------------
  # Load software and prepare isolated Python and R environments
  # ---------------------------------------------------------------------------

  load_rainplot_modules

  # Python virtual environment used by the workflow.
  venv_dir="$SRC_DIR/.rainplot_env"

  # User-local R package library used by the workflow.
  r_env_dir="$SRC_DIR/.r_library"

  # Create the Python virtual environment only when its Python executable is
  # missing.
  if [[ ! -x "$venv_dir/bin/python" ]]; then
    echo "[INFO] Creating virtual environment: $venv_dir"
    python3 -m venv --clear "$venv_dir"
  fi

  echo "[INFO] Activating virtual environment."
  source "$venv_dir/bin/activate"

  # Install or verify all Python packages required by the workflow.
  echo "[INFO] Installing Python requirements."
  python -m pip install -r "$SRC_DIR/requirements.txt" --quiet

  # Create the local R package directory when needed.
  if [[ ! -d "$r_env_dir" ]]; then
    echo "[INFO] Creating local R library: $r_env_dir"
    mkdir -p "$r_env_dir"
  fi

  # Tell R to install and load packages from the workflow-local library.
  export R_LIBS_USER="$r_env_dir"

  # Export the workflow root for child Python and R scripts.
  export WORKFLOW_ROOT="$WORKFLOW_ROOT"

  # ---------------------------------------------------------------------------
  # Initialize default W303 annotation and output paths
  # ---------------------------------------------------------------------------

  annotation_file="$(resolve_annotation_file "W303")"

  set_named_paths "$output_basename"

  # RFB BED output for the requested interval.
  local rfb_output="$BED_DIR/rfb_bases${start}_to_${end}.bed"

  # The original BrdU BED is used unless liftOver is requested.
  local bed_to_use="$OUTPUT"

  # RFB overlays are allowed until disabled by incompatible input modes.
  local use_rfb_overlay="yes"

  # Genomic features are generated by default.
  local generate_genomic_features="yes"

  # Metadata file that may be created during rainplot generation.
  local region_metadata_file="$RESULTS_DIR/rainplot_regions.tsv"

  # Used for lifted output names.
  local liftover_basename=""

  # Resolved liftOver chain path.
  local chain_path=""

  # ---------------------------------------------------------------------------
  # Report the selected job configuration
  # ---------------------------------------------------------------------------

  echo "[INFO] Phase selected: $phase_label"
  echo "[INFO] BrdU T-window mode selected: $brdu_window_mode"

  if [[ "$brdu_window_mode" == "binary" ]]; then
    echo "[INFO] Binary BrdU probability threshold selected: $binary_threshold"
  fi

  echo "[INFO] RFB plot mode selected: $rfb_plot_mode"
  echo "[INFO] RFB workflow enabled: $enable_rfb_workflow"
  echo "[INFO] BAM: $bam_path"

  if [[ "$uses_read_ids_file" == "yes" && -z "$chrom" && -z "$start" && -z "$end" ]]; then
    echo "[INFO] Region: all coordinates for reads listed in read IDs file"
  else
    echo "[INFO] Region: $chrom:$start-$end"
  fi

  echo "[INFO] Read IDs file: ${read_ids_file:-none}"
  echo "[INFO] Rainplot read count: ${max_reads:-all selected reads}"
  echo "[INFO] Output BED: $OUTPUT"

  # ---------------------------------------------------------------------------
  # Prepare the local R package environment
  # ---------------------------------------------------------------------------

  echo "[INFO] Ensuring local R library is ready at $R_LIBS_USER..."

  Rscript \
    "$SRC_DIR/ensure_r_environment.R" \
    "$R_LIBS_USER"

  echo "[INFO] Local R environment ready."

  # ---------------------------------------------------------------------------
  # Build the BrdU extraction command
  # ---------------------------------------------------------------------------
  #
  # raw_data_extraction_on_bam.py extracts read-level BrdU information from the
  # selected BAM and writes it to a BED file.
  #
  cmd=(
    python
    "$SRC_DIR/raw_data_extraction_on_bam.py"
    "$bam_path"
    -o
    "$OUTPUT"
  )

  # Add interval arguments when any interval field was supplied.
  if [[ -n "$chrom" || -n "$start" || -n "$end" ]]; then
    cmd+=(
      -c "$chrom"
      -s "$start"
      -e "$end"
    )
  fi

  # Add a single-read filter when supplied.
  if [[ -n "$read_id" ]]; then
    cmd+=(
      -r "$read_id"
    )
  fi

  # Add the read IDs file when supplied.
  if [[ -n "$read_ids_file" ]]; then
    cmd+=(
      --read_ids_file "$read_ids_file"
    )
  fi

  # ---------------------------------------------------------------------------
  # Extract BrdU data
  # ---------------------------------------------------------------------------

  echo "[INFO] Running BrdU extraction..."
  "${cmd[@]}"

  # Require a nonempty BED output before continuing.
  if [[ ! -s "$OUTPUT" ]]; then
    echo "[ERROR] Output file $OUTPUT is empty or was not created. Exiting."
    deactivate
    exit 1
  fi

  echo "[INFO] BrdU extraction complete: $OUTPUT"
  echo "[INFO] BrdU data will remain in the BAM coordinate system."
  echo "[INFO] For W303 BAMs, this means BrdU data remains in W303 coordinates."

  # ---------------------------------------------------------------------------
  # Optional W303-to-sacCer coordinate liftOver
  # ---------------------------------------------------------------------------

  if [[ "$do_liftover" == "yes" ]]; then
    # Download the UCSC liftOver executable when it is not already installed.
    if [[ ! -x "$LIFTOVER_BIN" ]]; then
      echo "[INFO] liftOver not found at $LIFTOVER_BIN"
      echo "[INFO] Downloading UCSC liftOver into $TOOLS_DIR ..."

      wget \
        -O "$LIFTOVER_BIN" \
        "$LIFTOVER_URL"

      chmod +x "$LIFTOVER_BIN"

      echo "[INFO] liftOver installed successfully at: $LIFTOVER_BIN"
    else
      echo "[INFO] Using existing liftOver binary: $LIFTOVER_BIN"
    fi

    # Resolve the selected chain and annotation files.
    #
    # "|| true" prevents set -e from exiting before a custom error message can
    # be displayed.
    chain_path="$(resolve_liftover_chain "$target_strain" || true)"
    annotation_file="$(resolve_annotation_file "$target_strain" || true)"
    annotation_label="$target_strain"

    if [[ -z "$chain_path" ]]; then
      echo "[ERROR] No chain file found for target strain $target_strain in $CHAIN_DIR"
      echo "[ERROR] Expected W303TosacCer${target_strain#sacCer}.over.chain.gz"
      deactivate
      exit 1
    fi

    if [[ -z "$annotation_file" || ! -s "$annotation_file" ]]; then
      echo "[ERROR] Annotation file not found for target strain $target_strain"
      deactivate
      exit 1
    fi

    # Build output names that identify both the source and target assemblies.
    liftover_basename="W303_to_${target_strain}_${OUTPUT_BASENAME}"

    LIFTOVER_OUTPUT="$BED_DIR/${liftover_basename}.bed"

    LIFTOVER_UNMAPPED_OUTPUT="$BED_DIR/${liftover_basename}_unmapped.bed"

    RAINPLOT_MANIFEST="$RESULTS_DIR/${liftover_basename}_rainplots_manifest.txt"

    GENOMIC_FEATURE_PNG="$RESULTS_DIR/genomic_feature_${liftover_basename}.png"

    echo "[INFO] Running chain-based liftOver on BrdU BED using W303 GenBank chromosome names..."

    python "$UTILS_DIR/liftover_brdu_bed.py" \
      "$OUTPUT" \
      --chain "$chain_path" \
      --mapped "$LIFTOVER_OUTPUT" \
      --unmapped "$LIFTOVER_UNMAPPED_OUTPUT" \
      --liftover_bin "$LIFTOVER_BIN"

    # Require at least one successfully lifted interval.
    if [[ ! -s "$LIFTOVER_OUTPUT" ]]; then
      echo "[ERROR] LiftOver output file is empty or was not created: $LIFTOVER_OUTPUT"
      deactivate
      exit 1
    fi

    # Use lifted BrdU coordinates for downstream rainplot generation.
    bed_to_use="$LIFTOVER_OUTPUT"

    # Disable RFB overlays because the RFB BED remains in the original W303
    # coordinate system.
    use_rfb_overlay="no"

    generate_genomic_features="yes"

    # The default W303 G4 BED cannot be overlaid on a lifted assembly.
    g4_bed_file=""

    echo "[INFO] Chain-based LiftOver complete: $LIFTOVER_OUTPUT"
    echo "[INFO] Unmapped LiftOver intervals written to: $LIFTOVER_UNMAPPED_OUTPUT"
    echo "[INFO] The workflow will use the lifted BrdU BED file for plotting."
    echo "[INFO] RFB overlay is being skipped because the RFB BED remains in the original coordinate system."
    echo "[INFO] Genomic feature annotation will use the selected $target_strain annotation file."
  else
    echo "[INFO] UCSC liftOver skipped."
    echo "[INFO] The workflow will continue in the original BAM/W303 coordinate system."
  fi

  # ---------------------------------------------------------------------------
  # Optional S-phase replication-fork-barrier motif extraction
  # ---------------------------------------------------------------------------

  if [[ "$enable_rfb_workflow" == "yes" ]]; then
    # A read IDs file may contain reads from multiple chromosomes or regions.
    #
    # RFB extraction requires one defined interval, so it is disabled for an
    # unrestricted read-list run.
    if [[ "$uses_read_ids_file" == "yes" && -z "$chrom" && -z "$start" && -z "$end" ]]; then
      echo "[INFO] Skipping RFB motif extraction because TSV/list reads are not restricted to one region."

      use_rfb_overlay="no"
      filter_rfb_reads="no"
      request_rfb_overlay="no"
    else
      echo "[INFO] Running RFB motif extraction for S phase..."

      python "$SRC_DIR/rfb_seq_matcher.py" \
        "$bam_path" \
        -c "$chrom" \
        -s "$start" \
        -e "$end" \
        -o "$rfb_output"

      if [[ ! -s "$rfb_output" ]]; then
        echo "[WARN] RFB output file is empty or was not created: $rfb_output"
      else
        echo "[INFO] RFB extraction complete: $rfb_output"
      fi
    fi
  else
    echo "[INFO] Skipping RFB motif extraction because this run is not S phase."
  fi

  # ---------------------------------------------------------------------------
  # Generate rainplots and write a manifest of produced images
  # ---------------------------------------------------------------------------

  echo "[INFO] Generating rain plots..."

  # Remove stale manifest and region metadata files from an earlier run.
  rm -f \
    "$RAINPLOT_MANIFEST" \
    "$region_metadata_file"

  # Use the liftOver basename when coordinates were converted.
  local rainplot_filename_prefix="${liftover_basename:-$OUTPUT_BASENAME}"

  # RFB overlay is disabled unless every required condition is satisfied.
  local show_rfb_overlay="no"

  if [[ "$enable_rfb_workflow" == "yes" \
        && "$request_rfb_overlay" == "yes" \
        && "$use_rfb_overlay" == "yes" ]]; then

    show_rfb_overlay="yes"

  elif [[ "$enable_rfb_workflow" == "yes" \
          && "$request_rfb_overlay" == "yes" \
          && "$use_rfb_overlay" != "yes" ]]; then

    echo "[WARN] RFB overlay was requested, but it is unavailable after liftOver because the RFB BED remains in the original coordinate system."
  fi

  # ---------------------------------------------------------------------------
  # FUNCTION: run_rainplot_generation
  # ---------------------------------------------------------------------------
  #
  # Purpose:
  #   Assemble and execute the rainplot_generation.py command.
  #
  # Arguments:
  #   $1
  #       Output manifest path.
  #
  #   $2
  #       Filename prefix for generated images.
  #
  #   $3
  #       Whether to enable the RFB overlay.
  #
  #   $4
  #       Whether to filter for reads that contain RFB coordinates.
  #
  run_rainplot_generation() {
    local manifest_path="$1"
    local filename_prefix="$2"
    local enable_rfb_overlay="$3"
    local filter_rfb_reads_arg="$4"

    # Base rainplot-generation command.
    local rainplot_cmd=(
      python
      "$SRC_DIR/rainplot_generation.py"
      "$bed_to_use"
      -o
      "$RESULTS_DIR"
      --phase
      "$phase_label"
      --filename_prefix
      "$filename_prefix"
      --output_manifest
      "$manifest_path"
      --brdu_window_mode
      "$brdu_window_mode"
      --binary_threshold
      "$binary_threshold"
    )

    # Include region boundaries when the run uses one defined interval.
    if [[ -n "$start" && -n "$end" ]]; then
      rainplot_cmd+=(
        --region_start "$start"
        --region_end "$end"
      )
    fi

    # Restrict the number of plotted reads when requested.
    if [[ -n "$max_reads" ]]; then
      rainplot_cmd+=(
        --max_reads "$max_reads"
      )
    fi

    # Preserve the order from a ranked read IDs file when applying max_reads.
    if [[ -n "$read_ids_file" ]]; then
      rainplot_cmd+=(
        --preserve_read_order_for_max_reads
      )
    fi

    # Supply the RFB BED directory only when RFB overlaying or filtering is
    # required.
    if [[ "$enable_rfb_workflow" == "yes" \
          && ( "$enable_rfb_overlay" == "yes" \
               || "$filter_rfb_reads_arg" == "yes" ) ]]; then

      rainplot_cmd+=(
        --rfb_dir "$BED_DIR"
      )
    fi

    # Restrict output to reads containing RFB coordinates.
    if [[ "$filter_rfb_reads_arg" == "yes" ]]; then
      rainplot_cmd+=(
        --filter_reads_with_rfb
      )
    fi

    # Draw RFB coordinates on the rainplot.
    if [[ "$enable_rfb_overlay" == "yes" ]]; then
      rainplot_cmd+=(
        --show_rfb_overlay
      )
    fi

    # Execute the complete command array.
    "${rainplot_cmd[@]}"
  }

  # Generate the selected rainplots.
  run_rainplot_generation \
    "$RAINPLOT_MANIFEST" \
    "$rainplot_filename_prefix" \
    "$show_rfb_overlay" \
    "$filter_rfb_reads"

  echo "[INFO] Rain plots complete. Saved to: $RESULTS_DIR"

  # The manifest is required by combine_rainplot_images.py.
  if [[ ! -s "$RAINPLOT_MANIFEST" ]]; then
    echo "[ERROR] Rain plot manifest is empty or was not created: $RAINPLOT_MANIFEST"
    deactivate
    exit 1
  fi

  # ---------------------------------------------------------------------------
  # Generate genomic-feature panels and combine them with rainplots
  # ---------------------------------------------------------------------------

  if [[ "$generate_genomic_features" == "yes" ]]; then
    # Require a nonempty annotation file.
    if [[ ! -s "$annotation_file" ]]; then
      echo "[ERROR] Annotation file not found:"
      echo "$annotation_file"
      deactivate
      exit 1
    fi

    echo "[INFO] Generating genomic feature plot from selected annotation source..."
    echo "[INFO] Using annotation file: $annotation_file"

    # G4_BED_FILE is passed as an environment variable. It is blank after
    # liftOver because the W303 G4 BED is not in the lifted coordinate system.
    G4_BED_FILE="$g4_bed_file" \
      Rscript "$SRC_DIR/genomic_feature_plot.R" \
        "${chrom:-read_ids_file}" \
        "${start:-0}" \
        "${end:-1}" \
        "$GENOMIC_FEATURE_PNG" \
        "$annotation_file" \
        "$annotation_label"

    if [[ ! -s "$GENOMIC_FEATURE_PNG" ]]; then
      echo "[ERROR] Genomic feature plot generation failed."
      deactivate
      exit 1
    fi

    echo "[INFO] Genomic feature plot complete: $GENOMIC_FEATURE_PNG"

    # Combine every rainplot listed in the manifest with the genomic-feature
    # panel. --delete_inputs removes the separate component images after the
    # combined image has been created.
    echo "[INFO] Combining rain plots with genomic feature plot..."

    python "$SRC_DIR/combine_rainplot_images.py" \
      --manifest "$RAINPLOT_MANIFEST" \
      --annotation_png "$GENOMIC_FEATURE_PNG" \
      --output_dir "$RESULTS_DIR" \
      --delete_inputs

    echo "[INFO] Combined rain plots saved to: $RESULTS_DIR"
  else
    echo "[INFO] Combined genomic feature plots were skipped for this run."
  fi

  # The temporary manifest is no longer needed after image combination.
  rm -f "$RAINPLOT_MANIFEST"

  # Exit the workflow-specific Python virtual environment.
  deactivate

  echo "[INFO] Virtual environment deactivated."
  echo "[INFO] Workflow complete."
  echo "[INFO] Rainplot workflow finished: $(date)"
}


# ==============================================================================
# FUNCTION: submit_workflow
# ==============================================================================
#
# Purpose:
#   Collect and validate user input, then submit the workflow to SLURM.
#
# This function:
#   - Lists available BAM files.
#   - Accepts interval-based or read-list-based selection.
#   - Prompts for binary or mean BrdU processing.
#   - Prompts for the maximum number of reads.
#   - Configures Mitosis or S Phase behavior.
#   - Configures RFB behavior for S Phase.
#   - Optionally configures UCSC liftOver.
#   - Submits this script back to SLURM using --run-job.
#
# Positional arguments:
#   $1
#       Optional BAM filename.
#
#   $2
#       Optional chromosome.
#
#   $3
#       Optional start coordinate.
#
#   $4
#       Optional end coordinate.
#
#   $5
#       Optional single read ID.
#
#   $6
#       Optional output BED filename.
#
#   $7
#       Optional read IDs TSV/list file.
#
submit_workflow() {
  # ---------------------------------------------------------------------------
  # Read optional positional arguments
  # ---------------------------------------------------------------------------

  local bam_input="${1:-}"
  local chrom="${2:-}"
  local start="${3:-}"
  local end="${4:-}"
  local read_id="${5:-}"
  local output_input="${6:-}"
  local read_ids_file_input="${7:-}"

  # ---------------------------------------------------------------------------
  # Declare remaining submission variables
  # ---------------------------------------------------------------------------

  local bam_path
  local read_ids_file=""
  local output_basename
  local default_output
  local phase_input
  local phase_label
  local brdu_window_mode_input
  local brdu_window_mode
  local binary_threshold_input
  local binary_threshold="0.5"
  local max_reads_input
  local max_reads
  local s_phase_plot_choice

  # Default RFB configuration.
  local rfb_plot_mode="none"
  local filter_rfb_reads="no"
  local request_rfb_overlay="yes"
  local enable_rfb_workflow="no"

  # liftOver configuration.
  local do_liftover
  local target_strain="W303"

  # SLURM job identifiers.
  local job_id
  local log_job_id

  # Tracks whether a read IDs file controls read selection.
  local uses_read_ids_file="no"

  # Create the expected workflow directory structure.
  mkdir -p \
    "$BAM_DIR" \
    "$BED_DIR" \
    "$CHAIN_DIR" \
    "$NCBI_ROOT" \
    "$TOOLS_DIR" \
    "$RESULTS_DIR" \
    "$LOG_DIR"

  # ---------------------------------------------------------------------------
  # Select and validate the BAM file
  # ---------------------------------------------------------------------------

  while true; do
    if [[ -z "$bam_input" ]]; then
      echo "[INFO] Available BAM files in $BAM_DIR:"

      find "$BAM_DIR" \
        -maxdepth 1 \
        -type f \
        -name "*.bam" \
        -printf "  %f\n" |
        sort

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

  # ---------------------------------------------------------------------------
  # Optionally select a TSV/list containing read IDs
  # ---------------------------------------------------------------------------

  while true; do
    if [[ -z "$read_ids_file_input" ]]; then
      read -r -p "Enter read IDs TSV/list file [blank for none]: " read_ids_file_input
    fi

    # A single read ID and a read IDs file cannot both control selection.
    if [[ -n "$read_id" && -n "$read_ids_file_input" ]]; then
      warn_unexpected_input "$read_ids_file_input" "blank because a single read ID was already supplied"
      read_ids_file_input=""
      continue
    fi

    if [[ -z "$read_ids_file_input" ]]; then
      break
    fi

    if [[ -f "$read_ids_file_input" ]]; then
      read_ids_file="$(absolute_existing_file "$read_ids_file_input")"
      uses_read_ids_file="yes"
      break
    fi

    warn_unexpected_input "$read_ids_file_input" "an existing TSV/list file path, or blank for none"
    read_ids_file_input=""
  done

  # ---------------------------------------------------------------------------
  # Select binary or mean BrdU T-window processing
  # ---------------------------------------------------------------------------

  while true; do
    echo
    echo "Choose BrdU T-window calculation mode:"
    echo "  1) binary - threshold probabilities, then average calls in each 100-T window"
    echo "  2) mean   - average raw probabilities in each 100-T window"

    read -r -p "Enter binary or mean [binary]: " brdu_window_mode_input

    if brdu_window_mode="$(
      normalize_brdu_window_mode "$brdu_window_mode_input"
    )"; then
      break
    fi

    warn_unexpected_input "$brdu_window_mode_input" "binary or mean, for example binary"
  done

  # Binary mode requires a BrdU probability threshold.
  if [[ "$brdu_window_mode" == "binary" ]]; then
    while true; do
      read -r -p "Enter binary BrdU probability threshold [0.5]: " binary_threshold_input

      binary_threshold_input="${binary_threshold_input:-0.5}"

      if binary_threshold="$(
        normalize_probability_threshold "$binary_threshold_input"
      )"; then
        break
      fi

      warn_unexpected_input "$binary_threshold_input" "a number from 0 to 1, for example 0.5"
    done
  fi

  # ---------------------------------------------------------------------------
  # Collect interval coordinates unless a read IDs file controls selection
  # ---------------------------------------------------------------------------

  if [[ "$uses_read_ids_file" == "no" ]]; then
    while [[ -z "$chrom" ]]; do
      read -r -p "Enter chromosome GenBank ID: " chrom
      if [[ -z "$chrom" ]]; then
        warn_unexpected_input "$chrom" "a chromosome ID, for example CM007964.1 or chr1"
      fi
    done

    while [[ -z "$start" || ! "$start" =~ ^[0-9]+$ ]]; do
      if [[ -n "$start" ]]; then
        warn_unexpected_input "$start" "an integer start coordinate, for example 10000"
      fi
      read -r -p "Enter start coordinate: " start
    done

    while [[ -z "$end" || ! "$end" =~ ^[0-9]+$ ]]; do
      if [[ -n "$end" ]]; then
        warn_unexpected_input "$end" "an integer end coordinate greater than the start, for example 50000"
      fi
      read -r -p "Enter end coordinate: " end
    done
  else
    # A read IDs file may contain reads from multiple intervals. Clear region
    # fields so downstream scripts process all coordinates associated with the
    # selected reads.
    chrom=""
    start=""
    end=""
  fi

  # Ask for an optional single-read filter only when neither selection method
  # has already supplied one.
  if [[ -z "$read_id" && -z "$read_ids_file" ]]; then
    read -r -p "Enter read ID filter [blank for all reads in region]: " read_id
  fi

  # ---------------------------------------------------------------------------
  # Select how many reads will be plotted
  # ---------------------------------------------------------------------------

  while true; do
    if [[ "$uses_read_ids_file" == "yes" ]]; then
      read -r -p "How many reads should be plotted? [blank for all reads in TSV/list]: " max_reads_input
    else
      read -r -p "How many reads should be plotted? [blank for all reads in region]: " max_reads_input
    fi

    if max_reads="$(normalize_max_reads "$max_reads_input")"; then
      break
    fi

    warn_unexpected_input "$max_reads_input" "a positive integer, for example 25, or blank for all reads"
  done

  # ---------------------------------------------------------------------------
  # Validate interval coordinates
  # ---------------------------------------------------------------------------

  if [[ "$uses_read_ids_file" == "no" ]]; then
    while (( start >= end )); do
      warn_unexpected_input "$start-$end" "a start coordinate less than the end coordinate, for example 10000 then 50000"
      read -r -p "Enter start coordinate: " start
      while [[ -z "$start" || ! "$start" =~ ^[0-9]+$ ]]; do
        warn_unexpected_input "$start" "an integer start coordinate, for example 10000"
        read -r -p "Enter start coordinate: " start
      done
      read -r -p "Enter end coordinate: " end
      while [[ -z "$end" || ! "$end" =~ ^[0-9]+$ ]]; do
        warn_unexpected_input "$end" "an integer end coordinate greater than the start, for example 50000"
        read -r -p "Enter end coordinate: " end
      done
    done
  fi

  # ---------------------------------------------------------------------------
  # Build or collect the output BED basename
  # ---------------------------------------------------------------------------

  if [[ "$uses_read_ids_file" == "yes" ]]; then
    default_output="$(
      append_brdu_mode_to_basename \
        "$(default_read_ids_output_basename "$bam_path" "$max_reads")" \
        "$brdu_window_mode"
    ).bed"
  else
    default_output="$(
      append_brdu_mode_to_basename \
        "$(default_output_basename "$bam_path" "$chrom" "$start" "$end")" \
        "$brdu_window_mode"
    ).bed"
  fi

  if [[ -z "$output_input" ]]; then
    read -r -p "Output BED filename [$default_output]: " output_input
  fi

  # Use the generated default when the user presses Enter.
  output_input="${output_input:-$default_output}"

  # Remove an optional .bed extension and ensure the BrdU mode is present.
  output_basename="$(
    append_brdu_mode_to_basename \
      "$(basename "$output_input" .bed)" \
      "$brdu_window_mode"
  )"

  # ---------------------------------------------------------------------------
  # Select the cell-cycle phase
  # ---------------------------------------------------------------------------

  while true; do
    read -r -p "Which phase is this run for? [M/S]: " phase_input

    if phase_label="$(normalize_phase_label "$phase_input")"; then
      break
    fi

    warn_unexpected_input "$phase_input" "M or S, for example S"
  done

  # ---------------------------------------------------------------------------
  # Configure RFB behavior
  # ---------------------------------------------------------------------------
  #
  # Mitosis:
  #   RFB processing is disabled.
  #
  # S Phase:
  #   The user can choose:
  #
  #     a. RFB-coordinate reads only.
  #     b. Rainplots without RFB processing.
  #     c. Standard rainplots with an RFB coordinate overlay.
  #
  if [[ "$phase_label" == "Mitosis" ]]; then
    rfb_plot_mode="none"
    filter_rfb_reads="no"
    request_rfb_overlay="no"
    enable_rfb_workflow="no"
  else
    enable_rfb_workflow="yes"

    while true; do
      read -r -p "For S Phase, choose rain plot mode: [a] RFB coords only, [b] without RFB, [c] with and without RFB coords: " s_phase_plot_choice

      case "$s_phase_plot_choice" in
        [Aa])
          rfb_plot_mode="rfb_only"
          filter_rfb_reads="yes"
          request_rfb_overlay="yes"
          break
          ;;
        [Bb])
          rfb_plot_mode="without_rfb"
          filter_rfb_reads="no"
          request_rfb_overlay="no"
          break
          ;;
        [Cc])
          rfb_plot_mode="mixed"
          filter_rfb_reads="no"
          request_rfb_overlay="yes"
          break
          ;;
        *)
          warn_unexpected_input "$s_phase_plot_choice" "a, b, or c"
          ;;
      esac
    done
  fi

  # ---------------------------------------------------------------------------
  # Configure optional UCSC liftOver
  # ---------------------------------------------------------------------------

  while true; do
    read -r -p "Would you like to do a UCSC liftOver on the BrdU BED? [y/n]: " do_liftover

    case "$do_liftover" in
      [Yy])
        do_liftover="yes"

        while true; do
          read -r -p "Which target yeast strain do you want to lift over to? [sacCer1/sacCer2/sacCer3]: " target_strain

          case "$target_strain" in
            sacCer1|sacCer2|sacCer3)
              # Validate that the required chain file exists.
              if ! resolve_liftover_chain "$target_strain" >/dev/null; then
                warn_unexpected_input "$target_strain" "a target with chain file W303TosacCer${target_strain#sacCer}.over.chain.gz in $CHAIN_DIR"
                continue
              fi

              # Validate the target annotation file.
              if [[ ! -s "$(resolve_annotation_file "$target_strain" || true)" ]]; then
                warn_unexpected_input "$target_strain" "a target strain with an available annotation file"
                continue
              fi

              break
              ;;
            *)
              warn_unexpected_input "$target_strain" "sacCer1, sacCer2, or sacCer3"
              ;;
          esac
        done

        break
        ;;
      [Nn])
        do_liftover="no"
        break
        ;;
      *)
        warn_unexpected_input "$do_liftover" "y or n"
        ;;
    esac
  done

  # ---------------------------------------------------------------------------
  # Display the final workflow configuration
  # ---------------------------------------------------------------------------

  echo
  echo "[INFO] Workflow root: $WORKFLOW_ROOT"
  echo "[INFO] BAM: $bam_path"

  if [[ "$uses_read_ids_file" == "yes" ]]; then
    echo "[INFO] Region: all coordinates for reads listed in read IDs file"
  else
    echo "[INFO] Region: $chrom:$start-$end"
  fi

  if [[ "$uses_read_ids_file" == "yes" ]]; then
    echo "[INFO] Read ID filter: exact reads from read IDs file"
  else
    echo "[INFO] Read ID filter: ${read_id:-all reads in region}"
  fi

  echo "[INFO] Read IDs file: ${read_ids_file:-none}"
  echo "[INFO] Rainplot read count: ${max_reads:-all selected reads}"
  echo "[INFO] Output BED: $BED_DIR/${output_basename}.bed"
  echo "[INFO] BrdU T-window mode: $brdu_window_mode"

  if [[ "$brdu_window_mode" == "binary" ]]; then
    echo "[INFO] Binary BrdU probability threshold: $binary_threshold"
  fi

  echo "[INFO] Phase label: $phase_label"
  echo "[INFO] RFB plot mode: $rfb_plot_mode"
  echo "[INFO] liftOver: $do_liftover"

  if [[ "$do_liftover" == "yes" ]]; then
    echo "[INFO] liftOver target: $target_strain"
  fi

  echo
  echo "Submitting rainplot workflow to SLURM..."

  # ---------------------------------------------------------------------------
  # Submit this script to SLURM in --run-job mode
  # ---------------------------------------------------------------------------
  #
  # --parsable
  #   Return the submitted job ID in a machine-readable form.
  #
  # --chdir
  #   Run the compute job from the workflow root.
  #
  # --export
  #   Export the current environment and the resolved workflow root.
  #
  job_id="$(sbatch --parsable \
    --job-name="rainplot_${output_basename}" \
    --chdir="$WORKFLOW_ROOT" \
    --output="$LOG_DIR/${output_basename}.%j.slurm.log" \
    --error="$LOG_DIR/${output_basename}.%j.slurm.err" \
    --export=ALL,RAINPLOT_WORKFLOW_ROOT="$WORKFLOW_ROOT" \
    "$SCRIPT_PATH" \
    --run-job \
    --workflow-root "$WORKFLOW_ROOT" \
    --bam "$bam_path" \
    --chrom "$chrom" \
    --start "$start" \
    --end "$end" \
    --read-id "$read_id" \
    --read-ids-file "$read_ids_file" \
    --output-basename "$output_basename" \
    --phase-label "$phase_label" \
    --brdu-window-mode "$brdu_window_mode" \
    --binary-threshold "$binary_threshold" \
    --max-reads "$max_reads" \
    --rfb-plot-mode "$rfb_plot_mode" \
    --filter-rfb-reads "$filter_rfb_reads" \
    --request-rfb-overlay "$request_rfb_overlay" \
    --enable-rfb-workflow "$enable_rfb_workflow" \
    --do-liftover "$do_liftover" \
    --target-strain "$target_strain")"

  # Some SLURM environments append cluster information after a semicolon.
  # Retain only the leading job-ID portion for constructing log paths.
  log_job_id="${job_id%%;*}"

  # ---------------------------------------------------------------------------
  # Report the submitted job and expected outputs
  # ---------------------------------------------------------------------------

  echo "[INFO] Submitted SLURM job: $job_id"
  echo "[INFO] Expected BED:         $BED_DIR/${output_basename}.bed"
  echo "[INFO] Expected results dir: $RESULTS_DIR"
  echo "[INFO] SLURM log:            $LOG_DIR/${output_basename}.${log_job_id}.slurm.log"
  echo "[INFO] SLURM err:            $LOG_DIR/${output_basename}.${log_job_id}.slurm.err"
  echo "[INFO] Workflow log:         $LOG_DIR/${output_basename}.${log_job_id}.workflow.log"
}


# ==============================================================================
# MAIN SCRIPT ROUTING
# ==============================================================================
#
# -h or --help:
#   Print usage information and exit.
#
# --run-job:
#   Execute the compute stage. This is normally passed internally by the
#   submitted SLURM job.
#
# Any other invocation:
#   Enter interactive submission mode.
#

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

if [[ "${1:-}" == "--run-job" ]]; then
  run_job "$@"
else
  submit_workflow "$@"
fi
