#!/bin/bash

# Stop the script immediately if:
# - any command exits with an error,
# - an unset variable is referenced,
# - any command in a pipeline fails.
set -euo pipefail

# Display the expected command-line arguments for this script.
usage() {
    echo "Usage: bash src/nanopore_sequence_workflow/dorado_basecall.sh --pod5 POD5_FILE_OR_DIR --output-dir DIR --log-file LOG --model MODEL --device DEVICE --min-qscore N [--output-format fastq|bam] [--reference FASTA] [--emit-moves yes|no] [--mm2-opts OPTS] [--kit-name KIT] [--barcode-mode none|demux] [--modified-bases MODEL_OR_CODE]"
}

# Initialize all input variables.
#
# Required inputs:
# - POD5
# - OUTPUT_DIR
# - LOG_FILE
#
# Optional inputs receive default values later in the script.
POD5=""
OUTPUT_DIR=""
LOG_FILE=""
DORADO_MODEL=""
DEVICE=""
MIN_QSCORE=""
OUTPUT_FORMAT="fastq"
REFERENCE=""
EMIT_MOVES="yes"
MM2_OPTS=""
KIT_NAME=""
BARCODE_MODE="none"
MODIFIED_BASES=""

# Process all command-line arguments supplied to the script.
while [[ $# -gt 0 ]]; do
    case "$1" in
        # Path to the input POD5 raw-signal file.
        --pod5)
            POD5="$2"
            shift 2
            ;;

        # Directory where the FASTQ output will be written.
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;

        # Path to the log file that will record run details and QC results.
        --log-file)
            LOG_FILE="$2"
            shift 2
            ;;

        # Dorado basecalling model.
        #
        # Both --model and --dorado-model are accepted as equivalent options.
        --model|--dorado-model)
            DORADO_MODEL="$2"
            shift 2
            ;;

        # Compute device used by Dorado, such as cuda:0 or cpu.
        --device)
            DEVICE="$2"
            shift 2
            ;;

        # Minimum read quality score required for reads to be emitted.
        --min-qscore)
            MIN_QSCORE="$2"
            shift 2
            ;;

        # Select whether Dorado writes FASTQ or aligned BAM output.
        --output-format)
            OUTPUT_FORMAT="$2"
            shift 2
            ;;

        # Reference FASTA used when Dorado produces aligned BAM output.
        --reference|--ref)
            REFERENCE="$2"
            shift 2
            ;;

        # Controls whether Dorado writes move-table information to BAM.
        --emit-moves)
            EMIT_MOVES="$2"
            shift 2
            ;;

        # Optional minimap2 option string passed through Dorado alignment.
        --mm2-opts)
            MM2_OPTS="$2"
            shift 2
            ;;

        # Sequencing-kit name used for barcode classification.
        --kit-name)
            KIT_NAME="$2"
            shift 2
            ;;

        # Controls whether barcode demultiplexing is performed.
        #
        # Accepted values:
        # - none
        # - demux
        --barcode-mode)
            BARCODE_MODE="$2"
            shift 2
            ;;

        # Optional Dorado modified-base model or modification code.
        --modified-bases)
            MODIFIED_BASES="$2"
            shift 2
            ;;

        # Display usage information and exit successfully.
        -h|--help)
            usage
            exit 0
            ;;

        # Stop when an unsupported command-line argument is encountered.
        *)
            echo "[ERROR] Unknown argument: $1"
            usage
            exit 1
            ;;
    esac
done

# Apply default settings when optional values were not supplied.
#
# Default model:
#   sup
#
# Default device:
#   the first visible CUDA GPU
#
# Default minimum quality score:
#   Q6
#
# Default barcode mode:
#   no demultiplexing
DORADO_MODEL="${DORADO_MODEL:-sup}"
DEVICE="${DEVICE:-cuda:0}"
MIN_QSCORE="${MIN_QSCORE:-6}"
OUTPUT_FORMAT="${OUTPUT_FORMAT:-fastq}"
EMIT_MOVES="${EMIT_MOVES:-yes}"
BARCODE_MODE="${BARCODE_MODE:-none}"

# Confirm that the required arguments were provided.
if [[ -z "$POD5" || -z "$OUTPUT_DIR" || -z "$LOG_FILE" ]]; then
    echo "[ERROR] Missing required Dorado argument."
    usage
    exit 1
fi

# Confirm that the input POD5 file or directory exists.
if [[ ! -f "$POD5" && ! -d "$POD5" ]]; then
    echo "[ERROR] POD5 input not found: $POD5"
    exit 1
fi

case "$OUTPUT_FORMAT" in
    fastq|bam) ;;
    *)
        echo "[ERROR] Invalid output format: $OUTPUT_FORMAT"
        echo "[ERROR] Expected fastq or bam."
        exit 1
        ;;
esac

case "$EMIT_MOVES" in
    yes|no) ;;
    *)
        echo "[ERROR] Invalid --emit-moves value: $EMIT_MOVES"
        echo "[ERROR] Expected yes or no."
        exit 1
        ;;
esac

if [[ "$OUTPUT_FORMAT" == "bam" && -z "$REFERENCE" ]]; then
    echo "[ERROR] --reference is required when --output-format bam is used."
    exit 1
fi

if [[ -n "$REFERENCE" && ! -f "$REFERENCE" ]]; then
    echo "[ERROR] Reference FASTA not found: $REFERENCE"
    exit 1
fi

if [[ -f "$POD5" && "$POD5" != *.pod5 ]]; then
    echo "[ERROR] POD5 input file must end in .pod5: $POD5"
    exit 1
fi

if [[ -d "$POD5" ]] && ! find "$POD5" -maxdepth 1 -type f -name "*.pod5" -print -quit | grep -q .; then
    echo "[ERROR] POD5 directory contains no .pod5 files: $POD5"
    exit 1
fi

# Validate the barcode mode.
#
# Only no demultiplexing or demultiplexing are accepted.
case "$BARCODE_MODE" in
    none|demux) ;;
    *)
        echo "[ERROR] Invalid barcode mode: $BARCODE_MODE"
        echo "[ERROR] Expected none or demux."
        exit 1
        ;;
esac

# A sequencing-kit name is required for Dorado to identify
# which barcode arrangement should be used during demultiplexing.
if [[ "$BARCODE_MODE" == "demux" && -z "$KIT_NAME" ]]; then
    echo "[ERROR] --kit-name is required when --barcode-mode demux is used."
    exit 1
fi

# Create the FASTQ output directory and the parent directory
# containing the log file if they do not already exist.
mkdir -p "$OUTPUT_DIR" "$(dirname "$LOG_FILE")"

# Derive the sample name from the input POD5 file or directory name.
#
# Example:
# sample.pod5 becomes sample
# sample_pod5_dir becomes sample_pod5_dir
POD5_BASENAME="$(basename "$POD5")"
POD5_PREFIX="${POD5_BASENAME%.pod5}"
JOB_ID="${SLURM_JOB_ID:-manual}"

# Construct the main FASTQ output path.
#
# Include the SLURM job ID so repeated runs on the same POD5 input
# do not overwrite one another when parameters differ.
FASTQ_FILE="$OUTPUT_DIR/${POD5_PREFIX}_${JOB_ID}.fastq"
BAM_FILE="$OUTPUT_DIR/${POD5_PREFIX}_${JOB_ID}.bam"
OUTPUT_FILE="$FASTQ_FILE"
OUTPUT_LABEL="FASTQ"
if [[ "$OUTPUT_FORMAT" == "bam" ]]; then
    OUTPUT_FILE="$BAM_FILE"
    OUTPUT_LABEL="BAM"
fi

# Construct the directory used for barcode-separated FASTQ files.
BARCODE_DIR="$OUTPUT_DIR/${POD5_PREFIX}_${JOB_ID}_barcodes"

# Write the basecalling configuration to the log file.
#
# The single > operator creates a new log file or replaces an
# existing file with the same name.
{
    echo "========================================="
    echo "  DORADO BASECALL"
    echo "========================================="
    echo "Started:        $(date)"
    echo "POD5:           $POD5"
    echo "Output format:  $OUTPUT_FORMAT"
    echo "Output $OUTPUT_LABEL:   $OUTPUT_FILE"
    echo "Reference:      ${REFERENCE:-none}"
    echo "Model:          $DORADO_MODEL"
    echo "Device:         $DEVICE"
    echo "Min q-score:    $MIN_QSCORE"
    echo "Emit moves:     $EMIT_MOVES"
    echo "MM2 opts:       ${MM2_OPTS:-none}"
    echo "Barcode mode:   $BARCODE_MODE"
    echo "Kit name:       ${KIT_NAME:-none}"
    echo "Modified bases: ${MODIFIED_BASES:-none}"
    echo "SLURM job ID:   ${SLURM_JOB_ID:-not_set}"
    echo "========================================="
} > "$LOG_FILE"

# Record information about the assigned NVIDIA GPU when nvidia-smi
# is available in the current environment.
#
# A failure from nvidia-smi does not stop the workflow because of || true.
if command -v nvidia-smi >/dev/null 2>&1; then
    {
        echo ""
        echo "========================================="
        echo "  GPU STATUS"
        echo "========================================="
        echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-not_set}"
        nvidia-smi
    } >> "$LOG_FILE" 2>&1 || true
fi

# Build the main Dorado basecaller command as a Bash array.
#
# basecaller:
#   Runs Dorado basecalling.
#
# DORADO_MODEL:
#   Specifies the neural-network basecalling model.
#
# POD5:
#   Supplies the raw nanopore signal.
#
# --device:
#   Selects the GPU or CPU device.
#
# --min-qscore:
#   Filters out reads below the selected quality threshold.
#
DORADO_ARGS=(
    basecaller
    "$DORADO_MODEL"
    "$POD5"
    --device "$DEVICE"
    --min-qscore "$MIN_QSCORE"
)

if [[ "$OUTPUT_FORMAT" == "fastq" ]]; then
    # --emit-fastq produces FASTQ-formatted output for the Minimap2 path.
    DORADO_ARGS+=(--emit-fastq)
else
    # --reference makes Dorado basecall and align in one step, producing BAM.
    DORADO_ARGS+=(--reference "$REFERENCE")

    if [[ "$EMIT_MOVES" == "yes" ]]; then
        DORADO_ARGS+=(--emit-moves)
    fi

    if [[ -n "$MM2_OPTS" ]]; then
        DORADO_ARGS+=(--mm2-opts "$MM2_OPTS")
    fi
fi

# Add modified-base detection to the Dorado command when a model
# or modification code was supplied.
if [[ -n "$MODIFIED_BASES" ]]; then
    DORADO_ARGS+=(--modified-bases "$MODIFIED_BASES")
fi

# When demultiplexing is requested, supply the sequencing-kit name.
#
# --no-trim tells Dorado not to remove barcode or adapter sequence
# during the initial basecalling step.
if [[ "$BARCODE_MODE" == "demux" ]]; then
    DORADO_ARGS+=(--kit-name "$KIT_NAME" --no-trim)
fi

# Record the exact Dorado command in the log file.
#
# %q writes each argument in a shell-escaped form so that the command
# can be copied and reproduced more reliably.
{
    echo ""
    echo "[INFO] Running Dorado basecaller."
    printf '[INFO] Command: dorado'
    printf ' %q' "${DORADO_ARGS[@]}"
    echo " > $OUTPUT_FILE"
} >> "$LOG_FILE"

# Run Dorado basecalling.
#
# Standard output is written to the selected output file.
# Standard error, including Dorado progress information, is appended
# to the log file.
dorado "${DORADO_ARGS[@]}" > "$OUTPUT_FILE" 2>> "$LOG_FILE"

# Confirm that Dorado produced a nonempty output file.
if [[ ! -s "$OUTPUT_FILE" ]]; then
    echo "[ERROR] $OUTPUT_LABEL file is empty or was not created: $OUTPUT_FILE" | tee -a "$LOG_FILE"
    exit 1
fi

# Perform barcode demultiplexing when requested.
if [[ "$BARCODE_MODE" == "demux" && "$OUTPUT_FORMAT" == "fastq" ]]; then
    # Create the directory that will contain barcode-specific FASTQ files.
    mkdir -p "$BARCODE_DIR"

    # Record the start of the demultiplexing stage.
    {
        echo ""
        echo "[INFO] Running Dorado demux."
        echo "[INFO] Barcode output directory: $BARCODE_DIR"
    } >> "$LOG_FILE"

    # Separate the basecalled reads by barcode.
    #
    # --kit-name:
    #   Tells Dorado which barcode kit was used.
    #
    # --output-dir:
    #   Specifies where the barcode FASTQ files are written.
    #
    # --emit-fastq:
    #   Produces FASTQ output for each barcode classification.
    dorado demux "$FASTQ_FILE" \
        --kit-name "$KIT_NAME" \
        --output-dir "$BARCODE_DIR" \
        --emit-fastq \
        >> "$LOG_FILE" 2>&1

    # Rename Dorado's barcode FASTQ outputs so every filename begins
    # with the original POD5 sample prefix and job ID.
    for file in "$BARCODE_DIR"/*.fastq; do
        # Skip the loop when no FASTQ files matched the wildcard.
        [[ -e "$file" ]] || continue

        basename_file="$(basename "$file")"

        # Rename classified barcode files.
        #
        # Example:
        # barcode01.fastq becomes sample_jobid_barcode01.fastq
        if [[ "$basename_file" =~ (barcode[0-9]+) ]]; then
            mv "$file" "$BARCODE_DIR/${POD5_PREFIX}_${JOB_ID}_${BASH_REMATCH[1]}.fastq"

        # Rename reads that could not be assigned to a barcode.
        elif [[ "$basename_file" == *unclassified* ]]; then
            mv "$file" "$BARCODE_DIR/${POD5_PREFIX}_${JOB_ID}_unclassified.fastq"
        fi
    done
fi

if [[ "$OUTPUT_FORMAT" == "fastq" ]]; then
    # Calculate and write post-basecalling quality-control information.
    {
    echo ""
    echo "========================================="
    echo "  BASECALL QC"
    echo "========================================="

    # Count FASTQ records.
    #
    # Each complete FASTQ record contains four lines, so the total
    # line count is divided by four.
    FASTQ_RECORDS=$(awk 'END {print NR / 4}' "$FASTQ_FILE")

    # Reads present in the output FASTQ are treated as passing reads
    # because Dorado applied the minimum quality-score filter.
    PASS_COUNT="$FASTQ_RECORDS"

    # Attempt to extract Dorado's failed-read count from the log.
    #
    # tail -1 uses the final reported failed count when the value
    # appears multiple times in the progress output.
    FAILED_COUNT=$(grep -oE 'failed=[0-9]+' "$LOG_FILE" | tail -1 | cut -d= -f2 || true)

    # Attempt to extract the total number of basecalled reads from the log.
    BASECALLED_COUNT=$(grep -oE 'basecalled=[0-9]+' "$LOG_FILE" | tail -1 | cut -d= -f2 || true)

    # Calculate read-count and read-length statistics from the FASTQ.
    #
    # FASTQ sequence lines occur where the line number modulo four equals two.
    #
    # The calculation reports:
    # - total reads,
    # - mean read length,
    # - minimum read length,
    # - maximum read length.
    READ_LENGTH_STATS=$(awk 'NR % 4 == 2 {n++; len=length($0); sum+=len; if (min == "" || len < min) min=len; if (len > max) max=len} END {if (n > 0) printf "Reads: %d\nMean read length: %.1f bp\nMin read length: %d bp\nMax read length: %d bp\n", n, sum/n, min, max; else print "Reads: 0"}' "$FASTQ_FILE")

    # Count how many FASTQ reads are at least 1,000 bases long.
    READS_OVER_1000=$(awk 'NR % 4 == 2 && length($0) >= 1000 {n++} END {print n + 0}' "$FASTQ_FILE")

    # Report passing, failed, and total basecalled read counts.
    echo "Passed reads in FASTQ: ${PASS_COUNT:-N/A}"
    echo "Failed reads from log: ${FAILED_COUNT:-N/A}"
    echo "Total basecalled:      ${BASECALLED_COUNT:-N/A}"

    # Calculate the percentage of total basecalled reads that passed
    # the minimum quality threshold when a total count was found.
    if [[ -n "${BASECALLED_COUNT:-}" && "$BASECALLED_COUNT" -gt 0 ]]; then
        PASS_PCT=$(awk "BEGIN {printf \"%.1f\", ($PASS_COUNT / $BASECALLED_COUNT) * 100}")
        echo "Pass rate:             ${PASS_PCT}%"
    fi

    # Report the read-length statistics.
    echo ""
    echo "FASTQ read length summary:"
    echo "$READ_LENGTH_STATS"
    echo "Reads >= 1000 bp: $READS_OVER_1000"

    # When demultiplexing was performed, report the number of reads
    # in each barcode-specific FASTQ file.
    if [[ "$BARCODE_MODE" == "demux" ]]; then
        echo ""
        echo "Barcode FASTQ counts:"

        for file in "$BARCODE_DIR"/*.fastq; do
            # Skip the loop if no barcode FASTQ files exist.
            [[ -e "$file" ]] || continue

            # Count FASTQ headers by counting lines beginning with @.
            read_count=$(grep -c "^@" "$file" || true)
            echo "$(basename "$file"): $read_count reads"
        done
    fi

    # Record the completion time.
    echo "Completed: $(date)"
    echo "========================================="
    } >> "$LOG_FILE"
else
    {
        echo ""
        echo "========================================="
        echo "  BASECALL QC"
        echo "========================================="
        echo "Dorado produced aligned BAM output."
        echo "Move table emitted: $EMIT_MOVES"
        echo "Completed: $(date)"
        echo "========================================="
    } >> "$LOG_FILE"
fi

# Print the primary output path to standard output.
#
# This allows a parent workflow or calling script to capture the
# generated filename.
echo "$OUTPUT_FILE"
