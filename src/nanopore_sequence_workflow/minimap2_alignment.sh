#!/bin/bash

# Stop the script immediately if:
# - any command exits with an error,
# - an unset variable is referenced,
# - any command within a pipeline fails.
set -euo pipefail

# Display the expected command-line arguments for the alignment script.
usage() {
    echo "Usage: bash src/nanopore_sequence_workflow/minimap2_alignment.sh --fastq FASTQ --reference FASTA --output-dir DIR --log-file LOG --preset PRESET --threads N --secondary yes|no --sort yes|no --index yes|no --min-mapq N --min-read-length N --primary-only yes|no"
}

# Initialize the required input paths.
FASTQ=""
REFERENCE=""
OUTPUT_DIR=""
LOG_FILE=""

# Initialize Minimap2, samtools, and filtering settings.
#
# These starting values also act as the normal defaults unless they
# are replaced by command-line arguments.
PRESET="map-ont"
THREADS="8"
SECONDARY="no"
SORT_BAM="yes"
INDEX_BAM="yes"
MIN_MAPQ="20"
MIN_READ_LENGTH="1000"
PRIMARY_ONLY="yes"

# Process all command-line arguments supplied to the script.
while [[ $# -gt 0 ]]; do
    case "$1" in
        # Path to the FASTQ file containing basecalled nanopore reads.
        --fastq)
            FASTQ="$2"
            shift 2
            ;;

        # Path to the reference genome FASTA.
        #
        # Both --reference and --ref are accepted.
        --reference|--ref)
            REFERENCE="$2"
            shift 2
            ;;

        # Directory where BAM files and the QC report will be written.
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;

        # Path to the log file for alignment messages and QC output.
        --log-file)
            LOG_FILE="$2"
            shift 2
            ;;

        # Minimap2 alignment preset.
        #
        # map-ont is intended for Oxford Nanopore reads.
        --preset)
            PRESET="$2"
            shift 2
            ;;

        # Number of CPU threads used by Minimap2 and samtools sort.
        --threads)
            THREADS="$2"
            shift 2
            ;;

        # Controls whether Minimap2 reports secondary alignments.
        #
        # Accepted values:
        # - yes
        # - no
        --secondary)
            SECONDARY="$2"
            shift 2
            ;;

        # Controls whether the BAM is coordinate-sorted.
        --sort)
            SORT_BAM="$2"
            shift 2
            ;;

        # Controls whether a BAM index is created.
        --index)
            INDEX_BAM="$2"
            shift 2
            ;;

        # Minimum mapping-quality threshold used during QC.
        --min-mapq)
            MIN_MAPQ="$2"
            shift 2
            ;;

        # Minimum read-length threshold used during QC.
        --min-read-length)
            MIN_READ_LENGTH="$2"
            shift 2
            ;;

        # Controls whether secondary, supplementary, and QC-failed
        # alignments are excluded while creating the raw BAM.
        --primary-only)
            PRIMARY_ONLY="$2"
            shift 2
            ;;

        # Display usage instructions and exit successfully.
        -h|--help)
            usage
            exit 0
            ;;

        # Stop if the script encounters an unsupported argument.
        *)
            echo "[ERROR] Unknown argument: $1"
            usage
            exit 1
            ;;
    esac
done

# Apply default values when any optional variable is empty.
PRESET="${PRESET:-map-ont}"
THREADS="${THREADS:-8}"
SECONDARY="${SECONDARY:-yes}"
SORT_BAM="${SORT_BAM:-yes}"
INDEX_BAM="${INDEX_BAM:-yes}"
MIN_MAPQ="${MIN_MAPQ:-20}"
MIN_READ_LENGTH="${MIN_READ_LENGTH:-1000}"
PRIMARY_ONLY="${PRIMARY_ONLY:-no}"

# Confirm that all required file and directory arguments were provided.
if [[ -z "$FASTQ" || -z "$REFERENCE" || -z "$OUTPUT_DIR" || -z "$LOG_FILE" ]]; then
    echo "[ERROR] Missing required Minimap2 argument."
    usage
    exit 1
fi

# Confirm that the FASTQ file exists and is not empty.
if [[ ! -s "$FASTQ" ]]; then
    echo "[ERROR] FASTQ file is empty or not found: $FASTQ"
    exit 1
fi

# Confirm that the reference FASTA exists.
if [[ ! -f "$REFERENCE" ]]; then
    echo "[ERROR] Reference FASTA not found: $REFERENCE"
    exit 1
fi

# Validate the value used to control secondary alignments.
case "$SECONDARY" in
    yes|no) ;;
    *)
        echo "[ERROR] Invalid --secondary value: $SECONDARY"
        echo "[ERROR] Expected yes or no."
        exit 1
        ;;
esac

# Validate the BAM sorting option.
case "$SORT_BAM" in
    yes|no) ;;
    *)
        echo "[ERROR] Invalid --sort value: $SORT_BAM"
        echo "[ERROR] Expected yes or no."
        exit 1
        ;;
esac

# Validate the BAM indexing option.
case "$INDEX_BAM" in
    yes|no) ;;
    *)
        echo "[ERROR] Invalid --index value: $INDEX_BAM"
        echo "[ERROR] Expected yes or no."
        exit 1
        ;;
esac

# Validate the primary-only filtering option.
case "$PRIMARY_ONLY" in
    yes|no) ;;
    *)
        echo "[ERROR] Invalid --primary-only value: $PRIMARY_ONLY"
        echo "[ERROR] Expected yes or no."
        exit 1
        ;;
esac

# Require both BAM sorting and indexing.
#
# The downstream coverage, mapping-quality, read-length, and chromosome
# 12/rDNA checks depend on a coordinate-sorted and indexed BAM.
if [[ "$SORT_BAM" != "yes" || "$INDEX_BAM" != "yes" ]]; then
    echo "[ERROR] This nanopore workflow requires --sort yes and --index yes."
    echo "[ERROR] The requested coverage, MAPQ/read-length, and chromosome 12/rDNA checks run on the sorted indexed BAM."
    exit 1
fi

# Create the output directory and the directory containing the log file.
mkdir -p "$OUTPUT_DIR" "$(dirname "$LOG_FILE")"

# Derive the sample prefix from the FASTQ filename.
#
# Examples:
# sample.fastq becomes sample
# sample.fq becomes sample
FASTQ_BASENAME="$(basename "$FASTQ")"
FASTQ_PREFIX="${FASTQ_BASENAME%.fastq}"
FASTQ_PREFIX="${FASTQ_PREFIX%.fq}"
JOB_ID="${SLURM_JOB_ID:-manual}"

# FASTQs produced by dorado_basecall.sh are named sample_jobid.fastq.
# Reuse that job ID when present; otherwise use the current job ID.
SAMPLE_PREFIX="$FASTQ_PREFIX"
OUTPUT_JOB_ID="$JOB_ID"
if [[ "$FASTQ_PREFIX" == *_"$JOB_ID" ]]; then
    SAMPLE_PREFIX="${FASTQ_PREFIX%_"$JOB_ID"}"
elif [[ "$FASTQ_PREFIX" =~ ^(.+)_([0-9]+)$ ]]; then
    SAMPLE_PREFIX="${BASH_REMATCH[1]}"
    OUTPUT_JOB_ID="${BASH_REMATCH[2]}"
fi

# Define the raw unsorted BAM output.
RAW_BAM="$OUTPUT_DIR/${SAMPLE_PREFIX}_${OUTPUT_JOB_ID}.bam"

# Define the final sorted and indexed BAM output.
SORTED_INDEXED_BAM="$OUTPUT_DIR/${SAMPLE_PREFIX}.sorted.indexed_${OUTPUT_JOB_ID}.bam"

# Store the final BAM path in a general variable used throughout QC.
BAM_FILE="$SORTED_INDEXED_BAM"

# Define the separate alignment QC report.
QC_REPORT="$OUTPUT_DIR/${SAMPLE_PREFIX}_${OUTPUT_JOB_ID}.alignment_qc.txt"

# Construct the optional Minimap2 flag that disables secondary alignments.
#
# When SECONDARY is yes, this variable stays empty and Minimap2 uses
# its normal secondary-alignment behavior.
SECONDARY_FLAG=""
if [[ "$SECONDARY" == "no" ]]; then
    SECONDARY_FLAG="--secondary=no"
fi

# Build the samtools view arguments used while converting SAM to BAM.
#
# -b:
#   Produce BAM output.
#
# -S:
#   Treat the input as SAM.
SAMTOOLS_VIEW_ARGS=(-bS)

# Set the initial filtering arguments used during QC.
#
# -F 4 excludes unmapped reads while retaining aligned primary,
# secondary, and supplementary records.
PRIMARY_FILTER_ARGS=(-F 4)

# When primary-only mode is enabled:
#
# -F 2308 excludes:
# - unmapped reads,
# - secondary alignments,
# - supplementary alignments,
# - QC-failed reads.
#
# The same filter is used both while creating the raw BAM and during
# the MAPQ/read-length QC calculations.
if [[ "$PRIMARY_ONLY" == "yes" ]]; then
    SAMTOOLS_VIEW_ARGS+=(-F 2308)
    PRIMARY_FILTER_ARGS=(-F 2308)
fi

# Write the complete alignment configuration to the main log file.
{
    echo "========================================="
    echo "  MINIMAP2 ALIGNMENT"
    echo "========================================="
    echo "Started:         $(date)"
    echo "FASTQ:           $FASTQ"
    echo "Reference:       $REFERENCE"
    echo "Raw BAM:         $RAW_BAM"
    echo "Final BAM:       $BAM_FILE"
    echo "Preset:          $PRESET"
    echo "Threads:         $THREADS"
    echo "Secondary:       $SECONDARY"
    echo "Sort BAM:        $SORT_BAM"
    echo "Index BAM:       $INDEX_BAM"
    echo "Min MAPQ:        $MIN_MAPQ"
    echo "Min read length: $MIN_READ_LENGTH"
    echo "Primary only:    $PRIMARY_ONLY"
    echo "SLURM job ID:    ${SLURM_JOB_ID:-not_set}"
    echo "Output job ID:   $OUTPUT_JOB_ID"
    echo "========================================="
} > "$LOG_FILE"

# Create the header of the separate QC report.
echo "=========================================" > "$QC_REPORT"
echo "  ALIGNMENT QC REPORT" >> "$QC_REPORT"
echo "  $(date)" >> "$QC_REPORT"
echo "  Sample: $SAMPLE_PREFIX" >> "$QC_REPORT"
echo "  Job ID: $OUTPUT_JOB_ID" >> "$QC_REPORT"
echo "=========================================" >> "$QC_REPORT"

# Build the Minimap2 command as a Bash array.
#
# -a:
#   Produce SAM-formatted alignments.
#
# -x:
#   Select the requested alignment preset.
#
# -t:
#   Set the number of CPU threads.
MINIMAP2_ARGS=(-ax "$PRESET" -t "$THREADS")

# Add --secondary=no when secondary alignments were disabled.
if [[ -n "$SECONDARY_FLAG" ]]; then
    MINIMAP2_ARGS+=("$SECONDARY_FLAG")
fi

# Add the reference FASTA and FASTQ input paths to the command.
MINIMAP2_ARGS+=("$REFERENCE" "$FASTQ")

# Record the exact Minimap2 and samtools conversion command in the log.
#
# %q prints each Minimap2 argument in a shell-escaped form.
{
    echo ""
    echo "[INFO] Running Minimap2 alignment."
    printf '[INFO] Command: minimap2'
    printf ' %q' "${MINIMAP2_ARGS[@]}"
    echo " | samtools view ${SAMTOOLS_VIEW_ARGS[*]} -o $RAW_BAM"
} >> "$LOG_FILE"

# Align the nanopore reads to the reference genome.
#
# Minimap2 writes SAM records to standard output.
# The SAM records are immediately piped into samtools view, which:
# - converts them to BAM,
# - optionally applies primary-only filtering,
# - writes the raw BAM file.
#
# Error messages from both commands are appended to the log file.
minimap2 "${MINIMAP2_ARGS[@]}" 2>> "$LOG_FILE" | samtools view "${SAMTOOLS_VIEW_ARGS[@]}" -o "$RAW_BAM" - 2>> "$LOG_FILE"

# Confirm that alignment produced a nonempty raw BAM file.
if [[ ! -s "$RAW_BAM" ]]; then
    echo "[ERROR] BAM file is empty or was not created: $RAW_BAM" | tee -a "$LOG_FILE"
    exit 1
fi

# Coordinate-sort the raw BAM.
#
# Sorting is required before indexing and before region-based queries.
echo "[INFO] Sorting BAM: $BAM_FILE" >> "$LOG_FILE"
samtools sort -@ "$THREADS" -o "$BAM_FILE" "$RAW_BAM" >> "$LOG_FILE" 2>&1

# Create the BAM index.
#
# The index allows samtools to quickly access individual chromosomes
# and genomic regions.
echo "[INFO] Indexing BAM: $BAM_FILE" >> "$LOG_FILE"
samtools index "$BAM_FILE" >> "$LOG_FILE" 2>&1

# Run alignment quality-control checks.
#
# Output from this block is:
# - appended to the QC report through tee,
# - also appended to the main log file.
{
    echo ""
    echo "========================================="
    echo "  1. AVERAGE COVERAGE"
    echo "========================================="

    # Count all mapped alignment records.
    #
    # -F 4 excludes unmapped records.
    ALIGNED_READS=$(samtools view -c -F 4 "$BAM_FILE")
    echo "Aligned reads: $ALIGNED_READS"

    # Use samtools coverage to calculate coverage statistics for each
    # reference sequence.
    #
    # Column 7 contains mean depth for each reported reference region.
    # The awk command calculates the arithmetic mean of those values.
    samtools coverage "$BAM_FILE" | awk 'NR>1 {sum+=$7; n++} END {if (n > 0) print "Mean coverage across all regions:", sum/n "x"; else print "Mean coverage across all regions: N/A"}'

    echo ""

    # Print the full samtools coverage table for every reference contig.
    samtools coverage "$BAM_FILE"

    echo ""
    echo "========================================="
    echo "  2. MAPQ AND READ LENGTH FILTER"
    echo "========================================="

    # Count aligned reads using the selected primary filtering behavior.
    #
    # With PRIMARY_ONLY=no, this uses -F 4.
    # With PRIMARY_ONLY=yes, this uses -F 2308.
    TOTAL_ALIGNED=$(samtools view -c "${PRIMARY_FILTER_ARGS[@]}" "$BAM_FILE")

    # Count reads that meet both requirements:
    # - mapping quality greater than or equal to MIN_MAPQ,
    # - sequence length greater than or equal to MIN_READ_LENGTH.
    #
    # Column 10 of a SAM record contains the read sequence.
    PASSING=$(samtools view -q "$MIN_MAPQ" "${PRIMARY_FILTER_ARGS[@]}" "$BAM_FILE" | awk -v min_len="$MIN_READ_LENGTH" 'length($10) >= min_len' | wc -l)

    echo "Total aligned reads:                         $TOTAL_ALIGNED"
    echo "Reads with MAPQ >= $MIN_MAPQ and length >= $MIN_READ_LENGTH bp: $PASSING"

    # Calculate the percentage of aligned reads passing the MAPQ and
    # read-length thresholds.
    if [[ "$TOTAL_ALIGNED" -gt 0 ]]; then
        PCT=$(awk "BEGIN {printf \"%.1f\", ($PASSING / $TOTAL_ALIGNED) * 100}")
        echo "Percentage:                                  ${PCT}%"
    else
        echo "Percentage:                                  N/A"
    fi

    echo ""
    echo "========================================="
    echo "  3. SAMTOOLS FLAGSTAT"
    echo "========================================="

    # Print standard alignment summary statistics, including:
    # - total records,
    # - primary and secondary alignments,
    # - mapped reads,
    # - supplementary alignments,
    # - paired-read statistics when applicable.
    samtools flagstat "$BAM_FILE"

    echo ""
    echo "========================================="
    echo "  4. rDNA / CHROMOSOME 12 READ CHECK"
    echo "========================================="

    # Search BAM index statistics for a chromosome 12 contig name.
    #
    # The first search checks several common chromosome 12 identifiers,
    # including the W303 GenBank accession CM007975.1 and the S288C
    # accession NC_001144.5.
    CHR12_NAME=$(samtools idxstats "$BAM_FILE" | awk '$1 ~ /^(chr12|chrXII|XII|12|CM007975\.1|NC_001144\.5)$/ {print $1; exit}')

    # If no exact chromosome 12 name was found, perform a broader,
    # case-insensitive search for names containing chromosome 12 patterns.
    if [[ -z "$CHR12_NAME" ]]; then
        CHR12_NAME=$(samtools idxstats "$BAM_FILE" | awk 'tolower($1) ~ /(chr)?xii|chromosome[_ -]?12|chr12/ {print $1; exit}')
    fi

    # Run chromosome 12/rDNA-specific QC when a matching contig was found.
    if [[ -n "$CHR12_NAME" ]]; then
        # Count all mapped records aligned to chromosome 12.
        CHR12_READS=$(samtools view -c -F 4 "$BAM_FILE" "$CHR12_NAME")

        # Count chromosome 12 reads that pass both the MAPQ and
        # minimum read-length thresholds.
        CHR12_MAPQ_LENGTH_READS=$(samtools view -q "$MIN_MAPQ" -F 4 "$BAM_FILE" "$CHR12_NAME" | awk -v min_len="$MIN_READ_LENGTH" 'length($10) >= min_len' | wc -l)

        # Calculate the mean sequence length of mapped chromosome 12 reads.
        CHR12_MEAN_LENGTH=$(samtools view -F 4 "$BAM_FILE" "$CHR12_NAME" | awk '{sum+=length($10); n++} END {if (n > 0) printf "%.1f", sum/n; else print "N/A"}')

        echo "Chromosome 12 contig used: $CHR12_NAME"
        echo "Aligned reads on chromosome 12/rDNA: $CHR12_READS"
        echo "Chromosome 12 reads with MAPQ >= $MIN_MAPQ and length >= $MIN_READ_LENGTH bp: $CHR12_MAPQ_LENGTH_READS"
        echo "Mean chromosome 12 read length: $CHR12_MEAN_LENGTH bp"

        # Calculate what percentage of all mapped alignment records
        # are aligned to chromosome 12.
        if [[ "$ALIGNED_READS" -gt 0 ]]; then
            CHR12_PCT=$(awk "BEGIN {printf \"%.1f\", ($CHR12_READS / $ALIGNED_READS) * 100}")
            echo "Percent of aligned reads on chromosome 12/rDNA: ${CHR12_PCT}%"
        fi
    else
        # Warn the user when the script cannot recognize a chromosome
        # 12 contig from the reference sequence names.
        echo "WARNING: Could not identify a chromosome 12/rDNA contig in BAM index stats."
        echo "Review reference contig names if rDNA recovery is expected."
    fi

    echo ""
    echo "========================================="
    echo "All checks complete."
    echo "Completed: $(date)"
    echo "========================================="

# tee appends the QC output to the QC report while also passing it
# onward to be appended to the main alignment log.
} | tee -a "$QC_REPORT" >> "$LOG_FILE"

# Record the QC report path in the main log.
echo "[INFO] QC report saved to: $QC_REPORT" >> "$LOG_FILE"

# Print the final sorted and indexed BAM path to standard output.
#
# A parent workflow can capture this output and pass the BAM path
# into later steps such as DNAscent.
echo "$BAM_FILE"
