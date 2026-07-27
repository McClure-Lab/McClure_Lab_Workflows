#!/bin/bash

# ==============================================================================
# BAM MERGE, SORT, AND INDEX WORKFLOW
# ==============================================================================
#
# Purpose
# -------
# This script interactively selects two BAM files, submits a SLURM job, merges
# the selected BAM files, coordinate-sorts the merged alignments, creates a BAM
# index, validates the result, and compares the final alignment count with the
# sum of the two input BAM alignment counts.
#
# The workflow is designed for BAM files stored under:
#
#   /beevol/home/pineirok/workflows/data/bam
#
# Basic execution
# ---------------
#
# Run the script from the workflow project with:
#
#   bash src/utils/merge_bams.sh
#
# The script may also be launched from another directory because it resolves
# its own absolute path before submitting itself to SLURM.
#
# Workflow modes
# --------------
#
# This script has two execution modes:
#
#   1. Interactive launcher mode
#
#      This is the default mode when MERGE_BATCH_JOB is unset or is not equal
#      to 1.
#
#      The launcher:
#
#        - Lists available BAM files.
#        - Prompts for two different BAM files.
#        - Prompts for the final output BAM filename.
#        - Prevents the output filename from matching either input.
#        - Warns before replacing an existing BAM or index.
#        - Displays the selected configuration.
#        - Submits this same script to SLURM.
#
#   2. Submitted batch-job mode
#
#      This mode is entered when the environment variable below is exported:
#
#        MERGE_BATCH_JOB=1
#
#      The submitted job receives:
#
#        BAM1
#        BAM2
#        OUTPUT_BAM
#
#      through sbatch --export.
#
#      The batch job:
#
#        - Loads Samtools.
#        - Validates both BAM files.
#        - Compares their @SQ reference headers.
#        - Merges the alignments.
#        - Coordinate-sorts the merged BAM.
#        - Indexes the final BAM.
#        - Validates the final BAM.
#        - Counts input and output alignments.
#        - Compares the final count with the expected sum.
#
# Input compatibility requirement
# -------------------------------
#
# The two BAM files must have identical @SQ header lines.
#
# @SQ lines describe the reference sequences used during alignment, including
# chromosome or contig names and lengths. If these lines differ, the BAM files
# may have been aligned against different references or differently ordered
# reference sequences.
#
# The script stops before merging when the @SQ headers do not match.
#
# Output files
# ------------
#
# Final coordinate-sorted BAM:
#
#   workflow_root/data/bam/<output_name>.bam
#
# BAM index:
#
#   workflow_root/data/bam/<output_name>.bam.bai
#
# SLURM standard-output log:
#
#   merge_bams_<job_id>.log
#
# SLURM standard-error log:
#
#   merge_bams_<job_id>.err
#
# Temporary files
# ---------------
#
# A temporary unsorted merged BAM is created:
#
#   <output_prefix>.temporary_merged.bam
#
# Samtools sort may also create temporary files beginning with:
#
#   <output_prefix>.sort_tmp
#
# These temporary files are removed automatically when the script exits.
#
# Alignment-count verification
# ----------------------------
#
# The script calculates:
#
#   expected_count = alignments_in_BAM1 + alignments_in_BAM2
#
# It then compares this value with the number of alignments in the final BAM.
#
# This count includes all alignment records returned by samtools view -c. It
# does not filter secondary, supplementary, unmapped, or duplicate alignments.
#
# A mismatch produces a warning rather than deleting the completed BAM.
#
# ==============================================================================


# ==============================================================================
# SLURM RESOURCE REQUESTS
# ==============================================================================
#
# Requested resources:
#
#   Job name:        merge_bams
#   Tasks:           1
#   CPU cores:       32
#   Memory:          32 GB
#   Runtime limit:   4 hours
#
# %j is replaced by the SLURM job ID in the output and error filenames.
#
#SBATCH --job-name=merge_bams
#SBATCH --output=merge_bams_%j.log
#SBATCH --error=merge_bams_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=04:00:00


# ==============================================================================
# BASH SAFETY SETTINGS
# ==============================================================================
#
# -e
#   Exit when an unhandled command returns a nonzero status.
#
# -u
#   Treat references to unset variables as errors.
#
# -o pipefail
#   Treat a pipeline as failed when any command within the pipeline fails.
#
set -euo pipefail


# ==============================================================================
# FIXED WORKFLOW PATHS
# ==============================================================================

# Root directory of the workflow project.
WORKFLOW_DIR="/beevol/home/pineirok/workflows"

# Directory containing input BAM files and the final merged BAM.
BAM_DIR="${WORKFLOW_DIR}/data/bam"

# Resolve the absolute path to this script.
#
# The script submits itself to SLURM later. Using an absolute path ensures that
# the submission works regardless of the user's current working directory.
SCRIPT_PATH="$(realpath "$0")"


# ==============================================================================
# INTERACTIVE LAUNCHER MODE
# ==============================================================================
#
# This section runs when MERGE_BATCH_JOB is:
#
#   - Unset
#   - Empty
#   - Any value other than 1
#
# The launcher collects input from the user and submits the actual merge job.
#
if [[ "${MERGE_BATCH_JOB:-0}" != "1" ]]; then

    # --------------------------------------------------------------------------
    # Validate the BAM directory
    # --------------------------------------------------------------------------

    if [[ ! -d "$BAM_DIR" ]]; then
        echo "[ERROR] BAM directory does not exist:"
        echo "        $BAM_DIR"
        exit 1
    fi

    # --------------------------------------------------------------------------
    # Collect available BAM filenames
    # --------------------------------------------------------------------------
    #
    # mapfile reads the sorted command output into the BAM_FILES Bash array.
    #
    # find options:
    #
    #   -maxdepth 1
    #       Search only files directly inside BAM_DIR.
    #
    #   -type f
    #       Return regular files only.
    #
    #   -name "*.bam"
    #       Return files ending in .bam.
    #
    #   -printf "%f\n"
    #       Print only the filename, not the complete path.
    #
    mapfile -t BAM_FILES < <(
        find "$BAM_DIR" \
            -maxdepth 1 \
            -type f \
            -name "*.bam" \
            -printf "%f\n" \
            | sort
    )

    # At least two BAM files are required for a merge.
    if [[ "${#BAM_FILES[@]}" -lt 2 ]]; then
        echo "[ERROR] At least two BAM files are required in:"
        echo "        $BAM_DIR"
        exit 1
    fi

    # --------------------------------------------------------------------------
    # Display available BAM files
    # --------------------------------------------------------------------------

    echo
    echo "Available BAM files in:"
    echo "$BAM_DIR"
    echo

    # Display a one-based numbered menu.
    for i in "${!BAM_FILES[@]}"; do
        printf "  %3d) %s\n" \
            "$((i + 1))" \
            "${BAM_FILES[$i]}"
    done

    echo

    # --------------------------------------------------------------------------
    # Select the first BAM file
    # --------------------------------------------------------------------------
    #
    # The prompt repeats until the user provides an integer corresponding to
    # one of the displayed BAM files.
    #
    while true; do
        read -r -p "Select the first BAM file by number: " BAM1_NUMBER

        if [[ "$BAM1_NUMBER" =~ ^[0-9]+$ ]] &&
           (( BAM1_NUMBER >= 1 &&
              BAM1_NUMBER <= ${#BAM_FILES[@]} )); then

            # Convert the one-based menu selection into a zero-based array
            # index and construct the complete BAM path.
            BAM1="${BAM_DIR}/${BAM_FILES[$((BAM1_NUMBER - 1))]}"
            break
        fi

        echo "[ERROR] Enter a number between 1 and ${#BAM_FILES[@]}."
    done

    # --------------------------------------------------------------------------
    # Select the second BAM file
    # --------------------------------------------------------------------------
    #
    # The second BAM must be a different file from the first BAM.
    #
    while true; do
        read -r -p "Select the second BAM file by number: " BAM2_NUMBER

        if [[ "$BAM2_NUMBER" =~ ^[0-9]+$ ]] &&
           (( BAM2_NUMBER >= 1 &&
              BAM2_NUMBER <= ${#BAM_FILES[@]} )); then

            BAM2="${BAM_DIR}/${BAM_FILES[$((BAM2_NUMBER - 1))]}"

            # Prevent the same BAM from being selected twice.
            if [[ "$BAM2" == "$BAM1" ]]; then
                echo "[ERROR] Select a different BAM file."
                continue
            fi

            break
        fi

        echo "[ERROR] Enter a number between 1 and ${#BAM_FILES[@]}."
    done

    # --------------------------------------------------------------------------
    # Collect and validate the output BAM filename
    # --------------------------------------------------------------------------

    while true; do
        echo
        read -r -p "Enter the final combined BAM filename: " OUTPUT_NAME

        # Remove any directory components entered by the user.
        #
        # This ensures the final BAM is always written under BAM_DIR.
        OUTPUT_NAME="$(basename "$OUTPUT_NAME")"

        # Reject a blank output filename.
        if [[ -z "$OUTPUT_NAME" ]]; then
            echo "[ERROR] Output filename cannot be empty."
            continue
        fi

        # Add the .bam extension when the user omits it.
        if [[ "$OUTPUT_NAME" != *.bam ]]; then
            OUTPUT_NAME="${OUTPUT_NAME}.bam"
        fi

        # Construct the complete output path.
        OUTPUT_BAM="${BAM_DIR}/${OUTPUT_NAME}"

        # Prevent either input BAM from being overwritten by the output.
        if [[ "$OUTPUT_BAM" == "$BAM1" ||
              "$OUTPUT_BAM" == "$BAM2" ]]; then
            echo "[ERROR] The output name cannot match either input BAM."
            continue
        fi

        break
    done

    # --------------------------------------------------------------------------
    # Display the selected merge configuration
    # --------------------------------------------------------------------------

    echo
    echo "Merge configuration"
    echo "=================================================="
    echo "First BAM:"
    echo "  $BAM1"
    echo
    echo "Second BAM:"
    echo "  $BAM2"
    echo
    echo "Final output BAM:"
    echo "  $OUTPUT_BAM"
    echo
    echo "Resources:"
    echo "  CPUs:   32"
    echo "  Memory: 32 GB"
    echo "=================================================="
    echo

    # --------------------------------------------------------------------------
    # Check for an existing output BAM or index
    # --------------------------------------------------------------------------
    #
    # Samtools index files may appear as either:
    #
    #   output.bam.bai
    #
    # or:
    #
    #   output.bai
    #
    # Both naming patterns are checked.
    #
    if [[ -e "$OUTPUT_BAM" ||
          -e "${OUTPUT_BAM}.bai" ||
          -e "${OUTPUT_BAM%.bam}.bai" ]]; then

        echo "[WARNING] The output BAM or index already exists:"
        echo "          $OUTPUT_BAM"
        echo

        read -r -p "Overwrite the existing output? [y/N]: " OVERWRITE

        # Only an explicit y or Y authorizes replacement.
        if [[ ! "$OVERWRITE" =~ ^[Yy]$ ]]; then
            echo "Job cancelled."
            exit 0
        fi
    fi

    # --------------------------------------------------------------------------
    # Final submission confirmation
    # --------------------------------------------------------------------------

    read -r -p "Submit this merge job? [y/N]: " CONFIRM

    if [[ ! "$CONFIRM" =~ ^[Yy]$ ]]; then
        echo "Job cancelled."
        exit 0
    fi

    # --------------------------------------------------------------------------
    # Submit this same script to SLURM
    # --------------------------------------------------------------------------
    #
    # --chdir
    #   Run the submitted job from WORKFLOW_DIR.
    #
    # --export
    #   Export the existing environment and explicitly pass:
    #
    #     MERGE_BATCH_JOB=1
    #         Selects batch-job mode.
    #
    #     BAM1
    #         Complete path to the first input BAM.
    #
    #     BAM2
    #         Complete path to the second input BAM.
    #
    #     OUTPUT_BAM
    #         Complete path to the final BAM.
    #
    sbatch \
        --chdir="$WORKFLOW_DIR" \
        --export=ALL,MERGE_BATCH_JOB=1,BAM1="$BAM1",BAM2="$BAM2",OUTPUT_BAM="$OUTPUT_BAM" \
        "$SCRIPT_PATH"

    # Stop the launcher after successful submission so it does not continue
    # into the batch-job section in the current shell.
    exit 0
fi


# ==============================================================================
# SUBMITTED SLURM JOB MODE
# ==============================================================================
#
# This section runs only when MERGE_BATCH_JOB=1.
#
# BAM1, BAM2, and OUTPUT_BAM are expected to have been exported by the launcher.
# The parameter-expansion checks below stop immediately with a descriptive
# message when one of the required variables is missing or empty.
#

: "${BAM1:?BAM1 was not passed to the SLURM job}"
: "${BAM2:?BAM2 was not passed to the SLURM job}"
: "${OUTPUT_BAM:?OUTPUT_BAM was not passed to the SLURM job}"


# ==============================================================================
# LOAD REQUIRED SOFTWARE
# ==============================================================================

# Load Samtools through the cluster module system.
module load samtools


# ==============================================================================
# JOB SETTINGS AND TEMPORARY PATHS
# ==============================================================================

# Use the CPU count assigned by SLURM.
#
# Fall back to 32 when SLURM_CPUS_PER_TASK is unavailable.
THREADS="${SLURM_CPUS_PER_TASK:-32}"

# Directory where the final BAM will be written.
OUTPUT_DIR="$(dirname "$OUTPUT_BAM")"

# Output filename without the .bam extension.
OUTPUT_BASE="$(basename "$OUTPUT_BAM" .bam)"

# Temporary unsorted BAM produced by samtools merge.
TEMP_MERGED_BAM="${OUTPUT_DIR}/${OUTPUT_BASE}.temporary_merged.bam"

# Prefix used by samtools sort for temporary sorting files.
TEMP_SORT_PREFIX="${OUTPUT_DIR}/${OUTPUT_BASE}.sort_tmp"


# ==============================================================================
# FUNCTION: cleanup
# ==============================================================================
#
# Purpose:
#   Remove temporary merge and sorting files.
#
# The function is registered with an EXIT trap, so cleanup is attempted after:
#
#   - Successful completion.
#   - Validation failure.
#   - Merge failure.
#   - Sorting failure.
#   - Indexing failure.
#
# The final output BAM and its index are not removed by this function.
#
cleanup() {
    # Remove the temporary unsorted merged BAM.
    rm -f "$TEMP_MERGED_BAM"

    # Remove any temporary BAM files generated by samtools sort using the
    # configured temporary prefix.
    rm -f "${TEMP_SORT_PREFIX}"*.bam
}

# Run cleanup whenever the script exits.
trap cleanup EXIT


# ==============================================================================
# REPORT THE BATCH-JOB CONFIGURATION
# ==============================================================================

echo "=================================================="
echo "BAM merge job"
echo "=================================================="
echo "[INFO] Job ID:      ${SLURM_JOB_ID:-unknown}"
echo "[INFO] First BAM:   $BAM1"
echo "[INFO] Second BAM:  $BAM2"
echo "[INFO] Output BAM:  $OUTPUT_BAM"
echo "[INFO] CPUs:        $THREADS"
echo "[INFO] Memory:      32 GB"
echo "=================================================="


# ==============================================================================
# VALIDATE INPUT FILE EXISTENCE
# ==============================================================================

if [[ ! -f "$BAM1" ]]; then
    echo "[ERROR] First BAM does not exist:"
    echo "        $BAM1"
    exit 1
fi

if [[ ! -f "$BAM2" ]]; then
    echo "[ERROR] Second BAM does not exist:"
    echo "        $BAM2"
    exit 1
fi

# Create the output directory when it does not already exist.
mkdir -p "$OUTPUT_DIR"


# ==============================================================================
# VALIDATE INPUT BAM STRUCTURE
# ==============================================================================
#
# samtools quickcheck performs a fast structural validation by checking the BAM
# header and end-of-file information. It does not inspect every alignment.
#
# Because set -e is enabled, the workflow stops if either command fails.
#

echo "[INFO] Validating input BAM files..."

samtools quickcheck -v "$BAM1"
samtools quickcheck -v "$BAM2"


# ==============================================================================
# COMPARE REFERENCE-SEQUENCE HEADERS
# ==============================================================================
#
# samtools view -H prints each BAM header.
#
# grep '^@SQ' retains only reference-sequence records.
#
# diff -q compares the @SQ records without printing their full differences.
#
# Matching @SQ records help confirm that:
#
#   - Both BAM files use the same chromosome/contig names.
#   - Both BAM files use the same reference lengths.
#   - The reference sequences appear in the same order.
#
# The script does not compare all other header records, such as:
#
#   @RG read groups
#   @PG program records
#   @CO comments
#
echo "[INFO] Comparing reference sequence headers..."

if ! diff -q \
    <(samtools view -H "$BAM1" | grep '^@SQ') \
    <(samtools view -H "$BAM2" | grep '^@SQ') \
    >/dev/null; then

    echo "[ERROR] The BAM files do not have matching @SQ headers."
    echo "[ERROR] They may have been aligned to different references."
    exit 1
fi

echo "[INFO] Reference sequence headers match."


# ==============================================================================
# REMOVE PREVIOUS OUTPUTS AND STALE TEMPORARY FILES
# ==============================================================================
#
# The interactive launcher already asked the user for overwrite confirmation
# when an existing output was detected.
#
# This removal occurs inside the compute job immediately before processing.
#

rm -f \
    "$OUTPUT_BAM" \
    "${OUTPUT_BAM}.bai" \
    "${OUTPUT_BAM%.bam}.bai" \
    "$TEMP_MERGED_BAM"


# ==============================================================================
# MERGE INPUT BAM FILES
# ==============================================================================
#
# samtools merge options:
#
#   -@
#       Number of worker threads.
#
#   -f
#       Overwrite the temporary merged output if it already exists.
#
# The merged BAM produced at this stage is not assumed to be coordinate sorted.
#
echo "[INFO] Merging BAM files..."

samtools merge \
    -@ "$THREADS" \
    -f \
    "$TEMP_MERGED_BAM" \
    "$BAM1" \
    "$BAM2"


# ==============================================================================
# COORDINATE-SORT THE MERGED BAM
# ==============================================================================
#
# samtools sort options:
#
#   -@
#       Number of worker threads.
#
#   -m 750M
#       Maximum memory used per sorting thread.
#
#       With 32 threads, Samtools may theoretically request approximately:
#
#         32 x 750 MB = 24,000 MB
#
#       Additional overhead remains within the requested 32 GB SLURM memory.
#
#   -T
#       Prefix for temporary sorting files.
#
#   -o
#       Final coordinate-sorted BAM path.
#
echo "[INFO] Sorting merged BAM..."

samtools sort \
    -@ "$THREADS" \
    -m 750M \
    -T "$TEMP_SORT_PREFIX" \
    -o "$OUTPUT_BAM" \
    "$TEMP_MERGED_BAM"


# ==============================================================================
# INDEX THE FINAL BAM
# ==============================================================================
#
# samtools index creates an index used for random genomic-region access.
#
# With the current command, the typical index path is:
#
#   <output>.bam.bai
#
echo "[INFO] Indexing final BAM..."

samtools index \
    -@ "$THREADS" \
    "$OUTPUT_BAM"


# ==============================================================================
# VALIDATE THE FINAL BAM
# ==============================================================================
#
# Confirm that the coordinate-sorted output has a readable BAM header and a
# valid end-of-file block.
#
echo "[INFO] Validating final BAM..."

samtools quickcheck -v "$OUTPUT_BAM"


# ==============================================================================
# COUNT INPUT AND OUTPUT ALIGNMENTS
# ==============================================================================
#
# samtools view -c counts alignment records.
#
# No flag filters are applied here. The counts therefore include every
# alignment record stored in each BAM, including records that may be:
#
#   - Unmapped
#   - Secondary
#   - Supplementary
#   - Duplicate-marked
#
# This verification checks whether the final BAM contains the expected total
# number of alignment records from both input BAMs.
#

BAM1_COUNT="$(samtools view -@ "$THREADS" -c "$BAM1")"
BAM2_COUNT="$(samtools view -@ "$THREADS" -c "$BAM2")"
OUTPUT_COUNT="$(samtools view -@ "$THREADS" -c "$OUTPUT_BAM")"

# Calculate the expected number of alignment records after concatenating the
# two input BAM alignment sets.
EXPECTED_COUNT=$((BAM1_COUNT + BAM2_COUNT))


# ==============================================================================
# REPORT ALIGNMENT COUNTS
# ==============================================================================

echo
echo "Alignment counts"
echo "=================================================="
echo "First BAM:      $BAM1_COUNT"
echo "Second BAM:     $BAM2_COUNT"
echo "Expected total: $EXPECTED_COUNT"
echo "Combined BAM:   $OUTPUT_COUNT"
echo "=================================================="


# ==============================================================================
# VERIFY THE FINAL ALIGNMENT COUNT
# ==============================================================================
#
# A successful merge should normally produce:
#
#   OUTPUT_COUNT == BAM1_COUNT + BAM2_COUNT
#
# A mismatch produces a warning but does not remove the completed output.
#

if [[ "$OUTPUT_COUNT" -eq "$EXPECTED_COUNT" ]]; then
    echo "[INFO] Combined alignment count matches the input sum."
else
    echo "[WARNING] Combined alignment count does not match the input sum."
fi


# ==============================================================================
# SUCCESSFUL COMPLETION
# ==============================================================================

echo
echo "[SUCCESS] Combined BAM created:"
echo "          $OUTPUT_BAM"
echo
echo "[SUCCESS] BAM index created:"
echo "          ${OUTPUT_BAM}.bai"