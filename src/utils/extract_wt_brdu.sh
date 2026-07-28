#!/bin/bash
#SBATCH --job-name=extract_WT_BrdU
#SBATCH --output=extract_WT_BrdU_%j.log
#SBATCH --error=extract_WT_BrdU_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=04:00:00

set -euo pipefail

module load modkit
module load samtools

BAM="MQNFGR_2_s_phase_020.Brdu_positive.threshold_0p5.bam"
OUTPUT="MQNFGR_2_s_phase_020.Brdu_positive.threshold_0p5.extract.tsv"

echo "Started: $(date)"
echo "Input BAM: $BAM"
echo "Output TSV: $OUTPUT"

/usr/bin/time -v modkit extract full \
    --mapped-only \
    --threads "${SLURM_CPUS_PER_TASK}" \
    "$BAM" \
    "$OUTPUT"

echo "Completed: $(date)"
echo "Output size:"
ls -lh "$OUTPUT"
