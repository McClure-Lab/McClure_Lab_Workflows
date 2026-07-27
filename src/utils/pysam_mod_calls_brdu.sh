#!/bin/bash
#SBATCH --job-name=pysam_BrdU_calls
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=12:00:00

set -euo pipefail


# =============================================================================
# PySAM BrdU-positive read extraction
#
# Run from the workflow root:
#
#   bash src/utils/pysam_mod_calls_brdu.sh
#
# Interactive stage:
#   1. Lists BAM files under <workflow_root>/data/bam.
#   2. Prompts for an input BAM.
#   3. Prompts for a BrdU modification threshold.
#   4. Submits this script as an SBATCH job.
#
# Analysis stage:
#   1. Streams through the BAM one record at a time using PySAM.
#   2. Parses the MM/ML tags using read.modified_bases.
#   3. Selects modification code "b", regardless of whether the BAM encodes
#      it as N+b or T+b.
#   4. Keeps a complete primary mapped read when at least one BrdU call has
#      an ML value greater than or equal to the selected threshold.
#   5. Sorts and indexes the resulting BAM.
#
# Your BAM uses:
#
#   MM:Z:N+b?,...;N+e?,...;
#
# Therefore, this script filters by modification code "b" and does not
# require the canonical-base component to be T.
#
# Example outputs:
#
#   MQNFGR_3_mitosis_WT_020.Brdu_positive.pysam.threshold_0p5.bam
#   MQNFGR_3_mitosis_WT_020.Brdu_positive.pysam.threshold_0p5.bam.bai
#   MQNFGR_3_mitosis_WT_020.Brdu_positive.pysam.threshold_0p5.read_ids.txt
#   MQNFGR_3_mitosis_WT_020.Brdu_positive.pysam.threshold_0p5.summary.tsv
# =============================================================================


# =============================================================================
# Message functions
# =============================================================================

print_info() {
    printf '[INFO] %s\n' "$1"
}


print_warning() {
    printf '[WARNING] %s\n' "$1" >&2
}


print_error() {
    printf '[ERROR] %s\n' "$1" >&2
}


# =============================================================================
# Workflow-root detection
# =============================================================================

determine_workflow_root() {
    local script_directory

    if [[ -d "$PWD/data/bam" && -d "$PWD/src/utils" ]]; then
        printf '%s\n' "$PWD"
        return 0
    fi

    script_directory="$(
        cd "$(dirname "${BASH_SOURCE[0]}")" &&
        pwd
    )"

    cd "${script_directory}/../.." && pwd
}


WORKFLOW_ROOT_ARGUMENT=""

for ((argument_index = 1; argument_index <= $#; argument_index++)); do
    if [[ "${!argument_index}" == "--workflow-root" ]]; then
        next_argument_index=$((argument_index + 1))
        WORKFLOW_ROOT_ARGUMENT="${!next_argument_index:-}"
        break
    fi
done


if [[ -n "$WORKFLOW_ROOT_ARGUMENT" ]]; then
    WORKFLOW_ROOT="$(
        cd "$WORKFLOW_ROOT_ARGUMENT" &&
        pwd
    )"
elif [[ -n "${PYSAM_BRDU_WORKFLOW_ROOT:-}" ]]; then
    WORKFLOW_ROOT="$PYSAM_BRDU_WORKFLOW_ROOT"
else
    WORKFLOW_ROOT="$(determine_workflow_root)"
fi


# =============================================================================
# Workflow paths
# =============================================================================

BAM_DIRECTORY="${WORKFLOW_ROOT}/data/bam"
LOG_DIRECTORY="${WORKFLOW_ROOT}/logs/pysam_mod_calls_brdu"
VENV_DIRECTORY="${WORKFLOW_ROOT}/.venv"
SCRIPT_PATH="${WORKFLOW_ROOT}/src/utils/pysam_mod_calls_brdu.sh"

JOB_MEMORY="32G"

mkdir -p "$LOG_DIRECTORY"


# =============================================================================
# Naming functions
# =============================================================================

strip_bam_suffixes() {
    local filename="$1"

    filename="${filename%.bam}"

    filename="${filename%.sorted.indexed.BrdU.detect}"
    filename="${filename%.sorted.indexed.Brdu.detect}"
    filename="${filename%.sorted.indexed.brdu.detect}"

    filename="${filename%.indexed.BrdU.detect}"
    filename="${filename%.indexed.Brdu.detect}"
    filename="${filename%.indexed.brdu.detect}"

    filename="${filename%.sorted.BrdU.detect}"
    filename="${filename%.sorted.Brdu.detect}"
    filename="${filename%.sorted.brdu.detect}"

    filename="${filename%.BrdU.detect}"
    filename="${filename%.Brdu.detect}"
    filename="${filename%.brdu.detect}"

    filename="${filename%.sorted.indexed}"
    filename="${filename%.indexed}"
    filename="${filename%.sorted}"

    printf '%s\n' "$filename"
}


normalize_threshold_for_filename() {
    local threshold="$1"

    # Examples:
    #
    #   0.50  -> 0p5
    #   0.60  -> 0p6
    #   0.625 -> 0p625

    printf '%s\n' "$threshold" |
        sed -E 's/0+$//; s/\.$//; s/\./p/g'
}


validate_threshold() {
    local threshold="$1"

    awk -v value="$threshold" 'BEGIN {
        if (value ~ /^[0-9]+([.][0-9]+)?$/ && value >= 0 && value <= 1) {
            exit 0
        }

        exit 1
    }'
}


check_required_command() {
    local command_name="$1"

    if ! command -v "$command_name" >/dev/null 2>&1; then
        print_error "Required command is unavailable: ${command_name}"
        exit 1
    fi
}


# =============================================================================
# Module and Python environment setup
# =============================================================================

load_required_modules() {
    if [[ -f /etc/profile.d/modules.sh ]]; then
        # shellcheck source=/dev/null
        source /etc/profile.d/modules.sh
    fi

    if command -v module >/dev/null 2>&1; then
        module load python
        module load samtools
    else
        print_warning "Environment modules are unavailable."
        print_warning "Continuing with Python and samtools from PATH."
    fi
}


activate_python_environment() {
    load_required_modules

    check_required_command python
    check_required_command samtools

    if [[ ! -d "$VENV_DIRECTORY" ]]; then
        print_info "Creating workflow Python virtual environment:"
        print_info "$VENV_DIRECTORY"

        python -m venv "$VENV_DIRECTORY"
    fi

    # shellcheck source=/dev/null
    source "${VENV_DIRECTORY}/bin/activate"

    if ! python -c 'import pysam' >/dev/null 2>&1; then
        print_info "Installing pysam into the workflow virtual environment."

        python -m pip install --upgrade pip
        python -m pip install pysam
    fi
}


# =============================================================================
# Usage
# =============================================================================

usage() {
    cat <<'EOF'
Usage:

  bash src/utils/pysam_mod_calls_brdu.sh

Optional noninteractive submission:

  bash src/utils/pysam_mod_calls_brdu.sh BAM_FILENAME MOD_THRESHOLD

Example:

  bash src/utils/pysam_mod_calls_brdu.sh \
      MQNFGR_3_mitosis_WT_020.sorted.indexed.BrdU.detect.bam \
      0.50
EOF
}


# =============================================================================
# Interactive submission mode
# =============================================================================

submit_workflow() {
    local bam_input="${1:-}"
    local mod_threshold="${2:-}"

    local bam_files=()
    local bam_selection
    local selected_bam
    local sample_prefix
    local threshold_tag

    local output_prefix
    local output_bam
    local output_index
    local output_read_ids
    local output_summary

    local existing_output_found
    local slurm_log_prefix
    local submission_output
    local job_id

    printf '\n'
    printf '%s\n' '============================================================'
    printf '%s\n' 'Extract BrdU-positive reads using PySAM'
    printf '%s\n' '============================================================'
    printf '\n'

    print_info "Workflow root: ${WORKFLOW_ROOT}"
    print_info "BAM directory: ${BAM_DIRECTORY}"

    if [[ ! -d "$BAM_DIRECTORY" ]]; then
        print_error "BAM directory does not exist:"
        print_error "$BAM_DIRECTORY"
        exit 1
    fi

    if [[ ! -f "$SCRIPT_PATH" ]]; then
        print_error "Could not locate this script:"
        print_error "$SCRIPT_PATH"
        exit 1
    fi

    if ! command -v sbatch >/dev/null 2>&1; then
        print_error "The sbatch command is unavailable."
        exit 1
    fi


    # =========================================================================
    # Prompt for input BAM
    # =========================================================================

    if [[ -z "$bam_input" ]]; then
        mapfile -t bam_files < <(
            find "$BAM_DIRECTORY" \
                -maxdepth 1 \
                -type f \
                -name '*.bam' \
                -printf '%f\n' |
            sort
        )

        if (( ${#bam_files[@]} == 0 )); then
            print_error "No BAM files were found under:"
            print_error "$BAM_DIRECTORY"
            exit 1
        fi

        printf '\n'
        printf '%s\n' 'Available BAM files:'
        printf '\n'

        for index in "${!bam_files[@]}"; do
            printf '  %3d) %s\n' \
                "$((index + 1))" \
                "${bam_files[$index]}"
        done

        printf '\n'

        while true; do
            read -r -p "Select a BAM file by number: " bam_selection

            if [[ ! "$bam_selection" =~ ^[0-9]+$ ]]; then
                print_warning "Enter a numeric selection."
                continue
            fi

            if (( bam_selection < 1 || bam_selection > ${#bam_files[@]} )); then
                print_warning \
                    "Selection must be between 1 and ${#bam_files[@]}."
                continue
            fi

            bam_input="${bam_files[$((bam_selection - 1))]}"
            break
        done
    fi

    selected_bam="$(basename "$bam_input")"

    if [[ ! -f "${BAM_DIRECTORY}/${selected_bam}" ]]; then
        print_error "Input BAM does not exist:"
        print_error "${BAM_DIRECTORY}/${selected_bam}"
        exit 1
    fi


    # =========================================================================
    # Prompt for threshold
    # =========================================================================

    if [[ -z "$mod_threshold" ]]; then
        printf '\n'

        while true; do
            read -r -p \
                "Enter the BrdU modification threshold [0.50]: " \
                mod_threshold

            mod_threshold="${mod_threshold:-0.50}"

            if validate_threshold "$mod_threshold"; then
                break
            fi

            print_warning \
                "The threshold must be a numeric value between 0 and 1."
        done
    elif ! validate_threshold "$mod_threshold"; then
        print_error "Invalid modification threshold: ${mod_threshold}"
        exit 1
    fi


    # =========================================================================
    # Build output paths
    # =========================================================================

    sample_prefix="$(strip_bam_suffixes "$selected_bam")"
    threshold_tag="$(normalize_threshold_for_filename "$mod_threshold")"

    output_prefix="${BAM_DIRECTORY}/${sample_prefix}.Brdu_positive.pysam.threshold_${threshold_tag}"

    output_bam="${output_prefix}.bam"
    output_index="${output_bam}.bai"
    output_read_ids="${output_prefix}.read_ids.txt"
    output_summary="${output_prefix}.summary.tsv"

    printf '\n'
    printf '%s\n' '------------------------------------------------------------'
    printf '%s\n' 'Submission summary'
    printf '%s\n' '------------------------------------------------------------'
    printf 'Input BAM:       %s\n' "$selected_bam"
    printf 'BrdU code:       b\n'
    printf 'Observed format: N+b\n'
    printf 'Mod threshold:   %s\n' "$mod_threshold"
    printf 'Output BAM:      %s\n' "$output_bam"
    printf 'Read-ID file:    %s\n' "$output_read_ids"
    printf 'Summary TSV:     %s\n' "$output_summary"
    printf '%s\n' '------------------------------------------------------------'
    printf '\n'


    # =========================================================================
    # Prevent accidental overwriting
    # =========================================================================

    existing_output_found=0

    for output_file in \
        "$output_bam" \
        "$output_index" \
        "$output_read_ids" \
        "$output_summary"
    do
        if [[ -e "$output_file" ]]; then
            print_error "Output already exists: ${output_file}"
            existing_output_found=1
        fi
    done

    if (( existing_output_found == 1 )); then
        print_error "Rename or remove existing outputs before rerunning."
        exit 1
    fi


    # =========================================================================
    # Submit job
    # =========================================================================

    slurm_log_prefix="${LOG_DIRECTORY}/${sample_prefix}.Brdu_positive.pysam.threshold_${threshold_tag}"

    submission_output="$(
        sbatch \
            --job-name="pysam_BrdU_${sample_prefix}" \
            --chdir="$WORKFLOW_ROOT" \
            --mem="$JOB_MEMORY" \
            --output="${slurm_log_prefix}.%j.log" \
            --error="${slurm_log_prefix}.%j.err" \
            --export=ALL,PYSAM_BRDU_WORKFLOW_ROOT="$WORKFLOW_ROOT" \
            "$SCRIPT_PATH" \
            --run-analysis \
            --workflow-root "$WORKFLOW_ROOT" \
            "$selected_bam" \
            "$mod_threshold"
    )"

    printf '%s\n' "$submission_output"

    job_id="$(
        printf '%s\n' "$submission_output" |
        awk '{print $NF}'
    )"

    printf '\n'
    print_info "Submitted job ID: ${job_id}"

    print_info "Expected output BAM:"
    printf '  %s\n' "$output_bam"

    print_info "Expected read-ID file:"
    printf '  %s\n' "$output_read_ids"

    print_info "Expected summary TSV:"
    printf '  %s\n' "$output_summary"

    print_info "SLURM logs:"
    printf '  %s.%s.log\n' "$slurm_log_prefix" "$job_id"
    printf '  %s.%s.err\n' "$slurm_log_prefix" "$job_id"
}


# =============================================================================
# Submitted analysis mode
# =============================================================================

run_analysis() {
    local selected_bam="$1"
    local mod_threshold="$2"

    local input_bam
    local sample_prefix
    local threshold_tag

    local output_prefix
    local output_bam
    local output_read_ids
    local output_summary

    local threads
    local temp_base_directory
    local temp_directory
    local temp_unsorted_bam
    local temp_sorted_bam
    local temp_read_ids
    local temp_summary

    local output_count
    local read_id_count

    if [[ -z "${SLURM_JOB_ID:-}" ]]; then
        print_error \
            "The --run-analysis mode must run through a submitted SLURM job."
        exit 1
    fi

    if ! validate_threshold "$mod_threshold"; then
        print_error "Invalid modification threshold: ${mod_threshold}"
        exit 1
    fi

    input_bam="${BAM_DIRECTORY}/${selected_bam}"

    if [[ ! -f "$input_bam" ]]; then
        print_error "Input BAM does not exist:"
        print_error "$input_bam"
        exit 1
    fi

    sample_prefix="$(strip_bam_suffixes "$selected_bam")"
    threshold_tag="$(normalize_threshold_for_filename "$mod_threshold")"

    output_prefix="${BAM_DIRECTORY}/${sample_prefix}.Brdu_positive.pysam.threshold_${threshold_tag}"

    output_bam="${output_prefix}.bam"
    output_read_ids="${output_prefix}.read_ids.txt"
    output_summary="${output_prefix}.summary.tsv"

    threads="${SLURM_CPUS_PER_TASK:-8}"

    temp_base_directory="${SLURM_TMPDIR:-/tmp}"

    temp_directory="$(
        mktemp \
            -d \
            "${temp_base_directory}/pysam_mod_calls_brdu.${SLURM_JOB_ID}.XXXXXX"
    )"

    temp_unsorted_bam="${temp_directory}/Brdu_positive.unsorted.bam"
    temp_sorted_bam="${temp_directory}/Brdu_positive.sorted.bam"
    temp_read_ids="${temp_directory}/Brdu_positive.read_ids.txt"
    temp_summary="${temp_directory}/Brdu_positive.summary.tsv"


    # =========================================================================
    # Cleanup
    # =========================================================================

    cleanup() {
        rm -rf "$temp_directory"
    }

    trap cleanup EXIT


    # =========================================================================
    # Environment
    # =========================================================================

    activate_python_environment

    check_required_command samtools
    check_required_command wc
    check_required_command mktemp


    # =========================================================================
    # Prevent overwrite
    # =========================================================================

    for output_file in \
        "$output_bam" \
        "${output_bam}.bai" \
        "$output_read_ids" \
        "$output_summary"
    do
        if [[ -e "$output_file" ]]; then
            print_error "Output already exists: ${output_file}"
            exit 1
        fi
    done


    # =========================================================================
    # Print job information
    # =========================================================================

    printf '\n'
    printf '%s\n' '============================================================'
    printf '%s\n' 'PySAM BrdU-positive read extraction'
    printf '%s\n' '============================================================'
    printf '\n'

    printf 'Analysis date:                 %s\n' \
        "$(date '+%Y-%m-%d %H:%M:%S')"

    printf 'SLURM job ID:                  %s\n' \
        "$SLURM_JOB_ID"

    printf 'Workflow root:                 %s\n' \
        "$WORKFLOW_ROOT"

    printf 'Input BAM:                     %s\n' \
        "$input_bam"

    printf 'Observed MM encoding:          N+b and N+e\n'
    printf 'Selected modification code:    b\n'

    printf 'BrdU modification threshold:   %s\n' \
        "$mod_threshold"

    printf 'Threads:                       %s\n' \
        "$threads"

    printf 'Requested memory:              %s\n' \
        "$JOB_MEMORY"

    printf 'Output BAM:                    %s\n' \
        "$output_bam"

    printf 'Read-ID file:                  %s\n' \
        "$output_read_ids"

    printf 'Summary TSV:                   %s\n' \
        "$output_summary"

    printf 'Python version:                %s\n' \
        "$(python --version 2>&1)"

    printf 'PySAM version:                 %s\n' \
        "$(python -c 'import pysam; print(pysam.__version__)')"

    printf 'Samtools version:              %s\n' \
        "$(samtools --version | head -n 1)"

    printf '\n'


    # =========================================================================
    # Stream the BAM one read at a time
    # =========================================================================

    python - \
        "$input_bam" \
        "$temp_unsorted_bam" \
        "$temp_read_ids" \
        "$temp_summary" \
        "$mod_threshold" \
        "$threads" <<'PYTHON'
#!/usr/bin/env python3

import math
import os
import sys
from typing import Any

import pysam


def normalize_modification_code(value: Any) -> str:
    """
    Normalize the modification code returned by PySAM.

    Depending on PySAM/HTSlib, code b may appear as:

        "b"
        b"b"
        98

    Integer 98 is the ASCII value for b.
    """
    if isinstance(value, bytes):
        return value.decode(
            "utf-8",
            errors="replace",
        )

    if isinstance(value, int):
        try:
            return chr(value)
        except (ValueError, OverflowError):
            return str(value)

    return str(value)


def normalize_canonical_base(value: Any) -> str:
    """
    Normalize the canonical-base component of a modified-base key.
    """
    if isinstance(value, bytes):
        return value.decode(
            "utf-8",
            errors="replace",
        ).upper()

    if isinstance(value, int):
        try:
            return chr(value).upper()
        except (ValueError, OverflowError):
            return str(value).upper()

    return str(value).upper()


def get_modified_bases(
    read: pysam.AlignedSegment,
) -> dict:
    """
    Return PySAM's parsed MM/ML dictionary.

    modified_bases is preferred. modified_bases_forward is used only as a
    fallback when modified_bases is unavailable or empty.
    """
    try:
        modified_bases = read.modified_bases
    except (KeyError, ValueError, TypeError):
        modified_bases = None

    if modified_bases:
        return modified_bases

    try:
        modified_bases_forward = read.modified_bases_forward
    except (KeyError, ValueError, TypeError):
        modified_bases_forward = None

    return modified_bases_forward or {}


def count_passing_brdu_calls(
    read: pysam.AlignedSegment,
    minimum_ml: int,
) -> tuple[int, list[str], bool]:
    """
    Count passing BrdU calls on one read.

    This intentionally mirrors the Modkit selection rule:

        call_code == "b"

    The canonical-base component is not required to be T because this BAM
    uses N+b rather than T+b.

    Returns:
        passing_call_count
        descriptions of observed modified-base keys
        whether any b modification key was observed
    """
    modified_bases = get_modified_bases(read)

    if not modified_bases:
        return 0, [], False

    passing_calls = 0
    observed_keys: list[str] = []
    observed_brdu_key = False

    for key, calls in modified_bases.items():
        if not isinstance(key, tuple) or len(key) != 3:
            observed_keys.append(
                f"unrecognized_key={key!r}"
            )
            continue

        canonical_base, strand, modification_code = key

        normalized_base = normalize_canonical_base(
            canonical_base
        )

        normalized_code = normalize_modification_code(
            modification_code
        ).lower()

        observed_keys.append(
            f"base={normalized_base},"
            f"strand={strand},"
            f"code={normalized_code},"
            f"raw_key={key!r}"
        )

        # Match Modkit filtering by modification code b.
        #
        # Do not require base == T because this BAM uses N+b.
        if normalized_code != "b":
            continue

        observed_brdu_key = True

        for call in calls:
            if not isinstance(call, tuple) or len(call) < 2:
                continue

            _query_position, ml_quality = call[:2]

            if ml_quality is None:
                continue

            try:
                ml_value = int(ml_quality)
            except (TypeError, ValueError):
                continue

            # PySAM may use -1 when no probability is available.
            if ml_value < 0:
                continue

            if ml_value >= minimum_ml:
                passing_calls += 1

    return passing_calls, observed_keys, observed_brdu_key


def main() -> None:
    if len(sys.argv) != 7:
        raise SystemExit(
            "Expected: input.bam output.bam read_ids.txt "
            "summary.tsv threshold threads"
        )

    input_bam = sys.argv[1]
    output_bam = sys.argv[2]
    read_ids_path = sys.argv[3]
    summary_path = sys.argv[4]
    threshold = float(sys.argv[5])
    threads = max(1, int(sys.argv[6]))

    if not 0.0 <= threshold <= 1.0:
        raise SystemExit(
            "Threshold must be between 0 and 1."
        )

    # Match the threshold conversion used by the existing workflow:
    #
    #   0.50 -> 128
    #   0.60 -> 153

    minimum_ml = math.ceil(
        threshold * 255
    )

    total_bam_records = 0
    primary_mapped_reads = 0

    reads_with_mm_tag = 0
    reads_with_ml_tag = 0
    reads_with_both_tags = 0

    reads_with_parsed_modifications = 0
    reads_with_brdu_key = 0

    brdu_positive_reads = 0
    total_passing_brdu_calls = 0

    modification_parse_errors = 0

    observed_key_descriptions: list[str] = []
    maximum_diagnostic_keys = 20

    with pysam.AlignmentFile(
        input_bam,
        "rb",
        threads=threads,
        check_sq=False,
    ) as input_handle:

        with pysam.AlignmentFile(
            output_bam,
            "wb",
            template=input_handle,
            threads=threads,
        ) as output_handle, open(
            read_ids_path,
            "w",
            encoding="utf-8",
        ) as read_ids_handle:

            for read in input_handle.fetch(
                until_eof=True
            ):
                total_bam_records += 1

                # Match samtools view -F 2308:
                #
                #   unmapped       0x4
                #   secondary      0x100
                #   supplementary  0x800

                if (
                    read.is_unmapped
                    or read.is_secondary
                    or read.is_supplementary
                ):
                    continue

                primary_mapped_reads += 1

                if primary_mapped_reads % 100_000 == 0:
                    print(
                        "[PROGRESS] "
                        f"Primary reads processed: "
                        f"{primary_mapped_reads:,}; "
                        f"BrdU-positive reads: "
                        f"{brdu_positive_reads:,}; "
                        f"passing BrdU calls: "
                        f"{total_passing_brdu_calls:,}",
                        flush=True,
                    )

                has_mm_tag = (
                    read.has_tag("MM")
                    or read.has_tag("Mm")
                )

                has_ml_tag = (
                    read.has_tag("ML")
                    or read.has_tag("Ml")
                )

                if has_mm_tag:
                    reads_with_mm_tag += 1

                if has_ml_tag:
                    reads_with_ml_tag += 1

                if has_mm_tag and has_ml_tag:
                    reads_with_both_tags += 1

                try:
                    modified_bases = get_modified_bases(
                        read
                    )
                except Exception:
                    modification_parse_errors += 1
                    continue

                if modified_bases:
                    reads_with_parsed_modifications += 1

                try:
                    (
                        passing_call_count,
                        observed_keys,
                        observed_brdu_key,
                    ) = count_passing_brdu_calls(
                        read=read,
                        minimum_ml=minimum_ml,
                    )
                except Exception:
                    modification_parse_errors += 1
                    continue

                if observed_brdu_key:
                    reads_with_brdu_key += 1

                for key_description in observed_keys:
                    if (
                        key_description
                        not in observed_key_descriptions
                        and len(observed_key_descriptions)
                        < maximum_diagnostic_keys
                    ):
                        observed_key_descriptions.append(
                            key_description
                        )

                if passing_call_count == 0:
                    continue

                output_handle.write(read)

                read_ids_handle.write(
                    f"{read.query_name}\n"
                )

                brdu_positive_reads += 1

                total_passing_brdu_calls += (
                    passing_call_count
                )

    percent_selected = (
        brdu_positive_reads
        / primary_mapped_reads
        * 100.0
        if primary_mapped_reads > 0
        else 0.0
    )

    average_calls_per_total_read = (
        total_passing_brdu_calls
        / primary_mapped_reads
        if primary_mapped_reads > 0
        else 0.0
    )

    average_calls_per_positive_read = (
        total_passing_brdu_calls
        / brdu_positive_reads
        if brdu_positive_reads > 0
        else 0.0
    )

    with open(
        summary_path,
        "w",
        encoding="utf-8",
    ) as summary_handle:
        summary_handle.write(
            "input_bam\t"
            "threshold\t"
            "minimum_ml\t"
            "total_bam_records\t"
            "primary_mapped_reads\t"
            "reads_with_mm_tag\t"
            "reads_with_ml_tag\t"
            "reads_with_both_tags\t"
            "reads_with_parsed_modifications\t"
            "reads_with_brdu_key\t"
            "brdu_positive_reads\t"
            "total_passing_brdu_calls\t"
            "percent_selected\t"
            "average_calls_per_total_read\t"
            "average_calls_per_positive_read\t"
            "modification_parse_errors\n"
        )

        summary_handle.write(
            f"{input_bam}\t"
            f"{threshold:.6f}\t"
            f"{minimum_ml}\t"
            f"{total_bam_records}\t"
            f"{primary_mapped_reads}\t"
            f"{reads_with_mm_tag}\t"
            f"{reads_with_ml_tag}\t"
            f"{reads_with_both_tags}\t"
            f"{reads_with_parsed_modifications}\t"
            f"{reads_with_brdu_key}\t"
            f"{brdu_positive_reads}\t"
            f"{total_passing_brdu_calls}\t"
            f"{percent_selected:.6f}\t"
            f"{average_calls_per_total_read:.6f}\t"
            f"{average_calls_per_positive_read:.6f}\t"
            f"{modification_parse_errors}\n"
        )

    print()
    print("=" * 60)
    print("PySAM BrdU-positive read extraction summary")
    print("=" * 60)

    print(
        f"Input BAM:                         "
        f"{input_bam}"
    )

    print(
        f"Modification encoding:            "
        f"N+b"
    )

    print(
        f"Selected modification code:       "
        f"b"
    )

    print(
        f"Modification threshold:           "
        f"{threshold:.6f}"
    )

    print(
        f"Minimum ML value:                 "
        f"{minimum_ml}"
    )

    print(
        f"Total BAM records:                "
        f"{total_bam_records:,}"
    )

    print(
        f"Primary mapped reads:             "
        f"{primary_mapped_reads:,}"
    )

    print(
        f"Reads with MM tags:               "
        f"{reads_with_mm_tag:,}"
    )

    print(
        f"Reads with ML tags:               "
        f"{reads_with_ml_tag:,}"
    )

    print(
        f"Reads with both MM and ML tags:   "
        f"{reads_with_both_tags:,}"
    )

    print(
        f"Reads with parsed modifications:  "
        f"{reads_with_parsed_modifications:,}"
    )

    print(
        f"Reads with a parsed BrdU b key:   "
        f"{reads_with_brdu_key:,}"
    )

    print(
        f"BrdU-positive reads written:      "
        f"{brdu_positive_reads:,}"
    )

    print(
        f"Total passing BrdU calls:         "
        f"{total_passing_brdu_calls:,}"
    )

    print(
        f"Percent of reads selected:        "
        f"{percent_selected:.6f}%"
    )

    print(
        f"Average calls per total read:     "
        f"{average_calls_per_total_read:.6f}"
    )

    print(
        f"Average calls per positive read:  "
        f"{average_calls_per_positive_read:.6f}"
    )

    print(
        f"Modification parse errors:        "
        f"{modification_parse_errors:,}"
    )

    print()
    print("Observed modified-base keys:")

    if observed_key_descriptions:
        for key_description in observed_key_descriptions:
            print(
                f"  {key_description}"
            )
    else:
        print(
            "  No parsed modified-base keys were observed."
        )

    print("=" * 60)

    if brdu_positive_reads == 0:
        raise SystemExit(
            "No primary mapped reads contained a passing "
            "modification-code b call. Inspect the observed "
            "modified-base keys in the log."
        )

    if not os.path.exists(output_bam):
        raise SystemExit(
            "The temporary output BAM was not created."
        )


if __name__ == "__main__":
    main()
PYTHON


    # =========================================================================
    # Coordinate-sort the filtered BAM
    # =========================================================================

    print_info "Coordinate-sorting the filtered BAM."

    samtools sort \
        -@ "$threads" \
        -o "$temp_sorted_bam" \
        "$temp_unsorted_bam"


    # =========================================================================
    # Move completed outputs into data/bam
    # =========================================================================

    mv "$temp_sorted_bam" "$output_bam"
    mv "$temp_read_ids" "$output_read_ids"
    mv "$temp_summary" "$output_summary"


    # =========================================================================
    # Index and validate
    # =========================================================================

    print_info "Indexing the output BAM."

    samtools index \
        -@ "$threads" \
        "$output_bam"

    print_info "Validating the output BAM."

    samtools quickcheck \
        -v \
        "$output_bam"


    # =========================================================================
    # Independent read-count verification
    # =========================================================================

    output_count="$(
        samtools view \
            -@ "$threads" \
            -c \
            -F 2308 \
            "$output_bam"
    )"

    read_id_count="$(
        wc -l < "$output_read_ids"
    )"


    # =========================================================================
    # Final summary
    # =========================================================================

    printf '\n'
    printf '%s\n' '============================================================'
    printf '%s\n' 'PySAM BrdU extraction complete'
    printf '%s\n' '============================================================'
    printf '\n'

    printf 'Output primary mapped reads:   %s\n' \
        "$output_count"

    printf 'Read IDs written:              %s\n' \
        "$read_id_count"

    if [[ "$output_count" -eq "$read_id_count" ]]; then
        printf 'Read-count agreement:          PASS\n'
    else
        printf 'Read-count agreement:          WARNING\n'
    fi

    printf '\n'

    printf 'Output BAM:                    %s\n' \
        "$output_bam"

    printf 'Output BAM index:              %s\n' \
        "${output_bam}.bai"

    printf 'Selected read IDs:             %s\n' \
        "$output_read_ids"

    printf 'Summary TSV:                   %s\n' \
        "$output_summary"

    printf '\n'
    printf '%s\n' '============================================================'

    trap - EXIT
    cleanup
}


# =============================================================================
# Command routing
# =============================================================================

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi


if [[ "${1:-}" == "--run-analysis" ]]; then
    shift

    if [[ "${1:-}" == "--workflow-root" ]]; then
        shift

        if (( $# < 1 )); then
            print_error "Missing path after --workflow-root."
            exit 1
        fi

        # WORKFLOW_ROOT was already parsed before command routing.
        shift
    fi

    if (( $# != 2 )); then
        print_error \
            "Expected: --run-analysis [--workflow-root PATH] BAM THRESHOLD"
        exit 1
    fi

    run_analysis "$1" "$2"
else
    submit_workflow "$@"
fi
