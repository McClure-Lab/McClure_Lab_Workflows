#!/bin/bash

# ==============================================================================
# BrdU READ-SUMMARY DASHBOARD WRAPPER
# ==============================================================================
#
# Purpose
# -------
# This Bash script is the SLURM launcher and compute-job wrapper for:
#
#   plot_brdu_read_summary_dashboard.py
#
# The Python script reads one or more BrdU read-percentage summary logs and
# creates graphical dashboard PNG files that compare the selected samples.
#
# This Bash wrapper is responsible for:
#
#   1. Resolving the workflow root.
#   2. Locating BrdU read-percentage summary logs.
#   3. Prompting the user to select one or more logs.
#   4. Creating or reusing a Python virtual environment.
#   5. Installing Matplotlib when needed.
#   6. Submitting this same Bash script to SLURM.
#   7. Running the Python dashboard generator inside the SLURM allocation.
#   8. Writing SLURM and workflow-specific logs.
#
# Execution modes
# ---------------
#
# Interactive submission mode:
#
#   bash src/utils/plot_brdu_read_summary_dashboard.sh
#
# In this mode, the script:
#
#   - Lists BrdU read-percentage summary logs.
#   - Asks how many logs should be included.
#   - Prompts for each log.
#   - Submits the plotting job to SLURM.
#
# Positional-input mode:
#
#   bash src/utils/plot_brdu_read_summary_dashboard.sh \
#       first_summary.log \
#       second_summary.log
#
# In this mode, the supplied logs are validated and submitted without the
# interactive selection prompts.
#
# Internal compute-job mode:
#
#   --run-job
#
# This mode is normally entered automatically by the submitted SLURM job.
# Users generally do not need to invoke it directly.
#
# Expected input files
# --------------------
#
# Summary logs are expected under:
#
#   workflow_root/logs/read_pct
#
# Valid summary-log names must:
#
#   - Contain:
#
#       .BrdU_read_percentage.
#
#   - End in:
#
#       .log
#
# Example:
#
#   sample.BrdU_read_percentage.threshold_0p6.20260727_101500.log
#
# Python dashboard script
# -----------------------
#
# The wrapper runs:
#
#   workflow_root/src/utils/plot_brdu_read_summary_dashboard.py
#
# The selected summary logs are passed to that script as positional inputs.
#
# Output files
# ------------
#
# Dashboard PNG files are written under:
#
#   workflow_root/results/read_pct_dashboard
#
# Expected dashboard naming:
#
#   brdu_read_summary_dashboard.threshold_<threshold>.png
#
# SLURM logs are written under:
#
#   workflow_root/logs/read_pct_dashboard
#
# Expected files:
#
#   brdu_read_pct_dashboard.<job_id>.slurm.log
#   brdu_read_pct_dashboard.<job_id>.slurm.err
#   brdu_read_pct_dashboard.<job_id>.workflow.log
#
# Python environment
# ------------------
#
# The workflow creates or reuses:
#
#   workflow_root/src/utils/.venv
#
# The environment must contain:
#
#   matplotlib
#
# Matplotlib is installed automatically when it cannot be imported.
#
# ==============================================================================


# ==============================================================================
# SLURM RESOURCE REQUESTS
# ==============================================================================
#
# Requested resources:
#
#   Job name:       brdu_pct_dashboard
#   CPU cores:      1
#   Memory:         4 GB
#   Runtime limit:  30 minutes
#   Partition:      normal
#
#SBATCH --job-name=brdu_pct_dashboard
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --partition=normal


# ==============================================================================
# BASH SAFETY SETTINGS
# ==============================================================================
#
# -e
#   Exit when an unhandled command fails.
#
# -u
#   Treat references to unset variables as errors.
#
# -o pipefail
#   Treat a pipeline as failed when any command within the pipeline fails.
#
set -euo pipefail


# ==============================================================================
# RESOLVE THE WORKFLOW ROOT
# ==============================================================================
#
# Resolution order:
#
#   1. A path supplied through --workflow-root.
#   2. BRDU_DASHBOARD_WORKFLOW_ROOT from the environment.
#   3. Two directories above this script.
#
# The third method assumes this script is stored at:
#
#   workflow_root/src/utils/plot_brdu_read_summary_dashboard.sh
#

WORKFLOW_ROOT_ARG=""

# Search all command-line arguments for an explicit --workflow-root value.
for ((arg_i = 1; arg_i <= $#; arg_i++)); do
    if [[ "${!arg_i}" == "--workflow-root" ]]; then
        next_arg_i=$((arg_i + 1))
        WORKFLOW_ROOT_ARG="${!next_arg_i:-}"
        break
    fi
done

# Select the final workflow root.
if [[ -n "$WORKFLOW_ROOT_ARG" ]]; then
    WORKFLOW_ROOT="$(cd "$WORKFLOW_ROOT_ARG" && pwd)"
elif [[ -n "${BRDU_DASHBOARD_WORKFLOW_ROOT:-}" ]]; then
    WORKFLOW_ROOT="$BRDU_DASHBOARD_WORKFLOW_ROOT"
else
    WORKFLOW_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
fi


# ==============================================================================
# WORKFLOW PATH CONFIGURATION
# ==============================================================================

# Absolute expected path to this Bash wrapper.
#
# The script submits itself back to SLURM using --run-job.
SCRIPT_PATH="$WORKFLOW_ROOT/src/utils/plot_brdu_read_summary_dashboard.sh"

# Python script that parses the summary logs and generates the PNG dashboards.
PYTHON_DASHBOARD_SCRIPT="$WORKFLOW_ROOT/src/utils/plot_brdu_read_summary_dashboard.py"

# Directory containing BrdU read-percentage analysis logs.
LOG_DIR="$WORKFLOW_ROOT/logs/read_pct"

# Directory containing SLURM and workflow logs for dashboard generation.
JOB_LOG_DIR="$WORKFLOW_ROOT/logs/read_pct_dashboard"

# Directory containing generated dashboard PNG files.
RESULT_DIR="$WORKFLOW_ROOT/results/read_pct_dashboard"

# Python virtual environment shared by the dashboard utilities.
VENV_DIR="$WORKFLOW_ROOT/src/utils/.venv"


# ==============================================================================
# FUNCTION: usage
# ==============================================================================
#
# Purpose:
#   Print command-line usage and a brief description of the workflow.
#
# Positional arguments:
#
#   summary_log ...
#       Zero or more BrdU read-percentage summary log files.
#
# When no logs are supplied, the script enters interactive selection mode.
#
usage() {
    echo "Usage: bash src/utils/plot_brdu_read_summary_dashboard.sh [summary_log ...]"
    echo
    echo "summary_log files are resolved by filename under logs/read_pct."
    echo "If arguments are omitted, the script prompts for summary logs before submitting a SLURM job."
    echo "Dashboard PNGs are written to results/read_pct_dashboard."
}


# ==============================================================================
# FUNCTION: print_error
# ==============================================================================
#
# Purpose:
#   Print a consistently formatted error message to standard error.
#
print_error() {
    echo "[ERROR] $*" >&2
}


# ==============================================================================
# FUNCTION: print_info
# ==============================================================================
#
# Purpose:
#   Print a consistently formatted informational message.
#
print_info() {
    echo "[INFO] $*"
}


# ==============================================================================
# FUNCTION: absolute_existing_file
# ==============================================================================
#
# Purpose:
#   Convert a file path into a normalized absolute path.
#
# Arguments:
#
#   $1
#       Existing file path.
#
# Output:
#   Absolute file path written to standard output.
#
absolute_existing_file() {
    local path="$1"
    local dir
    local file

    # Resolve the containing directory.
    dir="$(cd "$(dirname "$path")" && pwd)"

    # Preserve the original filename.
    file="$(basename "$path")"

    printf '%s/%s\n' "$dir" "$file"
}


# ==============================================================================
# FUNCTION: require_existing_dir
# ==============================================================================
#
# Purpose:
#   Verify that a required directory exists and is visible on the current node.
#
# Arguments:
#
#   $1
#       Human-readable directory name.
#
#   $2
#       Directory path.
#
# Behavior:
#   The script exits when the directory is missing.
#
require_existing_dir() {
    local name="$1"
    local path="$2"

    if [[ ! -d "$path" ]]; then
        print_error "$name does not exist: $path"
        print_error "Confirm the project path is mounted on the compute node."
        exit 1
    fi
}


# ==============================================================================
# FUNCTION: is_summary_log_name
# ==============================================================================
#
# Purpose:
#   Determine whether a filename matches the expected BrdU read-percentage log
#   naming convention.
#
# Requirements:
#
#   - The name contains:
#
#       .BrdU_read_percentage.
#
#   - The name ends in:
#
#       .log
#
# Arguments:
#
#   $1
#       Filename or path component to test.
#
# Return status:
#
#   0 when the name is valid.
#   1 when the name is invalid.
#
is_summary_log_name() {
    local name="$1"

    [[ "$name" == *".BrdU_read_percentage."* && "$name" == *.log ]]
}


# ==============================================================================
# FUNCTION: resolve_summary_log_file
# ==============================================================================
#
# Purpose:
#   Resolve a user-supplied summary-log filename or path.
#
# Resolution order:
#
#   1. Accept the supplied path directly when it exists and has a valid name.
#   2. Search LOG_DIR using the basename of the supplied value.
#
# Arguments:
#
#   $1
#       Summary-log filename or path.
#
# Output:
#   Absolute path when the log is resolved successfully.
#
# Return status:
#
#   0 when the log exists and matches the required naming convention.
#   1 when it cannot be resolved or is not a valid summary log.
#
resolve_summary_log_file() {
    local log_input="$1"
    local log_path

    # Accept a complete or relative file path.
    if [[ -f "$log_input" ]]; then
        if is_summary_log_name "$(basename "$log_input")"; then
            absolute_existing_file "$log_input"
            return 0
        fi

        return 1
    fi

    # Otherwise, search by basename under logs/read_pct.
    log_path="$LOG_DIR/$(basename "$log_input")"

    if [[ -f "$log_path" ]]; then
        if is_summary_log_name "$(basename "$log_path")"; then
            absolute_existing_file "$log_path"
            return 0
        fi
    fi

    return 1
}


# ==============================================================================
# FUNCTION: load_required_modules
# ==============================================================================
#
# Purpose:
#   Load Python through the cluster module system.
#
# Behavior:
#
#   - Sources the cluster module initialization script when available.
#   - Attempts to load Python 3.13.7.
#   - Falls back to the generic Python module.
#   - Falls back to Python already available through PATH when modules are
#     unavailable.
#   - Verifies that python3 exists before continuing.
#
load_required_modules() {
    # Initialize the module command when the cluster provides this file.
    if [[ -f /etc/profile.d/modules.sh ]]; then
        # shellcheck source=/dev/null
        source /etc/profile.d/modules.sh
    fi

    if command -v module >/dev/null 2>&1; then
        module load python/3.13.7 || module load python
    else
        echo "[WARN] Environment modules are not available in this shell."
        echo "[WARN] Continuing with python from PATH."
    fi

    # Confirm that Python is available after module loading.
    if ! command -v python3 >/dev/null 2>&1; then
        print_error "python3 was not found."
        exit 1
    fi
}


# ==============================================================================
# FUNCTION: prepare_python_environment
# ==============================================================================
#
# Purpose:
#   Create or reuse the dashboard Python virtual environment and ensure that
#   Matplotlib is installed.
#
# Environment:
#
#   workflow_root/src/utils/.venv
#
# Steps:
#
#   1. Create the virtual environment when its Python executable is missing.
#   2. Verify the environment Python exists.
#   3. Install pip through ensurepip when necessary.
#   4. Install Matplotlib when it cannot be imported.
#   5. Verify Matplotlib can be imported.
#
prepare_python_environment() {
    # Create the environment only when its Python executable is unavailable.
    if [[ ! -x "$VENV_DIR/bin/python" ]]; then
        print_info "Creating Python virtual environment: $VENV_DIR"
        python3 -m venv "$VENV_DIR"
    else
        print_info "Reusing Python virtual environment: $VENV_DIR"
    fi

    # Confirm that environment creation succeeded.
    if [[ ! -x "$VENV_DIR/bin/python" ]]; then
        print_error "Virtual-environment Python was not created: $VENV_DIR/bin/python"
        exit 1
    fi

    # Install pip when it is not available inside the environment.
    if ! "$VENV_DIR/bin/python" -m pip --version >/dev/null 2>&1; then
        print_info "Installing pip in the Python virtual environment."

        "$VENV_DIR/bin/python" \
            -m ensurepip \
            --upgrade
    fi

    # Install Matplotlib only when it cannot already be imported.
    if ! "$VENV_DIR/bin/python" -c "import matplotlib" >/dev/null 2>&1; then
        print_info "Installing matplotlib in the Python virtual environment."

        "$VENV_DIR/bin/python" \
            -m pip install \
            matplotlib
    else
        print_info "Required Python packages are already installed."
    fi

    # Perform a final import check.
    if ! "$VENV_DIR/bin/python" -c "import matplotlib" >/dev/null 2>&1; then
        print_error "matplotlib is not available in $VENV_DIR."
        exit 1
    fi
}


# ==============================================================================
# FUNCTION: run_job
# ==============================================================================
#
# Purpose:
#   Execute dashboard generation inside the submitted SLURM allocation.
#
# Internal arguments:
#
#   --workflow-root PATH
#       Already resolved near the beginning of this script. It is consumed here
#       so it is not treated as an unknown argument.
#
#   --output-dir PATH
#       Accepted for compatibility, but the wrapper uses RESULT_DIR.
#
#   --log-dir PATH
#       Accepted for compatibility, but the wrapper uses LOG_DIR.
#
#   --summary-log PATH
#       May be supplied multiple times. Each occurrence adds one input log.
#
# Major stages:
#
#   1. Parse internal job arguments.
#   2. Validate workflow directories and Python script.
#   3. Create result and job-log directories.
#   4. Configure a combined workflow log.
#   5. Load Python and prepare the virtual environment.
#   6. Run the Python dashboard generator.
#
run_job() {
    # Array containing all selected summary-log paths.
    local summary_logs=()

    # Combined workflow-log path.
    local workflow_log

    # Remove the leading --run-job argument.
    shift

    # --------------------------------------------------------------------------
    # Parse internal job arguments
    # --------------------------------------------------------------------------

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --workflow-root)
                shift 2
                ;;
            --output-dir)
                shift 2
                ;;
            --log-dir)
                shift 2
                ;;
            --summary-log)
                summary_logs+=("$2")
                shift 2
                ;;
            *)
                print_error "Unknown job argument: $1"
                exit 1
                ;;
        esac
    done

    # --------------------------------------------------------------------------
    # Validate required directories and output locations
    # --------------------------------------------------------------------------

    require_existing_dir "Workflow root" "$WORKFLOW_ROOT"
    require_existing_dir "read_pct log directory" "$LOG_DIR"

    mkdir -p \
        "$JOB_LOG_DIR" \
        "$RESULT_DIR"

    # Require the Python plotting implementation.
    if [[ ! -f "$PYTHON_DASHBOARD_SCRIPT" ]]; then
        print_error "Python dashboard script not found: $PYTHON_DASHBOARD_SCRIPT"
        exit 1
    fi

    # At least one selected summary log is required.
    if [[ ${#summary_logs[@]} -eq 0 ]]; then
        print_error "No summary log files were provided to the SLURM job."
        exit 1
    fi

    # --------------------------------------------------------------------------
    # Configure the combined workflow log
    # --------------------------------------------------------------------------
    #
    # All subsequent standard output and standard error are written both to the
    # active SLURM streams and the workflow-specific log through tee.
    #

    workflow_log="$JOB_LOG_DIR/brdu_read_pct_dashboard.${SLURM_JOB_ID:-manual}.workflow.log"

    exec > >(tee -a "$workflow_log") 2>&1

    print_info "BrdU read-summary dashboard started: $(date)"
    print_info "Workflow log: $workflow_log"

    # --------------------------------------------------------------------------
    # Prepare Python and generate the dashboard
    # --------------------------------------------------------------------------

    load_required_modules
    prepare_python_environment

    # Run the Python implementation using the prepared virtual environment.
    #
    # --run-job:
    #   Prevents the Python script from submitting another SLURM job.
    #
    # The summary logs are passed as positional arguments.
    #
    "$VENV_DIR/bin/python" \
        "$PYTHON_DASHBOARD_SCRIPT" \
        --run-job \
        --workflow-root "$WORKFLOW_ROOT" \
        --log-dir "$LOG_DIR" \
        --output-dir "$RESULT_DIR" \
        "${summary_logs[@]}"

    print_info "BrdU read-summary dashboard finished: $(date)"
}


# ==============================================================================
# FUNCTION: submit_workflow
# ==============================================================================
#
# Purpose:
#   Run the interactive selection stage and submit the dashboard job to SLURM.
#
# Positional arguments:
#
#   $@
#       Optional summary-log filenames or paths.
#
# Behavior when no logs are supplied:
#
#   - Lists matching logs under LOG_DIR.
#   - Prompts for the number of logs to include.
#   - Prompts for each individual log.
#
# Behavior when logs are supplied:
#
#   - Validates each supplied log.
#   - Resolves each log to an absolute path.
#
# All selected logs are converted to repeated --summary-log arguments before
# the job is submitted.
#
submit_workflow() {
    # Logs supplied directly as positional arguments.
    local summary_log_inputs=("$@")

    # Resolved absolute summary-log paths.
    local summary_logs=()

    # Available filenames displayed during interactive selection.
    local available_logs=()

    # Number of logs requested by the user.
    local log_count

    # Current interactive selection.
    local selection

    # Resolved path for the current selection.
    local resolved_log

    # Submitted SLURM job ID.
    local job_id
    local log_job_id

    # Repeated --summary-log arguments passed to run_job().
    local sbatch_args=()

    # --------------------------------------------------------------------------
    # Validate workflow directories
    # --------------------------------------------------------------------------

    require_existing_dir "Workflow root" "$WORKFLOW_ROOT"
    require_existing_dir "read_pct log directory" "$LOG_DIR"

    mkdir -p \
        "$JOB_LOG_DIR" \
        "$RESULT_DIR"

    # --------------------------------------------------------------------------
    # Interactively select logs when none were supplied
    # --------------------------------------------------------------------------

    if [[ ${#summary_log_inputs[@]} -eq 0 ]]; then
        # List matching summary logs directly under LOG_DIR.
        mapfile -t available_logs < <(
            find "$LOG_DIR" \
                -maxdepth 1 \
                -type f \
                -name "*.BrdU_read_percentage.*.log" \
                -printf "%f\n" |
            sort
        )

        # Stop when no valid summary logs are available.
        if [[ ${#available_logs[@]} -eq 0 ]]; then
            print_error "No BrdU read-percentage summary .log files were found in: $LOG_DIR"
            exit 1
        fi

        # Display a one-based numbered menu.
        echo "[INFO] Available BrdU read-percentage summary .log files in $LOG_DIR:"

        for index in "${!available_logs[@]}"; do
            printf "  %3d) %s\n" \
                "$((index + 1))" \
                "${available_logs[index]}"
        done

        echo

        # Ask how many logs should be included.
        while true; do
            read -r -p "How many summary .log files should be included? " log_count

            if [[ "$log_count" =~ ^[0-9]+$ && "$log_count" -ge 1 ]]; then
                break
            fi

            echo "Please enter a positive whole number."
        done

        # Prompt separately for each requested summary log.
        for ((log_i = 1; log_i <= log_count; log_i++)); do
            while true; do
                read -r -p "Select summary .log file ${log_i}/${log_count} by number or filename: " selection

                if [[ -z "$selection" ]]; then
                    echo "Please provide a selection."
                    continue
                fi

                # Numeric input is interpreted as a one-based menu selection.
                if [[ "$selection" =~ ^[0-9]+$ ]]; then
                    if (( selection >= 1 && selection <= ${#available_logs[@]} )); then
                        selection="${available_logs[selection - 1]}"
                    else
                        echo "Selection must be between 1 and ${#available_logs[@]}."
                        continue
                    fi
                fi

                # Resolve and validate the selected log.
                if resolved_log="$(resolve_summary_log_file "$selection")"; then
                    summary_logs+=("$resolved_log")
                    break
                fi

                echo "Summary .log file not found or not a BrdU read-percentage summary log: $selection"
            done
        done
    else
        # ----------------------------------------------------------------------
        # Validate logs supplied as positional arguments
        # ----------------------------------------------------------------------

        for summary_log_input in "${summary_log_inputs[@]}"; do
            if resolved_log="$(resolve_summary_log_file "$summary_log_input")"; then
                summary_logs+=("$resolved_log")
            else
                print_error "Summary .log file not found or invalid: $summary_log_input"
                exit 1
            fi
        done
    fi

    # --------------------------------------------------------------------------
    # Confirm SLURM is available
    # --------------------------------------------------------------------------

    if ! command -v sbatch >/dev/null 2>&1; then
        print_error "The sbatch command is unavailable."
        print_error "Run this script on the cluster login node."
        exit 1
    fi

    # --------------------------------------------------------------------------
    # Construct repeated --summary-log job arguments
    # --------------------------------------------------------------------------

    for summary_log in "${summary_logs[@]}"; do
        sbatch_args+=(
            --summary-log
            "$summary_log"
        )
    done

    # --------------------------------------------------------------------------
    # Display the final submission configuration
    # --------------------------------------------------------------------------

    echo

    print_info "Workflow root: $WORKFLOW_ROOT"
    print_info "Output directory: $RESULT_DIR"
    print_info "Selected summary log file(s):"

    for summary_log in "${summary_logs[@]}"; do
        print_info "  $summary_log"
    done

    echo
    echo "Submitting BrdU read-summary dashboard to SLURM..."

    # --------------------------------------------------------------------------
    # Submit this script back to SLURM
    # --------------------------------------------------------------------------
    #
    # --parsable
    #   Return a machine-readable job ID.
    #
    # --chdir
    #   Run the compute job from the workflow root.
    #
    # --output / --error
    #   Write SLURM output and error logs under JOB_LOG_DIR.
    #
    # --export
    #   Export the current environment and the resolved workflow root.
    #
    # --run-job
    #   Select the compute-job path in this Bash wrapper.
    #

    job_id="$(sbatch --parsable \
        --job-name="brdu_read_pct_dashboard" \
        --chdir="$WORKFLOW_ROOT" \
        --output="$JOB_LOG_DIR/brdu_read_pct_dashboard.%j.slurm.log" \
        --error="$JOB_LOG_DIR/brdu_read_pct_dashboard.%j.slurm.err" \
        --export=ALL,BRDU_DASHBOARD_WORKFLOW_ROOT="$WORKFLOW_ROOT" \
        "$SCRIPT_PATH" \
        --run-job \
        --workflow-root "$WORKFLOW_ROOT" \
        "${sbatch_args[@]}")"

    # Some SLURM systems append cluster information after a semicolon.
    log_job_id="${job_id%%;*}"

    # --------------------------------------------------------------------------
    # Report the submitted job and expected output paths
    # --------------------------------------------------------------------------

    print_info "Submitted SLURM job: $job_id"

    print_info "Expected SLURM log: $JOB_LOG_DIR/brdu_read_pct_dashboard.${log_job_id}.slurm.log"

    print_info "Expected SLURM err: $JOB_LOG_DIR/brdu_read_pct_dashboard.${log_job_id}.slurm.err"

    print_info "Expected workflow log: $JOB_LOG_DIR/brdu_read_pct_dashboard.${log_job_id}.workflow.log"

    print_info "Expected dashboard output:"

    print_info "  $RESULT_DIR/brdu_read_summary_dashboard.threshold_<threshold>.png"
}


# ==============================================================================
# MAIN SCRIPT ROUTING
# ==============================================================================
#
# -h or --help
#   Print usage information and exit.
#
# --run-job
#   Execute dashboard generation inside the submitted SLURM job.
#
# Any other invocation
#   Enter interactive or positional-input submission mode.
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