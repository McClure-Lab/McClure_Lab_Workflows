#!/usr/bin/env python3
#SBATCH --job-name=brdu_pct_dashboard
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --partition=normal
"""
Generate BrdU read-summary dashboard plots from read-percentage summary files.

Purpose
-------
This script creates one four-panel PNG dashboard for each BrdU modification
threshold represented in the selected summary files.

It accepts human-readable ``.log`` reports and machine-readable ``.tsv`` files
created by ``calculate_brdu_read_percentage.sh``. Records are parsed, cleaned,
deduplicated, grouped by threshold, sorted into a consistent biological sample
order, and rendered with Matplotlib.

Execution modes
---------------
Interactive launcher mode
    Running the script from a terminal prompts the user to choose one or more
    BrdU read-percentage summary files and then submits this same Python script
    to SLURM with ``--run-job``.

Submitted plotting mode
    ``--run-job`` performs parsing and plotting directly. Noninteractive SLURM
    execution is also detected automatically through ``SLURM_JOB_ID`` and
    ``sys.stdin.isatty()``.

Python environment
------------------
The script creates or reuses:

    src/utils/.venv

It ensures pip and Matplotlib are available, then restarts itself with the
virtual-environment Python interpreter when necessary.

Default inputs
--------------
Interactive selection scans:

    {workflow_root}/logs/read_pct

Supported files must:

* End in ``.log`` or ``.tsv``.
* Contain ``.BrdU_read_percentage.`` in the filename, case-insensitively.

Files and directories may also be supplied as positional arguments.

Default outputs
---------------
Dashboard PNG files are written under:

    {workflow_root}/results/read_pct_dashboard

One image is generated per threshold:

    brdu_read_summary_dashboard.threshold_<threshold>.png

Dashboard panels
----------------
Panel 1
    Stacked bars showing BrdU-positive reads and remaining reads.

Panel 2
    Percentage of reads containing at least one passing BrdU call.

Panel 3
    Average number of passing BrdU calls per total read.

Panel 4
    A summary table containing the values shown in the plots.

Important metric distinction
----------------------------
``percent_reads_with_brdu`` measures the percentage of reads containing at
least one passing BrdU call.

``avg_calls_per_total_read`` divides the total number of passing BrdU calls by
the total read count. One read may contribute multiple BrdU calls.

SLURM resources
---------------
The script requests one CPU, 4 GB RAM, a 30-minute runtime limit, and the
``normal`` partition.
"""

import argparse
import csv
import datetime as dt
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict


plt = None
FuncFormatter = None
SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
VENV_DIR = os.path.join(SCRIPT_DIR, ".venv")
REQUIRED_PYTHON_PACKAGES = {
    "matplotlib": "matplotlib",
}


# =============================================================================
# Helpers
# =============================================================================


def thousands_formatter(x, pos):
    """
    Format a Matplotlib tick value with comma-separated thousands.

    Parameters
    ----------
    x
        Numeric tick value supplied by Matplotlib.

    pos
        Tick position supplied by Matplotlib. It is required by the formatter
        API but is not otherwise used.

    Returns
    -------
    str
        Integer-formatted tick label, such as ``250,000``.
    """
    return f"{int(x):,}"


def require_matplotlib():
    """
    Import Matplotlib lazily and configure the noninteractive Agg backend.

    Matplotlib is delayed until plotting is required so launcher operations,
    argument parsing, and environment preparation can run before the package is
    imported. The Agg backend allows PNG creation on headless compute nodes.

    The imported pyplot module and FuncFormatter class are cached in module-
    level globals.
    """
    global plt
    global FuncFormatter

    if plt is not None and FuncFormatter is not None:
        return

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as pyplot
        from matplotlib.ticker import FuncFormatter as MatplotlibFuncFormatter
    except ModuleNotFoundError as exc:
        if exc.name != "matplotlib":
            raise

        print(
            "[ERROR] matplotlib is required to generate BrdU dashboard plots.\n"
            "Install it in your active Python environment, for example:\n"
            "  python3 -m pip install matplotlib",
            file=sys.stderr,
        )
        sys.exit(1)

    plt = pyplot
    FuncFormatter = MatplotlibFuncFormatter


def clean_int(value):
    """
    Convert a count containing optional commas and whitespace into an integer.

    Examples
    --------
    ``"228,797"`` becomes ``228797``.
    """
    return int(str(value).replace(",", "").strip())


def workflow_root_from_script():
    """
    Infer the workflow root from this script's expected ``src/utils`` location.
    """
    return os.path.abspath(
        os.path.join(
            os.path.dirname(__file__),
            "..",
            "..",
        )
    )


def absolute_existing_file(path):
    """
    Return an absolute path for a selected input or script path.

    File existence is validated by the calling selection or normalization
    routine.
    """
    return os.path.abspath(path)


def running_as_noninteractive_slurm_job():
    """
    Return True for a submitted SLURM process that cannot answer prompts.

    Interactive SLURM shells retain a terminal and therefore remain eligible
    for the normal interactive launcher workflow.
    """
    return bool(os.environ.get("SLURM_JOB_ID")) and not sys.stdin.isatty()


def venv_python_path():
    """
    Return the expected Python executable inside the local virtual environment.
    """
    return os.path.join(
        VENV_DIR,
        "bin",
        "python",
    )


def current_python_is_venv_python():
    """
    Determine whether the current interpreter is the prepared environment's
    Python executable.
    """
    expected = os.path.realpath(venv_python_path())
    current = os.path.realpath(sys.executable)

    return current == expected


def run_python_command(
    python_cmd,
    command_args,
    error_message,
):
    """
    Run a Python subprocess and convert command failure into a concise exit
    message.

    Parameters
    ----------
    python_cmd
        Python executable used to run the command.

    command_args
        Arguments passed after the Python executable.

    error_message
        Message shown when the command returns a nonzero exit status.
    """
    try:
        subprocess.run(
            [
                python_cmd,
                *command_args,
            ],
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        raise SystemExit(
            f"[ERROR] {error_message}"
        ) from exc


def package_is_importable(
    python_cmd,
    import_name,
):
    """
    Test whether a module can be imported by a specific Python interpreter.

    Parameters
    ----------
    python_cmd
        Python executable to test.

    import_name
        Module name used by the import statement.

    Returns
    -------
    bool
        True when the import succeeds.
    """
    result = subprocess.run(
        [
            python_cmd,
            "-c",
            f"import {import_name}",
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )

    return result.returncode == 0


def prepare_python_environment():
    """
    Create or reuse ``src/utils/.venv`` and install required packages.

    The environment is created with the interpreter currently running this
    script. pip is installed with ensurepip when necessary. Only missing
    packages are installed.

    Returns
    -------
    str
        Path to the virtual-environment Python executable.
    """
    python_cmd = venv_python_path()

    if not os.path.exists(python_cmd):
        print(
            f"[INFO] Creating Python virtual environment: "
            f"{VENV_DIR}"
        )

        run_python_command(
            sys.executable,
            [
                "-m",
                "venv",
                VENV_DIR,
            ],
            "Could not create the Python virtual environment.",
        )
    else:
        print(
            f"[INFO] Reusing Python virtual environment: "
            f"{VENV_DIR}"
        )

    if not os.path.exists(python_cmd):
        raise SystemExit(
            "[ERROR] Virtual-environment Python was not created: "
            f"{python_cmd}"
        )

    if not package_is_importable(
        python_cmd,
        "pip",
    ):
        print(
            "[INFO] Installing pip in the virtual environment "
            "with ensurepip."
        )

        run_python_command(
            python_cmd,
            [
                "-m",
                "ensurepip",
                "--upgrade",
            ],
            "Could not install pip in the Python virtual environment.",
        )

    missing_packages = [
        package_name
        for import_name, package_name in REQUIRED_PYTHON_PACKAGES.items()
        if not package_is_importable(
            python_cmd,
            import_name,
        )
    ]

    if missing_packages:
        print("[INFO] Installing missing Python package(s):")

        for package_name in missing_packages:
            print(f"[INFO]   {package_name}")

        run_python_command(
            python_cmd,
            [
                "-m",
                "pip",
                "install",
                *missing_packages,
            ],
            "Could not install required Python package(s) "
            "in the virtual environment.",
        )
    else:
        print(
            "[INFO] Required Python packages are already installed."
        )

    return python_cmd


def ensure_running_inside_prepared_environment():
    """
    Prepare dependencies and restart this script with the environment Python.

    ``os.execv`` replaces the current process rather than launching a nested
    long-lived process. Existing command-line arguments are preserved.
    """
    python_cmd = prepare_python_environment()

    if current_python_is_venv_python():
        return

    print(
        "[INFO] Restarting with virtual-environment Python: "
        f"{python_cmd}"
    )

    os.execv(
        python_cmd,
        [
            python_cmd,
            absolute_existing_file(__file__),
            *sys.argv[1:],
        ],
    )


def clean_sample_name(input_bam):
    """
    Convert an input BAM filename into a concise biological sample label.

    Known control and cell-cycle names are normalized. Unknown names retain
    their basename with underscores replaced by spaces.

    Parameters
    ----------
    input_bam
        BAM path or filename recorded in the summary file.

    Returns
    -------
    str
        Cleaned label used in the dashboard.
    """
    base = os.path.basename(
        input_bam.strip()
    )

    suffixes = [
        ".sorted.indexed.BrdU.detect.bam",
        ".BrdU.detect.bam",
        ".sorted.bam",
        ".bam",
    ]

    for suffix in suffixes:
        if base.endswith(suffix):
            base = base[: -len(suffix)]
            break

    lower = base.lower()

    if "negative" in lower:
        return "Negative control"

    if (
        "s_phase" in lower
        or "s-phase" in lower
        or "s phase" in lower
    ):
        return "S phase"

    if "mitosis" in lower and "wt" in lower:
        return "Mitosis WT"

    if "mitosis" in lower and "mms" in lower:
        return "Mitosis MMS"

    return base.replace("_", " ")


def sample_sort_key(label):
    """
    Return a stable biological display-order key for dashboard samples.

    Known samples are ordered as:

    1. Negative control
    2. S phase
    3. Mitosis WT
    4. Mitosis MMS

    Unrecognized labels follow the named samples and are sorted alphabetically.
    """
    order = {
        "Negative control": 0,
        "S phase": 1,
        "Mitosis WT": 2,
        "Mitosis MMS": 3,
    }

    return (
        order.get(label, 999),
        label,
    )


def threshold_label(threshold):
    """
    Convert a numeric threshold into a compact filename-safe label.

    Example
    -------
    ``0.6`` becomes ``0p6``.
    """
    s = str(threshold)

    if "." in s:
        s = s.rstrip("0").rstrip(".")

    return s.replace(".", "p")


def find_value(
    pattern,
    text,
    flags=0,
    default=None,
):
    """
    Return the first regular-expression capture group or a supplied default.

    Parameters
    ----------
    pattern
        Regular-expression pattern containing a capture group.

    text
        Text searched by the regular expression.

    flags
        Optional regular-expression flags.

    default
        Value returned when no match is found.
    """
    match = re.search(
        pattern,
        text,
        flags,
    )

    if match:
        return match.group(1).strip()

    return default


def is_summary_file(file_path):
    """
    Return True when a path has the expected read-percentage name and extension.

    Supported summary filenames must:

    - Contain ``.BrdU_read_percentage.`` case-insensitively.
    - End in ``.log`` or ``.tsv``.
    """
    name = os.path.basename(file_path)
    lower_name = name.lower()

    is_read_pct_output = (
        ".brdu_read_percentage." in lower_name
    )

    is_supported_ext = lower_name.endswith(
        (
            ".log",
            ".tsv",
        )
    )

    return (
        is_read_pct_output
        and is_supported_ext
    )


def list_summary_files(default_log_dir):
    """
    List supported BrdU summary files directly inside the default log directory.

    Directory scanning is nonrecursive.

    Parameters
    ----------
    default_log_dir
        Directory containing the summary files.

    Returns
    -------
    list[str]
        Sorted list of supported file paths.
    """
    if not os.path.isdir(default_log_dir):
        return []

    return [
        os.path.join(
            default_log_dir,
            name,
        )
        for name in sorted(
            os.listdir(default_log_dir)
        )
        if is_summary_file(name)
    ]


def prompt_for_positive_int(prompt):
    """
    Prompt repeatedly until the user enters an integer greater than zero.
    """
    while True:
        value = input(prompt).strip()

        if value.isdigit() and int(value) > 0:
            return int(value)

        print(
            "Please enter a positive whole number."
        )


def prompt_for_summary_file(
    default_log_dir,
    selected_files=None,
):
    """
    Prompt for one summary file using a number, filename, or complete path.

    Previously selected paths are hidden from the numbered menu to prevent
    duplicate interactive selections.

    Parameters
    ----------
    default_log_dir
        Default directory searched for summary files.

    selected_files
        Optional collection of files already selected.

    Returns
    -------
    str
        Selected summary-file path.
    """
    selected_files = set(
        selected_files or []
    )

    summary_files = [
        file_path
        for file_path in list_summary_files(
            default_log_dir
        )
        if absolute_existing_file(file_path)
        not in selected_files
    ]

    if summary_files:
        print(
            "[INFO] Available BrdU read-percentage summary files "
            f"in {default_log_dir}:"
        )

        for index, summary_file in enumerate(
            summary_files,
            start=1,
        ):
            print(
                f"  {index:3d}) "
                f"{os.path.basename(summary_file)}"
            )

        print()
    else:
        print(
            "[INFO] No BrdU read-percentage summary .log/.tsv "
            f"files found in {default_log_dir}."
        )

    while True:
        selection = input(
            "Enter a summary .log/.tsv file path "
            "or selection number: "
        ).strip()

        if not selection:
            print(
                "Please enter a summary .log/.tsv file path "
                "or selection number."
            )
            continue

        if selection.isdigit() and summary_files:
            selected_index = int(selection)

            if 1 <= selected_index <= len(summary_files):
                return summary_files[
                    selected_index - 1
                ]

            print(
                "Selection must be between 1 and "
                f"{len(summary_files)}."
            )
            continue

        expanded_selection = os.path.expanduser(
            selection
        )

        if os.path.isfile(expanded_selection):
            if is_summary_file(expanded_selection):
                return expanded_selection

            print(
                "Please choose a BrdU read-percentage "
                "summary .log or .tsv file."
            )
            continue

        candidate = os.path.join(
            default_log_dir,
            expanded_selection,
        )

        if os.path.isfile(candidate):
            if is_summary_file(candidate):
                return candidate

            print(
                "Please choose a BrdU read-percentage "
                "summary .log or .tsv file."
            )
            continue

        print(
            f"File not found: {selection}"
        )


def prompt_for_summary_files(default_log_dir):
    """
    Prompt for the number of files and collect that many unique summary paths.

    Parameters
    ----------
    default_log_dir
        Default directory containing summary files.

    Returns
    -------
    list[str]
        Absolute paths selected by the user.
    """
    file_count = prompt_for_positive_int(
        "How many BrdU read-percentage summary "
        ".log/.tsv files should be included? "
    )

    selected_files = []

    for file_index in range(
        1,
        file_count + 1,
    ):
        print(
            f"\n[INFO] Select summary file "
            f"{file_index} of {file_count}."
        )

        selected_file = prompt_for_summary_file(
            default_log_dir,
            selected_files,
        )

        selected_files.append(
            absolute_existing_file(
                selected_file
            )
        )

    return selected_files


def gather_input_files(
    paths,
    default_log_dir,
):
    """
    Expand files and directories into a list of supported summary files.

    Missing paths and unsupported files are skipped with warnings. Directory
    scanning is nonrecursive.

    Resolution order for each path
    ------------------------------
    1. Existing file as entered.
    2. Filename relative to the default log directory.
    3. Existing directory scanned for supported files.

    Parameters
    ----------
    paths
        File and directory paths supplied by the user.

    default_log_dir
        Default directory used to resolve filenames.

    Returns
    -------
    list[str]
        Supported input file paths.
    """
    files = []

    for path in paths:
        expanded_path = os.path.expanduser(path)

        if os.path.isfile(expanded_path):
            if is_summary_file(expanded_path):
                files.append(expanded_path)
            else:
                print(
                    "[WARNING] Skipping non-summary file: "
                    f"{path}"
                )

            continue

        candidate = os.path.join(
            default_log_dir,
            expanded_path,
        )

        if os.path.isfile(candidate):
            if is_summary_file(candidate):
                files.append(candidate)
            else:
                print(
                    "[WARNING] Skipping non-summary file: "
                    f"{path}"
                )

            continue

        if os.path.isdir(expanded_path):
            for name in sorted(
                os.listdir(expanded_path)
            ):
                candidate = os.path.join(
                    expanded_path,
                    name,
                )

                if is_summary_file(candidate):
                    files.append(
                        os.path.join(
                            expanded_path,
                            name,
                        )
                    )

            continue

        print(
            f"[WARNING] Skipping missing path: {path}"
        )

    return files


def normalize_existing_inputs(
    paths,
    default_log_dir,
):
    """
    Resolve supported inputs to absolute paths and require at least one file.

    Parameters
    ----------
    paths
        Files and directories requested by the user.

    default_log_dir
        Default directory used to resolve filenames.

    Returns
    -------
    list[str]
        Absolute summary-file paths.
    """
    input_files = gather_input_files(
        paths,
        default_log_dir,
    )

    if not input_files:
        raise SystemExit(
            "[ERROR] No BrdU read-percentage summary "
            ".log/.tsv files found."
        )

    return [
        absolute_existing_file(path)
        for path in input_files
    ]


# =============================================================================
# Parsing
# =============================================================================


def parse_log_text(
    text,
    source_name="<memory>",
):
    """
    Parse one text blob containing one or more human-readable analysis blocks.

    Each detected block is converted into the common internal record format.
    Incomplete blocks are skipped.

    Parameters
    ----------
    text
        Full text of the human-readable report or combined log.

    source_name
        Source path retained in each output record.

    Returns
    -------
    list[dict]
        Parsed summary records.
    """
    records = []

    blocks = re.split(
        r"BrdU read-percentage analysis",
        text,
    )

    for block in blocks[1:]:
        input_bam = find_value(
            r"Input BAM:\s*(.+)",
            block,
        )

        threshold = find_value(
            r"BrdU modification threshold:\s*([0-9.]+)",
            block,
        )

        total_reads = find_value(
            r"Total reads:\s*([\d,]+)",
            block,
        )

        brdu_reads = find_value(
            r"Reads with at least one BrdU call:\s*([\d,]+)",
            block,
        )

        total_calls = find_value(
            r"Total passing BrdU calls:\s*([\d,]+)",
            block,
        )

        if not all(
            [
                input_bam,
                threshold,
                total_reads,
                brdu_reads,
                total_calls,
            ]
        ):
            continue

        total_reads = clean_int(total_reads)
        brdu_reads = clean_int(brdu_reads)
        total_calls = clean_int(total_calls)
        threshold = float(threshold)

        percent_reads_with_brdu = (
            (brdu_reads / total_reads) * 100
            if total_reads > 0
            else 0.0
        )

        avg_calls_per_total_read = (
            total_calls / total_reads
            if total_reads > 0
            else 0.0
        )

        avg_calls_per_positive_read = (
            total_calls / brdu_reads
            if brdu_reads > 0
            else 0.0
        )

        records.append(
            {
                "source_name": source_name,
                "input_bam": input_bam,
                "sample": clean_sample_name(
                    input_bam
                ),
                "threshold": threshold,
                "total_reads": total_reads,
                "brdu_positive_reads": brdu_reads,
                "reads_without_brdu": (
                    total_reads - brdu_reads
                ),
                "total_passing_brdu_calls": (
                    total_calls
                ),
                "percent_reads_with_brdu": (
                    percent_reads_with_brdu
                ),
                "avg_calls_per_total_read": (
                    avg_calls_per_total_read
                ),
                "avg_calls_per_positive_read": (
                    avg_calls_per_positive_read
                ),
            }
        )

    return records


def parse_log_file(file_path):
    """
    Read and parse a human-readable BrdU read-percentage log file.

    Parameters
    ----------
    file_path
        Path to the log file.

    Returns
    -------
    list[dict]
        Parsed summary records.
    """
    with open(
        file_path,
        "r",
        encoding="utf-8",
        errors="replace",
    ) as handle:
        text = handle.read()

    return parse_log_text(
        text,
        source_name=file_path,
    )


def build_record(
    source_name,
    input_bam,
    threshold,
    total_reads,
    brdu_reads,
    total_calls,
):
    """
    Normalize raw summary fields and calculate all dashboard metrics.

    This helper creates the common internal record format used by both the log
    and TSV parsing paths.

    Parameters
    ----------
    source_name
        Source file retained for traceability.

    input_bam
        BAM filename or path from the summary.

    threshold
        BrdU modification threshold.

    total_reads
        Total mapped primary read count.

    brdu_reads
        Reads containing at least one passing BrdU call.

    total_calls
        Total number of passing BrdU calls.

    Returns
    -------
    dict
        Normalized dashboard record.
    """
    total_reads = clean_int(total_reads)
    brdu_reads = clean_int(brdu_reads)
    total_calls = clean_int(total_calls)
    threshold = float(threshold)

    percent_reads_with_brdu = (
        (brdu_reads / total_reads) * 100
        if total_reads > 0
        else 0.0
    )

    avg_calls_per_total_read = (
        total_calls / total_reads
        if total_reads > 0
        else 0.0
    )

    avg_calls_per_positive_read = (
        total_calls / brdu_reads
        if brdu_reads > 0
        else 0.0
    )

    return {
        "source_name": source_name,
        "input_bam": input_bam,
        "sample": clean_sample_name(
            input_bam
        ),
        "threshold": threshold,
        "total_reads": total_reads,
        "brdu_positive_reads": brdu_reads,
        "reads_without_brdu": (
            total_reads - brdu_reads
        ),
        "total_passing_brdu_calls": (
            total_calls
        ),
        "percent_reads_with_brdu": (
            percent_reads_with_brdu
        ),
        "avg_calls_per_total_read": (
            avg_calls_per_total_read
        ),
        "avg_calls_per_positive_read": (
            avg_calls_per_positive_read
        ),
    }


def parse_summary_tsv_file(file_path):
    """
    Parse a machine-readable read-percentage summary TSV.

    A whitespace-separated fallback is retained for older or malformed exports.
    The ``ample_bam`` fallback accommodates a historical truncated header.

    Parameters
    ----------
    file_path
        Path to the summary TSV.

    Returns
    -------
    list[dict]
        Parsed summary records.
    """
    records = []

    with open(
        file_path,
        "r",
        encoding="utf-8",
        errors="replace",
    ) as handle:
        lines = [
            line.strip()
            for line in handle
            if line.strip()
        ]

    if not lines:
        return records

    delimiter = (
        "\t"
        if "\t" in lines[0]
        else None
    )

    if delimiter == "\t":
        rows = csv.DictReader(
            lines,
            delimiter="\t",
        )
    else:
        header = re.split(
            r"\s+",
            lines[0],
        )

        rows = (
            dict(
                zip(
                    header,
                    re.split(
                        r"\s+",
                        line,
                    ),
                )
            )
            for line in lines[1:]
        )

    for row in rows:
        if not row:
            continue

        input_bam = (
            row.get("sample_bam")
            or row.get("ample_bam")
        )

        threshold = row.get(
            "mod_threshold"
        )

        total_reads = row.get(
            "total_reads"
        )

        brdu_reads = row.get(
            "reads_with_brdu_calls"
        )

        total_calls = row.get(
            "brdu_calls"
        )

        if not all(
            [
                input_bam,
                threshold,
                total_reads,
                brdu_reads,
                total_calls,
            ]
        ):
            continue

        records.append(
            build_record(
                source_name=file_path,
                input_bam=input_bam,
                threshold=threshold,
                total_reads=total_reads,
                brdu_reads=brdu_reads,
                total_calls=total_calls,
            )
        )

    return records


def parse_summary_file(file_path):
    """
    Dispatch a summary file to the TSV or human-readable log parser.
    """
    if file_path.lower().endswith(".tsv"):
        return parse_summary_tsv_file(
            file_path
        )

    return parse_log_file(
        file_path
    )


# =============================================================================
# Plotting
# =============================================================================


def plot_dashboard(
    records,
    threshold,
    output_path,
):
    """
    Create and save the four-panel dashboard for one threshold.

    Parameters
    ----------
    records
        Parsed records belonging to one modification threshold.

    threshold
        Numeric BrdU modification threshold represented by the records.

    output_path
        Destination PNG path.
    """
    require_matplotlib()

    records = sorted(
        records,
        key=lambda record: sample_sort_key(
            record["sample"]
        ),
    )

    sample_labels = [
        record["sample"]
        for record in records
    ]

    total_reads = [
        record["total_reads"]
        for record in records
    ]

    brdu_reads = [
        record["brdu_positive_reads"]
        for record in records
    ]

    non_brdu_reads = [
        record["reads_without_brdu"]
        for record in records
    ]

    brdu_percent = [
        record["percent_reads_with_brdu"]
        for record in records
    ]

    avg_calls_per_read = [
        record["avg_calls_per_total_read"]
        for record in records
    ]

    x = list(
        range(
            len(records)
        )
    )

    fig = plt.figure(
        figsize=(
            14,
            9,
        )
    )

    fig.suptitle(
        "BrdU Read Summary Dashboard "
        f"(mod threshold = {threshold:.2f})",
        fontsize=20,
        fontweight="bold",
        y=0.98,
    )

    grid_spec = fig.add_gridspec(
        2,
        2,
        height_ratios=[
            1,
            1,
        ],
    )

    ax1 = fig.add_subplot(
        grid_spec[0, 0]
    )

    ax2 = fig.add_subplot(
        grid_spec[0, 1]
    )

    ax3 = fig.add_subplot(
        grid_spec[1, 0]
    )

    ax4 = fig.add_subplot(
        grid_spec[1, 1]
    )

    # -------------------------------------------------------------------------
    # Panel 1: total reads and BrdU-positive reads
    # -------------------------------------------------------------------------

    ax1.bar(
        x,
        brdu_reads,
        label="BrdU+ Reads",
        color="tab:blue",
        width=0.6,
    )

    ax1.bar(
        x,
        non_brdu_reads,
        bottom=brdu_reads,
        label="Total Reads (minus BrdU+)",
        color="lightgray",
        edgecolor="gray",
        width=0.6,
    )

    ax1.set_title(
        "Total Reads and BrdU-Positive Reads",
        fontsize=16,
        fontweight="bold",
        pad=12,
    )

    ax1.set_ylabel(
        "Reads",
        fontsize=12,
    )

    ax1.set_xticks(x)

    ax1.set_xticklabels(
        sample_labels,
        rotation=0,
    )

    ax1.yaxis.set_major_formatter(
        FuncFormatter(
            thousands_formatter
        )
    )

    ax1.grid(
        axis="y",
        linestyle="--",
        alpha=0.5,
    )

    ax1.set_axisbelow(True)

    ax1.legend(
        loc="upper right",
        fontsize=10,
    )

    ymax_total = (
        max(total_reads)
        if total_reads
        else 1
    )

    ax1.set_ylim(
        0,
        ymax_total * 1.15,
    )

    # Add total-read labels above each stacked bar.
    for index, total in enumerate(
        total_reads
    ):
        ax1.text(
            index,
            total + ymax_total * 0.015,
            f"{total:,}",
            ha="center",
            va="bottom",
            fontsize=10,
            fontweight="bold",
        )

    # Place BrdU-positive labels inside large segments and above very small
    # segments so that the values remain readable.
    for index, value in enumerate(
        brdu_reads
    ):
        if value < ymax_total * 0.03:
            y_position = (
                value
                + ymax_total * 0.01
            )
            text_color = "black"
            vertical_alignment = "bottom"
        else:
            y_position = value / 2
            text_color = "white"
            vertical_alignment = "center"

        ax1.text(
            index,
            y_position,
            f"{value:,}",
            ha="center",
            va=vertical_alignment,
            fontsize=10,
            color=text_color,
            fontweight="bold",
        )

    # -------------------------------------------------------------------------
    # Panel 2: percentage of reads containing BrdU calls
    # -------------------------------------------------------------------------

    bars2 = ax2.bar(
        x,
        brdu_percent,
        color="tab:blue",
        width=0.5,
    )

    ax2.set_title(
        "% Reads with BrdU Calls",
        fontsize=16,
        fontweight="bold",
        pad=12,
    )

    ax2.set_ylabel(
        "Percent (%)",
        fontsize=12,
    )

    ax2.set_xticks(x)

    ax2.set_xticklabels(
        sample_labels,
        rotation=0,
    )

    ax2.grid(
        axis="y",
        linestyle="--",
        alpha=0.5,
    )

    ax2.set_axisbelow(True)

    ymax_pct = (
        max(brdu_percent)
        if brdu_percent
        else 1
    )

    ax2.set_ylim(
        0,
        max(
            25,
            ymax_pct * 1.25,
        ),
    )

    for bar, value in zip(
        bars2,
        brdu_percent,
    ):
        ax2.text(
            bar.get_x()
            + bar.get_width() / 2,
            value
            + max(
                0.2,
                ymax_pct * 0.02,
            ),
            f"{value:.2f}%",
            ha="center",
            va="bottom",
            fontsize=10,
            fontweight="bold",
        )

    # -------------------------------------------------------------------------
    # Panel 3: average passing BrdU calls per total read
    # -------------------------------------------------------------------------

    bars3 = ax3.bar(
        x,
        avg_calls_per_read,
        color="tab:green",
        width=0.55,
    )

    ax3.set_title(
        "Average BrdU Calls per Total Read",
        fontsize=16,
        fontweight="bold",
        pad=12,
    )

    ax3.set_ylabel(
        "Avg Calls/Read",
        fontsize=12,
    )

    ax3.set_xticks(x)

    ax3.set_xticklabels(
        sample_labels,
        rotation=0,
    )

    ax3.grid(
        axis="y",
        linestyle="--",
        alpha=0.5,
    )

    ax3.set_axisbelow(True)

    ymax_avg = (
        max(avg_calls_per_read)
        if avg_calls_per_read
        else 1
    )

    ax3.set_ylim(
        0,
        max(
            5,
            ymax_avg * 1.2,
        ),
    )

    for bar, value in zip(
        bars3,
        avg_calls_per_read,
    ):
        ax3.text(
            bar.get_x()
            + bar.get_width() / 2,
            value
            + max(
                0.05,
                ymax_avg * 0.03,
            ),
            f"{value:.2f}",
            ha="center",
            va="bottom",
            fontsize=10,
            fontweight="bold",
        )

    # -------------------------------------------------------------------------
    # Panel 4: summary table
    # -------------------------------------------------------------------------

    ax4.set_title(
        "Summary Table",
        fontsize=16,
        fontweight="bold",
        pad=12,
    )

    ax4.axis("off")

    column_labels = [
        "Sample",
        "Total Reads",
        "BrdU+ Reads",
        "BrdU+ %",
        "Avg Calls/Read",
    ]

    table_data = [
        [
            record["sample"],
            f"{record['total_reads']:,}",
            f"{record['brdu_positive_reads']:,}",
            (
                f"{record['percent_reads_with_brdu']:.2f}%"
            ),
            (
                f"{record['avg_calls_per_total_read']:.2f}"
            ),
        ]
        for record in records
    ]

    table = ax4.table(
        cellText=table_data,
        colLabels=column_labels,
        cellLoc="center",
        loc="center",
        bbox=[
            0.0,
            0.05,
            1.0,
            0.92,
        ],
    )

    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(
        1.05,
        1.4,
    )

    for (
        row,
        column,
    ), cell in table.get_celld().items():
        if row == 0:
            cell.set_text_props(
                weight="bold"
            )

    plt.tight_layout(
        rect=[
            0,
            0,
            1,
            0.96,
        ]
    )

    fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)


# =============================================================================
# Main
# =============================================================================


def parse_args():
    """
    Define and parse command-line arguments for launcher and plotting modes.

    Returns
    -------
    argparse.Namespace
        Parsed command-line settings.
    """
    workflow_root = workflow_root_from_script()

    default_log_dir = os.path.join(
        workflow_root,
        "logs",
        "read_pct",
    )

    parser = argparse.ArgumentParser(
        description=(
            "Submit a SLURM job to generate BrdU read-summary dashboards "
            "from read_pct summary .log/.tsv files."
        )
    )

    parser.add_argument(
        "inputs",
        nargs="*",
        help=(
            "Optional summary .log/.tsv file(s) or directories. "
            "If omitted, the script prompts for files from "
            f"{default_log_dir}."
        ),
    )

    parser.add_argument(
        "--run-job",
        action="store_true",
        help=(
            "Run plotting directly. Used internally by the "
            "submitted SLURM job."
        ),
    )

    parser.add_argument(
        "--workflow-root",
        default=workflow_root,
        help=(
            "Workflow root. Defaults to the repository root "
            "inferred from this script."
        ),
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        default=os.path.join(
            "results",
            "read_pct_dashboard",
        ),
        help=(
            "Directory to save dashboard PNG files."
        ),
    )

    parser.add_argument(
        "--log-dir",
        default=os.path.join(
            "logs",
            "read_pct",
        ),
        help=(
            "Default directory used for interactive "
            "summary file selection."
        ),
    )

    return parser.parse_args()


def run_dashboard(args):
    """
    Parse selected summary files and generate one dashboard per threshold.

    Duplicate records sharing sample, threshold, and input BAM are collapsed.
    The last parsed record for a duplicate key is retained.

    Parameters
    ----------
    args
        Parsed command-line arguments.
    """
    default_log_dir = os.path.abspath(
        args.log_dir
    )

    output_dir = os.path.abspath(
        args.output_dir
    )

    os.makedirs(
        output_dir,
        exist_ok=True,
    )

    # When no explicit input is provided in run mode, scan the default log
    # directory.
    input_paths = (
        args.inputs
        or [default_log_dir]
    )

    input_files = normalize_existing_inputs(
        input_paths,
        default_log_dir,
    )

    all_records = []

    # Parse every selected input independently so one malformed file does not
    # prevent other valid files from being processed.
    for file_path in input_files:
        try:
            records = parse_summary_file(
                file_path
            )

            if records:
                all_records.extend(records)
            else:
                print(
                    "[WARNING] No valid summary records found in: "
                    f"{file_path}"
                )
        except Exception as exc:
            print(
                f"[WARNING] Failed to parse {file_path}: "
                f"{exc}"
            )

    if not all_records:
        raise SystemExit(
            "[ERROR] No usable records were parsed "
            "from the provided files."
        )

    # Deduplicate records by biological sample, threshold, and original BAM.
    #
    # Dictionary assignment means the last parsed duplicate replaces earlier
    # duplicates with the same key.
    deduplicated_records = {}

    for record in all_records:
        key = (
            record["sample"],
            record["threshold"],
            record["input_bam"],
        )

        deduplicated_records[key] = record

    all_records = list(
        deduplicated_records.values()
    )

    # Group records so each threshold receives its own dashboard PNG.
    grouped = defaultdict(list)

    for record in all_records:
        grouped[
            record["threshold"]
        ].append(record)

    for threshold in sorted(
        grouped.keys()
    ):
        records = grouped[threshold]

        output_name = (
            "brdu_read_summary_dashboard."
            f"threshold_{threshold_label(threshold)}.png"
        )

        output_path = os.path.join(
            output_dir,
            output_name,
        )

        plot_dashboard(
            records,
            threshold,
            output_path,
        )

        print(
            f"[INFO] Wrote: {output_path}"
        )

    print("[INFO] Done.")


def submit_dashboard(args):
    """
    Collect interactive inputs when needed and submit the plotting job to SLURM.

    Parameters
    ----------
    args
        Parsed command-line arguments.
    """
    workflow_root = os.path.abspath(
        args.workflow_root
    )

    default_log_dir = os.path.abspath(
        args.log_dir
    )

    output_dir = os.path.abspath(
        args.output_dir
    )

    slurm_log_dir = os.path.join(
        workflow_root,
        "logs",
        "read_pct_dashboard",
    )

    os.makedirs(
        output_dir,
        exist_ok=True,
    )

    os.makedirs(
        slurm_log_dir,
        exist_ok=True,
    )

    # Prompt for summary files only when no positional inputs were supplied.
    input_paths = args.inputs

    if not input_paths:
        input_paths = prompt_for_summary_files(
            default_log_dir
        )

    input_files = normalize_existing_inputs(
        input_paths,
        default_log_dir,
    )

    # Locate sbatch through PATH rather than assuming a fixed executable path.
    sbatch_path = shutil.which(
        "sbatch"
    )

    if sbatch_path is None:
        raise SystemExit(
            "[ERROR] The sbatch command is unavailable. "
            "Run this script on the cluster login node, "
            "or run the plotting step directly with --run-job."
        )

    # Include a launcher timestamp in SLURM filenames so separate submissions
    # do not overwrite each other's logs.
    timestamp = dt.datetime.now().strftime(
        "%Y%m%d_%H%M%S"
    )

    script_path = absolute_existing_file(
        __file__
    )

    slurm_stdout = os.path.join(
        slurm_log_dir,
        (
            "brdu_read_pct_dashboard."
            f"{timestamp}.%j.slurm.log"
        ),
    )

    slurm_stderr = os.path.join(
        slurm_log_dir,
        (
            "brdu_read_pct_dashboard."
            f"{timestamp}.%j.slurm.err"
        ),
    )

    # Build the submitted command as a list so paths are passed safely without
    # shell quoting or string concatenation.
    command = [
        sbatch_path,
        "--parsable",
        "--job-name=brdu_read_pct_dashboard",
        f"--chdir={workflow_root}",
        f"--output={slurm_stdout}",
        f"--error={slurm_stderr}",
        script_path,
        "--run-job",
        "--workflow-root",
        workflow_root,
        "--log-dir",
        default_log_dir,
        "--output-dir",
        output_dir,
        *input_files,
    ]

    print(
        f"[INFO] Workflow root: {workflow_root}"
    )

    print(
        f"[INFO] Output directory: {output_dir}"
    )

    print("[INFO] Input log file(s):")

    for input_file in input_files:
        print(
            f"[INFO]   {input_file}"
        )

    print()

    print(
        "[INFO] Submitting BrdU read-summary dashboard "
        "job to SLURM..."
    )

    result = subprocess.run(
        command,
        check=True,
        text=True,
        capture_output=True,
    )

    job_id = result.stdout.strip()

    # Some SLURM installations append a cluster name after a semicolon.
    log_job_id = job_id.split(";")[0]

    print(
        f"[INFO] Submitted SLURM job: {job_id}"
    )

    print(
        "[INFO] Expected SLURM log: "
        f"{slurm_stdout.replace('%j', log_job_id)}"
    )

    print(
        "[INFO] Expected SLURM err: "
        f"{slurm_stderr.replace('%j', log_job_id)}"
    )

    print(
        "[INFO] Expected dashboard output pattern:"
    )

    print(
        "[INFO]   "
        f"{output_dir}/"
        "brdu_read_summary_dashboard."
        "threshold_<threshold>.png"
    )


def main():
    """
    Prepare the environment, normalize paths, and route to submit or run mode.
    """
    args = parse_args()

    # Create or reuse the local environment and restart this script with its
    # Python executable when necessary.
    ensure_running_inside_prepared_environment()

    args.workflow_root = os.path.abspath(
        args.workflow_root
    )

    # Resolve relative log and output paths against the selected workflow root.
    if not os.path.isabs(args.log_dir):
        args.log_dir = os.path.join(
            args.workflow_root,
            args.log_dir,
        )

    if not os.path.isabs(args.output_dir):
        args.output_dir = os.path.join(
            args.workflow_root,
            args.output_dir,
        )

    # Explicit --run-job and noninteractive SLURM execution both select direct
    # parsing and plotting. A normal terminal execution enters launcher mode.
    if (
        args.run_job
        or running_as_noninteractive_slurm_job()
    ):
        run_dashboard(args)
    else:
        submit_dashboard(args)


if __name__ == "__main__":
    main()