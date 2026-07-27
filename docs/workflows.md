# Workflows

## Overview

This repository provides workflow-based analysis utilities for the McClure Lab.

Current workflows:

- `Rain Plot Workflow`: extracts BrdU calls from BAM data and generates per-read rain plot outputs for a selected genomic region
- `Nanopore Sequence Workflow`: runs Dorado basecalling, minimap2/Dorado alignment QC, and DNAscent BrdU detection from nanopore POD5 data
- `Genome Browser Workflow`: extracts BrdU bedGraphs from modBAM files and generates smoothed or unsmoothed per-chromosome genome browser plots
- `Utility Workflow Helpers`: optional scripts for merging BAMs, summarizing BrdU signal, extracting BrdU-positive reads, lifting BrdU BEDs, and comparing points of interest

Each workflow has its own dedicated markdown file with workflow-specific usage details.

## Repository Setup

To use the workflows in this repository, first create a `workflows` directory. Then change into that directory and clone this repository there.

Example:

```bash
mkdir -p workflows
cd workflows
git clone <REPOSITORY_URL>
cd McClure_Lab_Workflows
```

This keeps the repository inside a `workflows` directory so the project layout stays organized and the workflow paths are easier to manage.

## Using the Repo

After cloning, run workflow commands from the repository root:

```bash
cd /path/to/workflows/McClure_Lab_Workflows
```

From there, users can follow the instructions in the workflow-specific markdown files under `docs/`.

Current workflow documentation:

- `docs/rain_plot_workflow.md`
- `docs/nanopore_sequence_workflow.md`
- `docs/genome_browser_workflow.md`
- `docs/utils_workflow_helpers.md`

Additional workflow documentation can be added here as new workflows are completed.
