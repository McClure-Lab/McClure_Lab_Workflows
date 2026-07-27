# Genome Browser Workflow

## Overview

The genome browser workflow generates per-chromosome BrdU pileup plots from a
BrdU-detected BAM. It extracts strand-specific BrdU bedGraph files with
`modkit pileup`, then creates either smoothed or unsmoothed genome-browser PNGs.

Main entry point:

```bash
bash src/genome_browser_workflow/genome_browser_workflow_script.sh
```

The script prompts for inputs, validates them, and submits a SLURM job. The
internal `--run-job` mode is used by the submitted SLURM job and should not be
called directly unless you are debugging the workflow.

Use this workflow after you already have a BAM containing BrdU modified-base
calls. A normal aligned BAM without MM/ML modified-base tags will fail preflight.
The workflow specifically looks for BrdU modification code `b`.

## Requirements

Run the workflow from the repository root, the directory containing `data/`,
`docs/`, `logs/`, `results/`, and `src/`.

Required tools:

- SLURM with `sbatch`
- `python3`
- `samtools`
- `modkit`

On systems with environment modules, the workflow tries to load:

```text
python/3.13.7
samtools
modkit
```

If modules are unavailable, the script continues and expects `python`,
`samtools`, and `modkit` to already be available on `PATH`.

The workflow creates and manages a local Python virtual environment:

```text
src/genome_browser_workflow/.genome_browser_env
```

The first run can take longer because Python requirements may need to be
installed.

## Default Paths

Input BAM files are selected from:

```text
data/bam/
```

Important: unlike the nanopore and DNAscent workflows, this wrapper only uses the
basename of the BAM you enter. That means the selected BAM must be directly under
`data/bam/`.

The default reference FASTA is:

```text
data/ncbi/W303/ncbi_dataset/GCA_002163515.1_ASM216351v1_genomic.fna
```

You can enter a different reference FASTA at the prompt. The script checks:

1. the exact path you entered
2. a path relative to the repository root
3. a path relative to `data/`
4. the first matching filename found anywhere under `data/`

Generated bedGraph and bedMethyl files are written to:

```text
data/bedgraph/
```

Prepared sorted/indexed BAMs are written to:

```text
data/sorted_bam/
```

Genome browser PNGs are written to:

```text
results/genome_browser_results/smoothed/
results/genome_browser_results/unsmoothed/
```

Logs are written to:

```text
logs/genome_browser_workflow/
```

## Optional Feature Tracks

The plotting scripts add these feature tracks when the files exist:

```text
data/bed/W303_g4_motifs.bed
data/bed/w303_te_and_ltrs.bed
data/bed/trna_coordinates.bed
```

If a feature BED is missing, the workflow still runs; that track is left empty.

Each plot includes signal tracks for coverage, BrdU count, BrdU percentage, G4
motifs, tRNAs, TEs, and subtelomeric regions.

## Basic Run

From the repository root:

```bash
bash src/genome_browser_workflow/genome_browser_workflow_script.sh
```

Example prompt answers:

```text
Enter the BAM filename from data/bam: sample.sorted.indexed.BrdU.detect.bam
Reference FASTA for modkit pileup [/path/to/repo/data/ncbi/W303/ncbi_dataset/GCA_002163515.1_ASM216351v1_genomic.fna]:
Modkit mod threshold [0.5]: 0.5

Choose cell-cycle phase for genome browser plot titles:
  M) Mitosis
  S) S Phase
Enter M or S: M

Choose genome browser output mode:
  1) smoothed
  2) unsmoothed
Enter smoothed or unsmoothed [smoothed]: smoothed
```

After submission, the script prints the SLURM job ID and expected output/log
paths.

## Optional Positional Arguments

You can pre-fill the BAM, output prefix, and Modkit threshold:

```bash
bash src/genome_browser_workflow/genome_browser_workflow_script.sh \
  BAM \
  output_prefix \
  mod_threshold
```

Argument meanings:

- `BAM`: BAM filename under `data/bam/`
- `output_prefix`: prefix used for bedGraph, bedMethyl, log, and PNG filenames
- `mod_threshold`: BrdU probability threshold for `modkit pileup`

Example:

```bash
bash src/genome_browser_workflow/genome_browser_workflow_script.sh \
  sample.sorted.indexed.BrdU.detect.bam \
  sample_browser \
  0.5
```

The script still prompts for the reference FASTA, phase, and plot mode.

## Prompt Reference

### BAM

Prompt:

```text
Enter the BAM filename from data/bam:
```

Enter the filename only. The workflow resolves it under:

```text
data/bam/
```

The BAM should be a BrdU-detected modBAM with MM/ML tags and BrdU code `b`, such
as a BAM produced by the DNAscent workflow.

### Reference FASTA

Prompt:

```text
Reference FASTA for modkit pileup [default_reference]:
```

Press Enter to use the default W303 reference. Enter another filename or path if
the BAM was aligned to a different reference.

The reference must match the BAM alignment reference closely enough for
`modkit pileup` to work correctly.

### Modkit Threshold

Prompt:

```text
Modkit mod threshold [0.5]:
```

Use a value from `0` to `1`, for example:

```text
0.5
0.6
B:0.5
```

The workflow normalizes values into Modkit syntax, such as `B:0.5`.

The default threshold `0.5` does not change output filenames. Nondefault
thresholds add a suffix to the output prefix:

```text
0.6  -> _06
0.25 -> _025
```

### Phase

Prompt:

```text
Enter M or S:
```

Choose:

- `M`: label plots as `Mitosis`
- `S`: label plots as `S Phase`

This affects plot titles only.

### Plot Mode

Prompt:

```text
Enter smoothed or unsmoothed [smoothed]:
```

Choose:

- `smoothed`: uses `genomic_browser_generation.py`; applies a centered rolling
  window to signal tracks
- `unsmoothed`: uses `genomic_browser_generation_unsmoothed.py`; plots raw
  signal without smoothing

Accepted values include `smoothed`, `smooth`, `s`, `1`, `unsmoothed`,
`unsmooth`, `raw`, `u`, and `2`.

## Outputs

For output prefix `sample_browser`, the workflow creates:

```text
data/bedgraph/sample_browser.full.bedmethyl
data/bedgraph/sample_browser.positive.bedgraph
data/bedgraph/sample_browser.negative.bedgraph
data/bedgraph/sample_browser.modkit.log
```

The prepared BAM used by Modkit is written under:

```text
data/sorted_bam/
```

The prepared BAM name is standardized as:

```text
<sample>.sorted.indexed.BrdU.detect.bam
<sample>.sorted.indexed.BrdU.detect.bam.bai
```

Smoothed genome browser images are written to:

```text
results/genome_browser_results/smoothed/
```

Unsmoothed genome browser images are written to:

```text
results/genome_browser_results/unsmoothed/
```

PNG filenames follow this pattern:

```text
<output_prefix>.chromosome_<chromosome>.<smoothed_or_unsmoothed>.png
```

Example:

```text
sample_browser.chromosome_12.smoothed.png
```

Workflow logs are written to:

```text
logs/genome_browser_workflow/
```

Common log names:

```text
sample_browser.<job_id>.slurm.log
sample_browser.<job_id>.slurm.err
sample_browser.<job_id>.workflow.log
```

The `.workflow.log` and `data/bedgraph/<output_prefix>.modkit.log` files are the
best places to check progress and errors.

## What Happens During a Run

The submitted SLURM job performs these steps:

1. Loads Python, samtools, and Modkit.
2. Creates or reuses `src/genome_browser_workflow/.genome_browser_env`.
3. Installs Python requirements.
4. Sorts the input BAM into `data/sorted_bam/`.
5. Indexes the prepared BAM.
6. Checks that the BAM contains MM/ML modified-base tags and BrdU code `b`.
7. Runs `modkit pileup` with the selected threshold.
8. Splits Modkit bedMethyl output into positive- and negative-strand bedGraphs.
9. Runs the selected plotting script.
10. Writes one genome-browser PNG per chromosome found in the bedGraphs.

## Troubleshooting

### No BAM files are listed

Place the BrdU-detected BAM directly under:

```text
data/bam/
```

Then rerun the script.

### The BAM fails modified-base preflight

The workflow requires MM/ML modified-base tags and BrdU modification code `b`.
A regular aligned BAM without modified-base calls is not enough. Use a
DNAscent/Dorado BrdU-detected modBAM or check which modification codes are
present in the BAM.

### Reference FASTA not found

Press Enter for the default W303 reference, or provide an explicit path. If you
provide only a filename, make sure that filename is unique under `data/`.

### Positive or negative bedGraph is empty

The workflow stops if either strand-specific bedGraph is missing or empty. Check:

- the BAM contains BrdU calls with modification code `b`
- the reference FASTA matches the BAM alignment
- the Modkit threshold is not too strict
- `data/bedgraph/<output_prefix>.modkit.log` for Modkit errors

### Missing feature tracks

Missing G4, TE, or tRNA BED files do not stop the workflow. The corresponding
plot tracks are empty.

### The first run is slow

This is expected if the local Python environment needs to be created or updated.
Later runs should reuse:

```text
src/genome_browser_workflow/.genome_browser_env
```

## Useful Commands

Print workflow help:

```bash
bash src/genome_browser_workflow/genome_browser_workflow_script.sh --help
```

View submitted or running SLURM jobs:

```bash
squeue -u "$USER"
```

Watch the workflow log:

```bash
tail -f logs/genome_browser_workflow/<output_prefix>.<job_id>.workflow.log
```

List generated smoothed plots:

```bash
ls results/genome_browser_results/smoothed/*.png
```

List generated unsmoothed plots:

```bash
ls results/genome_browser_results/unsmoothed/*.png
```
