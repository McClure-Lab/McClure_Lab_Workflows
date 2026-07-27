# Rain Plot Workflow

## Overview

The rain plot workflow extracts BrdU modified-base calls from a BAM file and
generates per-read rain plots with genomic annotation panels. It is designed to
run from the repository root and submit the compute-heavy work to SLURM.

The workflow can:

- extract all reads from a genomic interval
- extract one read from a genomic interval
- extract a ranked or predefined set of reads from a TSV/list file
- summarize BrdU signal in 100-thymidine windows using binary or mean mode
- limit how many reads are plotted
- generate S-phase RFB-focused plots and optional RFB overlays
- optionally lift W303 BrdU BED coordinates to `sacCer1`, `sacCer2`, or
  `sacCer3`
- generate genomic feature panels and stack them below each rain plot

Normal users should run:

```bash
bash src/rainplot_workflow/rainplot_workflow_script.sh
```

The script prompts for inputs, validates them, and submits a SLURM job. The
internal `--run-job` mode is used by the submitted SLURM job and should not be
called directly unless you are debugging the workflow.

## Requirements

Run the workflow from the repository root, which is the directory containing
`data/`, `docs/`, `results/`, and `src/`.

Required command-line tools:

- `sbatch` from SLURM
- `python3`
- `samtools`
- `Rscript`
- `modkit`
- `wget` only if you request liftOver and `tools/liftOver` is missing

On systems with environment modules, the script tries to load:

```text
python/3.13.7
modkit
samtools
R/4.5.2
```

If modules are not available, the script continues and expects `python3`,
`samtools`, `modkit`, and `Rscript` to already be available on `PATH`.

The workflow creates and manages these local environments:

- Python virtual environment: `src/rainplot_workflow/.rainplot_env`
- R package library: `src/rainplot_workflow/.r_library`

The first run can take longer because Python packages and R/Bioconductor
packages may need to be installed.

## Expected Repository Layout

```text
<repo-root>/
├── data/
│   ├── bam/
│   ├── bed/
│   ├── liftover_chains/
│   └── ncbi/
│       ├── sacCer1/
│       ├── sacCer2/
│       └── sacCer3/
├── docs/
├── logs/
│   └── rainplot_workflow/
├── results/
│   └── rainplot_results/
├── src/
│   ├── rainplot_workflow/
│   └── utils/
└── tools/
```

Input BAM files must be placed in:

```text
data/bam/
```

Pass only the BAM filename when prompted, for example:

```text
synthetic_brdu_5reads.bam
```

## Coordinate System

Use 0-based, half-open coordinates for `start` and `end`, the same style used by
BED files. For example, `0` to `50000` means positions `0` through `49999`.

For W303 BAM files, enter W303 GenBank chromosome IDs. Common chromosome IDs are:

```text
CM007964.1 = chr1    CM007965.1 = chr2    CM007966.1 = chr3
CM007967.1 = chr4    CM007968.1 = chr5    CM007969.1 = chr6
CM007970.1 = chr7    CM007971.1 = chr8    CM007972.1 = chr9
CM007973.1 = chr10   CM007974.1 = chr11   CM007975.1 = chr12
CM007976.1 = chr13   CM007977.1 = chr14   CM007978.1 = chr15
CM007979.1 = chr16   CM007980.1 = p2-micron
CM007981.1 = MT
```

The annotation plotting code also accepts common numeric, `chr1`, Roman-numeral,
and sacCer RefSeq-style chromosome names, but W303 GenBank IDs are the safest
choice for extracting from W303 BAMs.

## Quick Start

From the repository root:

```bash
bash src/rainplot_workflow/rainplot_workflow_script.sh
```

For a basic M-phase interval run, answer the prompts like this:

```text
Enter the BAM filename from data/bam: synthetic_brdu_5reads.bam
Enter read IDs TSV/list file [blank for none]:
Enter binary or mean [binary]: binary
Enter binary BrdU probability threshold [0.5]: 0.5
Enter chromosome GenBank ID: CM007964.1
Enter start coordinate: 0
Enter end coordinate: 50000
Enter read ID filter [blank for all reads in region]:
How many reads should be plotted? [blank for all reads in region]:
Output BED filename [synthetic_brdu_5reads.chr1.0.50000based.binary.bed]:
Which phase is this run for? [M/S]: M
Would you like to do a UCSC liftOver on the BrdU BED? [y/n]: n
```

After submission, the script prints the SLURM job ID and the expected log files.
Use those logs to follow the run.

## Optional Positional Arguments

You can pre-fill the first inputs on the command line:

```bash
bash src/rainplot_workflow/rainplot_workflow_script.sh \
  BAM \
  CHROM \
  START \
  END \
  READ_ID \
  OUTPUT_BED \
  READ_IDS_FILE
```

The script still asks the remaining workflow prompts before submitting the SLURM
job.

Argument meanings:

- `BAM`: filename under `data/bam/`
- `CHROM`: chromosome or contig for interval extraction
- `START`: 0-based start coordinate
- `END`: half-open end coordinate
- `READ_ID`: optional single-read substring filter; use `""` for all reads
- `OUTPUT_BED`: optional output BED filename; written under `data/bed/`
- `READ_IDS_FILE`: optional TSV/list file containing read IDs

Example:

```bash
bash src/rainplot_workflow/rainplot_workflow_script.sh \
  synthetic_brdu_5reads.bam \
  CM007964.1 \
  0 \
  50000 \
  "" \
  synthetic_chr1.binary.bed
```

## Prompt Reference

### Read IDs File

Prompt:

```text
Enter read IDs TSV/list file [blank for none]:
```

Leave this blank for a normal interval run.

Provide a TSV, CSV, or plain-text file to plot a specific list of reads. If the
file has a header with a `read_id` column, that column is used. Otherwise, the
first field on each non-empty line is treated as the read ID.

When a read IDs file is provided:

- the script clears `CHROM`, `START`, and `END`
- BrdU extraction scans all coordinates for the listed reads
- plot order follows the read order in the file
- `How many reads should be plotted?` limits the number of plottable reads
- RFB extraction is skipped because there is no single interval to search

Do not provide both a single `READ_ID` and a `READ_IDS_FILE`.

### BrdU T-Window Mode

Prompt:

```text
Enter binary or mean [binary]:
```

Choose one:

- `binary`: each BrdU probability is converted to `0` or `1` using the selected
  threshold, then calls are averaged in each 100-thymidine window
- `mean`: raw BrdU probabilities are averaged directly in each 100-thymidine
  window

Binary mode also asks:

```text
Enter binary BrdU probability threshold [0.5]:
```

Use a value from `0` to `1`. The default is `0.5`.

Output names include the selected mode:

```text
*.binary.bed
*.mean.bed
```

### Region

For interval runs, the workflow asks:

```text
Enter chromosome GenBank ID:
Enter start coordinate:
Enter end coordinate:
```

The start and end must be integers, and start must be less than end.

### Single Read Filter

Prompt:

```text
Enter read ID filter [blank for all reads in region]:
```

Leave blank to plot all reads in the interval. Enter a read ID, or a unique
substring of a read ID, to restrict extraction to matching reads.

### Read Count

Prompt:

```text
How many reads should be plotted? [blank for all reads in region]:
```

or, for read-list runs:

```text
How many reads should be plotted? [blank for all reads in TSV/list]:
```

Leave blank to plot every selected read. Enter a positive integer to cap the
number of rain plots.

For read-list runs, the first plottable reads from the file are used in order.
For interval runs, the plotting script samples up to the requested number of
reads from the extracted region.

### Output BED Filename

Prompt:

```text
Output BED filename [default_name.bed]:
```

Press Enter to use the default. The workflow writes the BED file to:

```text
data/bed/
```

If you provide a name without `.binary` or `.mean`, the workflow appends the
selected BrdU mode automatically.

### Phase

Prompt:

```text
Which phase is this run for? [M/S]:
```

Choose:

- `M`: Mitosis. RFB motif extraction is skipped.
- `S`: S phase. RFB options are enabled.

M-phase plotting skips reads aligned outside chromosomes 1-16, including
`LYZE01000019.1`, `LYZE01000020.1`, `LYZE01000021.1`, `CM007980.1`, and
`CM007981.1`.

### S-Phase RFB Mode

If you choose `S`, the workflow asks:

```text
For S Phase, choose rain plot mode: [a] RFB coords only, [b] without RFB, [c] with and without RFB coords:
```

Choose:

- `a`: plot only reads with detected RFB motif coordinates and draw RFB overlays
- `b`: generate rain plots without RFB filtering or overlays
- `c`: plot the full selected BrdU set and draw RFB overlays when available

RFB motif extraction searches for approximate matches to:

```text
TTTACCAAGAAAGATGTAAG
```

The RFB BED is written to:

```text
data/bed/rfb_basesSTART_to_END.bed
```

It is normal for an S-phase run to produce no RFB BED rows if no motif-positive
reads are found.

### LiftOver

Prompt:

```text
Would you like to do a UCSC liftOver on the BrdU BED? [y/n]:
```

Choose `n` for a standard W303-coordinate run.

Choose `y` to lift W303 BrdU BED coordinates to one of:

```text
sacCer1
sacCer2
sacCer3
```

Required chain files must exist in `data/liftover_chains/`:

```text
W303TosacCer1.over.chain.gz
W303TosacCer2.over.chain.gz
W303TosacCer3.over.chain.gz
```

If `tools/liftOver` does not exist, the script downloads the UCSC liftOver
binary into `tools/`.

Important liftOver behavior:

- BrdU calls are extracted from the original BAM first
- only BED coordinates are lifted
- lifted outputs are used for downstream rainplot generation
- RFB overlays are disabled after liftOver because the RFB BED remains in the
  original coordinate system
- W303 G4 motifs are not overlaid after liftOver

Annotation sources:

```text
W303    -> data/ncbi/sacCer3/genomic.gff
sacCer1 -> data/ncbi/sacCer1/sacCer1_features.bed
sacCer2 -> data/ncbi/sacCer2/sacCer2_features.bed
sacCer3 -> data/ncbi/sacCer3/genomic.gff
```

## Example Runs

### M-Phase Interval Run

```bash
bash src/rainplot_workflow/rainplot_workflow_script.sh \
  synthetic_brdu_5reads.bam \
  CM007964.1 \
  0 \
  50000 \
  "" \
  synthetic_chr1.binary.bed
```

Recommended prompt answers:

```text
Enter read IDs TSV/list file [blank for none]:
Enter binary or mean [binary]: binary
Enter binary BrdU probability threshold [0.5]: 0.5
How many reads should be plotted? [blank for all reads in region]:
Which phase is this run for? [M/S]: M
Would you like to do a UCSC liftOver on the BrdU BED? [y/n]: n
```

### S-Phase rDNA Run With RFB Overlay

```bash
bash src/rainplot_workflow/rainplot_workflow_script.sh \
  synthetic_brdu_5reads.bam \
  rDNA_2_repeats \
  0 \
  50000 \
  "" \
  synthetic_rDNA.binary.bed
```

Recommended prompt answers:

```text
Enter read IDs TSV/list file [blank for none]:
Enter binary or mean [binary]: binary
Enter binary BrdU probability threshold [0.5]: 0.5
How many reads should be plotted? [blank for all reads in region]:
Which phase is this run for? [M/S]: S
For S Phase, choose rain plot mode: [a] RFB coords only, [b] without RFB, [c] with and without RFB coords: c
Would you like to do a UCSC liftOver on the BrdU BED? [y/n]: n
```

Use `a` instead of `c` if you only want reads with detected RFB coordinates.

### Ranked Read-List Run

Use this when you already have a ranked TSV of reads and want the top plottable
reads in that order.

```bash
bash src/rainplot_workflow/rainplot_workflow_script.sh \
  sample.bam \
  "" \
  "" \
  "" \
  "" \
  sample_top100.binary.bed \
  results/rainplot_results/top100_ranked_reads.tsv
```

Recommended prompt answers:

```text
Enter binary or mean [binary]: binary
Enter binary BrdU probability threshold [0.5]: 0.5
How many reads should be plotted? [blank for all reads in TSV/list]: 100
Which phase is this run for? [M/S]: M
Would you like to do a UCSC liftOver on the BrdU BED? [y/n]: n
```

Example TSV format:

```text
rank	read_id	read_length	evaluated_brdu_positions	brdu_probability_sum	average_brdu_probability
1	29e07de0-de05-4017-85aa-2646eee6e040	1787	675	452.791017	0.670802
```

## Outputs

The workflow writes BED outputs to:

```text
data/bed/
```

Common BED outputs:

- extracted BrdU BED: `data/bed/<output_basename>.bed`
- S-phase RFB BED: `data/bed/rfb_basesSTART_to_END.bed`
- lifted BrdU BED: `data/bed/W303_to_sacCerX_<output_basename>.bed`
- unmapped liftOver intervals:
  `data/bed/W303_to_sacCerX_<output_basename>_unmapped.bed`

Final image outputs are written to:

```text
results/rainplot_results/
```

The final user-facing images start with:

```text
combined_
```

The workflow creates separate rainplot and annotation PNGs first, combines them,
and then deletes those intermediate PNGs. Temporary manifest files are also
removed after combination.

Log files are written to:

```text
logs/rainplot_workflow/
```

For each submitted job, the script prints paths like:

```text
<output_basename>.<job_id>.slurm.log
<output_basename>.<job_id>.slurm.err
<output_basename>.<job_id>.workflow.log
```

The `.workflow.log` file is usually the best place to see workflow progress and
errors.

## What Happens During a Run

The submitted SLURM job performs these steps:

1. Loads required modules or uses tools already on `PATH`.
2. Creates or reuses the Python virtual environment.
3. Installs Python requirements from `src/rainplot_workflow/requirements.txt`.
4. Creates or reuses the local R library.
5. Installs missing R/Bioconductor packages.
6. Sorts and indexes the BAM if needed.
7. Extracts BrdU calls into a BED-like file.
8. Optionally runs liftOver.
9. Optionally extracts S-phase RFB motif coordinates.
10. Generates per-read rain plots.
11. Generates genomic annotation panels.
12. Combines each rain plot with its matching annotation panel.

## Troubleshooting

### No BAM files are listed

Put the input BAM under:

```text
data/bam/
```

Then rerun the script and enter the filename only.

### The job submits but fails immediately

Check:

```text
logs/rainplot_workflow/<output_basename>.<job_id>.slurm.err
logs/rainplot_workflow/<output_basename>.<job_id>.workflow.log
```

Common causes are missing `samtools`, missing `Rscript`, missing annotation
files, or a BAM name that does not exist under `data/bam/`.

### The first run is slow

This is expected if the workflow needs to install Python or R packages. Later
runs should reuse:

```text
src/rainplot_workflow/.rainplot_env
src/rainplot_workflow/.r_library
```

### The BrdU BED is empty

The workflow stops if the extracted BrdU BED is empty. Check that:

- the chromosome name exists in the BAM
- `start` and `end` cover reads with BrdU calls
- the optional read ID filter matches a read in the selected region
- the BAM contains BrdU modified-base calls with modification code `b`

### No final combined images are produced

Check the workflow log for messages about skipped reads. A read must have enough
BrdU/T-window data to plot. The rainplot code uses 100 thymidines per window and
a 45-thymidine step, so very short or sparse reads can be skipped.

### RFB output is empty

This can be normal. RFB extraction only runs for S-phase runs and only reports
reads with an approximate match to the RFB motif. If no motif-positive reads are
found, RFB-specific plots may contain no filtered reads.

### liftOver fails

Check that the target chain and annotation files exist:

```text
data/liftover_chains/W303TosacCer1.over.chain.gz
data/liftover_chains/W303TosacCer2.over.chain.gz
data/liftover_chains/W303TosacCer3.over.chain.gz
data/ncbi/sacCer1/sacCer1_features.bed
data/ncbi/sacCer2/sacCer2_features.bed
data/ncbi/sacCer3/genomic.gff
```

Also check whether the compute node can run or download `tools/liftOver`.

## Useful Commands

View submitted or running jobs:

```bash
squeue -u "$USER"
```

Watch the workflow log after submission:

```bash
tail -f logs/rainplot_workflow/<output_basename>.<job_id>.workflow.log
```

List final images:

```bash
ls results/rainplot_results/combined_*.png
```

Print script usage:

```bash
bash src/rainplot_workflow/rainplot_workflow_script.sh --help
```
