# Utility Workflow Helpers

## Overview

The scripts under `src/utils/` are optional helper tools. They are not required
for every workflow run, but they help inspect the data, combine related files,
extract smaller analysis-ready BAMs, summarize BrdU signal, and narrow results
down to points of interest.

This page intentionally skips these lower-level conversion/parser utilities:

- `src/utils/convert_genbank_bed_to_ucsc.py`
- `src/utils/gff_feature_parser.py`
- `src/utils/liftover.py`

The remaining utilities are described below in the usual analysis order.

## 1. Merge BAMs

Script:

```text
src/utils/merge_bams.sh
```

Use this when two BAM files need to be combined into one BAM. This is commonly
useful after demultiplexing or barcode workflows, where reads for one biological
sample may be split across two barcode BAMs.

What it does:

- lists BAM files under `data/bam`
- prompts for two different BAM files
- checks that the BAMs have matching reference-sequence headers
- merges the two BAMs with `samtools merge`
- coordinate-sorts the merged BAM
- creates a BAM index
- validates the final BAM
- compares the final alignment count against the sum of both inputs

How it helps:

This creates one combined BAM that can be passed into downstream BrdU counting,
read extraction, genome browser, or rain plot workflows instead of analyzing
barcode BAMs separately.

## 2. BrdU Read-Percentage Summaries

Scripts:

```text
src/utils/count_brdu_read_stats.py
src/utils/calculate_brdu_read_percentage.sh
```

Use these to quantify how much BrdU signal is present in a BAM before moving
into plotting or more targeted analysis.

`count_brdu_read_stats.py` is the Python helper. It reads a BAM with `pysam`,
looks at MM/ML modified-base tags, and counts:

- total primary mapped reads
- reads with MM tags
- reads with ML tags
- reads with parsed modified-base calls
- reads with at least one passing BrdU call
- total passing BrdU calls
- modified-base parse errors

`calculate_brdu_read_percentage.sh` is the user-facing SLURM wrapper. It runs
the Python helper and writes both human-readable and machine-readable summary
files.

What it reports:

- percent of primary mapped reads with at least one BrdU call
- BrdU calls per 100 primary mapped reads
- threshold used for BrdU calling
- supporting MM/ML tag counts

How it helps:

These summaries are the input for the BrdU read-summary dashboard scripts. They
also give a quick sample-level check of whether a BAM has the expected BrdU
signal before spending time on heavier downstream workflows.

Typical outputs are written under:

```text
logs/read_pct/
```

## 3. BrdU Read-Summary Dashboards

Scripts:

```text
src/utils/plot_brdu_read_summary_dashboard.py
src/utils/plot_brdu_read_summary_dashboard.sh
```

Use these after running `calculate_brdu_read_percentage.sh` on one or more BAMs.
The dashboard scripts turn read-percentage summary logs into comparison plots.

`plot_brdu_read_summary_dashboard.py` parses `.log` or `.tsv` summary files and
generates one dashboard PNG per BrdU threshold.

`plot_brdu_read_summary_dashboard.sh` is the SLURM wrapper. It lists available
summary logs under `logs/read_pct`, lets you select one or more logs, prepares a
small Python environment, and submits the dashboard job.

The generated dashboard includes:

- stacked bars for BrdU-positive reads and remaining reads
- percentage of reads containing BrdU
- average passing BrdU calls per total read
- a table summarizing the plotted values

How it helps:

The dashboard makes it easier to compare samples, controls, phases, and
thresholds at a glance. This helps decide which BAMs are worth deeper inspection
with genome browser plots, rain plots, or points-of-interest analysis.

Dashboard outputs are written under:

```text
results/read_pct_dashboard/
```

## 4. Extract BrdU-Positive BAMs

Scripts:

```text
src/utils/mod_calls_brdu.sh
src/utils/pysam_mod_calls_brdu.sh
```

Use these when you need a smaller BAM containing only reads with passing BrdU
calls.

### `mod_calls_brdu.sh`

This is the Modkit-based extractor. It is useful when the BAM is smaller or when
Modkit can comfortably process the file.

What it does:

- lists BAM files under `data/bam`
- prompts for a BrdU modification threshold
- runs `modkit extract calls`
- identifies read IDs with passing BrdU calls
- extracts complete primary mapped alignments for those reads
- coordinate-sorts and indexes the output BAM
- writes a read-ID list for the selected BrdU-positive reads

Example output naming:

```text
<sample>.Brdu_positive.threshold_0p5.bam
<sample>.Brdu_positive.threshold_0p5.bam.bai
<sample>.Brdu_positive.threshold_0p5.read_ids.txt
```

### `pysam_mod_calls_brdu.sh`

This is the PySAM-based extractor. Use it when the BAM is larger and you still
need a BrdU-positive BAM, or when streaming through reads with PySAM is a better
fit than the Modkit extraction route.

What it does:

- streams through the BAM one read at a time
- parses MM/ML tags through `pysam`
- selects modification code `b`
- keeps primary mapped reads with at least one passing BrdU call
- writes a sorted/indexed BrdU-positive BAM
- writes the selected read IDs
- writes a summary TSV

Example output naming:

```text
<sample>.Brdu_positive.pysam.threshold_0p5.bam
<sample>.Brdu_positive.pysam.threshold_0p5.bam.bai
<sample>.Brdu_positive.pysam.threshold_0p5.read_ids.txt
<sample>.Brdu_positive.pysam.threshold_0p5.summary.tsv
```

How these help:

Both scripts reduce a large modBAM to the reads that actually carry BrdU signal.
That makes follow-up inspection, plotting, and targeted analysis faster and more
focused.

## 5. LiftOver BrdU BEDs for Rain Plots

Script:

```text
src/utils/liftover_brdu_bed.py
```

Use this when rain plot BrdU BED coordinates need to be lifted from W303 into a
target sacCer coordinate system.

What it does:

- reads a BrdU BED-like file with at least six columns
- creates a temporary BED4 file with synthetic row IDs
- runs UCSC `liftOver` with the selected chain file
- reconstructs the mapped output with lifted coordinates
- preserves the original BrdU columns such as read ID, base, and probability
- writes unmapped intervals separately

How it helps:

UCSC `liftOver` treats BED6 columns as special BED fields. BrdU BED files use
those columns differently, so this helper prevents BrdU read/probability columns
from being corrupted during coordinate conversion.

This script is called by the rain plot workflow when liftOver is requested.

## 6. Further Points-of-Interest Analysis

Scripts:

```text
src/utils/top_20_poi.sh
src/utils/compare_top_20_poi.sh
src/utils/venn_diagram_generation.sh
```

Use these after generating strand-specific BrdU bedGraph files when you want to
narrow the data down to high-interest genomic positions and compare those
positions across samples.

### `top_20_poi.sh`

This script combines positive- and negative-strand bedGraph files and selects up
to 20 high-ranking BrdU points of interest.

What it does:

- prompts for positive- and negative-strand bedGraph files
- combines strand counts at matching genomic coordinates
- calculates combined fraction modified and BrdU percentage
- ranks candidate positions by fraction modified, coverage, and modified count
- excludes nearby points on the same chromosome so selected POIs are spread out
- writes a ranked CSV under `data/csv`

How it helps:

It converts genome-wide BrdU pileup signal into a small ranked list of candidate
locations for manual review or comparison.

### `compare_top_20_poi.sh`

This script compares two POI CSV files.

What it does:

- prompts for two top-POI CSVs
- identifies POIs on the same chromosome within the configured distance
- uses one-to-one matching so one POI is not counted multiple times
- writes overlap, first-unmatched, second-unmatched, and summary CSV files

How it helps:

It shows which high-interest BrdU sites are shared between samples and which are
unique to each sample.

### `venn_diagram_generation.sh`

This script builds a Venn-style summary from two original POI CSVs and the
overlap CSV produced by `compare_top_20_poi.sh`.

What it does:

- prompts for the first original POI CSV
- prompts for the second original POI CSV
- prompts for the overlap CSV
- validates that overlap rows map back to the selected source files
- writes membership and summary CSV files
- generates a Venn-style PNG

How it helps:

It turns POI comparison results into a quick visual summary of shared and unique
candidate regions.

Typical outputs:

```text
data/csv/*_top20_poi.csv
data/csv/*_overlaps.csv
data/csv/*_first_unmatched.csv
data/csv/*_second_unmatched.csv
data/csv/*_summary.csv
data/csv/*.venn_membership.csv
data/csv/*.venn_summary.csv
results/venn_diagram/*.venn_diagram.png
```

## Suggested Optional Analysis Flow

1. Merge barcode BAMs with `merge_bams.sh` when one sample is split across two
   BAM files.
2. Quantify BrdU read-level signal with `calculate_brdu_read_percentage.sh`.
3. Compare those summaries with `plot_brdu_read_summary_dashboard.sh`.
4. Extract a BrdU-positive BAM with `mod_calls_brdu.sh` for smaller files or
   `pysam_mod_calls_brdu.sh` for larger files.
5. Use `liftover_brdu_bed.py` indirectly through the rain plot workflow when
   coordinate liftOver is needed.
6. Use `top_20_poi.sh`, `compare_top_20_poi.sh`, and
   `venn_diagram_generation.sh` for deeper candidate-region analysis.

These utilities are exploratory and optional. They help narrow broad BrdU signal
down to the samples, reads, and genomic positions that are most useful for the
next workflow step.
