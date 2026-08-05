# McClure Lab Workflows

This repository contains workflow-based source code for McClure Lab sequencing
and BrdU data analysis. Most code under `src/` is organized as runnable
workflows that collect user inputs, validate files and parameters, then submit
the compute-heavy work as SLURM jobs. Supporting Python, R, and shell scripts
handle extraction, plotting, alignment, summary analysis, liftOver, and related
utility tasks.

Run workflow commands from the repository root so each script can find the
expected `data/`, `logs/`, `results/`, `src/`, and `docs/` directories.

## Project Directory

```text
McClure_Lab_Workflows/
├── data/
│   ├── bam/                 BAM inputs and BAM outputs
│   ├── bed/                 BED files, annotations, and generated BED outputs
│   ├── bedgraph/            genome-browser bedGraph outputs
│   ├── fastq/               FASTQ outputs and reference FASTA files
│   ├── liftover_chains/     liftOver chain files
│   ├── ncbi/                reference genomes and genome annotations
│   └── pod5/                POD5 inputs and DNAscent index files
├── docs/                    detailed workflow and utility documentation
├── logs/                    SLURM logs and workflow logs
├── results/                 generated plots, dashboards, and analysis outputs
├── src/
│   ├── genome_browser_workflow/
│   ├── nanopore_sequence_workflow/
│   ├── rainplot_workflow/
│   ├── utils/
└── tools/                   locally installed helper executables
```

## Workflows

### Nanopore Sequence Workflow

The nanopore sequence workflow starts from raw POD5 signal and produces aligned,
sorted, indexed BAM files. `nanopore_workflow.sh` runs Dorado basecalling and
alignment/QC. `dnascent_workflow.sh` is run after alignment when you have the
BAM, BAM index, reference FASTA, DNAscent index, and original POD5 signal; it
runs DNAscent BrdU detection and writes a BrdU-detected BAM.

Run Dorado basecalling and alignment:

```bash
bash src/nanopore_sequence_workflow/nanopore_workflow.sh
```

Optional example with positional inputs:

```bash
bash src/nanopore_sequence_workflow/nanopore_workflow.sh \
  sample_run.pod5 \
  W303_reference.fna \
  data/fastq \
  data/bam
```

Run DNAscent BrdU detection after alignment:

```bash
bash src/nanopore_sequence_workflow/dnascent_workflow.sh
```

### Genome Browser Workflow

The genome browser workflow starts from a BrdU-detected BAM and generates
strand-specific BrdU bedGraph files plus smoothed or unsmoothed per-chromosome
genome-browser plots. It uses the W303 reference under `data/ncbi/W303` by
default and writes outputs under `data/bedgraph/` and
`results/genome_browser_results/`.

Run:

```bash
bash src/genome_browser_workflow/genome_browser_workflow_script.sh
```

Optional example with positional inputs:

```bash
bash src/genome_browser_workflow/genome_browser_workflow_script.sh \
  sample.sorted.indexed.BrdU.detect.bam \
  sample_browser \
  0.5
```

### Rain Plot Workflow

The rain plot workflow extracts BrdU modified-base calls from a selected BAM
region, read ID, or read-ID list, then generates per-read rain plots with
genomic annotation panels. It can summarize BrdU in binary or mean mode, limit
the number of reads plotted, run S-phase RFB-specific plotting, and optionally
perform W303-to-sacCer liftOver.

Run:

```bash
bash src/rainplot_workflow/rainplot_workflow_script.sh
```

Optional example with positional inputs:

```bash
bash src/rainplot_workflow/rainplot_workflow_script.sh \
  sample.sorted.indexed.BrdU.detect.bam \
  CM007964.1 \
  0 \
  50000 \
  "" \
  sample_chr1_0_50000.bed
```

### Utility Scripts

The `src/utils/` directory contains helper scripts used before, after, or
between the main workflows. These include scripts for merging BAMs, extracting
BrdU-positive reads, calculating BrdU read percentages, creating Venn diagrams,
running liftOver helpers, converting annotation formats, and building read
summary dashboards.


Example utility run:

```bash
bash src/utils/merge_bams.sh
```

## More Documentation

More detailed information for each workflow and the utility scripts is available
under `docs/`:

- `docs/workflows.md`
- `docs/genome_browser_workflow.md`
- `docs/rain_plot_workflow.md`
- `docs/nanopore_sequence_workflow.md`
- `docs/utils_workflow_helpers.md`
