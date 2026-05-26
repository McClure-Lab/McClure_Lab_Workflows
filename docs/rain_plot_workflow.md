# Rain Plot Workflow

## Overview

The rain plot workflow extracts BrdU calls from a BAM file for a user-selected genomic region and generates per-read rain plots. It can also:

- optionally search for RFB-associated reads in the selected region
- generate genomic annotation panels for the plotted region
- optionally lift W303 coordinates into `sacCer1`, `sacCer2`, or `sacCer3`

All rain plot extraction starts from the lab's W303-based data. If liftOver is used, the workflow keeps the BrdU signal from the original data and remaps the genomic coordinates into the selected target strain.

The main entry point is:

```bash
bash src/rainplot_workflow/rainplot_workflow_script.sh BAM CHROM START END [READ_ID] [OUTPUT_BED]
```

## Requirements

Before running the workflow, make sure:

- you are running from the repository root
- the repository contains `src/`, `data/`, `results/`, and `docs/`
- the BAM file you want to analyze is present in `data/bam/`
- the system has the required modules available:
  `python/3.13.7`, `modkit`, `samtools`, and `R/4.5.2`
- the machine can create a local Python virtual environment under `src/rainplot_workflow/.rainplot_env`
- the machine can create a local R library under `src/rainplot_workflow/.r_library`

The workflow script resolves the repository root from its own location, so the top-level folder name does not need to be `McClure_Lab_Workflows`. A root named `workflow`, `McClure_Lab_Workflows`, or something else is fine as long as the internal layout is the same.

Expected layout:

```text
<repo-root>/
├── data/
│   ├── bam/
│   ├── bed/
│   ├── sorted_bam/
│   ├── index_sorted_bam_bai/
│   ├── liftover_chains/
│   └── ncbi/
│       ├── sacCer1/
│       ├── sacCer2/
│       └── sacCer3/
├── docs/
├── results/
│   └── rainplot_results/
└── src/
    ├── rainplot_workflow/
    └── utils/
```

Important: run the workflow while your current directory is the repository root, the directory that contains `data/`, `docs/`, `results/`, and `src/`.

## Required Parameters

The workflow expects:

```bash
bash src/rainplot_workflow/rainplot_workflow_script.sh \
  BAM_FILE \
  CHROM \
  START \
  END \
  READ_ID \
  OUTPUT_BED
```

### 1. `BAM_FILE`

The BAM filename in:

```text
data/bam/
```

Example:

```bash
synthetic_brdu_5reads.bam
```

Pass the filename only, not the full path.

### 2. `CHROM`

The chromosome identifier to extract from the BAM.

Use the W303 GenBank chromosome ID as input.

Example:

```bash
CM007964.1
```

That example corresponds to chromosome 1 in the W303-based input data.

Below is a list of the GenBank chromosome ID's:

```text
"CM007964.1": "1",  "CM007965.1": "2",  "CM007966.1": "3",  "CM007967.1": "4",
"CM007968.1": "5",  "CM007969.1": "6",  "CM007970.1": "7",  "CM007971.1": "8",
"CM007972.1": "9",  "CM007973.1": "10", "CM007974.1": "11", "CM007975.1": "12",
"CM007976.1": "13", "CM007977.1": "14", "CM007978.1": "15", "CM007979.1": "16",
 "CM007980.1": "p2-micron", "CM007981.1": "MT"
```

### 3. `START`

The genomic start coordinate for extraction.

Example:

```bash
0
```

### 4. `END`

The genomic end coordinate for extraction.

Example:

```bash
50000
```

### 5. `READ_ID`

An optional read filter.

If you want all reads in the selected region, pass:

```bash
""
```

If you provide a read ID, only that read will be used.

### 6. `OUTPUT_BED`

The output BED filename written to:

```text
data/bed/
```

Recommended naming convention:

```bash
chr#_start_to_endbases.bed
```

Example:

```bash
chr1_0_to_50000bases.bed
```

If you do not provide an output BED name, the workflow will automatically use the default style:

```bash
chr#_start_to_endbases.bed
```

## Interactive Prompts

After the script starts, it asks the user several prompts.

### 1. Phase selection

The workflow first asks:

```text
Which phase is this run for? [M/S]:
```

#### If the user enters `M`

- the phase label becomes `Mitosis`
- the workflow skips RFB motif extraction
- the workflow does not request RFB overlay on the rain plots

This is the correct choice for M-phase or mitotic datasets.

#### If the user enters `S`

- the phase label becomes `S Phase`
- the workflow enables RFB motif extraction
- the workflow asks one additional S-phase-specific prompt

### 2. S-phase rain plot mode

If the user selected `S`, the workflow then asks:

```text
For S Phase, choose rain plot mode: [a] RFB coords only, [b] without RFB, [c] with and without RFB coords:
```

#### If the user enters `a`

- only reads with matching RFB coordinates are plotted
- RFB overlay markers are drawn on the rain plots

#### If the user enters `b`

- rain plots are generated without RFB overlay
- the workflow still treats the run as S phase, but the plot output does not include RFB coordinates

#### If the user enters `c`

- rain plots are generated from the full BrdU dataset
- RFB overlay markers are drawn when available

### 3. LiftOver prompt

After BrdU extraction, the workflow asks:

```text
Would you like to do a UCSC liftOver on the BrdU BED? [y/n]:
```

#### If the user enters `n`

- no liftOver is performed
- the workflow stays in W303 coordinates
- the workflow continues plotting in the original BAM/W303 coordinate system

#### If the user enters `y`

The workflow then asks:

```text
Which target yeast strain do you want to lift over to? [sacCer1/sacCer2/sacCer3]:
```

The user must enter one of:

```text
sacCer1
sacCer2
sacCer3
```

The workflow then automatically:

- selects the correct `W303TosacCerX.over.chain.gz` file from `data/liftover_chains/`
- lifts the W303 BrdU coordinates into the selected target strain
- selects the matching genomic annotation file under `data/ncbi/`

Current annotation sources:

- `W303` -> `data/ncbi/sacCer3/genomic.gff`
- `sacCer1` -> `data/ncbi/sacCer1/sacCer1_features.bed`
- `sacCer2` -> `data/ncbi/sacCer2/sacCer2_features.bed`
- `sacCer3` -> `data/ncbi/sacCer3/genomic.gff`

Important:

- the BrdU signal still comes from the original W303 read data
- only the genomic coordinates are remapped during liftOver
- the RFB BED remains in the original coordinate system, so RFB overlay is skipped after liftOver
- this allows the user to compare the original BrdU incorporation pattern against `sacCer1`, `sacCer2`, or `sacCer3` genomic annotations

## Synthetic Dataset Examples

The synthetic dataset is tailored to show the expected signal patterns, so it is the best starting point for users learning the workflow.

### Example 1: M-phase run on chromosome 1

Use this example to generate a mitotic chromosome 1 run from the synthetic dataset:

```bash
bash src/rainplot_workflow/rainplot_workflow_script.sh \
  synthetic_brdu_5reads.bam \
  CM007964.1 \
  0 \
  50000 \
  "" \
  syn_chr1_0_to_50000bases.bed
```

Recommended prompt selections for this example:

- phase: `M`
- liftOver: `n`

This run generates 4 plots.

### Example 2: S-phase run on synthetic rDNA with RFB

Use this example to generate an rDNA rain plot that includes RFB information:

```bash
bash src/rainplot_workflow/rainplot_workflow_script.sh \
  synthetic_brdu_5reads.bam \
  rDNA_2_repeats \
  0 \
  50000 \
  "" \
  syn_rDNA_0_to_50000bases.bed
```

Recommended prompt selections for this example:

- phase: `S`
- rain plot mode: `a` or `c`
- liftOver: `n`

Choose `a` if you want the RFB-focused plot only. Choose `c` if you want the full BrdU plot with RFB coordinates overlaid.

## RFB Functionality

The workflow includes an RFB finder for the selected region, but it now runs only for S-phase datasets.

Important:

- RFB detection is intended for S phase runs
- M phase or mitotic runs skip `src/rainplot_workflow/rfb_seq_matcher.py`
- it is normal for S-phase runs to produce an empty RFB BED if no motif-positive reads are found
- an empty RFB result does not automatically mean the workflow failed

When enabled, the workflow writes the RFB BED to:

```text
data/bed/rfb_basesSTART_to_END.bed
```

## Expected Output

The workflow can produce the following files.

### Always produced

- extracted BrdU BED in `data/bed/`
- rain plot PNGs in `results/rainplot_results/`
- genomic annotation PNGs in `results/rainplot_results/`

### Produced for S-phase runs when RFB workflow is enabled

- an RFB BED file in `data/bed/`
- RFB overlay on the rain plots when the selected S-phase mode includes it

### Produced when liftOver is requested

- lifted BrdU BED in `data/bed/`
- unmapped liftOver intervals in `data/bed/`
- rain plot outputs named with a `W303_to_sacCerX_` prefix
