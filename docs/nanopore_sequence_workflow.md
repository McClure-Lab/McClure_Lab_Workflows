# Nanopore Sequence Workflow

## Overview

The nanopore sequence workflow has two main user-facing entry points:

1. `src/nanopore_sequence_workflow/nanopore_workflow.sh`
2. `src/nanopore_sequence_workflow/dnascent_workflow.sh`

Run `nanopore_workflow.sh` first when starting from raw POD5 data. It submits a
GPU SLURM job that runs Dorado basecalling and then prepares an aligned,
coordinate-sorted, indexed BAM for downstream analysis. It can use either:

- Dorado BAM output, where Dorado basecalls and aligns directly to a reference
- Dorado FASTQ output, followed by minimap2 alignment

Run `dnascent_workflow.sh` after you have an aligned BAM, its BAM index, a
reference FASTA, a DNAscent index, and the original POD5 signal. It submits a GPU
SLURM job that runs `DNAscent detect` in a Singularity/Apptainer container and
writes a BrdU-detected BAM.

Normal users should run these scripts from the repository root:

```bash
bash src/nanopore_sequence_workflow/nanopore_workflow.sh
bash src/nanopore_sequence_workflow/dnascent_workflow.sh
```

Both scripts prompt for inputs, validate paths and settings, then submit the
compute work to SLURM. Their internal `--run-job` modes are used by the submitted
jobs and should not be called directly unless you are debugging.

## Default Paths and Explicit Paths

The workflows support default project locations, but you can also provide
explicit paths for large files.

Default project paths:

```text
data/pod5/   POD5 files, POD5 directories, and DNAscent index files
data/fastq/  FASTQ outputs and commonly used reference FASTA files
data/bam/    BAM outputs and BAM index files
logs/        SLURM and workflow logs
```

POD5 files are often large and may be impractical to copy into the repository.
For POD5 inputs, you can enter either:

- a filename or directory name under `data/pod5/`
- an explicit path to a `.pod5` file
- an explicit path to a directory containing `.pod5` files

For example, all of these are valid if the files exist:

```text
sample.pod5
sample_pod5_directory
/scratch/my_project/raw_pod5/sample.pod5
/scratch/my_project/raw_pod5/sample_run/
```

Output directories for `nanopore_workflow.sh` must resolve under the repository
`data/` directory. If you enter `fastq`, it resolves to `data/fastq`. If you
enter `data/custom_bam`, it resolves to `<repo-root>/data/custom_bam`.

## Requirements

Run from the repository root, the directory containing `data/`, `docs/`,
`logs/`, `results/`, and `src/`.

Required for `nanopore_workflow.sh`:

- SLURM with `sbatch`
- GPU access for Dorado
- `dorado`
- `samtools`
- `minimap2` when using the FASTQ plus minimap2 path

Required for `dnascent_workflow.sh`:

- SLURM with `sbatch`
- GPU access for DNAscent
- `samtools`
- `singularity` or `apptainer`
- a DNAscent container image

On systems with environment modules, the scripts try to load the needed modules.
If modules are unavailable, they continue and expect the required tools to
already be available on `PATH`.

## Repository Layout

```text
<repo-root>/
├── data/
│   ├── pod5/
│   ├── fastq/
│   └── bam/
├── logs/
│   ├── nanopore_sequence_workflow/
│   └── dnascent_workflow/
├── src/
│   └── nanopore_sequence_workflow/
│       ├── nanopore_workflow.sh
│       ├── dorado_basecall.sh
│       ├── minimap2_alignment.sh
│       └── dnascent_workflow.sh
└── docs/
```

## Step 1: Dorado and Alignment

Use `nanopore_workflow.sh` to go from POD5 signal to aligned BAM.

```bash
bash src/nanopore_sequence_workflow/nanopore_workflow.sh
```

Optional positional arguments:

```bash
bash src/nanopore_sequence_workflow/nanopore_workflow.sh \
  POD5 \
  REFERENCE_FASTA \
  FASTQ_OUTPUT_DIR \
  BAM_OUTPUT_DIR
```

Argument meanings:

- `POD5`: a `.pod5` file or directory of `.pod5` files; resolved directly first,
  then under `data/pod5/`
- `REFERENCE_FASTA`: resolved directly first, then under `data/fastq/`, then
  under `data/`, then by filename anywhere under `data/`
- `FASTQ_OUTPUT_DIR`: optional output directory for FASTQ mode; defaults to
  `data/fastq`
- `BAM_OUTPUT_DIR`: optional BAM output directory; defaults to `data/bam`

### Common Dorado BAM Run

This is the default and recommended path when you want an aligned BAM directly
from Dorado:

```bash
bash src/nanopore_sequence_workflow/nanopore_workflow.sh
```

Example prompt answers:

```text
Enter POD5 file/directory from data/pod5 or an explicit path: /scratch/project/run42_pod5/
Enter reference FASTA filename from data/fastq, data/, or an explicit path: W303.fasta
Should Dorado produce fastq or bam? [bam]: bam
Enter BAM output directory under data/ [/path/to/repo/data/bam]:
Run dorado demux after basecalling? (yes or no) [no]: no
Enter --dorado-model (sup, hac, fast, or explicit model path) [sup]: sup
Enter --device [cuda:0]: cuda:0
Enter --min-qscore [6]: 6
Enter --max-reads (optional; press Enter to basecall all reads):
Enter --emit-moves (yes or no; BAM output only) [yes]: yes
Enter --modified-bases (optional; press Enter to skip):
Enter --preset [map-ont]: map-ont
Enter --threads [8]: 8
Enter --secondary (yes or no) [no]: no
Enter --sort (yes required for QC) [yes]: yes
Enter --index (yes required for QC) [yes]: yes
Enter --min-mapq [20]: 20
Enter --min-read-length [1000]: 1000
Enter --primary-only (yes or no) [yes]: yes
```

This produces a Dorado-aligned BAM, then creates a filtered raw BAM, a final
sorted/indexed BAM, a BAM index, and alignment QC logs.

### FASTQ Plus Minimap2 Run

Choose FASTQ output if you want Dorado to write FASTQ first and then align with
minimap2:

```text
Should Dorado produce fastq or bam? [bam]: fastq
```

In this mode:

- Dorado writes FASTQ to `data/fastq` by default
- `--emit-fastq` is enabled
- `--emit-moves` is ignored because FASTQ cannot store move-table tags
- minimap2 aligns the FASTQ to the reference
- samtools sorts and indexes the final BAM

### Dorado Demultiplexing

Prompt:

```text
Run dorado demux after basecalling? (yes or no) [no]:
```

Choose `yes` only when the run has barcodes and you know the Dorado barcode kit
name. The script then asks:

```text
Enter barcode-kit / --kit-name:
```

Demultiplexed outputs are written under:

```text
data/bam/<pod5_prefix>_<job_id>/
```

Each barcode FASTQ or BAM is aligned/QC'd separately.

### Nanopore Output Files

For a POD5 file or directory with prefix `sample` and SLURM job ID `12345`, the
main outputs are:

```text
data/bam/sample_12345.bam
data/bam/sample_12345.dorado.raw.bam
data/bam/sample.sorted.indexed_12345.bam
data/bam/sample.sorted.indexed_12345.bam.bai
data/bam/sample_12345.alignment_qc.txt
```

For FASTQ mode, the FASTQ output is:

```text
data/fastq/sample_12345.fastq
```

Log files are written to:

```text
logs/nanopore_sequence_workflow/
```

Common log names:

```text
sample.12345.slurm.log
sample.12345.slurm.err
sample.12345.dorado.log
sample.12345.dorado_alignment.log
sample.12345.minimap2.log
sample.12345.barcode_alignment.log
```

## Step 2: DNAscent BrdU Detection

Use `dnascent_workflow.sh` after the first workflow has produced a sorted,
indexed BAM.

```bash
bash src/nanopore_sequence_workflow/dnascent_workflow.sh
```

Optional positional arguments:

```bash
bash src/nanopore_sequence_workflow/dnascent_workflow.sh \
  BAM \
  BAM_INDEX \
  REFERENCE_FASTA \
  DNASCENT_INDEX \
  POD5
```

Argument meanings:

- `BAM`: filename under `data/bam/` or an explicit path
- `BAM_INDEX`: filename under `data/bam/` or an explicit path; default prompt
  value is `<BAM>.bai`
- `REFERENCE_FASTA`: filename under `data/fastq/`, filename elsewhere under
  `data/`, or an explicit path
- `DNASCENT_INDEX`: filename under `data/pod5/` or an explicit path
- `POD5`: `.pod5` file or directory; resolved directly first, then under
  `data/pod5/`

Example run using the BAM from `nanopore_workflow.sh` and a POD5 directory that
stays on scratch storage:

```bash
bash src/nanopore_sequence_workflow/dnascent_workflow.sh \
  data/bam/sample.sorted.indexed_12345.bam \
  data/bam/sample.sorted.indexed_12345.bam.bai \
  data/fastq/W303.fasta \
  data/pod5/sample.dnascent.index \
  /scratch/project/run42_pod5/
```

The script will still prompt for DNAscent runtime parameters.

### DNAscent Prompts

If you do not provide positional arguments, the script lists available files and
asks for:

```text
Enter BAM filename from data/bam or an explicit path:
Enter BAM index filename from data/bam or an explicit path [<BAM>.bai]:
Enter reference FASTA filename from data/fastq, data/, or an explicit path:
Enter DNAscent index filename from data/pod5 or an explicit path:
Enter POD5 file/directory from data/pod5 or an explicit path:
```

Then it asks for DNAscent and SLURM settings:

```text
Enter DNAscent container image path [/cluster/singularity_images/DNAscent.sif]:
Enter DNAscent threads / SLURM cpus-per-task [12]:
Enter DNAscent --GPU value [0]:
Enter DNAscent -q minimum alignment/read quality [20]:
Enter DNAscent -l minimum read length [1000]:
Enter SLURM memory [36G]:
Enter SLURM time [12:00:00]:
```

The default container path can also be changed with:

```bash
export DNASCENT_IMAGE=/path/to/DNAscent.sif
```

### DNAscent Output Files

DNAscent writes output BAMs to:

```text
data/bam/
```

For input prefix `sample` and SLURM job ID `67890`, expected output is:

```text
data/bam/sample.sorted.indexed.BrdU.detect_67890.bam
```

DNAscent logs are written to:

```text
logs/dnascent_workflow/
```

Common log names:

```text
sample.67890.slurm.log
sample.67890.slurm.err
sample.67890.dnascent.log
```

The DNAscent job binds the workflow root plus the directories containing the BAM,
BAM index, reference FASTA, DNAscent index, POD5 input, output directory, and log
directory into the container. This is why explicit external POD5 paths can work
as long as the compute node and container runtime can access them.

## Recommended End-to-End Flow

1. Put the reference FASTA somewhere under `data/`, commonly `data/fastq/`, or
   keep it elsewhere and provide an explicit path.
2. Keep POD5 files in `data/pod5/` if convenient, or leave large POD5 files on a
   mounted scratch/project path and enter that explicit path.
3. Run `nanopore_workflow.sh`.
4. Wait for the final sorted/indexed BAM in `data/bam/`.
5. Confirm the matching `.bam.bai` exists.
6. Run `dnascent_workflow.sh` with the BAM, BAM index, reference FASTA, DNAscent
   index, and original POD5 input.
7. Use the DNAscent BrdU-detected BAM for downstream BrdU workflows.

## Troubleshooting

### POD5 input not found

The scripts check the exact path first, then `data/pod5/`. If your POD5 data is
on scratch storage, enter the full path to the `.pod5` file or the directory
containing `.pod5` files.

### POD5 directory contains no reads

The directory must contain `.pod5` files directly inside it. The current checks
look one directory level deep, not recursively through nested folders.

### Reference FASTA not found

Use an explicit FASTA path, or place the reference somewhere under `data/`.
The scripts search `data/fastq/` first and then search by exact filename under
`data/`.

### Sorting or indexing errors

The alignment QC requires sorted and indexed BAMs. Leave these prompts as `yes`:

```text
Enter --sort (yes required for QC) [yes]: yes
Enter --index (yes required for QC) [yes]: yes
```

### DNAscent container not found

Check the container path shown in the prompt. The default is:

```text
/cluster/singularity_images/DNAscent.sif
```

Use `DNASCENT_IMAGE` or enter an explicit path if your container is elsewhere.

### DNAscent cannot see external POD5 data

The workflow binds the POD5 directory into the container, but the path still must
exist on the compute node. Use storage that is mounted on the nodes where the GPU
job runs.

## Useful Commands

Print the nanopore workflow help:

```bash
bash src/nanopore_sequence_workflow/nanopore_workflow.sh --help
```

Print the DNAscent workflow help:

```bash
bash src/nanopore_sequence_workflow/dnascent_workflow.sh --help
```

View submitted or running SLURM jobs:

```bash
squeue -u "$USER"
```

Watch a Dorado log:

```bash
tail -f logs/nanopore_sequence_workflow/sample.12345.dorado.log
```

Watch a DNAscent log:

```bash
tail -f logs/dnascent_workflow/sample.67890.dnascent.log
```
