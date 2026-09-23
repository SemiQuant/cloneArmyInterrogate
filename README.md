# CloneArmy

CloneArmy is a modern Python package for analyzing haplotypes from Illumina paired-end and Oxford Nanopore (dorado SUP basecalled) amplicon sequencing data. It provides a streamlined workflow for processing FASTQ files, aligning reads, identifying sequence variants, and performing comparative analyses between samples.

## Features

- Fast paired-end read processing using BWA-MEM
- Oxford Nanopore support: dorado SUP basecalling from POD5, barcode demultiplexing, minimap2 alignment, and homopolymer-aware error handling
- Quality-based filtering of bases and alignments
- Haplotype identification and frequency analysis
- Statistical comparison between samples with FDR correction
- Interactive visualization of mutation frequencies
- Rich command-line interface with progress tracking and tabular output
- Comprehensive HTML reports
- Multi-threading support
- Support for full-length sequence analysis
- Real-time progress monitoring with progress bars
- Automatic downsampling of large FASTQ files
- Exportable results in multiple formats (CSV, JSON, Excel)

## Installation

### micromamba (recommended)

The included [`environment.yml`](environment.yml) installs CloneArmy together with all command-line tools (bwa, samtools, seqtk, minimap2):

```bash
git clone https://github.com/SemiQuant/cloneArmyInterrogate.git
cd cloneArmyInterrogate
micromamba create -f environment.yml
micromamba activate clonearmy
clonearmy --version
```

Or without cloning:

```bash
curl -LO https://raw.githubusercontent.com/SemiQuant/cloneArmyInterrogate/main/environment.yml
micromamba create -f environment.yml
micromamba activate clonearmy
```

To update CloneArmy in an existing environment:

```bash
micromamba activate clonearmy
pip install --upgrade --force-reinstall --no-deps git+https://github.com/SemiQuant/cloneArmyInterrogate.git
```

For development, install from your clone in editable mode instead:

```bash
micromamba activate clonearmy
pip install -e .
```

#### Dorado (Nanopore only)

dorado is not available from conda/bioconda. Download the build for your platform from the [dorado releases](https://github.com/nanoporetech/dorado#installation), unpack it, and put its `bin/` on `PATH` (or pass `--dorado-bin` to `clonearmy nanopore`). For example, on Apple Silicon:

```bash
DORADO_VERSION=x.y.z   # latest version from the dorado releases page
curl -LO https://cdn.oxfordnanoportal.com/software/analysis/dorado-${DORADO_VERSION}-osx-arm64.zip
unzip dorado-${DORADO_VERSION}-osx-arm64.zip -d ~/software
export PATH="$HOME/software/dorado-${DORADO_VERSION}-osx-arm64/bin:$PATH"
dorado --version
```

Use `linux-x64` for Linux with an NVIDIA GPU. SUP basecalling needs a GPU (Apple Silicon Metal or NVIDIA CUDA).

### pip

```bash
pip install git+https://github.com/SemiQuant/cloneArmyInterrogate.git
```

### Requirements

- Python ≥ 3.8
- BWA (must be installed and available in PATH)
- Samtools (must be installed and available in PATH)
- Seqtk (must be installed and available in PATH)
- minimap2 ≥ 2.27 and dorado (Nanopore only)

If you are not using `environment.yml`, install the tools with micromamba:
```bash
micromamba install -c conda-forge -c bioconda bwa samtools seqtk minimap2
```

## Usage

### Command Line Interface

#### Basic Analysis

```bash
# Basic usage with progress tracking
clonearmy run /path/to/fastq/directory reference.fasta

# With all options
clonearmy run /path/to/fastq/directory reference.fasta \
    --threads 8 \
    --output results \
    --min-base-quality 20 \
    --min-mapping-quality 30 \
    --min-read-count 10 \
    --max-file-size 100000000 \  # Target size for downsampling (100MB)
    --report  # Generate HTML report (default: true)
```

The `--max-file-size` option specifies the target size for downsampling large FASTQ files. If your input files are larger than this size, they will be automatically downsampled while maintaining paired-end relationships. This is useful for quick testing or when working with very large datasets. The size is specified in bytes (e.g., 100000000 for 100MB).

#### Comparative Analysis

```bash
# Compare two samples
clonearmy compare \
    /path/to/sample1/fastq \
    /path/to/sample2/fastq \
    reference.fasta \
    --threads 8 \
    --output comparison_results \
    --min-base-quality 20 \
    --min-mapping-quality 30 \
    --min-read-count 10 \
    --max-file-size 100000000 \  # Target size for downsampling (100MB)
    --full-length-only  # Only consider full-length sequences
```

#### Oxford Nanopore (dorado SUP)

The `nanopore` command takes raw POD5 data and basecalls it with dorado using the super-accurate (`sup`) model. dorado picks the newest SUP model that matches your flowcell and chemistry. The command can demultiplex by barcode kit, aligns the reads with minimap2, and then runs the same haplotype analysis and HTML report as the Illumina workflow.

Additional requirements: minimap2 (included in `environment.yml`) and [dorado](https://github.com/nanoporetech/dorado) on `PATH` or passed with `--dorado-bin` (see [Dorado installation](#dorado-nanopore-only)). SUP basecalling needs a GPU: Apple Silicon (Metal) or NVIDIA (CUDA). SUP on CPU is impractically slow.

```bash
# Barcoded run: demultiplex with the kit and rename barcodes to sample names
clonearmy nanopore /path/to/pod5 reference.fasta \
    --kit-name SQK-NBD114-24 \
    --barcode-map barcodes.csv \
    --threads 16 \
    --output results_ont

# Non-barcoded run: all reads form one sample named after the POD5 directory
clonearmy nanopore /path/to/pod5 reference.fasta -o results_ont
```

`barcodes.csv` has two columns, `barcode,sample`, for example `barcode01,T1`. Barcodes can also be given as numbers (`1,T1`).

Key options:

| Option | Default | Description |
|---|---|---|
| `--model` | `sup` | dorado model or model complex |
| `--kit-name` | none | Barcoding kit; omit for a single sample |
| `--device` | auto | dorado device (`metal`, `cuda:all`, `cpu`) |
| `--min-qscore` | 10 | Minimum mean read Q-score kept by dorado |
| `--min-base-quality` | 20 | Bases below this quality are masked to the reference |
| `--min-mapping-quality` | 20 | Minimum MAPQ |
| `--minimap2-preset` | `lr:hq` | minimap2 preset (`map-ont` for older, lower-accuracy data) |
| `--ignore-homopolymer-indels` / `--keep-homopolymer-indels` | ignore | Handling of 1-bp indels in reference homopolymers |
| `--homopolymer-min-length` | 3 | Minimum reference run length treated as a homopolymer |
| `--force-basecall` | off | Redo basecalling and demultiplexing even if outputs exist |
| `--full-length-tolerance` | 5 | Bases a read may miss at either reference end and still be full length |
| `--depth-thresholds` | `10,100,1000,10000` | Minimum read depths for the single-mutant coverage table |
| `--qc` / `--no-qc` | on | Read filtering, mutation load and depth-threshold tables |

Nanopore error handling: even SUP reads carry a few errors per read, mostly homopolymer length errors. To keep these from fragmenting haplotype counts:
- any base below `--min-base-quality` is replaced with the reference base;
- an insertion is kept only if all of its bases pass the quality threshold;
- 1-bp insertions and deletions inside (or extending) a reference homopolymer of at least `--homopolymer-min-length` bases are ignored, unless `--keep-homopolymer-indels` is set.

Basecalling output is written to `<output>/basecalling/` (`calls.bam`, `fastq/<sample>.fastq`), and re-running the command reuses it. `basecalling_summary.csv` holds reads, bases, mean length, mean Q and N50 per sample.

An existing Nanopore BAM (for example from dorado with `--reference`, sorted and indexed) can be analyzed directly. This writes the haplotype tables, QC CSVs and HTML report, and can be re-run later with `clonearmy report`:

```bash
clonearmy process-bam sample.bam reference.fasta --platform ont -o results_ont
```

#### QC and single-mutant coverage tables

`run`, `nanopore` and `process-bam` now also print and save (turn off with `--no-qc`):

- **Read filtering**: input reads, unmapped, mapped, supplementary (chimeras on Nanopore), low MAPQ, reads analysed, **full-length reads** (alignment spans the whole reference), full-length reads in haplotypes passing `--min-read-count`, read length, mean read Q, median alignment identity, and per-position depth (min/median/mean). Counts are reads for Nanopore and read pairs for Illumina.
- **Mutations per read**: % of reads with 0, 1, 2, 3, 4 and 5+ mutations, for full-length reads and all reads.
- **Single mutants by minimum read depth**: for each threshold in `--depth-thresholds` (default `10,100,1000,10000`), the number of haplotypes and reads kept, wild-type reads, and how many distinct single-mutation variants are supported by at least that many reads: SNVs (and % of the 3 × length possible), positions with an SNV, single deletions and insertions. For coding references (length divisible by 3, no internal stops) it adds synonymous SNVs, distinct missense amino-acid changes and nonsense codons, as % of those reachable by one nucleotide change. Computed separately for full-length reads and all reads.
- **Single variants**: every single-mutation variant with its read and full-length read counts, codon and amino-acid change. The console shows the top variants; the full table is in `qc_single_variants.csv` and the HTML report.

A read, or an Illumina read pair, is full length when its alignment starts within `--full-length-tolerance` bases (default 5) of the reference start and ends within that distance of the end, with no gap between mates. This matters for Nanopore: positions a truncated read does not cover are filled with the reference base, so without this check a partial read looks like a full-length wild-type or single-mutant read.

#### Reporting on already processed results

`clonearmy report` rebuilds the HTML report and all QC tables from an existing output directory (from `run`, `nanopore`, `process-bam` or `compare`) without re-aligning:

```bash
# Uses the settings stored by the run; per-sample BAMs are read for read statistics
clonearmy report results_ont reference.fasta

# Different depth thresholds, report written elsewhere
clonearmy report results_ont reference.fasta -d 10,100,1000,10000 -o report_v2

# Results from older CloneArmy versions: give the settings used for that run,
# haplotypes are rebuilt from the BAMs and cached next to them
clonearmy report results reference.fasta -q 25 -Q 20

# Fast: existing haplotype CSVs only, no BAM access
clonearmy report results reference.fasta --no-bam
```

Each run now writes `{sample}_haplotypes_all.csv.gz` (all haplotypes, no read-count filter) and `{sample}_parameters.json` (settings) next to the BAM. `report` reuses the unfiltered table when its settings match; if you change a setting that affects haplotype calling (base or mapping quality, homopolymer handling, indel size, full-length tolerance), or the table is missing, it rebuilds it from the BAM. Unset options fall back to the stored settings, then to platform defaults, and the source of each setting is shown in a table. Without a BAM (`--no-bam`, or older results whose BAM is gone), only the filtered `{sample}_haplotypes.csv` can be used, so per-read full-length status is unavailable and depth thresholds below the original `--min-read-count` are skipped.

| Option | Default | Description |
|---|---|---|
| `--depth-thresholds`, `-d` | `10,100,1000,10000` | Minimum read depths for the single-mutant table |
| `--sample`, `-s` | all | Only report these samples (repeatable) |
| `--bam` | none | Aligned BAMs stored outside the results directory (repeatable) |
| `--platform`, `-p` | stored / detected | `illumina` or `ont` |
| `--min-read-count`, `-r` | stored / 10 | Haplotype filter for the haplotype table and plots |
| `--reprocess` | off | Rebuild haplotype tables from BAMs even if cached tables match |
| `--no-bam` | off | Use existing haplotype CSVs only |
| `--no-read-stats` | off | Skip read length, Q, identity and depth statistics |

### Output Examples

#### Sample Analysis Results
```
╒════════════════╤══════════╤════════════╤══════════════╕
│ Sample         │ Reads    │ Haplotypes │ Mutations    │
╞════════════════╪══════════╪════════════╪══════════════╡
│ sample1        │ 10000    │ 45         │ 2.3 avg      │
│ sample2        │ 12000    │ 52         │ 1.8 avg      │
╘════════════════╧══════════╧════════════╧══════════════╛
```

#### Comparative Analysis Results
```
╒══════════╤════════════╤════════════╤═══════════╤═══════════╕
│ Position │ Sample 1 % │ Sample 2 % │ P-value   │ FDR       │
╞══════════╪════════════╪════════════╪═══════════╪═══════════╡
│ 123 A>T  │ 45.2      │ 12.3       │ 0.001     │ 0.003     │
│ 456 G>C  │ 33.1      │ 28.9       │ 0.042     │ 0.063     │
╘══════════╧════════════╧════════════╧═══════════╧═══════════╛
```

### Python API

```python
from pathlib import Path
from clone_army.processor import AmpliconProcessor
from clone_army.comparison import run_comparative_analysis

# Initialize processor with automatic downsampling
processor = AmpliconProcessor(
    reference_path="reference.fasta",
    min_base_quality=20,
    min_mapping_quality=30,
    min_read_count=10,
    max_file_size=100_000_000  # 100MB target size
)

# Process samples
results1 = processor.process_sample(
    fastq_r1="sample1_R1.fastq.gz",
    fastq_r2="sample1_R2.fastq.gz",
    output_dir="results/sample1",
    threads=4
)

results2 = processor.process_sample(
    fastq_r1="sample2_R1.fastq.gz",
    fastq_r2="sample2_R2.fastq.gz",
    output_dir="results/sample2",
    threads=4
)

# Perform comparative analysis
comparison_results = run_comparative_analysis(
    results1=results1,
    results2=results2,
    reference_seq=ref_seq,
    output_path="comparison_results.csv",
    full_length_only=False
)
```

## Output Files

### Single Sample Analysis
- Sorted BAM file with alignments
- `{sample}_haplotypes.csv` containing:
  - Sequence in reference coordinates (substitutions in lowercase, deletions as `-`)
  - Insertions relative to the reference, e.g. `250_251insTA` (inserted between reference positions 250 and 251; multiple separated by `;`). Insertions are part of the haplotype identity.
  - Read count
  - Frequency
  - Number of mutations
  - Full-length status (`is_full_length`, sequence based) and `full_length_count` (reads in the haplotype whose alignment spans the reference)
  - Quality metrics
- `{sample}_haplotypes_all.csv.gz`: the same columns for all haplotypes (no `--min-read-count` filter)
- `{sample}_parameters.json`: settings used, read by `clonearmy report`
- In the output directory:
  - `analysis_summary.csv`: the summary table
  - `qc_read_filtering.csv`, `qc_mutation_load.csv`, `qc_depth_thresholds.csv`, `qc_single_variants.csv` (see [QC tables](#qc-and-single-mutant-coverage-tables))
- Interactive HTML report with:
  - Settings, read filtering, mutations per read, and single-mutant coverage by minimum read depth (table and plot)
  - Summary statistics
  - Mutation frequency plots
  - Position-based mutation diversity plots
  - Mutation spectrum analysis
- Console output with summary statistics

### Comparative Analysis
- `comparison_results.csv` with statistical comparisons:
  - Mutation positions and types
  - Frequencies in each sample
  - Statistical significance (p-values)
  - FDR-corrected p-values
- Interactive HTML plots:
  - Mutation frequency comparison
  - Position-based mutation diversity
- Console output with significant mutations