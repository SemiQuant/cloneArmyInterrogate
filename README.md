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

```bash
pip install clonearmy
```

### Requirements

- Python ≥ 3.8
- BWA (must be installed and available in PATH)
- Samtools (must be installed and available in PATH)
- Seqtk (must be installed and available in PATH)

You can install the required tools using conda:
```bash
conda install -c bioconda bwa samtools seqtk
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

Additional requirements:

```bash
micromamba install -c bioconda minimap2 samtools
```

- [dorado](https://github.com/nanoporetech/dorado) must be on `PATH` (or passed with `--dorado-bin`). SUP basecalling needs a GPU: Apple Silicon (Metal) or NVIDIA (CUDA). SUP on CPU is impractically slow.

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

Nanopore error handling: even SUP reads carry a few errors per read, mostly homopolymer length errors. To keep these from fragmenting haplotype counts:
- any base below `--min-base-quality` is replaced with the reference base;
- an insertion is kept only if all of its bases pass the quality threshold;
- 1-bp insertions and deletions inside (or extending) a reference homopolymer of at least `--homopolymer-min-length` bases are ignored, unless `--keep-homopolymer-indels` is set.

Basecalling output is written to `<output>/basecalling/` (`calls.bam`, `fastq/<sample>.fastq`), and re-running the command reuses it. `basecalling_summary.csv` holds reads, bases, mean length, mean Q and N50 per sample.

An existing Nanopore BAM (for example from dorado with `--reference`, sorted and indexed) can be analyzed directly:

```bash
clonearmy process-bam sample.bam reference.fasta --platform ont
```

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
  - Full-length status
  - Quality metrics
- Interactive HTML report with:
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