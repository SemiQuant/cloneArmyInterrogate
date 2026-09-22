import sys
from pathlib import Path
from typing import Optional
import time
import logging

import click
from rich.console import Console
from rich.table import Table
from rich.progress import Progress, SpinnerColumn, TextColumn
from rich.panel import Panel
from rich import print as rprint
from Bio import SeqIO

from . import process_samples, summarize_results, validate_input, __version__
from .utils import process_nanopore_samples, validate_nanopore_input
from .report import generate_report
from .comparison import run_comparative_analysis
from .processor import AmpliconProcessor

console = Console()

def load_reference_sequence(reference_path: Path) -> str:
    """Load the reference sequence from a FASTA file."""
    try:
        with open(reference_path) as handle:
            record = next(SeqIO.parse(handle, "fasta"))
            return str(record.seq)
    except Exception as e:
        console.print(f"[bold red]Error loading reference sequence:[/] {str(e)}")
        sys.exit(1)

def print_dataframe(df, title: str):
    """Print a DataFrame as a rich table."""
    table = Table(title=title)
    for col in df.columns:
        table.add_column(str(col))
    for _, row in df.iterrows():
        table.add_row(*[str(x) for x in row])
    console.print(table)

def print_version(ctx, param, value):
    """Print version and exit."""
    if not value or ctx.resilient_parsing:
        return
    console.print(f"CloneArmy version [bold cyan]{__version__}[/]")
    ctx.exit()

@click.group()
@click.option('--version', is_flag=True, callback=print_version,
              expose_value=False, is_eager=True, help="Show version and exit.")
@click.option('--debug/--no-debug', default=False, help="Enable debug logging.")
def cli(debug: bool):
    """
    CloneArmy: Analyze haplotypes from Illumina paired-end or Oxford Nanopore
    amplicon sequencing.
    
    This tool processes FASTQ files (Illumina) or raw POD5 data basecalled with
    dorado SUP (Nanopore) to identify and quantify sequence variants and
    haplotypes in amplicon sequencing data.
    """
    if debug:
        logging.getLogger().setLevel(logging.DEBUG)

@cli.command()
@click.argument('fastq_dir', type=click.Path(exists=True))
@click.argument('reference', type=click.Path(exists=True))
@click.option('--threads', '-t', default=8, help='Number of threads to use')
@click.option('--output', '-o', type=click.Path(), help='Output directory')
@click.option('--min-base-quality', '-q', default=20, 
              help='Minimum base quality score')
@click.option('--min-mapping-quality', '-Q', default=30,
              help='Minimum mapping quality score')
@click.option('--min-read-count', '-r', default=10,
              help='Minimum number of reads to consider a haplotype')
@click.option('--max-file-size', '-m', default=10_000_000_000,
              help='Maximum file size in bytes (default: 10GB)')
@click.option('--report/--no-report', default=True,
              help='Generate HTML report')
@click.option('--bed', '-b', type=click.Path(exists=True),
              help='BED file for comparing indel positions')
@click.option('--max-indel-size', '-i', default=50,
              help='Maximum size of indels to consider as small indels')
def run(fastq_dir: str, reference: str, threads: int, output: Optional[str],
        min_base_quality: int, min_mapping_quality: int, min_read_count: int,
        max_file_size: int, report: bool, bed: Optional[str], max_indel_size: int):
    """Process FASTQ files and analyze mutations.

    FASTQ_DIR: Directory containing paired FASTQ files (_R1.fastq.gz and _R2.fastq.gz)
    REFERENCE: Reference sequence in FASTA format
    """
    start_time = time.time()
    
    # Validate input files
    fastq_dir = Path(fastq_dir)
    reference = Path(reference)
    output_dir = Path(output) if output else fastq_dir / 'results'
    bed_path = Path(bed) if bed else None
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console
    ) as progress:
        task = progress.add_task("Validating input files...", total=None)
        try:
            validate_input(fastq_dir, reference)
        except Exception as e:
            console.print(f"[bold red]Error:[/] {str(e)}")
            sys.exit(1)
        progress.update(task, completed=True)
        
        # Load reference sequence
        task = progress.add_task("Loading reference sequence...", total=None)
        ref_seq = load_reference_sequence(reference)
        progress.update(task, completed=True)
        
        # Process samples
        task = progress.add_task("Processing samples...", total=None)
        try:
            results, processor = process_samples(
                fastq_dir=fastq_dir,
                reference=reference,
                output_dir=output_dir,
                threads=threads,
                min_base_quality=min_base_quality,
                min_mapping_quality=min_mapping_quality,
                min_read_count=min_read_count,
                max_file_size=max_file_size,
                bed_path=bed_path,
                max_indel_size=max_indel_size
            )
        except Exception as e:
            console.print(f"[bold red]Error processing samples:[/] {str(e)}")
            sys.exit(1)
        progress.update(task, completed=True)
        
        # Generate summary
        task = progress.add_task("Generating summary...", total=None)
        summary = summarize_results(results, processor)
        progress.update(task, completed=True)
        
        # Generate report
        if report:
            task = progress.add_task("Generating report...", total=None)
            try:
                report_path = output_dir / 'report.html'
                generate_report(results, summary, report_path, ref_seq)
                progress.update(task, completed=True)
            except Exception as e:
                console.print(f"[bold red]Error generating report:[/] {str(e)}")
                sys.exit(1)
    
    # Print summary table
    if summary is not None and not summary.empty:
        table = Table(title="Analysis Summary")
        for col in summary.columns:
            table.add_column(col)
        for _, row in summary.iterrows():
            table.add_row(*[str(x) for x in row])
        console.print(table)
    
    elapsed = time.time() - start_time
    console.print(f"\n[green]Analysis completed in {elapsed:.1f} seconds[/]")

@cli.command()
@click.argument('pod5_dir', type=click.Path(exists=True, file_okay=False))
@click.argument('reference', type=click.Path(exists=True))
@click.option('--threads', '-t', default=8, help='Number of threads to use')
@click.option('--output', '-o', type=click.Path(), help='Output directory')
@click.option('--kit-name', '-k', default=None,
              help='Barcoding kit for demultiplexing (e.g. SQK-NBD114-24). Omit for a single sample')
@click.option('--model', default='sup', show_default=True,
              help='dorado model; "sup" selects the latest super-accurate model for the data')
@click.option('--device', '-x', default=None,
              help='dorado device (e.g. metal, cuda:all, cpu). Default: dorado auto-detect')
@click.option('--dorado-bin', type=click.Path(), default=None,
              help='Path to the dorado executable (default: dorado on PATH)')
@click.option('--models-dir', type=click.Path(), default=None,
              help='Directory to cache downloaded dorado models')
@click.option('--min-qscore', default=10, show_default=True,
              help='Minimum mean read Q-score kept by dorado')
@click.option('--barcode-map', type=click.Path(exists=True), default=None,
              help='CSV mapping barcode to sample name (columns: barcode,sample)')
@click.option('--force-basecall', is_flag=True,
              help='Redo basecalling and demultiplexing even if outputs exist')
@click.option('--min-base-quality', '-q', default=20, show_default=True,
              help='Bases below this quality are masked to the reference')
@click.option('--min-mapping-quality', '-Q', default=20, show_default=True,
              help='Minimum mapping quality score')
@click.option('--min-read-count', '-r', default=10,
              help='Minimum number of reads to consider a haplotype')
@click.option('--minimap2-preset', default='lr:hq', show_default=True,
              help='minimap2 preset (lr:hq for SUP reads, map-ont for older data)')
@click.option('--ignore-homopolymer-indels/--keep-homopolymer-indels', default=True, show_default=True,
              help='Ignore 1-bp indels in reference homopolymer runs (common ONT error)')
@click.option('--homopolymer-min-length', default=3, show_default=True,
              help='Minimum reference run length treated as a homopolymer')
@click.option('--report/--no-report', default=True,
              help='Generate HTML report')
@click.option('--bed', '-b', type=click.Path(exists=True),
              help='BED file for comparing indel positions')
@click.option('--max-indel-size', '-i', default=50,
              help='Maximum size of indels to consider as small indels')
def nanopore(pod5_dir: str, reference: str, threads: int, output: Optional[str],
             kit_name: Optional[str], model: str, device: Optional[str],
             dorado_bin: Optional[str], models_dir: Optional[str], min_qscore: int,
             barcode_map: Optional[str], force_basecall: bool, min_base_quality: int,
             min_mapping_quality: int, min_read_count: int, minimap2_preset: str,
             ignore_homopolymer_indels: bool, homopolymer_min_length: int,
             report: bool, bed: Optional[str], max_indel_size: int):
    """Basecall Nanopore POD5 data with dorado SUP and analyze haplotypes.

    POD5_DIR: Directory containing raw POD5 files (searched recursively)
    REFERENCE: Reference sequence in FASTA format
    """
    start_time = time.time()

    pod5_dir = Path(pod5_dir)
    reference = Path(reference)
    output_dir = Path(output) if output else pod5_dir / 'results'
    bed_path = Path(bed) if bed else None
    existing_calls = output_dir / 'basecalling' / 'calls.bam'
    resuming = existing_calls.exists() and not force_basecall

    warnings = validate_nanopore_input(pod5_dir, reference, dorado_bin,
                                       require_dorado=not resuming)
    if warnings:
        table = Table(title="Input Validation Problems")
        table.add_column("Problem", style="red")
        for w in warnings:
            table.add_row(w)
        console.print(table)
        if any(not w.startswith('Only FAST5') for w in warnings):
            sys.exit(1)

    params = Table(title="Nanopore Run Settings")
    params.add_column("Setting")
    params.add_column("Value")
    for key, value in [
        ('POD5 directory', pod5_dir), ('Reference', reference), ('Output', output_dir),
        ('Model', model), ('Kit', kit_name or 'none (single sample)'),
        ('Device', device or 'auto'), ('Min read Q', min_qscore),
        ('Min base Q (masking)', min_base_quality), ('Min MAPQ', min_mapping_quality),
        ('minimap2 preset', minimap2_preset),
        ('Homopolymer indels', f"ignored (run >= {homopolymer_min_length})"
                               if ignore_homopolymer_indels else 'kept'),
        ('Resume from existing basecalls', 'yes' if resuming else 'no'),
    ]:
        params.add_row(key, str(value))
    console.print(params)

    ref_seq = load_reference_sequence(reference)

    try:
        results, processor, basecall_stats = process_nanopore_samples(
            pod5_dir=pod5_dir,
            reference=reference,
            output_dir=output_dir,
            threads=threads,
            kit_name=kit_name,
            model=model,
            device=device,
            dorado_bin=dorado_bin,
            models_dir=models_dir,
            min_qscore=min_qscore,
            barcode_map=barcode_map,
            force_basecall=force_basecall,
            min_base_quality=min_base_quality,
            min_mapping_quality=min_mapping_quality,
            min_read_count=min_read_count,
            bed_path=bed_path,
            max_indel_size=max_indel_size,
            minimap2_preset=minimap2_preset,
            ignore_homopolymer_indels=ignore_homopolymer_indels,
            homopolymer_min_length=homopolymer_min_length
        )
    except Exception as e:
        console.print(f"[bold red]Error processing Nanopore samples:[/] {str(e)}")
        sys.exit(1)

    if not basecall_stats.empty:
        basecall_stats.to_csv(output_dir / 'basecalling_summary.csv', index=False)
        print_dataframe(basecall_stats, "Basecalling Summary")

    summary = summarize_results(results, processor)

    if report:
        try:
            report_path = output_dir / 'report.html'
            generate_report(results, summary, report_path, ref_seq)
            console.print(f"Report written to {report_path}")
        except Exception as e:
            console.print(f"[bold red]Error generating report:[/] {str(e)}")
            sys.exit(1)

    if summary is not None and not summary.empty:
        print_dataframe(summary, "Analysis Summary")

    elapsed = time.time() - start_time
    console.print(f"\n[green]Nanopore analysis completed in {elapsed:.1f} seconds[/]")

@cli.command()
@click.argument('fastq_dir1', type=click.Path(exists=True))
@click.argument('fastq_dir2', type=click.Path(exists=True))
@click.argument('reference', type=click.Path(exists=True))
@click.option('--threads', '-t', default=8, help='Number of threads to use')
@click.option('--output', '-o', type=click.Path(), help='Output directory')
@click.option('--min-base-quality', '-q', default=20, 
              help='Minimum base quality score')
@click.option('--min-mapping-quality', '-Q', default=30,
              help='Minimum mapping quality score')
@click.option('--min-read-count', '-r', default=10,
              help='Minimum number of reads to consider a haplotype')
@click.option('--max-file-size', '-m', default=10_000_000_000,
              help='Maximum file size in bytes (default: 10GB)')
@click.option('--full-length-only', '-f', is_flag=True,
              help='Only consider sequences that cover the entire reference')
@click.option('--bed', '-b', type=click.Path(exists=True),
              help='BED file for comparing indel positions')
@click.option('--max-indel-size', '-i', default=50,
              help='Maximum size of indels to consider as small indels')
def compare(fastq_dir1: str, fastq_dir2: str, reference: str, threads: int, 
           output: Optional[str], min_base_quality: int, min_mapping_quality: int,
           min_read_count: int, max_file_size: int, full_length_only: bool,
           bed: Optional[str], max_indel_size: int):
    """Compare two sets of FASTQ files and analyze differences.

    FASTQ_DIR1: First directory containing paired FASTQ files
    FASTQ_DIR2: Second directory containing paired FASTQ files
    REFERENCE: Reference sequence in FASTA format
    """
    start_time = time.time()
    
    # Validate input files
    fastq_dir1 = Path(fastq_dir1)
    fastq_dir2 = Path(fastq_dir2)
    reference = Path(reference)
    output_dir = Path(output) if output else Path('comparison_results')
    bed_path = Path(bed) if bed else None
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console
    ) as progress:
        task = progress.add_task("Validating input files...", total=None)
        try:
            validate_input(fastq_dir1, reference)
            validate_input(fastq_dir2, reference)
        except Exception as e:
            console.print(f"[bold red]Error:[/] {str(e)}")
            sys.exit(1)
        progress.update(task, completed=True)
        
        # Load reference sequence
        task = progress.add_task("Loading reference sequence...", total=None)
        ref_seq = load_reference_sequence(reference)
        progress.update(task, completed=True)
        
        # Process first set of samples
        task = progress.add_task("Processing first set of samples...", total=None)
        try:
            results1, processor1 = process_samples(
                fastq_dir=fastq_dir1,
                reference=reference,
                output_dir=output_dir / 'set1',
                threads=threads,
                min_base_quality=min_base_quality,
                min_mapping_quality=min_mapping_quality,
                min_read_count=min_read_count,
                max_file_size=max_file_size,
                bed_path=bed_path,
                max_indel_size=max_indel_size
            )
        except Exception as e:
            console.print(f"[bold red]Error processing first set of samples:[/] {str(e)}")
            sys.exit(1)
        progress.update(task, completed=True)
        
        # Process second set of samples
        task = progress.add_task("Processing second set of samples...", total=None)
        try:
            results2, processor2 = process_samples(
                fastq_dir=fastq_dir2,
                reference=reference,
                output_dir=output_dir / 'set2',
                threads=threads,
                min_base_quality=min_base_quality,
                min_mapping_quality=min_mapping_quality,
                min_read_count=min_read_count,
                max_file_size=max_file_size,
                bed_path=bed_path,
                max_indel_size=max_indel_size
            )
        except Exception as e:
            console.print(f"[bold red]Error processing second set of samples:[/] {str(e)}")
            sys.exit(1)
        progress.update(task, completed=True)
        
        # Run comparative analysis
        task = progress.add_task("Running comparative analysis...", total=None)
        try:
            comparison_results = run_comparative_analysis(
                results1, results2, output_dir,
                full_length_only=full_length_only
            )
        except Exception as e:
            console.print(f"[bold red]Error in comparative analysis:[/] {str(e)}")
            sys.exit(1)
        progress.update(task, completed=True)
        
        # Generate summaries
        task = progress.add_task("Generating summaries...", total=None)
        summary1 = summarize_results(results1, processor1)
        summary2 = summarize_results(results2, processor2)
        progress.update(task, completed=True)
    
    elapsed = time.time() - start_time
    console.print(f"\n[green]Comparison completed in {elapsed:.1f} seconds[/]")

@cli.command()
@click.argument('bam_file', type=click.Path(exists=True))
@click.argument('reference', type=click.Path(exists=True))
@click.option('--threads', '-t', default=8, help='Number of threads to use')
@click.option('--output', '-o', type=click.Path(), help='Output directory')
@click.option('--min-base-quality', '-q', default=20, 
             help='Minimum base quality score')
@click.option('--min-mapping-quality', '-Q', default=30,
             help='Minimum mapping quality score')
@click.option('--min-read-count', '-r', default=10,
             help='Minimum number of reads to consider a haplotype')
@click.option('--bed', '-b', type=click.Path(exists=True),
             help='BED file for comparing indel positions')
@click.option('--max-indel-size', '-i', default=50,
             help='Maximum size of indels to consider as small indels')
@click.option('--platform', '-p', type=click.Choice(['illumina', 'ont'], case_sensitive=False),
             default='illumina', show_default=True,
             help='Sequencing platform: illumina (paired-end) or ont (Nanopore single reads)')
@click.option('--ignore-homopolymer-indels/--keep-homopolymer-indels', default=True, show_default=True,
             help='(ont) Ignore 1-bp indels in reference homopolymer runs')
@click.option('--homopolymer-min-length', default=3, show_default=True,
             help='(ont) Minimum reference run length treated as a homopolymer')
def process_bam(bam_file: str, reference: str, threads: int, output: Optional[str],
                min_base_quality: int, min_mapping_quality: int, min_read_count: int,
                bed: Optional[str], max_indel_size: int, platform: str,
                ignore_homopolymer_indels: bool, homopolymer_min_length: int):
    """Process an existing BAM file for haplotype analysis.
    
    BAM_FILE: Path to the input BAM file
    REFERENCE: Path to reference FASTA file
    """
    console = Console()
    
    # Set up output directory
    bam_path = Path(bam_file)
    output_dir = Path(output) if output else bam_path.parent / 'results'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Set up bed path
    bed_path = Path(bed) if bed else None
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console
    ) as progress:
        task = progress.add_task("Validating input files...", total=None)
        try:
            # Check if BAM file exists and is valid
            if not bam_path.exists():
                raise FileNotFoundError(f"BAM file not found: {bam_path}")
            
            # Check if reference exists
            ref_path = Path(reference)
            if not ref_path.exists():
                raise FileNotFoundError(f"Reference file not found: {ref_path}")
                
        except Exception as e:
            console.print(f"[bold red]Error:[/] {str(e)}")
            sys.exit(1)
        progress.update(task, completed=True)
        
        # Initialize processor
        task = progress.add_task("Initializing processor...", total=None)
        try:
            processor = AmpliconProcessor(
                reference_path=reference,
                bed_path=bed_path,
                min_base_quality=min_base_quality,
                min_mapping_quality=min_mapping_quality,
                min_read_count=min_read_count,
                max_indel_size=max_indel_size,
                platform=platform,
                ignore_homopolymer_indels=ignore_homopolymer_indels,
                homopolymer_min_length=homopolymer_min_length
            )
        except Exception as e:
            console.print(f"[bold red]Error initializing processor:[/] {str(e)}")
            sys.exit(1)
        progress.update(task, completed=True)
        
        # Process BAM file
        task = progress.add_task("Processing BAM file...", total=None)
        try:
            results = processor.process_bam(
                bam_path=bam_path,
                output_dir=output_dir,
                threads=threads
            )
            
            if results.empty:
                console.print("[yellow]Warning:[/] No results found in BAM file")
            else:
                console.print(f"[green]Success![/] Results saved to {output_dir}")
                
        except Exception as e:
            console.print(f"[bold red]Error processing BAM file:[/] {str(e)}")
            sys.exit(1)
        progress.update(task, completed=True)

if __name__ == '__main__':
    cli()