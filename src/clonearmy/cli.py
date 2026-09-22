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

import math

import pandas as pd

from . import process_samples, summarize_results, validate_input, __version__
from .utils import process_nanopore_samples, validate_nanopore_input
from .report import generate_report
from .comparison import run_comparative_analysis
from .processor import AmpliconProcessor
from .qc import (DEFAULT_DEPTHS_STR, QCResult, build_qc, discover_samples, load_haplotype_table,
                 parse_depths, resolve_parameters)

console = Console()

def _depths_callback(ctx, param, value):
    return parse_depths(value)

def _format_cell(value) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return '–'
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, int) or (hasattr(value, 'dtype') and 'int' in str(value.dtype)):
        return f"{int(value):,}"
    if isinstance(value, float):
        return f"{value:,.2f}"
    return str(value)

SUMMARY_COLUMNS = {
    'sample': 'Sample', 'total_reads': 'Reads', 'unique_haplotypes': 'Haplotypes',
    'unique_single_mut_haplotypes': 'Single-mut haps',
    'full_length_reads': 'Full length', 'full_length_percent': 'FL %',
    'avg_mutations': 'Avg muts', 'max_frequency': 'Max freq %',
}

def print_table(df: pd.DataFrame, title: str, columns: Optional[dict] = None):
    """Print selected DataFrame columns ({column: header}) as a rich table with formatted numbers."""
    if df is None or df.empty:
        return
    columns = columns or {c: c for c in df.columns}
    columns = {c: h for c, h in columns.items() if c in df.columns}
    table = Table(title=title)
    for header in columns.values():
        table.add_column(header, justify='left' if header in ('Sample', 'Reference', 'Reads used') else 'right')
    for _, row in df.iterrows():
        table.add_row(*[_format_cell(row[c]) for c in columns])
    console.print(table)

def print_qc(qc: Optional[QCResult], min_read_count: int, platform: str = 'illumina',
             show_basecalling: bool = True):
    """Print read filtering, mutation load and depth-threshold tables."""
    if qc is None:
        return
    multi_ref = (not qc.depth_thresholds.empty and qc.depth_thresholds['reference'].nunique() > 1)
    ref_col = {'reference': 'Reference'} if multi_ref else {}

    unit = 'read pairs' if platform == 'illumina' else 'reads'
    if show_basecalling:
        print_table(qc.basecalling, "Basecalling Summary")
    print_table(qc.read_filtering, f"Read Filtering (counts are {unit})", {
        'sample': 'Sample', **ref_col, 'input_reads': 'Input', 'unmapped': 'Unmapped',
        'mapped': 'Mapped', 'supplementary': 'Suppl.', 'low_mapq': 'Low MAPQ',
        'analysed': 'Analysed', 'full_length': 'Full length', 'full_length_pct': 'FL %',
        f'full_length_in_haplotypes_ge{min_read_count}': f'FL in haps >={min_read_count}',
        'median_read_length': 'Median len', 'mean_read_q': 'Mean Q',
        'median_identity_pct': 'Identity %', 'depth_median': 'Depth med', 'depth_min': 'Depth min',
    })

    print_table(qc.mutation_load, "Mutations per Read (% of reads)", {
        'sample': 'Sample', **ref_col, 'read_set': 'Reads used', 'reads': 'Reads',
        'mean_mutations': 'Mean', 'pct_0_mut': '0', 'pct_1_mut': '1', 'pct_2_mut': '2',
        'pct_3_mut': '3', 'pct_4_mut': '4', 'pct_5plus_mut': '5+',
    })

    thresholds = qc.depth_thresholds
    if not thresholds.empty:
        has_fl = (thresholds['read_set'] == 'full_length').any()
        shown = thresholds[thresholds['read_set'] == ('full_length' if has_fl else 'all')]
        label = 'full-length reads' if has_fl else 'all reads'
        for (sample, ref), group in shown.groupby(['sample', 'reference'], sort=False):
            print_table(group, f"Single Mutants by Minimum Depth: {sample}"
                               f"{' / ' + ref if multi_ref else ''} ({label})", {
                'min_depth': 'Min depth', 'haplotypes': 'Haplotypes', 'reads_pct': 'Reads %',
                'wild_type_reads': 'WT reads',
                'single_mut_variants': 'Single muts', 'single_mut_reads': 'Single-mut reads',
                'single_mut_pct': 'Single-mut %',
                'single_snvs': 'SNVs', 'snv_pct_of_possible': 'SNV % of 3L',
                'snv_positions': 'Positions', 'position_pct': 'Pos %',
                'single_deletions': 'Del', 'single_insertions': 'Ins',
                'missense_aa_changes': 'Missense aa', 'missense_pct_of_reachable': 'Missense %',
                'nonsense_codons': 'Nonsense',
            })
        if has_fl and (thresholds['read_set'] == 'all').any():
            console.print("[dim]All-reads depth tables are in qc_depth_thresholds.csv and the HTML report.[/]")

    if not qc.single_variants.empty:
        sv = qc.single_variants
        sort_col = 'full_length_reads' if sv['full_length_reads'].notna().any() else 'reads'
        for sample, group in sv.groupby('sample', sort=False):
            top = group.sort_values(sort_col, ascending=False).head(15)
            print_table(top, f"Top single variants: {sample} (by {sort_col.replace('_', ' ')})", {
                'type': 'Type', 'position': 'Pos', 'ref': 'Ref', 'alt': 'Alt',
                'reads': 'Reads', 'full_length_reads': 'FL reads',
                'aa_change': 'AA', 'consequence': 'Consequence',
            })
        if any(len(g) > 15 for _, g in sv.groupby('sample')):
            console.print("[dim]All single variants are in qc_single_variants.csv.[/]")
    for note in qc.notes:
        console.print(f"[yellow]Note:[/] {note}")

def run_qc(output_dir: Path, processor: AmpliconProcessor, sample_names, depths,
           extra_bams=(), basecalling: Optional[pd.DataFrame] = None,
           read_stats: bool = True) -> Optional[QCResult]:
    """Build and write QC tables for samples just processed into output_dir."""
    samples = discover_samples(output_dir, extra_bams=extra_bams, only=sample_names)
    tables = {}
    for name, sample in samples.items():
        source = load_haplotype_table(sample, processor)
        if source is not None:
            tables[name] = source
    if not tables:
        return None
    qc = build_qc(samples, tables, processor, depths, read_stats=read_stats, basecalling=basecalling)
    for path in qc.write(output_dir):
        console.print(f"QC table written to {path}")
    return qc

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
@click.option('--full-length-tolerance', default=5, show_default=True,
              help='Bases a read pair may miss at either reference end and still count as full length')
@click.option('--depth-thresholds', '-d', default=DEFAULT_DEPTHS_STR, show_default=True,
              callback=_depths_callback,
              help='Comma-separated minimum read depths for the single-mutant coverage table')
@click.option('--qc/--no-qc', default=True, show_default=True,
              help='Compute read filtering, mutation load and depth-threshold tables')
def run(fastq_dir: str, reference: str, threads: int, output: Optional[str],
        min_base_quality: int, min_mapping_quality: int, min_read_count: int,
        max_file_size: int, report: bool, bed: Optional[str], max_indel_size: int,
        full_length_tolerance: int, depth_thresholds, qc: bool):
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
                max_indel_size=max_indel_size,
                full_length_tolerance=full_length_tolerance
            )
        except Exception as e:
            console.print(f"[bold red]Error processing samples:[/] {str(e)}")
            sys.exit(1)
        progress.update(task, completed=True)
        
        # Generate summary
        task = progress.add_task("Generating summary...", total=None)
        summary = summarize_results(results, processor)
        if summary is not None and not summary.empty:
            summary.to_csv(output_dir / 'analysis_summary.csv', index=False)
        progress.update(task, completed=True)

        qc_result = None
        if qc:
            task = progress.add_task("Computing QC tables...", total=None)
            qc_result = run_qc(output_dir, processor, list(results), depth_thresholds)
            progress.update(task, completed=True)
        
        # Generate report
        if report:
            task = progress.add_task("Generating report...", total=None)
            try:
                report_path = output_dir / 'report.html'
                generate_report(results, summary, report_path, ref_seq, qc=qc_result,
                                parameters=processor.processing_parameters())
                progress.update(task, completed=True)
            except Exception as e:
                console.print(f"[bold red]Error generating report:[/] {str(e)}")
                sys.exit(1)
    
    if summary is not None and not summary.empty:
        print_table(summary, "Analysis Summary", SUMMARY_COLUMNS)

    print_qc(qc_result, min_read_count)
    
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
@click.option('--full-length-tolerance', default=5, show_default=True,
              help='Bases a read may miss at either reference end and still count as full length')
@click.option('--depth-thresholds', '-d', default=DEFAULT_DEPTHS_STR, show_default=True,
              callback=_depths_callback,
              help='Comma-separated minimum read depths for the single-mutant coverage table')
@click.option('--qc/--no-qc', default=True, show_default=True,
              help='Compute read filtering, mutation load and depth-threshold tables')
def nanopore(pod5_dir: str, reference: str, threads: int, output: Optional[str],
             kit_name: Optional[str], model: str, device: Optional[str],
             dorado_bin: Optional[str], models_dir: Optional[str], min_qscore: int,
             barcode_map: Optional[str], force_basecall: bool, min_base_quality: int,
             min_mapping_quality: int, min_read_count: int, minimap2_preset: str,
             ignore_homopolymer_indels: bool, homopolymer_min_length: int,
             report: bool, bed: Optional[str], max_indel_size: int,
             full_length_tolerance: int, depth_thresholds, qc: bool):
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
        ('Full-length tolerance (bp)', full_length_tolerance),
        ('Depth thresholds', ','.join(str(d) for d in depth_thresholds)),
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
            homopolymer_min_length=homopolymer_min_length,
            full_length_tolerance=full_length_tolerance
        )
    except Exception as e:
        console.print(f"[bold red]Error processing Nanopore samples:[/] {str(e)}")
        sys.exit(1)

    if not basecall_stats.empty:
        basecall_stats.to_csv(output_dir / 'basecalling_summary.csv', index=False)
        print_dataframe(basecall_stats, "Basecalling Summary")

    summary = summarize_results(results, processor)
    if summary is not None and not summary.empty:
        summary.to_csv(output_dir / 'analysis_summary.csv', index=False)

    qc_result = (run_qc(output_dir, processor, list(results), depth_thresholds,
                        basecalling=basecall_stats) if qc else None)

    if report:
        try:
            report_path = output_dir / 'report.html'
            parameters = {**processor.processing_parameters(), 'dorado_model': model,
                          'kit_name': kit_name, 'min_qscore': min_qscore,
                          'minimap2_preset': minimap2_preset}
            generate_report(results, summary, report_path, ref_seq, qc=qc_result,
                            parameters=parameters)
            console.print(f"Report written to {report_path}")
        except Exception as e:
            console.print(f"[bold red]Error generating report:[/] {str(e)}")
            sys.exit(1)

    if summary is not None and not summary.empty:
        print_table(summary, "Analysis Summary", SUMMARY_COLUMNS)

    print_qc(qc_result, min_read_count, 'ont', show_basecalling=False)

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
@click.option('--full-length-tolerance', default=5, show_default=True,
             help='Bases a read may miss at either reference end and still count as full length')
@click.option('--depth-thresholds', '-d', default=DEFAULT_DEPTHS_STR, show_default=True,
             callback=_depths_callback,
             help='Comma-separated minimum read depths for the single-mutant coverage table')
@click.option('--qc/--no-qc', default=True, show_default=True,
              help='Compute read filtering, mutation load and depth-threshold tables')
@click.option('--report/--no-report', default=True,
              help='Generate HTML report')
def process_bam(bam_file: str, reference: str, threads: int, output: Optional[str],
                min_base_quality: int, min_mapping_quality: int, min_read_count: int,
                bed: Optional[str], max_indel_size: int, platform: str,
                ignore_homopolymer_indels: bool, homopolymer_min_length: int,
                full_length_tolerance: int, depth_thresholds, qc: bool, report: bool):
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
                homopolymer_min_length=homopolymer_min_length,
                full_length_tolerance=full_length_tolerance
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

        qc_result = None
        if qc and not results.empty:
            task = progress.add_task("Computing QC tables...", total=None)
            qc_result = run_qc(output_dir, processor, [bam_path.stem], depth_thresholds,
                               extra_bams=[bam_path])
            progress.update(task, completed=True)

        summary = None
        if not results.empty:
            summary = summarize_results({bam_path.stem: results}, processor)
            if summary is not None and not summary.empty:
                summary.to_csv(output_dir / 'analysis_summary.csv', index=False)

        if report and not results.empty:
            task = progress.add_task("Generating report...", total=None)
            try:
                report_path = output_dir / 'report.html'
                generate_report({bam_path.stem: results}, summary, report_path,
                                load_reference_sequence(Path(reference)), qc=qc_result,
                                parameters=processor.processing_parameters(bam_path=bam_path))
                console.print(f"Report written to {report_path}")
            except Exception as e:
                console.print(f"[bold red]Error generating report:[/] {str(e)}")
                sys.exit(1)
            progress.update(task, completed=True)

    if summary is not None and not summary.empty:
        print_table(summary, "Analysis Summary", SUMMARY_COLUMNS)
    print_qc(qc_result, min_read_count, platform.lower())

@cli.command()
@click.argument('results_dir', type=click.Path(exists=True, file_okay=False))
@click.argument('reference', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), default=None,
              help='Directory for the report and QC tables (default: RESULTS_DIR)')
@click.option('--bam', 'bams', multiple=True, type=click.Path(exists=True, dir_okay=False),
              help='Aligned BAM stored outside RESULTS_DIR, e.g. a process-bam input (repeatable)')
@click.option('--sample', '-s', 'sample_names', multiple=True,
              help='Only report these samples (repeatable; default: all found)')
@click.option('--platform', '-p', type=click.Choice(['illumina', 'ont'], case_sensitive=False),
              default=None, help='Default: stored run settings, else detected from the BAM header')
@click.option('--min-base-quality', '-q', type=int, default=None,
              help='Default: stored run settings, else 20')
@click.option('--min-mapping-quality', '-Q', type=int, default=None,
              help='Default: stored run settings, else 30 (illumina) / 20 (ont)')
@click.option('--min-read-count', '-r', type=int, default=None,
              help='Minimum reads per haplotype for the haplotype table and plots. Default: stored, else 10')
@click.option('--max-indel-size', '-i', type=int, default=None,
              help='Default: stored run settings, else 50')
@click.option('--ignore-homopolymer-indels/--keep-homopolymer-indels', default=None,
              help='(ont) Default: stored run settings, else ignore')
@click.option('--homopolymer-min-length', type=int, default=None,
              help='(ont) Default: stored run settings, else 3')
@click.option('--full-length-tolerance', type=int, default=None,
              help='Bases a read may miss at either reference end and still count as full length. '
                   'Default: stored run settings, else 5')
@click.option('--depth-thresholds', '-d', default=DEFAULT_DEPTHS_STR, show_default=True,
              callback=_depths_callback,
              help='Comma-separated minimum read depths for the single-mutant coverage table')
@click.option('--bed', '-b', type=click.Path(exists=True),
              help='BED file for comparing indel positions')
@click.option('--reprocess', is_flag=True,
              help='Rebuild haplotype tables from the BAMs even when cached tables match the settings')
@click.option('--no-bam', is_flag=True,
              help='Never read BAMs; use existing haplotype tables only (fast, fewer statistics)')
@click.option('--read-stats/--no-read-stats', default=True, show_default=True,
              help='Read length, Q, identity and depth statistics from the BAMs')
def report(results_dir: str, reference: str, output: Optional[str], bams, sample_names,
           platform: Optional[str], min_base_quality: Optional[int],
           min_mapping_quality: Optional[int], min_read_count: Optional[int],
           max_indel_size: Optional[int], ignore_homopolymer_indels: Optional[bool],
           homopolymer_min_length: Optional[int], full_length_tolerance: Optional[int],
           depth_thresholds, bed: Optional[str], reprocess: bool, no_bam: bool, read_stats: bool):
    """Generate the report and QC tables for already processed results.

    Works on output of `run`, `nanopore`, `process-bam` or `compare`. Unfiltered
    haplotype tables are reused when their settings match; otherwise (or for
    results from older versions) they are rebuilt from the sample BAMs and cached.

    \b
    RESULTS_DIR: Output directory of a previous run
    REFERENCE: Reference sequence in FASTA format used for that run
    """
    start_time = time.time()
    results_dir = Path(results_dir)
    output_dir = Path(output) if output else results_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    samples = discover_samples(results_dir, extra_bams=bams,
                               only=sample_names if sample_names else None)
    if not samples:
        console.print(f"[bold red]Error:[/] no BAMs or haplotype tables found in {results_dir}")
        sys.exit(1)

    params, origin = resolve_parameters({
        'platform': platform, 'min_base_quality': min_base_quality,
        'min_mapping_quality': min_mapping_quality, 'min_read_count': min_read_count,
        'max_indel_size': max_indel_size, 'ignore_homopolymer_indels': ignore_homopolymer_indels,
        'homopolymer_min_length': homopolymer_min_length,
        'full_length_tolerance': full_length_tolerance,
    }, samples)

    found = Table(title=f"Samples in {results_dir}")
    for col in ('Sample', 'BAM', 'Unfiltered table', 'Filtered table'):
        found.add_column(col)
    for name, s in samples.items():
        found.add_row(name, 'yes' if s.bam else '–', 'yes' if s.haplotypes_all else '–',
                      'yes' if s.haplotypes else '–')
    console.print(found)

    settings = Table(title="Report Settings")
    for col in ('Setting', 'Value', 'Source'):
        settings.add_column(col)
    for key, value in params.items():
        settings.add_row(key, str(value), origin[key])
    settings.add_row('depth_thresholds', ','.join(str(d) for d in depth_thresholds), 'command line')
    console.print(settings)

    try:
        processor = AmpliconProcessor(
            reference_path=reference,
            bed_path=Path(bed) if bed else None,
            min_base_quality=params['min_base_quality'],
            min_mapping_quality=params['min_mapping_quality'],
            min_read_count=params['min_read_count'],
            max_indel_size=params['max_indel_size'],
            platform=params['platform'],
            ignore_homopolymer_indels=params['ignore_homopolymer_indels'],
            homopolymer_min_length=params['homopolymer_min_length'],
            max_homopolymer_indel=params['max_homopolymer_indel'],
            full_length_tolerance=params['full_length_tolerance'],
            check_dependencies=False
        )
    except Exception as e:
        console.print(f"[bold red]Error initializing processor:[/] {str(e)}")
        sys.exit(1)

    tables = {}
    for name, sample in samples.items():
        try:
            source = load_haplotype_table(sample, processor, reprocess=reprocess, use_bam=not no_bam)
        except Exception as e:
            console.print(f"[bold red]Error loading {name}:[/] {str(e)}")
            continue
        if source is None:
            console.print(f"[yellow]Skipping {name}:[/] no haplotype table and no usable BAM")
            continue
        tables[name] = source

    if not tables:
        console.print("[bold red]Error:[/] no sample could be loaded")
        sys.exit(1)

    results = {name: processor.filter_haplotypes(source.table) for name, source in tables.items()}
    summary = summarize_results(results, processor)
    if summary is not None and not summary.empty:
        summary.to_csv(output_dir / 'analysis_summary.csv', index=False)

    basecalling = None
    basecalling_csv = results_dir / 'basecalling_summary.csv'
    if basecalling_csv.exists():
        basecalling = pd.read_csv(basecalling_csv)
        basecalling = basecalling[basecalling['sample'].isin(set(tables) | {s.stem for s in samples.values()})]

    qc_result = build_qc(samples, tables, processor, depth_thresholds,
                         read_stats=read_stats and not no_bam, basecalling=basecalling)
    for path in qc_result.write(output_dir):
        console.print(f"QC table written to {path}")

    report_path = output_dir / 'report.html'
    try:
        generate_report(results, summary, report_path, load_reference_sequence(Path(reference)),
                        qc=qc_result, parameters={**params,
                                                  'depth_thresholds': list(depth_thresholds)})
        console.print(f"Report written to {report_path}")
    except Exception as e:
        console.print(f"[bold red]Error generating report:[/] {str(e)}")
        sys.exit(1)

    print_table(qc_result.sources, "Haplotype Table Sources",
                {'sample': 'Sample', 'source': 'Source', 'bam': 'BAM'})
    if summary is not None and not summary.empty:
        print_table(summary, "Analysis Summary", SUMMARY_COLUMNS)
    print_qc(qc_result, params['min_read_count'], params['platform'])

    elapsed = time.time() - start_time
    console.print(f"\n[green]Report completed in {elapsed:.1f} seconds[/]")

if __name__ == '__main__':
    cli()