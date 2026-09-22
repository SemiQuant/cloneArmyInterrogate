from pathlib import Path
from typing import Union, Dict, List, Tuple, Optional
import logging
import pandas as pd
from Bio import SeqIO
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import partial
import multiprocessing
import click

from .processor import AmpliconProcessor, full_length_read_count
from . import nanopore

logger = logging.getLogger(__name__)

def process_samples(
    fastq_dir: Union[str, Path],
    reference: Union[str, Path],
    output_dir: Union[str, Path, None] = None,
    threads: int = 4,
    min_base_quality: int = 20,
    min_mapping_quality: int = 30,
    min_read_count: int = 10,
    max_file_size: int = 10_000_000_000,
    bed_path: Union[str, Path, None] = None,
    max_indel_size: int = 50,
    parallel_samples: int = None,
    full_length_tolerance: int = 5
) -> Tuple[Dict[str, pd.DataFrame], AmpliconProcessor]:
    """
    Process all samples in a directory.

    Args:
        fastq_dir: Directory containing FASTQ files
        reference: Path to reference FASTA file
        output_dir: Directory for output files (default: fastq_dir/results)
        threads: Number of threads to use per sample
        min_base_quality: Minimum base quality score
        min_mapping_quality: Minimum mapping quality score
        min_read_count: Minimum number of reads to consider a haplotype
        max_file_size: Maximum file size in bytes
        bed_path: Optional path to BED file for indel comparison
        max_indel_size: Maximum size of indels to consider as small indels
        parallel_samples: Number of samples to process in parallel (default: min(4, CPU count))
        full_length_tolerance: Bases a read pair may miss at either reference end and still be full length

    Returns:
        Tuple of (Dictionary mapping sample names to their results DataFrames, AmpliconProcessor)
    """
    fastq_dir = Path(fastq_dir)
    reference = Path(reference)
    output_dir = Path(output_dir) if output_dir else fastq_dir / 'results'
    output_dir.mkdir(parents=True, exist_ok=True)

    # Set parallel processing parameters
    if parallel_samples is None:
        parallel_samples = min(4, multiprocessing.cpu_count())
    
    # Adjust threads per sample based on parallel processing
    threads_per_sample = max(1, threads // parallel_samples)

    # Initialize processor
    processor = AmpliconProcessor(
        reference_path=reference,
        bed_path=bed_path,
        min_base_quality=min_base_quality,
        min_mapping_quality=min_mapping_quality,
        min_read_count=min_read_count,
        max_file_size=max_file_size,
        max_indel_size=max_indel_size,
        full_length_tolerance=full_length_tolerance
    )

    def process_single_sample(r1_file: Path, processor: AmpliconProcessor) -> Tuple[str, pd.DataFrame]:
        """Process a single sample and return its name and results."""
        try:
            r2_file = r1_file.parent / r1_file.name.replace('_R1_', '_R2_')
            if not r2_file.exists():
                logger.warning(f"No R2 file found for {r1_file.name}")
                return None

            sample_name = r1_file.name.split('_R1_')[0]
            logger.info(f"Processing sample: {sample_name}")

            result = processor.process_sample(
                fastq_r1=r1_file,
                fastq_r2=r2_file,
                output_dir=output_dir / sample_name,
                threads=threads_per_sample
            )
            return sample_name, result

        except Exception as e:
            logger.error(f"Error processing sample {r1_file.stem}: {str(e)}")
            return None

    # Get list of R1 files
    r1_files = sorted(fastq_dir.glob('*_R1_001.fastq*'))
    
    # Process samples in parallel
    results = {}
    with ThreadPoolExecutor(max_workers=parallel_samples) as executor:
        # Create futures for each sample
        future_to_sample = {
            executor.submit(process_single_sample, r1_file, processor): r1_file
            for r1_file in r1_files
        }

        # Process results as they complete
        with click.progressbar(length=len(r1_files), 
                             label='Processing samples') as bar:
            for future in as_completed(future_to_sample):
                result = future.result()
                if result is not None:
                    sample_name, df = result
                    results[sample_name] = df
                bar.update(1)

    return results, processor

def load_barcode_map(barcode_map: Union[str, Path, None]) -> Dict[str, str]:
    """Load a barcode,sample CSV (header optional) into {barcode: sample}."""
    if not barcode_map:
        return {}
    df = pd.read_csv(barcode_map, header=None, dtype=str, comment='#')
    if df.shape[1] < 2:
        raise ValueError(f"Barcode map must have two columns (barcode,sample): {barcode_map}")
    if df.iloc[0, 0].strip().lower() == 'barcode':
        df = df.iloc[1:]
    mapping = {}
    for barcode, sample in zip(df.iloc[:, 0], df.iloc[:, 1]):
        barcode = str(barcode).strip().lower()
        if barcode.isdigit():
            barcode = f"barcode{int(barcode):02d}"
        mapping[barcode] = str(sample).strip()
    return mapping

def process_nanopore_samples(
    pod5_dir: Union[str, Path],
    reference: Union[str, Path],
    output_dir: Union[str, Path, None] = None,
    threads: int = 4,
    kit_name: Optional[str] = None,
    model: str = 'sup',
    device: Optional[str] = None,
    dorado_bin: Union[str, Path, None] = None,
    models_dir: Union[str, Path, None] = None,
    min_qscore: int = 10,
    barcode_map: Union[str, Path, None] = None,
    force_basecall: bool = False,
    min_base_quality: int = 20,
    min_mapping_quality: int = 20,
    min_read_count: int = 10,
    bed_path: Union[str, Path, None] = None,
    max_indel_size: int = 50,
    minimap2_preset: str = 'lr:hq',
    ignore_homopolymer_indels: bool = True,
    homopolymer_min_length: int = 3,
    parallel_samples: int = None,
    full_length_tolerance: int = 5
) -> Tuple[Dict[str, pd.DataFrame], AmpliconProcessor, pd.DataFrame]:
    """
    Basecall raw Nanopore POD5 data with dorado (SUP by default) and analyze haplotypes.

    Args:
        pod5_dir: Directory containing POD5 files (searched recursively)
        reference: Path to reference FASTA file
        output_dir: Directory for output files (default: pod5_dir/results)
        threads: Number of threads to use
        kit_name: Barcoding kit (e.g. SQK-NBD114-24); if None, all reads form one sample
        model: dorado model or model complex (default: sup)
        device: dorado device string (e.g. metal, cuda:all); dorado default if None
        dorado_bin: Path to dorado executable
        models_dir: Directory to cache downloaded dorado models
        min_qscore: Minimum mean read Q-score kept by dorado
        barcode_map: Optional CSV mapping barcode to sample name
        force_basecall: Redo basecalling/demux even if outputs exist
        min_base_quality: Bases below this quality are masked to the reference
        min_mapping_quality: Minimum mapping quality score
        min_read_count: Minimum number of reads to consider a haplotype
        bed_path: Optional path to BED file for indel comparison
        max_indel_size: Maximum size of indels to consider as small indels
        minimap2_preset: minimap2 preset for alignment
        ignore_homopolymer_indels: Ignore 1-bp indels in reference homopolymer runs
        homopolymer_min_length: Minimum run length treated as a homopolymer
        parallel_samples: Number of samples to process in parallel (default: min(4, CPU count))
        full_length_tolerance: Bases a read may miss at either reference end and still be full length

    Returns:
        Tuple of (results per sample, AmpliconProcessor, basecalling summary DataFrame)
    """
    pod5_dir = Path(pod5_dir)
    reference = Path(reference)
    output_dir = Path(output_dir) if output_dir else pod5_dir / 'results'
    output_dir.mkdir(parents=True, exist_ok=True)
    basecall_dir = output_dir / 'basecalling'

    calls_bam = nanopore.basecall(
        pod5_dir=pod5_dir,
        out_bam=basecall_dir / 'calls.bam',
        model=model,
        kit_name=kit_name,
        device=device,
        min_qscore=min_qscore,
        models_dir=models_dir,
        dorado_bin=dorado_bin,
        force=force_basecall
    )

    if kit_name:
        fastqs = nanopore.demux(calls_bam, basecall_dir, dorado_bin=dorado_bin,
                                threads=threads, force=force_basecall)
    else:
        sample = pod5_dir.resolve().name
        fastqs = {sample: nanopore.bam_to_fastq(calls_bam, basecall_dir / 'fastq' / f"{sample}.fastq",
                                                threads=threads, force=force_basecall)}

    mapping = load_barcode_map(barcode_map)
    if mapping:
        missing = sorted(set(fastqs) - set(mapping))
        if missing:
            logger.warning(f"Barcodes not in barcode map (kept with barcode names): {', '.join(missing)}")
        fastqs = {mapping.get(bc, bc): fq for bc, fq in fastqs.items()}

    stats_df = nanopore.basecall_summary(fastqs)
    empty = set(stats_df.loc[stats_df['reads'] == 0, 'sample'])
    fastqs = {s: fq for s, fq in fastqs.items() if s not in empty}

    if parallel_samples is None:
        parallel_samples = min(4, multiprocessing.cpu_count())
    parallel_samples = max(1, min(parallel_samples, len(fastqs) or 1))
    threads_per_sample = max(1, threads // parallel_samples)

    processor = AmpliconProcessor(
        reference_path=reference,
        bed_path=bed_path,
        min_base_quality=min_base_quality,
        min_mapping_quality=min_mapping_quality,
        min_read_count=min_read_count,
        max_indel_size=max_indel_size,
        platform='ont',
        minimap2_preset=minimap2_preset,
        ignore_homopolymer_indels=ignore_homopolymer_indels,
        homopolymer_min_length=homopolymer_min_length,
        full_length_tolerance=full_length_tolerance
    )

    def process_single_sample(sample_name: str, fastq: Path) -> Optional[Tuple[str, pd.DataFrame]]:
        try:
            logger.info(f"Processing sample: {sample_name}")
            return sample_name, processor.process_sample_ont(
                fastq=fastq,
                output_dir=output_dir / sample_name,
                sample_name=sample_name,
                threads=threads_per_sample
            )
        except Exception as e:
            logger.error(f"Error processing sample {sample_name}: {str(e)}")
            return None

    results = {}
    with ThreadPoolExecutor(max_workers=parallel_samples) as executor:
        futures = [executor.submit(process_single_sample, s, fq) for s, fq in fastqs.items()]
        with click.progressbar(length=len(futures), label='Processing samples') as bar:
            for future in as_completed(futures):
                result = future.result()
                if result is not None:
                    sample_name, df = result
                    results[sample_name] = df
                bar.update(1)

    return results, processor, stats_df

def summarize_results(results: Dict[str, pd.DataFrame], processor: Optional['AmpliconProcessor'] = None) -> pd.DataFrame:
    """Create a summary of results across all samples."""
    if not results:
        return pd.DataFrame()
        
    summaries = []
    # Store all data for overall statistics
    all_data = []
    
    for sample, df in results.items():
        if df.empty:
            continue
            
        try:
            total_reads = df['count'].sum()
            
            # Calculate per-reference statistics
            ref_stats = []
            for ref in df['reference'].unique():
                ref_df = df[df['reference'] == ref]
                ref_reads = ref_df['count'].sum()
                
                # Get reference sequence
                ref_seq = processor.reference[ref] if processor else None
                
                # Get theoretical maximum SNPs for this reference
                theoretical_max_snps = ref_df['theoretical_max_snps'].iloc[0]
                
                # Count unique single SNP haplotypes
                single_snp_haplotypes = ref_df[ref_df['snp_count'] == 1]
                unique_snp_positions = set()
                
                # Only count SNPs from haplotypes with exactly one SNP
                for _, row in single_snp_haplotypes.iterrows():
                    haplotype = row['haplotype']
                    for pos, (ref, var) in enumerate(zip(ref_seq, haplotype)):
                        if var.islower() and ref != var.upper():
                            unique_snp_positions.add(pos + 1)
                
                unique_single_snp_count = len(unique_snp_positions)
                
                # Add warning if we exceed theoretical maximum
                if unique_single_snp_count > theoretical_max_snps:
                    logger.warning(f"Sample {sample} reference {ref} has more unique single SNP positions "
                                 f"than theoretically possible: found={unique_single_snp_count}, "
                                 f"max={theoretical_max_snps}")
                
                ref_stats.append({
                    'reference': ref,
                    'reads': ref_reads,
                    'unique_haplotypes': len(ref_df),
                    'unique_single_mut_haplotypes': len(ref_df[ref_df['mutations'] == 1]),
                    'unique_single_snp_haplotypes': unique_single_snp_count,
                    'unique_single_indel_haplotypes': len(ref_df[ref_df['indel_count'] == 1]),
                    'max_frequency': ref_df['frequency'].max(),
                    'avg_mutations': (ref_df['mutations'] * ref_df['count']).sum() / ref_reads if ref_reads > 0 else 0,
                    'avg_snps': (ref_df['snp_count'] * ref_df['count']).sum() / ref_reads if ref_reads > 0 else 0,
                    'avg_indels': (ref_df['indel_count'] * ref_df['count']).sum() / ref_reads if ref_reads > 0 else 0,
                    'full_length_reads': full_length_read_count(ref_df),
                    'full_length_percent': (full_length_read_count(ref_df) / ref_reads * 100) if ref_reads > 0 else 0,
                    'theoretical_max_snps': theoretical_max_snps
                })
            
            # Create sample summary
            summary = {
                'sample': sample,
                'total_reads': total_reads,
                'unique_haplotypes': len(df),
                'unique_single_mut_haplotypes': len(df[df['mutations'] == 1]),
                'unique_single_snp_haplotypes': len(df[df['snp_count'] == 1]),
                'unique_single_indel_haplotypes': len(df[df['indel_count'] == 1]),
                'max_frequency': df['frequency'].max(),
                'avg_mutations': (df['mutations'] * df['count']).sum() / total_reads if total_reads > 0 else 0,
                'avg_snps': (df['snp_count'] * df['count']).sum() / total_reads if total_reads > 0 else 0,
                'avg_indels': (df['indel_count'] * df['count']).sum() / total_reads if total_reads > 0 else 0,
                'full_length_reads': full_length_read_count(df),
                'full_length_percent': (full_length_read_count(df) / total_reads * 100) if total_reads > 0 else 0,
                'num_references': len(df['reference'].unique())
            }
            
            # Add single mutation stats if available
            if 'single_mutations' in df.columns:
                single_mut_reads = df[df['mutations'] == 1]['count'].sum()
                summary.update({
                    'single_mutation_reads': single_mut_reads,
                    'single_mutation_percent': (single_mut_reads / total_reads * 100) if total_reads > 0 else 0
                })
            
            summaries.append(summary)
            all_data.append(df)
            
            # Add reference-specific statistics
            for ref_stat in ref_stats:
                ref_summary = {
                    'sample': f"{sample}_{ref_stat['reference']}",
                    'total_reads': ref_stat['reads'],
                    'unique_haplotypes': ref_stat['unique_haplotypes'],
                    'unique_single_mut_haplotypes': ref_stat['unique_single_mut_haplotypes'],
                    'unique_single_snp_haplotypes': ref_stat['unique_single_snp_haplotypes'],
                    'unique_single_indel_haplotypes': ref_stat['unique_single_indel_haplotypes'],
                    'max_frequency': ref_stat['max_frequency'],
                    'avg_mutations': ref_stat['avg_mutations'],
                    'avg_snps': ref_stat['avg_snps'],
                    'avg_indels': ref_stat['avg_indels'],
                    'full_length_reads': ref_stat['full_length_reads'],
                    'full_length_percent': ref_stat['full_length_percent'],
                    'num_references': 1,
                    'theoretical_max_snps': ref_stat['theoretical_max_snps']
                }
                if 'single_mutations' in df.columns:
                    ref_single_mut_reads = df[(df['reference'] == ref_stat['reference']) & (df['mutations'] == 1)]['count'].sum()
                    ref_summary.update({
                        'single_mutation_reads': ref_single_mut_reads,
                        'single_mutation_percent': (ref_single_mut_reads / ref_stat['reads'] * 100) if ref_stat['reads'] > 0 else 0
                    })
                summaries.append(ref_summary)
                
        except Exception as e:
            logger.error(f"Error summarizing results for {sample}: {str(e)}")
            continue
    
    # Create per-sample summary DataFrame
    summary_df = pd.DataFrame(summaries)
    
    # Calculate overall statistics if we have data
    if all_data:
        combined_df = pd.concat(all_data, ignore_index=True)
        total_reads = combined_df['count'].sum()
        
        # Calculate per-reference overall statistics
        ref_overall_stats = []
        for ref in combined_df['reference'].unique():
            ref_df = combined_df[combined_df['reference'] == ref]
            ref_reads = ref_df['count'].sum()
            
            # Get reference sequence
            ref_seq = processor.reference[ref] if processor else None
            
            # Get theoretical maximum SNPs for this reference
            theoretical_max_snps = ref_df['theoretical_max_snps'].iloc[0]
            
            # Count unique single SNP haplotypes
            single_snp_haplotypes = ref_df[ref_df['snp_count'] == 1]
            unique_snp_positions = set()
            
            # Only count SNPs from haplotypes with exactly one SNP
            for _, row in single_snp_haplotypes.iterrows():
                haplotype = row['haplotype']
                for pos, (ref, var) in enumerate(zip(ref_seq, haplotype)):
                    if var.islower() and ref != var.upper():
                        unique_snp_positions.add(pos + 1)
            
            unique_single_snp_count = len(unique_snp_positions)
            
            # Add warning if we exceed theoretical maximum
            if unique_single_snp_count > theoretical_max_snps:
                logger.warning(f"Overall results for reference {ref} have more unique single SNP positions "
                             f"than theoretically possible: found={unique_single_snp_count}, "
                             f"max={theoretical_max_snps}")
            
            ref_overall_stats.append({
                'reference': ref,
                'reads': ref_reads,
                'unique_haplotypes': len(ref_df),
                'unique_single_mut_haplotypes': len(ref_df[ref_df['mutations'] == 1]),
                'unique_single_snp_haplotypes': unique_single_snp_count,
                'unique_single_indel_haplotypes': len(ref_df[ref_df['indel_count'] == 1]),
                'max_frequency': ref_df['frequency'].max(),
                'avg_mutations': (ref_df['mutations'] * ref_df['count']).sum() / ref_reads if ref_reads > 0 else 0,
                'avg_snps': (ref_df['snp_count'] * ref_df['count']).sum() / ref_reads if ref_reads > 0 else 0,
                'avg_indels': (ref_df['indel_count'] * ref_df['count']).sum() / ref_reads if ref_reads > 0 else 0,
                'full_length_reads': full_length_read_count(ref_df),
                'full_length_percent': (full_length_read_count(ref_df) / ref_reads * 100) if ref_reads > 0 else 0,
                'theoretical_max_snps': theoretical_max_snps
            })
        
        # Create overall summary
        overall_summary = {
            'sample': 'OVERALL',
            'total_reads': total_reads,
            'unique_haplotypes': len(combined_df),
            'unique_single_mut_haplotypes': len(combined_df[combined_df['mutations'] == 1]),
            'unique_single_snp_haplotypes': len(combined_df[combined_df['snp_count'] == 1]),
            'unique_single_indel_haplotypes': len(combined_df[combined_df['indel_count'] == 1]),
            'max_frequency': combined_df['frequency'].max(),
            'avg_mutations': (combined_df['mutations'] * combined_df['count']).sum() / total_reads if total_reads > 0 else 0,
            'avg_snps': (combined_df['snp_count'] * combined_df['count']).sum() / total_reads if total_reads > 0 else 0,
            'avg_indels': (combined_df['indel_count'] * combined_df['count']).sum() / total_reads if total_reads > 0 else 0,
            'full_length_reads': full_length_read_count(combined_df),
            'full_length_percent': (full_length_read_count(combined_df) / total_reads * 100) if total_reads > 0 else 0,
            'num_references': len(combined_df['reference'].unique())
        }
        
        # Add single mutation stats if available
        if 'single_mutations' in combined_df.columns:
            single_mut_reads = combined_df[combined_df['mutations'] == 1]['count'].sum()
            overall_summary.update({
                'single_mutation_reads': single_mut_reads,
                'single_mutation_percent': (single_mut_reads / total_reads * 100) if total_reads > 0 else 0
            })
        
        # Add overall summary as the first row
        summary_df = pd.concat([pd.DataFrame([overall_summary]), summary_df], ignore_index=True)
        
        # Add reference-specific overall statistics
        for ref_stat in ref_overall_stats:
            ref_overall_summary = {
                'sample': f"OVERALL_{ref_stat['reference']}",
                'total_reads': ref_stat['reads'],
                'unique_haplotypes': ref_stat['unique_haplotypes'],
                'unique_single_mut_haplotypes': ref_stat['unique_single_mut_haplotypes'],
                'unique_single_snp_haplotypes': ref_stat['unique_single_snp_haplotypes'],
                'unique_single_indel_haplotypes': ref_stat['unique_single_indel_haplotypes'],
                'max_frequency': ref_stat['max_frequency'],
                'avg_mutations': ref_stat['avg_mutations'],
                'avg_snps': ref_stat['avg_snps'],
                'avg_indels': ref_stat['avg_indels'],
                'full_length_reads': ref_stat['full_length_reads'],
                'full_length_percent': ref_stat['full_length_percent'],
                'num_references': 1,
                'theoretical_max_snps': ref_stat['theoretical_max_snps']
            }
            if 'single_mutations' in combined_df.columns:
                ref_single_mut_reads = combined_df[(combined_df['reference'] == ref_stat['reference']) & (combined_df['mutations'] == 1)]['count'].sum()
                ref_overall_summary.update({
                    'single_mutation_reads': ref_single_mut_reads,
                    'single_mutation_percent': (ref_single_mut_reads / ref_stat['reads'] * 100) if ref_stat['reads'] > 0 else 0
                })
            summary_df = pd.concat([pd.DataFrame([ref_overall_summary]), summary_df], ignore_index=True)
    
    # Format numeric columns
    if not summary_df.empty:
        # Round percentages to 2 decimal places
        for col in ['max_frequency', 'full_length_percent', 'single_mutation_percent']:
            if col in summary_df.columns:
                summary_df[col] = summary_df[col].round(2)
        
        # Round average mutations to 2 decimal places
        if 'avg_mutations' in summary_df.columns:
            summary_df['avg_mutations'] = summary_df['avg_mutations'].round(2)
    
    return summary_df

def validate_input(
    fastq_dir: Union[str, Path],
    reference: Union[str, Path]
) -> List[str]:
    """
    Validate input files and return any warnings.

    Args:
        fastq_dir: Directory containing FASTQ files
        reference: Path to reference FASTA file

    Returns:
        List of warning messages, empty if all valid
    """
    warnings = []
    
    # Check reference file
    ref_path = Path(reference)
    if not ref_path.exists():
        warnings.append(f"Reference file not found: {ref_path}")
    elif ref_path.stat().st_size == 0:
        warnings.append(f"Reference file is empty: {ref_path}")
    else:
        # Validate FASTA format
        try:
            with open(ref_path) as handle:
                records = list(SeqIO.parse(handle, "fasta"))
                if not records:
                    warnings.append(f"No valid FASTA sequences found in: {ref_path}")
        except Exception as e:
            warnings.append(f"Error reading reference file: {str(e)}")

    # Check BWA index files
    for ext in ['.amb', '.ann', '.bwt', '.pac', '.sa']:
        if not (ref_path.parent / (ref_path.name + ext)).exists():
            warnings.append(f"BWA index file missing: {ref_path}{ext}")

    # Check required executables
    for cmd in ['bwa', 'samtools']:
        if not shutil.which(cmd):
            warnings.append(f"Required program not found: {cmd}")

    # Check FASTQ directory
    fastq_dir = Path(fastq_dir)
    if not fastq_dir.is_dir():
        warnings.append(f"FASTQ directory not found: {fastq_dir}")
    else:
        # Look for both .fastq and .fastq.gz files
        r1_files = list(fastq_dir.glob('*_R1_001.fastq*'))
        if not r1_files:
            warnings.append(f"No R1 FASTQ files found in: {fastq_dir}")
        
        # Check for matching R2 files
        for r1 in r1_files:
            r2 = r1.parent / r1.name.replace('_R1_', '_R2_')
            if not r2.exists():
                warnings.append(f"No matching R2 file for: {r1.name}")
            
            # Check file sizes
            try:
                if r1.stat().st_size == 0:
                    warnings.append(f"Empty R1 file: {r1.name}")
                if r2.exists() and r2.stat().st_size == 0:
                    warnings.append(f"Empty R2 file: {r2.name}")
            except Exception as e:
                warnings.append(f"Error checking file sizes: {str(e)}")

    return warnings

def validate_nanopore_input(
    pod5_dir: Union[str, Path],
    reference: Union[str, Path],
    dorado_bin: Union[str, Path, None] = None,
    require_dorado: bool = True
) -> List[str]:
    """
    Validate Nanopore inputs and return any warnings.

    Args:
        pod5_dir: Directory containing POD5 files (searched recursively)
        reference: Path to reference FASTA file
        dorado_bin: Optional path to dorado executable
        require_dorado: Whether dorado must be available (False when resuming from existing basecalls)

    Returns:
        List of warning messages, empty if all valid
    """
    warnings = []

    ref_path = Path(reference)
    if not ref_path.exists():
        warnings.append(f"Reference file not found: {ref_path}")
    elif ref_path.stat().st_size == 0:
        warnings.append(f"Reference file is empty: {ref_path}")
    else:
        try:
            with open(ref_path) as handle:
                if not list(SeqIO.parse(handle, "fasta")):
                    warnings.append(f"No valid FASTA sequences found in: {ref_path}")
        except Exception as e:
            warnings.append(f"Error reading reference file: {str(e)}")

    for cmd in ['minimap2', 'samtools']:
        if not shutil.which(cmd):
            warnings.append(f"Required program not found: {cmd}")

    if require_dorado:
        try:
            nanopore.check_dorado(dorado_bin)
        except RuntimeError as e:
            warnings.append(str(e))

    pod5_dir = Path(pod5_dir)
    if not pod5_dir.is_dir():
        warnings.append(f"POD5 directory not found: {pod5_dir}")
    else:
        raw = nanopore.find_raw_files(pod5_dir)
        if raw['fast5'] and not raw['pod5']:
            warnings.append(f"Only FAST5 files found in {pod5_dir}; convert with `pod5 convert fast5` "
                            f"for best dorado performance")
        if not raw['pod5'] and not raw['fast5'] and require_dorado:
            warnings.append(f"No POD5 files found in: {pod5_dir}")

    return warnings

def load_results(results_dir: Union[str, Path]) -> Dict[str, pd.DataFrame]:
    """Load previously generated results from CSV files."""
    results_dir = Path(results_dir)
    results = {}
    
    try:
        # Look for results in sample subdirectories
        for csv_file in results_dir.rglob('*_haplotypes.csv'):
            sample_name = csv_file.name.replace('_haplotypes.csv', '')
            try:
                df = pd.read_csv(csv_file)
                
                # Ensure required columns are present
                required_cols = ['reference', 'haplotype', 'count', 'frequency', 'mutations']
                missing_cols = [col for col in required_cols if col not in df.columns]
                
                if missing_cols:
                    logger.warning(f"Results file {csv_file} missing columns: {missing_cols}")
                    continue
                    
                results[sample_name] = df
                logger.debug(f"Loaded results for sample: {sample_name}")
                
            except Exception as e:
                logger.error(f"Error loading results for {sample_name}: {str(e)}")
                continue
    except Exception as e:
        logger.error(f"Error scanning results directory: {str(e)}")
    
    return results