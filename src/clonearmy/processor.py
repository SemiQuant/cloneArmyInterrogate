from dataclasses import dataclass
from pathlib import Path
import subprocess
import tempfile
from typing import List, Dict, Generator, Tuple, Set, Union
import logging
from collections import Counter, defaultdict
import re
import shutil
import threading
import click

import pysam
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from rich.progress import track
import pyranges as pr

logger = logging.getLogger(__name__)

@dataclass
class AmpliconRead:
    """Represents a processed amplicon read pair."""
    sequence: str
    mutations: int
    quality: float
    indels: List[Dict]  # Add indels field
    insertions: str = ''

class AmpliconProcessor:
    """Process amplicon sequencing data."""
    
    def __init__(self, 
                 reference_path: Union[str, Path],
                 bed_path: Union[str, Path, None] = None,
                 min_base_quality: int = 20,
                 min_mapping_quality: int = 30,
                 min_read_count: int = 10,
                 max_file_size: int = 10_000_000_000,
                 max_indel_size: int = 50,  # Add max_indel_size parameter
                 platform: str = "illumina",
                 minimap2_preset: str = "lr:hq",
                 ignore_homopolymer_indels: bool = True,
                 homopolymer_min_length: int = 3,
                 max_homopolymer_indel: int = 1):
        """
        Initialize the processor.
        
        Args:
            reference_path: Path to reference FASTA file
            bed_path: Optional path to BED file for indel comparison
            min_base_quality: Minimum base quality score
            min_mapping_quality: Minimum mapping quality score
            min_read_count: Minimum number of reads to consider a haplotype
            max_file_size: Maximum file size in bytes
            max_indel_size: Maximum size of indels to consider as "small" indels
            platform: Sequencing platform, "illumina" (paired-end) or "ont" (Nanopore single reads)
            minimap2_preset: minimap2 preset used for ONT alignment
            ignore_homopolymer_indels: (ONT) ignore short indels in reference homopolymer runs
            homopolymer_min_length: (ONT) minimum run length to treat as a homopolymer
            max_homopolymer_indel: (ONT) maximum indel size ignored within homopolymers
        """
        platform = platform.lower()
        if platform not in ('illumina', 'ont'):
            raise ValueError(f"Unsupported platform: {platform} (expected 'illumina' or 'ont')")
        self.platform = platform
        self.minimap2_preset = minimap2_preset
        self.ignore_homopolymer_indels = ignore_homopolymer_indels
        self.homopolymer_min_length = homopolymer_min_length
        self.max_homopolymer_indel = max_homopolymer_indel
        self.reference_path = Path(reference_path)
        self.bed_path = Path(bed_path) if bed_path else None
        self.min_base_quality = min_base_quality
        self.min_mapping_quality = min_mapping_quality
        self.min_read_count = min_read_count
        self.max_file_size = max_file_size
        self.max_indel_size = max_indel_size
        
        # Load reference sequences
        self.reference = {}
        try:
            for record in SeqIO.parse(self.reference_path, "fasta"):
                self.reference[record.id] = str(record.seq)
        except Exception as e:
            logger.error(f"Error loading reference sequences: {str(e)}")
            raise
            
        # Load BED regions if provided
        self.bed_regions = None
        if self.bed_path:
            try:
                self.bed_regions = pr.read_bed(str(self.bed_path))
            except Exception as e:
                logger.error(f"Error loading BED file: {str(e)}")
                raise
        
        # Initialize aligner for indel detection
        self.aligner = PairwiseAligner()
        self.aligner.mode = 'global'
        self.aligner.match_score = 2
        self.aligner.mismatch_score = -1
        self.aligner.open_gap_score = -2
        self.aligner.extend_gap_score = -0.5
        self._alignment_cache: Dict[Tuple[str, str], str] = {}
        # The pairwise aligner is reconfigured per call, so it cannot be shared across threads
        self._aligner_lock = threading.Lock()
        
        # Check for required executables and index files
        self._check_dependencies()
        if self.platform == 'illumina':
            self._check_and_create_bwa_index()

    def _load_reference(self) -> Dict[str, str]:
        """Load reference sequences."""
        try:
            reference_dict = {}
            with open(self.reference_path) as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    reference_dict[record.id] = str(record.seq)
            if not reference_dict:
                raise ValueError(f"No sequences found in reference file: {self.reference_path}")
            return reference_dict
        except Exception as e:
            logger.error(f"Error loading reference sequence: {str(e)}")
            raise

    def _check_dependencies(self):
        """Check if required external programs are available."""
        required = ['minimap2', 'samtools'] if self.platform == 'ont' else ['bwa', 'samtools', 'seqtk']
        for cmd in required:
            if not shutil.which(cmd):
                raise RuntimeError(f"{cmd} not found in PATH. Please install {cmd}.")

    def _check_and_create_bwa_index(self):
        """Check if BWA index files exist, create them if they don't."""
        index_extensions = ['.amb', '.ann', '.bwt', '.pac', '.sa']
        missing_indices = [ext for ext in index_extensions 
                         if not (self.reference_path.parent / f"{self.reference_path.name}{ext}").exists()]
        
        if missing_indices:
            logger.info(f"Creating BWA index for {self.reference_path}")
            with click.progressbar(length=1, label='Indexing reference', show_eta=True) as bar:
                try:
                    # First check if reference file exists
                    if not self.reference_path.exists():
                        raise RuntimeError(f"Reference file not found: {self.reference_path}")
                    
                    # Check if reference file is empty
                    if self.reference_path.stat().st_size == 0:
                        raise RuntimeError(f"Reference file is empty: {self.reference_path}")
                    
                    # Run BWA index with detailed error capture
                    result = subprocess.run(
                        ['bwa', 'index', str(self.reference_path)],
                        check=True,
                        stderr=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        text=True
                    )
                    
                    # Verify index creation
                    still_missing = [ext for ext in index_extensions 
                                   if not (self.reference_path.parent / f"{self.reference_path.name}{ext}").exists()]
                    
                    if still_missing:
                        raise RuntimeError(f"BWA indexing failed to create files: {', '.join(still_missing)}")
                    
                    bar.update(1)
                    logger.info("BWA index created successfully")
                    
                except subprocess.CalledProcessError as e:
                    error_msg = f"BWA indexing failed: {e.stderr}"
                    logger.error(error_msg)
                    raise RuntimeError(error_msg)
                except Exception as e:
                    error_msg = f"Error during BWA indexing: {str(e)}"
                    logger.error(error_msg)
                    raise RuntimeError(error_msg)

    def _align_sequence_to_reference(self, sequence: str, ref_seq: str) -> str:
        """Align a sequence to reference to detect indels and mutations."""
        try:
            # Optimize aligner settings for better alignment
            self.aligner.mode = 'global'
            self.aligner.match_score = 2
            self.aligner.mismatch_score = -1  # Reduced penalty
            self.aligner.open_gap_score = -2  # Reduced penalty
            self.aligner.extend_gap_score = -0.5  # Reduced penalty
            self.aligner.target_internal_open_gap_score = -2
            self.aligner.target_internal_extend_gap_score = -0.5
            self.aligner.target_left_open_gap_score = -2
            self.aligner.target_left_extend_gap_score = -0.5
            self.aligner.target_right_open_gap_score = -2
            self.aligner.target_right_extend_gap_score = -0.5
            self.aligner.query_internal_open_gap_score = -2
            self.aligner.query_internal_extend_gap_score = -0.5
            self.aligner.query_left_open_gap_score = -2
            self.aligner.query_left_extend_gap_score = -0.5
            self.aligner.query_right_open_gap_score = -2
            self.aligner.query_right_extend_gap_score = -0.5
            
            # First try score_only mode to check if alignment is possible
            try:
                score = self.aligner.score(ref_seq, sequence)
                if score < -len(sequence):  # More lenient threshold
                    logger.warning("Poor alignment score, trying local alignment")
                    # If global alignment fails, try local alignment
                    self.aligner.mode = 'local'
                    self.aligner.mismatch_score = -1
                    self.aligner.open_gap_score = -2
                    
                    score = self.aligner.score(ref_seq, sequence)
                    if score < -len(sequence):
                        logger.warning("Poor alignment score even with local alignment, returning original sequence")
                        return sequence
                    
                    alignment = next(self.aligner.align(ref_seq, sequence))
                else:
                    alignment = next(self.aligner.align(ref_seq, sequence))
                
            except (StopIteration, ValueError) as e:
                logger.warning(f"Alignment failed: {str(e)}, returning original sequence")
                return sequence
            
            # Process the alignment to mark mutations and gaps
            result = []
            target_seq = str(alignment.target)
            query_seq = str(alignment.query)
            
            for t, q in zip(target_seq, query_seq):
                if t == '-':  # Insertion relative to reference
                    result.append(q.lower())
                elif q == '-':  # Deletion relative to reference
                    result.append('-')
                elif t.upper() != q.upper():  # Mismatch
                    result.append(q.lower())
                else:  # Match
                    result.append(q.upper())
            
            return ''.join(result)
            
        except Exception as e:
            logger.error(f"Error in sequence alignment: {str(e)}")
            return sequence

    def _reconstruct_sequence(self,
                            read1: pysam.AlignedSegment,
                            read2: pysam.AlignedSegment,
                            ref_seq: str) -> str:
        """Reconstruct the amplicon sequence from paired reads."""
        return self._reconstruct_from_reads((read1, read2), ref_seq)

    def _homopolymer_run_length(self, ref_seq: str, pos: int) -> int:
        """Length of the reference homopolymer run containing pos (0 if out of range)."""
        if pos < 0 or pos >= len(ref_seq):
            return 0
        base = ref_seq[pos].upper()
        start = pos
        while start > 0 and ref_seq[start - 1].upper() == base:
            start -= 1
        end = pos
        while end < len(ref_seq) - 1 and ref_seq[end + 1].upper() == base:
            end += 1
        return end - start + 1

    def _in_homopolymer(self, ref_seq: str, pos: int, min_len: int, base: str = None) -> bool:
        """True if pos lies in a reference homopolymer of at least min_len (optionally of base)."""
        if pos < 0 or pos >= len(ref_seq):
            return False
        if base is not None and ref_seq[pos].upper() != base.upper():
            return False
        return self._homopolymer_run_length(ref_seq, pos) >= min_len

    def _is_homopolymer_indel(self, ref_seq: str, ref_pos: int, op: int, bases: str) -> bool:
        """Whether an ONT indel looks like a homopolymer length error and should be ignored."""
        if not self.ignore_homopolymer_indels or len(bases) > self.max_homopolymer_indel:
            return False
        if len(set(bases.upper())) != 1:
            return False
        base = bases[0]
        min_len = self.homopolymer_min_length
        if op == 2:  # Deletion of ref_seq[ref_pos:ref_pos+len]
            return self._in_homopolymer(ref_seq, ref_pos, min_len, base)
        # Insertion between ref_pos-1 and ref_pos extends an adjacent run of the same base
        return (self._in_homopolymer(ref_seq, ref_pos - 1, min_len, base) or
                self._in_homopolymer(ref_seq, ref_pos, min_len, base))

    @staticmethod
    def format_insertions(insertions: Dict[int, str]) -> str:
        """Format {junction: bases} as 'P_P+1insBASES;...' (P = 1-based preceding reference base)."""
        return ';'.join(f"{pos}_{pos + 1}ins{bases.upper()}" for pos, bases in sorted(insertions.items()))

    @staticmethod
    def parse_insertions(insertions) -> List[Tuple[int, str]]:
        """Parse an insertions string into [(1-based preceding reference position, BASES), ...]."""
        if not isinstance(insertions, str) or not insertions:
            return []
        parsed = []
        for item in insertions.split(';'):
            match = re.fullmatch(r'(\d+)_\d+ins([A-Za-z]+)', item.strip())
            if match:
                parsed.append((int(match.group(1)), match.group(2).upper()))
        return parsed

    def _reconstruct_from_reads(self, reads, ref_seq: str) -> str:
        """Reconstruct the amplicon sequence from one or more reads aligned to ref_seq."""
        return self._reconstruct_with_insertions(reads, ref_seq)[0]

    def _reconstruct_with_insertions(self, reads, ref_seq: str) -> Tuple[str, str]:
        """Reconstruct the amplicon from reads aligned to ref_seq.

        Returns the reference-length haplotype (substitutions lowercase, deletions '-')
        and the insertions string. Insertions are kept out of the haplotype so positions
        stay in reference coordinates. Where reads overlap, the later read decides at any
        insertion junction it spans, as it does for substitutions.
        """
        sequence = list(ref_seq.upper())
        is_ont = self.platform == 'ont'
        # Keyed by junction j: inserted between 0-based reference positions j-1 and j
        insertions: Dict[int, str] = {}
        
        for read in reads:
            read_seq = read.query_sequence
            quals = read.query_qualities
            ref_pos = read.reference_start
            read_insertions: Dict[int, str] = {}
            uninformative = set()
            
            query_pos = 0
            for op, length in read.cigartuples:
                if op in (0, 7, 8):  # Match or mismatch (M, =, X)
                    for i in range(length):
                        if (quals[query_pos + i] >= self.min_base_quality and
                            ref_pos + i < len(sequence)):
                            base = read_seq[query_pos + i].upper()
                            if base != ref_seq[ref_pos + i].upper():
                                sequence[ref_pos + i] = base.lower()
                            else:
                                sequence[ref_pos + i] = base.upper()
                    query_pos += length
                    ref_pos += length
                elif op == 1:  # Insertion
                    inserted_bases = read_seq[query_pos:query_pos + length].upper()
                    passes_quality = all(q >= self.min_base_quality
                                         for q in quals[query_pos:query_pos + length])
                    if not passes_quality:
                        uninformative.add(ref_pos)
                    elif (0 <= ref_pos <= len(sequence) and
                          not (is_ont and self._is_homopolymer_indel(ref_seq, ref_pos, 1, inserted_bases))):
                        read_insertions[ref_pos] = inserted_bases
                    query_pos += length
                elif op == 2:  # Deletion
                    deleted = ref_seq[ref_pos:ref_pos + length]
                    if not (is_ont and self._is_homopolymer_indel(ref_seq, ref_pos, 2, deleted)):
                        # Mark deletion with '-'
                        for i in range(length):
                            if ref_pos + i < len(sequence):
                                sequence[ref_pos + i] = '-'
                    ref_pos += length
                elif op == 3:  # Skipped reference region
                    ref_pos += length
                elif op == 4:  # Soft clip
                    query_pos += length

            # Junctions strictly inside this read's alignment are covered by it
            for junction in [j for j in insertions
                             if read.reference_start < j < read.reference_end
                             and j not in read_insertions and j not in uninformative]:
                del insertions[junction]
            insertions.update(read_insertions)
        
        reconstructed = ''.join(sequence)
        cache_key = (reconstructed, ref_seq)
        with self._aligner_lock:
            cached = self._alignment_cache.get(cache_key)
            if cached is None:
                cached = self._align_sequence_to_reference(reconstructed, ref_seq)
                self._alignment_cache[cache_key] = cached
        return cached, self.format_insertions(insertions)

    def _is_full_length(self, sequence: str, ref_seq: str) -> bool:
        """Check if a sequence covers the full reference length.
        
        A sequence is considered full-length if:
        1. It has approximately the same length as the reference when accounting for indels
        2. It doesn't start or end with a deletion
        3. It doesn't start or end with an N
        """
        # Remove all deletion markers
        seq_no_indels = sequence.replace('-', '')
        ref_no_indels = ref_seq.replace('-', '')
        
        # Check if sequence starts/ends with N
        if sequence.startswith('N') or sequence.endswith('N'):
            return False
            
        # Check if sequence starts/ends with deletion
        if sequence.startswith('-') or sequence.endswith('-'):
            return False
            
        # For sequences with indels, we need to be more lenient
        # Calculate length difference accounting for indels
        len_diff = abs(len(seq_no_indels) - len(ref_no_indels))
        max_allowed_diff = min(self.max_indel_size, len(ref_no_indels) * 0.1)  # 10% of reference length or max_indel_size
        
        return len_diff <= max_allowed_diff

    def _get_read_pairs(self, 
                       bam: pysam.AlignmentFile,
                       ref_name: str) -> Generator[Tuple[pysam.AlignedSegment, pysam.AlignedSegment], None, None]:
        """Generate properly paired reads."""
        reads = {}
        for read in bam.fetch(ref_name):
            if (not read.is_proper_pair or 
                read.is_secondary or 
                read.is_supplementary or 
                read.mapping_quality < self.min_mapping_quality):
                continue
                
            qname = read.query_name
            if qname in reads:
                pair = reads.pop(qname)
                yield (read, pair) if read.is_read1 else (pair, read)
            else:
                reads[qname] = read

    def _align_reads(self, 
                    fastq_r1: Path,
                    fastq_r2: Path, 
                    temp_dir: Path,
                    output_dir: Path,
                    threads: int) -> Path:
        """Align reads using BWA-MEM and convert to sorted BAM."""
        sample_name = fastq_r1.stem.replace("_R1_001.fastq", "")
        temp_sam = temp_dir / f"{sample_name}.sam"
        temp_bam = temp_dir / f"{sample_name}.temp.bam"
        final_bam = output_dir / f"{sample_name}.bam"
        
        try:
            with click.progressbar(length=4, label='Aligning and processing reads') as bar:
                # Run BWA-MEM
                bwa_cmd = [
                    'bwa', 'mem',
                    '-t', str(threads),
                    str(self.reference_path),
                    str(fastq_r1),
                    str(fastq_r2)
                ]
                
                with open(temp_sam, 'w') as sam_out:
                    click.echo("\nRunning BWA alignment...")
                    logger.debug(f"Running BWA: {' '.join(bwa_cmd)}")
                    subprocess.run(
                        bwa_cmd,
                        stdout=sam_out,
                        stderr=subprocess.PIPE,
                        check=True
                    )
                bar.update(1)
                
                # Convert SAM to BAM
                click.echo("Converting SAM to BAM...")
                subprocess.run(
                    ['samtools', 'view', '-b', '-@', str(threads), '-o', str(temp_bam), str(temp_sam)],
                    check=True,
                    stderr=subprocess.PIPE
                )
                bar.update(1)
                
                # Sort BAM
                click.echo("Sorting BAM file...")
                subprocess.run(
                    [
                        'samtools', 'sort',
                        '-@', str(threads),
                        '-m', '1G',
                        '-T', str(temp_dir / f"{sample_name}.sort"),
                        '-o', str(final_bam),
                        str(temp_bam)
                    ],
                    check=True,
                    stderr=subprocess.PIPE
                )
                bar.update(1)
                
                # Index BAM
                click.echo("Indexing BAM file...")
                subprocess.run(
                    ['samtools', 'index', str(final_bam)],
                    check=True,
                    stderr=subprocess.PIPE
                )
                bar.update(1)
            
            return final_bam
            
        except subprocess.CalledProcessError as e:
            error_msg = e.stderr.decode() if e.stderr else str(e)
            raise RuntimeError(f"Alignment failed: {error_msg}")
        except Exception as e:
            raise RuntimeError(f"Alignment failed: {str(e)}")
        finally:
            for temp_file in [temp_sam, temp_bam]:
                try:
                    if temp_file.exists():
                        temp_file.unlink()
                except:
                    pass

    def _resolve_minimap2_preset(self) -> str:
        """Return the configured minimap2 preset, falling back to map-ont if lr:hq is unsupported."""
        preset = self.minimap2_preset
        if preset != 'lr:hq':
            return preset
        try:
            out = subprocess.run(['minimap2', '--version'], capture_output=True, text=True).stdout.strip()
            major, minor = (int(x) for x in out.split('-')[0].split('.')[:2])
            if (major, minor) < (2, 27):
                logger.warning(f"minimap2 {out} does not support lr:hq, using map-ont")
                return 'map-ont'
        except Exception as e:
            logger.debug(f"Could not determine minimap2 version: {e}")
        return preset

    def _align_reads_ont(self,
                         fastq: Path,
                         temp_dir: Path,
                         output_dir: Path,
                         threads: int,
                         sample_name: str) -> Path:
        """Align Nanopore reads with minimap2 and convert to sorted, indexed BAM."""
        temp_sam = temp_dir / f"{sample_name}.sam"
        final_bam = output_dir / f"{sample_name}.bam"
        preset = self._resolve_minimap2_preset()

        try:
            with click.progressbar(length=3, label='Aligning and processing reads') as bar:
                minimap2_cmd = [
                    'minimap2', '-ax', preset,
                    '-t', str(threads),
                    '--secondary=no',
                    str(self.reference_path),
                    str(fastq)
                ]
                with open(temp_sam, 'w') as sam_out:
                    click.echo(f"\nRunning minimap2 alignment ({preset})...")
                    logger.debug(f"Running minimap2: {' '.join(minimap2_cmd)}")
                    subprocess.run(minimap2_cmd, stdout=sam_out, stderr=subprocess.PIPE, check=True)
                bar.update(1)

                click.echo("Sorting BAM file...")
                subprocess.run(
                    [
                        'samtools', 'sort',
                        '-@', str(threads),
                        '-m', '1G',
                        '-T', str(temp_dir / f"{sample_name}.sort"),
                        '-o', str(final_bam),
                        str(temp_sam)
                    ],
                    check=True,
                    stderr=subprocess.PIPE
                )
                bar.update(1)

                click.echo("Indexing BAM file...")
                subprocess.run(['samtools', 'index', str(final_bam)], check=True, stderr=subprocess.PIPE)
                bar.update(1)

            return final_bam

        except subprocess.CalledProcessError as e:
            error_msg = e.stderr.decode() if e.stderr else str(e)
            raise RuntimeError(f"Alignment failed: {error_msg}")
        except Exception as e:
            raise RuntimeError(f"Alignment failed: {str(e)}")
        finally:
            try:
                if temp_sam.exists():
                    temp_sam.unlink()
            except OSError:
                pass

    def _is_usable_single_read(self, read: pysam.AlignedSegment) -> bool:
        """Filter for single-end (ONT) reads: mapped, primary, with base qualities and sufficient MAPQ."""
        return not (read.is_unmapped or
                    read.is_secondary or
                    read.is_supplementary or
                    read.query_qualities is None or
                    read.mapping_quality < self.min_mapping_quality)

    def _process_single_alignments(self,
                                   bam: pysam.AlignmentFile,
                                   ref_name: str,
                                   ref_seq: str) -> Generator[AmpliconRead, None, None]:
        """Process single-end (Nanopore) reads for a reference sequence."""
        mapped = next((s.mapped for s in bam.get_index_statistics() if s.contig == ref_name), 0)
        if mapped == 0:
            click.echo(f"No mapped reads found for {ref_name}")
            return

        with click.progressbar(length=mapped, label=f'Processing reads for {ref_name}') as bar:
            for read in bam.fetch(ref_name):
                bar.update(1)
                if not self._is_usable_single_read(read):
                    continue
                sequence, insertions = self._reconstruct_with_insertions((read,), ref_seq)
                mutations = (sum(1 for base in sequence if base.islower() or base == '-') +
                             len(self.parse_insertions(insertions)))
                yield AmpliconRead(sequence=sequence, mutations=mutations,
                                   quality=float(read.mapping_quality), indels=[],
                                   insertions=insertions)

    def _process_alignments(self, 
                          bam_path: Path,
                          ref_name: str) -> Generator[AmpliconRead, None, None]:
        """Process aligned reads for a reference sequence."""
        try:
            bam = pysam.AlignmentFile(bam_path, "rb")
            ref_seq = self.reference[ref_name]

            if self.platform == 'ont':
                yield from self._process_single_alignments(bam, ref_name, ref_seq)
                return
            
            # First pass: count total proper pairs
            read_pairs = {}
            total_pairs = 0
            
            for read in bam.fetch(ref_name):
                if (not read.is_proper_pair or 
                    read.is_secondary or 
                    read.is_supplementary or 
                    read.mapping_quality < self.min_mapping_quality):
                    continue
                
                qname = read.query_name
                if qname in read_pairs:
                    total_pairs += 1
                    read_pairs.pop(qname)
                else:
                    read_pairs[qname] = read
            
            # Reset file pointer and clear read_pairs
            bam.reset()
            read_pairs.clear()
            
            if total_pairs == 0:
                click.echo(f"No valid read pairs found for {ref_name}")
                return
            
            # Second pass: process read pairs
            with click.progressbar(length=total_pairs, 
                                 label=f'Processing reads for {ref_name}') as bar:
                
                for read in bam.fetch(ref_name):
                    if (not read.is_proper_pair or 
                        read.is_secondary or 
                        read.is_supplementary or 
                        read.mapping_quality < self.min_mapping_quality):
                        continue
                    
                    qname = read.query_name
                    if qname in read_pairs:
                        # Found a pair
                        pair = read_pairs.pop(qname)
                        read1, read2 = (read, pair) if read.is_read1 else (pair, read)
                        
                        sequence, insertions = self._reconstruct_with_insertions((read1, read2), ref_seq)
                        mutations = (sum(1 for base in sequence if base.islower() or base == '-') +
                                     len(self.parse_insertions(insertions)))
                        quality = (read1.mapping_quality + read2.mapping_quality) / 2
                        
                        bar.update(1)
                        
                        yield AmpliconRead(sequence=sequence, mutations=mutations, quality=quality,
                                           indels=[], insertions=insertions)
                    else:
                        read_pairs[qname] = read
            
        except Exception as e:
            logger.error(f"Error processing alignments for {ref_name}: {str(e)}")
            raise
        finally:
            if 'bam' in locals():
                bam.close()

    def _is_valid_snp(self, ref_base: str, alt_base: str) -> bool:
        """Validate if a mutation is a valid SNP.
        
        Args:
            ref_base: Reference base
            alt_base: Alternative base
            
        Returns:
            bool: True if the mutation is a valid SNP
        """
        valid_bases = {'A', 'C', 'G', 'T'}
        return (ref_base in valid_bases and 
                alt_base in valid_bases and 
                ref_base != alt_base)

    def _get_mutation_positions(self, haplotype: str, ref_seq: str, ref_name: str,
                                insertions: str = '') -> List[Dict[str, str]]:
        """Analyze mutation positions in a haplotype compared to reference.
        
        Args:
            haplotype: The haplotype sequence with mutations in lowercase
            ref_seq: The reference sequence
            ref_name: Name of the reference sequence
            insertions: Insertions string (see format_insertions); position is the preceding base
            
        Returns:
            List of dictionaries containing mutation information
        """
        mutations = []
        pos = 0  # 0-based position in reference
        hap_pos = 0  # 0-based position in haplotype
        
        # First validate sequence length (excluding indels)
        hap_no_indels = haplotype.replace('-', '')
        ref_no_indels = ref_seq.replace('-', '')
        
        # More lenient length validation
        len_diff = abs(len(hap_no_indels) - len(ref_no_indels))
        max_allowed_diff = min(self.max_indel_size, len(ref_no_indels) * 0.1)
        
        if len_diff > max_allowed_diff:
            logger.warning(f"Haplotype length mismatch for {ref_name}: "
                         f"reference={len(ref_no_indels)}, haplotype={len(hap_no_indels)}")
            return mutations
        
        while hap_pos < len(haplotype) and pos < len(ref_seq):
            if haplotype[hap_pos] == '-':  # Deletion
                # Calculate deletion size
                del_size = 1
                next_pos = hap_pos + 1
                while next_pos < len(haplotype) and haplotype[next_pos] == '-':
                    del_size += 1
                    next_pos += 1
                
                if del_size <= self.max_indel_size:
                    in_bed = self._is_indel_in_bed(ref_name, pos + 1, del_size)
                    mutations.append({
                        'position': pos + 1,  # Convert to 1-based
                        'ref': ref_seq[pos:pos+del_size].upper(),
                        'alt': '-',
                        'type': 'deletion',
                        'size': del_size,
                        'in_bed': in_bed,
                        'mutation_type': 'indel'
                    })
                pos += del_size
                hap_pos += del_size
            elif pos < len(ref_seq) and ref_seq[pos] == '-':  # Insertion
                # Calculate insertion size
                ins_size = 1
                next_pos = pos + 1
                while next_pos < len(ref_seq) and ref_seq[next_pos] == '-':
                    ins_size += 1
                    next_pos += 1
                
                if ins_size <= self.max_indel_size:
                    in_bed = self._is_indel_in_bed(ref_name, pos + 1, ins_size)
                    mutations.append({
                        'position': pos + 1,  # Convert to 1-based
                        'ref': '-',
                        'alt': haplotype[hap_pos:hap_pos+ins_size].upper(),
                        'type': 'insertion',
                        'size': ins_size,
                        'in_bed': in_bed,
                        'mutation_type': 'indel'
                    })
                hap_pos += ins_size
                pos += ins_size
            elif haplotype[hap_pos].islower():  # Substitution
                ref_base = ref_seq[pos].upper()
                alt_base = haplotype[hap_pos].upper()
                
                # Only count as SNP if it's a valid base substitution
                if self._is_valid_snp(ref_base, alt_base):
                    mutations.append({
                        'position': pos + 1,  # Convert to 1-based
                        'ref': ref_base,
                        'alt': alt_base,
                        'type': 'substitution',
                        'mutation_type': 'snp'
                    })
                pos += 1
                hap_pos += 1
            else:
                pos += 1
                hap_pos += 1

        for ins_pos, bases in self.parse_insertions(insertions):
            if len(bases) <= self.max_indel_size:
                mutations.append({
                    'position': ins_pos,
                    'ref': '-',
                    'alt': bases,
                    'type': 'insertion',
                    'size': len(bases),
                    'in_bed': self._is_indel_in_bed(ref_name, max(ins_pos, 1), 1),
                    'mutation_type': 'indel'
                })
        mutations.sort(key=lambda m: m['position'])
        
        return mutations

    def _analyze_amplicons(self,
                          amplicon_reads: List[AmpliconRead],
                          ref_name: str) -> List[Dict]:
        """Analyze processed amplicon reads."""
        results = []
        
        if not amplicon_reads:
            logger.warning(f"No valid reads found for reference {ref_name}")
            return results
        
        # Count total reads and calculate initial statistics
        haplotype_counts = Counter((read.sequence, read.insertions) for read in amplicon_reads)
        total_reads = sum(haplotype_counts.values())
        ref_seq = self.reference[ref_name].upper()
        
        logger.info(f"Found {len(haplotype_counts)} unique haplotypes from {total_reads} total reads for {ref_name}")
        
        # Filter by minimum read count
        filtered_haplotypes = {seq: count for seq, count in haplotype_counts.items() 
                             if count >= self.min_read_count}
        
        if not filtered_haplotypes:
            logger.warning(f"No haplotypes met minimum read count threshold ({self.min_read_count}) for {ref_name}")
            return results
            
        # Calculate statistics for filtered haplotypes
        filtered_total = sum(filtered_haplotypes.values())
        filtered_count = len(haplotype_counts) - len(filtered_haplotypes)
        filtered_reads = total_reads - filtered_total
        
        if filtered_count > 0:
            logger.info(
                f"Filtered out {filtered_count} low-count haplotypes "
                f"({filtered_reads:,} reads, {(filtered_reads/total_reads)*100:.1f}%) "
                f"for {ref_name}"
            )
        
        # Process each haplotype that passed the filter
        for (haplotype, insertions), count in sorted(filtered_haplotypes.items(), key=lambda x: x[1], reverse=True):
            # Calculate frequency based on total reads (not just filtered)
            frequency = (count / total_reads) * 100
            
            # Get mutation positions and types
            mutations = self._get_mutation_positions(haplotype, ref_seq, ref_name, insertions)
            
            # Count different types of mutations
            total_mutations = len(mutations)
            snp_count = sum(1 for m in mutations if m['mutation_type'] == 'snp')
            indel_count = sum(1 for m in mutations if m['mutation_type'] == 'indel')
            
            # Validate sequence length and composition
            is_full_length = self._is_full_length(haplotype, ref_seq)
            
            # Calculate theoretical maximum SNPs for this reference
            theoretical_max_snps = len(ref_seq) * 3  # 3 possible mutations per position
            
            # Add warning if we exceed theoretical maximum
            if snp_count > theoretical_max_snps:
                logger.warning(f"Haplotype has more SNPs than theoretically possible for {ref_name}: "
                             f"found={snp_count}, max={theoretical_max_snps}")
            
            results.append({
                'reference': ref_name,
                'haplotype': haplotype,
                'insertions': insertions,
                'count': count,
                'frequency': frequency,
                'mutations': total_mutations,
                'snp_count': snp_count,
                'indel_count': indel_count,
                'is_full_length': is_full_length,
                'theoretical_max_snps': theoretical_max_snps
            })
        
        return results

    def _is_indel_in_bed(self, ref_name: str, position: int, indel_size: int) -> bool:
        """Check if an indel overlaps with regions in the BED file.
        
        Args:
            ref_name: Reference sequence name
            position: 1-based position of the indel
            indel_size: Size of the indel
            
        Returns:
            bool: True if indel overlaps with BED regions
        """
        if not self.bed_regions:
            return False
            
        # Convert to 0-based position for pyranges
        pos_0based = position - 1
        
        # Create a small range for the indel
        indel_range = pr.PyRanges(
            chromosomes=[ref_name],
            starts=[pos_0based],
            ends=[pos_0based + abs(indel_size)]
        )
        
        # Check for overlap
        overlap = self.bed_regions.intersect(indel_range)
        return len(overlap) > 0

    def _downsample_fastq(self, 
                        input_fastq: Path, 
                        output_fastq: Path,
                        target_size: int) -> None:
        """Downsample a FASTQ file to approximately target size."""
        input_size = input_fastq.stat().st_size
        if input_size <= target_size:
            # If file is smaller than target, just create a symlink
            output_fastq.symlink_to(input_fastq)
            return
            
        # Calculate sampling fraction
        fraction = target_size / input_size
        
        logger.info(f"Downsampling {input_fastq.name} to {fraction:.2%} of original size")
        
        # Use seqtk to downsample
        try:
            seed = 100  # Fixed seed for reproducibility
            subprocess.run(
                [
                    'seqtk', 'sample',
                    '-s', str(seed),
                    str(input_fastq),
                    str(fraction)
                ],
                stdout=open(output_fastq, 'w'),
                stderr=subprocess.PIPE,
                check=True
            )
        except subprocess.CalledProcessError as e:
            error_msg = e.stderr.decode() if e.stderr else str(e)
            raise RuntimeError(f"Downsampling failed: {error_msg}")

    def process_sample(self, 
                      fastq_r1: Path,
                      fastq_r2: Path,
                      output_dir: Path,
                      threads: int = 4) -> pd.DataFrame:
        """Process a single sample's FASTQ files."""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir = Path(temp_dir)
            try:
                click.echo("Starting sample processing...")
                # Downsample FASTQ files if needed
                r1_size = fastq_r1.stat().st_size
                r2_size = fastq_r2.stat().st_size
                total_size = r1_size + r2_size
                
                if total_size > self.max_file_size:
                    click.echo(
                        f"Input files total size ({total_size:,} bytes) exceeds target size "
                        f"({self.max_file_size:,} bytes), downsampling..."
                    )
                    
                    # Calculate target size for each file proportionally
                    r1_target = int(self.max_file_size * (r1_size / total_size))
                    r2_target = int(self.max_file_size * (r2_size / total_size))
                    
                    # Create downsampled files
                    temp_r1 = temp_dir / "downsampled_R1.fastq"
                    temp_r2 = temp_dir / "downsampled_R2.fastq"
                    
                    with click.progressbar(length=2, label='Downsampling files') as bar:
                        self._downsample_fastq(fastq_r1, temp_r1, r1_target)
                        bar.update(1)
                        self._downsample_fastq(fastq_r2, temp_r2, r2_target)
                        bar.update(1)
                    
                    # Use downsampled files for processing
                    fastq_r1 = temp_r1
                    fastq_r2 = temp_r2
                
                click.echo("Aligning reads...")
                bam_path = self._align_reads(fastq_r1, fastq_r2, temp_dir, output_dir, threads)
                click.echo("Alignment complete. Processing references...")
                
                sample_name = fastq_r1.stem.replace("_R1_001.fastq", "")
                return self._results_from_bam(bam_path, output_dir, sample_name)
                
            except Exception as e:
                logger.error(f"Error processing sample: {str(e)}")
                return self._empty_results()

    def process_sample_ont(self,
                           fastq: Path,
                           output_dir: Path,
                           sample_name: str,
                           threads: int = 4) -> pd.DataFrame:
        """Process a single Nanopore sample's basecalled FASTQ file."""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir = Path(temp_dir)
            try:
                click.echo(f"Starting Nanopore sample processing: {sample_name}")
                bam_path = self._align_reads_ont(Path(fastq), temp_dir, output_dir, threads, sample_name)
                click.echo("Alignment complete. Processing references...")
                return self._results_from_bam(bam_path, output_dir, sample_name)
            except Exception as e:
                logger.error(f"Error processing Nanopore sample {sample_name}: {str(e)}")
                return self._empty_results()

    @staticmethod
    def _empty_results() -> pd.DataFrame:
        return pd.DataFrame(columns=['reference', 'haplotype', 'insertions', 'count', 'frequency',
                                     'mutations', 'snp_count', 'indel_count', 'is_full_length', 'theoretical_max_snps'])

    def _results_from_bam(self,
                          bam_path: Path,
                          output_dir: Path,
                          sample_name: str) -> pd.DataFrame:
        """Analyze haplotypes for every reference in an aligned BAM and write the results CSV."""
        results = []
        ref_count = len(self.reference)

        with click.progressbar(self.reference.items(),
                             length=ref_count,
                             label='Processing references') as refs:
            for ref_name, _ in refs:
                click.echo(f"\nProcessing reference: {ref_name}")
                amplicon_reads = list(self._process_alignments(bam_path, ref_name))
                if amplicon_reads:
                    click.echo(f"Found {len(amplicon_reads)} valid reads for {ref_name}")
                    results.extend(self._analyze_amplicons(amplicon_reads, ref_name))

        if not results:
            click.echo("No results found for any reference sequences")
            return self._empty_results()

        df = pd.DataFrame(results)

        # Filter by minimum read count
        df = df[df['count'] >= self.min_read_count].copy()

        if df.empty:
            click.echo(f"No haplotypes met minimum read count threshold ({self.min_read_count})")
            return self._empty_results()

        # Recalculate frequencies after filtering
        total_reads = df['count'].sum()
        df['frequency'] = (df['count'] / total_reads) * 100

        # Save results
        csv_path = output_dir / f"{sample_name}_haplotypes.csv"
        df.to_csv(csv_path, index=False)
        click.echo(f"Results saved to {csv_path}")

        return df

    def process_bam(self,
                   bam_path: Path,
                   output_dir: Path,
                   threads: int = 4) -> pd.DataFrame:
        """Process an existing BAM file.

        Args:
            bam_path: Path to the BAM file
            output_dir: Directory for output files
            threads: Number of threads to use

        Returns:
            pd.DataFrame: Results dataframe with haplotype information
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        bam_path = Path(bam_path)

        try:
            click.echo("Starting BAM processing...")
            return self._results_from_bam(bam_path, output_dir, bam_path.stem)

        except Exception as e:
            logger.error(f"Error processing BAM file: {str(e)}")
            return self._empty_results()
