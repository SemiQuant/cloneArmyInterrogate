"""Read-level QC and threshold tables for processed amplicon results.

Everything here works from files a run leaves behind (the aligned BAM and the
haplotype tables), so it can be regenerated later with `clonearmy report`.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union
import json
import logging

import click
import numpy as np
import pandas as pd
import pysam
from Bio.Seq import Seq

from .processor import (AmpliconProcessor, ALL_HAPLOTYPES_SUFFIX, HAPLOTYPES_SUFFIX,
                        PARAMS_SUFFIX, HAPLOTYPE_PARAMETERS)

logger = logging.getLogger(__name__)

DEFAULT_DEPTHS = (1, 2, 3, 5, 10, 20, 50, 100)
DEFAULT_DEPTHS_STR = ','.join(str(d) for d in DEFAULT_DEPTHS)

PLATFORM_DEFAULTS = {
    'illumina': {'min_base_quality': 20, 'min_mapping_quality': 30},
    'ont': {'min_base_quality': 20, 'min_mapping_quality': 20},
}
COMMON_DEFAULTS = {
    'min_read_count': 10,
    'max_indel_size': 50,
    'ignore_homopolymer_indels': True,
    'homopolymer_min_length': 3,
    'max_homopolymer_indel': 1,
    'full_length_tolerance': 5,
}

READ_SETS = ('full_length', 'all')
_ERROR_LUT = np.power(10.0, -np.arange(256) / 10.0)


def parse_depths(value: Union[str, Sequence[int], None]) -> List[int]:
    """Parse '1,5,10' (or a sequence) into sorted unique positive depth thresholds."""
    if value is None:
        return list(DEFAULT_DEPTHS)
    items = value.split(',') if isinstance(value, str) else value
    try:
        depths = sorted({int(str(d).strip()) for d in items if str(d).strip()})
    except ValueError:
        raise click.BadParameter(f"Depth thresholds must be comma-separated integers, got {value!r}")
    if not depths or depths[0] < 1:
        raise click.BadParameter("Depth thresholds must be integers >= 1")
    return depths


@dataclass
class SampleFiles:
    """Files belonging to one processed sample."""
    name: str
    stem: str
    directory: Path
    bam: Optional[Path] = None
    haplotypes: Optional[Path] = None
    haplotypes_all: Optional[Path] = None
    parameters: Optional[Path] = None

    def cached_parameters(self) -> Optional[dict]:
        if not self.parameters or not self.parameters.exists():
            return None
        try:
            with open(self.parameters) as fh:
                return json.load(fh)
        except (OSError, json.JSONDecodeError) as e:
            logger.warning(f"Could not read {self.parameters}: {e}")
            return None


def _is_aligned_bam(path: Path) -> bool:
    try:
        with pysam.AlignmentFile(str(path), 'rb', check_sq=False) as bam:
            return len(bam.references) > 0
    except (OSError, ValueError):
        return False


def discover_samples(results_dir: Union[str, Path],
                     extra_bams: Iterable[Union[str, Path]] = (),
                     only: Optional[Iterable[str]] = None) -> Dict[str, SampleFiles]:
    """Find per-sample BAMs and haplotype tables under a results directory.

    Files are grouped by directory and sample prefix. Basecalling output is skipped.
    extra_bams are aligned BAMs stored elsewhere (e.g. inputs of `process-bam`).
    """
    results_dir = Path(results_dir)
    groups: Dict[Tuple[Path, str], SampleFiles] = {}

    def entry(directory: Path, stem: str) -> SampleFiles:
        key = (directory.resolve(), stem)
        if key not in groups:
            groups[key] = SampleFiles(name=stem, stem=stem, directory=directory)
        return groups[key]

    for path in sorted(results_dir.rglob('*')):
        if not path.is_file() or 'basecalling' in path.relative_to(results_dir).parts:
            continue
        name = path.name
        if name.endswith(ALL_HAPLOTYPES_SUFFIX):
            entry(path.parent, name[:-len(ALL_HAPLOTYPES_SUFFIX)]).haplotypes_all = path
        elif name.endswith(HAPLOTYPES_SUFFIX):
            entry(path.parent, name[:-len(HAPLOTYPES_SUFFIX)]).haplotypes = path
        elif name.endswith(PARAMS_SUFFIX):
            entry(path.parent, name[:-len(PARAMS_SUFFIX)]).parameters = path
        elif name.endswith('.bam') and not name.endswith('.partial.bam') and _is_aligned_bam(path):
            entry(path.parent, path.stem).bam = path

    for bam in extra_bams:
        bam = Path(bam)
        match = next((g for g in groups.values() if g.stem == bam.stem and g.bam is None), None)
        if match is None:
            match = entry(results_dir, bam.stem)
        match.bam = bam

    samples = [g for g in groups.values() if g.bam or g.haplotypes or g.haplotypes_all]
    for g in samples:
        if g.bam is None:
            cached = g.cached_parameters() or {}
            stored = cached.get('bam')
            if stored:
                stored = Path(stored)
                if stored.exists() and _is_aligned_bam(stored):
                    g.bam = stored
    counts: Dict[str, int] = {}
    for g in samples:
        counts[g.stem] = counts.get(g.stem, 0) + 1
    for g in samples:
        if counts[g.stem] > 1:
            try:
                rel = g.directory.resolve().relative_to(results_dir.resolve()).as_posix()
            except ValueError:
                rel = g.directory.name
            g.name = f"{rel}/{g.stem}" if rel not in ('.', g.stem) else g.stem
    if only is not None:
        only = set(only)
        samples = [g for g in samples if g.name in only or g.stem in only]
    return {g.name: g for g in sorted(samples, key=lambda g: g.name)}


def detect_platform(bam_path: Union[str, Path]) -> Optional[str]:
    """Guess the platform from the aligner in the BAM header, then from read pairing."""
    try:
        with pysam.AlignmentFile(str(bam_path), 'rb', check_sq=False) as bam:
            for pg in bam.header.to_dict().get('PG', []):
                program = f"{pg.get('ID', '')} {pg.get('PN', '')}".lower()
                if 'minimap2' in program or 'dorado' in program:
                    return 'ont'
                if 'bwa' in program:
                    return 'illumina'
            for i, read in enumerate(bam.fetch(until_eof=True)):
                if read.is_paired:
                    return 'illumina'
                if i >= 1000:
                    return 'ont'
    except (OSError, ValueError) as e:
        logger.debug(f"Could not detect platform for {bam_path}: {e}")
    return None


def _parameters_match(cached: Optional[dict], current: dict) -> bool:
    if not cached:
        return False
    return all(cached.get(k) == current.get(k) for k in HAPLOTYPE_PARAMETERS)


def _read_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if 'insertions' not in df.columns:
        df['insertions'] = ''
    df['insertions'] = df['insertions'].fillna('')
    if 'is_full_length' in df.columns:
        df['is_full_length'] = df['is_full_length'].astype(bool)
    if 'full_length_count' not in df.columns:
        df['full_length_count'] = np.nan
    return df


@dataclass
class HaplotypeSource:
    """Unfiltered (when possible) haplotype table for one sample and where it came from."""
    table: pd.DataFrame
    source: str
    unfiltered: bool


def load_haplotype_table(sample: SampleFiles,
                         processor: AmpliconProcessor,
                         reprocess: bool = False,
                         use_bam: bool = True) -> Optional[HaplotypeSource]:
    """Load the unfiltered haplotype table, rebuilding it from the BAM when needed.

    Order: cached `_haplotypes_all.csv.gz` made with the same settings, else the BAM
    (the rebuilt table is cached), else a cached table made with other settings,
    else the filtered `_haplotypes.csv` (thresholds below its filter can't be evaluated).
    """
    current = processor.processing_parameters()
    cached_ok = _parameters_match(sample.cached_parameters(), current)

    if sample.haplotypes_all and cached_ok and not reprocess:
        return HaplotypeSource(_read_table(sample.haplotypes_all), 'cached unfiltered table', True)

    if use_bam and sample.bam:
        click.echo(f"\nRebuilding haplotypes for {sample.name} from {sample.bam}")
        table = processor.haplotype_table(sample.bam)
        path = processor.write_haplotype_table(table, sample.directory, sample.stem,
                                              bam_path=sample.bam)
        sample.haplotypes_all = path
        sample.parameters = sample.directory / f"{sample.stem}{PARAMS_SUFFIX}"
        return HaplotypeSource(_read_table(path), 'rebuilt from BAM', True)

    if sample.haplotypes_all:
        logger.warning(f"{sample.name}: cached unfiltered table was made with different settings "
                       f"and no BAM is available to rebuild it; using it as is")
        return HaplotypeSource(_read_table(sample.haplotypes_all),
                               'cached unfiltered table (different settings)', True)

    if sample.haplotypes:
        table = _read_table(sample.haplotypes)
        extras = []
        if not _has_full_length(table):
            extras.append("per-read full-length status needs the BAM")
        extras.append("depth thresholds below the table's minimum count are skipped")
        logger.warning(f"{sample.name}: only the filtered haplotype table is available; "
                       f"{'; '.join(extras)}")
        return HaplotypeSource(table, 'filtered haplotypes.csv', False)

    return None


def _has_full_length(df: pd.DataFrame) -> bool:
    return 'full_length_count' in df.columns and not df.empty and df['full_length_count'].notna().all()


def _read_set_counts(df: pd.DataFrame, read_set: str) -> pd.Series:
    return df['full_length_count'].astype(int) if read_set == 'full_length' else df['count'].astype(int)


def _available_read_sets(df: pd.DataFrame) -> List[str]:
    return list(READ_SETS) if _has_full_length(df) else ['all']


def is_coding(ref_seq: str) -> bool:
    """Treat the reference as an in-frame CDS if it is a whole number of codons without internal stops."""
    ref_seq = ref_seq.upper()
    if len(ref_seq) < 6 or len(ref_seq) % 3 or set(ref_seq) - set('ACGT'):
        return False
    return '*' not in str(Seq(ref_seq).translate())[:-1]


def _translate(codon: str) -> str:
    return str(Seq(codon).translate())


def reachable_protein_changes(ref_seq: str) -> Dict[str, int]:
    """Distinct missense substitutions and nonsense codons reachable by one nucleotide change."""
    ref_seq = ref_seq.upper()
    missense, nonsense = set(), set()
    for ci in range(len(ref_seq) // 3):
        codon = ref_seq[ci * 3:ci * 3 + 3]
        ref_aa = _translate(codon)
        for offset in range(3):
            for base in 'ACGT':
                if base == codon[offset]:
                    continue
                alt_aa = _translate(codon[:offset] + base + codon[offset + 1:])
                if alt_aa == ref_aa:
                    continue
                if alt_aa == '*':
                    nonsense.add(ci)
                else:
                    missense.add((ci, alt_aa))
    return {'missense': len(missense), 'nonsense': len(nonsense)}


def _annotate_snv(ref_seq: str, position: int, alt: str) -> Dict[str, str]:
    ci = (position - 1) // 3
    offset = (position - 1) % 3
    codon = ref_seq[ci * 3:ci * 3 + 3]
    alt_codon = codon[:offset] + alt + codon[offset + 1:]
    ref_aa, alt_aa = _translate(codon), _translate(alt_codon)
    if ref_aa == alt_aa:
        consequence = 'synonymous'
    elif alt_aa == '*':
        consequence = 'nonsense'
    elif ref_aa == '*':
        consequence = 'stop_lost'
    else:
        consequence = 'missense'
    return {'codon': ci + 1, 'ref_codon': codon, 'alt_codon': alt_codon,
            'aa_change': f"{ref_aa}{ci + 1}{alt_aa}", 'consequence': consequence}


def single_variant_table(df: pd.DataFrame, processor: AmpliconProcessor,
                         ref_name: str) -> pd.DataFrame:
    """One row per distinct single-mutation variant (haplotypes with exactly one mutation)."""
    columns = ['reference', 'type', 'position', 'ref', 'alt', 'reads', 'full_length_reads',
               'codon', 'ref_codon', 'alt_codon', 'aa_change', 'consequence']
    ref_seq = processor.reference[ref_name].upper()
    sub = df[(df['reference'] == ref_name) & (df['mutations'] == 1)]
    coding = is_coding(ref_seq)
    has_fl = _has_full_length(df)
    rows: Dict[Tuple, Dict] = {}
    for hap, ins, count, fl in zip(sub['haplotype'], sub['insertions'], sub['count'],
                                   sub['full_length_count']):
        muts = processor._get_mutation_positions(hap, ref_seq, ref_name, ins)
        if len(muts) != 1:
            continue
        m = muts[0]
        kind = 'snv' if m['mutation_type'] == 'snp' else m['type']
        key = (kind, m['position'], m['ref'], m['alt'])
        row = rows.get(key)
        if row is None:
            row = {'reference': ref_name, 'type': kind, 'position': m['position'],
                   'ref': m['ref'], 'alt': m['alt'], 'reads': 0,
                   'full_length_reads': 0 if has_fl else np.nan}
            if coding and kind == 'snv':
                row.update(_annotate_snv(ref_seq, m['position'], m['alt']))
            rows[key] = row
        row['reads'] += int(count)
        if has_fl:
            row['full_length_reads'] += int(fl)
    out = pd.DataFrame(list(rows.values()), columns=columns)
    return out.sort_values(['position', 'type', 'alt']).reset_index(drop=True)


def depth_threshold_table(df: pd.DataFrame, variants: pd.DataFrame, processor: AmpliconProcessor,
                          ref_name: str, depths: Sequence[int], sample: str,
                          count_floor: int = 1) -> pd.DataFrame:
    """Single-mutant coverage at each minimum read depth, for full-length and all reads.

    A haplotype or variant passes a threshold when at least that many reads support it.
    Thresholds below count_floor (the filter already applied to the table) are skipped.
    """
    ref_seq = processor.reference[ref_name].upper()
    ref_len = len(ref_seq)
    coding = is_coding(ref_seq)
    reachable = reachable_protein_changes(ref_seq) if coding else None
    ref_df = df[df['reference'] == ref_name]
    rows = []
    for read_set in _available_read_sets(ref_df):
        counts = _read_set_counts(ref_df, read_set)
        total = int(counts.sum())
        var_counts = variants['full_length_reads' if read_set == 'full_length' else 'reads']
        for depth in depths:
            if depth < count_floor:
                continue
            keep = counts >= depth
            reads = int(counts[keep].sum())
            v = variants[var_counts >= depth]
            snvs = v[v['type'] == 'snv']
            row = {
                'sample': sample, 'reference': ref_name, 'read_set': read_set, 'min_depth': depth,
                'total_reads': total,
                'haplotypes': int(keep.sum()),
                'reads': reads,
                'reads_pct': round(reads / total * 100, 2) if total else 0.0,
                'wild_type_reads': int(counts[keep & (ref_df['mutations'] == 0)].sum()),
                'single_mut_variants': len(v),
                'single_mut_reads': int(var_counts[var_counts >= depth].sum()),
                'single_mut_pct': round(int(var_counts[var_counts >= depth].sum()) / total * 100, 2)
                if total else 0.0,
                'single_snvs': len(snvs),
                'snv_pct_of_possible': round(len(snvs) / (3 * ref_len) * 100, 2) if ref_len else 0.0,
                'snv_positions': int(snvs['position'].nunique()),
                'position_pct': round(snvs['position'].nunique() / ref_len * 100, 2) if ref_len else 0.0,
                'single_deletions': int((v['type'] == 'deletion').sum()),
                'single_insertions': int((v['type'] == 'insertion').sum()),
            }
            if coding:
                missense = snvs[snvs['consequence'] == 'missense']
                nonsense = snvs[snvs['consequence'] == 'nonsense']
                n_missense = missense['aa_change'].nunique()
                n_nonsense = nonsense['codon'].nunique()
                row.update({
                    'synonymous_snvs': int((snvs['consequence'] == 'synonymous').sum()),
                    'missense_aa_changes': int(n_missense),
                    'missense_pct_of_reachable': round(n_missense / reachable['missense'] * 100, 2)
                    if reachable['missense'] else 0.0,
                    'nonsense_codons': int(n_nonsense),
                    'nonsense_pct_of_reachable': round(n_nonsense / reachable['nonsense'] * 100, 2)
                    if reachable['nonsense'] else 0.0,
                })
            rows.append(row)
    return pd.DataFrame(rows)


def mutation_load_table(df: pd.DataFrame, sample: str, ref_name: str,
                        max_bin: int = 5) -> pd.DataFrame:
    """Distribution of mutations per read, weighted by read counts."""
    ref_df = df[df['reference'] == ref_name]
    rows = []
    for read_set in _available_read_sets(ref_df):
        counts = _read_set_counts(ref_df, read_set)
        total = int(counts.sum())
        binned = counts.groupby(ref_df['mutations'].clip(upper=max_bin)).sum()
        row = {'sample': sample, 'reference': ref_name, 'read_set': read_set, 'reads': total,
               'mean_mutations': round(float((ref_df['mutations'] * counts).sum() / total), 3)
               if total else 0.0}
        for n in range(max_bin + 1):
            label = f"pct_{n}plus_mut" if n == max_bin else f"pct_{n}_mut"
            row[label] = round(float(binned.get(n, 0)) / total * 100, 2) if total else 0.0
        rows.append(row)
    return pd.DataFrame(rows)


def bam_read_stats(bam_path: Union[str, Path], processor: AmpliconProcessor,
                   label: str = '') -> Tuple[Dict[str, int], Dict[str, Dict]]:
    """Per-sample input/unmapped counts and per-reference read statistics from an aligned BAM.

    Counts are templates: reads for Nanopore, read pairs (read 1) for Illumina.
    """
    min_mapq = processor.min_mapping_quality
    sample_stats = {'input_reads': 0, 'unmapped': 0}
    per_ref = {ref: {'mapped': 0, 'supplementary': 0, 'low_mapq': 0, 'not_proper_pair': 0,
                     'lengths': [], 'read_q': [], 'identity': []}
               for ref in processor.reference}

    with pysam.AlignmentFile(str(bam_path), 'rb') as bam:
        try:
            total = sum(s.total for s in bam.get_index_statistics()) + bam.nocoordinate
        except (ValueError, AttributeError):
            total = None
        with click.progressbar(bam.fetch(until_eof=True), length=total,
                               label=f'Read statistics {label}'.strip()) as reads:
            for read in reads:
                if read.is_secondary:
                    continue
                if read.is_supplementary:
                    if read.reference_name in per_ref:
                        per_ref[read.reference_name]['supplementary'] += 1
                    continue
                if read.is_paired and not read.is_read1:
                    continue
                sample_stats['input_reads'] += 1
                if read.is_unmapped:
                    sample_stats['unmapped'] += 1
                    continue
                st = per_ref.get(read.reference_name)
                if st is None:
                    continue
                st['mapped'] += 1
                if read.mapping_quality < min_mapq:
                    st['low_mapq'] += 1
                if read.is_paired and not read.is_proper_pair:
                    st['not_proper_pair'] += 1
                st['lengths'].append(read.query_length)
                quals = read.query_qualities
                if quals is not None and len(quals):
                    err = _ERROR_LUT[np.frombuffer(quals, dtype=np.uint8)].mean()
                    st['read_q'].append(-10.0 * np.log10(err))
                if read.has_tag('NM') and read.cigartuples:
                    cols = sum(n for op, n in read.cigartuples if op in (0, 1, 2, 7, 8))
                    if cols:
                        st['identity'].append(1.0 - read.get_tag('NM') / cols)

        is_ont = processor.platform == 'ont'

        def passes(read) -> bool:
            return not (read.is_secondary or read.is_supplementary or read.is_unmapped or
                        read.mapping_quality < min_mapq or
                        (not is_ont and read.is_paired and not read.is_proper_pair))

        for ref, st in per_ref.items():
            st['depth'] = None
            if st['mapped'] == 0:
                continue
            try:
                cov = bam.count_coverage(ref, 0, len(processor.reference[ref]),
                                         quality_threshold=processor.min_base_quality,
                                         read_callback=passes)
                st['depth'] = np.sum(np.asarray(cov), axis=0)
            except ValueError as e:
                logger.debug(f"Could not compute depth for {ref} in {bam_path}: {e}")
    return sample_stats, per_ref


def _nan_stat(values, fn, digits: int = 1):
    return round(float(fn(values)), digits) if len(values) else np.nan


def read_filtering_rows(sample: str, df: pd.DataFrame, processor: AmpliconProcessor,
                        bam_path: Optional[Path]) -> List[Dict]:
    """Read filtering funnel per reference: input -> mapped -> analysed -> full length."""
    min_count = processor.min_read_count
    sample_stats, per_ref = (bam_read_stats(bam_path, processor, sample) if bam_path
                             else ({}, {}))
    has_fl = _has_full_length(df)
    refs = [r for r in processor.reference
            if (r in df['reference'].values) or per_ref.get(r, {}).get('mapped')]
    rows = []
    for ref in refs:
        ref_df = df[df['reference'] == ref]
        analysed = int(ref_df['count'].sum())
        full_length = int(ref_df['full_length_count'].sum()) if has_fl else np.nan
        fl_pass = (int(ref_df.loc[ref_df['count'] >= min_count, 'full_length_count'].sum())
                   if has_fl else np.nan)
        st = per_ref.get(ref, {})
        input_reads = sample_stats.get('input_reads', np.nan)
        row = {
            'sample': sample,
            'reference': ref,
            'input_reads': input_reads,
            'unmapped': sample_stats.get('unmapped', np.nan),
            'mapped': st.get('mapped', np.nan),
            'supplementary': st.get('supplementary', np.nan),
            'low_mapq': st.get('low_mapq', np.nan),
        }
        if processor.platform == 'illumina':
            row['not_proper_pair'] = st.get('not_proper_pair', np.nan)
        depth = st.get('depth')
        row.update({
            'ref_length': len(processor.reference[ref]),
            'analysed': analysed,
            'analysed_pct_of_input': round(analysed / input_reads * 100, 2)
            if isinstance(input_reads, (int, np.integer)) and input_reads else np.nan,
            'full_length': full_length,
            'full_length_pct': round(full_length / analysed * 100, 2)
            if has_fl and analysed else np.nan,
            f'full_length_in_haplotypes_ge{min_count}': fl_pass,
            'mean_read_length': _nan_stat(st.get('lengths', []), np.mean),
            'median_read_length': _nan_stat(st.get('lengths', []), np.median),
            'mean_read_q': _nan_stat(st.get('read_q', []), np.mean, 2),
            'median_identity_pct': round(float(np.median(st['identity'])) * 100, 2)
            if st.get('identity') else np.nan,
            'depth_min': int(depth.min()) if depth is not None else np.nan,
            'depth_median': int(np.median(depth)) if depth is not None else np.nan,
            'depth_mean': round(float(depth.mean()), 1) if depth is not None else np.nan,
        })
        rows.append(row)
    return rows


@dataclass
class QCResult:
    read_filtering: pd.DataFrame = field(default_factory=pd.DataFrame)
    mutation_load: pd.DataFrame = field(default_factory=pd.DataFrame)
    depth_thresholds: pd.DataFrame = field(default_factory=pd.DataFrame)
    single_variants: pd.DataFrame = field(default_factory=pd.DataFrame)
    basecalling: pd.DataFrame = field(default_factory=pd.DataFrame)
    sources: pd.DataFrame = field(default_factory=pd.DataFrame)
    notes: List[str] = field(default_factory=list)

    FILES = {
        'read_filtering': 'qc_read_filtering.csv',
        'mutation_load': 'qc_mutation_load.csv',
        'depth_thresholds': 'qc_depth_thresholds.csv',
        'single_variants': 'qc_single_variants.csv',
    }

    def write(self, output_dir: Union[str, Path]) -> List[Path]:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        written = []
        for attr, filename in self.FILES.items():
            df = getattr(self, attr)
            if not df.empty:
                path = output_dir / filename
                df.to_csv(path, index=False)
                written.append(path)
        return written


def build_qc(samples: Dict[str, SampleFiles],
             tables: Dict[str, HaplotypeSource],
             processor: AmpliconProcessor,
             depths: Sequence[int] = DEFAULT_DEPTHS,
             read_stats: bool = True,
             basecalling: Optional[pd.DataFrame] = None) -> QCResult:
    """Compute all QC tables for the samples that have a haplotype table."""
    qc = QCResult()
    filtering, load, thresholds, variants, sources = [], [], [], [], []
    for name, hs in tables.items():
        df = hs.table
        sf = samples.get(name)
        floor = 1 if hs.unfiltered else int(df['count'].min()) if not df.empty else 1
        sources.append({'sample': name, 'source': hs.source,
                        'bam': str(sf.bam) if sf and sf.bam else ''})
        if not hs.unfiltered and any(d < floor for d in depths):
            qc.notes.append(f"{name}: table was pre-filtered at >= {floor} reads, so depth "
                            f"thresholds below {floor} are not shown")
        if not _has_full_length(df):
            qc.notes.append(f"{name}: no per-read full-length status (needs the BAM); "
                            f"only the 'all reads' set is reported")
        filtering.extend(read_filtering_rows(name, df, processor,
                                             sf.bam if (read_stats and sf and sf.bam) else None))
        for ref in [r for r in processor.reference if r in df['reference'].values]:
            var = single_variant_table(df, processor, ref)
            if not var.empty:
                variants.append(var.assign(sample=name))
            thresholds.append(depth_threshold_table(df, var, processor, ref, depths, name, floor))
            load.append(mutation_load_table(df, name, ref))

    qc.read_filtering = pd.DataFrame(filtering)
    qc.mutation_load = pd.concat(load, ignore_index=True) if load else pd.DataFrame()
    qc.depth_thresholds = pd.concat(thresholds, ignore_index=True) if thresholds else pd.DataFrame()
    if variants:
        v = pd.concat(variants, ignore_index=True)
        qc.single_variants = v[['sample'] + [c for c in v.columns if c != 'sample']]
    qc.sources = pd.DataFrame(sources)
    if basecalling is not None:
        qc.basecalling = basecalling
    return qc


def resolve_parameters(explicit: Dict[str, object],
                       samples: Dict[str, SampleFiles]) -> Tuple[Dict[str, object], Dict[str, str]]:
    """Fill unset report parameters from the stored run parameters, then platform defaults.

    Returns (parameters, source of each parameter).
    """
    cached = next((p for p in (s.cached_parameters() for s in samples.values()) if p), None)
    params, origin = {}, {}

    platform = explicit.get('platform')
    if platform:
        origin['platform'] = 'command line'
    elif cached and cached.get('platform'):
        platform, origin['platform'] = cached['platform'], 'stored run parameters'
    else:
        detected = {detect_platform(s.bam) for s in samples.values() if s.bam} - {None}
        if len(detected) > 1:
            logger.warning(f"Samples look like different platforms ({', '.join(sorted(detected))}); "
                           f"pass --platform to choose")
        platform = sorted(detected)[0] if detected else 'illumina'
        origin['platform'] = 'detected from BAM' if detected else 'default'
    params['platform'] = platform.lower()

    defaults = {**COMMON_DEFAULTS, **PLATFORM_DEFAULTS[params['platform']]}
    for key, default in defaults.items():
        if explicit.get(key) is not None:
            params[key], origin[key] = explicit[key], 'command line'
        elif cached and cached.get(key) is not None:
            params[key], origin[key] = cached[key], 'stored run parameters'
        else:
            params[key], origin[key] = default, 'default'
    return params, origin
