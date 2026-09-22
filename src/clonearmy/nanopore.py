"""Oxford Nanopore support: Dorado basecalling, demultiplexing and read QC."""

from pathlib import Path
from typing import Dict, Optional, Union
import logging
import re
import shutil
import subprocess

import click
import numpy as np
import pandas as pd
import pysam

logger = logging.getLogger(__name__)

BARCODE_PATTERN = re.compile(r'(barcode\d+|unclassified)', re.IGNORECASE)
FASTQ_SUFFIXES = ('.fastq', '.fq', '.fastq.gz', '.fq.gz')


def check_dorado(dorado_bin: Union[str, Path, None] = None) -> str:
    """Resolve the dorado executable and return its path."""
    candidate = str(dorado_bin) if dorado_bin else 'dorado'
    resolved = shutil.which(candidate)
    if resolved is None and dorado_bin and Path(dorado_bin).is_file():
        resolved = str(Path(dorado_bin).resolve())
    if resolved is None:
        raise RuntimeError(
            "dorado not found. Install it from https://github.com/nanoporetech/dorado "
            "and add it to PATH, or pass --dorado-bin."
        )
    try:
        result = subprocess.run([resolved, '--version'], capture_output=True, text=True)
        version = (result.stdout or result.stderr).strip().splitlines()
        if version:
            logger.info(f"Using dorado {version[-1]} ({resolved})")
    except Exception as e:
        logger.debug(f"Could not determine dorado version: {e}")
    return resolved


def find_raw_files(pod5_dir: Union[str, Path]) -> Dict[str, list]:
    """Return raw signal files found recursively under pod5_dir."""
    pod5_dir = Path(pod5_dir)
    return {
        'pod5': sorted(pod5_dir.rglob('*.pod5')),
        'fast5': sorted(pod5_dir.rglob('*.fast5')),
    }


def basecall(pod5_dir: Union[str, Path],
             out_bam: Union[str, Path],
             model: str = 'sup',
             kit_name: Optional[str] = None,
             device: Optional[str] = None,
             min_qscore: int = 10,
             models_dir: Union[str, Path, None] = None,
             recursive: bool = True,
             dorado_bin: Union[str, Path, None] = None,
             force: bool = False) -> Path:
    """Basecall raw POD5 data with dorado and write an unaligned BAM."""
    out_bam = Path(out_bam)
    out_bam.parent.mkdir(parents=True, exist_ok=True)

    if out_bam.exists() and out_bam.stat().st_size > 0 and not force:
        click.echo(f"Found existing basecalls at {out_bam}, skipping basecalling "
                   f"(use --force-basecall to redo)")
        return out_bam

    dorado = check_dorado(dorado_bin)
    cmd = [dorado, 'basecaller', model, str(pod5_dir), '--min-qscore', str(min_qscore)]
    if recursive:
        cmd.append('--recursive')
    if kit_name:
        cmd.extend(['--kit-name', kit_name])
    if device:
        cmd.extend(['--device', device])
    if models_dir:
        Path(models_dir).mkdir(parents=True, exist_ok=True)
        cmd.extend(['--models-directory', str(models_dir)])

    click.echo(f"Running dorado basecaller ({model}): {' '.join(cmd)}")
    partial_bam = out_bam.with_suffix('.partial.bam')
    try:
        with open(partial_bam, 'wb') as bam_out:
            subprocess.run(cmd, stdout=bam_out, stderr=None, check=True)
    except subprocess.CalledProcessError as e:
        partial_bam.unlink(missing_ok=True)
        raise RuntimeError(f"dorado basecaller failed with exit code {e.returncode}")
    partial_bam.replace(out_bam)
    return out_bam


def _barcode_from_path(path: Path, root: Path) -> str:
    matches = BARCODE_PATTERN.findall(str(path.relative_to(root)))
    return matches[-1].lower() if matches else path.name.split('.')[0]


def demux(calls_bam: Union[str, Path],
          out_dir: Union[str, Path],
          dorado_bin: Union[str, Path, None] = None,
          threads: int = 8,
          force: bool = False) -> Dict[str, Path]:
    """Split basecalled reads by barcode and return {barcode: fastq}."""
    out_dir = Path(out_dir)
    raw_dir = out_dir / 'dorado_demux'
    fastq_dir = out_dir / 'fastq'

    existing = sorted(fastq_dir.glob('*.fastq')) if fastq_dir.exists() else []
    if existing and not force:
        click.echo(f"Found existing demultiplexed FASTQs in {fastq_dir}, skipping demux")
        return {f.stem: f for f in existing if f.stem != 'unclassified'}

    for d in (raw_dir, fastq_dir):
        if d.exists():
            shutil.rmtree(d)
        d.mkdir(parents=True)

    dorado = check_dorado(dorado_bin)
    cmd = [dorado, 'demux', '--no-classify', '--emit-fastq',
           '--threads', str(threads), '--output-dir', str(raw_dir), str(calls_bam)]
    click.echo(f"Running dorado demux: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, stderr=None, check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"dorado demux failed with exit code {e.returncode}")

    grouped: Dict[str, list] = {}
    for f in sorted(raw_dir.rglob('*')):
        if f.is_file() and f.name.endswith(FASTQ_SUFFIXES):
            grouped.setdefault(_barcode_from_path(f, raw_dir), []).append(f)

    if not grouped:
        raise RuntimeError(f"dorado demux produced no FASTQ files in {raw_dir}")

    barcodes = {}
    with click.progressbar(grouped.items(), length=len(grouped),
                           label='Collecting demultiplexed reads') as items:
        for barcode, files in items:
            merged = fastq_dir / f"{barcode}.fastq"
            with open(merged, 'w') as out:
                for f in files:
                    with pysam.FastxFile(str(f)) as fh:
                        for rec in fh:
                            out.write(str(rec) + '\n')
            if barcode != 'unclassified':
                barcodes[barcode] = merged
    return barcodes


def bam_to_fastq(calls_bam: Union[str, Path],
                 out_fastq: Union[str, Path],
                 threads: int = 4,
                 force: bool = False) -> Path:
    """Convert an unaligned dorado BAM to FASTQ, keeping base qualities."""
    out_fastq = Path(out_fastq)
    out_fastq.parent.mkdir(parents=True, exist_ok=True)
    if out_fastq.exists() and out_fastq.stat().st_size > 0 and not force:
        click.echo(f"Found existing FASTQ at {out_fastq}, skipping conversion")
        return out_fastq
    try:
        with open(out_fastq, 'w') as out:
            subprocess.run(['samtools', 'fastq', '-@', str(threads), str(calls_bam)],
                           stdout=out, stderr=subprocess.PIPE, check=True)
    except subprocess.CalledProcessError as e:
        error_msg = e.stderr.decode() if e.stderr else str(e)
        raise RuntimeError(f"samtools fastq failed: {error_msg}")
    return out_fastq


def fastq_stats(fastq: Union[str, Path]) -> Dict[str, float]:
    """Compute read count, bases, mean length, mean read Q and N50 for a FASTQ."""
    lengths = []
    read_qs = []
    with pysam.FastxFile(str(fastq)) as fh:
        for rec in fh:
            lengths.append(len(rec.sequence))
            quals = rec.get_quality_array() if rec.quality else None
            if quals is not None and len(quals):
                err = np.power(10.0, -np.asarray(quals, dtype=float) / 10.0).mean()
                read_qs.append(-10.0 * np.log10(err))

    if not lengths:
        return {'reads': 0, 'bases': 0, 'mean_length': 0.0, 'mean_q': 0.0, 'n50': 0}

    sorted_lengths = np.sort(np.asarray(lengths))[::-1]
    cumulative = np.cumsum(sorted_lengths)
    n50 = int(sorted_lengths[np.searchsorted(cumulative, cumulative[-1] / 2.0)])
    return {
        'reads': len(lengths),
        'bases': int(cumulative[-1]),
        'mean_length': round(float(np.mean(lengths)), 1),
        'mean_q': round(float(np.mean(read_qs)), 2) if read_qs else 0.0,
        'n50': n50,
    }


def basecall_summary(fastqs: Dict[str, Path]) -> pd.DataFrame:
    """Build a per-sample basecalling summary table."""
    rows = []
    with click.progressbar(fastqs.items(), length=len(fastqs),
                           label='Computing read statistics') as items:
        for sample, fastq in items:
            rows.append({'sample': sample, **fastq_stats(fastq)})
    return pd.DataFrame(rows, columns=['sample', 'reads', 'bases', 'mean_length', 'mean_q', 'n50'])
