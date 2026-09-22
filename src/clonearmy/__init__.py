"""
CloneArmy: A tool for analyzing haplotypes from Illumina paired-end and
Oxford Nanopore (dorado SUP) amplicon sequencing.
"""

from pathlib import Path
from typing import Union, Dict, List
import logging
from importlib.metadata import version

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

try:
    __version__ = version("clonearmy")
except Exception:
    __version__ = "0.3.0"

# Import after setting up logging and version
from .processor import AmpliconProcessor
from .utils import (validate_input, process_samples, summarize_results,
                    validate_nanopore_input, process_nanopore_samples)

__all__ = [
    'AmpliconProcessor',
    'process_samples',
    'process_nanopore_samples',
    'summarize_results',
    'validate_input',
    'validate_nanopore_input',
    '__version__'
]
