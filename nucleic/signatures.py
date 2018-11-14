from typing import List, Mapping, Optional

from nucleic import Spectrum

__all__ = [
    'assess_number_signatures',
    'cluster_spectrum',
    'deconstruct_sigs',
    'identify_signatures',
]


def assess_number_signatures(
    signatures: List[Spectrum],
    n: int,
    decomposition: str,
    decomposition_args: Optional[Mapping] = None,
    iterations: int = 1,
) -> None:
    """Asses the number of signatures."""
    raise NotImplementedError('Function placeholder.')


def cluster_spectrum(
    signatures: List[Spectrum],
    method: str = 'weighted',
    metric: str = 'cosine',
    optimal_order: bool = True,
) -> None:
    """Cluster spectrum."""
    raise NotImplementedError('Function placeholder.')


def deconstruct_sigs(spectrum: Spectrum, signatures: List[Spectrum], method: str) -> None:
    """Deconstruct spectrums into known signatures."""
    raise NotImplementedError('Function placeholder.')


def identify_signatures(signatures: List[Spectrum], n: int) -> None:
    """Identifiy de novo signatures from spectrums."""
    raise NotImplementedError('Function placeholder.')
