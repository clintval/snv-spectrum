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
    raise NotImplementedError('Function placeholder.')


def cluster_spectrum(
    signatures: List[Spectrum],
    method: str = 'weighted',
    metric: str = 'cosine',
    optimal_order: bool = True,
) -> None:
    raise NotImplementedError('Function placeholder.')


def deconstruct_sigs(spectrum: Spectrum, signatures: List[Spectrum], method: str) -> None:
    raise NotImplementedError('Function placeholder.')


def identify_signatures(signatures: List[Spectrum], n: int) -> None:
    raise NotImplementedError('Function placeholder.')
