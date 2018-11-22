from typing import List, Mapping, Optional

from nucleic import SnvSpectrum

__all__ = ['assess_number_signatures', 'deconstruct_into_signatures', 'identify_signatures']


def assess_number_signatures(
    signatures: List[SnvSpectrum],
    n: int,
    decomposition: str,
    decomposition_args: Optional[Mapping] = None,
    iterations: int = 1,
) -> None:
    """Asses the number of signatures."""
    raise NotImplementedError('Function placeholder.')


def deconstruct_into_signatures(
    spectrum: SnvSpectrum, signatures: List[SnvSpectrum], method: str
) -> None:
    """Deconstruct spectrums into known signatures."""
    raise NotImplementedError('Function placeholder.')


def identify_signatures(signatures: List[SnvSpectrum], n: int) -> None:
    """Identifiy *de novo* signatures from spectrums."""
    raise NotImplementedError('Function placeholder.')
