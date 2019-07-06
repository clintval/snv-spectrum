from pathlib import Path, PurePath
from typing import Union

from pyfaidx import Fasta

from nucleic import Dna


def subseq(reference: Union[Fasta, Path], contig: str, start: int, end: int) -> Dna:
    """Fetch a subsequence from a FASTA file.

    Args:
        reference: FASTA reference.
        contig: The reference sequence name containing the locus.
        start: The 0-based contig start position.
        end: The end exclusive contig end position.

    Notes:
        - The FASTA file will be indexed if it is not.

    """
    if isinstance(reference, (PurePath, str)):
        reference = Fasta(str(reference))
    return Dna(str(reference[str(contig)][int(start) : int(end)]), lowercase = True)


def centered_subseq(reference: Union[Fasta, Path], contig: str, position: int, k: int = 3) -> Dna:
    """Fetch a subsequence from a FASTA file centered on a position.

    Args:
        reference: FASTA reference.
        contig: The reference sequence name containing the position.
        position: The 0-based contig position that the Variant is centered on.
        k: The length of the context, must be positive and odd.

    Notes:
        - The FASTA file will be indexed if it is not.
        - The length of the returned subsequence will be odd.

    """
    if not isinstance(position, int) and position >= 0:
        raise TypeError('position must be a postitive integer')
    if not isinstance(k, int) and k % 2 != 1 and k > 0:
        raise TypeError('k must be a positive odd integer')

    flank_length = int((k - 1) / 2)
    start, end = position - flank_length - 1, position + flank_length
    subseq: Dna = subseq(reference, contig, start, end)
    return subseq
