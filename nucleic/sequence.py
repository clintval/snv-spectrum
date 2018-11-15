from itertools import combinations, product
from typing import Generator, List

from nucleic.constants import DNA_IUPAC_NONDEGENERATE

__all__ = ['dna_kmers', 'hamming_circle']


def dna_kmers(k: int = 3) -> Generator[str, None, None]:
    """Return the cartesian product of all DNA substrings of length `k`.

    Args:
        k: Length of of the DNA substring.

    Yields:
        Cartesian product of all DNA substrings of length `k`.

    Examples:
        >>> list(dna_kmers(1))
        ['A', 'C', 'G', 'T']
        >>> len(list(dna_kmers(3)))
        64

    """
    for parts in product(sorted(DNA_IUPAC_NONDEGENERATE), repeat=k):
        yield ''.join(parts)


def hamming_circle(string: str, n: int, alphabet: List[str]) -> Generator[str, None, None]:
    """Find strings, of a given alphabet, with a distance of `n` away from a string.

    Examples:
        >>> sorted(hamming_circle('abc', n=0, alphabet='abc'))
        ['abc']
        >>> sorted(hamming_circle('abc', n=1, alphabet='abc'))
        ['aac', 'aba', 'abb', 'acc', 'bbc', 'cbc']
        >>> sorted(hamming_circle('aaa', n=2, alphabet='ab'))
        ['abb', 'bab', 'bba']

    """
    for positions in combinations(range(len(string)), n):
        for replacements in product(range(len(alphabet)), repeat=n):
            skip = False
            cousin = list(string)

            for position, replacement in zip(positions, replacements):
                if cousin[position] == alphabet[replacement]:
                    skip = True
                else:
                    cousin[position] = alphabet[replacement]

            if skip is False:
                yield ''.join(cousin)
