import csv
import io
import urllib.request as request
from collections import defaultdict
from itertools import combinations, product
from typing import Dict, Generator, List, Mapping, Tuple

from nucleic.constants import DNA_IUPAC_NONDEGENERATE

__all__ = [
    'STRATTON_SNV_COLOR',
    'LONGFORM_LABEL',
    'COSMIC_SIGNATURE_URL',
    'dna_kmers',
    'fetch_cosmic_signatures',
    'hamming_circle',
]

STRATTON_SNV_COLOR: Mapping[str, str] = {
    'A→C': '#EDBFC2',
    'A→G': '#97D54C',
    'A→T': '#CBC9C8',
    'C→A': '#52C3F1',
    'C→G': '#231F20',
    'C→T': '#E62223',
    'G→A': '#E62223',
    'G→C': '#231F20',
    'G→T': '#52C3F1',
    'T→A': '#CBC9C8',
    'T→C': '#97D54C',
    'T→G': '#EDBFC2',
}

DEFAULT_SNV_COLOR: Mapping[str, str] = {
    'A→C': '#D53E4F',
    'A→G': '#FC8D59',
    'A→T': '#FEE08B',
    'C→A': '#3288BD',
    'C→G': '#99D594',
    'C→T': '#E6F598',
    'G→A': '#E6F598',
    'G→C': '#99D594',
    'G→T': '#3288BD',
    'T→A': '#FEE08B',
    'T→C': '#FC8D59',
    'T→G': '#D53E4F',
}

LONGFORM_LABEL: Mapping[str, str] = {
    'A→C': 'A:T→C:G',
    'A→G': 'A:T→G:C',
    'A→T': 'A:T→T:A',
    'C→A': 'C:G→A:T',
    'C→G': 'C:G→G:C',
    'C→T': 'C:G→T:A',
    'G→A': 'G:C→A:T',
    'G→C': 'G:C→C:G',
    'G→T': 'G:C→T:A',
    'T→A': 'A:T→T:A',
    'T→C': 'A:T→G:C',
    'T→G': 'A:T→C:G',
}

COSMIC_SIGNATURE_URL = (
    'http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt'
)


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


def fetch_cosmic_signatures() -> Dict:
    """Fetch the COSMIC published signatures from the following URL.

        - https://cancer.sanger.ac.uk/cosmic

    Returns:
        signatures: The probability masses of the COSMIC signatures.

    """
    from nucleic import Dna, Spectrum, Notation

    all_signatures: defaultdict = defaultdict(lambda: Spectrum(k=3, notation=Notation.pyrimidine))

    with request.urlopen(COSMIC_SIGNATURE_URL) as handle:
        reader = csv.reader(io.TextIOWrapper(handle), delimiter='\t')

        # First three columns are subtype, column, and label titles
        _, _, _, *signature_titles = list(filter(None, next(reader)))

        for line in reader:
            subtype, context, _, *points = list(filter(None, line))
            for title, point in zip(signature_titles, map(float, points)):
                # Split the subtype to get reference and alternate
                left, right = subtype.split('>')
                snv = Dna(left).to(right).within(context)
                all_signatures[title][snv] = point

    return dict(all_signatures)


def equal_partition(string: str) -> Tuple[str, str]:
    half, remainder = divmod(len(string), 2)
    return string[: half + remainder], string[half + remainder :]
