import csv
import io
import urllib.request as request

from collections import defaultdict
from itertools import product

from typing import Generator, Mapping, Tuple

__all__ = [
    'IUPAC_MAPPING',
    'PURINES',
    'PYRIMIDINES',
    'STRATTON_SNV_COLOR',
    'LONGFORM_LABEL',
    'COSMIC_SIGNATURE_URL',
    'dna_kmers',
    'fetch_cosmic_signatures',
]


IUPAC_MAPPING = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

PURINES: Tuple[str, str] = ('A', 'G')
PYRIMIDINES: Tuple[str, str] = ('C', 'T')

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


def complement(seq: str) -> str:
    translation_table = str.maketrans(IUPAC_MAPPING)
    return seq.translate(translation_table)


def reverse_complement(seq: str) -> str:
    return complement(seq)[::-1]


def dna_kmers(k: int = 3) -> Generator[str, None, None]:
    """Return the cartesian product of all DNA substrings of length k.

    Args:
        k: Length of of the DNA substring.

    Yields:
        Cartesian product of all DNA substrings of length k.

    Examples
        >>> list(dna_kmers(1))
        ['A', 'C', 'G', 'T']
        >>> len(list(dna_kmers(3)))
        64

    """
    for parts in product(sorted(IUPAC_MAPPING), repeat=k):
        yield ''.join(parts)


def fetch_cosmic_signatures() -> Mapping:
    """Fetch the COSMIC published signatures from:

        - https://cancer.sanger.ac.uk/cosmic

    Returns:
        cosmic_signatures: The probability masses of the COSMIC signatures.

    """
    from nucleic import Nt, Spectrum, Notation

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
                snv = Nt(left).to(right).within(context)
                all_signatures[title][snv] = point

    return dict(all_signatures)


# def kmer_frequencies_from_bed(
#     bed_file: Path,
#     reference_file: Path,
#     k: int = 3,
#     relative: bool = False,
#     reverse_complement: bool = None,
#     reference_notation: str = None,
# ) -> Counter:
#     """Count kmers in the sequences defined by intervals.

#     TODO: Use an interval tree to merge as default.

#     Args:
#         bed_file: Path to the bed intervals.
#         reference-file : Path to the reference genome.
#         k: The length of the kmer.
#         reverse_complement: To count on the reference or anti-reference strand.
#         reference_notation: ...

#     Returns:
#         kmer_dict: Mapping of kmers and their counts.

#     Note:
#         This function makes no attempt to merge overlapping intervals into one!

#     """
#     from skbio import DNA

#     if reverse_complement is not None and reference_notation is not None:
#         raise ValueError('reverse complement and reference notation can not both be set.')
#     kmer_dict = Counter()
#     reference = Fasta(str(reference_file), sequence_always_upper=True)

#     with open(bed_file) as handle:
#         for line in handle:
#             line = line.strip().split()
#             chrom, start, end, *_ = line

#             dna = DNA(str(reference[chrom][int(start) : int(end)]))
#             dna = dna.reverse_complement() if reverse_complement else dna

#             kmer_dict.update(dna.kmer_frequencies(k, relative=relative, overlap=True))

#     new = Counter()

#     for context, count in kmer_dict.items():
#         middle_base = str(context[int((k - 1) / 2)])

#         if (middle_base in PURINES and reference_notation == 'pyrimidine') or (
#             middle_base in PYRIMIDINES and reference_notation == 'purine'
#         ):
#             context = DNA(context).reverse_complement()
#         new[context] = count
#     return new


def equal_partition(string: str) -> Tuple[str, str]:
    half, remainder = divmod(len(string), 2)
    return string[: half + remainder], string[half + remainder :]
