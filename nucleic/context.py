from collections import Counter
from pathlib import Path
from typing import Iterable, Union

from nucleic.io.bed import Interval
from nucleic import ContextWeights, Dna, Notation
from nucleic.util import kmers


class ContextUtil(object):

    @staticmethod
    def counts_from_intervals(
        intervals: Iterable[Interval],
        reference: Union[Path, str],
        notation = Notation.none,
        k: int = 3,
    ):
        weights = Counter()
        weights.update({normalize_context(Dna(_), notation = notation): 0 for _ in kmers(k, alphabet = 'ACGT')})
        for interval in intervals:
            observed = interval.reference_seq(reference).kmer_frequencies(k, overlap = True)
            weights += {normalize_context(k, notation = notation): v for k, v in observed.items()}
        return weights

def normalize_context(context: Dna, notation: Notation = Notation.none):
    context: Dna = Dna(context)
    if len(context) % 2 != 1:
        raise ValueError('Context must be of odd length.')

    middle = context[int((len(context) - 1) / 2)]
    if middle.is_purine() and notation.name == 'pyrimidine':
        return context.reverse_complement()
    elif middle.is_pyrimidine() and notation.name == 'purine':
        return context.reverse_complement()
    else:
        return context