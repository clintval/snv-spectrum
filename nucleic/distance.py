from itertools import combinations, product
from typing import Callable, Generator, List, Optional, Union

import numpy as np

from scipy.spatial.distance import pdist

try:
    from fastcluster import linkage
except ImportError:
    from scipy.cluster.hierarchy import linkage

from nucleic import SnvSpectrum

__all__ = ['hamming_circle', 'hierarchy_cluster']


def hamming_circle(
    sequence: str, n: int, alphabet: Optional[str] = None
) -> Generator[str, None, None]:
    """Find strings, of a given character alphabet, with a distance of `n` away from a sequence.

    Args:
        sequence: A sequence of characters.
        n: The Hamming distance all returned items should be to *sequence*.
        alphabet: The alphabet to use when mutating *sequence*, optional.

    Yields:
        The next sequence found with a distance of `n`.

    Examples:
        >>> sorted(hamming_circle('abc', n=0))
        ['abc']
        >>> sorted(hamming_circle('abc', n=1, alphabet='abc'))
        ['aac', 'aba', 'abb', 'acc', 'bbc', 'cbc']
        >>> sorted(hamming_circle('aaa', n=2, alphabet='ab'))
        ['abb', 'bab', 'bba']

    """
    if not isinstance(sequence, str):
        raise TypeError('`sequence` must be of type `str`.')

    if alphabet is None:
        alphabet = ''.join(sorted(set(sequence)))

    for positions in combinations(range(len(sequence)), n):
        for replacements in product(range(len(alphabet)), repeat=n):
            skip = False
            cousin = list(sequence)

            for position, replacement in zip(positions, replacements):
                if cousin[position] == alphabet[replacement]:
                    skip = True
                else:
                    cousin[position] = alphabet[replacement]

            if skip is False:
                yield ''.join(cousin)


def hierarchy_cluster(
    spectrums: List[SnvSpectrum],
    method: str = 'weighted',
    metric: Union[str, Callable] = 'cosine',
    optimal_ordering: bool = True,
) -> np.ndarray:
    """Perform hierarchical or agglomerative clustering.

    Args:
        spectrums: The spectrums to cluster.
        method: The linkage algorithm to use.
        metric: The distance metric to use.

    Returns:
        The hierarchical clustering encoded as a linkage matrix.

    Notes:
        See SciPy's `Linkage Methods <https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html>`_
        section for a full description of available methods. Also see SciPy's
        :func:`pdist <scipy.spatial.distance.pdist>` function for a list of valid distance metrics.
        A custom distance metric can also be used for a metric.

    """
    observations = [spectrum.mass() for spectrum in spectrums]
    pairwise = pdist(observations, metric=metric)
    Z: np.ndarray = linkage(pairwise, method=method, metric=metric)
    if optimal_ordering is True:
        from polo import optimal_leaf_ordering

        Z = optimal_leaf_ordering(Z, pairwise)
    return Z
