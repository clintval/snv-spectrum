import json
import re
from heapq import nlargest
from itertools import product
from operator import itemgetter
from pathlib import Path
from pprint import pformat
from typing import Any, Dict, Generator, Iterable, List, Optional, Tuple, Union

import numpy as np

__all__ = [
    'DictMostCommonMixin',
    'DictNpArrayMixin',
    'DictPrettyReprMixin',
    'UnreachableException',
    'dataset',
    'kmers',
]


class UnreachableException(Exception):
    """Raise an exeption when a statement was meant to be unreachable."""

    pass


class DictMostCommonMixin(object):
    """Give any *dict*-like object a "most common" method.

    Examples:
        >>> class MyDict(DictMostCommonMixin, dict):
        ...     def __init__(self, *args, **kwargs):
        ...         super().__init__(*args, **kwargs)
        >>> mapping = MyDict({'sample-1': 2, 'sample-2': 10})
        >>> mapping.most_common()
        [('sample-2', 10), ('sample-1', 2)]
        >>> mapping.most_common(n=1)
        [('sample-2', 10)]

    """

    # Static typing has no features for mixins: https://github.com/python/typing/issues/246
    def most_common(self, n: Optional[int] = None) -> List[Tuple[Any, Any]]:
        """List the `n` most common items ordered by their values.

        This method returns dictionary items ordered from the most common to the least.
        If `n` is :class:`None`, then return all ordered items.

        Args:
            n: The `n` most common items to return, optional.

        """
        if n is None:
            return sorted(super().items(), key=itemgetter(1), reverse=True)  # type: ignore
        return nlargest(n, super().items(), key=itemgetter(1))  # type: ignore


class DictNpArrayMixin(object):
    """Make any *dict*-like object have methods that return :class:`numpy.ndarray` by default.

    Examples:
        >>> class MyDict(DictNpArrayMixin, dict):
        ...     def __init__(self, *args, **kwargs):
        ...         super().__init__(*args, **kwargs)
        >>> mapping = MyDict({'sample-1': 2})
        >>> mapping.keys()
        array(['sample-1'], dtype='<U8')
        >>> mapping.values()
        array([2])

    """

    # Static typing has no features for mixins: https://github.com/python/typing/issues/246
    def keys(self) -> np.ndarray:
        """Return this dictionary's keys as a :class:`numpy.ndarray`."""
        return np.array(list(super().keys()))  # type: ignore

    # Static typing has no features for mixins: https://github.com/python/typing/issues/246
    def values(self) -> np.ndarray:
        """Return this dictionary's values as a :class:`numpy.ndarray`."""
        return np.array(list(super().values()))  # type: ignore

    # Static typing has no features for mixins: https://github.com/python/typing/issues/246
    def items(self) -> np.ndarray:
        """Return this dictionary's items as a :class:`numpy.ndarray`."""
        return np.array(list(zip(self.keys(), self.values())))  # type: ignore


class DictPrettyReprMixin(object):
    """Return a pretty formatted *dict*-like object when called by :py:func:`repr`.

    Examples:
        >>> class AReallyLongDictName(DictPrettyReprMixin, dict):
        ...     def __init__(self, *args, **kwargs):
        ...         super().__init__(*args, **kwargs)
        >>> AReallyLongDictName({
        ...     'ScientificObservation1': 1,
        ...     'ScientificObservation2': 2,
        ...     'ScientificObservation3': 3,
        ...     'ScientificObservation4': 4})
        AReallyLongDictName({'ScientificObservation1': 1,
                             'ScientificObservation2': 2,
                             'ScientificObservation3': 3,
                             'ScientificObservation4': 4})

    """

    # Static typing has no features for mixins: https://github.com/python/typing/issues/246
    def __repr__(self) -> str:
        """Return a pretty formatted dictionary, but lead with the classname."""
        indent = len(self.__class__.__qualname__) + 2
        content = re.sub(
            r'{\s*', '{', pformat(dict(super().items()), indent=indent)  # type: ignore
        )
        return f'{self.__class__.__qualname__}({content})'


def dataset(identifier: str, database: str = 'published') -> Dict:
    """Return the filesystem path to the packaged data file.

    Args:
        identifier: The ID of the dataset, usually the PubMed ID.
        database: The database to search, defaults to "published".

    Return:
        The JSON serializable structure of data.

    Examples:
        >>> from pprint import pprint
        >>> from nucleic.util import dataset
        >>> pprint(dataset('28351974'))  # doctest:+ELLIPSIS
        [{'name': 'AFB1-gpt-10wk-exposure',
          'vector': [0.0034329041,
                     0.0137316166,
                     0.0344053282,
                     0.0188924926,
        ...

    """
    import nucleic

    module_file = Path(nucleic.__file__).expanduser().resolve()
    root_directory = module_file.parent.parent
    path = root_directory / 'data' / f'{database}.json'
    assert path.exists(), f'Database "{database}" does not exist!"'
    content: Dict = json.loads(path.read_text()).get(identifier, {})
    return content


def kmers(k: int, alphabet: Union[Iterable[str], str]) -> Generator[str, None, None]:
    """Return the cartesian product of all alphabet strings of length `k`.

    Args:
        k: Length of substring.
        alphabet: The characters to use for building the kmer set.

    Yields:
        Cartesian product of all alphabet strings of length `k`.

    Examples:
        >>> list(kmers(1, alphabet='ACGT'))
        ['A', 'C', 'G', 'T']
        >>> len(list(kmers(3, alphabet='ACGT')))
        64

    """
    yield from map(lambda _: ''.join(_), product(alphabet, repeat=k))
