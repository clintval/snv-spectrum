import re
from heapq import nlargest
from operator import itemgetter
from pprint import pformat
from typing import Any, List, Optional, Tuple

import numpy as np

__all__ = ['DictMostCommonMixin', 'DictNpArrayMixin', 'DictPrettyReprMixin']


class DictMostCommonMixin(object):
    """Give any *dict-like* object a most common method.

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
        """List the `n` most common elements and their counts.

        Method returns items from the most common to the least.
        If `n` is ``None``, then list all element counts.

        Args:
            n: The `n` most common items to return, optional.

        """
        if n is None:
            return sorted(super().items(), key=itemgetter(1), reverse=True)  # type: ignore
        return nlargest(n, super().items(), key=itemgetter(1))  # type: ignore


class DictNpArrayMixin(object):
    """Make any *dict-like* object methods return :class:`numpy.ndarray` by default.

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


class DictPrettyReprMixin(object):
    """Make any *dict-like* object pretty print when :meth:`DictPrettyReprMixin.__repr__` is called.

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
        """Pretty print a dictionary but lead with the classname."""
        indent = len(self.__class__.__qualname__) + 2
        content = re.sub(
            r'{\s*', '{', pformat(dict(super().items()), indent=indent)  # type: ignore
        )
        return f'{self.__class__.__qualname__}({content})'
