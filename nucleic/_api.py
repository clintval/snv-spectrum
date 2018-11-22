from collections import OrderedDict
from enum import Enum
from itertools import permutations
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional, Set, Tuple, Union

import numpy as np

from pyfaidx import Fasta

from skbio.sequence import GrammaredSequence
from skbio.sequence._nucleotide_mixin import NucleotideMixin
from skbio.util import classproperty

from nucleic.util import DictMostCommonMixin, DictNpArrayMixin, DictPrettyReprMixin, kmers


__all__ = ['Dna', 'Notation', 'Variant', 'SnvSpectrum']

#: The non-degenerate IUPAC DNA bases.
DNA_IUPAC_NONDEGENERATE: str = 'ACGT'


class Notation(Enum):
    """The spectrum notation types."""

    none: int = 0
    pyrimidine: int = 1
    purine: int = 2


class Dna(GrammaredSequence, NucleotideMixin):
    """Deoxyribonucleic acid composed of the following nucleotide sequences.

    ============ ==================================================== ==========
    String        Residue                                             Class
    ============ ==================================================== ==========
    :class:`A`   `Adenine <https://en.wikipedia.org/wiki/Adenine>`_   Purine
    :class:`C`   `Cytosine <https://en.wikipedia.org/wiki/Cytosine>`_ Pyrimidine
    :class:`G`   `Guanine <https://en.wikipedia.org/wiki/Guanine>`_   Purine
    :class:`T`   `Thymine <https://en.wikipedia.org/wiki/Thymine>`_   Pyrimidine
    ============ ==================================================== ==========

    Examples:
        >>> dna = Dna("A")
        >>> dna.is_purine()
        True
        >>> dna.complement()
        Dna("T")
        >>> Dna("T").to("A")
        Variant(ref=Dna("T"), alt=Dna("A"), context=Dna("T"))

    """

    @classproperty
    def degenerate_map(self) -> Dict[str, Set[str]]:
        """Return the mapping of none degenerate nucleotides."""
        return {'.': set(DNA_IUPAC_NONDEGENERATE)}

    @classproperty
    def complement_map(self) -> Dict[str, str]:
        """Return the complement mapping of none degenerate nucleotides."""
        return dict(zip(DNA_IUPAC_NONDEGENERATE, reversed(DNA_IUPAC_NONDEGENERATE)))

    @classproperty
    def definite_chars(self) -> Set[str]:
        """Return the definite characters of degenerate nucleotides."""
        return set(DNA_IUPAC_NONDEGENERATE)

    @classproperty
    def default_gap_char(self) -> str:
        """Return the default gap character."""
        return '-'

    @classproperty
    def gap_chars(self) -> Set[str]:
        """Return the set of allowed gap characters."""
        return set('-')

    def is_purine(self) -> bool:
        """Return if this sequence is a purine."""
        return str(self) in ('A', 'G')

    def is_pyrimidine(self) -> bool:
        """Return if this sequence is a pyrimdine."""
        return str(self) in ('C', 'T')

    def to(self, other: Union[str, 'Dna']) -> 'Variant':
        """Create a variant allele."""
        if isinstance(other, str):
            other = Dna(other)
        return Variant(self, other)

    def __hash__(self) -> int:
        return hash(repr(self))

    def __repr__(self) -> str:
        return f'{self.__class__.__qualname__}("{str(self)}\")'


class Variant(object):
    """A variant of DNA with both a reference and alternate allele."""

    def __init__(
        self,
        ref: Union[str, Dna],
        alt: Union[str, Dna],
        context: Optional[Union[str, Dna]] = None,
        data: Optional[Dict] = None,
    ) -> None:
        self.ref = ref
        self.alt = alt
        self.context = context
        self.data = data

    def color(self, scheme: str = 'default') -> Optional[str]:
        """Return the color representing this class of variant."""
        from nucleic.plotting.cmap import SnvCmap

        if self.is_snv():
            mapping: Dict[str, str] = getattr(SnvCmap, scheme)
            return mapping[self.label()]
        return None

    @property
    def context(self) -> Dna:
        """Return the context of this variant."""
        return self._context

    @context.setter
    def context(self, context: Optional[Union[Dna, str]]) -> None:
        """Verify that the context is valid.

        Args:
            context: The context of this variant, default to `ref`.

        """
        if context is None:
            self._context = self.ref
            return None
        elif isinstance(context, str):
            context = Dna(context)
        if len(context) % 2 != 1:
            raise ValueError(f'Context must be of odd length: {context}')

        middle = context[int(len(context) / 2)]
        if str(middle) != str(self.ref):
            raise ValueError(f'Middle of context must equal ref: {middle} != {self.ref}')

        self._context = context

    @property
    def ref(self) -> Dna:
        """Return the reference of this variant."""
        return self._ref

    @ref.setter
    def ref(self, ref: Dna) -> None:
        """Verify that the reference is of type :class:`nucleic.Dna`."""
        if isinstance(ref, str):
            ref = Dna(ref)
        self._ref = ref

    @property
    def alt(self) -> Dna:
        """Return the alternate of this variant."""
        return self._alt

    @alt.setter
    def alt(self, alt: Dna) -> None:
        """Verify that the alternate is of type :class:`nucleic.Dna`."""
        if isinstance(alt, str):
            alt = Dna(alt)
        self._alt = alt

    def complement(self) -> 'Variant':
        """Return the complement variant."""
        return self.copy(
            ref=self.ref.complement(), alt=self.alt.complement(), context=self.context.complement()
        )

    def reverse_complement(self) -> 'Variant':
        """Return the reverse complement of this variant."""
        return self.copy(
            ref=self.ref.reverse_complement(),
            alt=self.alt.reverse_complement(),
            context=self.context.reverse_complement(),
        )

    def is_null(self) -> bool:
        """Return if this variant has an alternate equivalent to its reference."""
        return self.ref == self.alt

    def is_snv(self) -> bool:
        """Return if this variant is a single nucleotide variant."""
        return not self.is_null() and len(self.ref) == len(self.alt) == 1

    def is_transition(self) -> bool:
        """Return if this variant is a transition."""
        return (
            self.is_snv()
            and (self.ref.is_pyrimidine() and self.alt.is_pyrimidine())
            or (self.ref.is_purine() and self.alt.is_purine())
        )

    def is_transversion(self) -> bool:
        """Return if this variant is a transversion."""
        return self.is_snv() and not self.is_transition()

    def is_mnv(self) -> bool:
        """Return if this variant is a multiple nucleotide variant."""
        return not self.is_snv() and len(self.ref) == len(self.alt)

    def is_insertion(self) -> bool:
        """Return if this variant is an insertion."""
        return len(self.ref) < len(self.alt)

    def is_deletion(self) -> bool:
        """Return if this variant is a deletion."""
        return len(self.ref) > len(self.alt)

    def is_indel(self) -> bool:
        """Return if this variant is an insertion or deletion."""
        return self.is_insertion() or self.is_deletion()

    def lseq(self) -> Dna:
        """Return the 5′ adjacent sequence to the variant."""
        return Dna(self.context[0 : int((len(self.context) - 1) / 2)])

    def rseq(self) -> Dna:
        """Return the 3′ adjacent sequence to the variant."""
        return Dna(self.context[int((len(self.context) - 1) / 2) + 1 :])

    def label(self) -> str:
        """Return a pretty representation of the variant."""
        return '→'.join(map(str, [self.ref, self.alt]))

    def copy(self, **kwargs: Any) -> 'Variant':
        """Make a copy of this variant."""
        kwargs = {} if kwargs is None else kwargs
        return Variant(
            ref=kwargs.pop('ref', self.ref),
            alt=kwargs.pop('alt', self.alt),
            context=kwargs.pop('context', self.context),
            data=kwargs.pop('data', self.data),
        )

    def within(self, context: Union[str, Dna]) -> 'Variant':
        """Set the context of this variant."""
        if isinstance(context, str):
            context = Dna(context)
        self.context = context
        return self

    def with_purine_ref(self) -> 'Variant':
        """Return this variant with its reference as a purine."""
        return self.reverse_complement() if self.ref.is_pyrimidine() else self

    def with_pyrimidine_ref(self) -> 'Variant':
        """Return this variant with its reference as a pyrimidine."""
        return self.reverse_complement() if self.ref.is_purine() else self

    def set_context_from_fasta(self, infile: Path, contig: str, position: int, k: int = 3) -> None:
        """Set the context by looking up a genomic loci from a FASTA.

        Args:
            infile: FASTA filepath.
            contig: The contig name with containing the locus.
            position: The 0-based contig position that the Variant is centered on.
            k: The length of the context, must be positive and odd.

        Notes:
            - The FASTA file will be indexed if it is not.
            - The length of the context must be odd so the context can be
              symmetrical in length about the the target position.

        """
        reference = Fasta(str(infile), sequence_always_upper=True)

        if not isinstance(position, int) and position >= 0:
            raise TypeError('position must be a postitive integer')
        if not isinstance(contig, str):
            raise TypeError('contig must be of type str')
        if not isinstance(k, int) and k % 2 != 1 and k > 0:
            raise TypeError('k must be a positive odd integer')

        flank_length = (k - 1) / 2
        start, end = position - flank_length - 1, position + flank_length
        self.context = Dna(reference[contig][int(start) : int(end)])

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, Variant):
            return NotImplemented
        return (
            self.ref == other.ref
            and self.alt == other.alt
            and str(self.context) == str(self.context)
        )

    def __hash__(self) -> int:
        return hash(f'{self.ref}{self.alt}{self.context}')

    def __len__(self) -> int:
        return len(self.context)

    def __repr__(self) -> str:
        return (
            f'{self.__class__.__qualname__}('
            f'ref={repr(self.ref)}, '
            f'alt={repr(self.alt)}, '
            f'context={repr(self.context)})'
        )

    def __str__(self) -> str:
        return f'{self.lseq()}[{self.label()}]{self.rseq()}'


class ContextWeights(DictPrettyReprMixin, DictMostCommonMixin, DictNpArrayMixin, OrderedDict):
    """A dictionary of sequences and their respective weights."""

    def __new__(cls, data: Optional[Any] = None, **kwargs: Mapping[Any, Any]) -> Any:
        return super().__new__(cls, data, **kwargs)


class Spectrum(DictMostCommonMixin, DictNpArrayMixin, OrderedDict):
    pass


class IndelSpectrum(Spectrum):
    """A spectrum of indel variants of various sizes.

    Warning:
        Not implemented.

    """

    def __init__(self) -> None:
        raise NotImplementedError('Class placeholder.')


class SnvSpectrum(Spectrum):
    """A spectrum of single nucleotide variants."""

    def __init__(self, k: int = 3, notation: Notation = Notation.none) -> None:
        self._initialized: bool = False

        if not isinstance(k, int) and k % 2 != 1 and k > 0:
            raise TypeError('`k` must be a positive odd integer')
        elif not isinstance(notation, Notation):
            raise TypeError('`notation` must be of type Notation')

        self.k = k
        self.notation = notation
        self.weights = ContextWeights()

        # Reverse the order of a `SnvSpectrum` built with purine notation.
        if notation.name == 'purine':
            codes = list(reversed(DNA_IUPAC_NONDEGENERATE))
        else:
            codes = list(DNA_IUPAC_NONDEGENERATE)

        for ref, alt in permutations(map(Dna, codes), 2):
            if (
                ref.is_purine()
                and notation.name == 'pyrimidine'
                or ref.is_pyrimidine()
                and notation.name == 'purine'
            ):
                continue
            for context in kmers(k, alphabet=sorted(DNA_IUPAC_NONDEGENERATE)):
                if context[int((k - 1) / 2)] != str(ref):
                    continue

                variant = ref.to(alt)
                self.weights[context] = 1
                self[variant.within(context)] = 0
        self._initialized = True

    @classmethod
    def from_iterable(
        cls, iterable: Iterable[int], k: int = 1, notation: Notation = Notation.none
    ) -> 'SnvSpectrum':
        """Build a spectrum from an interable of numbers."""
        spectrum = cls(k=k, notation=notation)
        iterable = list(iterable)

        if len(spectrum) != len(iterable):
            raise ValueError('`iterable` not of `len(SnvSpectrum())`')

        for snv, count in zip(spectrum.keys(), iterable):
            spectrum[snv] = count
        return spectrum

    def mass(self) -> np.ndarray:
        """Return the discrete probability mass of this spectrum.

        Raises:
            ValueError: if an observation is found with zero context weight.

        """
        # Normalized proportion for a variant within this spectrum.
        def norm_count(key: Variant) -> float:
            if self[key] != 0 and self.weights[str(key.context)] == 0:
                raise ValueError('Observations with no weight found: {self[key]}')
            return float(self[key] / self.weights[str(key.context)])

        array = np.array([norm_count(key) for key in self.keys()])
        return array / array.sum()

    def split_by_notation(self) -> Tuple['SnvSpectrum', 'SnvSpectrum']:
        """Split pyrimidine *vs* purine reference variants into seperate spectrum.

        Raises:
            ValueError: if the ``notation`` of this spectrum is not :class:`Notation.none`.

        Returns:
            spectrum_pu: A :class:`SnvSpectrum` holding purine reference variants.
            spectrum_py: A :class:`SnvSpectrum` holding pyrimidine reference variants.

        Note:
            - TODO: Return a collection holding the two spectrum, like ``namedtuple``.

        """
        if self.notation.name != 'none':
            raise ValueError('`SnvSpectrum` notation must be `Notation.none`')

        spectrum_pu = SnvSpectrum(self.k, Notation.purine)
        spectrum_py = SnvSpectrum(self.k, Notation.pyrimidine)

        for snv, count in self.items():
            if snv.ref.is_purine():
                spectrum_pu[snv] = count
                spectrum_pu.weights[str(snv.context)] = self.weights[str(snv.context)]
            elif snv.ref.is_pyrimidine():
                spectrum_py[snv] = count
                spectrum_py.weights[str(snv.context)] = self.weights[str(snv.context)]

        return spectrum_pu, spectrum_py

    def __setitem__(self, key: Variant, value: float) -> None:
        if not isinstance(key, Variant):
            raise TypeError(f'Key must be of type `Variant`: {key}')
        elif not key.is_snv():
            raise ValueError(f'`Variant` must be a single nucleotide variant: {key}')
        elif not isinstance(value, (int, float)):
            raise TypeError(f'Value must be a number: {value}')
        elif self._initialized and key not in self:
            raise KeyError(f'Variant not found in `SnvSpectrum.counts`: {key}')
        super().__setitem__(key, value)

    def __repr__(self) -> str:
        return f'{self.__class__.__qualname__}(k={self.k}, notation={self.notation})'
