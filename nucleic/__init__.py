import csv
import io
import warnings
from collections import OrderedDict, defaultdict
from enum import Enum
from itertools import permutations
from pathlib import Path
from typing import Any, Dict, Generator, Iterable, List, Mapping, Optional, Set, Tuple, Type, Union
import urllib.request as request

import numpy as np
from pyfaidx import Fasta
from skbio.util import classproperty
from skbio.sequence._nucleotide_mixin import NucleotideMixin
from skbio.sequence import GrammaredSequence

from nucleic.constants import DNA_IUPAC_NONDEGENERATE
from nucleic.sequence import dna_kmers
from nucleic.constants import DEFAULT_SNV_COLOR, STRATTON_SNV_COLOR
from nucleic.util import DictNpArrayMixin, DictMostCommonMixin, DictPrettyReprMixin


__all__ = ['DNA', 'Notation', 'Nt', 'Variant', 'Snv', 'SnvSpectrum', 'fetch_cosmic_signatures']

COSMIC_SIGNATURE_URL = (
    'http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt'
)


class Notation(Enum):
    none: int = 0
    pyrimidine: int = 1
    purine: int = 2


class DNA(GrammaredSequence, NucleotideMixin):
    """Deoxyribonucleic acid composed of the following nucleotide sequences:

    ============ ==================================================== ==========
    String        Residue                                             Class
    ============ ==================================================== ==========
    :class:`A`   `Adenine <https://en.wikipedia.org/wiki/Adenine>`_   Purine
    :class:`C`   `Cytosine <https://en.wikipedia.org/wiki/Cytosine>`_ Pyrimidine
    :class:`G`   `Guanine <https://en.wikipedia.org/wiki/Guanine>`_   Purine
    :class:`T`   `Thymine <https://en.wikipedia.org/wiki/Thymine>`_   Pyrimidine
    ============ ==================================================== ==========

    Examples:
        >>> dna = DNA("A")
        >>> dna.is_purine()
        True
        >>> dna.complement()
        DNA("T")
        >>> DNA("T").to("A")
        Variant(ref=DNA("T"), alt=DNA("A"), context=DNA("T"))

    """

    @classproperty
    def degenerate_map(self) -> Dict[str, Set[str]]:
        return {'.': set(DNA_IUPAC_NONDEGENERATE)}

    @classproperty
    def complement_map(self) -> Dict[str, str]:
        return dict(zip(DNA_IUPAC_NONDEGENERATE, reversed(DNA_IUPAC_NONDEGENERATE)))

    @classproperty
    def definite_chars(self) -> Set[str]:
        return set(DNA_IUPAC_NONDEGENERATE)

    @classproperty
    def default_gap_char(self) -> str:
        return '-'

    @classproperty
    def gap_chars(self) -> Set[str]:
        return set('-')

    def is_purine(self) -> bool:
        """Return if this sequence is a purine."""
        return str(self) in ('A', 'G')

    def is_pyrimidine(self) -> bool:
        """Return if this sequence is a pyrimdine."""
        return str(self) in ('C', 'T')

    def to(self, other: Union[str, 'DNA']) -> 'Variant':
        """Create a variant allele."""
        if isinstance(other, str):
            other = DNA(other)
        return Variant(self, other)

    def __hash__(self) -> int:
        return hash(repr(self))

    def __repr__(self) -> str:
        return f'{self.__class__.__qualname__}("{str(self)}\")'


def Nt(seq: Union[str, DNA]) -> DNA:
    """A single nucleotide of DNA.

    Warning:
        Will be deprecated in ``v0.7.0``. Use :class:`nucleic.DNA` instead.

    """
    warnings.warn('This function will be deprecated in v0.7.0. Use `nucleic.DNA`.')
    return DNA(seq)


class Variant(object):
    def __init__(
        self, ref: DNA, alt: DNA, context: Optional[DNA] = None, locus: Optional[str] = None
    ) -> None:
        if not isinstance(ref, DNA) or not isinstance(alt, DNA):
            raise TypeError('`ref` and `alt` must be of type `DNA`')
        if ref == alt:
            raise ValueError('`ref` and `alt` cannot be the same')
        if locus and not isinstance(locus, str):
            raise TypeError('`locus` must be of type `str`')
        if context and not isinstance(context, DNA):
            raise TypeError('`context` must be of type `DNA`')
        self.ref: DNA = ref
        self.alt: DNA = alt
        self.context = context
        self.locus = locus
        self.counts = Mapping[Variant, float]
        self.weights = Mapping[str, float]

    def color_default(self) -> str:
        """A neutral Hex color representing this class of variant."""
        return DEFAULT_SNV_COLOR[f'{self.ref}→{self.alt}']

    def color_stratton(self) -> str:
        """A Hex color representing this class of variant from Stratton's works."""
        return STRATTON_SNV_COLOR[f'{self.ref}→{self.alt}']

    @property
    def context(self) -> DNA:
        """Return the context of this variant."""
        return self._context

    @context.setter
    def context(self, context: Optional[DNA]) -> None:
        """Verify that the context is of odd length and its center base
        matches the reference.

        Args:
            context: The context of this variant, default to `ref`.

        """
        if context is None:
            self._context = self.ref
            return None

        if not isinstance(context, DNA):
            raise TypeError('`context` must be of type `DNA`')
        if len(context) % 2 != 1:
            raise ValueError(f'Context must be of odd length: {context}')
        if context[int(len(context) / 2)] != self.ref:
            raise ValueError(
                f'Middle of context must equal ref: '
                f'{context[int(len(context) / 2)]} != {self.ref}'
            )

        self._context = context

    @property
    def ref(self) -> DNA:
        """Return the reference of this variant."""
        return self._ref

    @ref.setter
    def ref(self, ref: DNA) -> None:
        """Verify that the reference is of type :class:`nucleic.DNA`."""
        if not isinstance(ref, DNA):
            raise TypeError('`ref` must be of type `DNA`')
        self._ref = ref

    @property
    def alt(self) -> DNA:
        """Return the alternate of this variant."""
        return self._alt

    @alt.setter
    def alt(self, alt: DNA) -> None:
        """Verify that the alternate is of type :class:`nucleic.DNA`."""
        if not isinstance(alt, DNA):
            raise TypeError('`alt` must be of type `DNA`')
        self._alt = alt

    def complement(self) -> 'Variant':
        """Return the complement variant."""
        return self.copy(
            ref=self.ref.complement(), alt=self.alt.complement(), context=self.context.complement()
        )

    def is_transition(self) -> bool:
        """Return if this variant is a transition."""
        return (self.ref.is_pyrimidine() and self.alt.is_pyrimidine()) or (
            self.ref.is_purine() and self.alt.is_purine()
        )

    def is_transversion(self) -> bool:
        """Return if this variant is a transversion."""
        return not self.is_transition()

    def lseq(self) -> DNA:
        """Return the 5′ adjacent sequence to the variant."""
        return DNA(self.context[0 : int((len(self.context) - 1) / 2)])

    def rseq(self) -> DNA:
        """Retrun the 3′ adjacent sequence to the variant."""
        return DNA(self.context[int((len(self.context) - 1) / 2) + 1 :])

    def label(self) -> str:
        """A pretty representation of the variants."""
        return '→'.join(map(str, [self.ref, self.alt]))

    def snv_label(self) -> str:
        """A pretty representation of the variant.

        Warning:
            Will be deprecated in ``v0.7.0``. Use :meth:`Variant.label` instead.

        """
        warnings.warn('This function will be deprecated in v0.7.0. Use `nucleic.DNA`.')
        return self.label()

    def copy(self, **kwargs: Any) -> 'Variant':
        """Make a deep copy of this variant."""
        kwargs = {} if kwargs is None else kwargs
        return Variant(
            ref=kwargs.pop('ref', self.ref),
            alt=kwargs.pop('alt', self.alt),
            context=kwargs.pop('context', self.context),
            locus=kwargs.pop('locus', self.locus),
        )

    def at(self, locus: str) -> 'Variant':
        self.locus = locus
        return self

    def reverse_complement(self) -> 'Variant':
        """Return the reverse complement of this variant."""
        return self.copy(
            ref=self.ref.reverse_complement(),
            alt=self.alt.reverse_complement(),
            context=self.context.reverse_complement(),
        )

    def within(self, context: Union[str, DNA]) -> 'Variant':
        """Set the context of this variant."""
        self.context = DNA(context)
        return self

    def with_purine_ref(self) -> 'Variant':
        """Return this variant with its reference as a purine."""
        return self.reverse_complement() if self.ref.is_pyrimidine() else self

    def with_pyrimidine_ref(self) -> 'Variant':
        """Return this variant with its reference as a pyrimidine."""
        return self.reverse_complement() if self.ref.is_purine() else self

    def set_context_from_fasta(self, infile: Path, contig: str, position: int, k: int = 3) -> str:
        """Set the context by looking up a genomic loci from a FASTA.

        Args:
            infile: FASTA filepath.
            contig: The contig name with containing the locus.
            position: The 0-based contig position that the Variant is centered on.
            k: The length of the context, must be positive and odd.

        Notes:
            - The FASTA file does not need to be indexed.
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
        context = DNA(reference[contig][int(start) : int(end)])
        self.context = context
        return context

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
        return f'{self.lseq()}[{self.ref}→{self.alt}]{self.rseq()   }'


def Snv(ref: DNA, alt: DNA, context: Optional[DNA] = None) -> 'Variant':
    """A single nucleotide variant of type :class:`Variant`.

    Warning:
        Will be deprecated in ``v0.7.0``. Use :class:`nucleic.Variant` instead.

    """
    warnings.warn('This function will be deprecated in v0.7.0. Use `nucleic.Variant`.')
    return Variant(ref, alt, context)


class ContextWeights(DictPrettyReprMixin, DictMostCommonMixin, DictNpArrayMixin, OrderedDict):
    """A dictionary of sequences and their respective weights."""

    def __new__(cls, data: Optional[Any] = None, **kwargs: Mapping[Any, Any]) -> Any:
        return super().__new__(cls, data, **kwargs)


class IndelSpectrum(object):
    """A spectrum of indel variants of various sizes.

    Warning:
        Not implemented.

    """

    def __init__(self) -> None:
        raise NotImplementedError('Class placeholder.')


class SnvSpectrum(DictMostCommonMixin, DictNpArrayMixin, OrderedDict):
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

        for ref, alt in permutations(map(DNA, codes), 2):
            if (
                ref.is_purine()
                and notation.name == 'pyrimidine'
                or ref.is_pyrimidine()
                and notation.name == 'purine'
            ):
                continue
            for context in dna_kmers(k):
                if context[int((k - 1) / 2)] != str(ref):
                    continue

                variant = ref.to(alt)
                self.weights[context] = 1
                self[variant.within(context)] = 0
        self._initialized = True

    @classmethod
    def from_iterable(
        cls, iterable: Iterable, k: int = 1, notation: Notation = Notation.none
    ) -> 'SnvSpectrum':
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

        def norm_count(snv: Variant) -> float:
            if self[snv] != 0 and self.weights[str(snv.context)] == 0:
                raise ValueError('Observations with no weight found: {self[snv]}')
            return float(self[snv] / self.weights[str(snv.context)])

        array = np.array([norm_count(snv) for snv in self.keys()])
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

    @property
    def counts(self) -> 'SnvSpectrum':
        """Return all single nucleotide variants and their counts.

        Warning:
            Will be deprecated in ``v0.7.0``. Use :class:`nucleic.SnvSpectrum` instead.

        """
        warnings.warn('This function will be deprecated in v0.7.0. Use `nucelic.SnvSpectrum`.')
        return self

    def counts_as_array(self) -> np.ndarray:
        """Return all counts as a :class:`numpy.ndarray`.
        
        Warning:
            Will be deprecated in ``v0.7.0``. Use :class:`nucleic.SnvSpectrum.values()` instead.

        """
        warnings.warn(
            'This function will be deprecated in v0.7.0. Use `nucleic.SnvSpectrum.values()`.'
        )
        return self.weights.keys()

    def weights_as_array(self) -> np.ndarray:
        """Return all weights as a :class:`numpy.ndarray`.

        Warning:
            Will be deprecated in ``v0.7.0``. Use :class:`nucleic.SnvSpectrum.weights.values` instead.

        """
        warnings.warn(
            'This function will be deprecated in v0.7.0. Use `nucleic.SnvSpectrum.weights.values`.'
        )
        return self.weights.values()

    def contexts(self) -> np.ndarray:
        """Return all :class:`Variant` key as a :class:`numpy.ndarray`.

        Warning:
            Will be deprecated in ``v0.7.0``. Use :class:`nucleic.SnvSpectrum.weights.keys` instead.

        """
        warnings.warn(
            'This function will be deprecated in v0.7.0. Use `nucleic.SnvSpectrum.weights.keys`.'
        )
        return self.weights.keys()

    def snvs(self) -> np.ndarray:
        """Return all :class:`Variant` key as a :class:`numpy.ndarray`.

        Warning:
            Will be deprecated in ``v0.7.0``. Use :class:`nucleic.SnvSpectrum.keys` instead.

        """
        warnings.warn(
            'This function will be deprecated in v0.7.0. Use `nucleic.SnvSpectrum.keys`.'
        )
        return self.keys()

    def __setitem__(self, key: Variant, value: float) -> None:
        if not isinstance(key, Variant):
            raise TypeError(f'Key must be of type `Variant`: {key}')
        elif not isinstance(value, (int, float)):
            raise TypeError(f'Value must be a number: {value}')
        elif self._initialized and key not in self:
            raise KeyError(f'Variant not found in `SnvSpectrum.counts`: {key}')
        super().__setitem__(key, value)

    def __repr__(self) -> str:
        return f'{self.__class__.__qualname__}(' f'k={self.k}, ' f'notation={self.notation})'


def Spectrum(k: int = 1, notation: Notation = Notation.none) -> SnvSpectrum:
    """A spectrum of single nucleotide variants.

    Warning:
        Will be deprecated in ``v0.7.0``. Use :class:`nucleic.SnvSpectrum` instead.

    """
    warnings.warn('This function will be deprecated in v0.7.0. Use `nucleic.SnvSpectrum`.')
    return SnvSpectrum(k, notation)


def fetch_cosmic_signatures() -> Dict:
    """Fetch the COSMIC published signatures from the following URL.

        - https://cancer.sanger.ac.uk/cosmic

    Returns:
        signatures: The probability masses of the COSMIC signatures.

    """
    from nucleic import DNA, SnvSpectrum, Notation

    all_signatures: defaultdict = defaultdict(
        lambda: SnvSpectrum(k=3, notation=Notation.pyrimidine)
    )

    with request.urlopen(COSMIC_SIGNATURE_URL) as handle:
        reader = csv.reader(io.TextIOWrapper(handle), delimiter='\t')

        # First three columns are subtype, column, and label titles
        _, _, _, *signature_titles = list(filter(None, next(reader)))

        for line in reader:
            subtype, context, _, *points = list(filter(None, line))
            for title, point in zip(signature_titles, map(float, points)):
                # Split the subtype to get reference and alternate
                left, right = subtype.split('>')
                snv = DNA(left).to(right).within(context)
                all_signatures[title][snv] = point

    return dict(all_signatures)
