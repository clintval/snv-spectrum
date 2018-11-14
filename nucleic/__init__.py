import warnings
from enum import Enum
from itertools import permutations
from pathlib import Path
from typing import Any, Dict, Generator, Iterable, List, Mapping, Optional, Set, Tuple, Union

import numpy as np
from pyfaidx import Fasta
from skbio.util import classproperty
from skbio.sequence._nucleotide_mixin import NucleotideMixin
from skbio.sequence import GrammaredSequence

from nucleic.util import dna_kmers
from nucleic.constants import DNA_IUPAC_NONDEGENERATE
from nucleic.util import DEFAULT_SNV_COLOR, STRATTON_SNV_COLOR

__all__ = ['Dna', 'Notation', 'Nt', 'Snv', 'Spectrum']


class Notation(Enum):
    none: int = 0
    pyrimidine: int = 1
    purine: int = 2


class Dna(GrammaredSequence, NucleotideMixin):
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
        >>> dna = Dna("A")
        >>> dna.is_purine()
        True
        >>> dna.complement()
        Dna("T")
        >>> Dna("T").to("A")
        Snv(ref=Dna("T"), alt=Dna("A"), context=Dna("T"))

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

    def to(self, other: Union[str, 'Dna']) -> 'Snv':
        """Create a variant allele."""
        if isinstance(other, str):
            other = Dna(other)
        return Snv(self, other)

    def __hash__(self) -> int:
        return hash(repr(self))

    def __repr__(self) -> str:
        return f'{self.__class__.__qualname__}("{str(self)}\")'


def Nt(seq: Union[str, Dna]) -> Dna:
    """A single nucleotide of DNA.

    Warning:
        Will be deprecated in ``v0.7.0``. Use :class:`nucleic.Dna` instead.

    """
    warnings.warn('This function will be deprecated in v0.7.0. Use `nucleic.Dna`.')
    return Dna(seq)


class Snv(object):
    def __init__(
        self, ref: Dna, alt: Dna, context: Optional[Dna] = None, locus: Optional[str] = None
    ) -> None:
        if not isinstance(ref, Dna) or not isinstance(alt, Dna):
            raise TypeError('`ref` and `alt` must be of type `Dna`')
        if ref == alt:
            raise ValueError('`ref` and `alt` cannot be the same')
        if locus and not isinstance(locus, str):
            raise TypeError('`locus` must be of type `str`')
        if context and not isinstance(context, Dna):
            raise TypeError('`context` must be of type `Dna`')
        self.ref: Dna = ref
        self.alt: Dna = alt
        self.context = context
        self.locus = locus
        self.counts = Mapping[Snv, float]
        self.weights = Mapping[str, float]

    def color_default(self) -> str:
        """A neutral Hex color representing this class of single nucleotide variant."""
        return DEFAULT_SNV_COLOR[f'{self.ref}→{self.alt}']

    def color_stratton(self) -> str:
        """A Hex color representing this class of single nucleotide variant from Stratton's works."""
        return STRATTON_SNV_COLOR[f'{self.ref}→{self.alt}']

    @property
    def context(self) -> Dna:
        """Return the context of this single nucleotide variant."""
        return self._context

    @context.setter
    def context(self, context: Optional[Dna]) -> None:
        """Verify that the context is of odd length and its center base
        matches the reference.

        Args:
            context: The context of this single nucleotide variant, default to `ref`.

        """
        if context is None:
            self._context = self.ref
            return None

        if not isinstance(context, Dna):
            raise TypeError('`context` must be of type `Dna`')
        if len(context) % 2 != 1:
            raise ValueError(f'Context must be of odd length: {context}')
        if context[int(len(context) / 2)] != self.ref:
            raise ValueError(
                f'Middle of context must equal ref: '
                f'{context[int(len(context) / 2)]} != {self.ref}'
            )

        self._context = context

    @property
    def ref(self) -> Dna:
        """Return the reference of this single nucleotide variant."""
        return self._ref

    @ref.setter
    def ref(self, ref: Dna) -> None:
        """Verify that the reference is of type :class:`nucleic.Dna`."""
        if not isinstance(ref, Dna):
            raise TypeError('`ref` must be of type `Dna`')
        self._ref = ref

    @property
    def alt(self) -> Dna:
        """Return the alternate of this single nucleotide variant."""
        return self._alt

    @alt.setter
    def alt(self, alt: Dna) -> None:
        """Verify that the alternate is of type :class:`nucleic.Dna`."""
        if not isinstance(alt, Dna):
            raise TypeError('`alt` must be of type `Dna`')
        self._alt = alt

    def complement(self) -> 'Snv':
        """Return the complement single nucleotide variant."""
        return self.copy(
            ref=self.ref.complement(), alt=self.alt.complement(), context=self.context.complement()
        )

    def is_transition(self) -> bool:
        """Return if this single nucleotide variant is a transition."""
        return (self.ref.is_pyrimidine() and self.alt.is_pyrimidine()) or (
            self.ref.is_purine() and self.alt.is_purine()
        )

    def is_transversion(self) -> bool:
        """Return if this single nucleotide variant is a transversion."""
        return not self.is_transition()

    def lseq(self) -> Dna:
        """Return the 5′ adjacent sequence to the single nucleotide variant."""
        return Dna(self.context[0 : int((len(self.context) - 1) / 2)])

    def rseq(self) -> Dna:
        """Retrun the 3′ adjacent sequence to the single nucleotide variant."""
        return Dna(self.context[int((len(self.context) - 1) / 2) + 1 :])

    def label(self) -> str:
        """A pretty representation of the single nucleotide variant."""
        return '>'.join(map(str, [self.ref, self.alt]))

    def snv_label(self) -> str:
        """A pretty representation of the single nucleotide variant.

        Warning:
            Will be deprecated in ``v0.7.0``. Use :meth:`Snv.label` instead.

        """
        warnings.warn('This function will be deprecated in v0.7.0. Use `nucleic.Dna`.')
        return self.label()

    def copy(self, **kwargs: Any) -> 'Snv':
        """Make a deep copy of this single nucleotide variant."""
        kwargs = {} if kwargs is None else kwargs
        return Snv(
            ref=kwargs.pop('ref', self.ref),
            alt=kwargs.pop('alt', self.alt),
            context=kwargs.pop('context', self.context),
            locus=kwargs.pop('locus', self.locus),
        )

    def at(self, locus: str) -> 'Snv':
        self.locus = locus
        return self

    def reverse_complement(self) -> 'Snv':
        """Return the reverse complement of this single nucleotide variant."""
        return self.copy(
            ref=self.ref.reverse_complement(),
            alt=self.alt.reverse_complement(),
            context=self.context.reverse_complement(),
        )

    def within(self, context: Union[str, Dna]) -> 'Snv':
        """Set the context of this single nucleotide variant."""
        self.context = Dna(context)
        return self

    def with_purine_ref(self) -> 'Snv':
        """Return this single nucleotide variant with its reference as a purine."""
        return self.reverse_complement() if self.ref.is_pyrimidine() else self

    def with_pyrimidine_ref(self) -> 'Snv':
        """Return this single nucleotide variant with its reference as a pyrimidine."""
        return self.reverse_complement() if self.ref.is_purine() else self

    def set_context_from_fasta(self, infile: Path, contig: str, position: int, k: int = 3) -> str:
        """Set the context by looking up a genomic loci from a FASTA.

        Args:
            infile: FASTA filepath.
            contig: The contig name with containing the locus.
            position: The 0-based contig position that the Snv is centered on.
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
        context = Dna(reference[contig][int(start) : int(end)])
        self.context = context
        return context

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, Snv):
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


class Spectrum(object):
    def __init__(self, k: int = 1, notation: Notation = Notation.none) -> None:
        if not isinstance(k, int) and k % 2 != 1 and k > 0:
            raise TypeError('`k` must be a positive odd integer')
        elif not isinstance(notation, Notation):
            raise TypeError('`notation` must be of type Notation')

        self.k = k
        self.notation = notation

        self.counts: Mapping[Snv, float] = {}
        self.weights: Mapping[str, float] = {}

        # Reverse the order of a `Spectrum` built with purine notation.
        if notation.value == 2:
            codes = list(reversed(DNA_IUPAC_NONDEGENERATE))
        else:
            codes = list(DNA_IUPAC_NONDEGENERATE)

        for ref, alt in permutations(map(Dna, codes), 2):
            if (
                ref.is_purine()
                and notation.value == 1
                or ref.is_pyrimidine()
                and notation.value == 2
            ):
                continue
            for context in filter(lambda kmer: kmer[int((k - 1) / 2)] == str(ref), dna_kmers(k)):
                snv = ref.to(alt).within(context)
                self.counts[snv] = 0
                self.weights[snv.context] = 1

    @property
    def contexts(self) -> List[Dna]:
        return [snv.context for snv in self.counts]

    @property
    def snvs(self) -> List[Snv]:
        return [snv for snv in self.counts]

    def mass(self) -> Mapping[Snv, float]:
        def norm_count(snv: Snv) -> float:
            return self.counts[snv] / self.weights[snv.context]

        proportions = {snv: norm_count(snv) for snv in self.counts}
        total = sum(proportions.values())
        if total == 0:
            return proportions
        else:
            return {snv: proportions[snv] / total for snv in self.counts}

    def counts_as_array(self) -> np.array:
        return np.array(list(self.counts.values()))

    def mass_as_array(self) -> np.array:
        return np.array(list(self.mass().values()))

    def weights_as_array(self) -> np.array:
        return np.array(list(self.weights.values()))

    def split_by_notation(self) -> Tuple['Spectrum', 'Spectrum']:
        if self.notation.value != 0:
            raise TypeError('`Spectrum` notation must be `Notation.none`')

        spectrum_pu = Spectrum(self.k, Notation.purine)
        spectrum_py = Spectrum(self.k, Notation.pyrimidine)

        for snv, count in self.counts.items():
            if snv.ref.is_purine:
                spectrum_pu.counts[snv] = count  # type: ignore
                spectrum_pu.weights[snv.context] = self.weights[snv.context]  # type: ignore
            elif snv.ref.is_pyrimidine():
                spectrum_py.counts[snv] = count  # type: ignore
                spectrum_py.weights[snv.context] = self.weights[snv.context]  # type: ignore

        return spectrum_pu, spectrum_py

    @classmethod
    def from_iterable(
        cls, iterable: Iterable, k: int = 1, notation: Notation = Notation.none
    ) -> 'Spectrum':
        spectrum = cls(k=k, notation=notation)
        iterable = list(iterable)  # Cast to list

        if len(spectrum) != len(iterable):
            raise ValueError('`iterable` not of `len(Spectrum())`')

        for snv, count in zip(spectrum.counts.keys(), iterable):
            spectrum.counts[snv] = count  # type: ignore # pylint: disable=E1101
        return spectrum

    def __iter__(self) -> Generator[Tuple[Snv, Union[float, int]], None, None]:
        yield from self.counts.items()

    def __len__(self) -> int:
        return len(self.counts)

    def __getitem__(self, key: Snv) -> float:
        return self.counts[key]

    def __setitem__(self, key: Snv, value: float) -> None:
        if not isinstance(key, Snv):
            raise TypeError(f'Key must be of type `Snv`: {key}')
        elif not isinstance(value, (int, float)):
            raise TypeError(f'Value must be a number: {value}')
        elif key not in self.counts:
            raise KeyError(f'Snv not found in `Spectrum.counts`: {key}')
        self.counts[key] = value  # type: ignore

    def __repr__(self) -> str:
        return f'{self.__class__.__qualname__}(' f'k={self.k}, ' f'notation={self.notation})'
