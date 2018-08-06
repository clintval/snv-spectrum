from collections import Counter, defaultdict
from enum import Enum
from itertools import permutations, product
from pathlib import Path
from typing import Generator, Iterable, List, Mapping, Optional, Tuple, Union

import numpy as np

from Bio.Seq import complement
from Bio.Seq import reverse_complement
from pyfaidx import Fasta

from snv_spectrum.util import dna_kmers
from snv_spectrum.util import NT_MAPPING, CONTEXT_TYPE
from snv_spectrum.util import PURINES, PYRIMIDINES, SNV_COLOR

__all__ = [
    'Notation',
    'Nt',
    'Snv',
    'Spectrum',
]


class Notation(Enum):
    none: int = 0
    pyrimidine: int = 1
    purine: int = 2


class Nt(object):
    def __init__(
        self,
        nt: Optional[str]='N'
    ) -> None:
        if nt not in NT_MAPPING:
            raise ValueError(f'"{nt}" not a valid IUPAC code')
        self._nt = nt.upper()

    def complement(self) -> 'Nt':
        # Use complement from skbio
        return Nt(complement(self._nt))

    @property
    def is_purine(self) -> bool:
        return self._nt in PURINES

    @property
    def is_pyrimidine(self) -> bool:
        return self._nt in PYRIMIDINES

    def to(self, other: 'Nt') -> 'Snv':
        return Snv(self, other)

    def __eq__(self, other: 'Nt') -> bool:
        return self._nt == other._nt

    def __len__(self) -> int:
        return 1

    def __repr__(self) -> str:
        return (
            f'{self.__class__.__qualname__}('
            f'"{self._nt}\")'
        )

    def __hash__(self) -> int:
        return hash(repr(self))

    def __str__(self) -> str:
        return self._nt


class Snv(object):
    def __init__(
        self,
        ref: Nt,
        alt: Nt,
        context: Optional[CONTEXT_TYPE]=None,
        locus: Optional[str]=None,
    ) -> None:
        if ref == alt:
            raise ValueError('`ref` and `alt` cannot be the same')
        self.ref = ref
        self.alt = alt
        self.context = context
        self.locus = locus

    @property
    def color(self) -> str:
        return SNV_COLOR[f'{self.ref}→{self.alt}']

    @property
    def context(self) -> CONTEXT_TYPE:
        """Return the context of this Snv."""
        return self._context

    @context.setter
    def context(self, context: Optional[CONTEXT_TYPE]) -> None:
        """Verify that the context is of odd length and its center base
        matches the reference.

        context: The context of this Snv, defaults to reference.

        """
        if context is None:
            self._context = CONTEXT_TYPE(self.ref)
            return None

        if len(context) % 2 != 1:
            raise ValueError(f'Context must be of odd length: {context}')
        if context[int(len(context) / 2)] != str(self.ref):
            raise ValueError(
                f'Middle of context must equal ref: '
                f'{context[int(len(context) / 2)]} != {self.ref}')

        self._context = context

    @property
    def is_transition(self):
        return (
            (self.ref.is_pyrimidine and self.alt.is_pyrimidine) or
            (self.ref.is_purine and self.alt.is_purine)
        )

    @property
    def is_transversion(self):
        return (
            (self.ref.is_pyrimidine and self.alt.is_purine) or
            (self.ref.is_purine and self.alt.is_pyrmidine)
        )

    @property
    def lseq(self) -> str:
        return self.context[0: int((len(self.context) - 1) / 2)]

    @property
    def rseq(self) -> str:
        return self.context[int((len(self.context) - 1) / 2) + 1:]

    @property
    def snv_label(self):
        return '>'.join(map(str, [self.ref, self.alt]))

    def copy(self, **kwargs: Mapping) -> 'Snv':
        kwargs = {} if kwargs is None else kwargs
        return Snv(
            ref=kwargs.pop('ref', self.ref),
            alt=kwargs.pop('alt', self.alt),
            context=kwargs.pop('context', self.context),
            locus=kwargs.pop('locus', self.locus)
        )

    def at(self, locus: str) -> 'Snv':
        self.locus = locus
        return self

    def complement(self) -> 'Snv':
        return self.copy(
            ref=self.ref.complement(),
            alt=self.alt.complement(),
            context=complement(self.context)
        )

    def reverse_complement(self) -> 'Snv':
        return self.copy(
            ref=self.ref.complement(),
            alt=self.alt.complement(),
            context=reverse_complement(self.context)
        )

    def within(self, context: Optional[CONTEXT_TYPE]) -> 'Snv':
        self.context = context
        return self

    def with_purine_ref(self) -> 'Snv':
        """Return this Snv with its reference as a purine."""
        return self.reverse_complement() if self.ref.is_pyrimidine else self

    def with_pyrimidine_ref(self) -> 'Snv':
        """Return this Snv with its reference as a pyrimidine."""
        return self.reverse_complement() if self.ref.is_purine else self

    def set_context_from_fasta_locus(
        self,
        infile: Path,
        contig: str,
        position: int,
        k: int=3
    ):
        """Set the context by looking up a genomic loci from a FASTA.

        The length of the context must be odd so the context can be symmetrical
        in length about the the target position.

        Args
            infile: Filepath location to the FASTA file (preferrably indexed).
            contig: The contig or chromosome name for the given locus.
            position: The 0-based position in the chromosome that the context
                will be centered on.
            k: The length of the context, must be positive and odd.

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
        context = str(reference[contig][int(start):int(end)])
        self.context = context
        return context

    def __eq__(self, other: 'Snv') -> bool:
        return (
            self.ref == other.ref and
            self.alt == other.alt and
            str(self.context) == str(self.context)
        )

    def __hash__(self) -> int:
        return hash(f'{self.ref}{self.alt}{self.context}')

    def __len__(self) -> int:
        return len(self.context)

    def __repr__(self) -> str:
        locus = "None" if self.locus is None else f'"{self.locus}"'
        return (
            f'{self.__class__.__qualname__}('
            f'ref={self.ref}, alt={self.alt}, '
            f'context="{self.context}")'
        )

    def __str__(self) -> str:
        return f'{self.lseq}[{self.ref}→{self.alt}]{self.rseq}'


class Spectrum(object):
    def __init__(self, k: int=1, notation: Notation=Notation.none):
        if not isinstance(k, int) and k % 2 != 1 and k > 0:
            raise TypeError('`k` must be a positive odd integer')
        elif not isinstance(notation, Notation):
            raise TypeError('`notation` must be of type Notation')

        self.k = k
        self.notation = notation

        self.counts: Mapping[Snv, float] = {}
        self.weights: Mapping[Snv, float] = {}

        # Reverse the order of a `Spectrum` built with purine notation.
        if notation.value == 2:
            codes = reversed(list(NT_MAPPING.keys()))
        else:
            codes = list(NT_MAPPING.keys())

        for ref, alt in permutations(map(Nt, codes), 2):
            if (
                ref.is_purine and notation.value == 1 or
                ref.is_pyrimidine and notation.value == 2
            ):
                continue
            for context in filter(
                lambda kmer: kmer[int((k - 1) / 2)] == str(ref),
                dna_kmers(k)
            ):
                snv = ref.to(alt).within(context)
                self.counts[snv] = 0
                self.weights[snv.context] = 1

    @property
    def contexts(self) -> List[CONTEXT_TYPE]:
        return [snv.context for snv in self.counts]

    @property
    def snvs(self) -> List[Snv]:
        return [snv for snv in self.counts]

    def mass(self) -> Mapping[Snv, float]:
        def norm_count(snv: Snv):
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
            context = snv.context
            if snv.ref.is_purine:
                spectrum_pu.counts[snv] = count
                spectrum_pu.weights[context] = self.weights[context]
            elif snv.ref.is_pyrimidine:
                spectrum_py.counts[snv] = count
                spectrum_py.weights[context] = self.weights[context]

        return spectrum_pu, spectrum_py

    @classmethod
    def from_iterable(
        cls,
        iterable: Iterable,
        k: int=1,
        notation: Notation=Notation.none,
    ) -> 'Spectrum':
        cls = cls(k=k, notation=notation)
        iterable = list(iterable)  # Cast to list

        if len(cls) != len(iterable):
            raise ValueError('`iterable` not of `len(Spectrum())`')

        for snv, count in zip(cls.counts.keys(), iterable):
            cls[snv] = count
        return cls

    def __iter__(self) -> Generator[Tuple[Snv, Union[float, int]], None, None]:
        yield from self.counts.items()

    def __len__(self) -> int:
        return len(self.counts)

    def __getitem__(self, key) -> Union[float, int]:
        return self.counts[key]

    def __setitem__(self, key, value) -> None:
        if not isinstance(key, Snv):
            raise TypeError(f'Key must be of type `Snv`: {key}')
        elif not isinstance(value, (int, float)):
            raise TypeError(f'Value must be a number: {value}')
        elif key not in self.counts:
            raise KeyError(f'Snv not found in `Spectrum.counts`: {key}')
        self.counts[key] = value

    def __repr__(self) -> str:
        return (
            f'{self.__class__.__qualname__}('
            f'k={self.k}, '
            f'notation={self.notation})')
