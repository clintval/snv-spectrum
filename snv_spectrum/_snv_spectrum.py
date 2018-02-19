import csv
import io
import urllib.request as request

from collections import Counter, defaultdict
from itertools import permutations, product

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np

from Bio.Seq import complement
from Bio.Seq import reverse_complement

__all__ = [
    'Snv',
    'Spectrum',
    'dna_kmers',
    'get_cosmic_signatures',
    'plot_spectrum',
    'purines',
    'pyrimidines',
    'signature_colors']


signature_colors = [
    '#52C3F1',
    '#231F20',
    '#E62223',
    '#CBC9C8',
    '#97D54C',
    '#EDBFC2']

signature_colors_by_subtype = {
    'A>C': '#EDBFC2', 'A>G': '#97D54C', 'A>T': '#CBC9C8',
    'C>A': '#52C3F1', 'C>G': '#231F20', 'C>T': '#E62223',
    'G>A': '#E62223', 'G>C': '#231F20', 'G>T': '#52C3F1',
    'T>A': '#CBC9C8', 'T>C': '#97D54C', 'T>G': '#EDBFC2'}

longform_subtype = {
    'A>C': 'A:T→C:G', 'A>G': 'A:T→G:C', 'A>T': 'A:T→T:A',
    'C>A': 'C:G→A:T', 'C>G': 'C:G→G:C', 'C>T': 'C:G→T:A',
    'G>A': 'G:C→A:T', 'G>C': 'G:C→C:G', 'G>T': 'G:C→T:A',
    'T>A': 'A:T→T:A', 'T>C': 'A:T→G:C', 'T>G': 'A:T→C:G'}

purines = {'A', 'G'}
pyrimidines = {'C', 'T'}
nucleotides = purines.union(pyrimidines)

COSMIC_SIGNATURE_URL = (
    'http://cancer.sanger.ac.uk/'
    'cancergenome/assets/signatures_probabilities.txt')


class Snv:
    def __init__(self, reference: str, alternate: str, context=None):
        if reference not in nucleotides:
            raise ValueError(f'Reference must be DNA nucleotide: {reference}')
        if alternate not in nucleotides:
            raise ValueError(f'Reference must be DNA nucleotide: {alternate}')

        self.reference = reference
        self.alternate = alternate
        self.context = reference if context is None else context

    @property
    def context(self):
        """Return the context of this Snv."""
        return self._context

    @context.setter
    def context(self, context: str):
        """Verify that the context is of odd length and its center base
        matches the reference.

        context : str
            The context of this Snv, defaults to reference.

        """
        if len(context) % 2 != 1:
            raise ValueError(f'Context must be of odd length: {context}')
        if context[int(len(context) / 2)] != self.reference:
            raise ValueError(
                f'Middle of context must equal reference: '
                f'{context[int(len(context) / 2)]} != {self.reference}')
        self._context = context

    @property
    def with_purine_reference(self):
        """Return this Snv with its reference as a purine."""
        if self.reference not in purines:
            return Snv(
                reference=complement(self.reference),
                alternate=complement(self.alternate),
                context=reverse_complement(self.context))
        else:
            return self.copy()

    @property
    def with_pyrimidine_reference(self):
        """Return this Snv with its reference as a pyrimidine."""
        if self.reference not in pyrimidines:
            return Snv(
                reference=complement(self.reference),
                alternate=complement(self.alternate),
                context=reverse_complement(self.context))
        else:
            return self.copy()

    def copy(self):
        """Make a deep copy of this object."""
        return Snv(self.reference, self.alternate, self.context)

    def __eq__(self, other):
        """Return if the reference, alternate, and context are equal."""
        return (
            self.reference == other.reference and
            self.alternate == other.alternate and
            self.context == other.context)

    def __hash__(self):
        return hash(repr(self))

    def __str__(self):
        """Return the substitution type as commonly seen in text."""
        return '>'.join([self.reference, self.alternate])

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'reference="{self.reference}", '
            f'alternate="{self.alternate}", '
            f'context="{self.context}")')


class Spectrum:
    def __init__(self, k=1, reference_notation=None):
        if not isinstance(k, int) and k % 2 != 1 and k > 0:
            raise ValueError('k must be a positive odd integer')
        if reference_notation not in (None, 'pyrimidine', 'purine'):
            raise ValueError(
                'reference_notation must be "pyrimidine", "purine", or None.')

        self.k = k
        self.reference_notation = reference_notation
        kmers = list(dna_kmers(k))

        self._context_weights = Counter()
        self._substitutions = dict()

        for reference, alternate in sorted(permutations(nucleotides, 2)):
            if (
                reference in purines and reference_notation == 'pyrimidine' or
                reference in pyrimidines and reference_notation == 'purine'
            ):
                continue
            for context in filter(
                lambda kmer: kmer[int((k - 1) / 2)] == reference,
                kmers
            ):
                snv = Snv(reference, alternate, context)
                self._context_weights[context] = 1
                self._substitutions[snv] = 0

    @property
    def as_array(self):
        return np.array(list(self.substitutions.values()))

    @property
    def density(self):
        if sum(self.substitutions.values()) == 0:
            return {snv: 0 for snv in self.substitutions}

        proportions = dict()
        for snv, count in self.substitutions.items():
            proportions[snv] = count / self.context_weights[snv.context]

        density = dict()
        for snv, proportion in proportions.items():
            density[snv] = proportion / sum(proportions.values())

        return density

    @property
    def density_as_array(self):
        return np.array(list(self.density.values()))

    @property
    def contexts(self):
        return [snv.context for snv in self.substitutions]

    @property
    def substitution_types(self):
        return sorted(set(map(str, self.substitutions)))

    @property
    def substitutions(self):
        return self._substitutions

    @property
    def context_weights(self):
        return self._context_weights

    def __iter__(self):
        yield from self.substitutions.items()

    def __len__(self):
        return len(self.substitutions)

    def __getitem__(self, key):
        return self._substitutions[key]

    def __setitem__(self, key, item):
        if key in self._substitutions:
            self._substitutions[key] = item
        else:
            raise KeyError(f'Snv not in spectrum: {key}')

    def __repr__(self):
        if self.reference_notation:
            reference_notation = f'"{self.reference_notation}"'
        else:
            reference_notation = 'None'

        return (
            f'{self.__class__.__name__}('
            f'k={self.k}, '
            f'reference_notation={reference_notation})')


def dna_kmers(k=3):
    """Return the cartesian product of all DNA substrings of length k.

    Parameters
    ----------
    k : int
        Length of of the DNA substring.

    Returns
    -------
    generator of str
        Cartesian product of all DNA substrings of length k.

    Examples
    --------
    >>> list(dna_kmers(1))
    ['A', 'C', 'G', 'T']
    >>> len(list(dna_kmers(3)))
    64

    """
    for parts in product(sorted(nucleotides), repeat=k):
        yield ''.join(parts)


def get_cosmic_signatures():
    """Download the COSMIC published signatures.

    Returns
    -------
    cosmic_signatures : defaultdict
        The probability densities of the COSMIC signatures.

    """
    all_signatures = defaultdict(
        lambda: Spectrum(k=3, reference_notation='pyrimidine'))

    with request.urlopen(COSMIC_SIGNATURE_URL) as handle:
        reader = csv.reader(io.TextIOWrapper(handle), delimiter='\t')
        # First three columns are subtype, column, and label titles
        _, _, _, *signature_titles = list(filter(None, next(reader)))

        for line in reader:
            subtype, context, _, *points = list(filter(None, line))
            for title, point in zip(signature_titles, map(float, points)):
                # Split the subtype to get reference and alternate
                snv = Snv(*subtype.split('>'), context)
                all_signatures[title][snv] = point

    return all_signatures


def plot_spectrum(
    spectrum,
    kind='density',
    bar_width=0.65,
    patch_padding=0.2,
    dpi=180
):
    if not isinstance(spectrum, Spectrum):
        raise ValueError('spectrum is not of class Spectrum')
    elif len(spectrum) != 96:
        raise ValueError('spectrum is not of length 96')
    elif kind not in ('count', 'density'):
        raise ValueError('kind must be "count" or "density"')

    xlim = (0 - bar_width, 96 - 1 + bar_width)

    fig, (ax_main, ax_cbar) = plt.subplots(
        nrows=2,
        ncols=1,
        dpi=dpi,
        figsize=(16, 16 / 3),
        gridspec_kw={
            'height_ratios': [28, 1],
            'hspace': 0.25,
            'wspace': 0.07})

    ax_main.set_axisbelow(True)
    ax_main.yaxis.grid(True, color='0.8', ls='-')

    if kind == 'density':
        vector = spectrum.density_as_array
    elif kind == 'count':
        vector = spectrum.as_array

    bars = ax_main.bar(range(96), vector, width=bar_width)

    for i, color in enumerate(signature_colors):
        for bar in bars[16 * i: 16 * i + 16]:
            bar.set_color(color)

        ax_cbar.add_patch(patches.Rectangle(
            xy=(16 * i - 1 + patch_padding / 2 + bar_width, 0),
            width=16 - patch_padding,
            height=1,
            color=color))

    ax_main.set_xticks(range(96))
    ax_main.set_xticklabels(
        spectrum.contexts,
        rotation=90,
        ha='center',
        family='monospace')

    for spine in ('top', 'right', 'bottom', 'left'):
        ax_cbar.spines[spine].set_visible(False)

    labels = [_.replace('>', '→') for _ in spectrum.substitution_types]

    ax_cbar.get_yaxis().set_visible(False)
    ax_cbar.set_xticks((np.linspace(*xlim, 7) + 16 / 2)[:6])
    ax_cbar.set_xticklabels(labels)

    ax_main.set_xlim(xlim)
    ax_cbar.set_xlim(xlim)

    return fig, (ax_main, ax_cbar)
