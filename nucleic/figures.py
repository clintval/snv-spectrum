from collections import namedtuple
from typing import Tuple

import numpy as np

from ordered_set import OrderedSet

from toyplot import Canvas
from toyplot.locator import Explicit as ExplicitLocator

from nucleic import SnvSpectrum

__all__ = ['GridSpec', 'plot_stratton_spectrum']


GridSpec = namedtuple('Grid', ('rows', 'columns', 'i', 'rowspan', 'j', 'columnspan'))


def plot_stratton_spectrum(
    spectrum: SnvSpectrum, kind: str = 'count', title: str = ''
) -> Tuple[Canvas, Tuple[Canvas.cartesian, Canvas.cartesian]]:
    """Plot the trinucleotide spectrum of mutation.

    Args:
        spectrum: single nucleotide variants in trinucleotide contexts.
        kind: whether to plot data as counts or as a probability mass.
        title: the plot title.

    Note:
        The spectrum must be of pyrimidine notation.

    """
    N = 96
    num_bins = 6
    if not isinstance(spectrum, SnvSpectrum):
        raise ValueError('`spectrum` is not of class `SnvSpectrum`')
    elif spectrum.notation.name != 'pyrimidine':
        raise ValueError('`spectrum.notation` must have been `Notation.pyrimidine`')
    elif len(spectrum) != N:
        raise ValueError(f'`spectrum` is not of length {N}')
    elif kind not in ('count', 'mass'):
        raise ValueError('`kind` must be "count" or "mass"')

    cmap = OrderedSet(snv.color_stratton() for snv in spectrum)

    canvas = Canvas(width=900, height=300)
    ax1_grid = GridSpec(rows=1.2, columns=1, i=0, rowspan=1, j=0, columnspan=1)
    ax2_grid = GridSpec(rows=4.8, columns=1, i=3, rowspan=1.775, j=0, columnspan=1)

    ax1 = canvas.cartesian(grid=ax1_grid, label=title)
    ax2 = canvas.cartesian(grid=ax2_grid)

    ax1.x.ticks.show = True
    ax1.x.ticks.labels.angle = 90
    ax1.x.ticks.labels.style.update({'font-family': 'monospace'})

    ax2.x.label.text = 'Substitution Type'
    ax2.y.show = False

    ax1.x.ticks.locator = ExplicitLocator(np.arange(N), [str(snv.context) for snv in spectrum])
    ax2.x.ticks.locator = ExplicitLocator(
        np.arange(num_bins), OrderedSet(snv.label() for snv in spectrum)
    )

    vector = spectrum.values() if kind == 'count' else spectrum.mass()
    for i, (data, color) in enumerate(zip(np.split(vector, num_bins), cmap)):
        start = i * (N / num_bins)
        x = np.arange(start, start + (N / num_bins))
        ax1.bars(x, data, color=color)

    ax2.bars(np.repeat(1, num_bins), color=cmap)

    return canvas, (ax1, ax2)
