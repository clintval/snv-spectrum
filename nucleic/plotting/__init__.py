from collections import namedtuple
from typing import Any, Tuple

import numpy as np
from ordered_set import OrderedSet

from toyplot import Canvas
from toyplot.locator import Explicit as ExplicitLocator

from nucleic import SnvSpectrum

__all__ = ['GridSpec', 'trinucleotide_spectrum']


GridSpec = namedtuple('GridSpec', ('rows', 'columns', 'i', 'rowspan', 'j', 'columnspan'))


def _toyplot_trinucleotide_spectrum(
    spectrum: SnvSpectrum, kind: str = 'count', cmap: str = 'stratton', title: str = ''
) -> Tuple[Canvas, Tuple[Canvas.cartesian, Canvas.cartesian]]:
    """Use the toyplot library to plot a trinucleotide spectrum."""
    N = 96
    num_bins = 6

    cset = OrderedSet(snv.color(cmap) for snv in spectrum)

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
    for i, (data, color) in enumerate(zip(np.split(vector, num_bins), cset)):
        start = i * (N / num_bins)
        x = np.arange(start, start + (N / num_bins))
        ax1.bars(x, data, color=color)

    ax2.bars(np.repeat(1, num_bins), color=cset)

    ax1.x.interactive.coordinates.show = False
    ax2.x.interactive.coordinates.show = False
    ax2.y.interactive.coordinates.show = False

    return canvas, (ax1, ax2)


def _matplotlib_trinucleotide_spectrum(
    spectrum: SnvSpectrum, kind: str = 'count', cmap: str = 'stratton', title: str = ''
) -> Any:
    """Use the matplotlib library to plot a trinucleotide spectrum."""
    import matplotlib.patches as patches
    import matplotlib.pyplot as plt

    N: int = 96
    bar_width: float = 0.6
    patch_padding: float = 0.6
    figsize: Tuple[float, float] = (10, 10 / 3)
    dpi: int = 180

    cset = OrderedSet(snv.color(cmap) for snv in spectrum)
    vector = spectrum.values() if kind == 'count' else spectrum.mass()

    fig, (ax1, ax2) = plt.subplots(
        nrows=2,
        ncols=1,
        dpi=dpi,
        figsize=figsize,
        gridspec_kw={'height_ratios': [28, 1], 'hspace': 0.3, 'wspace': 0.07},
    )

    bars = ax1.bar(x=range(N), height=vector, width=bar_width)
    for i, color in enumerate(cset):
        for bar in bars[16 * i : 16 * i + 16]:
            bar.set_color(color)

        ax2.add_patch(
            patches.Rectangle(
                xy=(16 * i - 1 + patch_padding / 2 + bar_width, 0),
                width=16 - patch_padding,
                height=1,
                color=color,
            )
        )

    for spine in ('top', 'right', 'bottom', 'left'):
        ax2.spines[spine].set_visible(False)
    for spine in ('top', 'right'):
        ax1.spines[spine].set_visible(False)
    xlim = (0 - bar_width, N - 1 + bar_width)

    ax2.get_yaxis().set_visible(False)
    ax2.set_xticks((np.linspace(*xlim, 7) + 16 / 2)[:6])
    ax2.set_xticklabels(OrderedSet(snv.label() for snv in spectrum))

    ax1.set_xticks(range(N))
    ax1.set_xticklabels(
        [str(snv.context) for snv in spectrum],
        rotation=90,
        ha='center',
        y=0.02,
        family='monospace',
    )

    ax1.set_xlim(xlim)
    ax2.set_xlim(xlim)

    return fig, (ax1, ax2)


def trinucleotide_spectrum(
    spectrum: SnvSpectrum,
    kind: str = 'count',
    cmap: str = 'stratton',
    title: str = '',
    plt_library: str = 'toyplot',
) -> Any:
    """Plot the trinucleotide spectrum of mutation.

    Args:
        spectrum: Single nucleotide variants in trinucleotide contexts.
        kind: Whether to plot data as counts or as a probability mass.
        cmap: The color map to use, one of `"stratton"` or `"default"`.
        title: The plot title.
        plt_library: Either plot with *"toyplot"* *"matplotlib"*.

    Returns:
        A canvas, and both Cartesian plotting axes.

    Note:
        The spectrum must be of *pyrimidine* notation.

    """
    N = 96
    if not isinstance(spectrum, SnvSpectrum):
        raise ValueError('`spectrum` is not of class `SnvSpectrum`')
    elif spectrum.notation.name != 'pyrimidine':
        raise ValueError('`spectrum.notation` must have been `Notation.pyrimidine`')
    elif len(spectrum) != N:
        raise ValueError(f'`spectrum` is not of length {N}')
    elif kind not in ('count', 'mass'):
        raise ValueError('`kind` must be "count" or "mass"')
    elif cmap not in ('default', 'stratton'):
        raise ValueError('`cmap` must be "default" or "stratton"')

    if plt_library == 'toyplot':
        return _toyplot_trinucleotide_spectrum(spectrum, kind, cmap, title)
    elif plt_library == 'matplotlib':
        return _matplotlib_trinucleotide_spectrum(spectrum, kind, cmap, title)
    else:
        raise ValueError('`plt_library` must be "toyplot" or "matplotlib".')
