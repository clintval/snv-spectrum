from collections import namedtuple
from typing import Any, Tuple

import palettable
import numpy as np
from ordered_set import OrderedSet

from toyplot import Canvas
from toyplot.locator import Explicit as ExplicitLocator

from nucleic import SnvSpectrum
from .seq_record import *

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

class TvTi:
    def __init__(
        self,
        labels,
        spectrums,
        colors=None,
        bar_width=0.925,
        TiTv=None,
        a_priori=None,
        cluster=True,
        optimal_order=False,
        figsize=(10, 7),
        dpi=180,
        gridspec_kw={},
        linkage_kw={},
        dendrogram_kw={},
        **kwargs
    ):
    
        import matplotlib.pyplot as plt
        from nucleic import Snv
        from nucleic.distance import hierarchy_cluster
        def ticks_off(ax, which='both'):
            assert which in ('x', 'y', 'both'), 'Which must be `x`, `y`, or `both`.'
            axis = ('x', 'y') if which == 'both' else (which, )

            for axe in axis:
                for tic in getattr(ax, f'{axe}axis').get_major_ticks():
                    tic.tick1On = tic.tick2On = False

            return ax
        def axis_off(ax, which='x'):
            """Turn off a specific axis in an ``ax``."""
            getattr(ax, f'get_{which}axis')().set_visible(False)
            return ax
        def despine(ax, top=True, left=True, bottom=True, right=True):
            """Selectively remove spines from an ``ax``."""
            for spine, on in zip(
                ('top', 'left', 'bottom', 'right'), (top, left, bottom, right)
            ):
                ax.spines[spine].set_visible(not on)
            return ax
        self.spectrums = spectrums
        self.labels = labels
        nrows = ncols = 1
        height_ratios, width_ratios = [1], [1]

        height_ratios = gridspec_kw.pop('height_ratios', height_ratios)
        width_ratios = gridspec_kw.pop('width_ratios', width_ratios)
        wspace = gridspec_kw.pop('wspace', 0.07)
        hspace = gridspec_kw.pop('hspace', 0.06)

        method = linkage_kw.pop('method', 'ward')
        metric = linkage_kw.pop('metric', 'cosine')
        color_threshold = dendrogram_kw.pop('color_threshold', 0)

        above_threshold_color = dendrogram_kw.pop(
            'above_threshold_color',
            '0.4')

        colors = palettable.colorbrewer.diverging.Spectral_6.hex_colors[::-1]

        if cluster:
            nrows += 1
            height_ratios.append(3)

        if TiTv is not None or a_priori is not None:
            ncols += 1
            width_ratios.append(20)

        fig, axes = plt.subplots(
            nrows=nrows,
            ncols=ncols,
            figsize=figsize,
            dpi=210,
            gridspec_kw={
                'height_ratios': height_ratios,
                'width_ratios': width_ratios,
                'wspace': wspace,
                'hspace': hspace})

        if nrows == 1 and ncols == 1:
            ax_observed = axes
        elif nrows == 1 and ncols == 2:
            ax_expected, ax_observed = axes
            self.ax_expected = ax_expected
        elif nrows == 2 and ncols == 1:
            ax_dendrogram, ax_observed = axes
            self.ax_dendrogram = ax_dendrogram
        if nrows == 2 and ncols == 2:
            empty, ax_dendrogram, ax_expected, ax_observed = axes.flatten()
            self.ax_dendrogram = ax_dendrogram
            self.ax_expected = ax_expected
            empty.axis('off')

        self.ax_observed = ax_observed

        if cluster is True:
            try:
                from fastcluster import linkage
            except ImportError:
                from scipy.cluster.hierarchy import linkage

            from scipy.spatial.distance import pdist
            from scipy.cluster.hierarchy import dendrogram

            Z = hierarchy_cluster(spectrums)

            self.dend = dendrogram(
                Z,
                above_threshold_color=above_threshold_color,
                color_threshold=color_threshold,
                ax=ax_dendrogram,
                **dendrogram_kw)

            self.ax_dendrogram.set_clip_on(False)

            # Redfine order of labels.
            labels = [labels[int(sample)] for sample in self.dend['ivl']]
            self.labels = labels
            spectrums = [spectrums[int(sample)] for sample in self.dend['ivl']]
            self.spectrums = spectrums

            # Format the dendrogram ax canvas.
            ax_dendrogram = axis_off(despine(ax_dendrogram), 'x')
            ax_dendrogram.spines['right'].set_visible(True)
            ax_dendrogram.yaxis.tick_right()

            for tick in ax_dendrogram.yaxis.get_major_ticks():
                tick.label.set_fontsize(8)

        substitution_types = spectrums[0].keys()

        bottom = np.repeat(0, len(labels))
        locations = np.arange(len(labels)) + bar_width / 2
        import operator
        for i, (subtype, color) in enumerate(zip(substitution_types, colors)):
            height = [s.mass()[i] for s in spectrums]
            ax_observed.bar(
                x=locations,
                height=height,
                bottom=bottom,
                width=bar_width,
                label=subtype.label().replace('→', '>'),
                color=color)
            bottom = list(map(operator.add, height, bottom))

        ax_observed.set_xticks(locations)
        ax_observed.set_xticklabels(
            labels,
            rotation=90,
            ha='center',
            family='monospace')

        ax_observed.set_xlim(0, len(labels) + bar_width - 1)
        ax_observed.set_ylim(0, 1)

        handles, labels = ax_observed.get_legend_handles_labels()
        self.ax_observed_legend = ax_observed.legend(
            handles[::-1],
            labels[::-1],
            frameon=False,
            bbox_to_anchor=(1, 1))

        if TiTv is not None:
            rates = [TiTv if s in transitions else 1
                     for s in substitution_types]
        elif a_priori is not None:
            rates = [a_priori[subtype] for subtype in substitution_types]

        if TiTv is not None or a_priori is not None:
            bottom = 0
            expected = [rate / sum(rates) for rate in rates]

            for tier, color in zip(expected, colors):
                ax_expected.bar(
                    x=[0],
                    height=tier,
                    bottom=bottom,
                    width=bar_width,
                    color=color)
                bottom = bottom + tier
                ax_observed.axhline(bottom, linestyle='--', alpha=0.4)
                ax_expected.axhline(bottom, linestyle='--', alpha=0.4)

            ax_expected.set_xticks([0])
            ax_expected.set_xticklabels(['Expected'])
            ax_expected.set_xlim(0 - bar_width / 2, 0 + bar_width / 2)
            ax_expected.set_ylim(0, 1)

            ticks_off(ax_expected, axis='y')

            ax_observed.set_ylabel('')
            ax_observed.set_yticklabels([])

        self.fig = fig
        self.axes = axes