import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np

from typing import Any, Tuple

from nucleic import Spectrum

__all__ = ['plot_spectrum']

signature_colors = ['#52C3F1', '#231F20', '#E62223', '#CBC9C8', '#97D54C', '#EDBFC2']

signature_cmap = {
    'A>C': '#EDBFC2',
    'A>G': '#97D54C',
    'A>T': '#CBC9C8',
    'C>A': '#52C3F1',
    'C>G': '#231F20',
    'C>T': '#E62223',
    'G>A': '#E62223',
    'G>C': '#231F20',
    'G>T': '#52C3F1',
    'T>A': '#CBC9C8',
    'T>C': '#97D54C',
    'T>G': '#EDBFC2',
}


def plot_spectrum(
    spectrum: Spectrum,
    kind: str = 'density',
    bar_width: float = 0.8,
    patch_padding: float = 0.2,
    figsize: Tuple[float, float] = (11, 11 / 3),
    dpi: int = 180,
) -> Any:
    N: int = 96

    if not isinstance(spectrum, Spectrum):
        raise ValueError('spectrum is not of class Spectrum')
    elif len(spectrum) != N:
        raise ValueError(f'spectrum is not of length {N}')
    elif kind not in ('count', 'density'):
        raise ValueError('kind must be "count" or "density"')

    fig, (ax_main, ax_cbar) = plt.subplots(
        nrows=2,
        ncols=1,
        dpi=dpi,
        figsize=figsize,
        gridspec_kw={'height_ratios': [28, 1], 'hspace': 0.25, 'wspace': 0.07},
    )

    fig.set_facecolor('white')
    ax_main.set_axisbelow(True)
    ax_main.yaxis.grid(True, color='0.8', ls='-')

    if kind == 'density':
        vector = spectrum.mass_as_array()
    elif kind == 'count':
        vector = spectrum.counts_as_array()

    bars = ax_main.bar(x=range(N), height=vector, width=bar_width)

    for i, color in enumerate(signature_colors):
        for bar in bars[16 * i : 16 * i + 16]:
            bar.set_color(color)

        rectangle = patches.Rectangle(
            xy=(16 * i - 1 + patch_padding / 2 + bar_width, 0),
            width=16 - patch_padding,
            height=1,
            color=color,
        )
        ax_cbar.add_patch(rectangle)

    ax_main.set_xticks(range(N))
    ax_main.set_xticklabels(
        spectrum.contexts, rotation=90, ha='center', y=0.02, family='monospace'
    )

    for spine in ('top', 'right', 'bottom', 'left'):
        ax_cbar.spines[spine].set_visible(False)

    xlim = (0 - bar_width, N - 1 + bar_width)
    # labels = map(lambda _: _.replace('>', ' to '), spectrum.substitution_types)

    ax_cbar.get_yaxis().set_visible(False)
    ax_cbar.set_xticks((np.linspace(*xlim, 7) + 16 / 2)[:6])
    # ax_cbar.set_xticklabels(labels)

    ax_main.set_xlim(xlim)
    ax_cbar.set_xlim(xlim)

    return fig, (ax_main, ax_cbar)


# def TiTv(
#     labels,
#     spectrums,
#     colors=None,
#     bar_width=0.925,
#     TiTv=None,
#     a_priori=None,
#     cluster=True,
#     optimal_order=False,
#     figsize=(10, 7),
#     dpi=180,
#     gridspec_kw=dict,
#     linkage_kw=dict,
#     dendrogram_kw=dict,
#     **kwargs,
# ):
#     spectrums = spectrums
#     labels = labels
#     nrows = ncols = 1
#     height_ratios, width_ratios = [1], [1]

#     height_ratios = gridspec_kw.pop('height_ratios', height_ratios)
#     width_ratios = gridspec_kw.pop('width_ratios', width_ratios)
#     wspace = gridspec_kw.pop('wspace', 0.07)
#     hspace = gridspec_kw.pop('hspace', 0.06)

#     method = linkage_kw.pop('method', 'ward')
#     metric = linkage_kw.pop('metric', 'cosine')
#     color_threshold = dendrogram_kw.pop('color_threshold', 0)
#     above_threshold_color = dendrogram_kw.pop('above_threshold_color', '0.4')

#     colors = palettable.colorbrewer.diverging.Spectral_6.hex_colors[::-1]

#     if cluster:
#         nrows += 1
#         height_ratios.append(3)

#     if TiTv is not None or a_priori is not None:
#         ncols += 1
#         width_ratios.append(20)

#     fig, axes = plt.subplots(
#         nrows=nrows,
#         ncols=ncols,
#         figsize=figsize,
#         dpi=210,
#         gridspec_kw={
#             'height_ratios': height_ratios,
#             'width_ratios': width_ratios,
#             'wspace': wspace,
#             'hspace': hspace,
#         },
#     )

#     if nrows == 1 and ncols == 1:
#         ax_observed = axes
#     elif nrows == 1 and ncols == 2:
#         ax_expected, ax_observed = axes
#         ax_expected = ax_expected
#     elif nrows == 2 and ncols == 1:
#         ax_dendrogram, ax_observed = axes
#         ax_dendrogram = ax_dendrogram
#     if nrows == 2 and ncols == 2:
#         empty, ax_dendrogram, ax_expected, ax_observed = axes.flatten()
#         ax_dendrogram = ax_dendrogram
#         ax_expected = ax_expected
#         empty.axis('off')

#     ax_observed = ax_observed

#     if cluster is True:
#         try:
#             from fastcluster import linkage
#         except ImportError:
#             from scipy.cluster.hierarchy import linkage

#         from scipy.spatial.distance import pdist
#         from scipy.cluster.hierarchy import dendrogram

#         vectors = [spectrum.counts_as_array() for spectrum in spectrums]

#         Z = linkage(X=vectors, method=method, metric=metric, **linkage_kw)

#         if optimal_order is True:
#             from polo import optimal_leaf_ordering

#             Z = optimal_leaf_ordering(Z, pdist(vectors, metric))

#         dend = dendrogram(
#             Z,
#             above_threshold_color=above_threshold_color,
#             color_threshold=color_threshold,
#             ax=ax_dendrogram,
#             **dendrogram_kw,
#         )

#         ax_dendrogram.set_clip_on(False)

#         # Redfine order of labels.
#         labels = [labels[int(sample)] for sample in dend['ivl']]
#         labels = labels
#         spectrums = [spectrums[int(sample)] for sample in dend['ivl']]
#         spectrums = spectrums

#         # Format the dendrogram ax canvas.
#         ax_dendrogram = axes_off(despine(ax_dendrogram), 'x')
#         ax_dendrogram.spines['right'].set_visible(True)
#         ax_dendrogram.yaxis.tick_right()

#         for tick in ax_dendrogram.yaxis.get_major_ticks():
#             tick.label.set_fontsize(8)

#     substitution_types = spectrums[0].substitution_types

#     bottom = np.repeat(0, len(labels))
#     locations = np.arange(len(labels)) + bar_width / 2

#     for subtype, color in zip(substitution_types, colors):
#         height = [s.density()[Snv(*subtype.split('>'))] for s in spectrums]
#         ax_observed.bar(
#             x=locations,
#             height=height,
#             bottom=bottom,
#             width=bar_width,
#             label=subtype.replace('>', ' to '),
#             color=color,
#         )
#         bottom = list(map(operator.add, height, bottom))

#     ax_observed.set_xticks(locations)
#     ax_observed.set_xticklabels(labels, rotation=90, ha='center', family='monospace')

#     ax_observed.set_xlim(0, len(labels) + bar_width - 1)
#     ax_observed.set_ylim(0, 1)

#     handles, labels = ax_observed.get_legend_handles_labels()
#     ax_observed.legend(handles[::-1], labels[::-1], frameon=False, bbox_to_anchor=(1, 1))

#     if TiTv is not None:
#         rates = [TiTv if s in transitions else 1 for s in substitution_types]
#     elif a_priori is not None:
#         rates = [a_priori[subtype] for subtype in substitution_types]

#     if TiTv is not None or a_priori is not None:
#         bottom = 0
#         expected = [rate / sum(rates) for rate in rates]

#         for tier, color in zip(expected, colors):
#             ax_expected.bar(x=[0], height=tier, bottom=bottom, width=bar_width, color=color)
#             bottom = bottom + tier
#             ax_observed.axhline(bottom, linestyle='--', alpha=0.4)
#             ax_expected.axhline(bottom, linestyle='--', alpha=0.4)

#         ax_expected.set_xticks([0])
#         ax_expected.set_xticklabels(['Expected'])
#         ax_expected.set_xlim(0 - bar_width / 2, 0 + bar_width / 2)
#         ax_expected.set_ylim(0, 1)

#         ticks_off(ax_expected, axis='y')

#         ax_observed.set_ylabel('')
#         ax_observed.set_yticklabels([])

#     return fig, axes
