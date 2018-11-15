from typing import Any, Tuple

import matplotlib.patches as patches
import matplotlib.pyplot as plt

import numpy as np

from nucleic import SnvSpectrum

__all__ = ['plot_spectrum']


signature_colors = ['#52C3F1', '#231F20', '#E62223', '#CBC9C8', '#97D54C', '#EDBFC2']

def plot_spectrum(
    spectrum: SnvSpectrum,
    kind: str = 'density',
    bar_width: float = 0.8,
    patch_padding: float = 0.2,
    figsize: Tuple[float, float] = (11, 3.7),
    dpi: int = 180,
) -> Any:
    """Plot the spectrum of mutation."""
    N: int = 96

    if not isinstance(spectrum, SnvSpectrum):
        raise ValueError('`spectrum` is not of class `SnvSpectrum`')
    elif len(spectrum) != N:
        raise ValueError(f'`spectrum` is not of length {N}')
    elif kind not in ('count', 'density'):
        raise ValueError('`kind` must be "count" or "density"')

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
        vector = spectrum.mass()
    elif kind == 'count':
        vector = spectrum.values()

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
    labels = sorted(set([snv.label().replace('>', ' to ') for snv in spectrum.keys()]))

    ax_cbar.get_yaxis().set_visible(False)
    ax_cbar.set_xticks((np.linspace(*xlim, 7) + 16 / 2)[:6])
    ax_cbar.set_xticklabels(labels)

    ax_main.set_xlim(xlim)
    ax_cbar.set_xlim(xlim)

    return fig, (ax_main, ax_cbar)
