from typing import Any, List, Optional

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Shadow

from palettable.cartocolors.qualitative import Prism_10

__all__ = ['plot_seq_record']


FACECOLOR = '0.978'
LABEL_SPACE = 4
TEXT_ZORDER = 20
SHADOW_PROPS = {'alpha': 0.075, 'color': 'k', 'ec': 'none'}

COLOR_MAP = {
    'dna binding': Prism_10.hex_colors[-2],
    'domain': Prism_10.hex_colors[6],
    'helix': Prism_10.hex_colors[1],
    'chain': '0.675',
}

ALPHA_MAP = {'dna binding': 1, 'helix': 1, 'chain': 1}

ZORDER_MAP = {'dna binding': 7, 'helix': 6, 'chain': 3}

HEIGHT_PADDING_MAP = {'dna binding': 0.2, 'helix': 0, 'chain': -0.75}

MOTIF_DEFINITION = {
    'disorder': 'Disordered region (Pfam/IUPred)',
    'low_complexity': 'Low complexity region (Pfam/SEG)',
    'sig_p': 'Signal peptide region (Pfam/Phobius)',
    'coiled_coil': 'Coiled-coil motif (Pfam/ncoils)',
    'transmembrane': 'Transmembrane region (Pfam/Phobius)',
}


def plot_seq_record(seq_record: Any, ax: Optional[Any] = None) -> Any:
    """Render a sequence schematic on a Matplotlib `ax`."""
    ax = ax or plt.gca()

    xticks: List[int] = []
    x_min = float('inf')
    x_max = float('-inf')

    for feature in sorted(
        seq_record.features, key=lambda _: _.location.end.position, reverse=True
    ):
        if feature.type.lower() not in ('dna binding', 'helix', 'chain', 'domain'):
            continue

        start = feature.location.start.position
        end = feature.location.end.position

        x_min, x_max = min(start, x_min), max(end, x_max)

        color = COLOR_MAP.get(feature.type.lower(), '0.7')
        zorder = ZORDER_MAP.get(feature.type.lower(), 6)
        alpha = ALPHA_MAP.get(feature.type.lower(), 1)
        height_padding = HEIGHT_PADDING_MAP.get(feature.type.lower(), 0)

        rectangle = Rectangle(
            xy=(start, -0.5 - (height_padding / 2)),
            width=end - start,
            height=1 + height_padding,
            color=color,
            zorder=zorder,
            alpha=alpha,
        )
        ax.add_patch(rectangle)

        shadow = Shadow(rectangle, ox=0.06, oy=-0.007, props=SHADOW_PROPS)
        shadow.set_zorder(1)

        # Add shadows to all features except chains.
        if feature.type.lower() != 'chain':
            ax.add_patch(shadow)

        # Add tick marks. Starts are guaranteed, ends only appear if they are
        # not within LABEL_SPACE of another start.
        if feature.type.lower() in ('helix', 'chain'):
            if len(xticks) != 0 and xticks[-1] - end < LABEL_SPACE:
                xticks.append(start)
            else:
                xticks.extend((end, start))

        if feature.type.lower() in ('helix', 'dna binding'):
            if feature.type == 'Helix':
                feature_type = r'$\alpha$ helix'
            elif feature.type == 'DNA binding':
                feature_type = 'DNA-binding\ndomain'
            else:
                feature_type = feature.type

            ax.annotate(
                s=feature_type,
                xy=(start + (end - start) / 2, 0),
                color='w',
                ha='center',
                va='center',
                fontsize=9,
                zorder=ZORDER_MAP.get(feature.type.lower(), 6) + 0.5,
            )

    ax.set_xticks(sorted(set(xticks + [77])))
    ax.set_xticklabels([_ + 1 for _ in sorted(set(xticks + [77]))], fontsize=8)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(-0.75, 0.75)

    for spine in ('top', 'left', 'right'):
        ax.spines[spine].set_visible(False)
    ax.yaxis.set_visible(False)

    return ax
