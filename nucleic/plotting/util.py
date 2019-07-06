import math
import operator

import matplotlib.ticker as mticker

from matplotlib.patches import Rectangle

from more_itertools import windowed

from nucleic.util import float_range, filter_every


def darken_rgb(rgb, p):
    """Darken an "rgb" value by `p` proportion."""
    assert 0 <= p <= 1, "Proportion must be [0, 1]"
    return [int(x * (1 - p)) for x in rgb]

def ticklabels_to_integer(ax, axis='y'):
    getattr(ax, f'{axis}axis').set_major_locator(
        mticker.MaxNLocator(integer=True))
    return ax


def ticklabels_to_percent(ax, axis='y', precision=1):
    getattr(ax, f'{axis}axis').set_major_formatter(
        mticker.FuncFormatter(
            lambda s, position: f'{{0:0.{precision}f}}%'.format(s * 100)))
    return ax


def ticklabels_to_scientific(ax, axis='y', precision=2):
    getattr(ax, f'{axis}axis').set_major_formatter(
        mticker.FuncFormatter(
            lambda s, position: f'{{:0.{precision}e}}'.format(s)))
    return ax

def ticklabels_to_fancy_science(ax, axis='y', precision=2):
    def SuperScriptinate(number):
        return number.replace('0','⁰').replace('1','¹').replace('2','²').replace('3','³').replace('4','⁴').replace('5','⁵').replace('6','⁶').replace('7','⁷').replace('8','⁸').replace('9','⁹').replace('-','⁻')
    def sci_notation(number, sig_fig=2):
        if number == 0.0:
            return 0
        ret_string = "{0:.{1:d}e}".format(number, sig_fig)
        a,b = ret_string.split("e")
        b = int(b)
        return a + "×10" + SuperScriptinate(str(b))
    getattr(ax, f'{axis}axis').set_major_formatter(
        mticker.FuncFormatter(
            lambda s, position: sci_notation(s, precision)))
    return ax

def ticklabels_to_thousands_sep(ax, axis='y'):
    getattr(ax, f'{axis}axis').set_major_formatter(
        mticker.FuncFormatter(lambda s, position: f'{int(s):,}'))
    return ax


def tile_highlight(
    ax,
    start=0,
    width=1,
    every=1,
    color='0.945',
):

    diff = math.ceil(abs(operator.sub(*ax.get_xlim())))
    ymin, ymax = ax.get_ylim()
    
    patches = []
    for x, _ in filter_every(windowed(float_range(start, start + diff, width), 2, 2), 2):
        patches.append(ax.add_patch(Rectangle(xy=(x, ymin), width=width, height=ymax - ymin, color=color)))
    return patches

def line_every(
    ax,
    start,
    every=1,
    color='black',
    width=1,
    style='-'
):
    diff = math.ceil(abs(operator.sub(*ax.get_xlim())))
    ymin, ymax = ax.get_ylim()

    line_locations = float_range(start, start + diff, jump = every)
    lines = []
    for x in line_locations:
        lines.append(ax.axvline(x, linewidth=width, color=color, linestyle=style))
    return lines