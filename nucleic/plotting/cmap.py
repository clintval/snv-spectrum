__all__ = ['SnvCmap']


class SnvCmap(object):
    """Colormaps for SNV variants."""

    #: The default color scheme for SNV variants.
    default = {
        'A→C': '#D53E4F',
        'A→G': '#FC8D59',
        'A→T': '#FEE08B',
        'C→A': '#3288BD',
        'C→G': '#99D594',
        'C→T': '#E6F598',
        'G→A': '#E6F598',
        'G→C': '#99D594',
        'G→T': '#3288BD',
        'T→A': '#FEE08B',
        'T→C': '#FC8D59',
        'T→G': '#D53E4F',
    }

    #: The Stratton et. al. stylized color scheme for SNV variants.
    stratton = {
        'A→C': '#EDBFC2',
        'A→G': '#97D54C',
        'A→T': '#CBC9C8',
        'C→A': '#52C3F1',
        'C→G': '#231F20',
        'C→T': '#E62223',
        'G→A': '#E62223',
        'G→C': '#231F20',
        'G→T': '#52C3F1',
        'T→A': '#CBC9C8',
        'T→C': '#97D54C',
        'T→G': '#EDBFC2',
    }
