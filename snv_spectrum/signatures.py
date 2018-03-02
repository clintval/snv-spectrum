__all__ = [
    'assess_number_signatures',
    'cluster_spectrum',
    'deconstruct_sigs',
    'identify_signatures'
]


def assess_number_signatures(
    signatures,
    n,
    decomposition,
    decomposition_args={},
    iterations=1,
):
    raise NotImplementedError('Function placeholder.')


def cluster_spectrum(
    signatures,
    method='weighted',
    metric='cosine',
    optimal_order=True,
):
    raise NotImplementedError('Function placeholder.')


def deconstruct_sigs(spectrum, signatures, method):
    raise NotImplementedError('Function placeholder.')


def identify_signatures(signatures, n):
    raise NotImplementedError('Function placeholder.')
