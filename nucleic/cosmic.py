import csv
import io
import urllib.request as request
from collections import defaultdict
from typing import Dict

from nucleic import Notation, Snv, SnvSpectrum

__all__ = ['fetch_cosmic_signatures']

#: Remote host for the COSMIC cancer signatures flatfile.
COSMIC_SIGNATURE_URL = (
    'http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt'
)


def fetch_cosmic_signatures() -> Dict[str, SnvSpectrum]:
    """Fetch the COSMIC published signatures from the following URL.

        - https://cancer.sanger.ac.uk/cosmic

    Returns:
        signatures: The probability masses of the COSMIC signatures.

    """
    all_signatures: defaultdict = defaultdict(
        lambda: SnvSpectrum(k=3, notation=Notation.pyrimidine)
    )

    with request.urlopen(COSMIC_SIGNATURE_URL) as handle:
        reader = csv.reader(io.TextIOWrapper(handle), delimiter='\t')

        # First three columns are subtype, column, and label titles
        _, _, _, *signature_titles = list(filter(None, next(reader)))

        for line in reader:
            subtype, context, _, *points = list(filter(None, line))
            for title, point in zip(signature_titles, map(float, points)):
                # Split the subtype to get reference and alternate
                left, right = subtype.split('>')
                snv = Snv(left, right, context=context)
                all_signatures[title][snv] = point

    return dict(all_signatures)
