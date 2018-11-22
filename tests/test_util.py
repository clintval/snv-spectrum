import pytest

from nucleic import SnvSpectrum
from nucleic.cosmic import fetch_cosmic_signatures
from nucleic.util import kmers


class TestUtil(object):
    """Unit tests for ``nucleic.util``."""

    def test_kmers(self):
        alphabet = ['A', 'C', 'G', 'T']
        assert list(kmers(1, alphabet)) == alphabet
        assert len(list(kmers(3, alphabet))) == 64

    @pytest.mark.xfail(reason='No internet')
    def test_fetch_cosmic_signatures(self):
        cosmic_signatures = fetch_cosmic_signatures()

        for key, value in cosmic_signatures.items():
            assert isinstance(key, str)
            assert isinstance(value, SnvSpectrum)
