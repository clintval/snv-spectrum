import pytest

from nucleic import SnvSpectrum, fetch_cosmic_signatures
from nucleic.sequence import dna_kmers


class TestUtil(object):
    """Unit tests for ``nucleic.util``."""

    def test_dna_kmers(self):
        assert list(dna_kmers(1)) == ['A', 'C', 'G', 'T']
        assert len(list(dna_kmers(3))) == 64

    @pytest.mark.xfail(reason='No internet')
    def test_fetch_cosmic_signatures(self):
        cosmic_signatures = fetch_cosmic_signatures()

        for key, value in cosmic_signatures.items():
            assert isinstance(key, str)
            assert isinstance(value, SnvSpectrum)
