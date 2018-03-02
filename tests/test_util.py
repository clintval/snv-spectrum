from nose.tools import assert_is_instance
from nose.tools import assert_list_equal
from nose.tools import eq_

from unittest import TestCase

from snv_spectrum import *  # Test import of __all__


class TestDnaKmers(TestCase):
    """Unit tests for ``dna_kmers``"""

    def test_dna_kmers(self):
        """Test ``dna_kmers()``"""
        assert_list_equal(list(dna_kmers(1)), ['A', 'C', 'G', 'T'])
        eq_(len(list(dna_kmers(3))), 64)


class TestGetCosmicSignatures(TestCase):
    """Unit tests for ``fetch_cosmic_signatures``"""

    def test_successful_download_and_class_types(self):
        """Test ``fetch_cosmic_signatures()`` keys and values"""
        cosmic_signatures = fetch_cosmic_signatures()

        for key, value in cosmic_signatures.items():
            assert_is_instance(key, str)
            assert_is_instance(value, Spectrum)
