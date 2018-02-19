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
