from nose.tools import assert_is_instance

from unittest import TestCase

from snv_spectrum import *  # Test import of __all__


class TestGetCosmicSignatures(TestCase):
    """Unit tests for ``get_cosmic_signatures``"""

    def test_successful_download_and_class_types(self):
        """Test ``get_cosmic_signatures()`` keys and values"""
        cosmic_signatures = get_cosmic_signatures()

        for key, value in cosmic_signatures.items():
            assert_is_instance(key, str)
            assert_is_instance(value, Spectrum)
