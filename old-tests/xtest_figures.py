import numpy as np

from nose.tools import assert_raises

from unittest import TestCase

from snvkit import *  # Test import of __all__
from snvkit.figures import *  # Test import of __all__

# Will need to follow instructions here for easy image comparison
# http://www.davidketcheson.info/2015/01/13/using_matplotlib_image_comparison.html

# from matplotlib.testing.decorators import image_comparison


class TestPlotSpectrum(TestCase):
    """Unit tests for ``plot_spectrum``"""

    def test_plot_spectrum_output(self):
        """Test ``plot_spectrum()`` image output"""
        np.random.seed(1)
        vector = np.random.randint(10, size=96)

        spectrum = Spectrum(k=3, reference_notation='pyrimidine')
        for snv, count in zip(spectrum.substitutions, vector):
            spectrum[snv] = count

        fig, (ax_main, ax_cbar) = plot_spectrum(spectrum, kind='count')
        fig, (ax_main, ax_cbar) = plot_spectrum(spectrum, kind='density')

    def test_init_with_invalid_argument_values(self):
        """Test for ``plot_spectrum`` with invalid argument values"""
        assert_raises(ValueError, plot_spectrum, ['list', 'of', 'values'])
        assert_raises(ValueError, plot_spectrum, Spectrum(k=1))
        assert_raises(
            ValueError,
            plot_spectrum,
            Spectrum(k=3, reference_notation='pyrimidine'),
            kind='boxplot',
        )


class TestTiTv(TestCase):
    """Unit tests for ``TiTv``"""

    def test_plot_spectrum_output(self):
        """Test ``TiTv()`` image output"""
        np.random.seed(1)
        vector = np.random.randint(10, size=96)

        spectrum = Spectrum(k=3, reference_notation='pyrimidine')
        for snv, count in zip(spectrum.substitutions, vector):
            spectrum[snv] = count
