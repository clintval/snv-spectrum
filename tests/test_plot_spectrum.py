from unittest import TestCase

import numpy as np

from snv_spectrum import *  # Test import of __all__

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
