import numpy as np

from nose.tools import assert_dict_equal
from nose.tools import assert_list_equal
from nose.tools import assert_raises
from nose.tools import assert_set_equal
from nose.tools import eq_

from numpy.testing import assert_almost_equal
from numpy.testing import assert_array_equal

from unittest import TestCase

from snv_spectrum import *  # Test import of __all__


class TestSpectrum(TestCase):
    """Unit tests for ``Snv``"""

    def test_error_checking_on_init(self):
        """Test init with invalid argument values"""
        assert_raises(ValueError, Spectrum, 3.2)
        assert_raises(TypeError, Spectrum, None, None)
        assert_raises(ValueError, Spectrum, reference_notation='transversion')
        assert_raises(ValueError, Spectrum, reference_notation=2)

    def test_init_with_valid_arguments(self):
        """Test default init with valid arguments"""
        spectrum = Spectrum()

        eq_(spectrum.k, 1)
        eq_(spectrum.reference_notation, None)

        assert_set_equal(set(spectrum.context_weights.values()), {1})
        assert_set_equal(set(spectrum.substitutions.values()), {0})

    def test_reference_notation_none(self):
        """Test ``reference`` == `None`"""
        spectrum = Spectrum()
        assert_set_equal(set(spectrum.contexts), {'A', 'C', 'G', 'T'})

    def test_substitution_types(self):
        """Test ``substitution_types()```"""
        assert_list_equal(
            Spectrum(reference_notation=None).substitution_types,
            ['A>C', 'A>G', 'A>T', 'C>A', 'C>G', 'C>T',
             'G>A', 'G>C', 'G>T', 'T>A', 'T>C', 'T>G'])

        assert_list_equal(
            Spectrum(reference_notation='pyrimidine').substitution_types,
            ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'])

        assert_list_equal(
            Spectrum(reference_notation='purine').substitution_types,
            ['A>C', 'A>G', 'A>T', 'G>A', 'G>C', 'G>T'])

    def test_context_weights(self):
        """Test ``context_weights`` for init and __set_item__"""
        spectrum = Spectrum()

        assert_dict_equal(
            dict(spectrum.context_weights),
            {'A': 1, 'C': 1, 'G': 1, 'T': 1})

        spectrum.context_weights['C'] += 1

        assert_dict_equal(
            dict(spectrum.context_weights),
            {'A': 1, 'C': 2, 'G': 1, 'T': 1})

    def test_density_with_uniform_weights(self):
        """Test ``density`` with default uniform weights set"""
        spectrum = Spectrum(reference_notation='pyrimidine')

        assert_list_equal(list(spectrum.density.values()), [0] * len(spectrum))

        spectrum[Snv('C', 'A')] += 2
        spectrum[Snv('C', 'G')] += 5
        spectrum[Snv('C', 'T')] += 0
        spectrum[Snv('T', 'A')] += 0
        spectrum[Snv('T', 'C')] += 20
        spectrum[Snv('T', 'G')] += 0
        total = 2 + 5 + 20

        assert_almost_equal(spectrum.density[Snv('C', 'T')], 0)
        assert_almost_equal(spectrum.density[Snv('T', 'A')], 0)
        assert_almost_equal(spectrum.density[Snv('T', 'G')], 0)

        assert_almost_equal(spectrum.density[Snv('C', 'A')], 2 / total)
        assert_almost_equal(spectrum.density[Snv('C', 'G')], 5 / total)
        assert_almost_equal(spectrum.density[Snv('T', 'C')], 20 / total)

    def test_density_array_with_uniform_weights(self):
        """Test ``density_as_array`` with default uniform weights set"""
        spectrum = Spectrum(reference_notation='pyrimidine')

        assert_array_equal(spectrum.density_as_array, np.zeros(len(spectrum)))

        spectrum[Snv('C', 'A')] += 2
        spectrum[Snv('C', 'G')] += 5
        spectrum[Snv('C', 'T')] += 0
        spectrum[Snv('T', 'A')] += 0
        spectrum[Snv('T', 'C')] += 20
        spectrum[Snv('T', 'G')] += 0
        total = 2 + 5 + 20

        assert_array_equal(
            spectrum.density_as_array,
            np.array([2 / total, 5 / total, 0, 0, 20 / total, 0]))

    def test_density_with_non_uniform_weights(self):
        """Test ``density`` with default random weights set"""
        pass  # TODO

    def test_density_array_with_non_uniform_weights(self):
        """Test ``density_as_array`` with random weights set"""
        pass  # TODO

    def test_len(self):
        """Test ``__len__()``"""
        eq_(len(Spectrum()), 12)
        eq_(len(Spectrum(reference_notation='purine')), 6)
        eq_(len(Spectrum(reference_notation='pyrimidine')), 6)

        eq_(len(Spectrum(3)), 96 * 2)
        eq_(len(Spectrum(3, reference_notation='purine')), 96)
        eq_(len(Spectrum(3, reference_notation='pyrimidine')), 96)

    def test_str(self):
        """Test ``__str__()`` to be the same as ``__repr__()``"""
        eq_(str(Spectrum()),
            repr(Spectrum()))
        eq_(str(Spectrum(3)),
            repr(Spectrum(3)))
        eq_(str(Spectrum(3, 'pyrimidine')),
            repr(Spectrum(3, 'pyrimidine')))
        eq_(str(Spectrum(3, 'purine')),
            repr(Spectrum(3, 'purine')))

    def test_repr(self):
        """Test ``__repr__()`` with different parameter combinations"""
        eq_(repr(Spectrum()),
            'Spectrum(k=1, reference_notation=None)')
        eq_(repr(Spectrum(3)),
            'Spectrum(k=3, reference_notation=None)')
        eq_(repr(Spectrum(3, 'pyrimidine')),
            'Spectrum(k=3, reference_notation="pyrimidine")')
        eq_(repr(Spectrum(3, 'purine')),
            'Spectrum(k=3, reference_notation="purine")')
