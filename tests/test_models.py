import pkgutil
import unittest

from unittest import TestCase
from tempfile import NamedTemporaryFile

from nose.tools import assert_dict_equal
from nose.tools import assert_list_equal
from nose.tools import assert_set_equal
from nose.tools import assert_is_none
from nose.tools import assert_not_equal
from nose.tools import assert_raises
from nose.tools import eq_

import numpy as np

from numpy.testing import assert_almost_equal
from numpy.testing import assert_array_equal

from snv_spectrum import *  # Test import of __all__


class TestSnv(TestCase):
    """Unit tests for ``Snv``"""

    def test_reference_alternate_init(self):
        """Test init with non-single nucleotide reference/alternate"""
        assert_raises(ValueError, Snv, 'Z', 'C')
        assert_raises(ValueError, Snv, 'C', 'Z')

    def test_default_context(self):
        """Test that context defaults to reference"""
        snv = Snv(reference='C', alternate='T')
        eq_(snv.context, 'C')

    def test_invalid_context_set(self):
        """Test setter property for a valid context"""
        snv = Snv(reference='C', alternate='T')

        assert_is_none(setattr(snv, 'context', 'ACA'))

    def test_invalid_context_with_reference(self):
        """Test setter property for an invalid context"""
        snv = Snv(reference='C', alternate='T')

        assert_raises(ValueError, lambda: setattr(snv, 'context', 'ATA'))

    def test_invalid_context_with_length(self):
        """Test setter property for an invalid context"""
        snv = Snv(reference='C', alternate='T')

        assert_raises(ValueError, lambda: setattr(snv, 'context', 'CCCC'))

    def test_invalid_context_with_empty(self):
        """Test setter property for an invalid context"""
        snv = Snv(reference='C', alternate='T')

        assert_raises(ValueError, lambda: setattr(snv, 'context', ''))

    def test_with_purine_reference_no_change(self):
        """Test ``with_purine_reference()`` with a purine reference"""
        snv = Snv(reference='G', alternate='T')

        eq_(snv.with_purine_reference(), snv)

    def test_with_purine_reference_change(self):
        """Test ``with_purine_reference()`` with a pyrimidine reference"""
        snv1 = Snv(reference='C', alternate='T')
        snv2 = Snv(reference='G', alternate='A')

        eq_(snv1.with_purine_reference(), snv2)

    def test_with_purine_reference_change_larger_context(self):
        """Test ``with_pyrimidine_reference()`` with a pyrimidine reference and a
        larger context

        """
        snv1 = Snv(reference='G', alternate='A', context='AGG')
        snv2 = Snv(reference='C', alternate='T', context='CCT')

        eq_(snv1.with_pyrimidine_reference(), snv2)

    def test_with_pyrimidine_reference_no_change(self):
        """Test ``with_pyrimidine_reference()`` with a pyrimidine reference"""
        snv = Snv(reference='C', alternate='A')

        eq_(snv.with_pyrimidine_reference(), snv)

    def test_with_pyrimidine_reference_change(self):
        """Test ``with_pyrimidine_reference()`` with a purine reference"""
        snv1 = Snv(reference='G', alternate='A')
        snv2 = Snv(reference='C', alternate='T')

        eq_(snv1.with_pyrimidine_reference(), snv2)

    def test_with_pyrimidine_reference_change_larger_context(self):
        """Test ``with_purine_reference()`` with a purine reference and a
        larger context

        """
        snv1 = Snv(reference='C', alternate='T', context='CCT')
        snv2 = Snv(reference='G', alternate='A', context='AGG')

        eq_(snv1.with_purine_reference(), snv2)

    def test_copy(self):
        """Test ``copy()`` for equality"""
        snv1 = Snv(reference='C', alternate='G')
        snv2 = snv1.copy()
        eq_(snv1, snv2)

    @unittest.skipIf(
        pkgutil.find_loader('pyfaidx') is None,
        'module pyfaidx could not be found')
    def def_set_context_from_fasta_locus(self):
        """Test ``set_context_from_fasta_locus`` to set context from file"""
        with NamedTemporaryFile('w') as infile:
            fn = infile.name

            infile.write('>seq\nAGTGAGGATGAGA')
            infile.flush()

            snv = Snv('G', 'A')
            snv.set_context_from_fasta_locus(fn, contig='seq', position=1)
            eq_(snv.context, 'AGT')

            snv.set_context_from_fasta_locus(fn, contig='seq', position=6, k=5)
            eq_(snv.context, 'AGGAT')

            # Test that setting the context also emits the context.
            context = snv.set_context_from_fasta_locus(
                fn,
                contig='seq',
                position=1)
            eq_(context, 'AGT')

            # Test argument validation on ``position`` and ``k``.
            snv = Snv('G', 'A')
            assert_raises(
                ValueError,
                lambda: snv.set_context_from_fasta_locus,
                fn, contig='seq', position=1, k=2)
            assert_raises(
                ValueError,
                lambda: snv.set_context_from_fasta_locus,
                fn, contig='seq', position=1, k=-1)
            assert_raises(
                ValueError,
                lambda: snv.set_context_from_fasta_locus,
                fn, contig='seq', position=1.0, k=3)
            assert_raises(
                ValueError,
                lambda: snv.set_context_from_fasta_locus,
                fn, contig='seq', position='200', k=3)

    def test_id(self):
        """Test ``.__hash__()`` for memory inequality on copy"""
        snv1 = Snv(reference='C', alternate='G')
        snv2 = snv1.copy()
        assert_not_equal(id(snv1), id(snv2))

    def test_str(self):
        """Test ``__str__()``"""
        eq_(str(Snv(reference='C', alternate='G')), 'C>G')

    def test_repr(self):
        """Test ``__repr__()``"""
        eq_(repr(Snv(reference='C', alternate='G')),
            'Snv(reference="C", alternate="G", context="C")')


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

    def test_invalid_settitem(self):
        """Test ``__setitem__()`` for not in keys"""
        spectrum = Spectrum(k=1, reference_notation='pyrimidine')

        assert_raises(
            KeyError,
            lambda snv: spectrum[snv],
            Snv(reference='G', alternate='A', context='G'))

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

        assert_list_equal(
            list(spectrum.density().values()),
            [0] * len(spectrum))

        spectrum[Snv('C', 'A')] += 2
        spectrum[Snv('C', 'G')] += 5
        spectrum[Snv('C', 'T')] += 0
        spectrum[Snv('T', 'A')] += 0
        spectrum[Snv('T', 'C')] += 20
        spectrum[Snv('T', 'G')] += 0
        total = 2 + 5 + 20

        assert_almost_equal(spectrum.density()[Snv('C', 'T')], 0)
        assert_almost_equal(spectrum.density()[Snv('T', 'A')], 0)
        assert_almost_equal(spectrum.density()[Snv('T', 'G')], 0)

        assert_almost_equal(spectrum.density()[Snv('C', 'A')], 2 / total)
        assert_almost_equal(spectrum.density()[Snv('C', 'G')], 5 / total)
        assert_almost_equal(spectrum.density()[Snv('T', 'C')], 20 / total)

    def test_density_array_with_uniform_weights(self):
        """Test ``density_as_array`` with default uniform weights set"""
        spectrum = Spectrum(reference_notation='pyrimidine')

        assert_array_equal(
            spectrum.density_as_array(),
            np.zeros(len(spectrum)))

        spectrum[Snv('C', 'A')] += 2
        spectrum[Snv('C', 'G')] += 5
        spectrum[Snv('C', 'T')] += 0
        spectrum[Snv('T', 'A')] += 0
        spectrum[Snv('T', 'C')] += 20
        spectrum[Snv('T', 'G')] += 0
        total = 2 + 5 + 20

        assert_array_equal(
            spectrum.density_as_array(),
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

    def test_iter(self):
        """Test ``__iter__()`` to iterate over substitution items"""
        spectrum = Spectrum()
        assert_list_equal(
            list(spectrum),
            list(spectrum.substitutions.items()))

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
