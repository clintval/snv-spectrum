import pkgutil
import unittest

from tempfile import NamedTemporaryFile

from nose.tools import assert_is_none
from nose.tools import assert_not_equal
from nose.tools import assert_raises
from nose.tools import eq_

from unittest import TestCase

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

        eq_(snv.with_purine_reference, snv)

    def test_with_purine_reference_change(self):
        """Test ``with_purine_reference()`` with a pyrimidine reference"""
        snv1 = Snv(reference='C', alternate='T')
        snv2 = Snv(reference='G', alternate='A')

        eq_(snv1.with_purine_reference, snv2)

    def test_with_purine_reference_change_larger_context(self):
        """Test ``with_purine_reference()`` with a pyrimidine reference and a
        larger context

        """
        snv1 = Snv(reference='G', alternate='A', context='AGG')
        snv2 = Snv(reference='C', alternate='T', context='CCT')

        eq_(snv1.with_pyrimidine_reference, snv2)

    def test_with_pyrimidine_reference_no_change(self):
        """Test ``with_pyrimidine_reference()`` with a pyrimidine reference"""
        snv = Snv(reference='C', alternate='A')

        eq_(snv.with_pyrimidine_reference, snv)

    def test_with_pyrimidine_reference_change(self):
        """Test ``with_pyrimidine_reference()`` with a purine reference"""
        snv1 = Snv(reference='G', alternate='A')
        snv2 = Snv(reference='C', alternate='T')

        eq_(snv1.with_pyrimidine_reference, snv2)

    def test_with_pyrimidine_reference_change_larger_context(self):
        """Test ``with_pyrimidine_reference()`` with a pyrimidine reference and
        a larger context

        """
        snv1 = Snv(reference='C', alternate='T', context='CCT')
        snv2 = Snv(reference='G', alternate='A', context='AGG')

        eq_(snv1.with_purine_reference, snv2)

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
