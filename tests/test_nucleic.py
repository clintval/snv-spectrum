import pytest

from nucleic import Dna, Variant, SnvSpectrum

from hypothesis import given
from hypothesis.strategies import sampled_from
from hypothesis import strategies as st

purines = PURINES = ('A', 'G')
pyrimidines = PYRIMIDINES = ('C', 'T')

nucleotides = sorted(Dna.complement_map.keys())


class TestDna(object):
    """Unit tests for ``nucleic.Dna``."""

    @given(sampled_from(nucleotides))
    def test_nt_init(self, nt):
        Dna(nt)

    @given(sampled_from(purines))
    def test_nt_is_purine(self, nt):
        assert Dna(nt).is_purine()
        assert not Dna(nt).is_pyrimidine()

    @given(sampled_from(pyrimidines))
    def test_nt_is_pyrimidine(self, nt):
        assert Dna(nt).is_pyrimidine()
        assert not Dna(nt).is_purine()

    @pytest.mark.parametrize('left', map(Dna, 'ACGT'))
    @pytest.mark.parametrize('right', map(Dna, 'ACGT'))
    def test_nt_to(self, left, right):
        if left == right:
            return
        assert left.to(right) == Variant(left, right)
        assert left.to(str(right)) == Variant(left, right)

    @pytest.mark.parametrize('nt', PURINES + PYRIMIDINES)
    def test_nt__repr__(self, nt):
        assert Dna(nt).__repr__() == f'Dna("{nt}")'


class TestVariant(object):
    """Unit tests for ``nucleic.Variant``."""

    @pytest.mark.parametrize('nt', map(Dna, (PURINES + PYRIMIDINES)))
    def test_variant_is_null(self, nt):
        assert Variant(nt, nt).is_null()

    @pytest.mark.parametrize(
        'snv,color_default,color_stratton',
        [
            [Variant(Dna('A'), Dna('C')), '#D53E4F', '#EDBFC2'],
            [Variant(Dna('T'), Dna('G')), '#D53E4F', '#EDBFC2'],
            [Variant(Dna('C'), Dna('A')), '#3288BD', '#52C3F1'],
        ],
    )
    def test_snv_color_spot_check(self, snv, color_default, color_stratton):
        assert snv.color('default') == color_default
        assert snv.color('stratton') == color_stratton

    @pytest.mark.parametrize('left', map(Dna, (PURINES + PYRIMIDINES)))
    @pytest.mark.parametrize('right', map(Dna, (PURINES + PYRIMIDINES)))
    def test_snv_color_complement_is_same_color(self, left, right):
        if left == right:
            return
        snv1 = left.to(right)
        snv2 = snv1.complement()
        assert snv1.color(scheme='default') == snv2.color(scheme='default')
        assert snv1.color(scheme='stratton') == snv2.color(scheme='stratton')

    @pytest.mark.parametrize('left', map(Dna, (PURINES + PYRIMIDINES)))
    @pytest.mark.parametrize('right', map(Dna, (PURINES + PYRIMIDINES)))
    def test_snv_context(self, left, right):
        if left == right:
            return
        snv = left.to(right)
        assert left == snv.context

    def test_snv_context_setter_illegal_values(self):
        snv = Dna('A').to('C')
        with pytest.raises(ValueError):
            snv.context = Dna('AC')
        with pytest.raises(ValueError):
            snv.context = Dna('CCC')

    @pytest.mark.parametrize('left', map(Dna, (PURINES + PYRIMIDINES)))
    @pytest.mark.parametrize('right', map(Dna, (PURINES + PYRIMIDINES)))
    def test_snv_is_transition_or_transversion(self, left, right):
        if left == right:
            return
        snv = left.to(right)
        if snv.ref.is_purine() and snv.alt.is_purine():
            assert snv.is_transition()
        if snv.ref.is_pyrimidine() and snv.alt.is_pyrimidine():
            assert snv.is_transition()
        if snv.ref.is_purine() and snv.alt.is_pyrimidine():
            assert snv.is_transversion()
        if snv.ref.is_pyrimidine() and snv.alt.is_purine():
            assert snv.is_transversion()

    @pytest.mark.parametrize('left', map(Dna, (PURINES + PYRIMIDINES)))
    @pytest.mark.parametrize('right', map(Dna, (PURINES + PYRIMIDINES)))
    @pytest.mark.parametrize('lseq', ['', 'A', 'AA', 'AAA'])
    @pytest.mark.parametrize('rseq', ['', 'A', 'AA', 'AAA'])
    def test_snv_within_context_lseq_and_rseq(self, left, right, lseq, rseq):
        if len(lseq) != len(rseq):
            return
        if left == right:
            return
        context = Dna(lseq + str(left) + rseq)
        snv = Variant(left, right).within(context)
        assert snv.lseq() == Dna(lseq)
        assert snv.rseq() == Dna(rseq)

    @pytest.mark.parametrize('snv,result', [(Dna('A').to('C'), 'A→C'), (Dna('T').to('A'), 'T→A')])
    def test_snv_for_snv_label(self, snv, result):
        assert snv.label() == result

    def test_snv_copy(self):
        snv1 = Dna('A').to('C')
        snv2 = snv1.copy()
        assert id(snv1) != id(snv2)

    @pytest.mark.parametrize(
        'left,right', list(zip(Dna.complement_map.keys(), Dna.complement_map.values()))
    )
    def test_snv_complement(self, left, right):
        snv = Dna(left).to(Dna(right))
        assert snv.complement() == Dna(right).to(Dna(left))
        assert snv.context == Dna(right).complement()

        snv.context = Dna(right + left + 'C')
        expected = Dna(right + left + 'C').complement()
        assert snv.complement().context == expected

    @pytest.mark.parametrize(
        'left,right', list(zip(Dna.complement_map.keys(), Dna.complement_map.values()))
    )
    def test_snv_reverse_complement(self, left, right):
        snv = Dna(left).to(Dna(right))
        assert snv.reverse_complement() == Dna(right).to(Dna(left))
        assert snv.context == Dna(right).reverse_complement()

        snv.context = Dna(right + left + 'C')
        expected = Dna(right + left + 'C').reverse_complement()
        assert snv.reverse_complement().context == expected

    @pytest.mark.parametrize('left', map(Dna, PURINES))
    @pytest.mark.parametrize('right', map(Dna, PURINES + PYRIMIDINES))
    def test_snv_with_purine_ref_as_purine_ref(self, left, right):
        if left == right:
            return
        snv = left.to(right)
        assert snv.with_purine_ref() == snv

    @pytest.mark.parametrize('left', map(Dna, PYRIMIDINES))
    @pytest.mark.parametrize('right', map(Dna, PURINES + PYRIMIDINES))
    def test_snv_with_purine_ref_as_pyrimidine_ref(self, left, right):
        if left == right:
            return
        snv = left.to(right)
        assert snv.with_purine_ref() == snv.complement()

    @pytest.mark.parametrize('left', map(Dna, PURINES))
    @pytest.mark.parametrize('right', map(Dna, PURINES + PYRIMIDINES))
    def test_snv_with_pyrimidine_ref_as_purine_ref(self, left, right):
        if left == right:
            return
        snv = left.to(right)
        assert snv.with_pyrimidine_ref() == snv.complement()

    @pytest.mark.parametrize('left', map(Dna, PYRIMIDINES))
    @pytest.mark.parametrize('right', map(Dna, PURINES + PYRIMIDINES))
    def test_snv_with_pyrimidinee_ref_as_pyrimidine_ref(self, left, right):
        if left == right:
            return
        snv = left.to(right)
        assert snv.with_pyrimidine_ref() == snv

    @pytest.mark.parametrize('lseq,rseq,length', [['', '', 1], ['A', 'A', 3], ['AA', 'AA', 5]])
    def test_snv__len__(self, lseq, rseq, length):
        snv = Dna('C').to('A').within(Dna(lseq + 'C' + rseq))
        assert len(snv) == length

    @pytest.mark.parametrize(
        'left,right,context,expected',
        [[Dna('T'), Dna('A'), Dna('T'), '[T→A]'], [Dna('C'), Dna('G'), Dna('ACC'), 'A[C→G]C']],
    )
    def test_snv__str__(self, left, right, context, expected):
        snv = left.to(right).within(context)
        assert str(snv) == expected

    @pytest.mark.parametrize('left,right', [[Dna('T'), Dna('A')], [Dna('C'), Dna('G')]])
    def test_snv__repr__(self, left, right):
        snv = left.to(right)
        assert (
            snv.__repr__() == f'Variant(ref={repr(left)}, alt={repr(right)}, context={repr(left)})'
        )
