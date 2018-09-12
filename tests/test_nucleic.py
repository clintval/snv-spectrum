import pytest

from Bio.Seq import complement, reverse_complement

from nucleic import Nt, Snv, Spectrum
from nucleic import IUPAC_MAPPING
from nucleic.util import PURINES, PYRIMIDINES


class TestNt(object):
    """Unit tests for ``nucleic.Nt``."""

    @pytest.mark.parametrize('nt', IUPAC_MAPPING.keys())
    def test_nt_init(self, nt):
        Nt(nt)

    @pytest.mark.parametrize('nt', ['X', 'a', 1])
    def test_nt_non_iupac_init(self, nt):
        with pytest.raises(ValueError):
            Nt(nt)

    @pytest.mark.parametrize('key,value', IUPAC_MAPPING.items())
    def test_nt_complement(self, key, value):
        assert Nt(key).complement() == Nt(value)

    @pytest.mark.parametrize('nt', PURINES)
    def test_nt_is_purine(self, nt):
        assert Nt(nt).is_purine
        assert not Nt(nt).is_pyrimidine

    @pytest.mark.parametrize('nt', PYRIMIDINES)
    def test_nt_is_pyrimidine(self, nt):
        assert Nt(nt).is_pyrimidine
        assert not Nt(nt).is_purine

    @pytest.mark.parametrize('left', map(Nt, 'ACGT'))
    @pytest.mark.parametrize('right', map(Nt, 'ACGT'))
    def test_nt_to(self, left, right):
        if left == right:
            return
        assert left.to(right) == Snv(left, right)
        assert left.to(str(right)) == Snv(left, right)

    @pytest.mark.parametrize('left,right', [[code] * 2 for code in 'ACGT'])
    def test_nt__eq__(self, left, right):
        assert Nt(left) == Nt(right)

    @pytest.mark.parametrize('nt', PURINES + PYRIMIDINES)
    def test_nt__len__(self, nt):
        assert len(Nt(nt)) == 1

    @pytest.mark.parametrize('nt', PURINES + PYRIMIDINES)
    def test_nt__repr__(self, nt):
        assert Nt(nt).__repr__() == f'Nt("{nt}")'

    @pytest.mark.parametrize('nt', PURINES + PYRIMIDINES)
    def test_nt__str__(self, nt):
        assert nt == str(Nt(nt))


class TestSnv(object):
    """Unit tests for ``nucleic.Snv``."""

    def test_snv_illegal_init(self):
        left, right = Nt('C'), Nt('T')
        with pytest.raises(TypeError):
            Snv(left, 'T')
        with pytest.raises(ValueError):
            Snv(left, left)
        with pytest.raises(TypeError):
            Snv(left, right, locus=2)
        with pytest.raises(TypeError):
            Snv(left, right, context=Nt('C'))

    @pytest.mark.parametrize('nt', map(Nt, (PURINES + PYRIMIDINES)))
    def test_snv_illegal_init_equal_ref_and_alt(self, nt):
        with pytest.raises(ValueError):
            Snv(nt, nt)

    @pytest.mark.parametrize(
        'snv,default_color,stratton_color',
        [
            [Snv(Nt('A'), Nt('C')), '#D53E4F', '#EDBFC2'],
            [Snv(Nt('T'), Nt('G')), '#D53E4F', '#EDBFC2'],
            [Snv(Nt('C'), Nt('A')), '#3288BD', '#52C3F1'],
        ],
    )
    def test_snv_color_spot_check(self, snv, default_color, stratton_color):
        assert snv.default_color == default_color
        assert snv.stratton_color == stratton_color

    @pytest.mark.parametrize('left', map(Nt, (PURINES + PYRIMIDINES)))
    @pytest.mark.parametrize('right', map(Nt, (PURINES + PYRIMIDINES)))
    def test_snv_color_complement_is_same_color(self, left, right):
        if left == right:
            return
        snv1 = left.to(right)
        snv2 = snv1.complement()
        assert snv1.default_color == snv2.default_color
        assert snv1.stratton_color == snv2.stratton_color

    @pytest.mark.parametrize('left', map(Nt, (PURINES + PYRIMIDINES)))
    @pytest.mark.parametrize('right', map(Nt, (PURINES + PYRIMIDINES)))
    def test_snv_context(self, left, right):
        if left == right:
            return
        snv = left.to(right)
        assert str(left) == snv.context

    def test_snv_context_setter_illegal_values(self):
        snv = Nt('A').to('C')
        with pytest.raises(ValueError):
            snv.context = 'AC'
        with pytest.raises(ValueError):
            snv.context = 'CCC'

    @pytest.mark.parametrize('left', map(Nt, (PURINES + PYRIMIDINES)))
    @pytest.mark.parametrize('right', map(Nt, (PURINES + PYRIMIDINES)))
    def test_snv_is_transition_or_transversion(self, left, right):
        if left == right:
            return
        snv = left.to(right)
        if snv.ref.is_purine and snv.alt.is_purine:
            assert snv.is_transition
        if snv.ref.is_pyrimidine and snv.alt.is_pyrimidine:
            assert snv.is_transition
        if snv.ref.is_purine and snv.alt.is_pyrimidine:
            assert snv.is_transversion
        if snv.ref.is_pyrimidine and snv.alt.is_purine:
            assert snv.is_transversion

    @pytest.mark.parametrize('left', map(Nt, (PURINES + PYRIMIDINES)))
    @pytest.mark.parametrize('right', map(Nt, (PURINES + PYRIMIDINES)))
    @pytest.mark.parametrize('lseq', ['', 'A', 'AA', 'AAA'])
    @pytest.mark.parametrize('rseq', ['', 'A', 'AA', 'AAA'])
    def test_snv_within_context_lseq_and_rseq(self, left, right, lseq, rseq):
        if len(lseq) != len(rseq):
            return
        if left == right:
            return
        context = lseq + str(left) + rseq
        snv = Snv(left, right).within(context)
        assert snv.lseq == lseq
        assert snv.rseq == rseq

    @pytest.mark.parametrize('snv,result', [(Nt('A').to('C'), 'A>C'), (Nt('T').to('A'), 'T>A')])
    def test_snv_for_snv_label(self, snv, result):
        assert snv.snv_label == result

    def test_snv_at_locus(self):
        locus = 'chr1:200'
        snv = Nt('A').to('C').at(locus)
        assert snv.locus == locus

    def test_snv_copy(self):
        snv1 = Nt('A').to('C')
        snv2 = snv1.copy()
        assert id(snv1) != id(snv2)

    @pytest.mark.parametrize('left,right', list(zip(IUPAC_MAPPING.keys(), IUPAC_MAPPING.values())))
    def test_snv_complement(self, left, right):
        snv = Nt(left).to(Nt(right))
        assert snv.complement() == Nt(right).to(Nt(left))
        assert snv.context == complement(right)

        snv.context = right + left + 'C'
        expected = complement(right + left + 'C')
        assert snv.complement().context == expected

    @pytest.mark.parametrize('left,right', list(zip(IUPAC_MAPPING.keys(), IUPAC_MAPPING.values())))
    def test_snv_reverse_complement(self, left, right):
        snv = Nt(left).to(Nt(right))
        assert snv.reverse_complement() == Nt(right).to(Nt(left))
        assert snv.context == reverse_complement(right)

        snv.context = right + left + 'C'
        expected = reverse_complement(right + left + 'C')
        assert snv.reverse_complement().context == expected

    @pytest.mark.parametrize('left', map(Nt, PURINES))
    @pytest.mark.parametrize('right', map(Nt, PURINES + PYRIMIDINES))
    def test_snv_with_purine_ref_as_purine_ref(self, left, right):
        if left == right:
            return
        snv = left.to(right)
        assert snv.with_purine_ref() == snv

    @pytest.mark.parametrize('left', map(Nt, PYRIMIDINES))
    @pytest.mark.parametrize('right', map(Nt, PURINES + PYRIMIDINES))
    def test_snv_with_purine_ref_as_pyrimidine_ref(self, left, right):
        if left == right:
            return
        snv = left.to(right)
        assert snv.with_purine_ref() == snv.complement()

    @pytest.mark.parametrize('left', map(Nt, PURINES))
    @pytest.mark.parametrize('right', map(Nt, PURINES + PYRIMIDINES))
    def test_snv_with_pyrimidine_ref_as_purine_ref(self, left, right):
        if left == right:
            return
        snv = left.to(right)
        assert snv.with_pyrimidine_ref() == snv.complement()

    @pytest.mark.parametrize('left', map(Nt, PYRIMIDINES))
    @pytest.mark.parametrize('right', map(Nt, PURINES + PYRIMIDINES))
    def test_snv_with_pyrimidinee_ref_as_pyrimidine_ref(self, left, right):
        if left == right:
            return
        snv = left.to(right)
        assert snv.with_pyrimidine_ref() == snv

    @pytest.mark.parametrize('lseq,rseq,length', [['', '', 1], ['A', 'A', 3], ['AA', 'AA', 5]])
    def test_snv__len__(self, lseq, rseq, length):
        snv = Nt('C').to('A').within(lseq + 'C' + rseq)
        assert len(snv) == length

    @pytest.mark.parametrize(
        'left,right,context', [[Nt('T'), Nt('A'), 'T'], [Nt('C'), Nt('G'), 'ACC']]
    )
    def test_snv__repr__(self, left, right, context):
        snv = left.to(right).within(context)
        assert snv.__repr__() == f'Snv(ref={left}, alt={right}, context="{context}")'

    @pytest.mark.parametrize(
        'left,right,context,expected',
        [[Nt('T'), Nt('A'), 'T', '[T→A]'], [Nt('C'), Nt('G'), 'ACC', 'A[C→G]C']],
    )
    def test_snv__repr__(self, left, right, context, expected):
        snv = left.to(right).within(context)
        assert str(snv) == expected
