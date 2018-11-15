import pytest

from nucleic import DNA, Variant, Spectrum

PURINES = ('A', 'G')
PYRIMIDINES = ('C', 'T')


class TestDNA(object):
    """Unit tests for ``nucleic.DNA``."""

    @pytest.mark.parametrize('nt', DNA.complement_map.keys())
    def test_nt_init(self, nt):
        DNA(nt)

    @pytest.mark.parametrize('nt', PURINES)
    def test_nt_is_purine(self, nt):
        assert DNA(nt).is_purine()
        assert not DNA(nt).is_pyrimidine()

    @pytest.mark.parametrize('nt', PYRIMIDINES)
    def test_nt_is_pyrimidine(self, nt):
        assert DNA(nt).is_pyrimidine()
        assert not DNA(nt).is_purine()

    @pytest.mark.parametrize('left', map(DNA, 'ACGT'))
    @pytest.mark.parametrize('right', map(DNA, 'ACGT'))
    def test_nt_to(self, left, right):
        if left == right:
            return
        assert left.to(right) == Variant(left, right)
        assert left.to(str(right)) == Variant(left, right)

    @pytest.mark.parametrize('nt', PURINES + PYRIMIDINES)
    def test_nt__repr__(self, nt):
        assert DNA(nt).__repr__() == f'DNA("{nt}")'


class TestVariant(object):
    """Unit tests for ``nucleic.Variant``."""

    def test_snv_illegal_init(self):
        left, right = DNA('C'), DNA('T')
        with pytest.raises(TypeError):
            Variant(left, 'T')
        with pytest.raises(ValueError):
            Variant(left, left)
        with pytest.raises(TypeError):
            Variant(left, right, locus=2)
        with pytest.raises(TypeError):
            Variant(left, right, context='C')

    @pytest.mark.parametrize('nt', map(DNA, (PURINES + PYRIMIDINES)))
    def test_snv_illegal_init_equal_ref_and_alt(self, nt):
        with pytest.raises(ValueError):
            Variant(nt, nt)

    @pytest.mark.parametrize(
        'snv,color_default,color_stratton',
        [
            [Variant(DNA('A'), DNA('C')), '#D53E4F', '#EDBFC2'],
            [Variant(DNA('T'), DNA('G')), '#D53E4F', '#EDBFC2'],
            [Variant(DNA('C'), DNA('A')), '#3288BD', '#52C3F1'],
        ],
    )
    def test_snv_color_spot_check(self, snv, color_default, color_stratton):
        assert snv.color_default() == color_default
        assert snv.color_stratton() == color_stratton

    @pytest.mark.parametrize('left', map(DNA, (PURINES + PYRIMIDINES)))
    @pytest.mark.parametrize('right', map(DNA, (PURINES + PYRIMIDINES)))
    def test_snv_color_complement_is_same_color(self, left, right):
        if left == right:
            return
        snv1 = left.to(right)
        snv2 = snv1.complement()
        assert snv1.color_default() == snv2.color_default()
        assert snv1.color_stratton() == snv2.color_stratton()

    @pytest.mark.parametrize('left', map(DNA, (PURINES + PYRIMIDINES)))
    @pytest.mark.parametrize('right', map(DNA, (PURINES + PYRIMIDINES)))
    def test_snv_context(self, left, right):
        if left == right:
            return
        snv = left.to(right)
        assert left == snv.context

    def test_snv_context_setter_illegal_values(self):
        snv = DNA('A').to('C')
        with pytest.raises(ValueError):
            snv.context = DNA('AC')
        with pytest.raises(ValueError):
            snv.context = DNA('CCC')

    @pytest.mark.parametrize('left', map(DNA, (PURINES + PYRIMIDINES)))
    @pytest.mark.parametrize('right', map(DNA, (PURINES + PYRIMIDINES)))
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

    @pytest.mark.parametrize('left', map(DNA, (PURINES + PYRIMIDINES)))
    @pytest.mark.parametrize('right', map(DNA, (PURINES + PYRIMIDINES)))
    @pytest.mark.parametrize('lseq', ['', 'A', 'AA', 'AAA'])
    @pytest.mark.parametrize('rseq', ['', 'A', 'AA', 'AAA'])
    def test_snv_within_context_lseq_and_rseq(self, left, right, lseq, rseq):
        if len(lseq) != len(rseq):
            return
        if left == right:
            return
        context = DNA(lseq + str(left) + rseq)
        snv = Variant(left, right).within(context)
        assert snv.lseq() == DNA(lseq)
        assert snv.rseq() == DNA(rseq)

    @pytest.mark.parametrize('snv,result', [(DNA('A').to('C'), 'A→C'), (DNA('T').to('A'), 'T→A')])
    def test_snv_for_snv_label(self, snv, result):
        assert snv.label() == result

    def test_snv_at_locus(self):
        locus = 'chr1:200'
        snv = DNA('A').to('C').at(locus)
        assert snv.locus == locus

    def test_snv_copy(self):
        snv1 = DNA('A').to('C')
        snv2 = snv1.copy()
        assert id(snv1) != id(snv2)

    @pytest.mark.parametrize(
        'left,right', list(zip(DNA.complement_map.keys(), DNA.complement_map.values()))
    )
    def test_snv_complement(self, left, right):
        snv = DNA(left).to(DNA(right))
        assert snv.complement() == DNA(right).to(DNA(left))
        assert snv.context == DNA(right).complement()

        snv.context = DNA(right + left + 'C')
        expected = DNA(right + left + 'C').complement()
        assert snv.complement().context == expected

    @pytest.mark.parametrize(
        'left,right', list(zip(DNA.complement_map.keys(), DNA.complement_map.values()))
    )
    def test_snv_reverse_complement(self, left, right):
        snv = DNA(left).to(DNA(right))
        assert snv.reverse_complement() == DNA(right).to(DNA(left))
        assert snv.context == DNA(right).reverse_complement()

        snv.context = DNA(right + left + 'C')
        expected = DNA(right + left + 'C').reverse_complement()
        assert snv.reverse_complement().context == expected

    @pytest.mark.parametrize('left', map(DNA, PURINES))
    @pytest.mark.parametrize('right', map(DNA, PURINES + PYRIMIDINES))
    def test_snv_with_purine_ref_as_purine_ref(self, left, right):
        if left == right:
            return
        snv = left.to(right)
        assert snv.with_purine_ref() == snv

    @pytest.mark.parametrize('left', map(DNA, PYRIMIDINES))
    @pytest.mark.parametrize('right', map(DNA, PURINES + PYRIMIDINES))
    def test_snv_with_purine_ref_as_pyrimidine_ref(self, left, right):
        if left == right:
            return
        snv = left.to(right)
        assert snv.with_purine_ref() == snv.complement()

    @pytest.mark.parametrize('left', map(DNA, PURINES))
    @pytest.mark.parametrize('right', map(DNA, PURINES + PYRIMIDINES))
    def test_snv_with_pyrimidine_ref_as_purine_ref(self, left, right):
        if left == right:
            return
        snv = left.to(right)
        assert snv.with_pyrimidine_ref() == snv.complement()

    @pytest.mark.parametrize('left', map(DNA, PYRIMIDINES))
    @pytest.mark.parametrize('right', map(DNA, PURINES + PYRIMIDINES))
    def test_snv_with_pyrimidinee_ref_as_pyrimidine_ref(self, left, right):
        if left == right:
            return
        snv = left.to(right)
        assert snv.with_pyrimidine_ref() == snv

    @pytest.mark.parametrize('lseq,rseq,length', [['', '', 1], ['A', 'A', 3], ['AA', 'AA', 5]])
    def test_snv__len__(self, lseq, rseq, length):
        snv = DNA('C').to('A').within(DNA(lseq + 'C' + rseq))
        assert len(snv) == length

    @pytest.mark.parametrize(
        'left,right,context,expected',
        [[DNA('T'), DNA('A'), DNA('T'), '[T→A]'], [DNA('C'), DNA('G'), DNA('ACC'), 'A[C→G]C']],
    )
    def test_snv__str__(self, left, right, context, expected):
        snv = left.to(right).within(context)
        assert str(snv) == expected

    @pytest.mark.parametrize('left,right', [[DNA('T'), DNA('A')], [DNA('C'), DNA('G')]])
    def test_snv__repr__(self, left, right):
        snv = left.to(right)
        assert (
            snv.__repr__() == f'Variant(ref={repr(left)}, alt={repr(right)}, context={repr(left)})'
        )
