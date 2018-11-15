from typing import Mapping

DNA_IUPAC_NONDEGENERATE: str = 'ACGT'


STRATTON_SNV_COLOR: Mapping[str, str] = {
    'A→C': '#EDBFC2',
    'A→G': '#97D54C',
    'A→T': '#CBC9C8',
    'C→A': '#52C3F1',
    'C→G': '#231F20',
    'C→T': '#E62223',
    'G→A': '#E62223',
    'G→C': '#231F20',
    'G→T': '#52C3F1',
    'T→A': '#CBC9C8',
    'T→C': '#97D54C',
    'T→G': '#EDBFC2',
}

DEFAULT_SNV_COLOR: Mapping[str, str] = {
    'A→C': '#D53E4F',
    'A→G': '#FC8D59',
    'A→T': '#FEE08B',
    'C→A': '#3288BD',
    'C→G': '#99D594',
    'C→T': '#E6F598',
    'G→A': '#E6F598',
    'G→C': '#99D594',
    'G→T': '#3288BD',
    'T→A': '#FEE08B',
    'T→C': '#FC8D59',
    'T→G': '#D53E4F',
}

LONGFORM_LABEL: Mapping[str, str] = {
    'A→C': 'A:T→C:G',
    'A→G': 'A:T→G:C',
    'A→T': 'A:T→T:A',
    'C→A': 'C:G→A:T',
    'C→G': 'C:G→G:C',
    'C→T': 'C:G→T:A',
    'G→A': 'G:C→A:T',
    'G→C': 'G:C→C:G',
    'G→T': 'G:C→T:A',
    'T→A': 'A:T→T:A',
    'T→C': 'A:T→G:C',
    'T→G': 'A:T→C:G',
}
