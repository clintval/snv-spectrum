from nucleic import Nt

__all__ = []

package_scope = locals()
for residue in ('adenine', 'cytosine', 'guanine', 'thymine'):
    code, *_ = residue.upper()
    package_scope[residue] = package_scope[code] = Nt(code)
    __all__.append(code)
