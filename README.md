# nucleic

[![Testing Status](https://travis-ci.org/clintval/nucleic.svg?branch=master)](https://travis-ci.org/clintval/nucleic)
[![codecov](https://codecov.io/gh/clintval/nucleic/branch/master/graph/badge.svg)](https://codecov.io/gh/clintval/nucleic)
[![Documentation Build Status](https://readthedocs.org/projects/nucleic/badge/?version=latest)](https://nucleic.readthedocs.io/en/latest/?badge=latest)
[![PyPi Release](https://badge.fury.io/py/nucleic.svg)](https://badge.fury.io/py/nucleic)
[![Python Versions](https://img.shields.io/pypi/pyversions/nucleic.svg)](https://pypi.python.org/pypi/nucleic/)
[![MyPy Checked](http://www.mypy-lang.org/static/mypy_badge.svg)](http://mypy-lang.org/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)

Analysis and plotting library for base substitution spectra and signatures.

```bash
❯ pip install nucleic
```

Features:

- Model DNA and variant alleles within their local context using an elegant API
- Combine single nucleotide variants into spectrums of mutagenesis
- Fetch COSMIC signatures of mutation, as well as other published signatures
- SVG plotting functions for displaying single nucleotide variant spectrums

Read the documentation at: [nucleic.readthedocs.io](http://nucleic.readthedocs.io/)

```python
from nucleic.cosmic import fetch_cosmic_signatures
from nucleic.plotting import trinucleotide_spectrum

signatures = fetch_cosmic_signatures()
canvas, (ax1, ax2) = trinucleotide_spectrum(signatures['Signature 24'])
```

![signature-24](https://raw.githubusercontent.com/clintval/nucleic/master/docs/img/signature-24.svg?sanitize=true "Signature 24")
