<h1 align="center">snv-spectrum</h2>

<p align="center">Analysis and plotting library for base substitution spectra and signatures</p>

<p align="center">
  <a href="#installation"><strong>Installation</strong></a>
  ·
  <a href="#tutorial"><strong>Tutorial</strong></a>
</p>

<p align="center">
  <a href="https://badge.fury.io/py/snv_spectrum"><img src="https://badge.fury.io/py/snv_spectrum.svg" alt="PyPI version"></img></a>
  <a href="https://travis-ci.org/clintval/snv-spectrum"><img src="https://travis-ci.org/clintval/snv-spectrum.svg?branch=master"></img></a>
  <a href="http://snv-spectrum.readthedocs.io/en/latest/?badge=latest"><img src="https://readthedocs.org/projects/snv-spectrum/badge/?version=latest"></img></a>
  <a href="https://codecov.io/gh/clintval/snv-spectrum"><img src="https://codecov.io/gh/clintval/snv-spectrum/branch/master/graph/badge.svg"></img></a>
  <a href="https://codeclimate.com/github/clintval/snv-spectrum/maintainability"><img src="https://api.codeclimate.com/v1/badges/7f6ce7780716a92c40b8/maintainability"></img></a>
  <a href="http://mypy-lang.org/"><img src="http://www.mypy-lang.org/static/mypy_badge.svg"></img></a>
  <a href="https://github.com/clintval/snv-spectrum/blob/master/LICENSE"><img src="https://img.shields.io/pypi/l/snv-spectrum.svg"></img></a>
</p>

<br>

A base substitution spectrum includes the frequency or probability of all transition and transversions in the DNA grammar. Including the local context of the base substitution aids in the interpretation of biochemical process and disease etiology.

This library aids in the analysis of spectrums with different local context sizes and their visualization.

<br>

<h3 align="center">Installation</h3>

```
❯ pip install snv_spectrum
```

<br>

<h3 align="center">Tutorial</h3>

### Nucleotides

The `Nt` class represents an IUPAC valid code for a non-degenerate DNA nucleotide.

```python
>>> from snv_spectrum import Nt
>>> Nt('A').is_purine
True
```

### SNV

The codes for the for residues adenosine, cytosine, guanine, and thymine are aliased for ease of creating compound objects like a single nucleotide variant (`Snv`):

```python
>>> from snv_spectrum.extra import A, C, G, T
>>> A == Nt('A')
True
>>> A.to(C)
Snv(ref=A, alt=C, context="A")
```

By default, the context of the variant is assigned to the reference base, although a larger context can be set.
The context must be symmetrical in length about the base substitution otherwise an error will be raised.

```python
>>> A.to(C).within('TAG')
Snv(ref=A, alt=C, context="TAG")
```

Unless the chemical process for the base substitution is known, it is useful to represent all base substitutions in a canonical form, with either a purine or pyrimidine as the reference base.

```python
>>> A.to(C).within('TAG').with_pyrimidine_ref()
Snv(ref=T, alt=G, context="CTA")
```

A complete example showing the creation of a notation normalized `Snv` object from strings only:

```python
>>> ref, alt, context = Nt('A'), Nt('C'), 'TAG'
>>> snv = ref.to(alt).within(context).with_pyrimidine_ref()
>>> snv.is_transversion
True
```

Each `Snv` has a color associated with it for a uniform color palette.

```python
>>> snv.color
'#EDBFC2'
```

An `Snv` can also hold a positional identifier in the `locus` property.

```python
>>> snv.at('chr3:2000')
>>> snv.locus
"ch3:2000"
```

### SNV Spectrums

A `Spectrum` can be initialized by specifying the size of the local context and the reference notation.

```python
>>> from snv_spectrum import Spectrum, Notation
>>> spectrum = Spectrum(k=3, notation=Notation.pyrimidine)
>>> spectrum.counts
"""
{Snv(ref=C, alt=A, context="ACA"): 0,
 Snv(ref=C, alt=A, context="ACC"): 0,
 Snv(ref=C, alt=A, context="ACG"): 0,
 Snv(ref=C, alt=A, context="ACT"): 0,
 Snv(ref=C, alt=A, context="CCA"): 0,
 Snv(ref=C, alt=A, context="CCC"): 0,
 ...
"""
```

Record observations by accessing the `Spectrum` like a Python dictionary.

```python
spectrum[snv] += 2
```

> *Note*: this is shorthand for `spectrum.counts[snv] += 2`.

If you have a vector of counts, or probabilities, then you can directly build a `Spectrum` as long as the data is listed in the correct alphabetic order of the `Spectrum` keys.

```python
>>> vector = [6, 5, 2, 2, 3, 8]
>>> Spectrum.from_iterable(vector, k=1, notation=Notation.pyrimidine).counts
"""
{Snv(ref=C, alt=A, context="C"): 6,
 Snv(ref=C, alt=G, context="C"): 5,
 Snv(ref=C, alt=T, context="C"): 2,
 Snv(ref=T, alt=A, context="T"): 2,
 Snv(ref=T, alt=C, context="T"): 3,
 Snv(ref=T, alt=G, context="T"): 8}
"""
```

### Working with Probability

Many spectra are produced from whole-genome or whole-exome sequencing experiments. Spectra must be normalized to the _kmer_ frequencies in the target study. Without normalization, no valid spectrum comparison can be made between data generated from different target territories or species.

By default each `Snv` is given a weight of 1 and calling `spectrum.mass()` will simply give the proportion of `Snv` counts in the `Spectrum`. After weights are set to the observed _kmer_ counts or frequency of the target territory, calling `spectrum.mass()` will compute a true normalized probability mass.

All weights can be set with assignment _e.g._: `spectrum.context_weights['ACA'] = 23420`.

```python
>>> spectrum.mass()
"""
{Snv(ref=C, alt=A, context="ACA"): 0.015677491601343786,
 Snv(ref=C, alt=A, context="ACC"): 0.007838745800671893,
 Snv(ref=C, alt=A, context="ACG"): 0.0011198208286674132,
 Snv(ref=C, alt=A, context="ACT"): 0.006718924972004479,
 Snv(ref=C, alt=A, context="CCA"): 0.010078387458006719,
 Snv(ref=C, alt=A, context="CCC"): 0.008958566629339306,
 ...
"""
```

_Kmer_ counts can be found with [`skbio.DNA.kmer_frequencies`](http://scikit-bio.org/docs/latest/generated/skbio.sequence.DNA.kmer_frequencies.html) for small targets and with [`jellyfish`](http://www.genome.umd.edu/jellyfish.html) for large targets.

### Fetching COSMIC Signatures

Download the published [COSMIC signatures](http://cancer.sanger.ac.uk/cosmic/signatures) of mutational processes in human cancer:

```python
>>> from snv_spectrum.util import fetch_cosmic_signatures
>>> fetch_cosmic_signatures()
"""
{'Signature 1': Spectrum(k=3, notation=Notation.pyrimidine),
 'Signature 2': Spectrum(k=3, notation=Notation.pyrimidine),
 'Signature 3': Spectrum(k=3, notation=Notation.pyrimidine),
 'Signature 4': Spectrum(k=3, notation=Notation.pyrimidine),
 'Signature 5': Spectrum(k=3, notation=Notation.pyrimidine),
...
"""
```

### Plotting Spectrums

Spectra with `k=3` in either `pyrimidine` or `purine` reference notation can be plotted using a style that was first used in Alexandrov _et. al._  in 2013 (PMID: [23945592](https://www.ncbi.nlm.nih.gov/pubmed/23945592)). Both `Snv` raw counts (`kind="count"`) or their probabilities (`kind="mass"`) can be plotted.

The figure and axes are returned to allow for custom formatting.

```python
from snv_spectrum import plot_spectrum

cosmic_signatures = fetch_cosmic_signatures()

fig, (ax_main, ax_cbar) = plot_spectrum(cosmic_signatures['Signature 1'], kind='mass')
fig, (ax_main, ax_cbar) = plot_spectrum(cosmic_signatures['Signature 14'], kind='mass')
```

![Signature 1][signature-1]
![Signature 14][signature-14]

[signature-1]: docs/img/signature-1.png "Signature 1"
[signature-14]: docs/img/signature-14.png "Signature 14"
