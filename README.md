<h1 align="center">snv-spectrum</h2>

<p align="center">Analysis and plotting library for base substitution spectra and signatures</p>

<p align="center">
  <a href="#installation"><strong>Installation</strong></a>
  ·
  <a href="#tutorial"><strong>Tutorial</strong></a>
  ·
  <a href="#contributing"><strong>Contributing</strong></a>
</p>

<p align="center">
  <a href="https://travis-ci.org/clintval/snv-spectrum">
    <img src="https://travis-ci.org/clintval/snv-spectrum.svg?branch=master"></img>
  </a>

  <a href="https://codecov.io/gh/clintval/snv-spectrum">
    <img src="https://codecov.io/gh/clintval/snv-spectrum/branch/master/graph/badge.svg"></img>
  </a>

  <a href="https://badge.fury.io/py/snv_spectrum">
    <img src="https://badge.fury.io/py/snv_spectrum.svg" alt="PyPI version"></img>
  </a>

  <a href="https://github.com/clintval/snv-spectrum/issues">
    <img src="https://img.shields.io/github/issues/clintval/snv-spectrum.svg"></img>
  </a>

  <a href="https://github.com/clintval/snv-spectrum/blob/master/LICENSE">
    <img src="https://img.shields.io/github/license/clintval/snv-spectrum.svg"></img>
  </a>
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

The base unit of the library is the `Snv` which represents a transition or transversion in a given local context.

The local context must be symmetrical in length about the base substitution.

```python
from snv_spectrum import Snv

snv = Snv(reference='G', alternate='T', context='AGC')
```

Unless the chemical process for the base substitution is specifically known it is useful to represent all base substitutions in a canonical form with either a pyrimidine or purine as the reference base.

```python
>>> snv.with_pyrimidine_reference
Snv(reference="C", alternate="A", context="GCT")
```

You can automatically generate a spectrum of `Snv` by specifying both the size of the local context and the reference notation.

```python
>>> from snv_spectrum import Spectrum
>>> spectrum = Spectrum(k=3, reference_notation='pyrimidine')
>>> list(Spectrum(k=3, reference_notation='pyrimidine'))
"""
[(Snv(reference="C", alternate="A", context="ACA"), 0),
 (Snv(reference="C", alternate="A", context="ACC"), 0),
 (Snv(reference="C", alternate="A", context="ACG"), 0),
 (Snv(reference="C", alternate="A", context="ACT"), 0),
 (Snv(reference="C", alternate="A", context="CCA"), 0),
 (Snv(reference="C", alternate="A", context="CCC"), 0),
 ...
"""
```

Begin to record observations by accessing the `Spectrum` like a Python dictionary.

```python
spectrum[snv] += 2
```

If you already have a vector of counts or probabilities then you can build a `Spectrum` quickly as long as the data is listed in the correct lexicographic order of the chosen reference notation.

```python
vector = [0, 2, 3, 4, ..., 95]
spectrum = Spectrum(k=3, reference_notation='pyrimidine')
for snv, count in zip(spectrum.substitutions, vector):
    spectrum[snv] = count
```

##### Working with Probability

Many spectra are produced from whole-genome or whole-exome sequencing experiments. Spectra must be normalized to the _kmer_ frequencies in the target study. Without normalization, no valid spectrum comparison can be made between data generated from different target territories or species.

By default each `Snv` is given a weight of 1 and calling `spectrum.density` will simply give the proportion of `Snv` counts in the `Spectrum`. After weights are set to the observed _kmer_ counts or frequency of the target territory, calling `spectrum.density` will compute a true normalized probability density.

All weights can be set with assignment _e.g._: `spectrum.weights['ACA'] = 23420`.

```python
>>> spectrum.density
"""
{Snv(reference="C", alternate="A", context="ACA"): 0.015677491601343786,
 Snv(reference="C", alternate="A", context="ACC"): 0.007838745800671893,
 Snv(reference="C", alternate="A", context="ACG"): 0.0011198208286674132,
 Snv(reference="C", alternate="A", context="ACT"): 0.006718924972004479,
 Snv(reference="C", alternate="A", context="CCA"): 0.010078387458006719,
 Snv(reference="C", alternate="A", context="CCC"): 0.008958566629339306,
 ...
"""
```

_Kmer_ counts can be found with [`skbio.DNA.kmer_frequencies`](http://scikit-bio.org/docs/latest/generated/skbio.sequence.DNA.kmer_frequencies.html) for small targets and with [`jellyfish`](http://www.genome.umd.edu/jellyfish.html) for large targets.

##### Plotting

Spectra with `k=3` in either `pyrimidine` or `purine` reference notation can be plotted using a style that was first used in Alexandrov _et. al._  in 2013 (PMID: [23945592](https://www.ncbi.nlm.nih.gov/pubmed/23945592)). Both `Snv` raw counts (`kind="count"`) or their probabilities (`kind="density"`) can be plotted.

The figure and axes are returned to allow for custom formatting.

```python
import numpy

from snv_spectrum import plot_spectrum

spectrum = Spectrum(k=3, reference_notation='pyrimidine')

for snv, count in zip(spectrum.substitutions, range(96)):
     spectrum[snv] = numpy.random.randint(20)

fig, (ax_main, ax_cbar) = plot_spectrum(spectrum, kind='density')
```

![Demo Plot Random][demo-plot-random]

##### Fetching COSMIC Signatures

You can also directly download and view the published [COSMIC signatures](http://cancer.sanger.ac.uk/cosmic/signatures) of mutational processes in human cancer.

```python
from snv_spectrum import get_cosmic_signatures

cosmic_signatures = get_cosmic_signatures()

fig, axes = plot_spectrum(cosmic_signatures['Signature 1'], kind='density')
fig, axes = plot_spectrum(cosmic_signatures['Signature 14'], kind='density')
```
![Signature 1][signature-1]
![Signature 14][signature-14]

<br>

<h3 align="center">Contributing</h3>

Pull requests, feature requests, and issues welcome!

To make a development install:

```bash
❯ git clone git@github.com:clintval/snv-spectrum.git
❯ pip install -e 'snv-spectrum[fancytest]'
```

To run the tests:

```bash
❯ ./snv-spectrum/tests/run-tests
```

[demo-plot-random]: docs/img/demo-plot-random.png "Demo Plot Random"
[signature-1]: docs/img/signature-1.png "Signature 1"
[signature-14]: docs/img/signature-14.png "Signature 14"
