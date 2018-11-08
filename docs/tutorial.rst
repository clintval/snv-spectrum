Tutorial
========

.. module:: nucleic

Nucleotides
~~~~~~~~~~~

The class :class:`Nt` represents an IUPAC valid code for a non-degenerate DNA nucleotide.

.. code-block:: python

    >>> from nucleic import Nt
    >>> Nt('A').is_purine()
    True

SNV
~~~

The codes for the for residues adenosine, cytosine, guanine, and thymine are aliased for ease of creating compound objects like a single nucleotide variant (:class:`Snv`):

.. code-block:: python

    >>> from nucleic.extra import A, C, G, T
    >>> A == Nt('A')
    True
    >>> A.to(C)
    Snv(ref=A, alt=C, context="A")

By default, the context of the variant is assigned to the reference base, although a larger context can be set.
The context must be symmetrical in length about the base substitution otherwise an error will be raised.

.. code-block:: python

    >>> A.to(C).within('TAG')
    Snv(ref=A, alt=C, context="TAG")

Unless the chemical process for the base substitution is known, it is useful to represent all base substitutions in a canonical form, with either a purine or pyrimidine as the reference base.

.. code-block:: python

    >>> A.to(C).within('TAG').with_pyrimidine_ref()
    Snv(ref=T, alt=G, context="CTA")

A complete example showing the creation of a notation-normalized :class:`Snv` from strings only:

.. code-block:: python

    >>> ref, alt, context = Nt('A'), Nt('C'), 'TAG'
    >>> snv = ref.to(alt).within(context).with_pyrimidine_ref()
    >>> snv.is_transversion()
    True

Each :class:`Snv` has a color associated with it for a uniform color palette.

.. code-block:: python

    >>> snv.stratton_color
    '#EDBFC2'

An :class:`Snv` can also hold a positional identifier in the :meth:`Snv.locus` property.

.. code-block:: python

    >>> snv = snv.at('chr3:2000')
    >>> snv.locus
    'chr3:2000'

SNV Spectrums
~~~~~~~~~~~~~

A :class:`Spectrum` can be initialized by specifying the size of the local context and the reference notation.

.. code-block:: python

    >>> from nucleic import Spectrum, Notation
    >>> spectrum = Spectrum(k=3, notation=Notation.pyrimidine)
    >>> # spectrum.counts

Record observations by accessing the :class:`Spectrum` like a Python dictionary.

.. code-block:: python

    spectrum[snv] += 2

> *Note*: this is shorthand for `spectrum.counts[snv] += 2`.

If you have a vector of counts, or probabilities, then you can directly build a :class:`Spectrum` as long as the data is listed in the correct alphabetic order of the :class:`Spectrum` keys.

.. code-block:: python

    >>> vector = [6, 5, 2, 2, 3, 8]
    >>> # Spectrum.from_iterable(vector, k=1, notation=Notation.pyrimidine).counts

Working with Probability
~~~~~~~~~~~~~~~~~~~~~~~~

Many spectra are produced from whole-genome or whole-exome sequencing experiments. Spectra must be normalized to the _kmer_ frequencies in the target study.
Without normalization, no valid spectrum comparison can be made between data generated from different target territories or species.

By default each `Snv` is given a weight of 1 and calling :meth:`Spectrum.mass_as_array()` will simply give the proportion of :class:`Snv` counts in the :class:`Spectrum`.
After weights are set to the observed _kmer_ counts or frequency of the target territory, calling `spectrum.mass()` will compute a true normalized probability mass.

All weights can be set with assignment _e.g._: `spectrum.context_weights['ACA'] = 23420`.

.. code-block:: python

    >>> # spectrum.mass()

_Kmer_ counts can be found with [`skbio.DNA.kmer_frequencies`](http://scikit-bio.org/docs/latest/generated/skbio.sequence.DNA.kmer_frequencies.html) for small targets and with [`jellyfish`](http://www.genome.umd.edu/jellyfish.html) for large targets.

Fetching COSMIC Signatures
~~~~~~~~~~~~~~~~~~~~~~~~~~

Download the published [COSMIC signatures](http://cancer.sanger.ac.uk/cosmic/signatures) of mutational processes in human cancer:

.. code-block:: python

    >>> from nucleic.util import fetch_cosmic_signatures
    >>> # fetch_cosmic_signatures()

Plotting Spectrums
~~~~~~~~~~~~~~~~~~

Spectra with `k=3` in either `pyrimidine` or `purine` reference notation can be plotted using a style that was first used in Alexandrov _et. al._  in 2013 (PMID: [23945592](https://www.ncbi.nlm.nih.gov/pubmed/23945592)). Both `Snv` raw counts (`kind="count"`) or their probabilities (`kind="mass"`) can be plotted.

The figure and axes are returned to allow for custom formatting.

.. code-block:: python

    # from nucleic import plot_spectrum

    # cosmic_signatures = fetch_cosmic_signatures()

    # fig, (ax_main, ax_cbar) = plot_spectrum(cosmic_signatures['Signature 1'], kind='mass')
    # fig, (ax_main, ax_cbar) = plot_spectrum(cosmic_signatures['Signature 14'], kind='mass')
