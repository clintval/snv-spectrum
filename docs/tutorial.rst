Tutorial
========

.. module:: nucleic

Nucleotides
~~~~~~~~~~~

The class :class:`DNA` is an IUPAC valid sequence of non-degenerate DNA nucleotides.
For the purposes of the tutorial we will assume single nucleotide sequences.

.. code-block:: python

    >>> from nucleic import DNA
    >>> DNA("A").is_purine()
    True

Creating Variant Alleles
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    >>> DNA("A").to("C")
    Variant(ref=DNA("A"), alt=DNA("C"), context=DNA("A"))

By default, the context of the variant is assigned to the reference base, although a larger context can be set.
The context must be symmetrical in length about the base substitution otherwise an error will be raised.

.. code-block:: python

    >>> DNA("A").to("C").within("TAG")
    Variant(ref=DNA("A"), alt=DNA("C"), context=DNA("TAG"))

Unless the chemical process for the base substitution is known, it is useful to represent all base substitutions in a canonical form, with either a purine or pyrimidine as the reference base.

.. code-block:: python

    >>> DNA("A").to("C").within("TAG").with_pyrimidine_ref()
    Variant(ref=DNA("T"), alt=DNA("G"), context=DNA("CTA"))

A complete example showing the creation of a notation-normalized :class:`Variant` from strings only:

.. code-block:: python

    >>> ref, alt, context = DNA("A"), DNA("C"), DNA("TAG")
    >>> snv = ref.to(alt).within(context).with_pyrimidine_ref()
    >>> snv.is_transversion()
    True

Each :class:`Variant` has a color associated with it for a uniform color palette.

.. code-block:: python

    >>> snv.color_stratton()
    '#EDBFC2'

Single Nucleotide Variant Spectrums
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A :class:`SnvSpectrum` can be initialized by specifying the size of the local context and the reference notation.

.. code-block:: python

    >>> from nucleic import SnvSpectrum, Notation
    >>> spectrum = SnvSpectrum(k=3, notation=Notation.pyrimidine)
    >>> spectrum
    SnvSpectrum(k=3, notation=Notation.pyrimidine)

Record observations by accessing the :class:`SnvSpectrum` like a Python dictionary.

.. code-block:: python

    spectrum[snv] += 2

*Note*: this is shorthand for ``spectrum.counts[snv] += 2``.

If you have a vector of counts, or probabilities, then you can directly build a :class:`SnvSpectrum` as long as the data is listed in the correct alphabetic order of the :class:`SnvSpectrum` keys.

.. code-block:: python

    >>> vector = [6, 5, 2, 2, 3, 8]
    >>> # SnvSpectrum.from_iterable(vector, k=1, notation=Notation.pyrimidine).counts

Working with Probability
~~~~~~~~~~~~~~~~~~~~~~~~

Many spectra are produced from whole-genome or whole-exome sequencing experiments. Spectra must be normalized to the _kmer_ frequencies in the target study.
Without normalization, no valid spectrum comparison can be made between data generated from different target territories or species.

By default each :class:`nucleic.Variant` is given a weight of 1 and calling :meth:`nucleic.SnvSpectrum.mass_as_array` will simply give the proportion of :class:`nucleic.Variant` counts in the :class:`nucleic.SnvSpectrum`.
After weights are set to the observed *k*-mer counts or frequency of the target territory, calling :meth:`SnvSpectrum.mass` will compute a true normalized probability mass.

All weights can be set with assignment *e.g.*: ``spectrum.context_weights["ACA"] = 23420``.

.. code-block:: python

    >>> # spectrum.mass()

*k*-mer counts can be found with :meth:`skbio.DNA.kmer_frequencies` for large targets.

Fetching COSMIC Signatures
~~~~~~~~~~~~~~~~~~~~~~~~~~

Download the published `COSMIC signatures <http://cancer.sanger.ac.uk/cosmic/signatures>`_ of mutational processes in human cancer:

.. code-block:: python

    >>> from nucleic import fetch_cosmic_signatures
    >>> signatures = fetch_cosmic_signatures()

Plotting Spectrums
~~~~~~~~~~~~~~~~~~

Spectra with ``k=3`` in either ``pyrimidine`` or ``purine`` reference notation can be plotted using a style that was first used in Alexandrov *et. al.*  in 2013 (PMID: `23945592 <https://www.ncbi.nlm.nih.gov/pubmed/23945592>`_). Both :class:`nucleic.Variant` raw counts (``kind="count"``) or their probabilities (``kind="mass"``) can be plotted.

The figure and axes are returned to allow for custom formatting.

.. code-block:: python

    from nucleic.figures import plot_stratton_spectrum

    cosmic_signatures = fetch_cosmic_signatures()

    fig, (ax_main, ax_cbar) = plot_stratton_spectrum(cosmic_signatures["Signature 1"], kind="mass")
    fig, (ax_main, ax_cbar) = plot_stratton_spectrum(cosmic_signatures["Signature 14"], kind="mass")
