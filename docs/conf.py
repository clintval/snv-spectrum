
# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

import os
import sys

sys.path.insert(0, os.path.abspath('../'))

# -- Project information -----------------------------------------------------

project = 'nucleic'
copyright = '2018, TwinStrand Biosciences'
author = 'clintval'

# The short X.Y version
version = '0.6'
# The full version, including alpha/beta/rc tags
release = '0.6.1'

# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.viewcode',
    'sphinx.ext.mathjax',
    'sphinx.ext.intersphinx',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
]

napoleon_include_private_with_doc = True
napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_use_param = False
napoleon_use_ivar = False
napoleon_use_rtype = False

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'biopython': ('https://biopython.readthedocs.io/en/latest/', None),
    'matplotlib': ('https://matplotlib.org/', None),
    'numpy': ('http://docs.scipy.org/doc/numpy/', None),
    'skbio': ('http://scikit-bio.org/docs/latest/', None),
}

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# Plugins for compiling alternative documentation formats
source_suffix = ['.rst', '.md']
source_parsers = {'.md': 'recommonmark.parser.CommonMarkParser'}

# The master toctree document.
master_doc = 'index'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = 'en'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path .
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'default'

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.
# See the documentation for a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"
# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
html_theme_options = {
    'canonical_url': '',
    'analytics_id': '',
    'logo_only': False,
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    # 'vcs_pageview_mode': '',
    # Toc options
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False,
}

# -- Options for HTMLHelp output ---------------------------------------------

# -- Options for LaTeX output ------------------------------------------------

# -- Options for manual page output ------------------------------------------

# -- Options for Texinfo output ----------------------------------------------
