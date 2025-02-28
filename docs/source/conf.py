# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'SqueezeMeta'
copyright = '2018, Javier Tamames & Fernando Puente-Sánchez'
author = 'Javier Tamames & Fernando Puente-Sánchez'

release = '1.7.0.alpha8'
version = '1.7.0'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
]

intersphinx_mapping = {
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'
html_logo  = '../resources/logo.svg'

# -- Options for EPUB output
epub_show_urls = 'footnote'
