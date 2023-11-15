# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys

project = 'PlanetProfile'
copyright = '2023, Steven D. Vance and Marshall J. Styczinski'
author = 'Steven D. Vance, Marshall J. Styczinski, and PlanetProfile collaborators'
release = 'v2.3.18'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.napoleon',
              'sphinxcontrib.apidoc',
              'sphinx.ext.autosummary',
              'sphinxcontrib.matlab',
              'myst_parser']
source_suffix = ['.rst', '.md']
sys.path.insert(0, os.path.abspath('../'))

apidoc_module_dir = '../'
apidoc_output_dir = 'stubs/'
apidoc_excluded_paths = ['configP*', 'setup.py']
apidoc_separate_modules = True
apidoc_module_first = True

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'configP*']


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'  # Install with pip install sphinx-rtd-theme
html_static_path = ['_static']
html_logo = '../misc/PPlogo.png'
html_favicon = '../misc/PPlogo.ico'

html_theme_options = {
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'style_nav_header_background': '#2980B9',  # Default is #2980B9
    'logo_only': True,
    # Toc options
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': -1
}
