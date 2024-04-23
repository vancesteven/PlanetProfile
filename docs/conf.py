# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys

sys.path.append(os.path.dirname(__file__))
from misc.upgreek import upgreekDefs

project = 'PlanetProfile'
copyright = '2024, Steven D. Vance and Marshall J. Styczinski'
author = 'Steven D. Vance, Marshall J. Styczinski, and PlanetProfile collaborators'
release = 'v2.5.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.napoleon',
              'sphinxcontrib.apidoc',
              'sphinxcontrib.matlab',
              'sphinx.ext.mathjax',
              'myst_parser']
source_suffix = ['.rst', '.md']
mathjaxVer = 2

_HERE = os.path.dirname(__file__)
_ROOT_DIR = os.path.abspath(os.path.join(_HERE, '..'))
_PACKAGE_DIR = os.path.abspath(os.path.join(_HERE, '../PlanetProfile'))
_TRAJEC_DIR = os.path.abspath(os.path.join(_HERE, '../PlanetProfile/TrajecAnalysis'))

sys.path.insert(0, _ROOT_DIR)
sys.path.insert(0, _PACKAGE_DIR)
sys.path.insert(0, _TRAJEC_DIR)

apidoc_module_dir = _ROOT_DIR
apidoc_output_dir = 'stubs'
apidoc_template_dir = 'templates'
apidoc_excluded_paths = ['configP*', 'setup.py']
apidoc_separate_modules = True
apidoc_module_first = True

templates_path = ['templates']
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


# -- Options for LaTeX math formatting-----------
# https://www.sphinx-doc.org/en/master/latex.html

# This block is used when the sphinx.ext.mathjax extension is loaded
latex_packages = ['upgreek', 'mhchem']
if mathjaxVer == 3:
    # Does not support STIX fonts as of v3.2.
    mathjax_path='https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js'
    mathjax3_config = {
      'loader': {'load': [f'[tex]/{pkg}' for pkg in latex_packages]},
      'tex': {'packages': {'[+]': [latex_packages]}}
    }
else:
    mathjax_path = 'https://cdn.jsdelivr.net/npm/mathjax@2/MathJax.js?config=TeX-AMS-MML_HTMLorMML'
    latex_packages.remove('upgreek')
    latex_packages += ['unicode']
    mathjax2_config = {
        'extensions': ['tex2jax.js'] , #+ [f'TeX/{pkg}.js' for pkg in latex_packages],
        'TeX': {'extensions': [f'{pkg}.js' for pkg in latex_packages], 'Macros': upgreekDefs},
        'HTML-CSS': {'fonts': ['STIX-Web']}
    }
