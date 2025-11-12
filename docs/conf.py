import os
import sys
sys.path.insert(0, os.path.abspath('../'))  # path to your package

project = 'gw_agn_watcher'
author = 'Hemanth Kumar'
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
]

html_theme = 'sphinx_rtd_theme'
