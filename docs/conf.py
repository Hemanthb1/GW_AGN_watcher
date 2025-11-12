import os
import sys
sys.path.insert(0, os.path.abspath('../GW_AGN_watcher/gw-agn-watcher'))  # path to your package

project = 'gw-agn-watcher'
author = 'Hemanth Kumar'
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
]

html_theme = 'sphinx_rtd_theme'
