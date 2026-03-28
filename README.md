# 🛰️ gw_agn_watcher

[![PyPI version](https://badge.fury.io/py/gw-agn-watcher.svg)](https://pypi.org/project/gw-agn-watcher/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Build](https://github.com/yourusername/gw_agn_watcher/actions/workflows/publish.yml/badge.svg)](https://github.com/yourusername/gw_agn_watcher/actions)

---

### Overview

**`gw_agn_watcher`** is a Python package for the **automated crossmatching of gravitational-wave (GW) sky maps** from the LIGO–Virgo–KAGRA (LVK) Collaboration with ** ZTF alerts and AGN catalogs** using ALeRCE infrastructure.  
It enables systematic searches for **electromagnetic counterparts** to compact binary mergers, with a particular focus on mergers that may occur in **active galactic nuclei (AGN) disks**.

---

### Key Features
- **Predicts the chirp mass of the GW superevent**
- 📡 **Ingest LVK skymaps** (`.fits`, HEALPix format) 
- 🌌 **Crossmatch ZTF alerts** with AGN catalogs (e.g., Milliquas)
- 🧠 **Apply ML-based filters** using ALeRCE classifiers, Pan-STARRS morphology, and Deep Real/Bogus scores
- 📅 **Temporal and spatial filtering** relative to the GW trigger time and sky localization
- 🎯 **Host-galaxy association** and ranking based on 2σ GW distance posteriors
- 🗺️ **Visualization tools** for probability maps, candidate locations, and sky coverage
- 🔧 **Modular and extensible** — suitable for ToO planning, multi-messenger analyses, and survey follow-up

---

### Installation

```bash
pip install gw_agn_watcher
```
### Citation

If you use this package, please cite:
@ARTICLE{2026arXiv260304342B,
       author = {{Bommireddy}, Hemanth and {Forster}, Francisco and {McMahon}, Isaac and {Pavez Herrera}, Manuel and {Cartier}, Regis and {Olivares Estay}, Felipe and {Hern{\'a}ndez Garc{\'\i}a}, Lorena and {Mart{\'\i}nez Aldama}, Mary Loli and {Mu{\~n}oz Arancibia}, Alejandra},
        title = "{A Broker Integrated Algorithm for Gravitational Wave - Electromagnetic Counterpart Searches in O4a and O4b Runs}",
      journal = {arXiv e-prints},
     keywords = {High Energy Astrophysical Phenomena, Instrumentation and Methods for Astrophysics},
         year = 2026,
        month = mar,
          eid = {arXiv:2603.04342},
        pages = {arXiv:2603.04342},
          doi = {10.48550/arXiv.2603.04342},
archivePrefix = {arXiv},
       eprint = {2603.04342},
 primaryClass = {astro-ph.HE},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2026arXiv260304342B},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}


