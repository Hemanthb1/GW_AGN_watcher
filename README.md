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
