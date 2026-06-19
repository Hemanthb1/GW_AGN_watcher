# 🛰️ gw_agn_watcher

[![PyPI version](https://badge.fury.io/py/gw-agn-watcher.svg)](https://pypi.org/project/gw-agn-watcher/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Build](https://github.com/yourusername/gw_agn_watcher/actions/workflows/publish.yml/badge.svg)]
[![Build Sphinx Docs](https://github.com/Hemanthb1/GW_AGN_watcher/actions/workflows/build_docs.yml/badge.svg)](https://github.com/Hemanthb1/GW_AGN_watcher/actions/workflows/build_docs.yml)

---

### Overview

**`gw_agn_watcher`** is a Python package for the **automated crossmatching of gravitational-wave (GW) sky maps** from the LIGO–Virgo–KAGRA (LVK) Collaboration with ** ZTF alerts and AGN catalogs** using ALeRCE infrastructure.  
It enables systematic searches for **electromagnetic counterparts** to compact binary mergers, with a particular focus on mergers that may occur in **active galactic nuclei (AGN) disks**.

---

### Key Features
- **Estimates the chirp mass of the GW superevent**
- 📡 **Ingest LVK skymaps** (`.fits`, HEALPix format) 
- 🌌 **Crossmatch ZTF alerts** with AGN catalogs (e.g., Milliquas)
- 🧠 **Apply ML-based filters** using ALeRCE classifiers, Pan-STARRS morphology, and Deep Real/Bogus scores
- 📅 **Temporal and spatial filtering** relative to the GW trigger time and sky localization
- 🎯 **Host-galaxy association** and ranking based on 2σ GW distance posteriors
- 🗺️ **Visualization tools** for probability maps, candidate locations, and sky coverage
- 🔧 **Modular and extensible** — suitable for ToO planning, multi-messenger analyses, and survey follow-up

---
## Outputs

The `run_pipeline()` function returns:

```python
final_cand, ra_deg, dec_deg, url, mjd_obs
```

### `final_cand`

A pandas DataFrame containing the final list of GW-associated AGN candidates after all filtering steps:

* Spatial localization within the GW skymap
* Crossmatch with the Milliquas AGN catalog
* Distance/redshift consistency with the GW event
* Stamp and light-curve classifier filtering
* Galactic plane and extinction cuts

This table represents the primary science output of the pipeline.

### `ra_deg`, `dec_deg`

Sky coordinates (in degrees) corresponding to the maximum-probability location of the gravitational-wave event.

```python
ra_deg   # Right Ascension (deg)
dec_deg  # Declination (deg)
```

These values can be used for visualization and follow-up observations.

### `url`

A pre-generated ALeRCE viewer URL containing all final candidate objects.

Opening this URL allows interactive inspection of the selected candidates through the ALeRCE portal.

Example:

```python
https://alerce.online/?oid=...
```

### `mjd_obs`

Modified Julian Date (MJD) corresponding to the gravitational-wave event time.

```python
mjd_obs
```

This value is used internally for temporal filtering and can be useful for downstream analyses.

---

## Intermediate Output Files

Several intermediate products are saved during execution for debugging, validation, and reproducibility.

| File                  | Description                                                               |
| --------------------- | ------------------------------------------------------------------------- |
| `redshift.csv`        | Objects passing the default redshift calculation step                     |
| `redshift_1sigma.csv` | Objects consistent with the GW distance at 1σ                             |
| `redshift_2sigma.csv` | Objects consistent with the GW distance at 2σ                             |
| `redshift_ksigma.csv` | Objects consistent with the user-defined kσ distance cut                  |
| `classifiers.csv`     | Results from ALeRCE stamp and light-curve classifiers                     |
| `final1.csv`          | Merged classifier and detection information prior to extinction filtering |

---

## Early Exit Conditions

The pipeline will terminate early and return an empty DataFrame if:

* No ALeRCE sources are found within the GW localization region.
* No Milliquas AGNs are spatially matched.
* No objects pass the selected redshift consistency criterion.
* No candidates survive extinction and Galactic plane filtering.

In these cases the returned candidate table will be empty, indicating no viable GW–AGN counterparts were identified.
---


### Installation

```bash
pip install gw_agn_watcher
```
### Citation

If you use this package, please cite:
```
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

```

