# GW_AGN_Watcher API Documentation

**Version:** 1.2.0  
**Author:** Hemanth Kumar  
**License:** MIT  
**Description:** Python tools for automated crossmatching of gravitational-wave (GW) sky maps with ZTF alerts and AGN catalogs.

---

## Installation

```bash
pip install gw_agn_watcher
```

### Dependencies

- numpy
- astropy
- matplotlib
- scipy
- requests
- pandas
- astroquery
- alphashape
- psycopg2-binary
- ephem
- dustmaps
- dust_extinction
- ligo.skymap

---

## Package Overview

```python
from gw_agn_watcher import (
    divide,
    main_pipeline,
    radecligo,
    findminclust,
    mainquery,
    match_milliquas,
    classifiers,
    db,
    detections,
    extinction,
    mass_estimation,
    redshift,
    get_public_alerts
)
```

---

## Core Modules

### 1. `radecligo` - Skymap Processing

Download and process LIGO/Virgo/KAGRA skymap FITS files.

#### `radecligo(url, credible_level=0.9, plot=False)`

Download and extract RA/Dec and probability information from a GW skymap.

**Parameters:**
- `url` (str): URL to the GW skymap FITS file
- `credible_level` (float): Cumulative probability cutoff (default: 0.9)
- `plot` (bool): If True, displays a scatter plot of the skymap

**Returns:**
- `skymap` (QTable): Original GW skymap table truncated to the credible region
- `df` (DataFrame): DataFrame with columns `meanra`, `meandec`, pixel index, and cumulative probability
- `ra_deg` (ndarray): RA values in degrees
- `dec_deg` (ndarray): Dec values in degrees
- `mjd_obs` (float): Observation MJD time from FITS header
- `event_name` (str): Extracted event name from the URL

**Example:**
```python
from gw_agn_watcher import radecligo

url = "https://example.com/skymap.fits"
skymap, df, ra_deg, dec_deg, mjd_obs, event_name = radecligo(url, credible_level=0.9)
```

---

### 2. `divide` - Sky Map Segmentation

Divide a GW skymap into spatial regions using K-Means clustering.

#### `dividemap(num_regions, df, plot=False, random_state=42)`

**Parameters:**
- `num_regions` (int): Number of clusters to create
- `df` (DataFrame): Input DataFrame with `meanra` and `meandec` columns
- `plot` (bool): If True, displays the clustered skymap
- `random_state` (int): Random seed for reproducibility (default: 42)

**Returns:**
- `df_out` (DataFrame): Input DataFrame with additional `cluster_label` column
- `kmeans` (KMeans): The fitted KMeans model

**Example:**
```python
from gw_agn_watcher import divide

df_clustered, kmeans = divide.dividemap(num_regions=50, df=skymap_df, plot=True)
```

---

### 3. `findminclust` - Adaptive Clustering

Find the minimum number of clusters that satisfy polygon constraints.

#### `find_min_clusters(df, max_vertices=100, max_diameter=25, max_try=200, random_state=42, plot=False)`

**Parameters:**
- `df` (DataFrame): DataFrame with `meanra` and `meandec` columns
- `max_vertices` (int): Maximum vertices allowed per polygon (default: 100)
- `max_diameter` (float): Maximum angular diameter in degrees (default: 25)
- `max_try` (int): Maximum number of cluster attempts (default: 200)
- `random_state` (int): Random seed (default: 42)
- `plot` (bool): If True, displays clustering results

**Returns:**
- `int`: Minimum valid number of clusters

#### `check_clusters(df, n_clusters, max_vertices=100, max_diameter=25, random_state=42)`

Check if all polygons from clustering satisfy vertex and diameter limits.

**Returns:** `bool` - True if all clusters satisfy constraints

#### `polygon_diameter(polygon)`

Return maximum angular separation (degrees) among polygon vertices.

**Parameters:** `polygon` (shapely geometry)  
**Returns:** `float`

---

### 4. `mainquery` - ALeRCE Database Queries

Query ALeRCE for objects within sky map regions.

#### `query_alerce_clusters(conn, skymap_df, time, ra, dec, ndays=200, alpha=0.01)`

Divide the sky map into alpha-shape polygons and query ALeRCE for objects.

**Parameters:**
- `conn`: Active ALeRCE database connection
- `skymap_df` (DataFrame): Sky map data with cluster labels
- `time` (float): Reference time (MJD)
- `ra` (float): Right ascension of the event
- `dec` (float): Declination of the event
- `ndays` (int): Time window in days (default: 200)
- `alpha` (float): Alpha shape parameter (default: 0.01)

**Returns:**
- `new_df` (DataFrame): Query results with object data

---

### 5. `detections` - Detection Queries

#### `query_detections(stamplc, conn)`

Query ALeRCE detections and PS1 matches for a list of object IDs.

**Parameters:**
- `stamplc` (DataFrame): DataFrame containing `oid` column
- `conn`: Active ALeRCE PostgreSQL database connection

**Returns:**
- `detections` (DataFrame): Filtered detections with PS1 metadata

---

### 6. `classifiers` - Object Classification

#### `query_classifiers(conn, new_df, batch_size=10000)`

Query ALeRCE for stamp_classifier and lc_classifier results.

**Parameters:**
- `conn`: Active ALeRCE database connection
- `new_df` (DataFrame): DataFrame with `oid` column (or index as `oid`)
- `batch_size` (int): Number of OIDs per batch (default: 10,000)

**Returns:**
- `candidates` (DataFrame): Combined SN/AGN/QSO/Blazar/SLSN sources with probabilities

---

### 7. `match_milliquas` - Catalog Crossmatching

#### `match_with_milliquas(cr, df1, output_csv='matched_milliquas.csv')`

Crossmatch candidate sources with the Milliquas AGN catalog.

**Parameters:**
- `cr` (DataFrame): Candidate objects with `meanra` and `meandec` columns
- `df1` (DataFrame): Milliquas catalog with `ra` and `dec` columns
- `output_csv` (str): Output file path (default: 'matched_milliquas.csv')

**Returns:**
- DataFrame: Matched sources with AGN name, redshift, and separation

**Example:**
```python
from gw_agn_watcher import match_milliquas

matched = match_milliquas.match_with_milliquas(candidates_df, milliquas_df)
```

---

### 8. `extinction` - Extinction Corrections

#### `compute_lat_extinction(final_df, rv=3.1, apply_cuts=True)`

Compute ecliptic latitude, galactic latitude, and A_g extinction for each object.

**Parameters:**
- `final_df` (DataFrame): Must contain `meanra`, `meandec`, `ndet` columns; index should be `oid`
- `rv` (float): R_V value for extinction law (default: 3.1)
- `apply_cuts` (bool): Apply astrophysical cuts (default: True)

**Returns:**
- `dfp` (DataFrame): Columns `oid`, `ecl_lat`, `gal_lat`, `gal_A_g`
- `candidates` (DataFrame): Filtered DataFrame (if apply_cuts=True)

---

### 9. `redshift` - Redshift Filtering

#### `compute_distance_redshift(url)`

Compute distance and redshift bounds from a GW skymap URL.

**Parameters:** `url` (str) - GW skymap URL  
**Returns:** `dict` with distance stats and redshift limits

#### `filter_agn_by_redshift(nagn, z_bounds)`

Filter crossmatched AGNs within given GW redshift ranges.

**Parameters:**
- `nagn` (DataFrame): AGN candidates with redshift data
- `z_bounds` (dict): Redshift bounds from `compute_distance_redshift`

**Returns:** `dict` with AGN subsets for three sigma ranges

---

### 10. `mass_estimation` - Mass Classification

#### Class: `MassEstimator`

Handle chirp mass classification using machine learning.

**Constructor:** `MassEstimator(data_npix, downsample_nside)`

**Methods:**

- `load_simulation_data(sim_dirs, det_data_file, inj_data_file, snr_threshold, overwrite)`
  - Load BAYESTAR simulation data for training

- `train(model_path, nbins, mchirp_max, save)`
  - Train Random Forest Classifier

- `predict_mass(downsample_skymap)`
  - Predict chirp mass PDF for a single event
  - **Returns:** `pdf` (ndarray), `midpoints` (ndarray)

#### `run_simulation(n_events, outdir, min_mass, max_mass, max_distance, snr_threshold, downsample_nside)`

Run BAYESTAR simulations for training data generation.

---

### 11. `db` - Database Connections

#### `get_alerce_connection()`

Get connection to ALeRCE database. Attempts remote credentials first, falls back to local parameters.

**Returns:** Database connection object (psycopg2)

**Example:**
```python
from gw_agn_watcher import db

conn = db.get_alerce_connection()
```

---

### 12. `main_pipeline` - End-to-End Pipeline

#### `run_pipeline(skymap_url, milliquas_csv, sigma_cut="2sigma")`

Execute the complete GW-AGN crossmatching pipeline.

**Parameters:**
- `skymap_url` (str): URL to the GW skymap FITS file
- `milliquas_csv` (str): Path to Milliquas catalog CSV
- `sigma_cut` (str): Sigma cut for filtering (default: "2sigma")

**Returns:**
- `candidates` (DataFrame): Final candidate AGN list
- `ra_deg` (float): Event RA in degrees
- `dec_deg` (float): Event Dec in degrees
- `url` (str): Skymap URL used

**Example:**
```python
from gw_agn_watcher import main_pipeline

candidates, ra, dec, url = main_pipeline.run_pipeline(
    skymap_url="https://example.com/skymap.fits",
    milliquas_csv="milliquas.csv"
)
```

---

### 13. `get_public_alerts` - LVK Alert Parsing

Utilities for parsing LVK public alerts from GCN.

**Functions:**
- `get_params_for_group(voevent_xml, name)` - Extract group parameters
- `get_params_for_param(voevent_xml, name)` - Extract specific parameter
- `get_skymap(url)` - Retrieve skymap from URL
- `get_skymap_stats(skymap)` - Compute skymap statistics
- `get_info(superevent)` - Extract superevent information

---

## Workflow Summary

```
1. Download skymap      → radecligo()
2. Segment sky map      → findminclust.find_min_clusters() + divide.dividemap()
3. Query candidates     → mainquery.query_alerce_clusters()
4. Get classifications  → classifiers.query_classifiers()
5. Get detections       → detections.query_detections()
6. Crossmatch AGN       → match_milliquas.match_with_milliquas()
7. Filter by redshift   → redshift.filter_agn_by_redshift()
8. Apply extinction     → extinction.compute_lat_extinction()
```

---

## Quick Start Example

```python
from gw_agn_watcher import main_pipeline

# Run the full pipeline
candidates, ra, dec, url = main_pipeline.run_pipeline(
    skymap_url="https://gracedb.ligo.org/superevents/S230529M/files/bayestar.fits",
    milliquas_csv="path/to/milliquas.csv",
    sigma_cut="2sigma"
)

print(f"Found {len(candidates)} candidate AGN counterparts")
print(candidates[['oid', 'agn_name', 'redshift']].head())
```

---

## Support

- **Repository:** https://github.com/Hemanthb1/GW_AGN_watcher
- **Issues:** https://github.com/Hemanthb1/GW_AGN_watcher/issues
