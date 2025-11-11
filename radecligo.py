"""
radecligo.py

Module to download and extract RA/Dec and probability information
from a LIGO/Virgo/KAGRA GW skymap FITS file.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.table import QTable
from astropy.utils.data import download_file
from astropy.io import fits
import astropy_healpix as ah
import astropy.units as u


def radecligo(url, credible_level=0.9, plot=False):
    """
    Download and process a LIGO/Virgo/KAGRA skymap FITS file.

    Parameters
    ----------
    url : str
        URL to the GW skymap FITS file.
    credible_level : float, optional
        Cumulative probability cutoff (default 0.9 for 90% region).
    plot : bool, optional
        If True, shows a scatter plot of RA vs Dec of selected pixels.

    Returns
    -------
    skymap : QTable
        Original GW skymap table truncated to the credible region.
    df : pandas.DataFrame
        DataFrame with RA, Dec, pixel index, and cumulative probability.
    ra_deg, dec_deg : ndarray
        RA and Dec of pixels (degrees).
    mjd_obs : float
        Observation MJD time from FITS header.
    event_name : str
        Extracted event name from the URL (between 'superevents' and 'files').
    """

    # --- Download and open skymap ---
    gw_file = download_file(url, cache=True)
    with fits.open(gw_file) as hdul:
        mjd_obs = hdul[1].header.get('MJD-OBS', np.nan)

    skymap = QTable.read(gw_file)
    skymap.sort('PROBDENSITY', reverse=True)

    # --- Compute cumulative probability ---
    level, ipix = ah.uniq_to_level_ipix(skymap['UNIQ'])
    nside = ah.level_to_nside(level)
    pixel_area = ah.nside_to_pixel_area(nside)
    prob = pixel_area * skymap['PROBDENSITY']
    cumprob = np.cumsum(prob)
    skymap['PROB'] = cumprob

    # --- Truncate to the credible region ---
    cutoff_index = np.searchsorted(cumprob, credible_level)
    skymap = skymap[:cutoff_index]

    # --- Convert to RA/Dec ---
    level, ipix = ah.uniq_to_level_ipix(skymap['UNIQ'])
    nside = ah.level_to_nside(level)
    ra, dec = ah.healpix_to_lonlat(ipix, nside, order='nested')
    ra_deg, dec_deg = ra.to(u.deg).value, dec.to(u.deg).value

    # --- Make DataFrame ---
    df = pd.DataFrame({
        'meanra': ra_deg,
        'meandec': dec_deg,
        'pixel_no': skymap['UNIQ'],
        'prob_contour': skymap['PROB']
    })

    # --- Extract event name from URL ---
    event_name = None
    parts = url.split('/')
    if 'superevents' in parts and 'files' in parts:
        i1, i2 = parts.index('superevents'), parts.index('files')
        event_name = parts[i1 + 1] if i2 > i1 + 1 else None

    # --- Optional plot ---
    if plot:
        plt.figure(figsize=(5, 4))
        plt.scatter(ra_deg, dec_deg, s=0.1, c='k')
        plt.xlabel("RA [deg]")
        plt.ylabel("Dec [deg]")
        plt.title(f"{event_name or 'GW event'}: {credible_level:.0%} credible region")
        plt.show()

    return skymap, df, ra_deg, dec_deg, mjd_obs, event_name


if __name__ == "__main__":
    # Example use (for testing)
    test_url = "https://gracedb.ligo.org/api/superevents/S230518h/files/bayestar.fits.gz"
    skymap, df, ra, dec, mjd, event = radecligo(test_url, plot=True)
    print(f"Event: {event}, MJD: {mjd}, Pixels: {len(df)}")
