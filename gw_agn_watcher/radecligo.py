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
    gw_skymap = download_file(url, cache=True)
    head=fits.open(gw_skymap)
    time=head[1].header['MJD-OBS']
    skymap = QTable.read(gw_skymap)

    skymap2 = QTable.read(gw_skymap)
    skymap.sort('PROBDENSITY', reverse=True)
    #skymap.sort('PROB', reverse=True)
    level, ipix = ah.uniq_to_level_ipix(skymap['UNIQ'])
    nside = ah.level_to_nside(level)
    pixel_area = ah.nside_to_pixel_area(ah.level_to_nside(level))
    prob = pixel_area * skymap['PROBDENSITY']
    cumprob = np.cumsum(prob)
    i = cumprob.searchsorted(0.9)
    skymap['PROB']=cumprob
    skymap = skymap[:i]

    skymap.sort('UNIQ')
    skymap = skymap['UNIQ','PROB']
    level, ipix = ah.uniq_to_level_ipix(skymap['UNIQ'])
    nside = ah.level_to_nside(level)

    # Ensure little-endian for safety
    def ensure_little_endian(data):
        if data.dtype.byteorder == '>':
            data = data.byteswap().newbyteorder('<')
        return data

    skymap= ensure_little_endian(np.array(skymap))

    # RA/Dec conversion
    ra, dec = ah.healpix_to_lonlat(ipix, nside, order='nested')
    df = pd.DataFrame(skymap)
    ra_deg = np.rad2deg(ra.value)
    dec_deg = np.rad2deg(dec.value)
    selected_elements=df['PROB']
    arr=df['UNIQ']

    # Build pandas DataFrame
    skymap1 = pd.DataFrame({
        'meanra': ra_deg.flatten(),
        'meandec': dec_deg.flatten(),
        'pixel_no': arr.values.flatten(),
        'prob_contour': selected_elements.values.flatten()
    })

    # Optional plotting
    if plot:
        plt.figure(figsize=(8,4))
        plt.scatter(ra_deg, dec_deg, s=0.5)
        plt.xlabel('RA [deg]')
        plt.ylabel('Dec [deg]')
        plt.title('GW Skymap Pixels')
        plt.show()

    # Extract event name from URL
    strings_list = url.split('/')
    start_index, end_index = -1, -1
    for i, string in enumerate(strings_list):
        if string == 'superevents':
            start_index = i
        elif string == 'files' and start_index != -1:
            end_index = i
            break
    event_name = 'unknown'
    if start_index != -1 and end_index != -1:
        event_name = ' '.join(strings_list[start_index+1:end_index]).strip()

    return skymap, skymap1, ra_deg, dec_deg, time, event_name

# Example usage if run as script
if __name__ == "__main__":
    test_url = "https://gracedb.ligo.org/api/superevents/S230518h/files/bayestar.fits.gz"
    skymap, df, ra, dec, mjd, event = radecligo(test_url, plot=True)
    print(f"Event: {event}, MJD: {mjd}, Pixels: {len(df)}")

