"""
Crossmatch candidates with the Milliquas AGN catalog.
Author: Hemanth Kumar
Date: 2025-11-11
"""

import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u


def match_with_milliquas(cr, df1, output_csv='matched_milliquas.csv'):
    """
    Crossmatch candidate sources with the Milliquas catalog.

    Parameters
    ----------
   cr : dataframe
        Dataframe containing candidate objects (must include 'meanra' and 'meandec' columns).
    milliquas_csv : dataframe
        Dataframe of Milliquas catalog file (must include 'ra' and 'dec' columns).
    output_csv : str, optional
        File to save the crossmatched results. Default is 'matched_milliquas.csv'.

    Returns
    -------
    pd.DataFrame
        DataFrame of matched sources with AGN name, redshift, and separation.
    """

    # === Load Data ===
   # cr = pd.read_csv(candidates_csv)
   # df1 = pd.read_csv(milliquas_csv)

    # === Extract coordinates ===
    ra1 = np.array(cr['meanra'])
    dec1 = np.array(cr['meandec'])
    ra2 = np.array(df1['ra'])
    dec2 = np.array(df1['dec'])

    # === Build SkyCoord objects ===
    candidates = SkyCoord(ra=ra1 * u.deg, dec=dec1 * u.deg)
    milliquas = SkyCoord(ra=ra2 * u.deg, dec=dec2 * u.deg)

    # === Match to nearest source ===
    idx, d2d, d3d = candidates.match_to_catalog_sky(milliquas)

    # === Extract AGN and redshift info ===
    milliquas_name_col = df1.columns[2]  # or replace with 'Name' if known
    milliquas_z_col = df1.columns[4]     # or replace with 'z' if known

    cr['agn'] = [df1.iloc[i][milliquas_name_col] for i in idx]
    cr['z'] = [df1.iloc[i][milliquas_z_col] for i in idx]
    cr['agnsep'] = d2d.degree

    # === Filter by close match (<= 0.0008 deg ≈ 2.88 arcsec) ===
    nagn = cr[cr['agnsep'] <= 0.0008]

    # === Save and report ===
    nagn.to_csv(output_csv, index=False)
    print(f"Matched {len(nagn)} candidates to Milliquas within 2.9 arcsec.")
    print(f"Results saved to {output_csv}")

    return nagn


if __name__ == "__main__":
    # Example usage — modify paths below
    candidates_csv = '/path/to/your_candidates.csv'
    milliquas_csv = '/Users/mayhem/Downloads/milqd.csv'

    match_with_milliquas(candidates_csv, milliquas_csv)
