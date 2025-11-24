# extinction_latitudes.py
"""
Module: extinction_latitudes
Purpose: Compute ecliptic latitude, galactic latitude, and galactic extinction (A_g)
         for a catalog of objects with mean RA/Dec coordinates, and apply spatial cuts.
"""

import os
import numpy as np
import pandas as pd
import ephem
import astropy.coordinates as coord
from astropy import units as u
from dustmaps.sfd import SFDQuery, fetch
from dustmaps.config import config
from astropy import coordinates
from dust_extinction.parameter_averages import F19


def alam_fromarrays(ebv, alam_ebv):
    """Compute extinction A_lambda = E(B-V) * (A_lambda / E(B-V))."""
    alam = np.outer(ebv, alam_ebv)
    if not isinstance(ebv, np.ndarray):
        alam = alam[0]
    return alam


def compute_lat_extinction(final_df, rv=3.1, apply_cuts=True):
    """
    Compute ecliptic latitude, galactic latitude, and A_g extinction
    for each object in the input DataFrame. Optionally apply astrophysical cuts.
    """

    # --- Ensure dustmaps data is available ---
    if not os.path.exists(os.path.join(config["data_dir"], "sfd")):
        print("📥 SFD maps missing, downloading...")
        fetch()  # will download SFD maps into config['data_dir']

    ecl_lat = {}
    gal_lat = {}

    for oid in final_df.index:
        try:
            ra = final_df.loc[oid].meanra
            dec = final_df.loc[oid].meandec

            # Ecliptic and galactic latitude
            ecl_lat[oid] = np.rad2deg(
                ephem.Ecliptic(ephem.Equatorial(f"{ra/15.}", f"{dec}", epoch=ephem.J2000)).lat
            )
            gal_lat[oid] = np.rad2deg(
                ephem.Galactic(ephem.Equatorial(f"{ra/15.}", f"{dec}", epoch=ephem.J2000)).lat
            )

        except Exception as e:
            print(f"⚠️ Error processing {oid}: {e}")
            ecl_lat[oid] = np.nan
            gal_lat[oid] = np.nan

    # --- Dust extinction via SFD + Fitzpatrick (2019) ---
    coords = coordinates.SkyCoord(
        ra=final_df["meanra"], dec=final_df["meandec"], unit=(u.deg, u.deg), frame="icrs"
    )
    sfd = SFDQuery()
    ebv = sfd(coords)

    ext_model = F19(Rv=rv)
    lam = np.array([4716.7, 6165.1])  # g and r band
    x_lam = 1e4 / lam
    alam_f19 = alam_fromarrays(ebv, ext_model(x_lam / u.micron)) * rv
    A_g = alam_f19[:, 0]

    dfp = pd.DataFrame.from_dict({"ecl_lat": ecl_lat, "gal_lat": gal_lat})
    dfp["gal_A_g"] = A_g
    dfp = dfp.reset_index().rename(columns={"index": "oid"})
    dfp['oid'] = dfp['oid'].astype(str)
    final_df['oid'] = final_df['oid'].astype(str)

    # Merge with input catalog
    candidates = final_df.merge(dfp, left_index=True, right_index=True)

    # Apply astrophysical filtering
    if apply_cuts:
        before = len(candidates)
        candidates = candidates[
            ((candidates["ndet"] > 1) | (candidates["ecl_lat"] > 20))
            & (candidates["gal_lat"].abs() > 20)
            & (candidates["gal_A_g"] < 1)
        ]
        after = len(candidates)
        print(f"✅ Applied sky-plane cuts: {before} → {after} candidates retained.")

    print(f"✅ Computed extinction and latitude for {len(dfp)} sources.")
    return dfp, candidates
