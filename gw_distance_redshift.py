# gw_distance_redshift.py
"""
Compute distance and redshift range from a LIGO/Virgo skymap URL.

Author: Hemanth Kumar
Date: 2025-11-11
"""

import re
import numpy as np
from astropy.utils.data import download_file
from ligo.skymap.io import fits
from astropy.coordinates import Distance
from ligo.skymap.distance import parameters_to_marginal_moments
from astropy import units as u
import astropy.cosmology.units as cu
from astropy.cosmology import WMAP9


def compute_distance_redshift(url):
    """
    Compute the distance and redshift bounds from a GW skymap URL.

    Parameters
    ----------
    url : str
        URL to the GW skymap FITS file.

    Returns
    -------
    dict
        Dictionary containing:
        - 'event_name': GW event identifier
        - 'distmean': mean distance (Mpc)
        - 'diststd': standard deviation (Mpc)
        - 'z_min', 'z_max' : 1.28σ redshift bounds
        - 'z_min1', 'z_max1' : 2σ redshift bounds
        - 'z_min2', 'z_max2' : kσ redshift bounds
    """

    # Extract event name
    strings_list = url.split('/')
    start_index, end_index = -1, -1
    for i, string in enumerate(strings_list):
        if string == 'superevents':
            start_index = i
        elif string == 'files':
            end_index = i
            break

    event_name = None
    if start_index != -1 and end_index != -1:
        event_name = ' '.join(strings_list[start_index + 1:end_index]).strip()

    # Download and read FITS skymap
    file = download_file(url, cache=True)
    skymap, metadata = fits.read_sky_map(file, nest=False, distances=True)

    map_struct = {
        "prob": skymap[0],
        "distmu": skymap[1],
        "distsigma": skymap[2],
        "distnorm": skymap[3],
    }

    distmean, diststd = parameters_to_marginal_moments(
        map_struct["prob"],
        map_struct["distmu"],
        map_struct["distsigma"]
    )

    distance = Distance(distmean * u.Mpc)
    sig = distmean / diststd
    k = 3 if sig > 3 else np.round(sig, 3)
    print(f"k = {k}")

    # Define distance bounds
    distance_lower = Distance((distmean - 1.28 * diststd) * u.Mpc)
    distance_upper = Distance((distmean + 1.28 * diststd) * u.Mpc)
    distance_lower1 = Distance((distmean - 2 * diststd) * u.Mpc)
    distance_upper1 = Distance((distmean + 2 * diststd) * u.Mpc)
    distance_lower2 = Distance(max(distmean - k * diststd, 0) * u.Mpc)
    distance_upper2 = Distance((distmean + k * diststd) * u.Mpc)

    # Convert to redshift bounds
    z_min = distance_lower.to(cu.redshift, cu.redshift_distance(WMAP9, kind="comoving"))
    z_max = distance_upper.to(cu.redshift, cu.redshift_distance(WMAP9, kind="comoving"))
    z_min1 = distance_lower1.to(cu.redshift, cu.redshift_distance(WMAP9, kind="comoving"))
    z_max1 = distance_upper1.to(cu.redshift, cu.redshift_distance(WMAP9, kind="comoving", zmax=3000))
    z_min2 = distance_lower2.to(cu.redshift, cu.redshift_distance(WMAP9, kind="comoving", zmin=1e-12))
    z_max2 = distance_upper2.to(cu.redshift, cu.redshift_distance(WMAP9, kind="comoving", zmax=20000))

    result = {
        "event_name": event_name,
        "distmean_Mpc": distmean,
        "diststd_Mpc": diststd,
        "z_min": z_min.value,
        "z_max": z_max.value,
        "z_min1": z_min1.value,
        "z_max1": z_max1.value,
        "z_min2": z_min2.value,
        "z_max2": z_max2.value,
    }

    print(f"Event: {event_name}")
    print(f"Mean distance: {distmean:.2f} ± {diststd:.2f} Mpc")
    print(f"Redshift range (1.28σ): {z_min.value:.4f} – {z_max.value:.4f}")
    print(f"Redshift range (2σ): {z_min1.value:.4f} – {z_max1.value:.4f}")
    print(f"Redshift range (kσ): {z_min2.value:.4f} – {z_max2.value:.4f}")

    return result


if __name__ == "__main__":
    # Example usage
    url = "https://gracedb.ligo.org/api/superevents/S230518h/files/bayestar.fits.gz"
    compute_distance_redshift(url)
