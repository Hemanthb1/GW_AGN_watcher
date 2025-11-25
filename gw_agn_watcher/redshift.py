import re
import numpy as np
from astropy.utils.data import download_file
from ligo.skymap.io import fits
from astropy.coordinates import Distance
from ligo.skymap.distance import parameters_to_marginal_moments
from astropy import units as u
import astropy.cosmology.units as cu
from astropy.cosmology import WMAP9
import pandas as pd


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
        Dictionary containing distance mean/std and redshift bounds.
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
    
    ZMIN_FLOOR = 1e-11
    DIST_FLOOR = 0.01 * u.Mpc  # 10 kpc — too small for reliable z inversion
    ZMAX1 = 3000
    ZMAX2 = 20000
    
    # --- Convert to redshift bounds (safe with logs) ---
    if distance_lower < DIST_FLOOR:
        print(f"⚠️ {event_name}: distance_lower below {DIST_FLOOR}, setting z_min=0.0")
        z_min = 0.0
    else:
        z_min = distance_lower.to(cu.redshift, cu.redshift_distance(WMAP9, kind='comoving')).value
    
    if distance_upper < DIST_FLOOR:
        print(f"⚠️ {event_name}: distance_upper below {DIST_FLOOR}, setting z_max=0.0")
        z_max = 0.0
    else:
        z_max = distance_upper.to(cu.redshift, cu.redshift_distance(WMAP9, kind='comoving')).value
    
    if distance_lower1 < DIST_FLOOR:
        print(f"⚠️ {event_name}: distance_lower1 below {DIST_FLOOR}, setting z_min1=0.0")
        z_min1 = 0.0
    else:
        z_min1 = distance_lower1.to(cu.redshift, cu.redshift_distance(WMAP9, kind='comoving')).value
    
    if distance_upper1 < DIST_FLOOR:
        print(f"⚠️ {event_name}: distance_upper1 below {DIST_FLOOR}, setting z_max1=0.0")
        z_max1 = 0.0
    else:
        z_max1 = distance_upper1.to(cu.redshift, cu.redshift_distance(WMAP9, kind='comoving', zmax=ZMAX1)).value
    
    if distance_lower2 < DIST_FLOOR:
        print(f"⚠️ {event_name}: distance_lower2 below {DIST_FLOOR}, setting z_min2=0.0")
        z_min2 = 0.0
    else:
        z_min2 = distance_lower2.to(cu.redshift, cu.redshift_distance(WMAP9, kind='comoving', zmin=ZMIN_FLOOR)).value
    
    if distance_upper2 < DIST_FLOOR:
        print(f"⚠️ {event_name}: distance_upper2 below {DIST_FLOOR}, setting z_max2=0.0")
        z_max2 = 0.0
    else:
        z_max2 = distance_upper2.to(cu.redshift, cu.redshift_distance(WMAP9, kind='comoving', zmax=ZMAX2)).value
    
    
    # --- Compile results ---
    result = {
        "event_name": event_name,
        "distmean_Mpc": distmean,
        "diststd_Mpc": diststd,
        "z_min": z_min,
        "z_max": z_max,
        "z_min1": z_min1,
        "z_max1": z_max1,
        "z_min2": z_min2,
        "z_max2": z_max2,
    }

    # Convert to redshift bounds
   

    print(f"Event: {event_name}")
    print(f"Mean distance: {distmean:.2f} ± {diststd:.2f} Mpc")
    print(f"Redshift range (1.28σ): {z_min.value:.4f} – {z_max.value:.4f}")
    print(f"Redshift range (2σ): {z_min1.value:.4f} – {z_max1.value:.4f}")
    print(f"Redshift range (kσ): {z_min2.value:.4f} – {z_max2.value:.4f}")

    return result


def filter_agn_by_redshift(nagn, z_bounds):
    """
    Filter crossmatched AGNs within given GW redshift ranges.

    Parameters
    ----------
    nagn : pandas.DataFrame
        DataFrame with at least a 'z' column.
    z_bounds : dict
        Dictionary from compute_distance_redshift().

    Returns
    -------
    dict
        Filtered AGN subsets for 1.28σ, 2σ, and kσ ranges.
    """
    z_min, z_max = z_bounds["z_min"], z_bounds["z_max"]
    z_min1, z_max1 = z_bounds["z_min1"], z_bounds["z_max1"]
    z_min2, z_max2 = z_bounds["z_min2"], z_bounds["z_max2"]

    final_1sigma = nagn[(nagn["z"] >= z_min) & (nagn["z"] < z_max)]
    final_2sigma = nagn[(nagn["z"] >= z_min1) & (nagn["z"] < z_max1)]
    final_ksigma = nagn[(nagn["z"] >= z_min2) & (nagn["z"] < z_max2)]

    print(f"1.28σ AGNs: {len(final_1sigma)} | 2σ AGNs: {len(final_2sigma)} | kσ AGNs: {len(final_ksigma)}")

    return {
        "final_1sigma": final_1sigma,
        "final_2sigma": final_2sigma,
        "final_ksigma": final_ksigma,
    }


if __name__ == "__main__":
    # Example usage
    url = "https://gracedb.ligo.org/api/superevents/S230518h/files/bayestar.fits.gz"
    z_bounds = compute_distance_redshift(url)

    # Example: Load crossmatched AGN file
    nagn = pd.read_csv("crossmatched_agn.csv")  # Must contain 'z' column

    filtered = filter_agn_by_redshift(nagn, z_bounds)

    # Save outputs
    for key, df in filtered.items():
        df.to_csv(f"{key}.csv", index=False)
    print("✅ Saved filtered AGN subsets by GW redshift bounds.")
