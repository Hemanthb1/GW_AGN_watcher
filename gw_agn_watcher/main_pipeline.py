import os
import pandas as pd
import matplotlib.pyplot as plt
import importlib

from . import radecligo, findminclust, divide, mainquery, match_milliquas
from . import redshift, classifiers, detections, extinction
from .db import get_alerce_connection


def run_pipeline(skymap_url, milliquas_csv):
    print("🚀 Starting GW–AGN crossmatching pipeline...")
    print(f"🔗 Skymap: {skymap_url}")
    print(f"📂 Milliquas catalog: {milliquas_csv}\n")

    # --- Step 1: Download and process skymap ---
    skymap, skymap1, ra_deg, dec_deg, mjd_obs, event_name = radecligo.radecligo(skymap_url)
    print(f"✅ Loaded skymap '{event_name}' with {len(skymap1)} pixels in 90% region.\n")

    # --- Step 2: Find clusters in the skymap ---
    num = findminclust.find_min_clusters(skymap1)
    df_out, kmeans = divide.dividemap(num, skymap1)
    print(f"✅ Divided into {len(df_out)} clusters (k={num}).\n")

    # --- Step 3: Query ALeRCE clusters ---
    conn = get_alerce_connection()
    new_df = mainquery.query_alerce_clusters(conn, df_out, mjd_obs, ra_deg, dec_deg)
    print(f"✅ Queried ALeRCE: {len(new_df)} sources retrieved from cluster regions.\n")

    if new_df.empty:
        print("⚠️ No ALeRCE sources found near GW localization — stopping early.")
        return pd.DataFrame(), ra_deg, dec_deg, None

    # --- Step 4: Match with Milliquas ---
    agn = pd.read_csv(milliquas_csv)
    nagn = match_milliquas.match_with_milliquas(new_df, agn)
    print(f"✅ Matched with Milliquas: {len(nagn)} candidate AGNs after spatial crossmatch.\n")

    if nagn.empty:
        print("⚠️ No Milliquas matches found — stopping early.")
        return nagn, ra_deg, dec_deg, None

    # --- Step 5: Redshift filtering ---
    res = redshift.compute_distance_redshift(skymap_url)
    res1 = redshift.filter_agn_by_redshift(nagn, res)
    res1["final_2sigma"].to_csv("redshift.csv", index=False)
    print(f"✅ Redshift filtering complete: {len(res1)} objects remain within 2σ distance.\n")

    if res1.empty or "final_2sigma" not in res1:
        print("⚠️ No AGNs passed redshift cut — stopping early.")
        return res1, ra_deg, dec_deg, None

    # --- Step 6: Query classifiers and detections ---
    conn = get_alerce_connection()
    cand = classifiers.query_classifiers(conn, res1["final_2sigma"])
    print(f"✅ Classifiers queried: {len(cand)} objects classified (stamp/lc).\n")

    cand.to_csv("classifiers.csv", index=False)

    det = detections.query_detections(cand, conn)
    print(f"✅ Detections queried: {len(det)} rows retrieved from database.\n")

    # --- Step 7: Merge and compute extinction ---
    final1 = pd.merge(cand, det, on=["oid"], how="inner")
    final1["event_id"] = event_name
    final1.to_csv("final1.csv", index=False)
    print(f"✅ Merged classifiers + detections: {len(final1)} objects.\n")

    if final1.empty:
        print("⚠️ No valid objects for extinction step — stopping early.")
        return final1, ra_deg, dec_deg, None

    importlib.reload(extinction)
    dust, candidates = extinction.compute_lat_extinction(final1, apply_cuts=True)
    print(f"✅ Extinction computed for {len(dust)} sources.")
    print(f"✅ {len(candidates)} sources remain after sky-plane & dust cuts.\n")

    if candidates.empty:
        print("⚠️ No candidates remain after extinction filtering. Returning empty set.")
        return candidates, ra_deg, dec_deg, None

    # --- Step 8: Generate ALeRCE viewer URL ---
    suffix = "&count=true&page=1&perPage=1000&sortDesc=true&selectedClassifier=stamp_classifier"
    url = "https://alerce.online/?" + "&".join(f"oid={i}" for i in candidates.oid_x) + suffix
    print(f"🔗 Final ALeRCE viewer link generated.\n")

    print("🏁 Pipeline completed successfully.")
    return candidates, ra_deg, dec_deg, url
