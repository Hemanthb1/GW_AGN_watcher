# treasure_trove/pipeline.py

import os
import pandas as pd
import matplotlib.pyplot as plt
import importlib

from . import radecligo, findminclust, divide, mainquery, match_milliquas
from . import redshift, classifiers, detections, extinction
from .db import get_alerce_connection




def run_pipeline(skymap_url, milliquas_csv):
    # Download and process skymap
    skymap, skymap1, ra_deg, dec_deg, mjd_obs, event_name = radecligo.radecligo(skymap_url)
    
    # Find clusters
    num = findminclust.find_min_clusters(skymap1)
    df_out, kmeans = divide.dividemap(num, skymap1)
    
    # Query ALeRCE clusters
    conn = get_alerce_connection()
    new_df = mainquery.query_alerce_clusters(conn, df_out, mjd_obs, ra_deg, dec_deg)
    
    # Match with Milliquas
    agn = pd.read_csv(milliquas_csv)
    nagn = match_milliquas.match_with_milliquas(new_df, agn)
    
    # Redshift filtering
    res = redshift.compute_distance_redshift(skymap_url)
    res1 = redshift.filter_agn_by_redshift(nagn, res)
    res1["final_2sigma"].to_csv('redshift.csv', index=False)
    
    # Query classifiers and detections (requires database connection)
    conn = get_alerce_connection()
    
    cand = classifiers.query_classifiers(conn, res1['final_2sigma']) 
    cand.to_csv('classifiers.csv', index=False)# Pass actual conn if available
    det = detections.query_detections(cand, conn)
    
    
    # Merge and compute extinction
    final1 = pd.merge(cand, det, on=['oid'])
    final1['event_id'] = event_name
    final1.to_csv('final1.csv',index=False)
    importlib.reload(extinction)
    dust, candidates = extinction.compute_lat_extinction(final1, apply_cuts=True)
    suffix = "&count=true&page=1&perPage=1000&sortDesc=true&selectedClassifier=stamp_classifier"
    url = "https://alerce.online/?" + "&".join("oid=%s" % i for i in candidates.oid) + suffix

    
    
    return candidates, ra_deg, dec_deg, url
