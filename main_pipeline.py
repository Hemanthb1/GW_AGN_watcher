# main_pipeline.py
from input_config import url, milquas_path, db_credentials
from database_connection import connect_to_db

from radecligo import radecligo
from divide import divide_map
from findminclusters import find_min_clusters
from divide_query import query_alerce_clusters
from stamplc_classify import query_stamp_and_l_classifiers
from ps1_query import query_detections
from extinction import compute_lat_detections
from gw_distance_redshift import compute_distance_redshift
from match_milliquas import match_with_milliquas

def main():
    
    # Step 2: read & divide skymap
    skymap, df, ra_deg, dec_deg, mjd_obs, event_name = read_skymap(url)
    df_out,kmeans = divide_map(df)

    # Step 3: minimum number of clusters
    n = find_min_clusters(df_out, plot=True)

   #step4: divide_query
    # Step 3: query ALeRCE for each region
    all_candidates = []
    for region in regions:
        df = query_alerce_region(conn, region)
        all_candidates.append(df)
    import pandas as pd
    candidates = pd.concat(all_candidates, ignore_index=True)
    #step 4: agn crossmatch
    final = crossmatch_milliquas(filtered, milquas_path)
    #step 4: stamp lc cut
    stamplc= query_stamp_and_lc_classifiers(conn, new_df, batch_size=10000)
    #step 5: ps1 cut
    detections=query_detections(stamplc, conn)
    # Step 6: extinction cuts
    dfp,candidates=compute_lat_extinction(final_df, rv=3.1, apply_cuts=True)
    #step 7: gw redshift cut
    result=compute_distance_redshift(url)
 
    # Step 8: save output
    final.to_csv("final_candidates.csv", index=False)
    print("✅ Pipeline complete — results saved to final_candidates.csv")

if __name__ == "__main__":
    main()
