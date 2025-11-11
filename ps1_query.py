# query_detections.py
"""
Module: query_detections
Purpose: Query ALeRCE detection and PS1 data for a given list of oids,
         applying sgscore, distance, and DRB cuts.
"""

import pandas as pd

def query_detections(stamplc, conn):
    """
    Query detections and PS1 matches for a given set of object IDs (oids).

    Parameters
    ----------
    stamplc : pandas.DataFrame
        DataFrame containing an 'oid' column with object IDs.
    conn : psycopg2 connection
        Active database connection to the ALeRCE PostgreSQL database.

    Returns
    -------
    detections : pandas.DataFrame
        Filtered detections joined with PS1 metadata.
    """

    if stamplc.empty:
        print("⚠️ No oids provided — returning empty DataFrame.")
        return pd.DataFrame()

    # Build comma-separated list of OIDs for the query
    oid_list = ",".join([f"'{x}'" for x in stamplc["oid"].unique()])

    query = f"""
        SELECT
            det.oid, det.drb, det.fid,
            det.mjd, det.magpsf, det.sigmapsf,
            det.has_stamp,
            ps1.sgscore1, ps1.distpsnr1
        FROM
            (SELECT *
             FROM detection
             WHERE oid IN ({oid_list})
            ) AS det
        INNER JOIN
            (SELECT *
             FROM ps1_ztf
             WHERE oid IN ({oid_list})
            ) AS ps1
        ON det.oid = ps1.oid
        WHERE
            (ps1.sgscore1 < 0.5 OR ps1.distpsnr1 > 1)
            AND det.drb > 0.5;
    """

    # Execute SQL and return DataFrame
    detections = pd.read_sql_query(query, conn)

    # Drop duplicate OIDs to keep one per source
    detections = detections.drop_duplicates(subset="oid", keep="first")

    print(f"✅ Retrieved {len(detections)} detections after filtering.")
    return detections
