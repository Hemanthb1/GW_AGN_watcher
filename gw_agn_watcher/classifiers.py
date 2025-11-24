# query_classifiers.py
import pandas as pd
import numpy as np

def query_classifiers(conn, new_df, batch_size=10000):
    """
    Query ALeRCE for both stamp_classifier and lc_classifier results,
    filter high-probability sources, and merge them into one DataFrame.

    Parameters
    ----------
    conn : psycopg2 connection
        Active connection to the ALeRCE database.
    new_df : pd.DataFrame
        DataFrame with 'oid' (or index = 'oid') to query.
    batch_size : int, optional
        Number of OIDs per batch (default 10,000).

    Returns
    -------
    candidates : pd.DataFrame
        Combined SN/AGN/QSO/Blazar/SLSN sources with probabilities.
    """

    if 'oid' in new_df.columns:
        n = new_df.set_index('oid')
    else:
        n = new_df

    n_batches = int(np.ceil(n.shape[0] / batch_size))
    dfs = np.array_split(n, n_batches)

    stamp_class = pd.DataFrame()
    lc_class = pd.DataFrame()

    print(f"🔍 Running {n_batches} batches of up to {batch_size} OIDs each...")

    for i, batch in enumerate(dfs):
        oids_str = ",".join(["'%s'" % oid for oid in batch.index])

        # 1️⃣ Stamp classifier query
        query_stamp = f"""
            SELECT
                object.oid, object.meanra, object.meandec, object.firstmjd,
                object.ndet, probability.probability, probability.class_name,
                probability.classifier_name
            FROM 
                object 
            INNER JOIN
                probability
            ON 
                object.oid = probability.oid
            WHERE 
                object.oid IN ({oids_str})
                AND probability.classifier_name = 'stamp_classifier'
                AND probability.class_name IN ('SN','AGN')
            GROUP BY
                object.oid, object.meanra, object.meandec, object.firstmjd,
                object.ndet, probability.classifier_name,
                probability.probability, probability.class_name
            HAVING
                SUM(
                    CASE WHEN probability.class_name = 'SN' THEN probability.probability ELSE 0 END +
                    CASE WHEN probability.class_name = 'AGN' THEN probability.probability ELSE 0 END
                ) > 0.5;
        """
        sn = pd.read_sql_query(query_stamp, conn)
        stamp_class = pd.concat([stamp_class, sn], ignore_index=True)

        # 2️⃣ Light-curve classifier query
        query_lc = f"""
            SELECT
                object.oid, object.meanra, object.meandec, object.firstmjd,
                object.ndet, probability.probability, probability.class_name,
                probability.classifier_name
            FROM 
                object 
            INNER JOIN
                probability
            ON 
                object.oid = probability.oid
            WHERE 
                object.oid IN ({oids_str})
                AND probability.classifier_name = 'lc_classifier'
                AND probability.class_name IN ('AGN','QSO','Blazar','SLSN','SNII','SNIbc','SNIa')
                AND probability.ranking = 1;
        """
        sn1 = pd.read_sql_query(query_lc, conn)
        lc_class = pd.concat([lc_class, sn1], ignore_index=True)

        print(f"✅ Batch {i+1}/{n_batches}: stamp={sn.shape[0]}, lc={sn1.shape[0]}")

    # Drop duplicates
    stamp_class.drop_duplicates(subset='oid', inplace=True)
    lc_class.drop_duplicates(subset='oid', inplace=True)

    # Filter lc_class to only high-probability sources
    lc_class = lc_class.groupby('oid').filter(lambda g: g['probability'].sum() > 0.5)

    # Combine both results (avoid overlap)
    unique_to_lc = lc_class[~lc_class['oid'].isin(stamp_class['oid'])]
    candidates = pd.concat([stamp_class, unique_to_lc], ignore_index=True)

    print(f"\n🏁 Final combined sample: {candidates.shape[0]} objects.")
    return candidates
