import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import alphashape
import psycopg2
import requests
import warnings
from astropy.time import Time

warnings.simplefilter(action='ignore', category=UserWarning)


def connect_alerce():
    """Connect to the ALeRCE PostgreSQL database."""
    url = "https://raw.githubusercontent.com/alercebroker/usecases/master/alercereaduser_v4.json"
    params = requests.get(url).json()['params']

    conn = psycopg2.connect(
        dbname=params['dbname'],
        user=params['user'],
        host=params['host'],
        password=params['password']
    )
    return conn


def query_alerce_clusters(skymap_df, time, ndays=200, alpha=0.01):
    """
    Divide the sky map into alpha-shape polygons by cluster label,
    query ALeRCE for objects inside each polygon and within [time, time+ndays].
    """
    conn = connect_alerce()
    new_df = pd.DataFrame()
    n_clusters = len(skymap_df['cluster_label'].unique())

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='astro hours mollweide')

    for i in range(n_clusters):
        cluster_data = skymap_df[skymap_df['cluster_label'] == i].reset_index(drop=True)
        ra_vals = cluster_data['meanra'].to_numpy()
        dec_vals = cluster_data['meandec'].to_numpy()

        if len(ra_vals) < 3:
            continue

        points_2d = np.column_stack((ra_vals, dec_vals))
        alpha_shape = alphashape.alphashape(points_2d, alpha)

        if alpha_shape.geom_type == 'Polygon':
            print(f"Processing cluster {i}")
            x = np.array(alpha_shape.exterior.coords.xy[0])
            y = np.array(alpha_shape.exterior.coords.xy[1])
            ax.plot(x, y, 'g', linewidth=1, transform=ax.get_transform('world'))

            coords = np.column_stack((x, y)).ravel().tolist()

            mjd_first = int(time)
            mjd_last = int(time) + ndays

            query = f"""
                SELECT
                    object.oid, object.meanra, object.meandec, object.firstmjd,
                    object.stellar, object.ndet
                FROM object
                WHERE q3c_poly_query(meanra, meandec, ARRAY[{','.join(map(str, coords))}])
                    AND object.firstmjd >= {mjd_first}
                    AND object.firstmjd <= {mjd_last};
            """

            try:
                results = pd.read_sql_query(query, conn)
                new_df = pd.concat([new_df, results], ignore_index=True)
                ax.scatter(results['meanra'], results['meandec'], s=1, alpha=0.1,
                           transform=ax.get_transform('world'))
            except Exception as e:
                print(f"⚠️ Query failed for cluster {i}: {e}")

    conn.close()
    plt.show()
    plt.close()
    return new_df


if __name__ == "__main__":
    # Example usage
    # Load your sky map DataFrame first (must include 'meanra', 'meandec', 'cluster_label')
    # Example:
    # skymap = pd.read_csv('skymap_clusters.csv')
    # time = 60250.0  # example MJD

    print("Module loaded: define your skymap and run query_alerce_clusters(skymap, time)")
