import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import alphashape
import psycopg2
import requests
import warnings
from astropy.time import Time

warnings.simplefilter(action='ignore', category=UserWarning)




def query_alerce_clusters(conn,skymap_df, time,ra,dec, ndays=200, alpha=0.01):
    """
    Divide the sky map into alpha-shape polygons by cluster label,
    query ALeRCE for objects inside each polygon and within [time, time+ndays].
    """
    conn = connect_alerce()
    new_df = pd.DataFrame()
    n_clusters = len(skymap_df['cluster_label'].unique())

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='mollweide')

    for i in range(n_clusters):
        cluster_data = skymap_df[skymap_df['cluster_label'] == i]
    #plt.scatter(cluster_data['meanra'],cluster_data['meandec'])
        alt = cluster_data
        alt=alt.reset_index(drop=True)
        rag=alt['meanra'].to_numpy()
        decg=alt['meandec'].to_numpy()
        combined_array3= np.concatenate([rag[:, np.newaxis], decg[:, np.newaxis]], axis=1).ravel()

        dd = combined_array3.reshape(int(combined_array3.shape[0]/2),2)
        points_2d = [(x, y) for x, y in zip(dd[:, 0], dd[:, 1])]
        alpha_shape = alphashape.alphashape(points_2d,0.01)
        
        if alpha_shape.geom_type == 'Polygon':
        # Process the single Polygon
            print("querying cluster:", i)
            x = np.array(alpha_shape.exterior.coords.xy[0])
            y = np.array(alpha_shape.exterior.coords.xy[1])
            ax.plot(x, y,'g',linewidth=1)#S,transform=ax.get_transform('world'))
            result = []
            for l in range(len(x)):
                result.append(x[l])
                result.append(y[l])

            ndays = 200
            #mjd_last = Time(datetime.utcnow(), scale='utc').mjd - ndays
            mjd_last = int(time) + ndays
            mjd_first= int(time)


            query = f"""
            SELECT
                object.oid, object.meanra, object.meandec, object.firstmjd, object.stellar,
                object.ndet
            FROM 
                object 
            WHERE q3c_poly_query(meanra, meandec,ARRAY[{','.join([str(coord) for coord in result])}])
                AND object.firstMJD >= %s
                AND object.firstMJD <= %s;;
            """%(mjd_first,mjd_last)


            try:
                results = pd.read_sql_query(query, conn)
                new_df = pd.concat([new_df, results], ignore_index=True)
                ax.scatter(new_df['meanra'], new_df['meandec'], s=1, alpha=0.1
                )#transform=ax.get_transform('world'))
            except Exception as e:
                print(f"⚠️ Query failed for cluster {i}: {e}")

        ax.scatter(ra,dec,s=10)#transform=ax.get_transform('world'))

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
