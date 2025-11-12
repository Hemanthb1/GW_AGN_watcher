"""
findclusters.py

Adaptive sky map segmentation using K-Means and alpha shapes.
Finds the minimum number of clusters (regions) such that
each region satisfies angular and geometric constraints.

Dependencies:
    pandas, numpy, matplotlib, scikit-learn, alphashape, astropy
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import alphashape
from astropy.coordinates import SkyCoord
from astropy import units as u


def polygon_diameter(polygon):
    """
    Return maximum angular separation (deg) among polygon vertices.

    Parameters
    ----------
    polygon : shapely.geometry.Polygon
        Polygon describing the region boundary.

    Returns
    -------
    float
        Maximum angular separation in degrees between any two vertices.
    """
    coords = np.array(polygon.exterior.coords)
    sky = SkyCoord(coords[:, 0], coords[:, 1], unit="deg")
    sep = sky[:, None].separation(sky[None, :])
    return sep.max().deg


def check_clusters(df, n_clusters, max_vertices=100, max_diameter=25, random_state=42):
    """
    Cluster points and check if all polygons satisfy vertex and diameter limits.

    Parameters
    ----------
    df : pandas.DataFrame
        Must contain 'meanra' and 'meandec' columns.
    n_clusters : int
        Number of clusters for K-Means.
    max_vertices : int
        Maximum allowed number of polygon vertices per cluster.
    max_diameter : float
        Maximum allowed angular diameter (deg) per cluster.
    random_state : int
        Random seed for reproducibility.

    Returns
    -------
    bool
        True if all clusters satisfy constraints, else False.
    """
    kmeans = KMeans(n_clusters=n_clusters, n_init=10, random_state=random_state)
    df['cluster_label'] = kmeans.fit_predict(df[['meanra', 'meandec']])

    for label in df['cluster_label'].unique():
        cluster_points = df[df['cluster_label'] == label][['meanra', 'meandec']].to_numpy()
        if len(cluster_points) < 3:
            continue

        try:
            alpha_shape = alphashape.alphashape(cluster_points, 0.01)
        except Exception:
            continue

        if alpha_shape.geom_type != "Polygon":
            continue

        verts = len(alpha_shape.exterior.coords)
        size = polygon_diameter(alpha_shape)

        if verts >= max_vertices or size >= max_diameter:
            return False
    return True


def find_min_clusters(df, max_vertices=100, max_diameter=25, max_try=200,
                      random_state=42, plot=False):
    """
    Find minimum n_clusters that satisfies polygon vertex and diameter constraints.

    Parameters
    ----------
    df : pandas.DataFrame
        Must contain 'meanra' and 'meandec' columns.
    max_vertices : int
        Maximum vertices allowed per alpha shape polygon.
    max_diameter : float
        Maximum angular diameter allowed per region (deg).
    max_try : int
        Maximum number of clusters to try.
    random_state : int
        Random seed for reproducibility.
    plot : bool
        If True, plots the final valid clustering.

    Returns
    -------
    int
        Minimum valid number of clusters.
    """
    for n_clusters in range(1, max_try + 1):
        if check_clusters(df.copy(), n_clusters, max_vertices, max_diameter,
                          random_state=random_state):
            if plot:
                kmeans = KMeans(n_clusters=n_clusters, n_init=10, random_state=random_state)
                df['cluster_label'] = kmeans.fit_predict(df[['meanra', 'meandec']])
                plt.scatter(df['meanra'], df['meandec'], c=df['cluster_label'],
                            cmap='tab20', alpha=0.4, s=5)
                plt.xlabel('RA [deg]')
                plt.ylabel('Dec [deg]')
                plt.title(f'Valid segmentation with n_clusters={n_clusters}')
                plt.show()
            return n_clusters

    raise RuntimeError("No valid clustering found within max_try clusters")


if __name__ == "__main__":
    # Example usage (standalone test)
    ra = np.random.uniform(0, 360, 2000)
    dec = np.random.uniform(-60, 60, 2000)
    df_mock = pd.DataFrame({'meanra': ra, 'meandec': dec})

    n = find_min_clusters(df_mock, plot=True)
    print(f"Minimum valid number of clusters: {n}")
