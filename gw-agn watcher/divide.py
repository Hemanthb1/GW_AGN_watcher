"""
dividemap.py

Divide a GW skymap (RA, Dec positions) into spatial regions using K-Means clustering.
"""

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans


def dividemap(num_regions, df, plot=True, random_state=42):
    """
    Divide a LIGO/Virgo/KAGRA skymap DataFrame into clusters (spatial regions)
    using K-Means clustering on RA and Dec.

    Parameters
    ----------
    num_regions : int
        Number of regions (clusters) to divide the map into.
    df : pandas.DataFrame
        DataFrame containing at least the columns 'meanra' and 'meandec'.
    plot : bool, optional
        If True, shows a scatter plot of RA/Dec colored by cluster label.
    random_state : int, optional
        Random seed for KMeans reproducibility.

    Returns
    -------
    df_out : pandas.DataFrame
        Input DataFrame with an additional column 'cluster_label' (int).
    kmeans : sklearn.cluster.KMeans
        The fitted KMeans model.
    """

    if not {'meanra', 'meandec'}.issubset(df.columns):
        raise ValueError("DataFrame must contain columns 'meanra' and 'meandec'")

    # Fit a K-Means model
    kmeans = KMeans(n_clusters=num_regions, random_state=random_state)
    kmeans.fit(df[['meanra', 'meandec']])

    # Add cluster labels to the DataFrame
    df_out = df.copy()
    df_out['cluster_label'] = kmeans.labels_

    # Optional plot
    if plot:
        plt.figure(figsize=(6, 5))
        plt.scatter(df_out['meanra'], df_out['meandec'],
                    c=df_out['cluster_label'], cmap='tab20', s=2, alpha=0.6)
        plt.xlabel('RA [deg]')
        plt.ylabel('Dec [deg]')
        plt.title(f'K-Means Clustering into {num_regions} Regions')
        plt.show()

    return df_out, kmeans


if __name__ == "__main__":
    # Example usage
    # Generate a mock dataset for testing
    import numpy as np
    ra = np.random.uniform(0, 360, 5000)
    dec = np.random.uniform(-60, 60, 5000)
    df_mock = pd.DataFrame({'meanra': ra, 'meandec': dec})

    df_clustered, model = dividemap(10, df_mock)
    print(df_clustered.head())
