import pandas as pd
import numpy as np
import pytest

from gw_agn_watcher import main_pipeline


def test_run_pipeline(monkeypatch, tmp_path):
    # --- Mock radecligo ---
    def mock_radecligo(url):
        skymap = "mock_skymap"
        skymap1 = pd.DataFrame({
            "meanra": [10, 20],
            "meandec": [-5, 15],
            "pixel_no": [1, 2],
            "prob_contour": [0.6, 0.4]
        })
        return skymap, skymap1, np.array([10, 20]), np.array([-5, 15]), 60000.0, "S000001"
    monkeypatch.setattr(main_pipeline.radecligo, "radecligo", mock_radecligo)

    # --- Mock findminclust ---
    monkeypatch.setattr(main_pipeline.findminclust, "find_min_clusters", lambda df: 2)

    # --- Mock divide ---
    monkeypatch.setattr(main_pipeline.divide, "dividemap",
                        lambda num, df: (pd.DataFrame({"meanra": [10, 20]}), "mock_kmeans"))

    # --- Mock DB connection ---
    monkeypatch.setattr(main_pipeline, "get_alerce_connection", lambda: "mock_conn")

    # --- Mock mainquery ---
    monkeypatch.setattr(main_pipeline.mainquery, "query_alerce_clusters",
                        lambda conn, df_out, mjd, ra, dec: pd.DataFrame({"oid": ["a1", "a2"], "meanra": [10, 20]}))

    # --- Mock match_milliquas ---
    monkeypatch.setattr(main_pipeline.match_milliquas, "match_with_milliquas",
                        lambda df, agn: df.assign(matched=True))

    # --- Mock redshift ---
    monkeypatch.setattr(main_pipeline.redshift, "compute_distance_redshift",
                        lambda url: pd.DataFrame({"dist": [1]}))
    monkeypatch.setattr(main_pipeline.redshift, "filter_agn_by_redshift",
                        lambda nagn, res: pd.DataFrame({"final_2sigma": pd.Series(["a1", "a2"])}))

    # --- Mock classifiers and detections ---
    monkeypatch.setattr(main_pipeline.classifiers, "query_classifiers",
                        lambda conn, df: pd.DataFrame({"oid": ["a1", "a2"], "cls": ["AGN", "SN"]}))
    monkeypatch.setattr(main_pipeline.detections, "query_detections",
                        lambda df, conn: pd.DataFrame({"oid": ["a1", "a2"], "ndet": [5, 10]}))

    # --- Mock extinction ---
    monkeypatch.setattr(main_pipeline.extinction, "compute_lat_extinction",
                        lambda df, apply_cuts=True: (
                            pd.DataFrame({"ecl_lat": [10, 20], "gal_lat": [30, 40], "gal_A_g": [0.1, 0.2]}),
                            df
                        ))

    # --- Create dummy Milliquas CSV ---
    milliquas_path = tmp_path / "milliquas.csv"
    pd.DataFrame({"ra": [10], "dec": [-5]}).to_csv(milliquas_path, index=False)

    # --- Run the pipeline ---
    candidates, ra, dec, url = main_pipeline.run_pipeline(
        "https://mock_url", str(milliquas_path)
    )

    # --- Validate outputs ---
    assert isinstance(candidates, pd.DataFrame)
    assert "oid" in candidates.columns
    assert isinstance(ra, np.ndarray)
    assert isinstance(dec, np.ndarray)
    assert isinstance(url, str)
    assert "https://alerce.online/?" in url
    assert "&oid=" in url
