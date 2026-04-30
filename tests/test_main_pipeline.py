def test_run_pipeline(monkeypatch, tmp_path):

    def mock_radecligo(url):
        return (
            "mock_skymap",
            pd.DataFrame({
                "meanra": [10, 20],
                "meandec": [-5, 15],
                "pixel_no": [1, 2],
                "prob_contour": [0.6, 0.4]
            }),
            np.array([10, 20]),
            np.array([-5, 15]),
            60000.0,
            "S000001"
        )

    # Patch correctly (IMPORTANT)
    monkeypatch.setattr(main_pipeline, "radecligo", mock_radecligo)
    monkeypatch.setattr(main_pipeline, "get_alerce_connection", lambda: "mock_conn")

    monkeypatch.setattr(main_pipeline, "find_min_clusters", lambda df: 2)
    monkeypatch.setattr(main_pipeline, "dividemap",
                        lambda num, df: (pd.DataFrame({"meanra": [10, 20]}), "mock_kmeans"))

    monkeypatch.setattr(main_pipeline, "query_alerce_clusters",
                        lambda conn, df_out, mjd, ra, dec:
                        pd.DataFrame({"oid": ["a1", "a2"], "meanra": [10, 20]}))

    monkeypatch.setattr(main_pipeline, "match_with_milliquas",
                        lambda df, agn: df.assign(matched=True))

    monkeypatch.setattr(main_pipeline, "compute_distance_redshift",
                        lambda url: pd.DataFrame({"dist": [1]}))

    monkeypatch.setattr(main_pipeline, "filter_agn_by_redshift",
                        lambda nagn, res: pd.DataFrame({"final_2sigma": ["a1", "a2"]}))

    monkeypatch.setattr(main_pipeline, "query_classifiers",
                        lambda conn, df:
                        pd.DataFrame({"oid": ["a1", "a2"], "cls": ["AGN", "SN"]}))

    monkeypatch.setattr(main_pipeline, "query_detections",
                        lambda df, conn:
                        pd.DataFrame({"oid": ["a1", "a2"], "ndet": [5, 10]}))

    monkeypatch.setattr(main_pipeline, "compute_lat_extinction",
                        lambda df, apply_cuts=True:
                        (pd.DataFrame({"gal_A_g": [0.1, 0.2]}), df))

    # Dummy file
    path = tmp_path / "milliquas.csv"
    pd.DataFrame({"ra": [10], "dec": [-5]}).to_csv(path, index=False)

    candidates, ra, dec, url = main_pipeline.run_pipeline(
        "https://mock_url", str(path)
    )

    # Strong assertions
    assert len(candidates) == 2
    assert set(candidates["oid"]) == {"a1", "a2"}
    assert isinstance(ra, np.ndarray)
    assert isinstance(dec, np.ndarray)

    assert url.startswith("https://alerce.online/?")
    assert "oid=" in url
