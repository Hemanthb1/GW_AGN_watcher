import numpy as np
import pandas as pd
import pytest
from astropy.table import QTable

from gw_agn_watcher.radecligo import radecligo


def test_radecligo_basic(monkeypatch, tmp_path):
    # --- Create a dummy FITS file path ---
    dummy_file = tmp_path / "dummy_skymap.fits"

    # --- Mock: download_file simply returns dummy file path ---
    monkeypatch.setattr(
        "gw_agn_watcher.radecligo.download_file", lambda url, cache=True: str(dummy_file)
    )

    # --- Mock: FITS header reading ---
    class DummyHeader:
        header = {"MJD-OBS": 60000.0}

    monkeypatch.setattr(
        "gw_agn_watcher.radecligo.fits.open", lambda path: [None, DummyHeader()]
    )

    # --- Mock: QTable.read returns minimal fake data ---
    dummy_qtable = QTable({
        "UNIQ": np.array([1, 2, 3, 4], dtype=np.int64),
        "PROBDENSITY": np.array([0.4, 0.3, 0.2, 0.1])
    })
    monkeypatch.setattr("gw_agn_watcher.radecligo.QTable.read", lambda path: dummy_qtable)

    # --- Mock: healpix conversions ---
    monkeypatch.setattr("gw_agn_watcher.radecligo.ah.uniq_to_level_ipix",
                        lambda uniq: (np.array([12, 12, 12, 12]), np.arange(len(uniq))))
    monkeypatch.setattr("gw_agn_watcher.radecligo.ah.level_to_nside", lambda lvl: 8)
    monkeypatch.setattr("gw_agn_watcher.radecligo.ah.nside_to_pixel_area", lambda ns: 1.0)
    monkeypatch.setattr("gw_agn_watcher.radecligo.ah.healpix_to_lonlat",
                        lambda ipix, nside, order="nested": (np.linspace(0,1,len(ipix))*u.rad, np.linspace(0,1,len(ipix))*u.rad))

    # --- Run the function ---
    url = "https://gracedb.ligo.org/api/superevents/S999999x/files/bayestar.fits.gz"
    skymap, df, ra, dec, mjd, event = radecligo(url, credible_level=0.9, plot=False)

    # --- Verify output types ---
    assert isinstance(skymap, QTable)
    assert isinstance(df, pd.DataFrame)
    assert isinstance(ra, np.ndarray)
    assert isinstance(dec, np.ndarray)
    assert isinstance(mjd, float)
    assert isinstance(event, str)

    # --- Check DataFrame columns ---
    for col in ["meanra", "meandec", "pixel_no", "prob_contour"]:
        assert col in df.columns

    # --- Check event name extraction ---
    assert event == "S999999x"

    # --- Probability and shape sanity checks ---
    assert len(df) <= len(skymap)
    assert np.all((ra >= 0) & (ra <= 360))
