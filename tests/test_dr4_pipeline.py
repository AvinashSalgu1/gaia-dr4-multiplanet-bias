"""The DR4 input pipeline on REAL pre-release epoch astrometry (12 sources, June 2026)."""
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
CANDIDATES = [ROOT / "data" / "dr4_prerelease" / "GAIA_DR4_PRERELEASE_EPOCH_ASTROMETRY_RAW.xml",
              Path(__file__).parent / "data" / "GAIA_DR4_PRERELEASE_EPOCH_ASTROMETRY_RAW.xml"]
DATA = next((p for p in CANDIDATES if p.exists()), None)
pytestmark = pytest.mark.skipif(DATA is None, reason="pre-release file not present")
GAIA4, LSPM = 1457486023639239296, 1663617687609809280


def test_per_transit_reduction():
    """Transits combine the used CCDs; a single star keeps its parallax and a good fit."""
    from gaiabias import dr4
    from gaiabias.fitting import fit_5par
    d = dr4.load_ccd_table(DATA)
    assert len(dr4.source_ids(d)) == 12
    ep, y, err, n, G = dr4.source_transits(d, LSPM)
    assert len(ep) > 60 and np.median(n) >= 8
    assert np.all(err < 0.5) and np.all(np.diff(np.sort(ep.t_yr)) >= 0)
    chi2, mu = fit_5par(ep, y, err)
    assert mu[4] == pytest.approx(10.1, abs=0.1)
    assert chi2 / (len(y) - 5) == pytest.approx(1.0, abs=0.35)
    ep_x, _, err_x, _, _ = dr4.source_transits(d, LSPM, excess=True)
    assert np.all(err_x >= err - 1e-12)


def test_orbit_from_transits_matches_published():
    """Gaia-4 b from the per-transit data: P ~ 570-580 d, ~11 M_J (as from the per-CCD data)."""
    from gaiabias import dr4
    from gaiabias.fitting import fit_orbit_trend
    d = dr4.load_ccd_table(DATA)
    ep, y, err, _, _ = dr4.source_transits(d, GAIA4)
    f = fit_orbit_trend(ep, y, err, 0.64, n_trend=0, P_max=2e4)
    assert f.P_days == pytest.approx(579, rel=0.02)
    assert 9 < f.m_mjup < 13


def test_astropy_formats_read_the_same():
    """ECSV written from the VOTable gives the same CCD table (the path the bulk DR4 files will take)."""
    pytest.importorskip("astropy")
    from astropy.table import Table
    from gaiabias import dr4
    from gaiabias.dr4_epoch import read_votable_binary2
    raw = read_votable_binary2(DATA)
    keep = ["source_id", "transit_id", "centroid_pos_al", "obs_time_tcb", "obs_time_bary_corr", "scan_pos_angle",
            "parallax_factor_al", "centroid_pos_error_al", "ipd_error_al", "agis_source_excess_noise",
            "used_by_agis_al", "g_mag", "ra0", "dec0"]
    sel = [i for i, s in enumerate(raw["source_id"]) if s == LSPM]
    t = Table({k: [raw[k][i] for i in sel] for k in keep})
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / "lspm.ecsv"
        t.write(p, format="ascii.ecsv")
        a = dr4.load_ccd_table(p)
    b = dr4.load_ccd_table(DATA)
    m = b["source_id"] == LSPM
    assert len(a["w"]) == m.sum()
    assert np.allclose(np.sort(a["w"]), np.sort(b["w"][m]), equal_nan=True)


def test_fit_source_adopts_one_orbit_for_gaia4():
    from gaiabias import gaiamock_bridge as gb
    try:
        gb.gaiamock()
    except Exception:
        pytest.skip("gaiamock not installed")
    from gaiabias import dr4
    d = dr4.load_ccd_table(DATA)
    ep, y, err, _, _ = dr4.source_transits(d, GAIA4)
    r = dr4.fit_source(ep, y, err, 0.64)
    assert r["model"] == "one planet"
    assert r["n_planets"] == 1 and r["planet1_P"] == pytest.approx(579, rel=0.03)
