"""Regression tests on REAL Gaia DR4 epoch astrometry (pre-release, 12 sources).

If a code change breaks reading or fitting real data, these fail.
"""
from pathlib import Path

import numpy as np
import pytest

DATA = Path(__file__).parent / "data" / "GAIA_DR4_PRERELEASE_EPOCH_ASTROMETRY_RAW.xml"
pytestmark = pytest.mark.skipif(not DATA.exists(), reason="pre-release file not present")


def _source(d, sid):
    from gedi.epochs import Epochs
    t_yr = d["t_ns"] / 1e9 / 86400 / 365.25 - 7.5
    m = (d["source_id"] == sid) & d["used"] & np.isfinite(d["w"]) & np.isfinite(t_yr)
    return Epochs(t_yr[m], np.radians(d["psi_deg"][m]), d["pf_al"][m]), d["w"][m], d["err"][m], d["transit"][m]


def _table():
    from gedi.dr4_epoch import ccd_table
    return ccd_table(DATA)


def test_reader_and_single_star_parallax():
    from gedi.fitting import fit_5par
    d = _table()
    assert len(np.unique(d["source_id"])) == 12
    ep, w, err, _ = _source(d, 1663617687609809280)          # LSPM J1324+6124, parallax ~10.1 mas
    chi2, mu = fit_5par(ep, w, err)
    assert mu[4] == pytest.approx(10.1, abs=0.1)
    assert chi2 / (len(w) - 5) == pytest.approx(1.0, abs=0.2)


def test_ccds_within_a_transit_are_independent():
    from gedi.fitting import _design_5par, fit_5par
    d = _table()
    z2 = []
    for sid in (4181040337841125632, 3926186255616949504, 435469040545191680, 1663617687609809280,
                2309425390592896, 20694084440761600, 2237987199365376, 60730287810150016, 10973744521070720):
        ep, w, err, tr = _source(d, sid)
        _, mu = fit_5par(ep, w, err)
        r = w - _design_5par(ep) @ mu
        for x in np.unique(tr):
            k = tr == x
            if k.sum() >= 4:
                wt = 1 / err[k] ** 2
                z2.append((np.sum(wt * r[k]) / wt.sum()) ** 2 * wt.sum())
    assert np.mean(z2) == pytest.approx(1.0, abs=0.12)       # independent CCDs -> 1


def test_known_orbits_recovered():
    from gedi.fitting import fit_orbit_trend
    d = _table()
    ep, w, err, _ = _source(d, 4318465066420528000)          # Gaia BH3
    f = fit_orbit_trend(ep, w, err, 0.76, n_trend=0, P_max=2e4)
    assert f.P_days == pytest.approx(4200, rel=0.05) and f.ecc == pytest.approx(0.73, abs=0.03)
    ep, w, err, _ = _source(d, 1457486023639239296)          # Gaia-4
    f = fit_orbit_trend(ep, w, err, 0.64, n_trend=0)
    assert f.P_days == pytest.approx(575, rel=0.03) and 8 < f.m_mjup < 14
    assert f.significance > 5
