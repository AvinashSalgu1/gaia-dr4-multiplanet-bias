"""Checks against gaiamock itself. Skipped unless gaiamock + scan files are installed."""
import numpy as np
import pytest

from gedi import gaiamock_bridge as gb

pytestmark = pytest.mark.skipif(not gb.available() or gb.n_scan_files() < 49152,
                                reason="gaiamock (compiled) and full healpix_scans/ required")

RA, DEC = 120.0, -35.0


def test_signal_matches_gaiamock_a0_predictor(monkeypatch):
    """Our reflex signal == gaiamock's, on identical epochs, with noise switched off."""
    from gedi.epochs import Epochs
    from gedi.orbits import thiele_innes, sky_offsets

    gm = gb.gaiamock()
    monkeypatch.setattr(gm, "al_uncertainty_per_ccd_interp", lambda G: 0.0 * np.asarray(G))
    a0, P, Tp, ecc, Om, inc, w = 0.8, 700.0, 100.0, 0.3, 0.7, 1.1, 2.0
    plx, pmra, pmdec = 12.0, -20.0, 7.0
    t, psi, pf, lam, _ = gm.predict_astrometry_binary_in_terms_of_a0(
        RA, DEC, plx, pmra, pmdec, P, Tp, ecc, Om, inc, w, a0, 15.0, "dr4", gb.c_funcs())
    five = pmra * t * np.sin(psi) + pmdec * t * np.cos(psi) + plx * pf
    A, B, F, G = thiele_innes(a0, inc, Om, w)
    dra, ddec = sky_offsets(t, A, B, F, G, P, ecc, 2 * np.pi * Tp / P)
    ours = dra * np.sin(psi) + ddec * np.cos(psi)
    assert np.allclose(lam - five, ours, atol=1e-9)


def test_single_planet_recovered():
    from gedi import epochs as E
    from gedi.fitting import fit_one_companion, lammers_winn_detected
    from gedi.orbits import Companion
    from gedi.simulate import Star, simulate

    rng = np.random.default_rng(3)
    ep = E.from_gaiamock(RA, DEC, "dr4", rng=rng)
    star = Star(parallax_mas=40.0, m_star_msun=0.7, G_mag=10.5, pmra=50.0, pmdec=-80.0)
    c = Companion(m_mjup=5.0, P_days=900.0, ecc=0.2, inc=1.0, Omega=0.5, omega=1.5, Tp_days=200.0)
    sim = simulate(ep, star, [c], rng=rng)
    fit = fit_one_companion(ep, sim.obs, sim.err, star.m_star_msun)
    assert fit.P_days == pytest.approx(900.0, rel=0.03)
    assert fit.m_mjup == pytest.approx(5.0, rel=0.1)
    assert lammers_winn_detected(fit)


def test_no_planet_not_detected():
    from gedi import epochs as E
    from gedi.fitting import fit_one_companion, lammers_winn_detected
    from gedi.simulate import Star, simulate

    rng = np.random.default_rng(4)
    ep = E.from_gaiamock(RA, DEC, "dr4", rng=rng)
    star = Star(parallax_mas=40.0, m_star_msun=0.7, G_mag=14.0)
    sim = simulate(ep, star, [], rng=rng)
    fit = fit_one_companion(ep, sim.obs, sim.err, star.m_star_msun)
    assert fit.delta_chi2 < 50
    assert not lammers_winn_detected(fit)


def test_f2_matches_gaiamock():
    from gedi.fitting import f2_statistic
    gm = gb.gaiamock()
    for n, red in [(60, 1.0), (80, 3.0), (150, 12.0)]:
        ours = f2_statistic(red * (n - 12), n)
        theirs = gm.predict_F2_unbinned_data(red, 12, n, Nbin=8)
        assert ours == pytest.approx(theirs, rel=1e-12)


def test_numpy_optimiser_matches_gaiamock():
    """Our orbit-only optimiser (used for the orbit+trend comparison) reaches
    gaiamock's chi2 minimum (or better) on detectable planets."""
    from gedi import epochs as E
    from gedi.fitting import fit_one_companion, fit_orbit_trend
    from gedi.orbits import Companion
    from gedi.simulate import Star, simulate

    rng = np.random.default_rng(21)
    worse = 0
    for _ in range(8):
        ep = E.from_gaiamock(rng.uniform(0, 360), rng.uniform(-60, 60), "dr4", rng=rng)
        star = Star(parallax_mas=30.0, m_star_msun=0.8, G_mag=11.0, pmra=20.0, pmdec=-10.0)
        P = 10 ** rng.uniform(2.3, 3.2)
        c = Companion(5.0, P, rng.uniform(0, 0.5), 1.0, 0.5, 1.5, rng.uniform(0, P))
        sim = simulate(ep, star, [c], rng=rng)
        gm_fit = fit_one_companion(ep, sim.obs, sim.err, star.m_star_msun)
        ours = fit_orbit_trend(ep, sim.obs, sim.err, star.m_star_msun, n_trend=0)
        worse += ours.chi2 > gm_fit.chi2 + 1.0
        assert ours.sigma_P == pytest.approx(gm_fit.sigma_P, rel=0.3)
        assert ours.sigma_a0 == pytest.approx(gm_fit.sigma_a0, rel=0.3)      # significance a0/sigma_a0
    assert worse <= 1


def test_cascade_outcomes():
    from gedi import epochs as E
    from gedi.fitting import run_cascade
    from gedi.orbits import Companion
    from gedi.simulate import Star, simulate

    rng = np.random.default_rng(22)
    ep = E.from_gaiamock(RA, DEC, "dr4", rng=rng)
    star = Star(parallax_mas=30.0, m_star_msun=0.8, G_mag=12.0)
    assert run_cascade(ep, *_oe(simulate(ep, star, [], rng=rng)))["kind"] == "none"
    out = run_cascade(ep, *_oe(simulate(ep, star, [Companion(8.0, 700.0, 0.1, 1.0, 0.3, 1.0, 50.0)], rng=rng)))
    assert out["kind"] == "orbit" and out["P"] == pytest.approx(700.0, rel=0.05)


def _oe(sim):
    return sim.obs, sim.err


def test_ruwe_formula_matches_gaiamock_cascade():
    """Our RUWE (from the 5-parameter chi2) vs gaiamock's cascade value on the same data."""
    from gedi import epochs as E
    from gedi.fitting import fit_5par, ruwe_from_chi2, run_cascade
    from gedi.orbits import Companion
    from gedi.simulate import Star, simulate

    rng = np.random.default_rng(31)
    for comps in ([], [Companion(1.0, 900.0, 0.1, 1.0, 0.3, 1.0, 50.0)]):
        ep = E.from_gaiamock(RA, DEC, "dr4", rng=rng)
        sim = simulate(ep, Star(parallax_mas=20.0, m_star_msun=0.8, G_mag=14.0), comps, rng=rng)
        res = run_cascade(ep, sim.obs, sim.err, ruwe_min=50.0)      # high threshold: always stops at 5-par
        ours = ruwe_from_chi2(fit_5par(ep, sim.obs, sim.err)[0], len(ep))
        assert ours == pytest.approx(res["ruwe"], rel=0.02)
