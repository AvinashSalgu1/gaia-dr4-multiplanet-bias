"""Pure-numpy checks of the orbit and fitting code (no gaiamock needed)."""
import numpy as np
import pytest

from gaiabias import epochs as E
from gaiabias.fitting import chi2_at, fit_5par
from gaiabias.orbits import (Companion, campbell_from_thiele_innes, companion_mass_from_a0,
                             reflex_a0_mas, reflex_along_scan, relative_semimajor_axis_au,
                             sky_offsets, solve_kepler, thiele_innes)
from gaiabias.simulate import Star, sigma_al_per_transit, simulate

RNG = np.random.default_rng(1)


@pytest.mark.parametrize("ecc", [0.0, 0.3, 0.7, 0.95])
def test_kepler_solver(ecc):
    M = np.linspace(-10, 10, 2001)
    Ea = solve_kepler(M, ecc)
    assert np.allclose(Ea - ecc * np.sin(Ea), np.mod(M, 2 * np.pi), atol=1e-10)


def test_jupiter_sun():
    a = relative_semimajor_axis_au(4332.59, 1.0 + 9.547919e-4)
    assert a == pytest.approx(5.2034, rel=2e-3)
    # Sun's reflex for Jupiter seen from 10 pc: ~0.497 mas
    a0 = reflex_a0_mas(1.0, 1.0, 4332.59, parallax_mas=100.0)
    assert a0 == pytest.approx(0.4966, rel=3e-3)


def test_mass_inversion():
    a0 = reflex_a0_mas(0.8, 3.3, 900.0, 25.0)
    assert companion_mass_from_a0(a0, 900.0, 0.8, 25.0) == pytest.approx(3.3, rel=1e-8)


def test_thiele_innes_roundtrip():
    for _ in range(500):
        a0 = RNG.uniform(0.01, 5)
        inc = np.arccos(RNG.uniform(-1, 1))
        Om, om = RNG.uniform(0, np.pi), RNG.uniform(0, 2 * np.pi)
        A, B, F, G = thiele_innes(a0, inc, Om, om)
        a0r, incr, Omr, omr = campbell_from_thiele_innes(A, B, F, G)
        assert a0r == pytest.approx(a0, rel=1e-9)
        assert incr == pytest.approx(inc, abs=1e-6)
        # cos i = (AG - BF)/a0^2 (Halbwachs et al. 2023; used by gaiamock)
        assert np.cos(incr) == pytest.approx((A * G - B * F) / a0**2, abs=1e-9)
        # (Omega, omega) recovered up to the usual simultaneous pi ambiguity
        assert np.allclose(thiele_innes(a0r, incr, Omr, omr), (A, B, F, G), atol=1e-9 * a0)


def test_face_on_circular_has_constant_separation():
    A, B, F, G = thiele_innes(1.7, 0.0, 0.4, 1.1)
    t = np.linspace(-3, 3, 400)
    dra, ddec = sky_offsets(t, A, B, F, G, P_days=700.0, ecc=0.0, phi_p=0.3)
    assert np.allclose(np.hypot(dra, ddec), 1.7)


def test_two_companions_superpose():
    ep = E.synthetic(rng=RNG)
    m_star, plx = 0.9, 20.0
    c1 = Companion(2.0, 400.0, 0.1, 0.7, 0.3, 1.0, 50.0)
    c2 = Companion(5.0, 1500.0, 0.3, 1.2, 2.0, 4.0, 300.0)
    both, a0s = reflex_along_scan(ep.t_yr, ep.psi, [c2, c1], m_star, plx)  # order should not matter
    s1, _ = reflex_along_scan(ep.t_yr, ep.psi, [c1], m_star, plx)
    # outer companion orbits star + inner companion (Jacobi)
    m_in = 2.0 * 9.547919e-4
    a0_2 = reflex_a0_mas(m_star, 5.0, 1500.0, plx, m_inner_msun=m_in)
    A, B, F, G = thiele_innes(a0_2, c2.inc, c2.Omega, c2.omega)
    dra, ddec = sky_offsets(ep.t_yr, A, B, F, G, c2.P_days, c2.ecc, c2.phi_p)
    s2 = dra * np.sin(ep.psi) + ddec * np.cos(ep.psi)
    assert np.allclose(both, s1 + s2)
    assert a0s[1] == pytest.approx(a0_2)


def test_noise_free_single_companion_is_fit_exactly():
    """Simulator and the 9-parameter linear model share conventions (gaiamock's)."""
    ep = E.synthetic(rng=RNG)
    star = Star(parallax_mas=15.0, m_star_msun=1.1, G_mag=11.0, pmra=-30.0, pmdec=12.0,
                ra_off=0.2, dec_off=-0.1)
    c = Companion(4.0, 800.0, 0.35, 1.0, 0.8, 2.5, 120.0)
    sim = simulate(ep, star, [c], add_noise=False)
    chi2, mu = chi2_at(ep, sim.obs, sim.err, c.P_days, c.phi_p, c.ecc)
    assert chi2 < 1e-12
    ra_off, pmra, dec_off, pmdec, plx, B, G, A, F = mu
    assert [ra_off, pmra, dec_off, pmdec, plx] == pytest.approx(
        [star.ra_off, star.pmra, star.dec_off, star.pmdec, star.parallax_mas], abs=1e-9)
    assert np.allclose([A, B, F, G], thiele_innes(sim.a0_mas[0], c.inc, c.Omega, c.omega), atol=1e-9)


def test_null_chi2_distribution():
    ep = E.synthetic(rng=RNG)
    star = Star(parallax_mas=10.0, m_star_msun=1.0, G_mag=14.0)
    chis = [fit_5par(ep, *_obs_err(simulate(ep, star, [], rng=RNG)))[0] for _ in range(400)]
    assert np.mean(chis) == pytest.approx(len(ep) - 5, rel=0.06)


def test_sigma_al():
    assert sigma_al_per_transit(12.0) == pytest.approx(0.125 / np.sqrt(8))


def _obs_err(sim):
    return sim.obs, sim.err


def test_hill_stability():
    from gaiabias.orbits import hill_stable
    # L&W Sec. 6.2: 5 M_J at 3 au around 0.8 Msun needs an outer 5 M_J at >~5.3 au
    P = lambda a_au: 365.25 * np.sqrt(a_au**3 / 0.8)
    inner = Companion(5.0, P(3.0))
    assert not hill_stable(0.8, inner, Companion(5.0, P(5.0)))
    assert hill_stable(0.8, inner, Companion(5.0, P(5.6)))
    # eccentricity makes it harder
    assert not hill_stable(0.8, inner, Companion(5.0, P(5.6), ecc=0.3))


def test_f2_good_fit_near_zero():
    from gaiabias.fitting import f2_statistic
    n = 80
    chis = RNG.chisquare(n - 12, 4000)
    f2 = np.array([f2_statistic(c, n) for c in chis])
    assert abs(np.median(f2)) < 0.5
    # Binned reduced chi2 of 5 is diluted over 8x more per-CCD dof: F2 ~ 7, NOT flagged
    assert 5 < f2_statistic(5.0 * (n - 12), n) < 9
    assert f2_statistic(25.0 * (n - 12), n) > 25


def test_orbit_plus_acceleration_recovered():
    """Noise-free orbit + acceleration: orbit+trend fit recovers P; orbit-only can't fit it."""
    from gaiabias.fitting import fit_orbit_trend
    rng = np.random.default_rng(5)
    ep = E.synthetic(n_transits=80, rng=rng)
    star = Star(parallax_mas=25.0, m_star_msun=0.9, G_mag=12.0, pmra=10.0, pmdec=-5.0)
    c = Companion(4.0, 700.0, 0.25, 0.9, 1.0, 2.0, 100.0)
    sim = simulate(ep, star, [c], add_noise=False)
    acc = 0.5 * ep.t_yr**2 * (0.08 * np.sin(ep.psi) - 0.05 * np.cos(ep.psi))   # mas/yr^2 terms
    obs = sim.obs + acc
    with_trend = fit_orbit_trend(ep, obs, sim.err, star.m_star_msun, n_trend=2)
    no_trend = fit_orbit_trend(ep, obs, sim.err, star.m_star_msun, n_trend=0)
    assert with_trend.P_days == pytest.approx(700.0, rel=0.01)
    assert with_trend.chi2 < 1e-3
    assert list(with_trend.mu[9:11]) == pytest.approx([0.08, -0.05], abs=1e-4)
    assert no_trend.chi2 > 100 * max(with_trend.chi2, 1e-6)


def test_orbit_only_optimiser_finds_true_minimum():
    """With noise, the numpy optimiser should reach at least the chi2 of the true orbit."""
    from gaiabias.fitting import chi2_at, fit_orbit_trend
    rng = np.random.default_rng(8)
    worse = 0
    for _ in range(8):
        ep = E.synthetic(n_transits=80, rng=rng)
        star = Star(parallax_mas=30.0, m_star_msun=0.8, G_mag=11.0)
        P = 10 ** rng.uniform(2.2, 3.2)
        c = Companion(3.0, P, rng.uniform(0, 0.5), 1.0, 0.5, 1.5, rng.uniform(0, P))
        sim = simulate(ep, star, [c], rng=rng, bright_star_jitter=False)
        f = fit_orbit_trend(ep, sim.obs, sim.err, star.m_star_msun, n_trend=0)
        chi2_true, _ = chi2_at(ep, sim.obs, sim.err, c.P_days, c.phi_p, c.ecc)
        worse += f.chi2 > chi2_true + 1.0
        assert np.isfinite(f.sigma_P) and f.sigma_P > 0
    assert worse == 0


def test_mutual_inclination_geometry():
    from gaiabias.orbits import mutual_inclination, orientation_relative_to
    rng = np.random.default_rng(3)
    for _ in range(50):
        i0, O0 = np.arccos(rng.uniform(-1, 1)), rng.uniform(0, 2 * np.pi)
        d, az = rng.uniform(0, np.pi), rng.uniform(0, 2 * np.pi)
        i1, O1 = orientation_relative_to(i0, O0, d, az)
        assert mutual_inclination(i0, O0, i1, O1) == pytest.approx(d, abs=1e-9)
    i1, O1 = orientation_relative_to(1.0, 0.7, 0.0, 2.0)          # coplanar: same plane
    assert [i1, O1] == pytest.approx([1.0, 0.7], abs=1e-12)
    assert mutual_inclination(0.3, 1.0, np.pi - 0.3, 1.0) == pytest.approx(np.pi - 0.6)


def test_Tp_for_mean_anomaly():
    from gaiabias.orbits import Tp_for_mean_anomaly
    P, t, M = 900.0, -56.875, 1.2
    Tp = Tp_for_mean_anomaly(M, t, P)
    assert (2 * np.pi * (t - Tp) / P) == pytest.approx(M)


def test_two_planet_fit_recovers_both():
    from gaiabias.fitting import fit_two_planets
    rng = np.random.default_rng(5)
    ep = E.synthetic(n_transits=90, rng=rng)
    star = Star(parallax_mas=40.0, m_star_msun=0.9, G_mag=11.0)
    t = Companion(8.0, 600.0, 0.15, 1.0, 0.5, 1.5, 100.0)
    c = Companion(10.0, 1500.0, 0.2, 0.7, 2.0, 0.8, 300.0)
    sim = simulate(ep, star, [t, c], rng=rng, bright_star_jitter=False)
    f = fit_two_planets(ep, sim.obs, sim.err, star.m_star_msun)
    Ps = sorted(f.P)
    assert Ps[0] == pytest.approx(600.0, rel=0.05)
    assert Ps[1] == pytest.approx(1500.0, rel=0.08)
    assert f.gain_over_one > 100 and min(f.drop_dchi2) > 100
    assert f.chi2 < f.one.chi2


def test_two_planet_fit_single_planet_gives_small_gain():
    from gaiabias.fitting import fit_two_planets
    rng = np.random.default_rng(6)
    ep = E.synthetic(n_transits=90, rng=rng)
    star = Star(parallax_mas=40.0, m_star_msun=0.9, G_mag=11.0)
    sim = simulate(ep, star, [Companion(8.0, 600.0, 0.15, 1.0, 0.5, 1.5, 100.0)], rng=rng,
                   bright_star_jitter=False)
    f = fit_two_planets(ep, sim.obs, sim.err, star.m_star_msun)
    assert f.gain_over_one < 40          # only noise left to fit (7 free parameters)
    assert f.P[0] == pytest.approx(600.0, rel=0.03) or f.P[1] == pytest.approx(600.0, rel=0.03)


def test_step2_outcome_rules():
    """Target / companion / wrong classification and the two-planet 'if justified' rule."""
    from gaiabias.step2.outcomes import judge
    nan = np.nan
    # system 0: target found correctly; 1: wrong (mass 50% off); 2: companion found; 3: nothing
    d = dict(P_t=np.array([500.0] * 4), m_t=np.array([3.0] * 4), P_c=np.array([3000.0] * 4),
             m_c=np.array([5.0] * 4),
             one_w_det=np.array([1, 1, 1, 0.0]), one_w_P=np.array([505, 500, 3050, nan]),
             one_w_m=np.array([3.1, 4.5, 5.2, nan]),
             one_a_det=np.array([1, 1, 1, 1.0]), one_a_P=np.array([500.0] * 4), one_a_m=np.array([3.0] * 4))
    o = judge(d, "one")
    assert list(o["target"]) == [True, False, False, False]
    assert list(o["wrong"]) == [False, True, False, False]
    assert list(o["companion"]) == [False, False, True, False]
    assert list(o["none"]) == [False, False, False, True]
    # two-planet: system 0 second planet significant -> both reported; system 1 not -> one-planet answer
    two = dict(d, **{f"two_{v}_{q}": x for v in "wa" for q, x in [
        ("P1", np.array([500.0, 500, 500, 500])), ("m1", np.array([3.0] * 4)), ("d1", np.array([300.0] * 4)),
        ("s1", np.array([5.0] * 4)), ("P2", np.array([3000.0, 3000, 3000, 3000])),
        ("m2", np.array([5.0] * 4)), ("d2", np.array([200.0, 10, 10, 10]) if v == "w" else np.array([1.0] * 4)),
        ("s2", np.array([100.0] * 4))]},
        **{f"accel_{v}_{q}": x for v in "wa" for q, x in [
            ("o_det", np.array([1, 1, 0, 0.0])), ("o_P", np.array([500.0, 500, nan, nan])),
            ("o_m", np.array([3.0, 4.5, nan, nan]) if v == "w" else np.array([3.0, 3.0, nan, nan]))]})
    o2 = judge(two, "two")
    assert list(o2["n_reported"]) == [2, 1, 0, 0]
    assert o2["target"][0] and o2["companion"][0]
    assert o2["wrong"][1] and not o2["target"][1]
    assert not o2["false_alarm"].any()


def test_bic_choice():
    from gaiabias.step2.outcomes import bic_choice
    n = np.array([80.0, 80.0, 80.0])          # ln 80 = 4.38: each extra parameter costs 4.38 in chi2
    d = {"n_transits": n,
         "accel_w_c5": np.array([1000.0, 100.0, 400.0]), "accel_w_c7": np.array([990.0, 95.0, 300.0]),
         "accel_w_c0": np.array([90.0, 90.0, 200.0]), "accel_w_c2": np.array([88.0, 89.0, 150.0]),
         "two_w_chi2": np.array([85.0, 88.0, 60.0])}
    # 0: one planet (others gain too little); 1: no planet; 2: two planets (jerk fits absent -> ignored)
    assert list(bic_choice(d, "w")) == [3, 0, 6]
    d.update(jerk_w_c9=np.array([985.0, 94.0, 300.0]), jerk_w_cJ=np.array([30.0, 89.0, 150.0]))
    assert list(bic_choice(d, "w")) == [5, 0, 6]        # a big jerk improvement wins for system 0


def test_ruwe_from_chi2_is_one_for_noise():
    from gaiabias.fitting import fit_5par, ruwe_from_chi2
    rng = np.random.default_rng(9)
    vals = []
    for _ in range(200):
        ep = E.synthetic(n_transits=80, rng=rng)
        sim = simulate(ep, Star(parallax_mas=20.0, m_star_msun=0.8, G_mag=14.0), [], rng=rng)
        vals.append(ruwe_from_chi2(fit_5par(ep, sim.obs, sim.err)[0], len(ep)))
    assert np.mean(vals) == pytest.approx(1.0, abs=0.01)
    assert np.std(vals) < 0.02


def test_significance_matches_monte_carlo():
    """Linearised a0 error vs the scatter of a0 over noise (robust width; a few blended fits excluded)."""
    from gaiabias.fitting import fit_orbit_trend
    rng = np.random.default_rng(4)
    ep = E.synthetic(n_transits=80, rng=rng)
    star = Star(parallax_mas=30.0, m_star_msun=0.9, G_mag=12.0)
    c = Companion(3.0, 700.0, 0.2, 1.0, 0.5, 1.5, 100.0)
    a0, s = [], []
    for _ in range(40):
        sim = simulate(ep, star, [c], rng=rng, bright_star_jitter=False)
        f = fit_orbit_trend(ep, sim.obs, sim.err, 0.9, n_trend=0)
        if abs(np.log(f.P_days / 700.0)) < 0.05:
            a0.append(f.a0_mas)
            s.append(f.sigma_a0)
    robust = 1.4826 * np.median(np.abs(np.array(a0) - np.median(a0)))
    assert np.median(s) == pytest.approx(robust, rel=0.35)
