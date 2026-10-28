"""Tests of the pure-numpy parts of scripts/extended (no gaiamock, rebound or dynesty needed)."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "extended"))


def test_linear_evidence_matches_numerical_integral():
    import g2_evidence as E
    rng = np.random.default_rng(0)
    n = 30
    X = rng.normal(size=(n, 1))
    y = 0.7 * X[:, 0] + 0.5 * rng.normal(size=n)
    err, s = np.full(n, 0.5), np.array([2.0])
    b = np.linspace(-10, 10, 200001)
    ll = -0.5 * np.sum((y[None, :] - b[:, None] * X[None, :, 0]) ** 2 / err ** 2, axis=1) - 0.5 * n * np.log(
        2 * np.pi * 0.25)
    lp = -0.5 * b ** 2 / s ** 2 - 0.5 * np.log(2 * np.pi * s ** 2)
    m = np.max(ll + lp)
    brute = m + np.log(np.sum(np.exp(ll + lp - m)) * (b[1] - b[0]))
    assert abs(E.log_evidence_linear(X, y, err, s) - brute) < 1e-6


def test_evidence_prefers_true_number_of_orbits_on_clear_signal():
    import g2_evidence as E
    from gedi import epochs as Ep
    rng = np.random.default_rng(3)
    ep = Ep.synthetic(80, rng=rng)
    nl = [700.0, 1.0, 0.1]
    X1 = E.design(ep, "one orbit", nl)
    beta = np.array([0, 0, 0, 0, 10.0, 0.8, -0.5, 0.3, 0.6])
    err = np.full(len(ep), 0.05)
    y = X1 @ beta + err * rng.standard_normal(len(ep))
    z1 = E.log_evidence_linear(X1, y, err, E.prior_widths("one orbit", 1.0))
    z0 = E.log_evidence_linear(E.design(ep, "no planet", []), y, err, E.prior_widths("no planet", 1.0))
    assert z1 > z0 + 50


def test_rv_campaign_confirms_easy_case_and_fails_hard_case():
    import g1_rv_grid as G
    rng = np.random.default_rng(1)
    inc = np.radians(60)
    easy = np.mean([G.campaign(rng, 3.0, 1.0, 9.0, 2.0, "HARPS", 2.0, [20], "extremes", inc)[0] for _ in range(30)])
    hard = np.mean([G.campaign(rng, 1.0, 5.0, 12.0, 10.0, "CORALIE", 2.0, [10], "random", inc)[0] for _ in range(30)])
    assert easy > 0.9 and hard < 0.2


def test_observable_nights_fraction_reasonable():
    import g1_rv_grid as G
    n = G.observable_nights(150.0, 2.0, np.random.default_rng(0))
    assert 0.2 * 730 < len(n) < 0.9 * 730


def test_census_expected_detections_scale_with_sensitivity():
    import g2_census as Cn
    full = np.ones((len(Cn.M_GRID), len(Cn.P_GRID)))
    none = np.zeros_like(full)
    mu_full = Cn.expected(full, 1.0)
    assert Cn.expected(none, 1.0) == 0.0
    assert 0.05 < mu_full < 0.4          # giants per star at 0.3-13 M_J, 0.3-20 yr (field model)
    assert Cn.expected(0.5 * full, 1.0) == np.float64(0.5 * mu_full) or abs(Cn.expected(0.5 * full, 1.0) - 0.5 * mu_full) < 1e-9


def test_binned_rate_recovers_input_when_response_is_identity():
    import g1_rate_binned as B
    R = np.eye(4) * 0.5
    truth = np.array([0.01, 0.02, 0.0, 0.03])
    n_star = 1e5
    nk = n_star * R @ truth
    r = B.fit_binned(nk, R, n_star, np.array([True, True, False, True]))
    assert np.allclose(r[[0, 1, 3]], truth[[0, 1, 3]], rtol=1e-3)


def test_photocentre_factor():
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "extended"))
    import numpy as np
    import g2_mixed_flags as M
    f = M.photocentre_factor(np.array([0.2, 0.8, 1.0]), 1.0)
    assert f[0] > 0.95 and 0.3 < f[1] < 0.8 and abs(f[2]) < 1e-9     # twins: photocentre does not move
    m = M.effective_mass_mj(1.0, 300.0, 500.0, 0.5)
    from gedi.orbits import reflex_a0_mas
    assert abs(reflex_a0_mas(1.0, m, 500.0, 1.0) / reflex_a0_mas(1.0, 300.0, 500.0, 1.0) - 0.5) < 1e-6
