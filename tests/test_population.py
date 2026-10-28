"""Step 3 building blocks: occurrence model, mock stars, planet draws, rate estimators."""
import numpy as np
import pytest

from gedi import population as P
from gedi import rates as R
from gedi.orbits import hill_stable


def test_occurrence_matches_cls():
    # L&W Eq. 23 with C = 0.15 per ln a per ln m over 0.3-13 M_J -> ~15 giants
    # per 100 stars at 2-8 au (CLS: 14.1 +2.0/-1.8, Fulton et al. 2021)
    assert P.rate_in_box(0.3, 13, 2, 8) == pytest.approx(0.154, abs=0.002)
    assert P.lambda_total(0.9) == pytest.approx(0.367, abs=0.002)
    assert P.lambda_total(1.8) == pytest.approx(2 * P.lambda_total(0.9))


def test_a_sampler_follows_shape():
    rng = np.random.default_rng(0)
    a = P.draw_a_au(rng, 400_000)
    frac = np.mean((a > 2) & (a < 8))
    assert frac == pytest.approx(P.rate_in_box(0.3, 13, 2, 8) / P.lambda_total(0.9), rel=0.02)


def test_clustered_host_fraction_reproduces_bryan():
    h = P.host_fraction_clustered()
    lam_box = P.rate_in_box(1, 13, 5, 20)
    assert 1 - np.exp(-lam_box / h) == pytest.approx(P.BRYAN_COMP_FRAC, rel=1e-9)
    assert 0.05 < h < 0.3


def test_stars_respect_limits():
    st = P.draw_stars(20_000, np.random.default_rng(1))
    d = 1000 / st.parallax
    assert np.all(st.G <= P.G_MAX + 1e-9)
    assert np.all(d <= P.r_max_pc(st.m_star) * (1 + 1e-9))
    assert np.all((st.m_star >= 0.3) & (st.m_star <= 1.4))
    assert np.all((st.dec >= -90) & (st.dec <= 90))


@pytest.mark.parametrize("scenario", ["independent", "clustered"])
def test_planets_same_mean_rate_and_stable(scenario):
    st = P.draw_stars(30_000, np.random.default_rng(2))
    pl = P.draw_planets(st, np.random.default_rng(3), scenario)
    expected = P.lambda_total(st.m_star).sum()
    # same mean occurrence in both scenarios (minus the few planets dropped as unstable)
    assert (len(pl) + pl.n_dropped) == pytest.approx(expected, rel=0.08)
    for s in np.nonzero(np.bincount(pl.star) >= 2)[0][:300]:
        ii = np.nonzero(pl.star == s)[0]
        for i in ii:
            for j in ii:
                if i < j:
                    assert hill_stable(st.m_star[s], pl.companion(i), pl.companion(j))


def test_clustered_has_more_multis():
    st = P.draw_stars(30_000, np.random.default_rng(4))
    frac = {}
    for sc in ("independent", "clustered"):
        n = np.bincount(P.draw_planets(st, np.random.default_rng(5), sc).star)
        frac[sc] = (n >= 2).sum() / (n >= 1).sum()
    assert frac["clustered"] > 3 * frac["independent"]


def _power_law_sample(rng, n_star, C, alpha, beta, m_edges, p_edges):
    """Draw planets from Gamma = C m^alpha P^beta per ln m per ln P (P in yr)."""
    def draw(k, lo, hi, size):
        u = rng.uniform(0, 1, size)
        if abs(k) < 1e-9:
            return lo * (hi / lo) ** u
        return (lo**k + u * (hi**k - lo**k)) ** (1 / k)
    n_exp = n_star * C * float(R._int_pow(alpha, m_edges[0], m_edges[-1])) * float(R._int_pow(beta, p_edges[0], p_edges[-1]))
    n = rng.poisson(n_exp)
    return draw(alpha, m_edges[0], m_edges[-1], n), draw(beta, p_edges[0], p_edges[-1], n)


def test_power_law_fit_recovers_truth_with_completeness():
    rng = np.random.default_rng(6)
    me, pe = np.geomspace(0.3, 13, 9), np.geomspace(0.1, 8, 13)
    m, Pyr = _power_law_sample(rng, 1_000_000, 0.02, -0.3, 0.4, me, pe)
    # a completeness that depends on both m and P, applied by thinning
    eta_of = lambda m_, p_: np.clip(0.2 * (m_ / 1.0) ** 0.8 * (p_ / 2.0) ** -0.5, 0, 1)   # noqa: E731
    keep = rng.uniform(0, 1, len(m)) < eta_of(m, Pyr)
    # eta per cell from injections, as in the real analysis
    mi = np.exp(rng.uniform(np.log(0.3), np.log(13), 400_000))
    pi = np.exp(rng.uniform(np.log(0.1), np.log(8), 400_000))
    eta, _, _ = R.completeness(mi, pi, rng.uniform(0, 1, len(mi)) < eta_of(mi, pi), me, pe)
    f = R.fit_power_law(m[keep], Pyr[keep], 1_000_000, me, pe, eta)
    t = R.fit_power_law(m, Pyr, 1_000_000, me, pe)
    # with ~1e5 detections the statistical errors are tiny; what remains is the
    # finite-cell completeness approximation, which must stay at the ~1% level
    assert abs(f.alpha + 0.3) < 0.03
    assert abs(f.beta - 0.4) < 0.03
    assert f.rate == pytest.approx(t.rate, rel=0.03)
    assert t.rate == pytest.approx(len(m) / 1_000_000, rel=1e-3)


def test_binned_rate_is_unbiased():
    rng = np.random.default_rng(7)
    me, pe = np.array([0.3, 1, 3, 13]), np.array([0.1, 0.3, 1, 2, 4, 8])
    m, Pyr = _power_law_sample(rng, 2_000_000, 0.01, 0.0, 0.0, me, pe)   # log-uniform
    det = rng.uniform(0, 1, len(m)) < 0.3
    mi = np.exp(rng.uniform(np.log(0.3), np.log(13), 300_000))
    pi = np.exp(rng.uniform(np.log(0.1), np.log(8), 300_000))
    eta, eerr, _ = R.completeness(mi, pi, rng.uniform(0, 1, len(mi)) < 0.3, me, pe)
    r = R.binned_rate(m[det], Pyr[det], 2_000_000, eta, eerr, me, pe)
    truth = R.true_rate(m, Pyr, 2_000_000, me, pe)
    pull = (r.f - truth) / r.err
    assert np.all(np.abs(pull) < 4)


def test_deep_sample_reference_planet_at_threshold():
    from gedi.orbits import reflex_a0_mas
    from gedi.simulate import sigma_al_per_transit
    for m in (0.35, 0.8, 1.3):
        r = float(P.r_ref_pc(m)[0])
        G = P.abs_G(m) + 5 * np.log10(r / 10)
        snr = reflex_a0_mas(m, P.M_REF_MJ, P.T_DR4_DAYS, 1000 / r) / sigma_al_per_transit(G)
        assert snr == pytest.approx(P.SNR_REF, rel=0.03)
    # deep sample is much closer than the L&W one
    assert np.all(P.r_max_pc([0.3, 1.0, 1.4], "deep") < 0.3 * P.r_max_pc([0.3, 1.0, 1.4], "lw"))


def test_gaia_selection_cuts():
    """L&W main-sequence and distance cuts, 250 pc, RUWE and mass assignment on a toy table."""
    from gedi import gaia_stars as GS
    # star 0: Sun-like at 50 pc, kept; 1: too red for the main sequence; 2: RUWE 2; 3: M_G 13 (too faint)
    plx = np.array([20.0, 20.0, 20.0, 20.0])
    MG = np.array([4.7, 4.7, 4.7, 13.0])
    raw = dict(source_id=np.arange(4), ra=np.zeros(4), dec=np.zeros(4), parallax=plx, parallax_error=plx / 50,
               pmra=np.zeros(4), pmdec=np.zeros(4), phot_g_mean_mag=MG - 5 * np.log10(plx / 100),
               bp_rp=np.array([0.8, 3.0, 0.8, 3.0]), ruwe=np.array([1.0, 1.0, 2.0, 1.0]),
               mass_flame=np.array([np.nan, np.nan, np.nan, np.nan]))
    sel, counts = GS.select(raw, ruwe_max=1.4)
    assert list(sel["source_id"]) == [0]
    assert sel["mass"][0] == pytest.approx(1.0, abs=0.05)
    assert sel["distance_pc"][0] == pytest.approx(50.0)
    sel2, _ = GS.select(raw, ruwe_max=None)
    assert list(sel2["source_id"]) == [0, 2]
    assert GS.mass_from_MG(np.array([11.5]))[0] == pytest.approx(0.23)


def test_lw_distance_limit_peaks_at_MG_2p8():
    """L&W Eqs. 36-37: 1600 pc at M_G = 2.8, closer for brighter AND for fainter stars."""
    r = P.r_max_lw(np.array([1.0, 2.8, 4.7, 10.0]))
    assert r[1] == pytest.approx(1600.0)
    assert r[0] < r[1] and r[2] < r[1] and r[3] < r[2]
    assert r[3] == pytest.approx(1600.0 * 10 ** (-0.093 * 7.2))      # ~340 pc for an early M dwarf


def test_gaia_sky_tiles_cover_the_sky_by_source_id():
    """Whole-sky layout: 768 HEALPix level-3 tiles as source_id ranges that tile the id space,
    match the HEALPix pixel of the position, and give a random-order subset for any sky fraction."""
    from gedi import gaia_stars as GS
    ranges = [GS.tile_source_id_range(p) for p in range(GS.N_TILES)]
    assert ranges[0][0] == 0 and all(a[1] == b[0] for a, b in zip(ranges, ranges[1:]))
    assert ranges[-1][1] == 12 * 4 ** 12 * 2 ** 35            # every level-12 pixel, nothing beyond
    assert GS.n_tiles_for(1.0) == 768 and GS.n_tiles_for(0.25) == 192
    with pytest.raises(SystemExit):
        GS.n_tiles_for(1.5)
    assert sorted(GS.tile_order()) == list(range(GS.N_TILES)) and np.array_equal(GS.tile_order(), GS.tile_order())
    assert GS.merged_path("tiles") != GS.merged_path("circles")
    for m in GS.MIRRORS.values():
        q = m["tile_query"].format(lo="0", hi=str(2 ** 53), plx_min="4.000000")
        assert ">= 0" in q and f"< {2 ** 53}" in q and "4.000000" in q and "{" not in q and "CIRCLE" not in q
        assert q.lstrip().upper().startswith("SELECT")
    hp = pytest.importorskip("healpy")
    rng = np.random.default_rng(1)
    theta, phi = np.arccos(rng.uniform(-1, 1, 500)), rng.uniform(0, 2 * np.pi, 500)
    sid = hp.ang2pix(4096, theta, phi, nest=True).astype(np.int64) << 35      # Gaia's source_id encoding
    assert np.array_equal(sid >> GS.TILE_SHIFT, hp.ang2pix(8, theta, phi, nest=True))


def test_gaia_download_queries_and_columns():
    """Both archives' queries format, and archive answers map to the Gaia column names whatever the case."""
    from gedi import gaia_stars as GS
    for m in GS.MIRRORS.values():
        q = m["query"].format(ra0="12.5", dec0="-33.1", r="3", plx_min="4.000000")
        assert "CIRCLE('ICRS', 12.5, -33.1, 3)" in q and "4.000000" in q and "{" not in q
    assert '"BP-RP" AS bp_rp' in GS.VIZIER_QUERY and 'FROM "I/355/gaiadr3"' in GS.VIZIER_QUERY
    # sky circles: never overlapping, fixed random order, sky fraction -> number of circles
    ra, dec = GS.cone_centres()
    v = np.c_[np.cos(np.radians(dec)) * np.cos(np.radians(ra)), np.cos(np.radians(dec)) * np.sin(np.radians(ra)),
              np.sin(np.radians(dec))]
    sep = np.degrees(np.arccos(np.clip(v @ v.T, -1, 1)))
    np.fill_diagonal(sep, 180)
    assert sep.min() > 2 * GS.CONE_R_DEG
    assert np.array_equal(GS.cone_order(), GS.cone_order()) and sorted(GS.cone_order()) == list(range(GS.N_CONES))
    assert GS.n_cones_for(0.25) == 365
    with pytest.raises(SystemExit):
        GS.n_cones_for(0.6)

    class Tab(dict):
        colnames = property(lambda self: list(self))
    tab = Tab({c.upper(): np.arange(3) * 1.0 for c in GS.COLS})
    arr = GS._table_to_arrays(tab)
    assert set(arr) == set(GS.COLS) and arr["source_id"].dtype == np.int64
    with pytest.raises(KeyError):
        GS._table_to_arrays(Tab({"ra": np.zeros(2)}))


def _toy_survey(rng, m, P):
    """Detect planets on a noisy signal; fitted m and P scatter (Eddington-like migration)."""
    snr_true = m * P ** (2 / 3)
    snr_obs = snr_true * np.exp(0.35 * rng.standard_normal(len(m)))
    det = snr_obs > 3.0
    m_fit = m * snr_obs / snr_true
    P_fit = P * np.exp(0.15 * rng.standard_normal(len(m)))
    return det, m_fit, P_fit


def test_migration_aware_rates_close():
    """Counting detections at FITTED values with completeness at TRUE values is biased; the fitted-space
    completeness and the forward-modelled power law recover the truth."""
    from gedi import rates as R
    rng = np.random.default_rng(5)
    mE, pE = np.array([0.3, 1.0, 3.0, 13.0]), np.array([0.1, 0.3, 1.0, 2.0, 4.0, 8.0])
    mT, pT = np.geomspace(0.3, 13, 7), np.r_[0.05, np.geomspace(0.1, 8, 9), 16.0]
    # injections: log-uniform
    ni = 400_000
    mi, Pi = np.exp(rng.uniform(np.log(0.3), np.log(13), ni)), np.exp(rng.uniform(np.log(0.05), np.log(16), ni))
    di, mfi, Pfi = _toy_survey(rng, mi, Pi)
    # population: alpha = 0, beta = 0.4 over 0.05-16 yr, 0.2 planets per star in total
    n_star, beta = 300_000, 0.4
    npl = rng.poisson(0.2 * n_star)
    m = np.exp(rng.uniform(np.log(0.3), np.log(13), npl))
    u = rng.uniform(size=npl)
    P = (0.05 ** beta + u * (16 ** beta - 0.05 ** beta)) ** (1 / beta)
    d, mf, Pf = _toy_survey(rng, m, P)
    truth = R.true_rate(m, P, n_star, mE, pE)
    eta_t, eta_te, _ = R.completeness(mi, Pi, di, mE, pE)
    eta_f, eta_fe, _ = R.completeness_fitted(mi, Pi, di, mfi[di], Pfi[di], mE, pE)
    naive = R.binned_rate(mf[d], Pf[d], n_star, eta_t, eta_te, mE, pE).f
    fixed = R.binned_rate(mf[d], Pf[d], n_star, eta_f, eta_fe, mE, pE).f
    ok = (eta_t > 0.05) & (truth > 0)
    assert np.max(np.abs(naive[ok] / truth[ok] - 1)) > 0.15          # the old way is visibly off
    assert np.max(np.abs(fixed[ok] / truth[ok] - 1)) < 0.2            # the new way is close in every bin
    assert fixed[ok].sum() / truth[ok].sum() == pytest.approx(1.0, abs=0.03)
    Rm = R.response_matrix(mi, Pi, mi[di], Pi[di], mfi[di], Pfi[di], mT, pT, mT, pT[1:-1])
    f = R.fit_power_law_forward(mf[d], Pf[d], n_star, Rm, mT, pT, mT, pT[1:-1])
    t = R.fit_power_law(m, P, n_star, mT, pT[1:-1])
    assert abs(f.alpha - 0.0) < 3 * np.sqrt(f.cov[1, 1]) + 0.05
    assert abs(f.beta - beta) < 3 * np.sqrt(f.cov[2, 2]) + 0.05
    assert f.rate == pytest.approx(t.rate, rel=0.08)


def test_multiplicity_reweighting():
    """Zero-truncated Poisson multiplicity: f <-> mu, and host weights reproduce the requested f."""
    from gedi import multiplicity as MU
    for f in (0.0, 0.12, 0.59, 0.8):
        assert MU.f_multi(MU.mu_for(f)) == pytest.approx(f, abs=1e-6) if f > 0 else MU.mu_for(f) < 1e-6
    rng = np.random.default_rng(0)
    k = np.r_[np.ones(3000, int), rng.integers(2, 6, 3000)]
    pool = MU.Pool(k=k, **{c: np.zeros(0) for c in ("t_host", "t_m", "t_P", "a_host", "a_m", "a_P", "a_cls",
                                                    "b_host", "b_m", "b_P", "b_cls")}, m_star=np.ones(len(k)))
    w, mu, lost, f_real = MU.host_weights_by_mass(pool, 0.4, np.zeros(len(k), int))
    assert w.sum() == pytest.approx(1.0)
    assert np.sum(w[k >= 2]) == pytest.approx(f_real, abs=1e-9)
    assert abs(f_real - 0.4) < 0.02 and lost < 0.02


def test_rv_helpers():
    """Jupiter's K on the Sun, V from G, precision and epoch rules behave as documented."""
    from gedi import rv as RV
    assert RV.semi_amplitude(1.0, 11.86, 1.0) == pytest.approx(12.5, rel=0.02)    # Jupiter: 12.5 m/s
    assert RV.v_mag(10.0, 0.82) == pytest.approx(10.17, abs=0.03)                # Sun-like colour: V - G ~ 0.17
    s8, s12 = RV.rv_precision(8.0, 1.0), RV.rv_precision(12.0, 1.0)
    assert 2.0 < s8 < s12                                                        # jitter floor, fainter is worse
    assert RV.rv_precision(8.0, 1.0, "FEROS") > s8
    assert RV.epochs_known_period(1000.0, 3.0) == 6                              # huge signal: minimum epochs
    assert RV.epochs_known_period(15.0, 3.0) == 6 and RV.epochs_known_period(3.0, 3.0) == 50
    # a slow companion is easier to catch over a longer programme
    assert RV.epochs_companion(30.0, 20.0, 3.0, 3.0) < RV.epochs_companion(30.0, 20.0, 3.0, 1.0)


def test_population_variants():
    """Variants change only multi-planet systems: near-coplanar orbits, low eccentricities; unknown names fail;
    the default draws are unchanged by the refactor (same rng sequence)."""
    from gedi.orbits import mutual_inclination
    st = P.draw_stars(8_000, np.random.default_rng(6))
    base = P.draw_planets(st, np.random.default_rng(7), "clustered")
    again = P.draw_planets(st, np.random.default_rng(7), "clustered", variant=())
    assert np.array_equal(base.P, again.P) and np.array_equal(base.inc, again.inc)
    cop = P.draw_planets(st, np.random.default_rng(7), "clustered", variant=("coplanar",))
    low = P.draw_planets(st, np.random.default_rng(7), "clustered", variant=("lowecc",))
    for pl, check in ((cop, "inc"), (low, "e")):
        n = np.bincount(pl.star)
        multi = n[pl.star] >= 2
        assert multi.sum() > 100
        if check == "e":
            assert np.median(pl.e[multi]) < 0.15 and np.median(pl.e[~multi]) > 0.15
        else:
            mi = []
            for s in np.nonzero(n >= 2)[0][:300]:
                ii = np.nonzero(pl.star == s)[0]
                mi.append(np.degrees(mutual_inclination(pl.inc[ii[0]], pl.Omega[ii[0]], pl.inc[ii[1]], pl.Omega[ii[1]])))
            assert np.median(mi) < 5 and max(mi) < 15
    with pytest.raises(ValueError):
        P.draw_planets(st, np.random.default_rng(7), "clustered", variant=("nonsense",))


def test_nbody_stability_screen():
    """A widely separated pair survives; a pair inside each other's Hill sphere does not."""
    pytest.importorskip("rebound")
    wide = [dict(m=1.0, P=365.25, e=0.0, inc=0.1, Omega=0.0, omega=0.0, Tp=0.0),
            dict(m=1.0, P=365.25 * 8, e=0.0, inc=0.1, Omega=0.0, omega=0.0, Tp=100.0)]
    tight = [dict(m=10.0, P=365.25, e=0.3, inc=0.1, Omega=0.0, omega=0.0, Tp=0.0),
             dict(m=10.0, P=365.25 * 1.15, e=0.3, inc=0.1, Omega=0.0, omega=3.0, Tp=100.0)]
    assert P.nbody_stable(1.0, wide, n_orbits=2000)
    assert not P.nbody_stable(1.0, tight, n_orbits=2000)
