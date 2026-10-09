"""Fitting simulated epoch astrometry with 0- and 1-companion models.

* The 5-parameter (no companion) fit and the 9-parameter linear fit at fixed
  (P, phi_p, e) are pure numpy and mirror gaiamock's C code
  (``check_ruwe`` / ``get_chi2_astrometry``), column for column.
* The search over (P, phi_p, e) uses gaiamock's compiled optimiser.

The detection test follows Lammers & Winn (2025, arXiv:2511.04673):
Delta chi^2 > 50 between the one-companion and no-companion models, best-fit
mass < 13 M_Jup, and a constrained period (P_84/P_16 < 1.5). The period
criterion is approximated here from the linearised period uncertainty; check
against an MCMC posterior for borderline cases (see docs/METHODS_AND_RESULTS.md).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .epochs import Epochs
from .orbits import campbell_from_thiele_innes, companion_mass_from_a0, orbit_XY


def _design_5par(ep: Epochs) -> np.ndarray:
    s, c = np.sin(ep.psi), np.cos(ep.psi)
    return np.vstack([s, ep.t_yr * s, c, ep.t_yr * c, ep.plx_factor]).T


def _design_9par(ep: Epochs, P_days: float, phi_p: float, ecc: float) -> np.ndarray:
    s, c = np.sin(ep.psi), np.cos(ep.psi)
    X, Y = orbit_XY(ep.t_yr, P_days, ecc, phi_p)
    # same column order as gaiamock's C code: ra, pmra, dec, pmdec, plx, B, G, A, F
    return np.vstack([s, ep.t_yr * s, c, ep.t_yr * c, ep.plx_factor,
                      X * s, Y * s, X * c, Y * c]).T


def _weighted_lstsq(M: np.ndarray, obs: np.ndarray, err: np.ndarray):
    w = 1.0 / err
    Mw, yw = M * w[:, None], obs * w
    try:   # normal equations: ~10x faster than SVD for these small, well-scaled systems
        mu = np.linalg.solve(Mw.T @ Mw, Mw.T @ yw)
    except np.linalg.LinAlgError:
        mu, *_ = np.linalg.lstsq(Mw, yw, rcond=None)
    resid = yw - Mw @ mu
    return float(resid @ resid), mu


def fit_5par(ep: Epochs, obs, err):
    """No-companion fit. Returns (chi2, [ra_off, pmra, dec_off, pmdec, plx])."""
    return _weighted_lstsq(_design_5par(ep), obs, err)


def chi2_at(ep: Epochs, obs, err, P_days: float, phi_p: float, ecc: float):
    """One-companion chi2 at fixed nonlinear parameters.

    Returns (chi2, mu) with mu = [ra_off, pmra, dec_off, pmdec, plx, B, G, A, F].
    """
    return _weighted_lstsq(_design_9par(ep, P_days, phi_p, ecc), obs, err)


@dataclass
class OneCompanionFit:
    P_days: float
    phi_p: float
    ecc: float
    A: float
    B: float
    F: float
    G: float
    parallax: float
    chi2: float
    chi2_null: float
    a0_mas: float
    inc: float
    sigma_P: float = np.nan
    sigma_a0: float = np.nan
    m_mjup: float = np.nan   # implied dark-companion mass

    @property
    def delta_chi2(self) -> float:
        return self.chi2_null - self.chi2


def fit_one_companion(ep: Epochs, obs, err, m_star_msun: float,
                      P_min: float = 10.0, P_max: float = 1e4,
                      optimizer: str = "periodogram_multistart",
                      uncertainties: bool = True) -> OneCompanionFit:
    """Best-fit single-companion orbit using gaiamock's optimiser."""
    from .gaiamock_bridge import c_funcs, gaiamock

    gm, cf = gaiamock(), c_funcs()
    t = np.ascontiguousarray(ep.t_yr, float)
    psi = np.ascontiguousarray(ep.psi, float)
    pf = np.ascontiguousarray(ep.plx_factor, float)
    obs = np.ascontiguousarray(obs, float)
    err = np.ascontiguousarray(err, float)

    P, phi_p, ecc = gm.fit_orbital_solution_nonlinear(
        t, psi, pf, obs, err, cf,
        L=np.array([P_min, 0.0, 0.0]), U=np.array([P_max, 2 * np.pi, 0.99]),
        optimizer=optimizer)
    chi2, mu = chi2_at(ep, obs, err, P, phi_p, ecc)
    chi2_null, _ = fit_5par(ep, obs, err)
    ra_off, pmra, dec_off, pmdec, plx, B, G, A, F = mu
    a0, inc, _, _ = campbell_from_thiele_innes(A, B, F, G)

    sigma_P = sigma_a0 = np.nan
    if uncertainties:
        p0 = [ra_off, dec_off, plx, pmra, pmdec, P, ecc, phi_p, A, B, F, G]
        errors, a0_gm, sigma_a0, _ = gm.get_uncertainties_at_best_fit_binary_solution(
            t, psi, pf, obs, err, p0, cf, binned=True)
        sigma_P = errors[5]

    m = companion_mass_from_a0(a0, P, m_star_msun, plx) if plx > 0 else np.nan
    return OneCompanionFit(P, phi_p, ecc, A, B, F, G, plx, chi2, chi2_null,
                           a0, inc, sigma_P, sigma_a0, m)


def lammers_winn_detected(fit: OneCompanionFit, dchi2_min: float = 50.0,
                          m_max_mjup: float = 13.0, P_ratio_max: float = 1.5) -> bool:
    """Lammers & Winn (2026) detection criteria (period test approximated)."""
    if not (fit.delta_chi2 > dchi2_min and fit.m_mjup < m_max_mjup):
        return False
    if not np.isfinite(fit.sigma_P):
        return False
    lo = fit.P_days - fit.sigma_P
    return lo > 0 and (fit.P_days + fit.sigma_P) / lo < P_ratio_max


def f2_statistic(chi2: float, n_obs: int, n_param: int = 12, n_bin: int = 8) -> float:
    """Gaia-style F2 goodness of fit for a fit to transit-binned data.

    Same as gaiamock.predict_F2_unbinned_data: the binned reduced chi2 is
    converted to its unbinned (per-CCD) equivalent, then Wilson-Hilferty
    transformed. F2 ~ N(0, 1) for a good fit; Gaia DR3 orbital solutions
    required F2 < 25.
    """
    chi2_red_binned = chi2 / (n_obs - n_param)
    nu_unbinned = n_obs * n_bin - n_param
    chi2_red_unbinned = (n_obs * n_bin - n_obs + chi2_red_binned * (n_obs - n_param)) / nu_unbinned
    return float(np.sqrt(9 * nu_unbinned / 2) * (chi2_red_unbinned ** (1 / 3) + 2 / (9 * nu_unbinned) - 1))


def ruwe_from_chi2(chi2_5par: float, n_obs: int, n_bin: int = 8) -> float:
    """RUWE a star would get from a single-star (5-parameter) fit.

    Same per-CCD conversion as f2_statistic: chi2 of the 5-parameter fit to the
    transit-binned data -> unbinned chi2 (each transit = n_bin CCD measurements,
    whose scatter about their own mean adds n_bin - 1 per transit), then the
    unit-weight error sqrt(chi2_unbinned / (n_bin N - 5)). Equal to 1 on average
    for a star with nothing but noise (the "renormalised" in RUWE is then 1).
    """
    return float(np.sqrt((n_obs * (n_bin - 1) + chi2_5par) / (n_obs * n_bin - 5)))


# ---------------------------------------------------------------------------
# One companion + polynomial trend (acceleration / jerk), pure numpy.
#
# Gaia's DR3 pipeline (and gaiamock's cascade) fits EITHER an acceleration
# solution OR an orbit, never both. A planet with a long-period companion is
# better described by orbit + trend; these functions fit that model. The
# trend terms are linear, so at fixed (P, phi_p, e) the whole model is a linear
# least-squares problem, exactly like gaiamock's 9-parameter solve.
# ---------------------------------------------------------------------------

def _trend_columns(ep: Epochs, n_trend: int) -> list:
    """n_trend = 0 (none), 2 (acceleration), 4 (acceleration + jerk); same
    parametrisation as gaiamock's 7- and 9-parameter models."""
    s, c, t = np.sin(ep.psi), np.cos(ep.psi), ep.t_yr
    cols = []
    if n_trend >= 2:
        cols += [0.5 * t**2 * s, 0.5 * t**2 * c]
    if n_trend >= 4:
        cols += [t**3 / 6 * s, t**3 / 6 * c]
    return cols


def fit_trend_only(ep: Epochs, obs, err, n_trend: int):
    """No-companion model with an acceleration (n_trend=2) or jerk (4) term."""
    M = np.column_stack([_design_5par(ep)] + _trend_columns(ep, n_trend)) if n_trend else _design_5par(ep)
    return _weighted_lstsq(M, obs, err)


def _design_orbit_trend(ep, P, phi_p, ecc, n_trend):
    base = _design_9par(ep, P, phi_p, ecc)
    extra = _trend_columns(ep, n_trend)
    return np.column_stack([base] + extra) if extra else base


def chi2_at_trend(ep, obs, err, P, phi_p, ecc, n_trend):
    return _weighted_lstsq(_design_orbit_trend(ep, P, phi_p, ecc, n_trend), obs, err)


@dataclass
class OrbitTrendFit:
    n_trend: int
    P_days: float
    phi_p: float
    ecc: float
    mu: np.ndarray          # ra, pmra, dec, pmdec, plx, B, G, A, F, [trend terms]
    chi2: float
    chi2_null: float        # trend-only model with the same n_trend
    a0_mas: float
    inc: float
    sigma_P: float
    m_mjup: float
    sigma_a0: float = np.nan     # uncertainty of a0 (linearised, inflated like sigma_P)

    @property
    def significance(self) -> float:
        """a0 / sigma_a0, the Gaia NSS 'significance' of an orbit (DR3 definition)."""
        return self.a0_mas / self.sigma_a0 if self.sigma_a0 > 0 else np.nan

    @property
    def delta_chi2(self) -> float:
        """Evidence for the companion over the matching no-companion model."""
        return self.chi2_null - self.chi2

    @property
    def parallax(self) -> float:
        return float(self.mu[4])


def _a0_sigma(bgaf, cov4):
    """Uncertainty of the photocentre semi-major axis a0 from the covariance of the
    Thiele-Innes constants, given in design-column order (B, G, A, F)."""
    B, G, A, F = bgaf

    def a0(x):
        b, g, a, f = x
        u = 0.5 * (a * a + b * b + f * f + g * g)
        v = a * g - b * f
        return np.sqrt(u + np.sqrt(max(u * u - v * v, 0.0)))

    x0 = np.array([B, G, A, F], float)
    h = 1e-6 * max(np.abs(x0).max(), 1e-6)
    grad = np.array([(a0(x0 + h * np.eye(4)[k]) - a0(x0 - h * np.eye(4)[k])) / (2 * h) for k in range(4)])
    var = float(grad @ cov4 @ grad)
    return np.sqrt(var) if var > 0 else np.inf


def _linearised_sigma_P(ep, obs, err, P, phi_p, ecc, mu, n_trend):
    return _linearised_errors(ep, obs, err, P, phi_p, ecc, mu, n_trend)[0]


def _linearised_errors(ep, obs, err, P, phi_p, ecc, mu, n_trend):
    """sigma_P from the Jacobian of the normalised residuals at the best fit,
    inflated by the goodness of fit exactly as in gaiamock's
    get_uncertainties_at_best_fit_binary_solution."""
    def model(P_, ph_, e_):
        return _design_orbit_trend(ep, P_, ph_, e_, n_trend) @ mu

    J_lin = _design_orbit_trend(ep, P, phi_p, ecc, n_trend) / err[:, None]
    cols = []
    for k, (val, step) in enumerate([(P, 1e-6 * P), (phi_p, 1e-6), (ecc, 1e-6)]):
        up = [P, phi_p, ecc]
        dn = [P, phi_p, ecc]
        up[k] = val + step
        dn[k] = max(val - step, 0.0) if k == 2 else val - step
        cols.append((model(*up) - model(*dn)) / (up[k] - dn[k]) / err)
    J = np.column_stack([J_lin] + cols)
    try:
        cov = np.linalg.inv(J.T @ J)
    except np.linalg.LinAlgError:
        return np.inf, np.inf
    n, k = len(obs), J.shape[1]
    resid = (obs - model(P, phi_p, ecc)) / err
    chi2_red_binned = float(resid @ resid) / (n - k)
    nu_unb = n * 8 - k
    chi2_red_unb = (n * 8 - n + chi2_red_binned * (n - k)) / nu_unb
    cc = np.sqrt(chi2_red_unb / (1 - 2 / (9 * nu_unb)) ** 3)
    var = cov[J_lin.shape[1], J_lin.shape[1]]
    sigma_P = float(cc * np.sqrt(var)) if var > 0 else np.inf
    sigma_a0 = cc * _a0_sigma(mu[5:9], cov[5:9, 5:9])
    return sigma_P, float(sigma_a0)


def _circular_periodogram(ep, obs, err, freqs, n_trend, extra=None):
    """chi2 of circular orbit (+ trend) at every trial frequency, batched.

    `extra` (n x k) adds fixed columns to every trial model, e.g. a first
    planet's orbit when searching for a second signal.
    """
    w = 1.0 / err
    cols = [_design_5par(ep)] + _trend_columns(ep, n_trend) + ([extra] if extra is not None else [])
    base = np.column_stack(cols) * w[:, None]
    s, c = np.sin(ep.psi) * w, np.cos(ep.psi) * w
    arg = 2 * np.pi * np.outer(freqs, ep.t_yr * 365.25)            # (Nf, n)
    X, Y = np.cos(arg), np.sin(arg)
    orb = np.stack([X * s, Y * s, X * c, Y * c], axis=-1)            # (Nf, n, 4)
    M = np.concatenate([np.broadcast_to(base, (len(freqs),) + base.shape), orb], axis=-1)
    yw = obs * w
    A = np.einsum("fni,fnj->fij", M, M)
    b = np.einsum("fni,n->fi", M, yw)
    mu = np.linalg.solve(A, b[..., None])[..., 0]
    resid = yw[None, :] - np.einsum("fni,fi->fn", M, mu)
    return np.sum(resid**2, axis=1)


def fit_orbit_trend(ep: Epochs, obs, err, m_star_msun: float, n_trend: int = 2,
                    P_min: float = 10.0, P_max: float = 1e4, oversample: float = 5.0,
                    n_peaks: int = 4) -> OrbitTrendFit:
    """One companion + polynomial trend, fitted with a numpy optimiser.

    1. Circular-orbit periodogram on a uniform frequency grid (for e = 0 the
       orbital phase is absorbed by the linear Thiele-Innes terms, so each
       trial period is an exact linear solve).
    2. Nelder-Mead polish of (log P, phi_p, e) from the best few peaks and
       several starting phases.
    With n_trend = 0 this is the same model as gaiamock's orbit fit, which is
    how the optimiser is validated (see tests/).
    """
    from scipy.optimize import minimize

    obs = np.asarray(obs, float)
    err = np.asarray(err, float)
    span = (ep.t_yr.max() - ep.t_yr.min()) * 365.25
    df = 1.0 / (oversample * span)
    freqs = np.arange(1.0 / P_max, 1.0 / P_min, df)
    chis = _circular_periodogram(ep, obs, err, freqs, n_trend)

    # distinct local minima, best first
    order = np.argsort(chis)
    peaks = []
    for i in order:
        if all(abs(i - j) > 2 for j in peaks):
            peaks.append(i)
        if len(peaks) >= n_peaks:
            break

    lnPmin, lnPmax = np.log(P_min), np.log(P_max)

    def objective(x):
        lnP, ph, e = x
        if not (lnPmin <= lnP <= lnPmax) or not (0.0 <= e < 0.99):
            return 1e30
        return chi2_at_trend(ep, obs, err, np.exp(lnP), ph % (2 * np.pi), e, n_trend)[0]

    best = (np.inf, None)
    for i in peaks:
        for ph0 in (0.0, np.pi):
            for e0 in (0.1, 0.5):
                r = minimize(objective, [np.log(1.0 / freqs[i]), ph0, e0], method="Nelder-Mead",
                             options=dict(xatol=1e-5, fatol=1e-4, maxiter=600))
                if r.fun < best[0]:
                    best = (r.fun, r.x)

    lnP, ph, e = best[1]
    P, ph = float(np.exp(lnP)), float(ph % (2 * np.pi))
    chi2, mu = chi2_at_trend(ep, obs, err, P, ph, e, n_trend)
    chi2_null, _ = fit_trend_only(ep, obs, err, n_trend)
    B, G, A, F = mu[5:9]
    a0, inc, _, _ = campbell_from_thiele_innes(A, B, F, G)
    sigma_P, sigma_a0 = _linearised_errors(ep, obs, err, P, ph, e, mu, n_trend)
    plx = mu[4]
    m = companion_mass_from_a0(a0, P, m_star_msun, plx) if plx > 0 else np.nan
    return OrbitTrendFit(n_trend, P, ph, float(e), mu, chi2, chi2_null, a0, inc, sigma_P, m, sigma_a0)


def detected_lw(delta_chi2: float, m_mjup: float, P_days: float, sigma_P: float,
                dchi2_min: float = 50.0) -> bool:
    """Lammers & Winn detection criteria on plain numbers (any model)."""
    if not (delta_chi2 > dchi2_min and np.isfinite(m_mjup) and m_mjup < 13.0):
        return False
    if not np.isfinite(sigma_P):
        return False
    lo = P_days - sigma_P
    return lo > 0 and (P_days + sigma_P) / lo < 1.5


def run_cascade(ep: Epochs, obs, err, ruwe_min: float = 1.4) -> dict:
    """What a DR3-style Gaia pipeline would publish (gaiamock's cascade).

    The cascade tries, in order: a 5-parameter single-star solution (accepted
    if RUWE < ruwe_min), then 9- and 7-parameter ACCELERATION solutions
    (accepted if they fit well and are significant - in which case no orbit
    is ever tried), and only then a one-companion orbit.

    Returns dict(kind=..., and orbit parameters if kind == "orbit"), where
    kind is one of: "insufficient", "none" (5-par), "accel9", "accel7", "orbit".
    """
    from .gaiamock_bridge import c_funcs, gaiamock

    gm = gaiamock()
    res = gm.fit_full_astrometric_cascade(
        np.ascontiguousarray(ep.t_yr, float), np.ascontiguousarray(ep.psi, float),
        np.ascontiguousarray(ep.plx_factor, float), np.ascontiguousarray(obs, float),
        np.ascontiguousarray(err, float), c_funcs(), ruwe_min=ruwe_min, skip_acceleration=False)
    head = res[0]
    if head == 0 and res[1] == 0:
        return dict(kind="insufficient")
    if head == -1:
        return dict(kind="none", ruwe=res[1])
    if head == -9:
        return dict(kind="accel9")
    if head == -7:
        return dict(kind="accel7")
    return dict(kind="orbit", parallax=res[0], P=res[10], sigma_P=res[11], phi_p=res[12],
                ecc=res[14], a0=res[17], sigma_a0=res[18], F2=res[21], ruwe=res[22])


# ----------------------------------------------------------------------------- two planets
def _design_two(ep, P1, ph1, e1, P2, ph2, e2):
    """5-parameter motion + two Keplerian orbits (13 linear columns)."""
    s, c = np.sin(ep.psi), np.cos(ep.psi)
    X, Y = orbit_XY(ep.t_yr, P2, e2, ph2)
    return np.column_stack([_design_9par(ep, P1, ph1, e1), X * s, Y * s, X * c, Y * c])


def _linearised_sigma_two(ep, obs, err, nl, mu):
    """sigma_P of both planets from the Jacobian at the best fit, inflated for
    goodness of fit as in _linearised_sigma_P. Returns (sigma_P1, sigma_P2)."""
    def model(x):
        return _design_two(ep, *x) @ mu

    J_lin = _design_two(ep, *nl) / err[:, None]
    cols = []
    for k, val in enumerate(nl):
        step = 1e-6 * val if k in (0, 3) else 1e-6
        up, dn = list(nl), list(nl)
        up[k] = val + step
        dn[k] = max(val - step, 0.0) if k in (2, 5) else val - step
        cols.append((model(up) - model(dn)) / (up[k] - dn[k]) / err)
    J = np.column_stack([J_lin] + cols)
    try:
        cov = np.linalg.inv(J.T @ J)
    except np.linalg.LinAlgError:
        return np.inf, np.inf, np.inf, np.inf
    n, k = len(obs), J.shape[1]
    resid = (obs - model(nl)) / err
    chi2_red_binned = float(resid @ resid) / max(n - k, 1)
    nu_unb = n * 8 - k
    chi2_red_unb = (n * 8 - n + chi2_red_binned * (n - k)) / nu_unb
    cc = np.sqrt(chi2_red_unb / (1 - 2 / (9 * nu_unb)) ** 3)
    out = []
    for idx in (J_lin.shape[1], J_lin.shape[1] + 3):
        v = cov[idx, idx]
        out.append(float(cc * np.sqrt(v)) if v > 0 else np.inf)
    sa = [float(cc * _a0_sigma(mu[c0:c0 + 4], cov[c0:c0 + 4, c0:c0 + 4])) for c0 in (5, 9)]
    return out[0], out[1], sa[0], sa[1]


@dataclass
class TwoPlanetFit:
    """Sequential two-planet fit: planet 1 = best single orbit, planet 2 = best
    extra signal found in what planet 1 leaves behind, then a joint refit."""
    one: OrbitTrendFit              # the one-planet fit (starting point)
    P: tuple                        # (P1, P2) days
    phi_p: tuple
    ecc: tuple
    mu: np.ndarray                  # 5 motion + 4 + 4 Thiele-Innes terms
    chi2: float                     # joint two-planet chi2
    a0: tuple
    m_mjup: tuple
    sigma_P: tuple
    drop_dchi2: tuple               # chi2 increase when each planet is removed (others re-solved)
    sigma_a0: tuple = (np.nan, np.nan)

    @property
    def significance(self) -> tuple:
        """a0 / sigma_a0 for each planet (Gaia NSS 'significance', DR3 definition)."""
        return tuple(a / s if s > 0 else np.nan for a, s in zip(self.a0, self.sigma_a0))

    @property
    def gain_over_one(self) -> float:
        """How much the second planet improves on the one-planet fit (delta chi2)."""
        return self.one.chi2 - self.chi2


def fit_two_planets(ep: Epochs, obs, err, m_star_msun: float, one: OrbitTrendFit | None = None,
                    P_min: float = 10.0, P_max: float = 1e4, oversample: float = 5.0,
                    n_peaks: int = 4) -> TwoPlanetFit:
    """Search for a second planet after modelling the first, then fit both jointly.

    1. Fit one planet (``fit_orbit_trend`` with no trend) unless given.
    2. Circular-orbit periodogram of the data with planet 1 held at its best
       fit (planet 1's four linear terms re-solved at every trial period).
    3. Nelder-Mead on all six orbital parameters (P, phi_p, e for each planet)
       from the best few periodogram peaks; the linear terms (5 motion + 4 + 4)
       are solved exactly each time.
    The joint solution is returned whether or not the second planet is
    significant; decide with ``gain_over_one`` and ``drop_dchi2``.
    """
    from scipy.optimize import minimize

    obs, err = np.asarray(obs, float), np.asarray(err, float)
    if one is None:
        one = fit_orbit_trend(ep, obs, err, m_star_msun, n_trend=0)
    first_cols = _design_9par(ep, one.P_days, one.phi_p, one.ecc)[:, 5:]

    span = (ep.t_yr.max() - ep.t_yr.min()) * 365.25
    freqs = np.arange(1.0 / P_max, 1.0 / P_min, 1.0 / (oversample * span))
    chis = _circular_periodogram(ep, obs, err, freqs, 0, extra=first_cols)
    peaks = []
    for i in np.argsort(chis):
        if all(abs(i - j) > 2 for j in peaks):
            peaks.append(i)
        if len(peaks) >= n_peaks:
            break

    lo, hi = np.log(P_min), np.log(P_max)

    def objective(x):
        l1, p1, e1, l2, p2, e2 = x
        if not (lo <= l1 <= hi and lo <= l2 <= hi) or not (0 <= e1 < 0.99 and 0 <= e2 < 0.99):
            return 1e30
        return _weighted_lstsq(_design_two(ep, np.exp(l1), p1 % (2 * np.pi), e1,
                                           np.exp(l2), p2 % (2 * np.pi), e2), obs, err)[0]

    best = (np.inf, None)
    for i in peaks:
        for ph0 in (0.0, np.pi):
            x0 = [np.log(one.P_days), one.phi_p, one.ecc, np.log(1.0 / freqs[i]), ph0, 0.1]
            r = minimize(objective, x0, method="Nelder-Mead",
                         options=dict(xatol=1e-5, fatol=1e-4, maxiter=1500, maxfev=1800))
            if r.fun < best[0]:
                best = (r.fun, r.x)
    x = best[1]
    nl = [float(np.exp(x[0])), float(x[1] % (2 * np.pi)), float(x[2]),
          float(np.exp(x[3])), float(x[4] % (2 * np.pi)), float(x[5])]
    M = _design_two(ep, *nl)
    chi2, mu = _weighted_lstsq(M, obs, err)
    a0s, ms, drops = [], [], []
    for k, (c0, P) in enumerate(((5, nl[0]), (9, nl[3]))):
        B, G, A, F = mu[c0:c0 + 4]
        a0, _, _, _ = campbell_from_thiele_innes(A, B, F, G)
        a0s.append(a0)
        ms.append(companion_mass_from_a0(a0, P, m_star_msun, mu[4]) if mu[4] > 0 else np.nan)
        keep = [j for j in range(13) if not (c0 <= j < c0 + 4)]
        drops.append(_weighted_lstsq(M[:, keep], obs, err)[0] - chi2)
    sP1, sP2, sa1, sa2 = _linearised_sigma_two(ep, obs, err, nl, mu)
    sP = (sP1, sP2)
    return TwoPlanetFit(one, (nl[0], nl[3]), (nl[1], nl[4]), (nl[2], nl[5]), mu, chi2,
                        tuple(a0s), tuple(ms), tuple(sP), tuple(drops), (sa1, sa2))
