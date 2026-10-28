"""Occurrence-rate inference as a Gaia DR4 user would do it (step 3).

Two estimators, both using a completeness map eta(m, P) measured from
SINGLE-planet injection-recovery on the same stellar sample:

* binned:  f_bin = N_det,bin / (N_star * eta_bin), planets per star in a
  (mass, period) bin, where N_det,bin counts detections whose FITTED mass and
  period fall in the bin and eta_bin is the probability that a planet whose
  TRUE parameters are in the bin is detected (Fulton et al. 2021 style
  "inverse detection efficiency"; it assumes planets are log-uniform within a
  bin, as the injections are);
* Poisson likelihood: an inhomogeneous Poisson process with occurrence
  density Gamma = C (m/M_J)^alpha (P/yr)^beta per ln m per ln P inside a
  window, and expected detections N_star * sum_cells eta_cell * int_cell Gamma
  (e.g. Youdin 2011; Foreman-Mackey et al. 2014).

The same estimators applied to the TRUE planets with eta = 1 give the truth.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


# ---------------------------------------------------------------- binned
def _cell(x, edges):
    i = np.digitize(x, edges) - 1
    return np.where((i >= 0) & (i < len(edges) - 1), i, -1)


def counts_2d(m, P, m_edges, P_edges, weights=None):
    """Counts in (mass, period) cells; points outside the grid are ignored."""
    m, P = np.asarray(m, float), np.asarray(P, float)
    im, ip = _cell(m, m_edges), _cell(P, P_edges)
    ok = (im >= 0) & (ip >= 0) & np.isfinite(m) & np.isfinite(P)
    out = np.zeros((len(m_edges) - 1, len(P_edges) - 1))
    w = np.ones(len(m)) if weights is None else np.asarray(weights, float)
    np.add.at(out, (im[ok], ip[ok]), w[ok])
    return out


def completeness(m_true, P_true, det, m_edges, P_edges):
    """eta per cell = detected / injected (by TRUE parameters), and its binomial error."""
    n_inj = counts_2d(m_true, P_true, m_edges, P_edges)
    n_det = counts_2d(m_true, P_true, m_edges, P_edges, weights=np.asarray(det, float))
    with np.errstate(invalid="ignore", divide="ignore"):
        eta = n_det / n_inj
        # Wilson-like floor so that a cell with 0 detections still has an error
        err = np.sqrt(np.maximum(eta * (1 - eta), 1.0 / np.maximum(n_inj, 1)) / n_inj)
    return eta, err, n_inj


@dataclass
class BinnedRate:
    f: np.ndarray          # planets per star
    err: np.ndarray        # 1-sigma (Poisson on N_det + completeness error)
    n_det: np.ndarray


def binned_rate(m_fit, P_fit, n_star, eta, eta_err, m_edges, P_edges) -> BinnedRate:
    n = counts_2d(m_fit, P_fit, m_edges, P_edges)
    with np.errstate(invalid="ignore", divide="ignore"):
        f = n / (n_star * eta)
        rel = np.sqrt(1.0 / np.maximum(n, 1) + (eta_err / eta) ** 2)
    return BinnedRate(f, f * rel, n)


def true_rate(m, P, n_star, m_edges, P_edges):
    """Planets per star in each cell, from the injected planets."""
    return counts_2d(m, P, m_edges, P_edges) / n_star


# ---------------------------------------------------------------- Poisson likelihood
def _int_pow(k, lo, hi):
    """int_lo^hi x^k dln x."""
    k = np.asarray(k, float)
    return np.where(np.abs(k) < 1e-9, np.log(hi / lo), (hi ** k - lo ** k) / np.where(np.abs(k) < 1e-9, 1, k))


@dataclass
class PowerLawFit:
    lnC: float
    alpha: float           # mass slope (per ln m)
    beta: float            # period slope (per ln P, P in years)
    cov: np.ndarray        # covariance of (lnC, alpha, beta)
    n: int
    rate: float            # planets per star in the whole window
    rate_err: float

    def summary(self) -> str:
        e = np.sqrt(np.diag(self.cov))
        return (f"alpha {self.alpha:+.2f}+/-{e[1]:.2f}  beta {self.beta:+.2f}+/-{e[2]:.2f}  "
                f"rate in window {self.rate:.4f}+/-{self.rate_err:.4f}  (n={self.n})")


def fit_power_law(m, P_yr, n_star, m_edges, P_edges_yr, eta=None) -> PowerLawFit:
    """Maximum-likelihood power law in the window spanned by the edges.

    ``eta`` (cells x cells, same edges) is the completeness; None means eta = 1
    (the truth). Cells with unknown completeness (NaN) are treated as eta = 0
    and detections falling in them are discarded.
    """
    from scipy.optimize import minimize

    m, P_yr = np.asarray(m, float), np.asarray(P_yr, float)
    shape = (len(m_edges) - 1, len(P_edges_yr) - 1)
    eta = np.ones(shape) if eta is None else np.nan_to_num(np.asarray(eta, float), nan=0.0)
    im, ip = _cell(m, m_edges), _cell(P_yr, P_edges_yr)
    ok = (im >= 0) & (ip >= 0)
    ok[ok] &= eta[im[ok], ip[ok]] > 0
    lm, lp, n = np.log(m[ok]), np.log(P_yr[ok]), int(ok.sum())
    slm, slp = lm.sum(), lp.sum()
    m_lo, m_hi = np.asarray(m_edges[:-1]), np.asarray(m_edges[1:])
    p_lo, p_hi = np.asarray(P_edges_yr[:-1]), np.asarray(P_edges_yr[1:])

    def expected_over_C(a, b, e):
        return n_star * float(_int_pow(a, m_lo, m_hi) @ e @ _int_pow(b, p_lo, p_hi))

    def nll(x):
        lnC, a, b = x
        return -(n * lnC + a * slm + b * slp - np.exp(lnC) * expected_over_C(a, b, eta))

    def nll_profiled(ab):        # C at its maximum for given slopes
        mu = expected_over_C(*ab, eta)
        return -(n * np.log(n / mu) + ab[0] * slm + ab[1] * slp - n)

    res = minimize(nll_profiled, [0.0, 0.0], method="Nelder-Mead",
                   options=dict(xatol=1e-7, fatol=1e-9, maxiter=4000))
    a, b = res.x
    lnC = np.log(n / expected_over_C(a, b, eta))
    x0 = np.array([lnC, a, b])

    # covariance from the numerical Hessian of the full negative log-likelihood
    h = np.array([1e-4, 1e-4, 1e-4])
    H = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            ei, ej = np.eye(3)[i] * h[i], np.eye(3)[j] * h[j]
            H[i, j] = (nll(x0 + ei + ej) - nll(x0 + ei - ej) - nll(x0 - ei + ej) + nll(x0 - ei - ej)) / (4 * h[i] * h[j])
    cov = np.linalg.inv(H)

    def rate_of(x):
        return np.exp(x[0]) * expected_over_C(x[1], x[2], np.ones(shape)) / n_star

    r = rate_of(x0)
    g = np.array([(rate_of(x0 + np.eye(3)[k] * 1e-5) - rate_of(x0 - np.eye(3)[k] * 1e-5)) / 2e-5 for k in range(3)])
    return PowerLawFit(float(lnC), float(a), float(b), cov, n, float(r), float(np.sqrt(g @ cov @ g)))


# ---------------------------------------------------------------- migration-aware versions
# Near the detection threshold the FITTED mass and period of a planet scatter
# away from the true ones (and detected planets are preferentially those whose
# signal fluctuated up), so detections migrate between bins. Measuring the
# completeness by TRUE parameters but counting detections by FITTED parameters
# then biases every bin (e.g. 1-3 M_J deficits and 3-13 M_J excesses). The
# injections record both, so a DR4 user can correct for it:
#  * binned: eta_fit[k] = (injections detected with FITTED parameters in cell k)
#                          / (injections with TRUE parameters in cell k)
#    (bin-by-bin correction factor; exact if the population has the injection
#    shape inside each bin, i.e. log-uniform);
#  * power law: forward model through the response matrix
#    R[k, j] = P(detected with fitted cell k | true cell j), Poisson likelihood
#    in fitted cells; true cells extend beyond the window so planets migrating
#    in from outside are modelled.

def completeness_fitted(m_true, P_true, det, m_fit, P_fit, m_edges, P_edges, w_inj=None):
    """eta_fit per cell (detections counted at FITTED parameters, injections at TRUE ones) and its error.

    det, m_fit, P_fit are per detection (an injection that yields two reported
    planets appears twice); m_true/P_true are all injections. w_inj: optional
    weights per injection (bootstrap) and per detection as a tuple (w_all, w_det).
    """
    w_all, w_det = (None, None) if w_inj is None else w_inj
    n_inj = counts_2d(m_true, P_true, m_edges, P_edges, weights=w_all)
    n_det = counts_2d(m_fit, P_fit, m_edges, P_edges, weights=w_det)
    with np.errstate(invalid="ignore", divide="ignore"):
        eta = n_det / n_inj
        err = np.sqrt(np.maximum(n_det, 1.0)) / n_inj
    return eta, err, n_inj


def response_matrix(m_true_all, P_true_all, det_true_m, det_true_P, det_fit_m, det_fit_P,
                    m_true_edges, P_true_edges, m_fit_edges, P_fit_edges, w_all=None, w_det=None):
    """R[k, j]: expected number of reported planets in fitted cell k per planet in true cell j."""
    n_inj = counts_2d(m_true_all, P_true_all, m_true_edges, P_true_edges, weights=w_all).ravel()
    jm, jp = _cell(np.asarray(det_true_m, float), m_true_edges), _cell(np.asarray(det_true_P, float), P_true_edges)
    km, kp = _cell(np.asarray(det_fit_m, float), m_fit_edges), _cell(np.asarray(det_fit_P, float), P_fit_edges)
    ok = (jm >= 0) & (jp >= 0) & (km >= 0) & (kp >= 0)
    nPt, nPf = len(P_true_edges) - 1, len(P_fit_edges) - 1
    R = np.zeros(((len(m_fit_edges) - 1) * nPf, len(n_inj)))
    wd = np.ones(len(jm)) if w_det is None else np.asarray(w_det, float)
    np.add.at(R, (km[ok] * nPf + kp[ok], jm[ok] * nPt + jp[ok]), wd[ok])
    with np.errstate(invalid="ignore", divide="ignore"):
        R = np.where(n_inj > 0, R / n_inj, 0.0)
    return R


def fit_power_law_forward(m_fit, P_fit_yr, n_star, R, m_true_edges, P_true_edges, m_fit_edges, P_fit_edges,
                          m_win=None, P_win=None, weights=None, errors=True) -> PowerLawFit:
    """Power law Gamma = C m^alpha P^beta (per ln m per ln P) fitted to the counts in FITTED cells,
    with expected counts N_star * R @ (integral of Gamma over each TRUE cell).
    The reported rate is the integral over the window (m_win, P_win; default: the fitted-cell grid)."""
    from scipy.optimize import minimize

    n = counts_2d(m_fit, P_fit_yr, m_fit_edges, P_fit_edges, weights=weights).ravel()
    mlo, mhi = np.asarray(m_true_edges[:-1]), np.asarray(m_true_edges[1:])
    plo, phi = np.asarray(P_true_edges[:-1]), np.asarray(P_true_edges[1:])
    m_win = (m_fit_edges[0], m_fit_edges[-1]) if m_win is None else m_win
    P_win = (P_fit_edges[0], P_fit_edges[-1]) if P_win is None else P_win
    N = n.sum()

    def cells(a, b):
        return np.outer(_int_pow(a, mlo, mhi), _int_pow(b, plo, phi)).ravel()

    def mu_over_C(a, b):
        return n_star * (R @ cells(a, b))

    def nll(x):
        lnC, a, b = x
        mu = np.exp(lnC) * mu_over_C(a, b)
        return -(np.sum(n * np.log(np.maximum(mu, 1e-300))) - mu.sum())

    def nll_prof(ab):
        u = mu_over_C(*ab)
        C = N / u.sum()
        return nll([np.log(C), *ab])

    res = minimize(nll_prof, [0.0, 0.0], method="Nelder-Mead", options=dict(xatol=1e-7, fatol=1e-9, maxiter=4000))
    a, b = res.x
    x0 = np.array([np.log(N / mu_over_C(a, b).sum()), a, b])

    def rate_of(x):
        return float(np.exp(x[0]) * _int_pow(x[1], *m_win) * _int_pow(x[2], *P_win))

    if not errors:
        return PowerLawFit(float(x0[0]), float(a), float(b), np.full((3, 3), np.nan), int(round(N)), rate_of(x0),
                           np.nan)
    h = 1e-4
    H = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            ei, ej = np.eye(3)[i] * h, np.eye(3)[j] * h
            H[i, j] = (nll(x0 + ei + ej) - nll(x0 + ei - ej) - nll(x0 - ei + ej) + nll(x0 - ei - ej)) / (4 * h * h)
    cov = np.linalg.inv(H)
    r = rate_of(x0)
    g = np.array([(rate_of(x0 + np.eye(3)[k] * 1e-5) - rate_of(x0 - np.eye(3)[k] * 1e-5)) / 2e-5 for k in range(3)])
    return PowerLawFit(float(x0[0]), float(a), float(b), cov, int(round(N)), r, float(np.sqrt(max(g @ cov @ g, 0))))
