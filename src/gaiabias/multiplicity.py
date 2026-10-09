"""Multiplicity scan for step 3, by reweighting systems that were already fitted.

Each star's Gaia result depends only on that star and its own planets, so a
population with ANY giant-planet multiplicity can be assembled from the host
stars simulated in the two step-3 scenarios (independent + clustered pooled),
without new fits:

* the number of giants k on a host follows a zero-truncated Poisson
  distribution ZTP(k; mu); mu is set by the target fraction f of giant hosts
  that have two or more giants, f = 1 - mu e^-mu / (1 - e^-mu). The
  independent scenario has f ~ 0.12, the clustered (Bryan et al. 2016) one
  f ~ 0.59;
* host stars with k giants are taken from the pool of simulated hosts with k
  giants, each weighted by ZTP(k; mu) / N_k(pool);
* the fraction of stars that host giants is then set so that the number of
  giants per star INSIDE the detectable window (0.3-13 M_J, 0.1-8 yr) equals
  the CLS-calibrated value of the occurrence model (Lammers & Winn Eq. 23):
  every multiplicity level has exactly the same true occurrence rate there.
  This is done separately in bins of stellar mass, so the occurrence keeps its
  stellar-mass dependence (more giants around more massive stars); otherwise
  changing the multiplicity would also change which stars host the planets,
  and with it their detectability.
  (This also removes the slightly lower in-window rate of the original
  clustered draw, caused by Hill-stability redraws in crowded systems.)

Errors: bootstrap with Poisson(1) weights on the hosts and on the injections.

Assumptions: planets on a k-giant host are distributed as in the simulation
(the occurrence model's period-mass distribution, Hill stable); stars without
giants never yield a detection (the pre-fit cut check found none below the
cut); stellar and brown-dwarf companions are not included.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def f_multi(mu):
    """Fraction of hosts with >= 2 giants for a zero-truncated Poisson with mean parameter mu."""
    mu = np.asarray(mu, float)
    return 1 - mu * np.exp(-mu) / -np.expm1(-mu)


def mu_for(f):
    """Inverse of f_multi (bisection); f = 0 gives mu -> 0 (every host has exactly one giant)."""
    if f <= 1e-9:
        return 1e-9
    lo, hi = 1e-9, 50.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if f_multi(mid) < f else (lo, mid)
    return 0.5 * (lo + hi)


def ztp_pmf(k, mu):
    from scipy.stats import poisson
    k = np.asarray(k)
    return np.where(k >= 1, poisson.pmf(k, mu) / -np.expm1(-mu), 0.0)


@dataclass
class Pool:
    """All simulated host stars (both scenarios), with their true planets and reported ones."""
    k: np.ndarray            # giants per host
    t_host: np.ndarray       # true planets
    t_m: np.ndarray
    t_P: np.ndarray
    a_host: np.ndarray       # reported planets, all planets on the star (A)
    a_m: np.ndarray
    a_P: np.ndarray
    a_cls: np.ndarray
    b_host: np.ndarray       # reported planets, each planet alone (B)
    b_m: np.ndarray
    b_P: np.ndarray
    b_cls: np.ndarray
    m_star: np.ndarray       # host star mass


def build_pool(scen_results: dict) -> Pool:
    """scen_results: {scenario: hosts dict from step3.analyse_one}."""
    parts = {k: [] for k in Pool.__dataclass_fields__}
    off = 0
    for scen, h in scen_results.items():
        npl, pst = h["npl"], h["pst"]
        hosts = np.nonzero(npl > 0)[0]
        hid = -np.ones(len(npl), int)
        hid[hosts] = off + np.arange(len(hosts))
        parts["k"].append(npl[hosts])
        parts["m_star"].append(h["m_star"][hosts])
        parts["t_host"].append(hid[pst]), parts["t_m"].append(h["pm"]), parts["t_P"].append(h["pP"])
        parts["a_host"].append(hid[h["A_star"]]), parts["a_m"].append(h["A_m"]), parts["a_P"].append(h["A_P"])
        parts["a_cls"].append(h["A_cls"])
        parts["b_host"].append(hid[pst[h["B_pl"]]]), parts["b_m"].append(h["B_m"]), parts["b_P"].append(h["B_P"])
        parts["b_cls"].append(h["B_cls"])
        off += len(hosts)
    return Pool(**{k: np.concatenate(v) for k, v in parts.items()})


def host_weights(pool: Pool, f: float, k_min_pool: int = 20):
    """Relative weight per host (sums to 1) for multiplicity f, and the probability mass lost
    to k values with fewer than k_min_pool hosts in the pool."""
    mu = mu_for(f)
    kmax = int(pool.k.max())
    N_k = np.bincount(pool.k, minlength=kmax + 1).astype(float)
    q = ztp_pmf(np.arange(kmax + 1), mu)
    usable = N_k >= k_min_pool
    lost = float(q[~usable].sum() + max(0.0, 1 - q.sum()))
    q = np.where(usable, q, 0.0)
    q /= q.sum()
    with np.errstate(invalid="ignore", divide="ignore"):
        per = np.where(N_k > 0, q / N_k, 0.0)
    return per[pool.k], mu, lost


def host_weights_by_mass(pool: Pool, f: float, mbin_host, k_min_pool: int = 10):
    """Like host_weights, but the multiplicity distribution is imposed within each stellar-mass bin.
    Returns per-host weights normalised to 1 within each bin, mu, the lost tail mass (max over bins) and the
    multiplicity actually realised after dropping that tail (host-weighted mean over bins)."""
    mu = mu_for(f)
    kmax = int(pool.k.max())
    q_all = ztp_pmf(np.arange(kmax + 1), mu)
    w = np.zeros(len(pool.k))
    lost, f_real, n_tot = 0.0, 0.0, 0
    for j in np.unique(mbin_host):
        sel = mbin_host == j
        N_k = np.bincount(pool.k[sel], minlength=kmax + 1).astype(float)
        usable = N_k >= k_min_pool
        lost = max(lost, float(q_all[~usable].sum() + max(0.0, 1 - q_all.sum())))
        q = np.where(usable, q_all, 0.0)
        q /= q.sum()
        f_real += float(q[2:].sum()) * sel.sum()
        n_tot += int(sel.sum())
        with np.errstate(invalid="ignore", divide="ignore"):
            per = np.where(N_k > 0, q / N_k, 0.0)
        w[sel] = per[pool.k[sel]]
    return w, mu, lost, f_real / max(n_tot, 1)


def in_window(m, P, m_win, P_win):
    return (m >= m_win[0]) & (m < m_win[1]) & (P >= P_win[0]) & (P < P_win[1])
