"""Simulated two-planet systems for step 2.

The draw (``draw_system``) and data generation (``make_data``) are copied
unchanged from the original step-2 script: the random numbers are consumed in
exactly the same order, so a given seed gives the same system, epochs and
noise as before and earlier results can be reused. Do not reorder rng calls.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

T_DR4_DAYS = 2458868.5 - 2456891.5     # 1977 d = 5.41 yr
F2_CUT = 25.0                           # Gaia DR3 goodness-of-fit cut
LOCK_TOL = 0.1                          # |ln(P_fit / P)| < 0.1 -> "locked" on that planet
SEED0 = 7                               # master seed of the step-2 population
RESULTS = Path(__file__).resolve().parents[3] / "results"


def population_seeds(n: int, seed0: int = SEED0) -> list[int]:
    """Seeds of the first n systems (a prefix of the full 20,000)."""
    return [int(s) for s in np.random.SeedSequence(seed0).generate_state(max(n, 20000), dtype=np.uint32)[:n]]


def _random_orientation(rng):
    return dict(inc=np.arccos(rng.uniform(-1, 1)), Omega=rng.uniform(0, np.pi),
                omega=rng.uniform(0, 2 * np.pi))


def draw_system(rng):
    """Target in the detectable range + a stable companion anywhere."""
    from ..orbits import Companion, companion_mass_from_a0, draw_ecc_kipping, hill_stable
    from ..simulate import sigma_al_per_transit

    while True:
        snr1 = 10 ** rng.uniform(np.log10(1.5), np.log10(10.0))
        P_t = 10 ** rng.uniform(np.log10(100.0), np.log10(0.7 * T_DR4_DAYS))
        G = rng.uniform(8.0, 15.0)
        m_star = rng.uniform(0.4, 1.3)
        parallax = 10 ** rng.uniform(np.log10(2.0), np.log10(100.0))
        m_t = companion_mass_from_a0(snr1 * float(sigma_al_per_transit(G)), P_t, m_star, parallax)
        if not 0.3 < m_t < 13.0:
            continue
        target = Companion(m_t, P_t, float(draw_ecc_kipping(rng)), Tp_days=rng.uniform(0, P_t),
                           **_random_orientation(rng))
        a0_t = snr1 * float(sigma_al_per_transit(G))
        for _ in range(500):
            P_c = 10 ** rng.uniform(np.log10(30.0), np.log10(10.0 * T_DR4_DAYS))
            rho = 10 ** rng.uniform(-2.0, 2.0)
            m_c = companion_mass_from_a0(rho * a0_t, P_c, m_star, parallax)
            if not 0.3 < m_c < 13.0:
                continue
            comp = Companion(m_c, P_c, float(draw_ecc_kipping(rng)), Tp_days=rng.uniform(0, P_c),
                             **_random_orientation(rng))
            if hill_stable(m_star, target, comp):
                break
        else:
            continue
        sky = dict(ra=rng.uniform(0, 360), dec=np.degrees(np.arcsin(rng.uniform(-1, 1))))
        return dict(G=G, m_star=m_star, parallax=parallax, pmra=rng.normal(0, 40),
                    pmdec=rng.normal(0, 40), target=target, comp=comp, **sky)


def make_data(seed: int):
    """One system's data: (s, ep, single, both, noise).

    single / both are noise-free simulations (target alone / target +
    companion); the SAME noise vector is added to both.
    """
    from .. import epochs as E
    from ..simulate import Star, sigma_al_per_transit, simulate

    rng = np.random.default_rng(seed)
    s = draw_system(rng)
    ep = E.from_gaiamock(s["ra"], s["dec"], "dr4", rng=rng)
    star = Star(parallax_mas=s["parallax"], m_star_msun=s["m_star"], G_mag=s["G"],
                pmra=s["pmra"], pmdec=s["pmdec"])
    t, c = s["target"], s["comp"]
    single = simulate(ep, star, [t], add_noise=False)
    both = simulate(ep, star, [t, c], add_noise=False)
    sigma = float(sigma_al_per_transit(s["G"]))
    noise = sigma * rng.standard_normal(len(ep))
    if s["G"] < 13:
        noise += rng.uniform(0, 0.04) * rng.standard_normal(len(ep))
    return s, ep, single, both, noise


META_COLS = ["seed", "n_transits", "G", "m_star", "parallax",
             "m_t", "P_t", "e_t", "a0_t", "snr1_t", "m_c", "P_c", "e_c", "a0_c", "snr1_c", "rho",
             "inc_t", "inc_c", "mutual_inc", "M_mid_c", "peri_in_window",
             "dchi2_seen_t", "dchi2_seen_c", "dchi2_left_c_accel"]


def describe(seed, s, ep, single, both) -> dict:
    """True parameters plus derived quantities used by the 'what else matters' stage."""
    from ..fitting import fit_5par, fit_trend_only
    from ..orbits import companions_ordered, mutual_inclination
    from ..simulate import snr1

    t, c = s["target"], s["comp"]
    a0_t = single.a0_mas[0]
    a0_c = [a for a, cc in zip(both.a0_mas, companions_ordered([t, c])) if cc is c][0]
    err = single.err
    sig_c = both.signal - single.signal
    t_lo, t_hi = ep.t_yr.min() * 365.25, ep.t_yr.max() * 365.25
    t_mid = 0.5 * (t_lo + t_hi)
    M_mid = (2 * np.pi * (t_mid - c.Tp_days) / c.P_days + np.pi) % (2 * np.pi) - np.pi
    k_lo, k_hi = np.ceil((t_lo - c.Tp_days) / c.P_days), np.floor((t_hi - c.Tp_days) / c.P_days)
    return dict(
        seed=seed, n_transits=len(ep), G=s["G"], m_star=s["m_star"], parallax=s["parallax"],
        m_t=t.m_mjup, P_t=t.P_days, e_t=t.ecc, a0_t=a0_t, snr1_t=snr1(a0_t, s["G"]),
        m_c=c.m_mjup, P_c=c.P_days, e_c=c.ecc, a0_c=a0_c, snr1_c=snr1(a0_c, s["G"]), rho=a0_c / a0_t,
        inc_t=t.inc, inc_c=c.inc, mutual_inc=mutual_inclination(t.inc, t.Omega, c.inc, c.Omega),
        M_mid_c=abs(M_mid), peri_in_window=int(k_lo <= k_hi),
        # chi2 of each planet's noise-free wobble left after removing position, proper
        # motion and parallax (= how much signal Gaia actually "sees"), and of the
        # companion's wobble left after also removing an acceleration
        dchi2_seen_t=fit_5par(ep, single.signal, err)[0],
        dchi2_seen_c=fit_5par(ep, sig_c, err)[0],
        dchi2_left_c_accel=fit_trend_only(ep, sig_c, err, 2)[0],
    )
