"""Keplerian reflex motion of a star with one or more dark companions.

Conventions follow gaiamock (El-Badry et al. 2024) exactly, so that data
simulated here can be passed straight into gaiamock's orbit fitter:

* time ``t_yr`` is in years relative to the data-release reference epoch;
* periods ``P`` and periastron times ``Tp`` are in days;
* mean anomaly  M = 2*pi*(t_yr*365.25 - Tp)/P  =  2*pi*t_yr*365.25/P - phi_p,
  with phi_p = 2*pi*Tp/P;
* Thiele-Innes constants (A, B, F, G) in mas, with
  delta_ra* = B*X + G*Y  and  delta_dec = A*X + F*Y,
  X = cos E - e,  Y = sqrt(1 - e^2) sin E;
* along-scan displacement  delta_eta = delta_ra* sin(psi) + delta_dec cos(psi).

Angles: ``inc`` (0 = face-on), ``Omega`` (node) and ``omega`` (argument of
periastron *of the star's reflex orbit*) are in radians.

Multi-companion systems are modelled as a superposition of Keplerian reflex
orbits (Jacobi-like hierarchy, no planet-planet interactions). That is adequate
for the ~5.5 yr DR4 baseline unless the planets are close to resonance or
strongly interacting.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .constants import AU_M, G_SI, MJUP_MSUN, MSUN_KG, DAY_S


@dataclass(frozen=True)
class Companion:
    """A dark companion (planet) on a Keplerian orbit around the star."""

    m_mjup: float      # companion mass [M_Jup]
    P_days: float      # orbital period [d]
    ecc: float = 0.0
    inc: float = 0.0   # [rad], 0 = face-on
    Omega: float = 0.0  # [rad]
    omega: float = 0.0  # [rad], of the star's reflex orbit
    Tp_days: float = 0.0  # time of periastron [d], relative to reference epoch

    @property
    def phi_p(self) -> float:
        return 2.0 * np.pi * self.Tp_days / self.P_days


def solve_kepler(M: np.ndarray, ecc: float, tol: float = 1e-12, maxiter: int = 50) -> np.ndarray:
    """Eccentric anomaly E from mean anomaly M (vectorised Newton iteration)."""
    M = np.mod(np.asarray(M, dtype=float), 2.0 * np.pi)
    E = M + ecc * np.sin(M) if ecc < 0.8 else np.full_like(M, np.pi)
    for _ in range(maxiter):
        f = E - ecc * np.sin(E) - M
        dE = f / (1.0 - ecc * np.cos(E))
        E = E - dE
        if np.max(np.abs(dE)) < tol:
            break
    return E


def relative_semimajor_axis_au(P_days: float, m_total_msun: float) -> float:
    """Kepler's third law: semi-major axis of the relative orbit [AU]."""
    P_s = P_days * DAY_S
    return (G_SI * m_total_msun * MSUN_KG * P_s**2 / (4.0 * np.pi**2)) ** (1.0 / 3.0) / AU_M


def reflex_a0_mas(m_star_msun: float, m_comp_mjup: float, P_days: float,
                  parallax_mas: float, m_inner_msun: float = 0.0) -> float:
    """Angular semi-major axis of the star's reflex orbit [mas].

    ``m_inner_msun`` is the mass of any companions on smaller orbits, which in
    a Jacobi hierarchy orbit together with the star.
    """
    m_c = m_comp_mjup * MJUP_MSUN
    m_central = m_star_msun + m_inner_msun
    a_rel = relative_semimajor_axis_au(P_days, m_central + m_c)
    return a_rel * m_c / (m_central + m_c) * parallax_mas


def thiele_innes(a0_mas: float, inc: float, Omega: float, omega: float):
    """Thiele-Innes constants (A, B, F, G) in mas (gaiamock convention)."""
    ci = np.cos(inc)
    cO, sO, cw, sw = np.cos(Omega), np.sin(Omega), np.cos(omega), np.sin(omega)
    A = a0_mas * (cw * cO - sw * sO * ci)
    B = a0_mas * (cw * sO + sw * cO * ci)
    F = -a0_mas * (sw * cO + cw * sO * ci)
    G = -a0_mas * (sw * sO - cw * cO * ci)
    return A, B, F, G


def orbit_XY(t_yr: np.ndarray, P_days: float, ecc: float, phi_p: float):
    """Elliptical rectangular coordinates X, Y of the unit orbit."""
    M = 2.0 * np.pi * np.asarray(t_yr) * 365.25 / P_days - phi_p
    E = solve_kepler(M, ecc)
    return np.cos(E) - ecc, np.sqrt(1.0 - ecc**2) * np.sin(E)


def sky_offsets(t_yr: np.ndarray, A, B, F, G, P_days, ecc, phi_p):
    """(delta_ra*, delta_dec) of the photocentre in mas for one orbit."""
    X, Y = orbit_XY(t_yr, P_days, ecc, phi_p)
    return B * X + G * Y, A * X + F * Y


def companions_ordered(companions):
    """Companions sorted by period (inner first)."""
    return sorted(companions, key=lambda c: c.P_days)


def reflex_along_scan(t_yr, psi, companions, m_star_msun: float, parallax_mas: float):
    """Total along-scan reflex displacement [mas] from all companions.

    Returns (delta_eta, a0_list) where a0_list gives each companion's reflex
    a0 in mas, in order of increasing period.
    """
    t_yr = np.asarray(t_yr, dtype=float)
    psi = np.asarray(psi, dtype=float)
    total = np.zeros_like(t_yr)
    a0s = []
    m_inner = 0.0
    for c in companions_ordered(companions):
        a0 = reflex_a0_mas(m_star_msun, c.m_mjup, c.P_days, parallax_mas, m_inner)
        A, B, F, G = thiele_innes(a0, c.inc, c.Omega, c.omega)
        dra, ddec = sky_offsets(t_yr, A, B, F, G, c.P_days, c.ecc, c.phi_p)
        total += dra * np.sin(psi) + ddec * np.cos(psi)
        a0s.append(a0)
        m_inner += c.m_mjup * MJUP_MSUN
    return total, a0s


def campbell_from_thiele_innes(A, B, F, G):
    """Convert (A, B, F, G) to (a0, inc, Omega, omega), angles in rad.

    Standard inversion (e.g. Halbwachs et al. 2023). Omega is returned in
    [0, pi) with omega adjusted accordingly; inc in [0, pi].
    """
    u = 0.5 * (A**2 + B**2 + F**2 + G**2)
    v = A * G - B * F
    a0 = np.sqrt(u + np.sqrt(np.maximum((u + v) * (u - v), 0.0)))
    w_plus_O = np.arctan2(B - F, A + G)
    w_minus_O = np.arctan2(-B - F, A - G)
    omega = 0.5 * (w_plus_O + w_minus_O)
    Omega = 0.5 * (w_plus_O - w_minus_O)
    if Omega < 0:
        Omega += np.pi
        omega += np.pi
    Omega = np.mod(Omega, np.pi)
    omega = np.mod(omega, 2 * np.pi)
    # (A+G, B-F) = a0 (1+cos i) (cos, sin)(omega+Omega)
    # (A-G, -(B+F)) = a0 (1-cos i) (cos, sin)(omega-Omega)
    # => tan^2(i/2) = |(A-G, B+F)| / |(A+G, B-F)|
    num = np.hypot(A - G, B + F)
    den = np.hypot(A + G, B - F)
    inc = 2.0 * np.arctan(np.sqrt(num / den)) if den > 0 else np.pi / 2
    return a0, inc, Omega, omega


def companion_mass_from_a0(a0_mas: float, P_days: float, m_star_msun: float,
                           parallax_mas: float) -> float:
    """Invert reflex_a0_mas for the companion mass [M_Jup] (dark companion)."""
    from scipy.optimize import brentq

    if a0_mas <= 0:
        return 0.0

    def f(m_mjup):
        return reflex_a0_mas(m_star_msun, m_mjup, P_days, parallax_mas) - a0_mas

    hi = 1.0
    while f(hi) < 0 and hi < 1e6:
        hi *= 2.0
    if f(hi) < 0:
        return np.inf
    return brentq(f, 1e-8, hi, xtol=1e-10)


def hill_stable(m_star_msun: float, c1: Companion, c2: Companion) -> bool:
    """Conservative two-planet stability test.

    Requires the gap between the inner planet's apocentre and the outer
    planet's pericentre to exceed 2*sqrt(3) mutual Hill radii (Gladman 1993
    criterion for circular orbits, applied to the closest approach, so
    eccentric pairs need wider spacing - as Lammers & Winn 2025 note).
    """
    inner, outer = companions_ordered([c1, c2])
    m_in = inner.m_mjup * MJUP_MSUN
    m_out = outer.m_mjup * MJUP_MSUN
    a_in = relative_semimajor_axis_au(inner.P_days, m_star_msun + m_in)
    a_out = relative_semimajor_axis_au(outer.P_days, m_star_msun + m_in + m_out)
    r_hill = ((m_in + m_out) / (3 * m_star_msun)) ** (1 / 3) * 0.5 * (a_in + a_out)
    gap = a_out * (1 - outer.ecc) - a_in * (1 + inner.ecc)
    return gap > 2 * np.sqrt(3) * r_hill


def draw_ecc_kipping(rng, size=None):
    """Eccentricity from Beta(0.867, 3.03) (Kipping 2013), as in Lammers & Winn."""
    return np.minimum(rng.beta(0.867, 3.03, size), 0.95)


# ----------------------------------------------------------------------------- geometry helpers
def _orbit_normal(inc: float, Omega: float) -> np.ndarray:
    """Unit orbit normal for (inc, Omega); any fixed sky convention works for angles between orbits."""
    return np.array([np.sin(inc) * np.sin(Omega), -np.sin(inc) * np.cos(Omega), np.cos(inc)])


def mutual_inclination(inc1: float, Omega1: float, inc2: float, Omega2: float) -> float:
    """Angle between two orbital planes [rad]:
    cos(Delta) = cos i1 cos i2 + sin i1 sin i2 cos(Omega1 - Omega2).

    Note: astrometry alone cannot tell (Omega, omega) from (Omega + pi, omega + pi),
    so the same pair of sky motions is consistent with two values of Delta.
    """
    c = (np.cos(inc1) * np.cos(inc2)
         + np.sin(inc1) * np.sin(inc2) * np.cos(Omega1 - Omega2))
    return float(np.arccos(np.clip(c, -1.0, 1.0)))


def orientation_relative_to(inc_ref: float, Omega_ref: float, mutual_inc: float,
                            azimuth: float) -> tuple[float, float]:
    """(inc, Omega) of an orbit tilted by `mutual_inc` from the reference orbit.

    `azimuth` (rad) says in which direction around the reference orbit's
    normal the tilt is applied. Returns inc in [0, pi], Omega in [0, 2 pi).
    """
    n_ref = _orbit_normal(inc_ref, Omega_ref)
    e1 = np.cross([0.0, 0.0, 1.0], n_ref)
    if np.linalg.norm(e1) < 1e-12:            # reference orbit face-on
        e1 = np.array([1.0, 0.0, 0.0])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(n_ref, e1)
    n = (np.cos(mutual_inc) * n_ref
         + np.sin(mutual_inc) * (np.cos(azimuth) * e1 + np.sin(azimuth) * e2))
    inc = float(np.arccos(np.clip(n[2], -1.0, 1.0)))
    Omega = float(np.arctan2(n[0], -n[1]) % (2 * np.pi))
    return inc, Omega


def Tp_for_mean_anomaly(M_rad: float, t_days: float, P_days: float) -> float:
    """Periastron time [d, relative to the reference epoch] that puts the orbit at
    mean anomaly `M_rad` at time `t_days` (same convention as orbit_XY)."""
    return float(t_days - M_rad / (2 * np.pi) * P_days)
