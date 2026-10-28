"""Synthetic stars and giant-planet populations for step 3 (occurrence rates).

Stars
-----
A mock main-sequence sample. Two distance limits are available:

* "lw": L&W's (AJ 171:18, Sec. 2.1) M_G-dependent r_max(M_G), within which a
  13 M_J planet on a 9.5-yr orbit is detectable at N_sigma = 0.5. Very
  permissive: the median star is ~1.6 kpc away and only ~1 in 6000 stars has
  a planet Gaia could detect, so L&W need tens of millions of real stars.
* "deep" (default for step 3): additionally, a REFERENCE planet of M_REF_MJ
  (1 M_J) on a P = T_DR4 orbit must give SNR_1 >= 1. This keeps the stars
  around which typical giants are detectable, so far fewer stars (and fits)
  are needed. Occurrence rates are per star of whatever sample is defined,
  and the completeness is measured on the same stars, so the choice changes
  the precision, not the validity, of the comparison. It does shift the
  detections towards higher SNR (fewer near threshold) than a real
  magnitude-limited DR4 sample.

Both also apply G <= 19. Unlike L&W (real Gaia DR3 stars) the stars are synthetic:

* masses 0.3-1.4 M_sun from a Kroupa (2001) IMF;
* absolute G from an approximate main-sequence mass-M_G relation
  (Pecaut & Mamajek 2013 M_V, converted to G; good to ~0.3 mag);
* positions uniform in space with an exponential disc (scale height 300 pc,
  Sun in the mid-plane); no extinction;
* proper motions N(0, 40 mas/yr) per axis, as in steps 1-2.

Planets
-------
Occurrence follows L&W Eq. 23 (Fulton et al. 2021 shape):

    dN / (dN_star dln a dln m) = C (M_star / 0.9 M_sun) (a/au)^beta [1 - exp(-(a/a0)^gamma)]

with C = 0.15, beta = -0.86, a0 = 3.6 au, gamma = 1.59, over 0.3-13 M_J and
0.1-30 au, uniform in ln m. This gives 15.4 giants per 100 stars (0.9 M_sun)
at 2-8 au, matching the CLS value of 14.1 (+2.0/-1.8), and 0.367 giants per
star over the full range.

Two multiplicity scenarios share this SAME mean occurrence (same planets per
star), and differ only in how planets are shared between stars:

* "independent": each star's number of giants is Poisson(lambda(M_star));
  planets are drawn independently. Multiplicity is the minimum implied by the
  occurrence rate.
* "clustered": only a fraction h of stars can host giants, and those host
  Poisson(lambda/h). h is set so that a star with a giant has the outer
  companion fraction measured by Bryan et al. (2016): 52 +/- 5 % for 1-20 M_J
  at 5-20 au, i.e. 44.5 % for 1-13 M_J when scaled with the ln-uniform mass
  distribution used here.

Every system with two or more planets must pass the pairwise Hill stability
test of step 2 (orbits.hill_stable). A planet that would make its system
unstable is redrawn (up to N_REDRAW times) and dropped if it never fits; the
"true" occurrence used to judge any inferred rate is always the realised set
of planets, so this does not bias the comparison. In the clustered scenario
(~3 giants per host) this drops ~5% of planets and preferentially redraws
close pairs, so its realised rate is slightly lower and its period
distribution slightly wider than the independent one; each scenario is always
compared with its own realised truth.

Variants (draw_planets(..., variant=...)) change only systems with two or more giants:

* "coplanar": the first planet's orientation is random; every further planet is
  tilted from it by a mutual inclination drawn from a Rayleigh distribution with
  sigma = 2 deg (the scale measured for Kepler multi-planet systems, Fabrycky et
  al. 2014), in a random direction. By default orientations are independent.
* "stable": every multi-planet system is integrated with rebound (WHFast, MEGNO)
  for STABLE_ORBITS orbits of its innermost planet; a close encounter (one mutual
  Hill radius), an ejection or MEGNO > 4 marks it unstable, and the whole system
  is redrawn (same number of planets) up to STABLE_ROUNDS times. A system that
  is still unstable keeps only its first planet. Needs rebound.
* "lowecc": eccentricities of planets in multi-planet systems from a Rayleigh
  distribution with sigma = 0.1 instead of Kipping (2013); a bounding case for
  dynamically cold systems (single planets keep Kipping).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .constants import MJUP_MSUN
from .orbits import Companion, draw_ecc_kipping, hill_stable, reflex_a0_mas

# ---------------------------------------------------------------- occurrence model
LW_C, LW_BETA, LW_A0_AU, LW_GAMMA = 0.15, -0.86, 3.6, 1.59
M_RANGE_MJ = (0.3, 13.0)
A_RANGE_AU = (0.1, 30.0)
M_STAR_NORM = 0.9
BRYAN_COMP_FRAC = 0.52 * np.log(13.0) / np.log(20.0)   # 1-13 M_J instead of 1-20 M_J
BRYAN_BOX = dict(m=(1.0, 13.0), a=(5.0, 20.0))
N_REDRAW = 20
MUTUAL_INC_SIGMA = np.radians(2.0)
LOWECC_SIGMA = 0.1
STABLE_ORBITS = 1e4
STABLE_ROUNDS = 6
VARIANTS = ("coplanar", "stable", "lowecc")

_LNA = np.linspace(np.log(A_RANGE_AU[0]), np.log(A_RANGE_AU[1]), 4001)


def a_shape(a_au):
    """Semi-major-axis shape of L&W Eq. 23 (per ln a, without C)."""
    a_au = np.asarray(a_au, float)
    return a_au ** LW_BETA * (1.0 - np.exp(-(a_au / LW_A0_AU) ** LW_GAMMA))


def _cum_shape():
    y = a_shape(np.exp(_LNA))
    c = np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(_LNA))])
    return c


_CUM = _cum_shape()


def rate_in_box(m_lo, m_hi, a_lo, a_hi, m_star=M_STAR_NORM) -> float:
    """Expected number of giants per star with m in [m_lo, m_hi] M_J and a in [a_lo, a_hi] au."""
    m_lo, m_hi = max(m_lo, M_RANGE_MJ[0]), min(m_hi, M_RANGE_MJ[1])
    if m_hi <= m_lo:
        return 0.0
    ia = np.interp(np.log(np.clip([a_lo, a_hi], *A_RANGE_AU)), _LNA, _CUM)
    return LW_C * (m_star / M_STAR_NORM) * np.log(m_hi / m_lo) * float(ia[1] - ia[0])


def lambda_total(m_star) -> np.ndarray:
    """Expected number of giants per star (0.3-13 M_J, 0.1-30 au)."""
    return LW_C * (np.asarray(m_star, float) / M_STAR_NORM) * np.log(M_RANGE_MJ[1] / M_RANGE_MJ[0]) * _CUM[-1]


def draw_a_au(rng, size):
    return np.exp(np.interp(rng.uniform(0, _CUM[-1], size), _CUM, _LNA))


def draw_m_mjup(rng, size):
    return np.exp(rng.uniform(np.log(M_RANGE_MJ[0]), np.log(M_RANGE_MJ[1]), size))


def host_fraction_clustered(p_comp: float = BRYAN_COMP_FRAC) -> float:
    """h such that P(>=1 companion in the Bryan box | star hosts a giant) = p_comp.

    Within a host, planets are Poisson with rate lambda/h, so
    P(companion in box | host) = 1 - exp(-lambda_box / h), independent of
    the inner planet. Calibrated at M_star = 0.9 M_sun.
    """
    lam_box = rate_in_box(*BRYAN_BOX["m"], *BRYAN_BOX["a"])
    return lam_box / -np.log(1.0 - p_comp)


def period_days(a_au, m_star, m_mjup):
    """Kepler's third law for the relative orbit."""
    m_tot = np.asarray(m_star) + np.asarray(m_mjup) * MJUP_MSUN
    return 365.25 * np.sqrt(np.asarray(a_au) ** 3 / m_tot)


# ---------------------------------------------------------------- stars
# Approximate main-sequence (mass [M_sun], M_G) relation: Pecaut & Mamajek
# (2013, updated table) M_V converted with G - V from the same table.
_MS_MASS = np.array([0.23, 0.37, 0.44, 0.50, 0.57, 0.64, 0.73, 0.82, 0.88, 0.94, 1.00, 1.08, 1.18, 1.33, 1.45])
_MS_MG = np.array([11.5, 10.3, 9.4, 8.85, 8.3, 7.6, 6.95, 6.1, 5.65, 5.3, 4.67, 4.25, 3.9, 3.33, 2.9])
STAR_MASS_RANGE = (0.3, 1.4)
G_MAX = 19.0
Z_SCALE_PC = 300.0
T_DR4_DAYS = 2458868.5 - 2456891.5   # 1977 d
M_REF_MJ, SNR_REF = 1.0, 1.0           # "deep" sample: 1 M_J at P = T_DR4 gives SNR_1 >= 1


def abs_G(m_star):
    return np.interp(m_star, _MS_MASS, _MS_MG)


def r_max_lw(M_G):
    """L&W Sec. 4.1 (Eqs. 36-37) distance limit [pc] for their stellar sample.

    The limit peaks at 1600 pc for M_G = 2.8 and falls on both sides: for
    brighter stars it "decreases with luminosity" (Eq. 36), for fainter stars
    it "increases with luminosity" (Eq. 37), i.e. fainter = closer.
    (Corrected 2026-10-07: both exponents had the wrong sign.)
    """
    M_G = np.asarray(M_G, float)
    return np.where(M_G <= 2.8, 1600.0 * 10 ** (0.069 * (M_G - 2.8)),
                    1600.0 * 10 ** (-0.093 * (M_G - 2.8)))


def r_ref_pc(m_star, m_ref=M_REF_MJ, P_ref=T_DR4_DAYS, snr_ref=SNR_REF):
    """Largest distance at which the reference planet has SNR_1 >= snr_ref.

    SNR_1 = a0 / sigma_fov(G) falls with distance (a0 ~ 1/d, sigma grows as
    the star gets fainter), so scan a fine distance grid from far to near.
    """
    from .orbits import reflex_a0_mas
    from .simulate import sigma_al_per_transit

    m_star = np.atleast_1d(np.asarray(m_star, float))
    r = np.geomspace(1.0, 2e4, 3000)                      # pc
    out = np.empty(len(m_star))
    # evaluate on a mass grid and interpolate (fast for millions of stars)
    mg = np.linspace(STAR_MASS_RANGE[0] - 0.01, STAR_MASS_RANGE[1] + 0.01, 120)
    rl = np.empty(len(mg))
    for k, m in enumerate(mg):
        snr = reflex_a0_mas(m, m_ref, P_ref, 1000.0 / r) / sigma_al_per_transit(abs_G(m) + 5 * np.log10(r / 10))
        ok = np.nonzero(snr >= snr_ref)[0]
        rl[k] = r[ok[-1]] if len(ok) else 0.0
    out[:] = np.interp(m_star, mg, rl)
    return out


def r_max_pc(m_star, sample: str = "deep"):
    """Distance limit: L&W's r_max and G <= G_MAX, plus the reference-planet limit for "deep"."""
    M_G = abs_G(m_star)
    r = np.minimum(r_max_lw(M_G), 10 ** ((G_MAX - M_G) / 5 + 1))
    if sample == "deep":
        r = np.minimum(r, r_ref_pc(m_star))
    elif sample != "lw":
        raise ValueError(sample)
    return r


def draw_imf(rng, size, lo=STAR_MASS_RANGE[0], hi=STAR_MASS_RANGE[1]):
    """Kroupa (2001): dN/dm ~ m^-1.3 below 0.5 M_sun, m^-2.3 above (continuous)."""
    grid = np.linspace(lo, hi, 2001)
    pdf = np.where(grid < 0.5, grid ** -1.3, 0.5 * grid ** -2.3)
    cdf = np.concatenate([[0.0], np.cumsum(0.5 * (pdf[1:] + pdf[:-1]) * np.diff(grid))])
    return np.interp(rng.uniform(0, cdf[-1], size), cdf, grid)


# ICRS -> Galactic rotation matrix (Hipparcos convention); we need its inverse.
_R_GAL = np.array([[-0.0548755604, -0.8734370902, -0.4838350155],
                   [+0.4941094279, -0.4448296300, +0.7469822445],
                   [-0.8676661490, -0.1980763734, +0.4559837762]])


@dataclass
class StarSample:
    m_star: np.ndarray
    G: np.ndarray
    parallax: np.ndarray   # mas
    ra: np.ndarray         # deg
    dec: np.ndarray        # deg
    pmra: np.ndarray
    pmdec: np.ndarray

    def __len__(self):
        return len(self.m_star)

    def take(self, idx):
        return StarSample(*(getattr(self, k)[idx] for k in self.__dataclass_fields__))


def draw_stars(n: int, rng, sample: str = "deep") -> StarSample:
    """n stars of the mock sample (see module docstring)."""
    r_big = float(r_max_pc(np.linspace(*STAR_MASS_RANGE, 200), sample).max())
    out = []
    got = 0
    while got < n:
        k = max(4 * (n - got), 10000)
        m = draw_imf(rng, k)
        # uniform in a sphere of radius r_big, then disc and distance limits
        xyz = rng.normal(size=(k, 3))
        xyz /= np.linalg.norm(xyz, axis=1)[:, None]
        r = r_big * rng.uniform(0, 1, k) ** (1 / 3)
        gal = xyz * r[:, None]                          # galactic cartesian, pc
        keep = (r < r_max_pc(m, sample)) & (rng.uniform(0, 1, k) < np.exp(-np.abs(gal[:, 2]) / Z_SCALE_PC))
        m, gal, r = m[keep], gal[keep], r[keep]
        icrs = gal @ _R_GAL                             # inverse rotation
        ra = np.degrees(np.arctan2(icrs[:, 1], icrs[:, 0])) % 360
        dec = np.degrees(np.arcsin(icrs[:, 2] / r))
        out.append((m, r, ra, dec))
        got += len(m)
    m, r, ra, dec = (np.concatenate(x)[:n] for x in zip(*out))
    G = abs_G(m) + 5 * np.log10(r / 10)
    return StarSample(m, G, 1000.0 / r, ra, dec, rng.normal(0, 40, n), rng.normal(0, 40, n))


# ---------------------------------------------------------------- planets
@dataclass
class PlanetTable:
    """All planets of a population, one row per planet."""
    star: np.ndarray       # index into the StarSample
    m: np.ndarray          # M_J
    P: np.ndarray          # days
    a: np.ndarray          # au (relative orbit)
    e: np.ndarray
    inc: np.ndarray
    Omega: np.ndarray
    omega: np.ndarray
    Tp: np.ndarray         # days
    n_dropped: int = 0     # planets dropped for instability
    n_redrawn: int = 0

    def __len__(self):
        return len(self.star)

    def companion(self, i) -> Companion:
        return Companion(float(self.m[i]), float(self.P[i]), float(self.e[i]), float(self.inc[i]),
                         float(self.Omega[i]), float(self.omega[i]), float(self.Tp[i]))


def _one_planet(rng, m_star):
    a = float(draw_a_au(rng, 1)[0])
    m = float(draw_m_mjup(rng, 1)[0])
    P = float(period_days(a, m_star, m))
    return dict(m=m, P=P, a=a, e=float(draw_ecc_kipping(rng)), inc=float(np.arccos(rng.uniform(-1, 1))),
                Omega=float(rng.uniform(0, np.pi)), omega=float(rng.uniform(0, 2 * np.pi)),
                Tp=float(rng.uniform(0, P)))


def _as_comp(p):
    return Companion(p["m"], p["P"], p["e"], p["inc"], p["Omega"], p["omega"], p["Tp"])


def _draw_system(rng, ms, n, variant=()):
    """n planets on one star (Hill-stable pairs; variants as in the module docstring).
    Returns (planets, n_dropped, n_redrawn)."""
    from .orbits import orientation_relative_to
    sysp, dropped, redrawn = [], 0, 0
    for _ in range(n):
        for k in range(N_REDRAW + 1):
            p = _one_planet(rng, ms)
            if n >= 2 and "lowecc" in variant:
                p["e"] = float(min(rng.rayleigh(LOWECC_SIGMA), 0.95))
            if sysp and "coplanar" in variant:
                p["inc"], p["Omega"] = orientation_relative_to(sysp[0]["inc"], sysp[0]["Omega"],
                                                               float(rng.rayleigh(MUTUAL_INC_SIGMA)),
                                                               float(rng.uniform(0, 2 * np.pi)))
            c = _as_comp(p)
            if all(hill_stable(ms, c, _as_comp(q)) for q in sysp):
                sysp.append(p)
                redrawn += k > 0
                break
        else:
            dropped += 1
    return sysp, dropped, redrawn


def nbody_stable(ms, planets, n_orbits=STABLE_ORBITS) -> bool:
    """True if the system survives n_orbits of its innermost planet (rebound WHFast + MEGNO)."""
    import rebound
    if len(planets) < 2:
        return True
    sim = rebound.Simulation()
    sim.units = ("yr", "AU", "Msun")
    sim.add(m=ms)
    for p in sorted(planets, key=lambda q: q["P"]):
        P_yr = p["P"] / 365.25
        sim.add(m=p["m"] * MJUP_MSUN, P=P_yr, e=p["e"], inc=p["inc"], Omega=p["Omega"],
                omega=p["omega"] + np.pi,                      # planet = star's reflex omega + pi
                M=(-2 * np.pi * (p["Tp"] / 365.25) / P_yr) % (2 * np.pi))
    sim.move_to_com()
    ps = sim.particles
    r_hill = min(((ps[i].m + ps[i + 1].m) / (3 * ms)) ** (1 / 3) * 0.5 * (ps[i].a + ps[i + 1].a)
                 for i in range(1, len(ps) - 1))
    sim.integrator = "whfast"
    sim.dt = ps[1].P / 30.0
    sim.exit_max_distance = 100 * ps[len(ps) - 1].a
    sim.exit_min_distance = r_hill
    sim.init_megno()
    try:
        sim.integrate(n_orbits * ps[1].P, exact_finish_time=0)
        return bool(sim.megno() <= 4)
    except (rebound.Escape, rebound.Encounter):
        return False
    except Exception:
        return False


def draw_planets(stars: StarSample, rng, scenario: str = "independent", variant=(), jobs: int = 1,
                 info: dict | None = None) -> PlanetTable:
    """Planets for every star under the given multiplicity scenario and variants (module docstring).
    `info`, if given, receives counts of the stability screening."""
    variant = tuple(variant)
    bad = set(variant) - set(VARIANTS)
    if bad:
        raise ValueError(f"unknown variant(s) {sorted(bad)}; choose from {VARIANTS}")
    lam = lambda_total(stars.m_star)
    if scenario == "independent":
        n_pl = rng.poisson(lam)
    elif scenario == "clustered":
        h = host_fraction_clustered()
        host = rng.uniform(0, 1, len(stars)) < h
        n_pl = np.where(host, rng.poisson(lam / h), 0)
    else:
        raise ValueError(scenario)

    systems = {}
    dropped = redrawn = 0
    for s in np.nonzero(n_pl)[0]:
        sysp, d_, r_ = _draw_system(rng, float(stars.m_star[s]), int(n_pl[s]), variant)
        systems[s] = sysp
        dropped, redrawn = dropped + d_, redrawn + r_
    if "stable" in variant:
        from joblib import Parallel, delayed
        pending = [s for s, p in systems.items() if len(p) >= 2]
        n_multi, n_first_unstable, rounds = len(pending), None, 0
        for rounds in range(STABLE_ROUNDS + 1):
            if not pending:
                break
            ok = Parallel(n_jobs=jobs, batch_size=64)(delayed(nbody_stable)(float(stars.m_star[s]), systems[s])
                                                      for s in pending)
            pending = [s for s, good in zip(pending, ok) if not good]
            if n_first_unstable is None:
                n_first_unstable = len(pending)
            if rounds == STABLE_ROUNDS:
                break
            for s in pending:                                  # redraw the whole system
                sysp, d_, r_ = _draw_system(rng, float(stars.m_star[s]), len(systems[s]), variant)
                systems[s] = sysp
        for s in pending:                                      # still unstable: keep the first planet only
            dropped += len(systems[s]) - 1
            systems[s] = systems[s][:1]
        if info is not None:
            info.update(n_multi=n_multi, n_unstable_first=n_first_unstable or 0, n_reduced=len(pending),
                        stable_orbits=STABLE_ORBITS)
    rows = [(s, p) for s, sysp in systems.items() for p in sysp]

    def col(k):
        return np.array([p[k] for _, p in rows], float)

    return PlanetTable(np.array([s for s, _ in rows], int), col("m"), col("P"), col("a"), col("e"),
                       col("inc"), col("Omega"), col("omega"), col("Tp"), dropped, redrawn)


def planet_snr1(stars: StarSample, planets: PlanetTable):
    """Per-planet SNR_1 = a0 / sigma_fov (no Jacobi correction; used only for the pre-fit cut)."""
    from .simulate import sigma_al_per_transit

    s = planets.star
    a0 = reflex_a0_mas(stars.m_star[s], planets.m, planets.P, stars.parallax[s])
    return a0 / sigma_al_per_transit(stars.G[s])
