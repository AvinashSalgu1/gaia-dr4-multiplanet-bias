"""Simulated Gaia epoch astrometry for a star with any number of companions.

Generalises gaiamock.predict_astrometry_binary_in_terms_of_a0 (single dark
companion) to N companions, with the same noise model so that results with
one companion are directly comparable with gaiamock and Lammers & Winn (2026).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .constants import N_CCD_PER_TRANSIT
from .epochs import Epochs
from .orbits import Companion, reflex_along_scan

# Per-CCD along-scan uncertainty [mas] vs G: the "EDR3 adjusted" curve from
# Fig. 3 of arXiv:2206.05439 (Gaia DR3 astrometric orbit paper, Holl et al.
# 2023a, the source Lammers & Winn 2026 cite for their Fig. 4), as tabulated in
# gaiamock.al_uncertainty_per_ccd_interp (MIT licence, (c) 2024 K. El-Badry).
_G_VALS = np.array([4, 5, 6, 7, 8.2, 8.4, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20])
_SIGMA_CCD = np.array([0.4, 0.35, 0.15, 0.17, 0.23, 0.13, 0.13, 0.135, 0.125, 0.13,
                       0.15, 0.23, 0.36, 0.63, 1.05, 2.05, 4.1])


def sigma_al_per_transit(G_mag) -> np.ndarray:
    """Along-scan uncertainty per field-of-view transit [mas] (8 CCDs averaged)."""
    return np.interp(G_mag, _G_VALS, _SIGMA_CCD) / np.sqrt(N_CCD_PER_TRANSIT)


@dataclass
class Star:
    parallax_mas: float
    m_star_msun: float
    G_mag: float
    pmra: float = 0.0   # mas/yr
    pmdec: float = 0.0  # mas/yr
    ra_off: float = 0.0  # mas, offset from reference position
    dec_off: float = 0.0


@dataclass
class SimulatedAstrometry:
    epochs: Epochs
    obs: np.ndarray          # along-scan positions [mas]
    err: np.ndarray          # reported uncertainties [mas]
    signal: np.ndarray       # noise-free reflex part [mas]
    a0_mas: list = field(default_factory=list)  # reflex a0 per companion (inner first)
    jitter_mas: float = 0.0  # unmodelled extra noise actually added


def simulate(epochs: Epochs, star: Star, companions: list[Companion],
             rng: np.random.Generator | None = None, add_noise: bool = True,
             bright_star_jitter: bool = True) -> SimulatedAstrometry:
    """Epoch astrometry Lambda = 5-parameter motion + sum of reflex orbits + noise.

    Noise follows gaiamock: Gaussian with the per-transit sigma, plus for
    G < 13 an unmodelled jitter with sigma drawn from U(0, 0.04) mas that is
    *not* included in the reported uncertainties.
    """
    rng = np.random.default_rng() if rng is None else rng
    t, psi, pf = epochs.t_yr, epochs.psi, epochs.plx_factor
    sin_psi, cos_psi = np.sin(psi), np.cos(psi)

    five_par = (star.ra_off * sin_psi + star.dec_off * cos_psi
                + star.pmra * t * sin_psi + star.pmdec * t * cos_psi
                + star.parallax_mas * pf)
    signal, a0s = reflex_along_scan(t, psi, companions, star.m_star_msun, star.parallax_mas)

    sigma = float(sigma_al_per_transit(star.G_mag))
    err = np.full(len(t), sigma)
    obs = five_par + signal
    jitter = 0.0
    if add_noise:
        obs = obs + sigma * rng.standard_normal(len(t))
        if bright_star_jitter and star.G_mag < 13:
            jitter = rng.uniform(0, 0.04)
            obs = obs + jitter * rng.standard_normal(len(t))
    return SimulatedAstrometry(epochs, obs, err, signal, a0s, jitter)


def snr1(a0_mas: float, G_mag: float) -> float:
    """Lammers & Winn (2026) per-transit signal-to-noise, SNR_1 = alpha / sigma_fov."""
    return a0_mas / float(sigma_al_per_transit(G_mag))
