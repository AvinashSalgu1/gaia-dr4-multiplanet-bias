"""Gaia observation epochs (time, scan angle, along-scan parallax factor)."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .constants import DR_REF_JD


@dataclass
class Epochs:
    t_yr: np.ndarray        # years relative to the DR reference epoch
    psi: np.ndarray         # scan angle [rad]
    plx_factor: np.ndarray  # along-scan parallax factor

    def __len__(self) -> int:
        return len(self.t_yr)

    def subset(self, mask) -> "Epochs":
        return Epochs(self.t_yr[mask], self.psi[mask], self.plx_factor[mask])


def from_gaiamock(ra_deg: float, dec_deg: float, data_release: str = "dr4",
                  reject_frac: float = 0.1, rng: np.random.Generator | None = None) -> Epochs:
    """Real Gaia scanning law at (ra, dec) via gaiamock's pre-computed GOST files.

    Like gaiamock, a random ``reject_frac`` of field-of-view transits is dropped.
    Requires gaiamock and its healpix_scans/ directory (see README).
    """
    from .gaiamock_bridge import gaiamock

    gm = gaiamock()
    tab = gm.get_gost_one_position(ra_deg, dec_deg, data_release=data_release)
    rng = np.random.default_rng() if rng is None else rng
    if reject_frac > 0:
        tab = tab[rng.uniform(0, 1, len(tab)) > reject_frac]
    psi, plx_factor, jd = gm.fetch_table_element(
        ["scanAngle[rad]", "parallaxFactorAlongScan",
         "ObservationTimeAtBarycentre[BarycentricJulianDateInTCB]"], tab)
    t_yr = (np.asarray(jd) - DR_REF_JD[data_release]) / 365.25
    return Epochs(np.asarray(t_yr, float), np.asarray(psi, float), np.asarray(plx_factor, float))


def synthetic(n_transits: int = 70, baseline_yr: float = 5.5,
              rng: np.random.Generator | None = None) -> Epochs:
    """Crude stand-in for the scanning law, FOR UNIT TESTS ONLY.

    Transits come in pairs ~2 h apart within ~30 visibility periods, with random
    scan angles and a parallax factor that varies annually. Do not use it for
    science: real results must use :func:`from_gaiamock`.
    """
    rng = np.random.default_rng() if rng is None else rng
    n_vis = max(n_transits // 2, 13)
    centres = np.sort(rng.uniform(0, baseline_yr, n_vis)) - baseline_yr / 2
    t = np.repeat(centres, 2)[:n_transits]
    t = t + np.tile([0.0, 2.0 / 24 / 365.25], n_vis)[:n_transits]
    psi = np.repeat(rng.uniform(0, 2 * np.pi, n_vis), 2)[:n_transits]
    plx_factor = np.sin(2 * np.pi * t + psi)
    return Epochs(t, psi, plx_factor)
