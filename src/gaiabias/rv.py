"""Radial-velocity follow-up of Gaia planets (step 4): signal sizes, instrument precision, epochs needed.

All instrument numbers are deliberately simple, documented assumptions (change them here):

* RV semi-amplitude K = 28.4329 m/s * m sin i / M_J * (P / 1 yr)^(-1/3) * (M_star / M_sun)^(-2/3)
  / sqrt(1 - e^2)  (planet mass << star mass).
* V magnitude from Gaia G and BP-RP with the Gaia DR3 relation (Riello et al. 2021):
  G - V = -0.02704 + 0.01424 x - 0.2156 x^2 + 0.01426 x^3, x = BP - RP. BP - RP is estimated from
  M_G with the main-sequence locus of the step 3 selection, M_G - 2.9 (BP - RP) = 2.45 (the
  middle of L&W's Eq. 38 band), when the colour is not available.
* Photon-noise precision per exposure scales as 10^(0.2 (V - 8)) (constant exposure time):
    HARPS (3.6 m, La Silla):  1.0 m/s at V = 8 in 900 s, instrumental floor 0.5 m/s;
    FEROS (2.2 m, La Silla):  5.0 m/s at V = 8 in 1200 s, instrumental floor 4 m/s.
  Stellar jitter added in quadrature: 2 m/s (FGK, M > 0.6 M_sun), 3 m/s (M dwarfs).
  Each epoch costs exposure + 300 s overhead.
* Epochs needed for a 5-sigma result with N evenly spread epochs of precision sigma:
    known period (confirm a Gaia orbit): sigma_K = sigma sqrt(2/N)       -> N = 2 (5 sigma / K)^2
    orbit with period shorter than the programme: same, plus N >= 3 per free parameter (>= 8)
    slow companion seen as a trend over a baseline B: sigma_(Delta v) = sigma sqrt(12/N)
                                                     -> N = 12 (5 sigma / Delta v)^2,
    with Delta v = min(2K, 2 pi K B / P) x 2/pi (the average change of a sinusoid over B).
  At least 6 epochs are always required.
"""
from __future__ import annotations

import numpy as np

K_CONST = 28.4329          # m/s for 1 M_J, 1 yr, 1 M_sun
INSTRUMENTS = {
    "HARPS": dict(sigma8=1.0, floor=0.5, t_exp=900.0),
    "FEROS": dict(sigma8=5.0, floor=4.0, t_exp=1200.0),
}
OVERHEAD_S = 300.0
N_MIN = 6
SNR_REQ = 5.0


def semi_amplitude(m_mj, P_yr, m_star, e=0.0, sin_i=1.0):
    """RV semi-amplitude in m/s."""
    m_mj, P_yr, m_star, e, sin_i = (np.asarray(x, float) for x in (m_mj, P_yr, m_star, e, sin_i))
    return K_CONST * m_mj * sin_i * P_yr ** (-1 / 3) * m_star ** (-2 / 3) / np.sqrt(1 - np.minimum(e, 0.99) ** 2)


def bp_rp_from_MG(M_G):
    return (np.asarray(M_G, float) - 2.45) / 2.9


def v_mag(G, bp_rp):
    x = np.asarray(bp_rp, float)
    return np.asarray(G, float) - (-0.02704 + 0.01424 * x - 0.2156 * x ** 2 + 0.01426 * x ** 3)


def rv_precision(V, m_star, instrument="HARPS"):
    """Per-epoch RV uncertainty in m/s, including jitter."""
    p = INSTRUMENTS[instrument]
    phot = np.hypot(p["sigma8"] * 10 ** (0.2 * (np.asarray(V, float) - 8)), p["floor"])
    jit = np.where(np.asarray(m_star, float) > 0.6, 2.0, 3.0)
    return np.hypot(phot, jit)


def epochs_known_period(K, sigma):
    """Epochs to measure K at SNR_REQ when the period is known (confirming a Gaia orbit)."""
    with np.errstate(divide="ignore"):
        return np.maximum(N_MIN, np.ceil(2 * (SNR_REQ * sigma / np.asarray(K, float)) ** 2))


def trend_change(K, P_yr, baseline_yr):
    """Typical RV change over the programme for a companion of semi-amplitude K and period P."""
    K, P_yr = np.asarray(K, float), np.asarray(P_yr, float)
    return np.minimum(2 * K, 2 * np.pi * K * baseline_yr / P_yr) * 2 / np.pi


def epochs_companion(K, P_yr, sigma, baseline_yr):
    """Epochs to detect a companion: as an orbit if its period fits in the programme, else as a trend."""
    K, P_yr = np.asarray(K, float), np.asarray(P_yr, float)
    with np.errstate(divide="ignore"):
        orbit = np.maximum(8, np.ceil(2 * (SNR_REQ * sigma / K) ** 2))
        dv = trend_change(K, P_yr, baseline_yr)
        trend = np.maximum(N_MIN, np.ceil(12 * (SNR_REQ * sigma / dv) ** 2))
    return np.where(P_yr <= baseline_yr, orbit, trend)


def hours(n_epochs, instrument="HARPS"):
    return np.asarray(n_epochs, float) * (INSTRUMENTS[instrument]["t_exp"] + OVERHEAD_S) / 3600.0
