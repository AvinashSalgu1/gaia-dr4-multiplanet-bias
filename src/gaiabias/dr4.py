"""Real Gaia DR4 epoch astrometry as input to the project's orbit fits.

The DR4 epoch-astrometry table has one row per field-of-view transit with per-CCD arrays
(DR4 draft data model, Sec. 3.1; checked on the June 2026 pre-release). The column names
below are the data-model ones. Any file that provides them works: the pre-release VOTable
(read without astropy by dr4_epoch.py), or ECSV, CSV or FITS files read with astropy (the
format of the bulk DR4 files is not yet published; only `_read_astropy` may need changing).

The fits and the simulations of this project work per transit: the CCD measurements of a
transit that AGIS used (used_by_agis_al, which excludes the SkyMapper CCD) are combined into
their inverse-variance weighted mean, with sigma_transit = (sum 1/sigma_CCD^2)^(-1/2). The
pre-release showed that CCDs within a transit are independent (correlation -0.001 +/- 0.007),
which justifies this. sigma_CCD is centroid_pos_error_al; for the single-star model the data
model adds agis_source_excess_noise in quadrature (excess=True); orbit fits use the formal
errors (excess=False), since the excess noise belongs to the single-star solution.

Times are TCB nanoseconds since J2010.0 plus the barycentric correction; they are returned in
years from the DR4 reference epoch J2017.5, the convention of gaiabias.epochs.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .epochs import Epochs

NS_PER_YR = 1e9 * 86400 * 365.25
J2010_TO_J2017P5_YR = 7.5


def load_ccd_table(path) -> dict:
    """One entry per CCD measurement (dict of numpy arrays), from any supported file."""
    path = Path(path)
    if path.suffix.lower() in (".xml", ".vot", ".votable"):
        from .dr4_epoch import ccd_table
        return ccd_table(path)
    return _read_astropy(path)


def _read_astropy(path: Path) -> dict:
    from astropy.table import Table
    t = Table.read(path)

    def arr(row, name):
        v = row[name]
        if isinstance(v, str):                       # CSV/ECSV array columns may come as "[a, b, ...]"
            v = np.array([float(x) for x in v.strip("[]()").replace(",", " ").split()]) if v.strip("[]() ") else []
        return np.atleast_1d(np.asarray(v))

    rows = {k: [] for k in ("source_id", "transit", "ccd", "t_ns", "psi_deg", "pf_al", "w", "err", "ipd_err",
                            "excess_src", "used", "g_mag", "ra0", "dec0")}
    for r in t:
        w = arr(r, "centroid_pos_al").astype(float)
        tt, psi = arr(r, "obs_time_tcb").astype(float), arr(r, "scan_pos_angle").astype(float)
        err, ipd = arr(r, "centroid_pos_error_al").astype(float), arr(r, "ipd_error_al").astype(float)
        used = arr(r, "used_by_agis_al")
        for c in range(len(w)):
            rows["source_id"].append(int(r["source_id"]))
            rows["transit"].append(int(r["transit_id"]))
            rows["ccd"].append(c)
            rows["t_ns"].append(tt[c] + float(r["obs_time_bary_corr"]))
            rows["psi_deg"].append(psi[c])
            rows["pf_al"].append(float(r["parallax_factor_al"]))
            rows["w"].append(w[c])
            rows["err"].append(err[c])
            rows["ipd_err"].append(ipd[c] if c < len(ipd) else np.nan)
            rows["excess_src"].append(float(r["agis_source_excess_noise"]))
            rows["used"].append(bool(used[c]) if c < len(used) else False)
            rows["g_mag"].append(float(r["g_mag"]))
            rows["ra0"].append(float(r["ra0"]))
            rows["dec0"].append(float(r["dec0"]))
    return {k: np.array(v) for k, v in rows.items()}


def source_ids(d: dict) -> np.ndarray:
    return np.unique(d["source_id"])


def source_ccds(d: dict, sid: int, excess: bool = False):
    """(Epochs, along-scan positions [mas], errors [mas], transit ids) of the used CCDs of one source."""
    t_yr = d["t_ns"] / NS_PER_YR - J2010_TO_J2017P5_YR
    m = (d["source_id"] == sid) & d["used"].astype(bool) & np.isfinite(d["w"]) & np.isfinite(t_yr) & (d["err"] > 0)
    err = d["err"][m]
    if excess:
        err = np.sqrt(err ** 2 + np.nan_to_num(d["excess_src"][m]) ** 2)
    return Epochs(t_yr[m], np.radians(d["psi_deg"][m]), d["pf_al"][m]), d["w"][m], err, d["transit"][m]


def source_transits(d: dict, sid: int, excess: bool = False, min_ccd: int = 1):
    """Per-transit (Epochs, along-scan position [mas], error [mas], number of CCDs, G) of one source.

    Each transit is the inverse-variance weighted mean of its used CCDs; time, scan angle and
    parallax factor are those of the same weighted mean (they change by < 1 s and negligibly in
    angle across the CCDs of a transit)."""
    ep, w, err, tr = source_ccds(d, sid, excess)
    ids, inv = np.unique(tr, return_inverse=True)
    wt = 1.0 / err ** 2
    S = np.bincount(inv, weights=wt)
    mean = lambda x: np.bincount(inv, weights=wt * x) / S  # noqa: E731
    # scan angle: weighted circular mean
    psi = np.arctan2(np.bincount(inv, weights=wt * np.sin(ep.psi)), np.bincount(inv, weights=wt * np.cos(ep.psi)))
    n = np.bincount(inv)
    k = n >= min_ccd
    g = d["g_mag"][d["source_id"] == sid]
    G = float(np.nanmedian(g)) if len(g) else np.nan
    return (Epochs(mean(ep.t_yr)[k], psi[k], mean(ep.plx_factor)[k]), mean(w)[k], (1.0 / np.sqrt(S))[k], n[k], G)


def fit_source(ep: Epochs, y, err, m_star: float) -> dict:
    """All DR4 orbit models of step 2 on one source, the BIC choice, and the reported planets.

    Returns flat columns: n_transits, chi2 of every model (as in step 2), the adopted model
    (index and name in outcomes.BIC_MODELS) and, for each planet the adopted model reports, its
    period [d], mass [M_J, dark-companion formula for m_star] and significance a0/sigma_a0."""
    from .step2 import analyses as A
    from .step2 import outcomes as O
    out = dict(n_transits=len(ep))
    out.update(A.run(ep, y, err, m_star, "w", set(A.FITTED)))
    d = {k: np.array([v]) for k, v in out.items()}
    k = int(O.bic_choice(d, "w")[0])
    out["orbit_P"], out["orbit_m"], out["orbit_e"] = out.get("one_w_P"), out.get("one_w_m"), out.get("one_w_e")
    out["model_index"] = k
    out["model"] = O.BIC_MODELS[k] if k >= 0 else "none"
    planets = [(P[0], m[0], s[0]) for mask, P, m, s in O.reported(d, "best", "w") if mask[0]]
    out["n_planets"] = len(planets)
    for i, (P, m, s) in enumerate(sorted(planets)[:2], 1):
        out[f"planet{i}_P"], out[f"planet{i}_m"], out[f"planet{i}_sig"] = P, m, s
    return out
