"""Check our assumptions against REAL Gaia DR4 epoch astrometry (June 2026 pre-release, 12 sources).

1. Noise level: the formal per-CCD along-scan errors vs the gaiamock / Holl et
   al. (2023a) per-CCD table used in all our simulations.
2. Are the CCD measurements within one field-of-view transit independent?
   Our simulations (like gaiamock and Lammers & Winn, footnote 5) assume so,
   i.e. sigma_transit = sigma_CCD / sqrt(N_CCD). For each single star / quasar we
   fit the 5-parameter model to all CCDs, then compare the scatter of each
   transit's mean residual with what independent CCD errors predict
   (z = mean residual / its formal error; rms(z) = 1 if independent,
   rms(z)^2 = 1 + (n - 1) rho for a within-transit correlation rho).
3. Our orbit fitter on real data: one-companion fits to the three sources with
   published orbits (HD 114762, Gaia-4, Gaia BH3), and the acceleration and
   two-planet models of step 2 on Gaia-4.

Weights follow the DR4 data model (Sec. 3.1): CCDs with used_by_agis_al only
(this excludes the SkyMapper CCD, whose errors are underestimated);
sigma = centroid_pos_error_al for source models other than the single-star
model (agis_source_excess_noise belongs to the single-star AGIS solution).

Input: the pre-release XML (data/dr4_prerelease/*EPOCH_ASTROMETRY*.xml), from
https://www.cosmos.esa.int/web/gaia/dr4-prerelease
Run:   python scripts/check_dr4_prerelease.py      (seconds; no gaiamock needed)
"""
from __future__ import annotations

import glob
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
OUT = ROOT / "results" / "dr4_prerelease"

SOURCES = {   # source_id: (label, kind, host mass [M_sun] for orbit fits)
    4181040337841125632: ("HD 183633 (G 9)", "single", None),
    3926186255616949504: ("G 19 star", "single", None),
    435469040545191680: ("G 14 star", "single", None),
    1663617687609809280: ("LSPM J1324+6124", "single", None),
    2309425390592896: ("G 14 star b", "single", None),
    20694084440761600: ("G 14 star c", "single", None),
    2237987199365376: ("quasar 1", "single", None),
    60730287810150016: ("quasar 2", "single", None),
    10973744521070720: ("quasar 3", "single", None),
    3937211745905473024: ("HD 114762", "orbit", 0.83),
    1457486023639239296: ("Gaia-4", "orbit", 0.64),
    4318465066420528000: ("Gaia BH3", "orbit", 0.76),
}


def find_input():
    hits = sorted(glob.glob(str(ROOT / "data" / "dr4_prerelease" / "**" / "*EPOCH_ASTROMETRY*.xml"), recursive=True))
    if len(sys.argv) > 1:
        hits = [sys.argv[1]]
    if not hits:
        raise SystemExit("put the pre-release XML in data/dr4_prerelease/ (or pass its path)")
    return hits[0]


def per_source(d, s):
    from gaiabias.epochs import Epochs
    t_yr = d["t_ns"] / 1e9 / 86400 / 365.25 - 7.5        # ns since J2010.0 (TCB) -> years from J2017.5
    m = (d["source_id"] == s) & d["used"] & np.isfinite(d["w"]) & np.isfinite(t_yr)
    return Epochs(t_yr[m], np.radians(d["psi_deg"][m]), d["pf_al"][m]), d["w"][m], d["err"][m], \
        d["transit"][m], float(np.nanmedian(d["g_mag"][d["source_id"] == s]))


def transit_test(ep, w, err, tr):
    """rms of standardised transit-mean residuals after a 5-parameter fit, and within-transit chi2/dof."""
    from gaiabias.fitting import _design_5par, fit_5par
    chi2, mu = fit_5par(ep, w, err)
    r = w - _design_5par(ep) @ mu
    z, within, nccd = [], [], []
    for x in np.unique(tr):
        k = tr == x
        if k.sum() < 4:
            continue
        wt = 1 / err[k] ** 2
        mean = np.sum(wt * r[k]) / wt.sum()
        z.append(mean * np.sqrt(wt.sum()))
        within.append(np.sum(wt * (r[k] - mean) ** 2) / (k.sum() - 1))
        nccd.append(k.sum())
    z = np.array(z)
    return dict(chi2_dof=chi2 / (len(w) - 5), z2=np.mean(z ** 2), z2_err=np.std(z ** 2) / np.sqrt(len(z)),
                within=np.median(within), n_tr=len(z), n_ccd=float(np.mean(nccd)), plx=mu[4])


def main():
    from gaiabias.dr4_epoch import ccd_table
    from gaiabias.fitting import fit_orbit_trend, fit_two_planets
    from gaiabias.simulate import _G_VALS, _SIGMA_CCD

    path = find_input()
    d = ccd_table(path)
    OUT.mkdir(parents=True, exist_ok=True)
    L = ["=== DR4 PRE-RELEASE CHECK ===", f"input: {Path(path).name}, "
         f"{len(np.unique(d['source_id']))} sources, {len(d['w'])} CCD measurements", "",
         "[1] noise and [2] within-transit independence (single stars and quasars, 5-parameter fit):",
         f"  {'source':<20}{'G':>6}{'CCD err':>9}{'gaiamock':>10}{'chi2/dof':>10}{'within':>8}"
         f"{'<z^2>':>12}{'parallax':>10}"]
    rows = []
    for s, (lab, kind, _) in SOURCES.items():
        if kind != "single" or s not in d["source_id"]:
            continue
        ep, w, err, tr, G = per_source(d, s)
        r = transit_test(ep, w, err, tr)
        gm = float(np.interp(G, _G_VALS, _SIGMA_CCD))
        rows.append((G, float(np.median(err)), gm, r))
        L.append(f"  {lab:<20}{G:6.2f}{np.median(err):9.3f}{gm:10.3f}{r['chi2_dof']:10.2f}{r['within']:8.2f}"
                 f"{r['z2']:7.2f}±{r['z2_err']:.2f}{r['plx']:10.3f}")
    z2 = np.array([r[3]["z2"] for r in rows])
    z2e = np.array([r[3]["z2_err"] for r in rows])
    w_ = 1 / z2e ** 2
    z2m, z2me = np.sum(w_ * z2) / w_.sum(), 1 / np.sqrt(w_.sum())
    nbar = np.mean([r[3]["n_ccd"] for r in rows])
    rho, rho_e = (z2m - 1) / (nbar - 1), z2me / (nbar - 1)
    L += [f"  pooled <z^2> = {z2m:.3f} ± {z2me:.3f}  ->  within-transit error correlation rho = {rho:+.3f} ± {rho_e:.3f}"
          f" (n = {nbar:.1f} CCDs per transit)",
          f"  => sigma_transit / (sigma_CCD/sqrt(n)) = {np.sqrt(z2m):.3f};  2-sigma upper limit on rho: {rho + 2 * rho_e:.3f}",
          "  (rho = 0: CCDs independent, as gaiamock and L&W assume; rho = 1 would make transits sqrt(n) ~ 3x noisier)",
          "", "[3] orbit fits on real data (one companion, all CCDs, numpy optimiser of step 2):"]
    orbit_rows = []
    for s, (lab, kind, mstar) in SOURCES.items():
        if kind != "orbit" or s not in d["source_id"]:
            continue
        ep, w, err, tr, G = per_source(d, s)
        f = fit_orbit_trend(ep, w, err, mstar, n_trend=0, P_max=2e4)
        orbit_rows.append((lab, ep, w, err, f, tr))
        L.append(f"  {lab:<11} P {f.P_days:8.1f} ± {f.sigma_P:6.1f} d   e {f.ecc:.3f}   a0 {f.a0_mas:7.3f} mas   "
                 f"inc {np.degrees(f.inc):5.1f} deg   parallax {f.parallax:7.3f}   companion {f.m_mjup:9.2f} M_J "
                 f"(dark-companion formula, M* = {mstar})   chi2/dof {f.chi2 / (len(w) - 12):.2f}")
        if lab == "Gaia-4":
            f2 = fit_orbit_trend(ep, w, err, mstar, n_trend=2)
            t = fit_two_planets(ep, w, err, mstar, one=f)
            lnN = np.log(len(w))
            L.append(f"     Gaia-4, step-2 models: + acceleration improves chi2 by {f.chi2 - f2.chi2:.1f} (BIC needs "
                     f"> {2 * lnN:.1f}); two planets by {t.gain_over_one:.1f} (BIC needs > {7 * lnN:.1f}); "
                     f"second-planet candidate P {t.P[1]:.1f} d, {t.m_mjup[1]:.1f} M_J")
    L.append("  published: HD 114762 b P ~ 84 d; Gaia-4 b ~ 570 d, ~12 M_J; Gaia BH3 ~ 11.6 yr, e ~ 0.73, ~33 M_sun "
             "(check against the papers before quoting)")
    print("\n".join(L))
    (OUT / "summary.txt").write_text("\n".join(L) + "\n")
    plot(rows, orbit_rows, _G_VALS, _SIGMA_CCD)


def plot(rows, orbit_rows, G_vals, sig_ccd):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from gaiabias import plotstyle as ps
    from gaiabias.fitting import _design_9par

    ps.apply()
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.8), constrained_layout=True)
    a = ax[0]
    g = np.linspace(6, 20, 200)
    a.plot(g, np.interp(g, G_vals, sig_ccd), color=ps.MUTED, lw=2, label="per-CCD error assumed in our simulations")
    a.scatter([r[0] for r in rows], [r[1] for r in rows], color=ps.C_TARGET, zorder=3,
              label="real DR4 formal per-CCD error (median)")
    a.set_yscale("log")
    a.set_xlabel("G magnitude")
    a.set_ylabel("along-scan error per CCD (mas)")
    a.set_title("(a) Noise level: real DR4 vs our model")
    a.legend(fontsize=8.5)
    a = ax[1]
    y = [r[3]["z2"] for r in rows]
    e = [r[3]["z2_err"] for r in rows]
    a.errorbar(range(len(rows)), y, e, fmt="o", color=ps.C_TARGET, capsize=3)
    a.axhline(1, color=ps.INK, lw=1, label="CCDs independent (our assumption)")
    a.axhline(1 + 7 * 0.1, color=ps.C_WRONG, ls="--", lw=1, label="10% correlation within a transit")
    a.set_xticks(range(len(rows)), [f"G {r[0]:.0f}" for r in rows], fontsize=8)
    a.set_ylabel("scatter of transit means ÷ expected\n(mean z², 1 = independent)")
    a.set_title("(b) Are the CCDs in one transit independent?")
    a.set_ylim(0.5, 1.85)
    a.legend(fontsize=8.5, loc="lower right")
    a = ax[2]
    for lab, ep, w, err, f, tr in orbit_rows:
        if lab != "Gaia-4":
            continue
        M = _design_9par(ep, f.P_days, f.phi_p, f.ecc)
        pred = M[:, 5:] @ f.mu[5:9]                      # fitted orbit, projected on the scan direction
        meas = w - M[:, :5] @ f.mu[:5]                   # data minus position, motion and parallax
        xs, ys, es = [], [], []
        for x in np.unique(tr):                          # one point per transit (weighted mean of its CCDs)
            k = tr == x
            wt = 1 / err[k] ** 2
            xs.append(np.mean(pred[k]))
            ys.append(np.sum(wt * meas[k]) / wt.sum())
            es.append(1 / np.sqrt(wt.sum()))
        a.errorbar(xs, ys, es, fmt="o", ms=3.5, color=ps.C_TARGET, ecolor=ps.C_NONE, lw=0.8,
                   label=f"one point per transit ({len(xs)})")
        lim = 1.15 * max(np.max(np.abs(xs)), np.max(np.abs(ys)))
        a.plot([-lim, lim], [-lim, lim], color=ps.INK, lw=1, label="perfect agreement")
        a.set_xlim(-lim, lim)
        a.set_ylim(-lim, lim)
        a.set_title(f"(c) Gaia-4: measured wobble vs our fitted orbit\n(P = {f.P_days:.0f} d, {f.m_mjup:.1f} M_J)")
    a.set_xlabel("wobble predicted by the fitted orbit (mas, along scan)")
    a.set_ylabel("measured wobble (mas, along scan)")
    a.legend(fontsize=8.5, loc="upper left")
    fig.suptitle("Checking the simulation assumptions on real Gaia DR4 epoch astrometry (pre-release, 12 sources)",
                 fontsize=12)
    out = OUT / "dr4_prerelease_check.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")


if __name__ == "__main__":
    main()
