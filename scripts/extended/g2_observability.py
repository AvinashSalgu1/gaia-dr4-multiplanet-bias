"""Which companions DR4 measures as full orbits, which only as accelerations, and where RVs help.

Needs gaiamock (real scanning law); about 15-30 minutes on a 12-core machine with the defaults:

    python scripts/extended/g2_observability.py            # [--n-inj 4] [-j -1]
    python scripts/extended/g2_observability.py --plot-only

For a Sun-like star (1.0 M_sun) and an M dwarf (0.5 M_sun) at 25, 50 and 100 pc, companions
on a grid of mass (0.1-30 M_J) and semi-major axis (0.3-30 au) are injected into simulated
DR4 epoch astrometry (random sky position, orientation, phase and Kipping eccentricity) and
fitted with the DR4 models. Each injection is classed as:

* orbit: an orbit (alone or with an acceleration) passes the Lammers & Winn detection test
  (delta chi2 > 50 and a constrained period; without their 13 M_J limit, since the grid
  extends to brown dwarfs) with significance >= 5 - the inclination and the node of the companion's orbit are measured;
* poorly constrained orbit: the same with significance < 5;
* acceleration only: no orbit, but an acceleration or jerk model is preferred over a single
  star by the BIC - no inclination; mass and separation are degenerate;
* not detected.

Overlaid (analytic): RV semi-amplitude of 3 and 10 m/s for sin i = 0.87, and the 20-yr
period that long RV surveys cover. Hipparcos data are not yet in the joint fit; their reach
will be added when the joint DR4 + Hipparcos + RV fit is implemented (Study 2).
Output in results/extended/observability/: summary.txt, observability.csv, observability.png.
"""
from __future__ import annotations

import argparse
import csv

import numpy as np

import common as C

STARS = [(1.0, 25.0), (1.0, 50.0), (1.0, 100.0), (0.5, 25.0), (0.5, 50.0), (0.5, 100.0)]
M_GRID = np.geomspace(0.1, 30.0, 10)
A_GRID = np.geomspace(0.3, 30.0, 12)
CATS = ["orbit", "poorly constrained orbit", "acceleration only", "not detected"]
SNR_FIT = 0.3


def one(args):
    m_star, dist, m, a_au, k, seed, synthetic = args
    from gaiabias import epochs as E
    from gaiabias import population as POP
    from gaiabias.fitting import fit_5par, fit_orbit_trend, fit_trend_only
    from gaiabias.orbits import Companion, draw_ecc_kipping, reflex_a0_mas
    from gaiabias.simulate import Star, sigma_al_per_transit, simulate

    rng = np.random.default_rng([seed, int(1000 * m_star), int(dist), int(1000 * m), int(1000 * a_au), k])
    G = float(POP.abs_G(np.array([m_star]))[0]) + 5 * np.log10(dist / 10.0)
    P = float(POP.period_days(a_au, m_star, m))
    if synthetic:
        ep = E.synthetic(80, rng=rng)
    else:
        ra, dec = rng.uniform(0, 360), np.degrees(np.arcsin(rng.uniform(-1, 1)))
        ep = E.from_gaiamock(ra, dec, "dr4", rng=rng)
    plx = 1000.0 / dist
    star = Star(parallax_mas=plx, m_star_msun=m_star, G_mag=G, pmra=rng.normal(0, 40), pmdec=rng.normal(0, 40))
    c = Companion(m, P, float(draw_ecc_kipping(rng)), np.arccos(rng.uniform(-1, 1)), rng.uniform(0, np.pi),
                  rng.uniform(0, 2 * np.pi), rng.uniform(0, P))
    sim = simulate(ep, star, [c], add_noise=False)
    sig = float(sigma_al_per_transit(G))
    obs = sim.obs + sig * rng.standard_normal(len(ep))
    if G < 13:
        obs += rng.uniform(0, 0.04) * rng.standard_normal(len(ep))
    n = len(ep)
    c5 = fit_5par(ep, obs, sim.err)[0]
    c7 = fit_trend_only(ep, obs, sim.err, 2)[0]
    c9 = fit_trend_only(ep, obs, sim.err, 4)[0]
    trend = (c5 - c7 > 2 * np.log(n)) or (c5 - c9 > 4 * np.log(n))
    if reflex_a0_mas(m_star, m, P, plx) / sig >= SNR_FIT:
        best = None
        for nt in (0, 2):
            try:
                f = fit_orbit_trend(ep, obs, sim.err, m_star, n_trend=nt)
            except Exception:
                continue
            # Lammers & Winn test without the 13 M_J mass limit (the grid includes brown dwarfs)
            sP = f.sigma_P
            if f.delta_chi2 > 50 and np.isfinite(sP) and 0 < sP < f.P_days and \
                    (f.P_days + sP) / (f.P_days - sP) < 1.5:
                if best is None or f.significance > best.significance:
                    best = f
        if best is not None:
            return 0 if best.significance >= 5 else 1
    return 2 if trend else 3


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--n-inj", type=int, default=4)
    p.add_argument("-j", "--jobs", type=int, default=-1)
    p.add_argument("--seed", type=int, default=31)
    p.add_argument("--plot-only", action="store_true")
    p.add_argument("--synthetic", action="store_true", help="synthetic epochs (code test only)")
    p.add_argument("--quick", action="store_true", help="coarse grid, 1 injection (test)")
    a = p.parse_args()
    out = C.out_dir("observability")
    path = out / "observability.csv"
    mg, ag = (M_GRID[::3], A_GRID[::3]) if a.quick else (M_GRID, A_GRID)
    if not a.plot_only:
        from joblib import Parallel, delayed
        n_inj = 1 if a.quick else a.n_inj
        tasks = [(ms, d, m, aa, k, a.seed, a.synthetic) for ms, d in STARS for m in mg for aa in ag
                 for k in range(n_inj)]
        print(f"{len(tasks)} injections ...", flush=True)
        t = C.Timer()
        res = Parallel(n_jobs=a.jobs, verbose=2)(delayed(one)(x) for x in tasks)
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["m_star", "dist_pc", "m_MJ", "a_au", "k", "category"])
            for x, r in zip(tasks, res):
                w.writerow([x[0], x[1], f"{x[2]:.4g}", f"{x[3]:.4g}", x[4], CATS[r]])
        print(f"done in {t}")
    rows = list(csv.DictReader(open(path)))
    mg = np.unique([float(r["m_MJ"]) for r in rows])
    ag = np.unique([float(r["a_au"]) for r in rows])
    frac = {}
    for ms, d in STARS:
        F = np.zeros((len(CATS), len(mg), len(ag)))
        for r in rows:
            if float(r["m_star"]) == ms and float(r["dist_pc"]) == d:
                i, j = np.searchsorted(mg, float(r["m_MJ"])), np.searchsorted(ag, float(r["a_au"]))
                F[CATS.index(r["category"]), i, j] += 1
        frac[(ms, d)] = F / np.maximum(F.sum(0, keepdims=True), 1)
    L = ["=== WHAT DR4 MEASURES: FULL ORBIT, ACCELERATION, OR NOTHING ===",
         f"{len(rows)} injections; per star and distance: smallest mass with >= 50% full-orbit (significance >= 5) "
         "detections, by semi-major axis.", ""]
    for (ms, d), F in frac.items():
        cells = []
        for j, aa in enumerate(ag):
            ok = np.nonzero(F[0, :, j] >= 0.5)[0]
            acc = np.nonzero(F[0, :, j] + F[1, :, j] + F[2, :, j] >= 0.5)[0]
            cells.append(f"{aa:.2g} au: orbit {mg[ok[0]]:.2g}" if len(ok) else f"{aa:.2g} au: orbit -")
            cells[-1] += f" / any {mg[acc[0]]:.2g}" if len(acc) else " / any -"
        L.append(f"  {ms:g} M_sun at {d:g} pc: " + "; ".join(cells))
    L += ["", "'orbit' = full orbit with significance >= 5 (inclination and node measured); 'any' = orbit, poorly "
          "constrained orbit or acceleration-only detection (M_J)."]
    C.write_summary(out / "summary.txt", L)

    plt, ps = C.setup_plots()
    from gaiabias import population as POP
    fig, axes = plt.subplots(2, 3, figsize=(16, 9.5), sharex=True, sharey=True, constrained_layout=True)
    for k, ((ms, d), F) in enumerate(frac.items()):
        ax = axes[k // 3, k % 3]
        im = ax.pcolormesh(ag, mg, F[0], cmap="Blues", vmin=0, vmax=1, shading="nearest")
        any_det = F[0] + F[1] + F[2]
        ax.contour(ag, mg, any_det, levels=[0.5], colors=[ps.C_ACCEL], linewidths=2)
        ax.contour(ag, mg, F[0] + F[1], levels=[0.5], colors=[ps.C_WRONG], linewidths=1.5, linestyles="--")
        aa = np.geomspace(ag[0], ag[-1], 100)
        for K, ls in ((3.0, ":"), (10.0, "-.")):
            P_yr = np.sqrt(aa ** 3 / ms)
            ax.plot(aa, K / (28.4329 * 0.87) * P_yr ** (1 / 3) * ms ** (2 / 3), ls, color=ps.INK, lw=1,
                    label=f"RV K = {K:g} m/s" if k == 0 else None)
        a20 = (20.0 ** 2 * ms) ** (1 / 3)
        ax.axvline(a20, color="0.5", lw=1, ls="--")
        ax.text(a20 * 1.05, mg[0] * 1.3, "P = 20 yr", fontsize=8, color="0.4", rotation=90)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(f"{ms:g} M_sun star at {d:g} pc (G = "
                     f"{float(POP.abs_G(np.array([ms]))[0]) + 5 * np.log10(d / 10):.1f})", fontsize=10.5)
        if k // 3 == 1:
            ax.set_xlabel("companion semi-major axis (au)")
        if k % 3 == 0:
            ax.set_ylabel("companion mass (M_J)")
    fig.colorbar(im, ax=axes, label="share measured as a full orbit (significance >= 5)", shrink=0.6)
    h = [plt.Line2D([], [], color=ps.C_ACCEL, lw=2), plt.Line2D([], [], color=ps.C_WRONG, lw=1.5, ls="--")]
    hl, ll = axes[0, 0].get_legend_handles_labels()
    axes[0, 0].legend(h + hl, ["50%: any DR4 detection (orbit or acceleration)", "50%: any orbit"] + ll,
                      fontsize=7.5, loc="lower right")
    fig.suptitle("DR4 alone: full orbits (inclination measured) vs accelerations only; RV lines for sin i = 0.87",
                 fontsize=12)
    fig.savefig(out / "observability.png", dpi=150, bbox_inches="tight")
    print(f"saved {out / 'observability.png'}")


if __name__ == "__main__":
    main()
