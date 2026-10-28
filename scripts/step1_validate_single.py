"""Step 1: validate single-planet detection against Lammers & Winn (2026).

Injects ONE planet per star over a range of SNR_1 = a0/sigma_fov and orbital
period, fits it with gaiamock's one-companion optimiser, and measures

  (a) detection probability vs SNR_1, for several P/T bins;
  (b) detection probability vs P/T, for several SNR_1 bins;
  (c) Delta chi^2 vs SNR_1^2  (L&W find Delta chi^2 ~ 52 SNR_1^2);
  (d) how well P and m are recovered for detected planets.

T is the DR4 observing window (5.41 yr). L&W reference points, to compare with:
  - Delta chi^2 ~= 52 SNR_1^2  (detection cut Delta chi^2 > 50 -> SNR_1 ~ 1)
  - longest detectable periods ~4.0 yr in DR4 (~72% of the baseline)

Run:
    python scripts/step1_validate_single.py            # 10000 injections (~3 core-hours)
    python scripts/step1_validate_single.py -n 2000    # quick look
    python scripts/step1_validate_single.py --plot-only   # re-plot saved results

Outputs: results/step1_single.csv, results/step1_single.png, and a short
summary printed at the end.
"""
from __future__ import annotations

import argparse
import csv
import time
from pathlib import Path

import numpy as np

RESULTS = Path(__file__).resolve().parents[1] / "results"
CSV_PATH = RESULTS / "step1_single.csv"
PNG_PATH = RESULTS / "step1_single.png"
T_DR4_DAYS = 2458868.5 - 2456891.5   # gaiamock's DR4 window: 1977 d = 5.41 yr
MJ_MAX = 13.0

COLUMNS = ["seed", "ra", "dec", "parallax", "G", "m_star", "m_true", "P_true", "ecc_true",
           "inc_true", "a0_true", "snr1", "n_transits", "dchi2_true", "dchi2_fit",
           "P_fit", "sigma_P", "ecc_fit", "a0_fit", "m_fit", "detected"]


# ----------------------------------------------------------------------------- draw
def draw_system(rng: np.random.Generator):
    """One random star + planet with prescribed (SNR_1, P) and m_true < 13 M_J.

    SNR_1 and P are drawn log-uniform so every part of the plane is sampled;
    the planet mass is then whatever is needed at the drawn distance.
    """
    from gedi.orbits import companion_mass_from_a0
    from gedi.simulate import sigma_al_per_transit

    while True:
        snr1 = 10 ** rng.uniform(np.log10(0.1), np.log10(5.0))
        P = 10 ** rng.uniform(np.log10(50.0), np.log10(3.0 * T_DR4_DAYS))
        G = rng.uniform(8.0, 15.0)
        m_star = rng.uniform(0.4, 1.3)
        parallax = 10 ** rng.uniform(np.log10(2.0), np.log10(100.0))  # 10-500 pc
        a0 = snr1 * float(sigma_al_per_transit(G))
        m = companion_mass_from_a0(a0, P, m_star, parallax)
        if 0.1 < m < MJ_MAX:
            break
    ra = rng.uniform(0, 360)
    dec = np.degrees(np.arcsin(rng.uniform(-1, 1)))
    return dict(ra=ra, dec=dec, parallax=parallax, G=G, m_star=m_star, m_true=m,
                P_true=P, ecc_true=rng.uniform(0.0, 0.6), inc_true=np.arccos(rng.uniform(-1, 1)),
                Omega=rng.uniform(0, np.pi), omega=rng.uniform(0, 2 * np.pi),
                Tp=rng.uniform(0, P), pmra=rng.normal(0, 40), pmdec=rng.normal(0, 40), snr1=snr1)


# ----------------------------------------------------------------------------- one run
def run_one(seed: int) -> dict:
    from gedi import epochs as E
    from gedi.fitting import chi2_at, fit_one_companion, fit_5par, lammers_winn_detected
    from gedi.orbits import Companion
    from gedi.simulate import Star, simulate

    rng = np.random.default_rng(seed)
    s = draw_system(rng)
    ep = E.from_gaiamock(s["ra"], s["dec"], "dr4", rng=rng)
    star = Star(parallax_mas=s["parallax"], m_star_msun=s["m_star"], G_mag=s["G"],
                pmra=s["pmra"], pmdec=s["pmdec"])
    comp = Companion(s["m_true"], s["P_true"], s["ecc_true"], s["inc_true"],
                     s["Omega"], s["omega"], s["Tp"])
    sim = simulate(ep, star, [comp], rng=rng)

    chi2_null, _ = fit_5par(ep, sim.obs, sim.err)
    chi2_true, _ = chi2_at(ep, sim.obs, sim.err, comp.P_days, comp.phi_p, comp.ecc)
    try:
        fit = fit_one_companion(ep, sim.obs, sim.err, s["m_star"])
        out = dict(dchi2_fit=fit.delta_chi2, P_fit=fit.P_days, sigma_P=fit.sigma_P,
                   ecc_fit=fit.ecc, a0_fit=fit.a0_mas, m_fit=fit.m_mjup,
                   detected=int(lammers_winn_detected(fit)))
    except Exception:  # keep going; failed fits count as non-detections
        out = dict(dchi2_fit=np.nan, P_fit=np.nan, sigma_P=np.nan, ecc_fit=np.nan,
                   a0_fit=np.nan, m_fit=np.nan, detected=0)

    return dict(seed=seed, ra=s["ra"], dec=s["dec"], parallax=s["parallax"], G=s["G"],
                m_star=s["m_star"], m_true=s["m_true"], P_true=s["P_true"],
                ecc_true=s["ecc_true"], inc_true=s["inc_true"], a0_true=sim.a0_mas[0],
                snr1=s["snr1"], n_transits=len(ep), dchi2_true=chi2_null - chi2_true, **out)


def run_all(n: int, n_jobs: int, seed0: int):
    from joblib import Parallel, delayed

    seeds = np.random.SeedSequence(seed0).generate_state(n, dtype=np.uint32)
    t0 = time.perf_counter()
    rows = Parallel(n_jobs=n_jobs, verbose=5)(delayed(run_one)(int(s)) for s in seeds)
    print(f"{n} injections in {(time.perf_counter() - t0) / 60:.1f} min")
    RESULTS.mkdir(exist_ok=True)
    with open(CSV_PATH, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in COLUMNS})
    print(f"saved {CSV_PATH}")


# ----------------------------------------------------------------------------- analysis
def load():
    with open(CSV_PATH) as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) for r in rows]) for k in COLUMNS}


def binned_fraction(x, flag, edges):
    idx = np.digitize(x, edges) - 1
    frac, err, n = [], [], []
    for i in range(len(edges) - 1):
        m = idx == i
        k, N = flag[m].sum(), m.sum()
        p = k / N if N else np.nan
        frac.append(p)
        err.append(np.sqrt(p * (1 - p) / N) if N else np.nan)
        n.append(N)
    return np.array(frac), np.array(err), np.array(n)


def snr_at_half(snr, det, mask):
    """SNR_1 where detection probability crosses 50% (interpolated)."""
    edges = np.logspace(-1, np.log10(5), 21)
    f, _, n = binned_fraction(snr[mask], det[mask], edges)
    c = np.sqrt(edges[1:] * edges[:-1])
    ok = n >= 10
    c, f = c[ok], f[ok]
    above = np.where(f >= 0.5)[0]
    if len(above) == 0 or above[0] == 0:
        return np.nan
    i = above[0]
    return 10 ** np.interp(0.5, [f[i - 1], f[i]], np.log10([c[i - 1], c[i]]))


def analyse(make_plot: bool = True):
    d = load()
    det = d["detected"].astype(bool)
    x = d["P_true"] / T_DR4_DAYS
    snr = d["snr1"]
    n = len(snr)

    # ---- numbers to compare with Lammers & Winn
    short = x < 0.5
    strong = short & (snr > 1.0)   # signal-dominated, so noise doesn't inflate the ratio
    ratio_true = np.median(d["dchi2_true"][strong] / snr[strong] ** 2)
    good = strong & np.isfinite(d["dchi2_fit"])
    ratio_fit = np.median(d["dchi2_fit"][good] / snr[good] ** 2)
    s50_short = snr_at_half(snr, det, short)

    xe = np.array([0.0, 0.25, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0, 3.0])
    hi = snr > 3
    f_hi, _, n_hi = binned_fraction(x[hi], det[hi], xe)
    xc = 0.5 * (xe[1:] + xe[:-1])
    ok = (n_hi >= 10) & np.isfinite(f_hi)
    below = np.where(ok & (f_hi < 0.5) & (xc > 0.5))[0]
    if len(below):
        j = below[0]
        x50 = np.interp(0.5, [f_hi[j], f_hi[j - 1]], [xc[j], xc[j - 1]])
    else:
        x50 = np.nan

    # Lammers & Winn's own definition of the maximum period: where the MEAN
    # relative period error of recovered orbits reaches 20% (their Fig. 10).
    # Use all strong-signal fits, detected or not, as an injection-recovery test.
    fin = np.isfinite(d["P_fit"]) & (snr > 3)
    relerr = np.abs(d["P_fit"] / d["P_true"] - 1)
    pe = np.array([0.0, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.7, 2.0, 3.0])
    pc = 0.5 * (pe[1:] + pe[:-1])
    mean_err, med_err = [], []
    for lo, hi_ in zip(pe[:-1], pe[1:]):
        m = fin & (x >= lo) & (x < hi_)
        mean_err.append(np.mean(relerr[m]) if m.sum() >= 10 else np.nan)
        med_err.append(np.median(relerr[m]) if m.sum() >= 10 else np.nan)
    mean_err, med_err = np.array(mean_err), np.array(med_err)

    def first_cross(curve, level=0.2):
        idx = np.where(np.isfinite(curve) & (curve >= level))[0]
        if len(idx) == 0 or idx[0] == 0:
            return np.nan
        j = idx[0]
        return np.interp(level, [curve[j - 1], curve[j]], [pc[j - 1], pc[j]])

    x20_mean, x20_med = first_cross(mean_err), first_cross(med_err)

    rec = det & short
    dP = np.abs(d["P_fit"][rec] / d["P_true"][rec] - 1)
    dm = d["m_fit"][rec] / d["m_true"][rec] - 1
    opt_miss = np.mean((d["dchi2_true"] - d["dchi2_fit"] > 10)[np.isfinite(d["dchi2_fit"])])

    lines = [
        "=== STEP 1 SUMMARY ===",
        f"injections: {n}, detected: {det.sum()} ({det.mean():.1%}), "
        f"median transits/star: {np.median(d['n_transits']):.0f}",
        f"Delta chi2 / SNR1^2 (P < 0.5T, SNR1 > 1): true-orbit {ratio_true:.1f}, fitted {ratio_fit:.1f}   "
        f"[L&W DR4, Eq. 31: ~22]",
        f"SNR1 at 50% detection (P < 0.5T): {s50_short:.2f}",
        f"P/T where detection drops below 50% (SNR1 > 3): {x50:.2f} = {x50 * T_DR4_DAYS / 365.25:.1f} yr"
        f"   [our metric; not L&W's definition]",
        f"P/T where MEAN period error reaches 20% (SNR1 > 3): {x20_mean:.2f} = "
        f"{x20_mean * T_DR4_DAYS / 365.25:.1f} yr   (median error: {x20_med:.2f} = "
        f"{x20_med * T_DR4_DAYS / 365.25:.1f} yr)   [L&W definition: ~4.0 yr]",
        f"detected, P < 0.5T: median |dP/P| {np.median(dP):.3f}, "
        f"median dm/m {np.median(dm):+.3f}, 16-84% dm/m [{np.percentile(dm, 16):+.2f}, "
        f"{np.percentile(dm, 84):+.2f}]",
        f"optimiser missed true orbit by >10 in chi2: {opt_miss:.1%} of fits",
    ]
    print("\n".join(lines))
    (RESULTS / "step1_summary.txt").write_text("\n".join(lines) + "\n")

    if make_plot:
        plot(d, det, x, snr)


def plot(d, det, x, snr):
    import matplotlib
    import matplotlib.ticker
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    blues = plt.get_cmap("Blues")
    ink, muted = "#222222", "#888888"
    plt.rcParams.update({"axes.edgecolor": muted, "axes.labelcolor": ink, "xtick.color": ink,
                         "ytick.color": ink, "axes.spines.top": False, "axes.spines.right": False,
                         "font.size": 10, "legend.frameon": False})
    fig, ax = plt.subplots(2, 2, figsize=(11, 8.5), constrained_layout=True)

    # (a) p_det vs SNR1 in P/T bins
    se = np.logspace(-1, np.log10(5), 16)
    sc = np.sqrt(se[1:] * se[:-1])
    pbins = [(0, 0.5), (0.5, 0.8), (0.8, 1.2), (1.2, 3.0)]
    for i, (lo, hi) in enumerate(pbins):
        m = (x >= lo) & (x < hi)
        f, e, n = binned_fraction(snr[m], det[m], se)
        ok = n >= 5
        ax[0, 0].errorbar(sc[ok], f[ok], e[ok], marker="o", ms=4, lw=2, capsize=0,
                          color=blues(0.4 + 0.18 * i),
                          label=f"period {lo * T_DR4_DAYS / 365.25:.1f}-{hi * T_DR4_DAYS / 365.25:.1f} yr")
    ax[0, 0].axvline(1.0, color=muted, lw=1, ls=":")
    ax[0, 0].set(xscale="log", xlabel=r"SNR$_1$ = $a_0/\sigma_{\rm fov}$", ylabel="detection probability",
                 title="(a) detection vs signal", ylim=(-0.03, 1.03))
    ax[0, 0].legend(loc="upper left")

    # (b) p_det vs orbital period (years) in SNR1 bins
    from gedi import plotstyle as ps
    P_yr = x * T_DR4_DAYS / 365.25
    pe = np.geomspace(50 / 365.25, 3 * T_DR4_DAYS / 365.25, 16)
    pcen = np.sqrt(pe[1:] * pe[:-1])
    for i, (lo, hi) in enumerate([(0.5, 1), (1, 2), (2, 5)]):
        m = (snr >= lo) & (snr < hi)
        f, e, n = binned_fraction(P_yr[m], det[m], pe)
        ok = n >= 5
        ax[0, 1].errorbar(pcen[ok], f[ok], e[ok], marker="o", ms=4, lw=2, capsize=0,
                          color=blues(0.45 + 0.25 * i), label=f"{lo} ≤ SNR$_1$ < {hi}")
    ps.years_axis(ax[0, 1], "x", lo=pe[0], hi=pe[-1], label="planet's orbital period")
    ps.mark_dr4(ax[0, 1])
    ax[0, 1].axvline(4.0, color=muted, lw=1, ls=":")
    ax[0, 1].text(3.9, 0.03, "L&W 4 yr", color=muted, fontsize=8, ha="right", va="bottom")
    ax[0, 1].set(ylabel="detection probability", title="(b) detection vs orbital period", ylim=(-0.03, 1.03))
    ax[0, 1].legend(loc="upper left")

    # (c) Delta chi2 vs SNR1^2
    short = x < 0.5
    ax[1, 0].scatter(snr[short] ** 2, d["dchi2_true"][short], s=8, color=blues(0.75),
                     alpha=0.5, lw=0, label="period < 2.7 yr (true orbit)")
    xx = np.logspace(-2, np.log10(25), 50)
    ax[1, 0].plot(xx, 22 * xx, color=ink, lw=2, label=r"L&W DR4 (their Eq. 31): 22 SNR$_1^2$")
    ax[1, 0].plot(xx, 52 * xx, color=muted, lw=1.2, ls="--", label=r"L&W DR5 (Eq. 26): 52 SNR$_1^2$")
    ax[1, 0].axhline(50, color=muted, lw=1, ls=":")
    ax[1, 0].set(xscale="log", yscale="log", xlabel=r"SNR$_1^2$", ylabel=r"$\Delta\chi^2$",
                 title=r"(c) $\Delta\chi^2$ scaling", ylim=(1, None))
    ax[1, 0].legend(loc="upper left")

    # (d) mass recovery for detections
    rec = det & short & np.isfinite(d["m_fit"])
    ax[1, 1].scatter(snr[rec], d["m_fit"][rec] / d["m_true"][rec], s=8, color=blues(0.75),
                     alpha=0.5, lw=0)
    ax[1, 1].axhline(1, color=ink, lw=1)
    ax[1, 1].set(xscale="log", yscale="log", xlabel=r"SNR$_1$", ylabel=r"$m_{\rm fit}/m_{\rm true}$",
                 title="(d) mass recovery, detected, period < 2.7 yr", ylim=(0.3, 3))
    for axis in (ax[1, 1].xaxis, ax[1, 1].yaxis):
        axis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
        axis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax[1, 1].set_yticks([0.3, 0.5, 1, 2, 3])
    ax[1, 1].set_xticks([1, 2, 3, 5])

    fig.suptitle("Step 1 - single-planet validation against Lammers & Winn (2026)", color=ink)
    fig.savefig(PNG_PATH, dpi=150)
    print(f"saved {PNG_PATH}")


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("-n", type=int, default=10000, help="number of injections")
    p.add_argument("-j", "--jobs", type=int, default=-1, help="parallel workers (-1 = all)")
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--plot-only", action="store_true")
    a = p.parse_args()
    if not a.plot_only:
        run_all(a.n, a.jobs, a.seed)
    analyse()


if __name__ == "__main__":
    main()
