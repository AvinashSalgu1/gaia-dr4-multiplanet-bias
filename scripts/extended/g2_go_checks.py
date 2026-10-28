"""Comparison with Greenfield & Owen (2026, arXiv:2610.07149): near-2:1 pairs and evidence thresholds.

Needs gaiamock and dynesty. Three parts, each resumable (one CSV row per simulated system):

    python scripts/extended/g2_go_checks.py resoecc   [--dr dr4] [-j 64]   # minutes
    python scripts/extended/g2_go_checks.py single    [-j 64]              # ~1-2 h on 64 cores
    python scripts/extended/g2_go_checks.py second    [-j 64]              # ~3-5 h on 64 cores
    python scripts/extended/g2_go_checks.py <part> --analyse              # summary and figure only

Set-up as in Greenfield & Owen (GO) where they state it: a 1 M_sun star at 16 pc, proper motion
(5, -24) mas/yr, sky position (45 deg, 45 deg), Gaussian noise of 30 micro-arcsec per field-of-view
transit, circular face-on orbits (random inclinations in part 2 as a second case). GO take epochs from
GOST; here they come from gaiamock's GOST-based scanning law at the same position (10% of transits
dropped, as everywhere in this project). Planet-planet interactions are ignored (Keplerian sum).
S/N is GO's per-transit signal-to-noise: reflex semi-major axis a0 / 30 micro-arcsec.

resoecc  Circular, coplanar, face-on pairs at exactly P_out / P_in = 2 with mass ratios
         m_out / m_in = 0.2, 0.4, 0.8, 1, 1.5, 2.5, 5, outer periods 0.6-8 yr, the stronger signal at
         S/N 5 or 10. Each is fitted with one orbit (and with two orbits; BIC choice between them).
         Compared with the analytic resonance-eccentricity mapping e_eff = 2^(1/3) m_in / m_out
         (GO Eq. 22; Yahalomi et al. 2026) and with GO's quoted biases (e ~0.8 and mass +100% at
         m_out/m_in = 0.8; e ~0.4 and mass +10% at 2.5).
single   One planet, e = 0.01, P 0.3-3 yr, S/N log-uniform 0.3-3, face-on and random inclination.
         Detected if log10 of the Bayes factor (one orbit vs no planet) > 2 ("decisive", as GO).
         The 50% S/N threshold is compared with GO's DR4 values (0.89 face-on, 1.01 random), and with
         the Lammers & Winn rule used in this project (Delta chi2 > 50 plus mass and period checks).
second   Two circular face-on planets at period ratios 1.7, 1.9, 1.95, 2.0, 2.05, 2.1, 2.3 and mass
         ratios 0.2, 0.4, 0.8, 1, 1.5, 2.5, 5 (as GO), inner period 0.3-2 yr, S/N of the WEAKER planet
         log-uniform 0.2-5. The second planet is detected if log10 BF(two orbits vs one orbit) > 2.
         The 50% threshold of the weaker planet per configuration is compared with GO (0.5-2, higher
         near 2:1). The BIC choice and the project's two-orbit detection rule are recorded too.
Output: results/extended/go_checks/<part>/ (results.csv, summary.txt, figure).
Evidence: linear parameters integrated analytically with the Gaussian priors of g2_evidence.py
(Thiele-Innes constants 1 mas x --prior-scale); period (10-10,000 d, log-uniform), phase and
eccentricity sampled by dynesty. GO use log-uniform priors on the orbit size, which favour small
orbits much more than a 1-mas Gaussian does, so the evidence thresholds depend on the prior: parts
single and second are therefore run with --prior-scale 1 and 0.1 (folders <part> and <part>_prior0.1).
"""
from __future__ import annotations

import argparse
import csv

import numpy as np

import common as C

SIG = 0.030                          # mas per transit
M_STAR, PLX, PMRA, PMDEC = 1.0, 62.5, 5.0, -24.0
RA, DEC = 45.0, 45.0
LN10 = np.log(10.0)
Q_GRID = [0.2, 0.4, 0.8, 1.0, 1.5, 2.5, 5.0]               # m_out / m_in
PR_GRID = [1.7, 1.9, 1.95, 2.0, 2.05, 2.1, 2.3]            # P_out / P_in
GO_SINGLE = {"face-on": 0.89, "random": 1.01}              # GO DR4 single-planet S/N thresholds


# ---------------------------------------------------------------- simulation helpers
def epochs(dr, seed):
    from gaiabias import epochs as E
    return E.from_gaiamock(RA, DEC, dr, rng=np.random.default_rng([seed, 1]))


def mass_for_snr(snr, P_days, m_inner=0.0):
    """Planet mass [M_J] whose reflex a0 is snr * SIG (a0 is linear in m for m << M_star)."""
    from gaiabias.orbits import reflex_a0_mas
    return snr * SIG / reflex_a0_mas(M_STAR, 1.0, P_days, PLX, m_inner)


def observe(ep, comps, rng):
    from gaiabias.orbits import reflex_along_scan
    s, c = np.sin(ep.psi), np.cos(ep.psi)
    sig, a0 = reflex_along_scan(ep.t_yr, ep.psi, comps, M_STAR, PLX)
    y = PMRA * ep.t_yr * s + PMDEC * ep.t_yr * c + PLX * ep.plx_factor + sig + SIG * rng.standard_normal(len(ep))
    return y, np.full(len(ep), SIG), a0


def companion(m, P, rng, inc=0.0, Omega=0.0, ecc=0.0):
    from gaiabias.orbits import Companion
    return Companion(m_mjup=m, P_days=P, ecc=ecc, inc=inc, Omega=Omega, omega=rng.uniform(0, 2 * np.pi),
                     Tp_days=rng.uniform(0, P))


def lnz(ep, y, err, model, rng, nlive, scale=1.0):
    """(ln Z, its error); scale multiplies the orbit and trend prior widths (1 mas)."""
    import g2_evidence as GE
    return GE.evidence(ep, y, err, model, scale, nlive, rng)[:2]


def bic(chi2, k, n):
    return chi2 + k * np.log(n)


# ---------------------------------------------------------------- part 1: resonance-eccentricity
def tasks_resoecc(dr):
    P_out = np.geomspace(0.6, 8.0, 7) * 365.25
    out, s = [], 0
    for snr in (5.0, 10.0):
        for q in Q_GRID:
            for P in P_out:
                for rep in range(12):
                    out.append(dict(part="resoecc", dr=dr, seed=s, q=q, P_out=P, snr_dom=snr, rep=rep))
                    s += 1
    return out


def run_resoecc(t):
    from gaiabias.fitting import fit_one_companion, fit_two_planets
    rng = np.random.default_rng([t["seed"], 11])
    ep = epochs(t["dr"], t["seed"])
    P_out, q = t["P_out"], t["q"]
    P_in = P_out / 2
    ratio = q * 2 ** (2 / 3)                       # S/N_out / S/N_in
    snr_out, snr_in = (t["snr_dom"], t["snr_dom"] / ratio) if ratio >= 1 else (t["snr_dom"] * ratio, t["snr_dom"])
    m_in = mass_for_snr(snr_in, P_in)
    m_out = mass_for_snr(snr_out, P_out, m_in * 9.546e-4)
    comps = [companion(m_in, P_in, rng), companion(m_out, P_out, rng)]
    y, err, a0 = observe(ep, comps, rng)
    one = fit_one_companion(ep, y, err, M_STAR, uncertainties=False)
    two = fit_two_planets(ep, y, err, M_STAR)
    n = len(ep)
    return dict(t, n_obs=n, m_in=m_in, m_out=m_out, snr_in=a0[0] / SIG, snr_out=a0[1] / SIG,
                P1=one.P_days, e1=one.ecc, m1=one.m_mjup, chi2_1=one.chi2,
                P2a=two.P[0], P2b=two.P[1], m2a=two.m_mjup[0], m2b=two.m_mjup[1], chi2_2=two.chi2,
                two_by_bic=bool(bic(two.chi2, 19, n) < bic(one.chi2, 12, n)))


def analyse_resoecc(rows, out):
    plt, ps = C.setup_plots()
    L = ["=== RESONANCE-ECCENTRICITY CHECK (circular, coplanar, face-on 2:1 pairs; one-orbit fits) ===",
         "Columns per mass ratio q = m_out/m_in: the one-orbit fit locks on P_out (within 10%) or P_in; median",
         "fitted e [16-84%] when locked on P_out; analytic e_eff = 1.26/q (valid for small e); median fitted mass",
         "/ m_out - 1; share of systems where the BIC prefers two orbits. GO quote e ~0.8 and mass +100% at",
         "q = 0.8, e ~0.4 and mass +10% at q = 2.5.", ""]
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    for dr in sorted({r["dr"] for r in rows}):
        for snr in sorted({float(r["snr_dom"]) for r in rows}):
            sel = [r for r in rows if r["dr"] == dr and float(r["snr_dom"]) == snr]
            L.append(f"--- {dr.upper()}, stronger planet at S/N {snr:g} ({len(sel)} systems) ---")
            L.append("    q     on P_out  on P_in   e (on P_out)        e_eff   mass bias   BIC two orbits")
            qs, med, lo, hi, mb = [], [], [], [], []
            for q in Q_GRID:
                s = [r for r in sel if float(r["q"]) == q]
                if not s:
                    continue
                P1, Po = np.array([float(r["P1"]) for r in s]), np.array([float(r["P_out"]) for r in s])
                on_out = np.abs(P1 / Po - 1) < 0.1
                on_in = np.abs(P1 / (Po / 2) - 1) < 0.1
                e = np.array([float(r["e1"]) for r in s])[on_out]
                mbias = (np.array([float(r["m1"]) for r in s]) / np.array([float(r["m_out"]) for r in s]) - 1)[on_out]
                two = np.mean([r["two_by_bic"] in (True, "True") for r in s])
                eq = np.percentile(e, [50, 16, 84]) if len(e) else [np.nan] * 3
                L.append(f"  {q:4.1f}   {100 * on_out.mean():5.0f}%   {100 * on_in.mean():5.0f}%   "
                         f"{eq[0]:.2f} [{eq[1]:.2f}-{eq[2]:.2f}]   {min(1.26 / q, 9.99):5.2f}   "
                         f"{100 * np.median(mbias) if len(mbias) else np.nan:+6.0f}%      {100 * two:4.0f}%")
                qs.append(q), med.append(eq[0]), lo.append(eq[1]), hi.append(eq[2])
                mb.append(np.median(mbias) if len(mbias) else np.nan)
            lab = f"{dr.upper()}, S/N {snr:g}"
            axs[0].errorbar(qs, med, yerr=[np.subtract(med, lo), np.subtract(hi, med)], marker="o", capsize=3,
                            label=lab)
            axs[1].plot(qs, 100 * np.array(mb), marker="o", label=lab)
            L.append("")
    qq = np.geomspace(0.2, 5, 100)
    axs[0].plot(qq, np.minimum(1.26 / qq, 1), "k--", lw=1, label="e_eff = 1.26 m_in/m_out")
    axs[0].scatter([0.8, 2.5], [0.8, 0.4], marker="*", s=140, color="k", zorder=5, label="GO quoted")
    axs[1].scatter([0.8, 2.5], [100, 10], marker="*", s=140, color="k", zorder=5, label="GO quoted")
    from matplotlib.ticker import NullFormatter
    for ax in axs:
        ax.set_xscale("log")
        ax.set_xticks(Q_GRID, [f"{q:g}" for q in Q_GRID])
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_xlabel("mass ratio m_out / m_in")
        ax.legend(frameon=False, fontsize=8)
    axs[0].set_ylabel("eccentricity of the one-orbit fit")
    axs[1].set_ylabel("one-orbit mass / outer mass - 1 (%)")
    axs[0].set_ylim(0, 1)
    fig.savefig(out / "resoecc.png", dpi=150)
    return L


# ---------------------------------------------------------------- part 2: single-planet threshold
def tasks_single(n=260):
    out = []
    for k, geo in enumerate(("face-on", "random")):
        for i in range(n):
            out.append(dict(part="single", geo=geo, seed=100000 * (k + 1) + i))
    return out


def run_single(t):
    from gaiabias.fitting import fit_one_companion, lammers_winn_detected
    rng = np.random.default_rng([t["seed"], 22])
    ep = epochs("dr4", t["seed"])
    P = np.exp(rng.uniform(np.log(0.3), np.log(3.0))) * 365.25
    snr = np.exp(rng.uniform(np.log(0.3), np.log(3.0)))
    inc = 0.0 if t["geo"] == "face-on" else float(np.arccos(rng.uniform(-1, 1)))
    c = companion(mass_for_snr(snr, P), P, rng, inc=inc, Omega=rng.uniform(0, np.pi), ecc=0.01)
    y, err, a0 = observe(ep, [c], rng)
    sc = float(t.get("prior_scale", 1.0))
    z0, _ = lnz(ep, y, err, "no planet", rng, 400, sc)
    z1, z1e = lnz(ep, y, err, "one orbit", rng, 400, sc)
    fit = fit_one_companion(ep, y, err, M_STAR)
    return dict(t, P=P, snr=a0[0] / SIG, inc=inc, n_obs=len(ep), lnZ0=z0, lnZ1=z1, lnZ1err=z1e,
                log10BF=(z1 - z0) / LN10, dchi2=fit.delta_chi2, lw=bool(lammers_winn_detected(fit)))


def threshold(snr, det):
    """S/N at which the detection probability is 50% (logistic in ln S/N), with a bootstrap range."""
    from scipy.optimize import minimize
    x, y = np.log(np.asarray(snr, float)), np.asarray(det, float)

    def nll(p, xx, yy):
        pr = 1 / (1 + np.exp(np.clip(-(xx - p[0]) / max(abs(p[1]), 1e-3), -50, 50)))
        pr = np.clip(pr, 1e-9, 1 - 1e-9)
        return -np.sum(yy * np.log(pr) + (1 - yy) * np.log(1 - pr))

    def fit(xx, yy):
        if yy.all() or not yy.any():
            return np.nan
        return float(np.exp(minimize(nll, [np.median(xx), 0.2], args=(xx, yy), method="Nelder-Mead").x[0]))

    best = fit(x, y)
    rng = np.random.default_rng(3)
    bs = [fit(x[i], y[i]) for i in (rng.integers(0, len(x), len(x)) for _ in range(200))]
    lo, hi = np.nanpercentile(bs, [16, 84]) if np.isfinite(bs).any() else (np.nan, np.nan)
    return best, lo, hi


def analyse_single(rows, out):
    plt, ps = C.setup_plots()
    L = ["=== SINGLE-PLANET DETECTION THRESHOLD (DR4, 30 micro-arcsec per transit) ===",
         "S/N = a0 / 30 micro-arcsec. 50% thresholds [16-84% bootstrap]:", ""]
    fig, ax = plt.subplots(figsize=(6.4, 4.4), constrained_layout=True)
    bins = np.geomspace(0.3, 3, 9)
    for geo in ("face-on", "random"):
        s = [r for r in rows if r["geo"] == geo]
        if not s:
            continue
        snr = np.array([float(r["snr"]) for r in s])
        bf = np.array([float(r["log10BF"]) for r in s]) > 2
        lw = np.array([r["lw"] in (True, "True") for r in s])
        t_bf, t_lw = threshold(snr, bf), threshold(snr, lw)
        L.append(f"  {geo:8s} ({len(s)} planets): evidence (log10 BF > 2) {t_bf[0]:.2f} [{t_bf[1]:.2f}-{t_bf[2]:.2f}]"
                 f"   GO: {GO_SINGLE[geo]:.2f}   |   Lammers & Winn rule {t_lw[0]:.2f} [{t_lw[1]:.2f}-{t_lw[2]:.2f}]")
        for det, ls, lab in ((bf, "-", "evidence"), (lw, "--", "L&W rule")):
            idx = np.digitize(snr, bins) - 1
            frac = [det[idx == k].mean() if np.any(idx == k) else np.nan for k in range(len(bins) - 1)]
            ax.plot(np.sqrt(bins[:-1] * bins[1:]), frac, ls, marker="o", label=f"{geo}, {lab}")
        ax.axvline(GO_SINGLE[geo], color="0.5", lw=0.8, ls=":")
    L += ["", "The evidence threshold depends on the prior on the orbit size: rerun with --prior-scale to see it.",
          "The L&W rule (Delta chi2 > 50 with mass and period checks) is shown for the same planets; step 1's",
          "SNR_1 ~ 1.8 uses the gaiamock noise table (with bright-star jitter) and periods up to 3 x the time span."]
    from matplotlib.ticker import NullFormatter
    ax.set_xscale("log")
    ax.set_xticks([0.3, 0.5, 1, 2, 3], ["0.3", "0.5", "1", "2", "3"])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel("S/N per transit (a0 / sigma)")
    ax.set_ylabel("detected fraction")
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(out / "single.png", dpi=150)
    return L


# ---------------------------------------------------------------- part 3: second-planet threshold
def tasks_second(n=30):
    out, s = [], 0
    for pr in PR_GRID:
        for q in Q_GRID:
            for i in range(n):
                out.append(dict(part="second", seed=500000 + s, pr=pr, q=q))
                s += 1
    return out


def run_second(t):
    from gaiabias.fitting import fit_orbit_trend, fit_two_planets
    rng = np.random.default_rng([t["seed"], 33])
    ep = epochs("dr4", t["seed"])
    P_in = np.exp(rng.uniform(np.log(0.3), np.log(2.0))) * 365.25
    P_out = P_in * t["pr"]
    ratio = t["q"] * t["pr"] ** (2 / 3)            # S/N_out / S/N_in
    weak = np.exp(rng.uniform(np.log(0.2), np.log(5.0)))
    snr_in, snr_out = (weak, weak * ratio) if ratio >= 1 else (weak / ratio, weak)
    m_in = mass_for_snr(snr_in, P_in)
    m_out = mass_for_snr(snr_out, P_out, m_in * 9.546e-4)
    comps = [companion(m_in, P_in, rng), companion(m_out, P_out, rng)]
    y, err, a0 = observe(ep, comps, rng)
    sc = float(t.get("prior_scale", 1.0))
    z1, _ = lnz(ep, y, err, "one orbit", rng, 400, sc)
    z2, z2e = lnz(ep, y, err, "two orbits", rng, 400, sc)
    one = fit_orbit_trend(ep, y, err, M_STAR, n_trend=0)
    two = fit_two_planets(ep, y, err, M_STAR, one=one)
    n = len(ep)
    found = sorted(two.P)
    both_right = (abs(found[0] / P_in - 1) < 0.1) and (abs(found[1] / P_out - 1) < 0.1)
    return dict(t, P_in=P_in, P_out=P_out, m_in=m_in, m_out=m_out, snr_in=a0[0] / SIG, snr_out=a0[1] / SIG,
                snr_weak=min(a0) / SIG, n_obs=n, lnZ1=z1, lnZ2=z2, lnZ2err=z2e, log10BF=(z2 - z1) / LN10,
                bic_two=bool(bic(two.chi2, 19, n) < bic(one.chi2, 12, n)),
                rule_two=bool(min(two.drop_dchi2) > 50), both_periods_right=bool(both_right))


def analyse_second(rows, out):
    plt, ps = C.setup_plots()
    L = ["=== SECOND-PLANET DETECTION THRESHOLD (DR4; circular, face-on) ===",
         "50% threshold on the S/N of the WEAKER planet, per period ratio (rows) and mass ratio m_out/m_in (columns).",
         "Detected = log10 BF(two orbits vs one) > 2. GO find 0.5-2, higher near 2:1. '-' = no threshold in range.", ""]
    L.append("   P ratio " + "".join(f"{q:>8.1f}" for q in Q_GRID) + "    all")
    grid = np.full((len(PR_GRID), len(Q_GRID)), np.nan)
    for i, pr in enumerate(PR_GRID):
        line = f"   {pr:6.2f}  "
        for j, q in enumerate(Q_GRID):
            s = [r for r in rows if float(r["pr"]) == pr and float(r["q"]) == q]
            if len(s) < 8:
                line += "       ."
                continue
            th = threshold([float(r["snr_weak"]) for r in s], [float(r["log10BF"]) > 2 for r in s])[0]
            grid[i, j] = th
            line += f"{th:8.2f}" if np.isfinite(th) else "       -"
        s = [r for r in rows if float(r["pr"]) == pr]
        th = threshold([float(r["snr_weak"]) for r in s], [float(r["log10BF"]) > 2 for r in s]) if s else (np.nan,)
        L.append(line + (f"  {th[0]:6.2f}" if s else ""))
    near = [r for r in rows if 1.95 <= float(r["pr"]) <= 2.05]
    far = [r for r in rows if float(r["pr"]) in (1.7, 2.3)]
    if near and far:
        tn = threshold([float(r["snr_weak"]) for r in near], [float(r["log10BF"]) > 2 for r in near])
        tf = threshold([float(r["snr_weak"]) for r in far], [float(r["log10BF"]) > 2 for r in far])
        L += ["", f"  near 2:1 (1.95-2.05): {tn[0]:.2f} [{tn[1]:.2f}-{tn[2]:.2f}]   "
                  f"far (1.7 and 2.3): {tf[0]:.2f} [{tf[1]:.2f}-{tf[2]:.2f}]"]
        for lab, key in (("BIC prefers two orbits", "bic_two"), ("project two-orbit rule (both drop dchi2 > 50)",
                                                                 "rule_two")):
            a = threshold([float(r["snr_weak"]) for r in near], [r[key] in (True, "True") for r in near])
            b = threshold([float(r["snr_weak"]) for r in far], [r[key] in (True, "True") for r in far])
            L.append(f"  {lab}: near {a[0]:.2f}, far {b[0]:.2f}")
        det = [r for r in rows if float(r["log10BF"]) > 2]
        if det:
            L.append(f"  of systems where evidence prefers two orbits, both periods found within 10% by the two-orbit "
                     f"fit: {100 * np.mean([r['both_periods_right'] in (True, 'True') for r in det]):.0f}%")
    fig, ax = plt.subplots(figsize=(6.4, 4.6), constrained_layout=True)
    im = ax.imshow(grid, origin="lower", cmap="viridis", aspect="auto", vmin=0, vmax=3)
    ax.set_xticks(range(len(Q_GRID)), [f"{q:g}" for q in Q_GRID])
    ax.set_yticks(range(len(PR_GRID)), [f"{p:g}" for p in PR_GRID])
    for i in range(len(PR_GRID)):
        for j in range(len(Q_GRID)):
            if np.isfinite(grid[i, j]):
                ax.text(j, i, f"{grid[i, j]:.1f}", ha="center", va="center", fontsize=7, color="w")
    ax.set_xlabel("mass ratio m_out / m_in")
    ax.set_ylabel("period ratio P_out / P_in")
    fig.colorbar(im, ax=ax, label="50% S/N threshold of the weaker planet")
    fig.savefig(out / "second.png", dpi=150)
    return L


PARTS = dict(resoecc=(tasks_resoecc, run_resoecc, analyse_resoecc),
             single=(tasks_single, run_single, analyse_single),
             second=(tasks_second, run_second, analyse_second))


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("part", choices=list(PARTS))
    p.add_argument("--dr", default="dr4", help="resoecc only: data release of the epochs (dr4 or dr5)")
    p.add_argument("-j", "--jobs", type=int, default=-1)
    p.add_argument("--analyse", action="store_true")
    p.add_argument("--quick", action="store_true", help="a handful of systems (code test)")
    p.add_argument("--prior-scale", type=float, default=1.0, help="single/second: prior widths x this")
    p.add_argument("-n", type=int, default=None, help="systems per configuration (single: per geometry)")
    a = p.parse_args()
    from joblib import Parallel, delayed
    make, run, analyse = PARTS[a.part]
    tag = "" if a.prior_scale == 1 or a.part == "resoecc" else f"_prior{a.prior_scale:g}"
    out = C.out_dir(f"go_checks/{a.part}{tag}" + ("_quick" if a.quick else ""))
    path = out / "results.csv"
    tasks = make(a.dr) if a.part == "resoecc" else (make(a.n) if a.n else make())
    for t in tasks:
        t["prior_scale"] = a.prior_scale
    if a.quick:
        rng = np.random.default_rng(0)
        tasks = [tasks[i] for i in rng.choice(len(tasks), 6, replace=False)]
    if not a.analyse:
        done = set()
        if path.exists():
            done = {(r.get("dr", ""), r["seed"]) for r in csv.DictReader(open(path))}
        todo = [t for t in tasks if (t.get("dr", ""), str(t["seed"])) not in done]
        print(f"{len(tasks) - len(todo)} systems done, {len(todo)} to go", flush=True)
        if todo:
            t0 = C.Timer()
            gen = Parallel(n_jobs=a.jobs, return_as="generator_unordered")(delayed(run)(t) for t in todo)
            w, fh = None, None
            for k, row in enumerate(gen, 1):
                if w is None:
                    new = not path.exists()
                    fh = open(path, "a", newline="")
                    w = csv.DictWriter(fh, fieldnames=list(row), extrasaction="ignore")
                    if new:
                        w.writeheader()
                w.writerow(row)
                fh.flush()
                if k % 25 == 0 or k == len(todo):
                    print(f"  {k}/{len(todo)} ({t0})", flush=True)
            fh.close()
    rows = list(csv.DictReader(open(path)))
    L = analyse(rows, out)
    L += ["", "Reference: Greenfield, E. J. & Owen, J. E. 2026, arXiv:2610.07149 (submitted to MNRAS)."]
    C.write_summary(out / "summary.txt", L)


if __name__ == "__main__":
    main()
