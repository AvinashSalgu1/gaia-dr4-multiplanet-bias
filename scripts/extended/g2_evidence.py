"""Does Bayesian evidence pick the same DR4 model as the BIC?

Needs gaiamock (scanning law) and `dynesty` (pip install dynesty). Resumable; about 30-60 min
on a 12-core machine with the defaults:

    python scripts/extended/g2_evidence.py           # [--per-group 50] [-j -1]
    python scripts/extended/g2_evidence.py --analyse

Systems from the step 2 experiment (same seeds as step 2) in four groups: weak target
(per-transit signal-to-noise below 2.5), eccentric (either planet e > 0.5), near 2:1 (period
ratio 1.8-2.25) and random. For each system, with both planets present:

* the BIC choice among the seven models of step 2 (no planet, acceleration, jerk, one orbit,
  orbit + acceleration, orbit + jerk, two orbits), from the step 2 fits;
* the Bayesian evidence of the same seven models. The position, proper motion, parallax,
  Thiele-Innes and trend terms are linear and are integrated analytically with Gaussian
  priors (widths below); the period, phase and eccentricity of each orbit are sampled with
  nested sampling (dynesty): log-uniform period 10-10,000 d, uniform phase, Kipping (2013)
  eccentricity prior.

The evidence choice is the model with the highest evidence, and the summary compares it
with the BIC choice, overall and per group. The nested-sampling error of every ln Z is kept;
systems where it exceeds 1 are counted in the summary (rerun them with more --nlive).

Linear prior widths (1 sigma): positions 100 mas, proper motions 100 mas/yr, parallax 500 mas,
Thiele-Innes constants 1 mas (x --prior-scale), acceleration 1 mas/yr^2, jerk 1 mas/yr^3.
The evidence depends on these widths; --prior-scale 3 repeats the analysis with 3x wider orbit
and trend priors (written to a separate folder).
"""
from __future__ import annotations

import argparse
import csv

import numpy as np

import common as C

MODELS = ["no planet", "acceleration only", "jerk only", "one orbit", "orbit + acceleration", "orbit + jerk",
          "two orbits"]
NL = {"no planet": 0, "acceleration only": 0, "jerk only": 0, "one orbit": 1, "orbit + acceleration": 1,
      "orbit + jerk": 1, "two orbits": 2}
TREND = {"no planet": 0, "acceleration only": 2, "jerk only": 4, "one orbit": 0, "orbit + acceleration": 2,
         "orbit + jerk": 4, "two orbits": 0}
BIC_K = [5, 7, 9, 12, 14, 16, 19]
LNP = (np.log(10.0), np.log(1e4))


def log_evidence_linear(X, y, err, s):
    """ln of the marginal likelihood with a Gaussian prior N(0, diag(s^2)) on the linear parameters."""
    w = 1.0 / err ** 2
    A = X.T @ (X * w[:, None]) + np.diag(1.0 / s ** 2)
    b = X.T @ (y * w)
    try:
        L = np.linalg.cholesky(A)
    except np.linalg.LinAlgError:
        return -np.inf
    z = np.linalg.solve(L, b)
    return float(-0.5 * np.sum(y ** 2 * w) + 0.5 * z @ z - np.sum(np.log(np.diag(L))) - np.sum(np.log(s))
                 - 0.5 * np.sum(np.log(2 * np.pi / w)))


def design(ep, model, nl):
    from gaiabias.fitting import _design_5par, _trend_columns
    from gaiabias.orbits import orbit_XY
    cols = [_design_5par(ep)] + _trend_columns(ep, TREND[model])
    s, c = np.sin(ep.psi), np.cos(ep.psi)
    for k in range(NL[model]):
        P, ph, e = nl[3 * k:3 * k + 3]
        X, Y = orbit_XY(ep.t_yr, P, e, ph)
        cols.append(np.column_stack([X * s, Y * s, X * c, Y * c]))
    return np.column_stack(cols)


def prior_widths(model, scale):
    s = [100.0, 100.0, 100.0, 100.0, 500.0]
    if TREND[model] >= 2:
        s += [1.0 * scale] * 2
    if TREND[model] >= 4:
        s += [1.0 * scale] * 2
    s += [1.0 * scale] * (4 * NL[model])
    return np.array(s)


def evidence(ep, y, err, model, scale, nlive, rng):
    """(ln Z, its error, highest ln marginal likelihood seen by the sampler) for one model."""
    s = prior_widths(model, scale)
    k = NL[model]
    if k == 0:
        return log_evidence_linear(design(ep, model, []), y, err, s), 0.0, np.nan
    import dynesty
    from scipy.stats import beta

    def ptf(u):
        x = np.empty_like(u)
        for j in range(k):
            x[3 * j] = np.exp(LNP[0] + u[3 * j] * (LNP[1] - LNP[0]))
            x[3 * j + 1] = 2 * np.pi * u[3 * j + 1]
            x[3 * j + 2] = min(beta.ppf(u[3 * j + 2], 0.867, 3.03), 0.95)
        return x

    def lnl(x):
        return log_evidence_linear(design(ep, model, x), y, err, s)

    ns = dynesty.NestedSampler(lnl, ptf, 3 * k, nlive=nlive * k, sample="rwalk", rstate=rng)
    ns.run_nested(dlogz=0.5, print_progress=False)
    r = ns.results
    return float(r.logz[-1]), float(r.logzerr[-1]), float(np.max(r.logl))


def select_seeds(per_group):
    """Seeds from the step 2 population for the four groups (system draws only, no epochs needed)."""
    from gaiabias.step2.sim import draw_system, population_seeds
    from gaiabias.simulate import snr1
    groups = {"weak": [], "eccentric": [], "near 2:1": [], "random": []}
    for s in population_seeds(20000):
        sy = draw_system(np.random.default_rng(s))
        t, c = sy["target"], sy["comp"]
        from gaiabias.orbits import reflex_a0_mas
        snr_t = snr1(reflex_a0_mas(sy["m_star"], t.m_mjup, t.P_days, sy["parallax"]), sy["G"])
        pr = max(t.P_days, c.P_days) / min(t.P_days, c.P_days)
        for g, ok in (("weak", snr_t < 2.5), ("eccentric", max(t.ecc, c.ecc) > 0.5),
                      ("near 2:1", 1.8 <= pr <= 2.25), ("random", True)):
            if ok and len(groups[g]) < per_group:
                groups[g].append(s)
                break
        if all(len(v) >= per_group for v in groups.values()):
            break
    return [(g, s) for g, v in groups.items() for s in v]


def run_one(group, seed, scale, nlive, synthetic):
    import sys
    sys.path.insert(0, str(C.REPO / "scripts" / "extended"))
    from g2_noise_robustness import make_data
    from gaiabias.step2 import analyses as A
    s, ep, single, both, noise = make_data(seed, "baseline", synthetic)
    y, err = both.obs + noise, both.err
    fit = A.run(ep, y, err, s["m_star"], "w", set(A.FITTED))
    chi2 = [fit.get(k, np.nan) for k in ("accel_w_c5", "accel_w_c7", "jerk_w_c9", "accel_w_c0", "accel_w_c2",
                                          "jerk_w_cJ", "two_w_chi2")]
    bic = np.array(chi2) + np.array(BIC_K) * np.log(len(ep))
    bic_choice = MODELS[int(np.nanargmin(bic))] if np.isfinite(bic).any() else "none"
    rng = np.random.default_rng([seed, 5])
    row = dict(group=group, seed=seed, n_transits=len(ep), bic_choice=bic_choice)
    lnz = {}
    for i, m in enumerate(MODELS):
        z, ze, lmax = evidence(ep, y, err, m, scale, nlive, rng)
        lnz[m] = z
        row[f"lnZ_{m}"], row[f"lnZerr_{m}"] = z, ze
        row[f"chi2opt_{m}"] = chi2[i]
        row[f"lmax_{m}"] = lmax
    row["evidence_choice"] = max(lnz, key=lnz.get)
    return row


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--per-group", type=int, default=50)
    p.add_argument("--nlive", type=int, default=400, help="live points per orbit")
    p.add_argument("--prior-scale", type=float, default=1.0)
    p.add_argument("-j", "--jobs", type=int, default=-1)
    p.add_argument("--analyse", action="store_true")
    p.add_argument("--quick", action="store_true", help="2 systems per group, 100 live points (test)")
    p.add_argument("--synthetic", action="store_true", help="synthetic epochs (code test only)")
    a = p.parse_args()
    from joblib import Parallel, delayed
    name = "evidence" if a.prior_scale == 1 else f"evidence_prior{a.prior_scale:g}"
    out = C.out_dir(name + ("_quick" if a.quick else ""))
    if a.quick:
        a.per_group, a.nlive = 2, 100
    path = out / "evidence.csv"
    cols = ["group", "seed", "n_transits", "bic_choice", "evidence_choice"] + \
        [f"{q}_{m}" for m in MODELS for q in ("lnZ", "lnZerr", "chi2opt", "lmax")]
    if not a.analyse:
        done = set()
        if path.exists():
            done = {int(r["seed"]) for r in csv.DictReader(open(path))}
        else:
            with open(path, "w", newline="") as f:
                csv.writer(f).writerow(cols)
        todo = [(g, s) for g, s in select_seeds(a.per_group) if s not in done]
        print(f"{len(done)} systems done, {len(todo)} to go", flush=True)
        t = C.Timer()
        tasks = (delayed(run_one)(g, s, a.prior_scale, a.nlive, a.synthetic) for g, s in todo)
        # every worker stays busy: rows are written as each system finishes (no waiting for the
        # slowest system of a batch), so an interrupted run loses at most the systems in progress
        with open(path, "a", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
            for k, r in enumerate(Parallel(n_jobs=a.jobs, return_as="generator_unordered")(tasks), 1):
                w.writerow(r)
                f.flush()
                if k % 24 == 0 or k == len(todo):
                    print(f"  {k}/{len(todo)} ({t})", flush=True)
    rows = list(csv.DictReader(open(path)))
    L = ["=== BAYESIAN EVIDENCE VS THE BIC CHOICE ===",
         f"{len(rows)} step-2 systems (both planets present); linear prior scale x{a.prior_scale:g}.",
         "Rows: BIC choice; columns: evidence choice (counts).", ""]
    short = ["none", "acc", "jerk", "orbit", "orb+acc", "orb+jerk", "two"]
    for g in ["all", "weak", "eccentric", "near 2:1", "random"]:
        sel = rows if g == "all" else [r for r in rows if r["group"] == g]
        if not sel:
            continue
        agree = np.mean([r["bic_choice"] == r["evidence_choice"] for r in sel])
        n_planets = lambda m: NL.get(m, 0)  # noqa: E731
        same_n = np.mean([n_planets(r["bic_choice"]) == n_planets(r["evidence_choice"]) for r in sel])
        L.append(f"--- {g} (n={len(sel)}): same model {100 * agree:.0f}%; same number of orbits {100 * same_n:.0f}% ---")
        L.append("              " + "".join(f"{s:>9}" for s in short))
        for i, mb in enumerate(MODELS):
            cnt = [sum(r["bic_choice"] == mb and r["evidence_choice"] == me for r in sel) for me in MODELS]
            if sum(cnt):
                L.append(f"   {short[i]:<10} " + "".join(f"{c:>9}" for c in cnt))
        L.append("")
    L.append("Sampler check: the nested-sampling error of every ln Z is in evidence.csv (lnZerr_*).")
    big = sum(any(float(r[f"lnZerr_{m}"]) > 1 for m in MODELS) for r in rows)
    L.append(f"systems with any lnZ error > 1: {big}")
    C.write_summary(out / "summary.txt", L)
    plt, ps = C.setup_plots()
    fig, ax = plt.subplots(figsize=(6.4, 5.4), constrained_layout=True)
    M = np.array([[sum(r["bic_choice"] == mb and r["evidence_choice"] == me for r in rows) for me in MODELS]
                  for mb in MODELS], float)
    Mn = M / np.maximum(M.sum(1, keepdims=True), 1)
    ax.imshow(Mn, cmap="Blues", vmin=0, vmax=1)
    for i in range(len(MODELS)):
        for j in range(len(MODELS)):
            if M[i, j]:
                ax.text(j, i, f"{int(M[i, j])}", ha="center", va="center", fontsize=8,
                        color="white" if Mn[i, j] > 0.5 else ps.INK)
    ax.set_xticks(range(len(MODELS)), short, rotation=45)
    ax.set_yticks(range(len(MODELS)), short)
    ax.set_xlabel("model with the highest Bayesian evidence")
    ax.set_ylabel("model with the lowest BIC")
    ax.set_title(f"BIC vs evidence ({len(rows)} systems)", fontsize=10.5)
    fig.savefig(out / "evidence.png", dpi=150, bbox_inches="tight")
    print(f"saved {out / 'evidence.png'}")


if __name__ == "__main__":
    main()
