"""Do the step 2 results survive a less ideal noise model?

Needs gaiamock (scanning law). Resumable; about 1.5-3 hours on a 12-core machine with the defaults:

    python scripts/extended/g2_noise_robustness.py              # [-n 1000] [-j -1]
    python scripts/extended/g2_noise_robustness.py --analyse

The first n systems of the step 2 experiment (same seeds, so the same planets, epochs and
baseline noise as step 2) are refitted with the DR4 models under five noise models. The
fits always use the REPORTED uncertainties of the reference model, as a DR4 user would:

* baseline     - as in step 2 (per-CCD errors / sqrt(8), extra unreported jitter for G < 13);
* ccd_corr     - CCD measurements in a transit correlated with rho = 0.05, so the true
                 per-transit error is sqrt(1 + 7 rho) = 1.16 times the reported one
                 (the DR4 pre-release gives rho = -0.001 +/- 0.007);
* transit_corr - an extra error shared by all transits within 4 days (one visibility
                 period), with 0.5 times the per-transit error;
* excess       - an extra unreported error of 0.03 mas, growing by 50% per magnitude fainter
                 than G = 13;
* underest     - all errors 30% larger than reported.

Output in results/extended/noise_robustness/: one CSV per model, summary.txt (the step 2
headline numbers under each model) and noise_robustness.png.
"""
from __future__ import annotations

import argparse
import csv

import numpy as np

import common as C

MODELS = ["baseline", "ccd_corr", "transit_corr", "excess", "underest"]
KEEP = ["seed", "n_transits", "G", "P_t", "m_t", "P_c", "m_c", "rho"]


def make_data(seed, model, synthetic=False):
    """Step 2 system, epochs and noise (identical to step2.sim.make_data for 'baseline')."""
    from gedi import epochs as E
    from gedi.simulate import Star, sigma_al_per_transit, simulate
    from gedi.step2.sim import draw_system

    rng = np.random.default_rng(seed)
    s = draw_system(rng)
    ep = E.synthetic(80, rng=rng) if synthetic else E.from_gaiamock(s["ra"], s["dec"], "dr4", rng=rng)
    star = Star(parallax_mas=s["parallax"], m_star_msun=s["m_star"], G_mag=s["G"], pmra=s["pmra"], pmdec=s["pmdec"])
    t, c = s["target"], s["comp"]
    single = simulate(ep, star, [t], add_noise=False)
    both = simulate(ep, star, [t, c], add_noise=False)
    sigma = float(sigma_al_per_transit(s["G"]))
    gauss = sigma * rng.standard_normal(len(ep))
    jit = rng.uniform(0, 0.04) * rng.standard_normal(len(ep)) if s["G"] < 13 else np.zeros(len(ep))
    r2 = np.random.default_rng([seed, 77])
    if model == "ccd_corr":
        gauss = gauss * np.sqrt(1 + 7 * 0.05)
    elif model == "underest":
        gauss = gauss * 1.3
        jit = jit * 1.3
    noise = gauss + jit
    if model == "transit_corr":
        order = np.argsort(ep.t_yr)
        grp = np.zeros(len(ep), int)
        grp[order[1:]] = np.cumsum(np.diff(ep.t_yr[order]) * 365.25 > 4.0)
        noise = noise + (0.5 * sigma * r2.standard_normal(grp.max() + 1))[grp]
    elif model == "excess":
        ex = 0.03 * 1.5 ** max(0.0, s["G"] - 13.0)
        noise = noise + ex * r2.standard_normal(len(ep))
    return s, ep, single, both, noise


def run_one(seed, model, synthetic):
    from gedi.step2 import analyses as A
    s, ep, single, both, noise = make_data(seed, model, synthetic)
    row = dict(seed=seed, n_transits=len(ep), G=s["G"], P_t=s["target"].P_days, m_t=s["target"].m_mjup,
               P_c=s["comp"].P_days, m_c=s["comp"].m_mjup, rho=np.nan)
    for v, obs in (("w", both.obs + noise), ("a", single.obs + noise)):
        row.update(A.run(ep, obs, both.err, s["m_star"], v, set(A.FITTED)))
    return row


def run(model, n, jobs, synthetic, out):
    from joblib import Parallel, delayed
    from gedi.step2 import analyses as A
    from gedi.step2.sim import population_seeds
    cols = KEEP + A.all_cols()
    path = out / f"{model}.csv"
    done = set()
    if path.exists():
        done = {int(float(r["seed"])) for r in csv.DictReader(open(path))}
    else:
        with open(path, "w", newline="") as f:
            csv.writer(f).writerow(cols)
    todo = [s for s in population_seeds(n) if s not in done]
    print(f"{model}: {len(done)} systems done, {len(todo)} to go", flush=True)
    t = C.Timer()
    for i in range(0, len(todo), 120):
        rows = Parallel(n_jobs=jobs)(delayed(run_one)(s, model, synthetic) for s in todo[i:i + 120])
        with open(path, "a", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
            for r in rows:
                w.writerow({k: r.get(k, "") for k in cols})
        print(f"  {model}: {i + len(rows)}/{len(todo)} ({t})", flush=True)


def load(path):
    from gedi.step2 import analyses as A
    rows = list(csv.DictReader(open(path)))
    d = {}
    for k in rows[0]:
        d[k] = np.array([r[k] if k in A.TEXT_FIELDS else (float(r[k]) if r[k] not in ("", "nan") else np.nan)
                         for r in rows])
    return d


def analyse(out, n):
    from gedi.step2.outcomes import judge
    res = {}
    for model in MODELS:
        path = out / f"{model}.csv"
        if not path.exists():
            continue
        d = load(path)
        if n:
            keep = np.argsort(d["seed"])[:n]
            d = {k: v[keep] for k, v in d.items()}
        r = {"n": len(d["seed"]), "alone": np.nanmean(d["one_a_det"] == 1), "with": np.nanmean(d["one_w_det"] == 1)}
        for an in ("one", "best"):
            o = judge(d, an)
            r[an] = {k: float(np.mean(o[k])) for k in ("target", "companion", "wrong", "none", "false_alarm")}
        res[model] = r
    common_n = min(r["n"] for r in res.values())
    L = ["=== STEP 2 HEADLINE NUMBERS UNDER DIFFERENT NOISE MODELS ===",
         f"first {common_n} step-2 systems (same planets and epochs); fits use the reported (baseline) errors.",
         "Percent of all systems: target detected alone -> with the companion (one-orbit fit); then for the one-orbit",
         "fit and the BIC choice: target correct / companion correct / wrong planet / nothing [false alarm].", ""]
    for model, r in res.items():
        L.append(f"  {model:<13} (n={r['n']}): detected {100 * r['alone']:.1f}% -> {100 * r['with']:.1f}%")
        for an, lab in (("one", "one orbit"), ("best", "best by BIC")):
            o = r[an]
            L.append(f"      {lab:<12} {100 * o['target']:5.1f} / {100 * o['companion']:4.1f} / {100 * o['wrong']:5.1f} / "
                     f"{100 * o['none']:5.1f}   [{100 * o['false_alarm']:.1f}]")
    C.write_summary(out / "summary.txt", L)
    plt, ps = C.setup_plots()
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6), constrained_layout=True)
    x = np.arange(len(res))
    for k, (an, lab) in enumerate((("one", "one-orbit fit"), ("best", "best model by BIC"))):
        for j, (key, col, name) in enumerate((("target", ps.C_TARGET, "target correct"),
                                              ("wrong", ps.C_WRONG, "wrong planet"),
                                              ("false_alarm", ps.C_NONE, "false alarm (companion removed)"))):
            ax[k].bar(x + (j - 1) * 0.26, [res[m][an][key] for m in res], 0.25, color=col, label=name)
        ax[k].set_xticks(x, list(res), rotation=20)
        ps.percent_axis(ax[k], "y", "share of all systems")
        ax[k].set_title(lab, fontsize=10.5)
        ax[k].legend(fontsize=8)
    fig.suptitle("Step 2 results under less ideal noise models (fits use the reported errors)", fontsize=12)
    fig.savefig(out / "noise_robustness.png", dpi=150, bbox_inches="tight")
    print(f"saved {out / 'noise_robustness.png'}")


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("-n", type=int, default=1000)
    p.add_argument("-j", "--jobs", type=int, default=-1)
    p.add_argument("--models", nargs="+", default=MODELS, choices=MODELS)
    p.add_argument("--analyse", action="store_true")
    p.add_argument("--synthetic", action="store_true", help="synthetic epochs (code test only)")
    a = p.parse_args()
    out = C.out_dir("noise_robustness")
    if not a.analyse:
        for m in a.models:
            run(m, a.n, a.jobs, a.synthetic, out)
    analyse(out, a.n)


if __name__ == "__main__":
    main()
