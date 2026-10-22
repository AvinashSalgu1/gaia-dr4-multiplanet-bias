"""Warning flags chosen on one half of the stars and tested on the other half.

Uses only the step 3 results (no new fits, about a minute):

    python scripts/extended/g1_flags_heldout.py [--root results/step3_gaia]

The step 4 flags (adopted model not a plain orbit; significance < 5; period > 0.7 x time
span; eccentricity > 0.5) were measured on the same simulated planets used to choose
them. Here the stars are split at random into a training half and a test half:

* the thresholds of the three continuous flags are chosen on the training half
  (maximum of recall - false-positive rate for "the reported planet is not accurate");
* the default rule and the tuned rule are then scored on the test half only:
  recall (share of inaccurate planets flagged), precision (share of flagged planets
  that are inaccurate), false-positive rate (share of accurate planets flagged) and
  the share of planets flagged, with bootstrap 16-84% ranges over test stars;
* calibration: the probability that a planet is inaccurate given how many flags it
  raises, estimated on the training half, compared with the test half.

Planet-only populations (independent and clustered multiplicity). Stellar and
brown-dwarf companions are not in these simulations, so the numbers apply to
planet-only systems.
"""
from __future__ import annotations

import argparse
import itertools

import numpy as np

import common as C

T_YR = 1977.0 / 365.25
SIG_GRID = [3, 4, 5, 6, 8, 10]
PF_GRID = [0.5, 0.6, 0.7, 0.8, 0.9]
E_GRID = [0.3, 0.4, 0.5, 0.6, 0.7]
DEFAULT = (5.0, 0.7, 0.5)


def flags(t, sig_cut, pf_cut, e_cut):
    """Boolean array per flag (n_flags, n_planets)."""
    return np.array([t["bic"] != 3,
                     np.nan_to_num(t["sig"], nan=0) < sig_cut,
                     t["P"] > pf_cut * T_YR,
                     np.nan_to_num(t["e1"], nan=0) > e_cut])


def scores(flagged, bad):
    tp, fp = np.sum(flagged & bad), np.sum(flagged & ~bad)
    fn, tn = np.sum(~flagged & bad), np.sum(~flagged & ~bad)
    return dict(recall=tp / max(tp + fn, 1), precision=tp / max(tp + fp, 1), fpr=fp / max(fp + tn, 1),
                share=flagged.mean(), n=len(bad))


def boot(star, flagged, bad, n_boot=300, seed=1):
    rng = np.random.default_rng(seed)
    us = np.unique(star)
    idx_of = {s: np.nonzero(star == s)[0] for s in us}
    out = {k: [] for k in ("recall", "precision", "fpr", "share")}
    for _ in range(n_boot):
        pick = np.concatenate([idx_of[s] for s in rng.choice(us, len(us))])
        sc = scores(flagged[pick], bad[pick])
        for k in out:
            out[k].append(sc[k])
    return {k: np.percentile(v, [16, 84]) for k, v in out.items()}


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--root", type=C.Path, default=C.STEP3)
    p.add_argument("--seed", type=int, default=2027)
    a = p.parse_args()
    import step4 as S4
    import step3 as S3

    z, meta, jobs, d = C.load_step3(a.root)
    out = C.out_dir("flags_heldout")
    L = ["=== WARNING FLAGS ON HELD-OUT STARS ===",
         f"step 3 folder: {a.root}; best model by BIC (DR4-like); stars split 50/50 at random (seed {a.seed}).",
         "'Not accurate' = the reported planet is imprecise or spurious; flags as in step 4.",
         "recall = share of inaccurate planets flagged; precision = share of flagged planets that are "
         "inaccurate; FPR = share of accurate planets flagged.", ""]
    res = {}
    for scen in S3.SCENARIOS:
        t = S4.planet_table(z, meta, jobs, d, "best", scen, 2.0)
        bad = t["cls"] >= 1
        rng = np.random.default_rng([a.seed, S3.SCENARIOS.index(scen)])
        us = np.unique(t["star"])
        train_stars = set(rng.choice(us, len(us) // 2, replace=False).tolist())
        tr = np.array([s in train_stars for s in t["star"]])
        te = ~tr
        best, best_j = None, -1
        for cut in itertools.product(SIG_GRID, PF_GRID, E_GRID):
            f = flags(t, *cut).any(axis=0)
            sc = scores(f[tr], bad[tr])
            j = sc["recall"] - sc["fpr"]
            if j > best_j:
                best, best_j = cut, j
        rows = {}
        for name, cut in (("default", DEFAULT), ("tuned on training half", best)):
            f = flags(t, *cut).any(axis=0)
            rows[name] = dict(train=scores(f[tr], bad[tr]), test=scores(f[te], bad[te]),
                              ci=boot(t["star"][te], f[te], bad[te]), cut=cut)
        # calibration by number of flags raised (default thresholds)
        nflag = flags(t, *DEFAULT).sum(axis=0)
        cal = []
        for k in range(5):
            sel = nflag == k if k < 4 else nflag >= 4
            cal.append((k, bad[tr & sel].mean() if (tr & sel).sum() >= 10 else np.nan,
                        bad[te & sel].mean() if (te & sel).sum() >= 10 else np.nan, int((te & sel).sum())))
        # hidden giant as the target (default rule)
        hid = t["n_hidden"] > 0
        f = flags(t, *DEFAULT).any(axis=0)
        rows["default, target = star hides another giant"] = dict(
            train=scores(f[tr], hid[tr]), test=scores(f[te], hid[te]), ci=boot(t["star"][te], f[te], hid[te]),
            cut=DEFAULT)
        res[scen] = dict(rows=rows, cal=cal, base=bad[te].mean())
        L.append(f"--- {scen} multiplicity: {len(t['P'])} reported planets, {len(us)} stars; inaccurate on "
                 f"test half {100 * bad[te].mean():.0f}% ---")
        for name, r in rows.items():
            c = r["cut"]
            L.append(f"  {name}: significance < {c[0]:g}, period > {c[1]:g} x span, eccentricity > {c[2]:g}")
            for half in ("train", "test"):
                s = r[half]
                extra = ""
                if half == "test":
                    ci = r["ci"]
                    extra = (f"   [16-84%: recall {100 * ci['recall'][0]:.0f}-{100 * ci['recall'][1]:.0f}, "
                             f"precision {100 * ci['precision'][0]:.0f}-{100 * ci['precision'][1]:.0f}, "
                             f"FPR {100 * ci['fpr'][0]:.0f}-{100 * ci['fpr'][1]:.0f}]")
                L.append(f"    {half:<5} recall {100 * s['recall']:5.1f}%  precision {100 * s['precision']:5.1f}%  "
                         f"FPR {100 * s['fpr']:5.1f}%  flagged {100 * s['share']:5.1f}%  (n={s['n']}){extra}")
        L.append("  calibration (default thresholds): P(inaccurate | number of flags) train -> test (n test)")
        for k, ptr, pte, n in cal:
            L.append(f"    {k if k < 4 else '4':>2}{'+' if k == 4 else ' '} flags: {ptr:.2f} -> {pte:.2f}  ({n})")
        L.append("")
    L.append("Planet-only simulations: no stellar or brown-dwarf companions are included.")
    C.write_summary(out / "summary.txt", L)

    plt, ps = C.setup_plots()
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8), constrained_layout=True)
    metrics = ["recall", "precision", "fpr"]
    x = np.arange(len(metrics))
    for i, scen in enumerate(S3.SCENARIOS):
        for j, name in enumerate(["default", "tuned on training half"]):
            r = res[scen]["rows"][name]
            v = [r["test"][m] for m in metrics]
            lo = [r["test"][m] - r["ci"][m][0] for m in metrics]
            hi = [r["ci"][m][1] - r["test"][m] for m in metrics]
            off = (2 * i + j - 1.5) * 0.19
            ax[0].bar(x + off, v, 0.18, yerr=[lo, hi], capsize=2,
                      color=[ps.C_TARGET, "#8fb8f0", ps.C_WRONG, "#f3b98c"][2 * i + j],
                      label=f"{scen}, {name}")
            tv = [r["train"][m] for m in metrics]
            ax[0].scatter(x + off, tv, marker="_", s=120, color=ps.INK, zorder=3,
                          label="same rule on training half" if (i, j) == (0, 0) else None)
    ax[0].set_xticks(x, ["recall:\ninaccurate\nplanets flagged", "precision:\nflagged planets\nthat are inaccurate",
                         "false-positive rate:\naccurate planets\nflagged"], fontsize=9)
    ps.percent_axis(ax[0], "y", "test half")
    ax[0].set_ylim(0, 1)
    ax[0].legend(fontsize=8, loc="upper right")
    ax[0].set_title("Suspect rule scored on stars not used to choose it", fontsize=10.5)
    for i, scen in enumerate(S3.SCENARIOS):
        cal = res[scen]["cal"]
        k = [c[0] for c in cal]
        ax[1].plot(k, [c[1] for c in cal], "o--", color=[ps.C_TARGET, ps.C_WRONG][i], mfc="white",
                   label=f"{scen}: training half")
        ax[1].plot(k, [c[2] for c in cal], "o-", color=[ps.C_TARGET, ps.C_WRONG][i], label=f"{scen}: test half")
    ax[1].set_xticks(range(5), ["0", "1", "2", "3", "4"])
    ax[1].set_xlabel("number of warning flags raised")
    ps.percent_axis(ax[1], "y", "share of planets that are not accurate")
    ax[1].set_ylim(0, 1)
    ax[1].legend(fontsize=8)
    ax[1].set_title("Calibration: chance of an inaccurate orbit vs flags raised", fontsize=10.5)
    fig.savefig(out / "flags_heldout.png", dpi=150, bbox_inches="tight")
    print(f"saved {out / 'flags_heldout.png'}")


if __name__ == "__main__":
    main()
