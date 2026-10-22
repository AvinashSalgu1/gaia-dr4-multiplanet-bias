"""Can DR4 itself tell us how often giant-planet hosts have a second giant?

Uses only the step 3 results (no new fits, about two minutes):

    python scripts/extended/g1_multiplicity_signal.py [--root results/step3_gaia]

The multiplicity scan of step 3 reweights the fitted host stars to populations in which
a fraction f of giant-planet hosts have two or more giants (same true occurrence rate).
Here the same reweighting is used to ask what a DR4 user could MEASURE about f: among
the stars with a reported planet, the share whose adopted solution (lowest BIC) is not
a plain orbit (one orbit + acceleration, + jerk, or two orbits), and the share with
two reported orbits. If these shares change with f, the observed mix of solution types
in DR4 constrains f, and with it the multiplicity correction of the occurrence rate.

For each f: the shares with bootstrap 16-84% ranges, and the precision on f that a
catalogue of N reported planets would give from the binomial error of the share alone
(sigma_f = sigma_share / |d share / d f|), for N = 1,000-10,000. This ignores the
model uncertainty of the simulation itself (the BIC rule, the noise model, the planet
population), which is why the joint rate-multiplicity fit is tested on simulated
catalogues before it is applied to DR4.
"""
from __future__ import annotations

import argparse

import numpy as np

import common as C

F_GRID = np.round(np.arange(0.0, 0.81, 0.1), 2)
N_CAT = [1000, 3000, 10000]


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--root", type=C.Path, default=C.STEP3)
    p.add_argument("--boot", type=int, default=200)
    a = p.parse_args()
    import step3 as S3
    from gaiabias import multiplicity as MU
    from gaiabias.step2.outcomes import bic_choice

    z, meta, jobs, d = C.load_step3(a.root)
    out = C.out_dir("multiplicity_signal")
    res = S3.analyse_one(z, meta, jobs, d, "best")
    ch = bic_choice(d, "w")
    pos = {j: k for k, j in enumerate(d["job"])}
    # per pool host (both scenarios, in build_pool order): adopted model and number of reported planets
    hosts = {sc: res["scen"][sc]["hosts"] for sc in S3.SCENARIOS}
    pool = MU.build_pool(hosts)
    model, n_rep = [], []
    for sc in S3.SCENARIOS:
        npl = hosts[sc]["npl"]
        hs = np.nonzero(npl > 0)[0]
        job_of = {int(j["star"]): j["job"] for j in jobs if j["kind"] == "A" and j["scen"] == sc}
        nr = np.bincount(hosts[sc]["A_star"], minlength=len(npl))
        for s in hs:
            k = pos.get(job_of.get(int(s), ""), None)
            model.append(int(ch[k]) if k is not None else -1)
            n_rep.append(int(nr[s]))
    model, n_rep = np.array(model), np.array(n_rep)
    det = n_rep > 0
    alt = np.isin(model, [4, 5, 6])          # one orbit + accel, + jerk, two orbits
    two = model == 6
    edges = np.quantile(pool.m_star, np.linspace(0, 1, S3.N_MBIN + 1))
    edges[-1] += 1e-9
    hbin = np.clip(np.digitize(pool.m_star, edges) - 1, 0, S3.N_MBIN - 1)
    rng = np.random.default_rng(3)
    shares = {k: np.full((len(F_GRID), a.boot + 1), np.nan) for k in ("alt", "two", "multi_det")}
    f_real = np.zeros(len(F_GRID))
    for i, f in enumerate(F_GRID):
        w, mu, lost, fr = MU.host_weights_by_mass(pool, f, hbin)
        f_real[i] = fr
        for b in range(a.boot + 1):
            g = np.ones(len(w)) if b == 0 else rng.poisson(1.0, len(w))
            ww = w * g * det
            if ww.sum() <= 0:
                continue
            shares["alt"][i, b] = np.sum(ww * alt) / ww.sum()
            shares["two"][i, b] = np.sum(ww * two) / ww.sum()
            shares["multi_det"][i, b] = np.sum(ww * (pool.k >= 2)) / ww.sum()
    L = ["=== DOES THE MIX OF DR4 SOLUTION TYPES MEASURE MULTIPLICITY? ===",
         f"step 3 folder: {a.root}; best model by BIC; hosts of both scenarios reweighted to a fraction f of giant "
         "hosts with 2+ giants (same true occurrence rate).",
         "share_alt = stars with a reported planet whose adopted model is orbit + acceleration, + jerk or two orbits;",
         "share_two = ... whose adopted model has two orbits. Ranges: bootstrap 16-84% over hosts.", "",
         "    f    share_alt              share_two              (reported hosts with 2+ giants)"]
    for i in range(len(F_GRID)):
        r = lambda k: (f"{shares[k][i, 0]:.3f} [{np.nanpercentile(shares[k][i, 1:], 16):.3f},"  # noqa: E731
                       f"{np.nanpercentile(shares[k][i, 1:], 84):.3f}]")
        L.append(f"   {f_real[i]:.2f}  {r('alt'):<22} {r('two'):<22} {shares['multi_det'][i, 0]:.2f}")
    sens = {}
    for key in ("alt", "two"):
        s0 = shares[key][:, 0]
        slope = np.gradient(s0, f_real)
        sens[key] = slope
        L.append("")
        L.append(f"precision on f from share_{key} alone (binomial error, N reported planets in the catalogue):")
        for i in range(0, len(F_GRID), 2):
            parts = []
            for n in N_CAT:
                sig = np.sqrt(s0[i] * (1 - s0[i]) / n) / max(abs(slope[i]), 1e-9)
                parts.append(f"N={n:,}: +/-{sig:.3f}")
            L.append(f"   f = {f_real[i]:.2f}: slope {slope[i]:+.3f};  " + ";  ".join(parts))
    L += ["", "The simulation's own uncertainty (model selection rule, noise model, planet population) is not in these",
          "numbers; it is what the joint rate-multiplicity fit on simulated catalogues has to measure."]
    C.write_summary(out / "summary.txt", L)

    plt, ps = C.setup_plots()
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8), constrained_layout=True)
    for key, col, lab in (("alt", ps.C_WRONG, "adopted model is not a plain orbit"),
                          ("two", ps.C_TARGET, "adopted model has two orbits")):
        lo = np.nanpercentile(shares[key][:, 1:], 16, axis=1)
        hi = np.nanpercentile(shares[key][:, 1:], 84, axis=1)
        ax[0].fill_between(100 * f_real, lo, hi, color=col, alpha=0.2, lw=0)
        ax[0].plot(100 * f_real, shares[key][:, 0], "o-", color=col, label=lab)
    ax[0].set_xlabel("giant-planet hosts with a second giant (%)")
    ps.percent_axis(ax[0], "y", "share of stars with a reported planet")
    ax[0].set_ylim(0, None)
    ax[0].legend(fontsize=8.5)
    ax[0].set_title("What a DR4 user sees as multiplicity changes", fontsize=10.5)
    nn = np.geomspace(300, 30000, 40)
    for i, ls in ((1, "--"), (6, "-")):
        s0, sl = shares["alt"][i, 0], sens["alt"][i]
        ax[1].plot(nn, np.sqrt(s0 * (1 - s0) / nn) / abs(sl), ls, color=ps.C_WRONG,
                   label=f"true f = {f_real[i]:.2f}")
    ax[1].set_xscale("log")
    ax[1].set_yscale("log")
    ax[1].set_xlabel("number of reported planets in the catalogue")
    ax[1].set_ylabel("statistical uncertainty on f (from the solution mix)")
    ax[1].legend(fontsize=8.5)
    ax[1].set_title("Statistical precision only (no model uncertainty)", fontsize=10.5)
    fig.savefig(out / "multiplicity_signal.png", dpi=150, bbox_inches="tight")
    print(f"saved {out / 'multiplicity_signal.png'}")


if __name__ == "__main__":
    main()
