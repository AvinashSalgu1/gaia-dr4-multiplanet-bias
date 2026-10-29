"""Population-mix test: do the per-range corrections depend on the assumed mass and period mix?

Uses only the step 3 results (no new fits):

    python scripts/extended/g1_mix_test.py [--root results/step3_allsky]

The step 3 simulation draws its giants from the California Legacy Survey (CLS) distribution, which
fixes both the overall rate and the mix of masses and periods. The multiplicity scan holds the CLS
in-window rate fixed and reports inferred / true rates; g1_scan_cells.py splits these by mass and
period. Those ratios are the corrections a DR4 user would apply, but they are averages over the
assumed mix inside each range and they feed a power law fitted over the whole window.

This script tilts the simulated population and recomputes the same ratios. Each host gets a factor

    t = prod over its true planets of (m / m0)^dalpha * (P / P0)^dbeta,   m0 = 1.97 M_J, P0 = 0.89 yr

(the geometric centres of the window), so dalpha > 0 gives more massive giants and dbeta > 0 more
long-period giants. The factor is renormalised within every (stellar-mass bin, number of giants)
class, so the multiplicity distribution imposed at each f is unchanged, and the in-window true rate
is then fixed per stellar-mass bin as in step 3. Same fits, same detections, same completeness;
only the weights change.

If inferred / true in each range stays close to the CLS case, the corrections are a property of
Gaia and the fitting, not of the assumed population. If it moves, the mix has to be fitted
together with the rate and the multiplicity from the DR4 data.
Output: results/extended/mix_test/ (summary.txt, mix_test.csv, mix_test_mass.png, mix_test_period.png).
"""
from __future__ import annotations

import argparse
import csv

import numpy as np

import common as C

M_LAB = ["0.3–1 M_J", "1–3 M_J", "3–13 M_J"]
P_LAB = ["0.1–0.3 yr", "0.3–1 yr", "1–2 yr", "2–4 yr", "4–8 yr"]
M0, P0 = np.sqrt(0.3 * 13.0), np.sqrt(0.1 * 8.0)
# (name, dalpha, dbeta); the CLS case first
TILTS = [("CLS", 0.0, 0.0),
         ("mass -0.6", -0.6, 0.0), ("mass -0.3", -0.3, 0.0), ("mass +0.3", 0.3, 0.0), ("mass +0.6", 0.6, 0.0),
         ("period -0.6", 0.0, -0.6), ("period -0.3", 0.0, -0.3), ("period +0.3", 0.0, 0.3),
         ("period +0.6", 0.0, 0.6)]
WELL = (3.0, 1.0)            # well-measured block: m >= 3 M_J and P >= 1 yr


def _int_pow(k, lo, hi):
    return np.log(hi / lo) if abs(k) < 1e-9 else (hi ** k - lo ** k) / k


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--root", default=str(C.STEP3))
    a = p.parse_args()
    import step3 as S3
    from gaiabias import multiplicity as MU
    from gaiabias import rates as R
    from pathlib import Path

    z, meta, jobs, d = C.load_step3(Path(a.root))
    out = C.out_dir("mix_test")
    n_star = meta["n_stars"]
    m_win, P_win = (S3.M_FINE[0], S3.M_FINE[-1]), (S3.P_FINE[0], S3.P_FINE[-1])
    ME, PE = S3.M_EDGES, S3.P_EDGES
    nm, nP = len(ME) - 1, len(PE) - 1
    well_cells = (ME[:-1, None] >= WELL[0]) & (PE[None, :-1] >= WELL[1])
    F = S3.SCAN_F
    nt = len(TILTS)
    res, share = {}, {}
    for an in S3.SCAN_METHODS:
        print(f"analysis {an} ...", flush=True)
        r = S3.analyse_one(z, meta, jobs, d, an)
        hosts = {sc: r["scen"][sc]["hosts"] for sc in S3.SCENARIOS}
        pool = MU.build_pool(hosts)
        h0 = hosts["independent"]
        n_hosts = len(pool.k)
        ms_all = h0["m_star"]
        edges = np.quantile(ms_all, np.linspace(0, 1, S3.N_MBIN + 1))
        edges[-1] += 1e-9
        star_bin = np.clip(np.digitize(ms_all, edges) - 1, 0, S3.N_MBIN - 1)
        frac_bin = np.bincount(star_bin, minlength=S3.N_MBIN) / n_star
        host_bin = np.clip(np.digitize(pool.m_star, edges) - 1, 0, S3.N_MBIN - 1)
        win0 = MU.in_window(h0["pm"], h0["pP"], m_win, P_win)
        lam_bin = np.bincount(star_bin[h0["pst"][win0]], minlength=S3.N_MBIN) / np.maximum(
            np.bincount(star_bin, minlength=S3.N_MBIN), 1)
        t_in = MU.in_window(pool.t_m, pool.t_P, m_win, P_win)
        t_win = np.bincount(pool.t_host, weights=t_in, minlength=n_hosts)
        eta = r["eta"]
        region = eta >= S3.MIN_ETA
        inj = r["inj"]
        Rm = R.response_matrix(inj["m_all"], inj["P_all"], inj["t_m"], inj["t_P"], inj["f_m"], inj["f_P"],
                               S3.M_FINE, S3.P_TRUE_EDGES, S3.M_FINE, S3.P_FINE)
        cls_id = host_bin * (int(pool.k.max()) + 1) + pool.k             # (stellar-mass bin, k) class
        q = {k: np.full((nt, len(F), nm, nP), np.nan) for k in ("plA", "plB", "binA", "true")}
        sh = np.zeros((nt, nm, nP))
        for it, (name, da, db) in enumerate(TILTS):
            lt = da * np.log(pool.t_m / M0) + db * np.log(pool.t_P / P0)
            tilt = np.exp(np.bincount(pool.t_host, weights=lt, minlength=n_hosts))
            mean_c = np.bincount(cls_id, weights=tilt) / np.maximum(np.bincount(cls_id), 1)
            tilt = tilt / mean_c[cls_id]
            print(f"  {name} ...", flush=True)
            for i, f in enumerate(F):
                v0, _, _, _ = MU.host_weights_by_mass(pool, f, host_bin)
                v0 = v0 * tilt
                ew = np.bincount(host_bin, weights=v0 * t_win, minlength=S3.N_MBIN)
                with np.errstate(invalid="ignore", divide="ignore"):
                    Hj = np.where(ew > 0, lam_bin / ew, 0.0) * frac_bin
                v = v0 * Hj[host_bin]
                q["true"][it, i] = R.counts_2d(pool.t_m, pool.t_P, ME, PE, weights=v[pool.t_host])
                if i == 0:
                    sh[it] = q["true"][it, 0] / np.nansum(q["true"][it, 0])
                for tag, mm, PP, hh in (("A", pool.a_m, pool.a_P, pool.a_host),
                                        ("B", pool.b_m, pool.b_P, pool.b_host)):
                    if tag == "A":
                        c = R.counts_2d(mm, PP, ME, PE, weights=v[hh])
                        with np.errstate(invalid="ignore", divide="ignore"):
                            q["binA"][it, i] = np.where(region, c / eta, np.nan)
                    try:
                        fit = R.fit_power_law_forward(mm, PP, n_star, Rm, S3.M_FINE, S3.P_TRUE_EDGES, S3.M_FINE,
                                                      S3.P_FINE, weights=v[hh] * n_star, errors=False)
                        Cc = np.exp(fit.lnC)
                        for jm in range(nm):
                            for jp in range(nP):
                                q[f"pl{tag}"][it, i, jm, jp] = Cc * _int_pow(fit.alpha, ME[jm], ME[jm + 1]) * \
                                    _int_pow(fit.beta, PE[jp], PE[jp + 1])
                    except Exception as e:                           # noqa: BLE001
                        print(f"  power law failed: {name}, f={f} ({tag}): {e}", flush=True)
        q["region"] = region
        res[an], share[an] = q, sh

    # ratios per range: by mass, by period, whole window, well-measured block
    def ratios(q):
        reg = q["region"]
        o = {}
        sums = {"mass": lambda x: np.nansum(x, axis=-1), "period": lambda x: np.nansum(x, axis=-2),
                "window": lambda x: np.nansum(x, axis=(-1, -2))[..., None],
                "well": lambda x: np.nansum(np.where(well_cells, x, 0.0), axis=(-1, -2))[..., None]}
        for dim, fn in sums.items():
            tr = fn(q["true"])
            with np.errstate(invalid="ignore", divide="ignore"):
                g = {"plA": fn(q["plA"]) / tr, "plB": fn(q["plB"]) / tr}
                g["AB"] = g["plA"] / g["plB"]
                num = fn(np.where(reg, q["binA"], 0.0))
                den = fn(np.where(reg, q["true"], 0.0))
                g["binA"] = np.where(den > 0, num / den, np.nan)
            o[dim] = g
        return o

    RT = {an: ratios(q) for an, q in res.items()}
    LABS = {"mass": M_LAB, "period": P_LAB, "window": ["whole window"], "well": ["3–13 M_J, 1–8 yr"]}
    iF = [int(np.argmin(abs(F - x))) for x in (0.0, 0.4, 0.8)]
    fmt = lambda x: "   -" if not np.isfinite(x) else f"{x:4.2f}"  # noqa: E731
    L = ["=== POPULATION-MIX TEST: INFERRED / TRUE RATE UNDER TILTED MASS AND PERIOD MIXES ===",
         f"step 3 results: {a.root} ({n_star:,} stars); no new fits. Tilt factor per host = prod over its planets of "
         f"(m/{M0:.2f} M_J)^dalpha (P/{P0:.2f} yr)^dbeta, renormalised within each (stellar-mass bin, k) class so the "
         "multiplicity is unchanged; in-window true rate fixed per stellar-mass bin. A = real systems, "
         "B = planets alone, power law unless 'binned' (completeness >= 5% cells only).", "",
         "--- share of true in-window planets per range (f = 0) ---"]
    for an in res:
        sh = share[an]
        L.append(f"  {S3.LABEL[an]}")
        L.append("   tilt          " + "  ".join(f"{s:>10}" for s in M_LAB + P_LAB))
        for it, (name, _, _) in enumerate(TILTS):
            vals = list(sh[it].sum(axis=1)) + list(sh[it].sum(axis=0))
            L.append(f"   {name:<13} " + "  ".join(f"{v:10.3f}" for v in vals))
        break                                                        # same truth for both analyses
    L.append("")
    rows = []
    for an in res:
        for dim in ("mass", "period", "window", "well"):
            g = RT[an][dim]
            for j, lab in enumerate(LABS[dim]):
                L.append(f"--- {S3.LABEL[an]}, {lab} ---")
                L.append("   tilt         " + "  ".join(f"A/true f={F[i]:.1f}" for i in iF) + "   "
                         + "  ".join(f"A/B f={F[i]:.1f}" for i in iF) + "   " + "  bin A/true f=0.8")
                for it, (name, _, _) in enumerate(TILTS):
                    L.append(f"   {name:<12} " + "  ".join(f"{fmt(g['plA'][it, i, j]):>12}" for i in iF) + "   "
                             + "  ".join(f"{fmt(g['AB'][it, i, j]):>9}" for i in iF) + "   "
                             + f"{fmt(g['binA'][it, iF[2], j]):>16}")
                    for i, f in enumerate(F):
                        rows.append(dict(analysis=an, tilt=name, dalpha=TILTS[it][1], dbeta=TILTS[it][2], by=dim,
                                         range=lab, f=f, pl_A_over_true=g["plA"][it, i, j],
                                         pl_B_over_true=g["plB"][it, i, j], pl_A_over_B=g["AB"][it, i, j],
                                         bin_A_over_true=g["binA"][it, i, j]))
                # spread over tilts at f = 0.8, relative to CLS
                k = iF[2]
                with np.errstate(invalid="ignore", divide="ignore"):
                    rel = g["plA"][:, k, j] / g["plA"][0, k, j]
                    relb = g["binA"][:, k, j] / g["binA"][0, k, j]
                L.append(f"   spread at f=0.8 relative to CLS: power law {np.nanmin(rel):.2f}–{np.nanmax(rel):.2f}"
                         + (f", binned {np.nanmin(relb):.2f}–{np.nanmax(relb):.2f}" if np.isfinite(relb).any() else ""))
                L.append("")
    with open(out / "mix_test.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    C.write_summary(out / "summary.txt", L)

    plt, ps = C.setup_plots()
    cmap = plt.get_cmap("coolwarm")
    for kind, dim, labs, sel in (("mass", "mass", M_LAB, [0, 1, 2, 3, 4]), ("period", "period", P_LAB, [0, 5, 6, 7, 8])):
        ncol = len(labs) + 2
        fig, ax = plt.subplots(len(res), ncol, figsize=(3.6 * ncol, 4.2 * len(res) + 0.8), sharex=True,
                               squeeze=False, constrained_layout=True)
        for ia, an in enumerate(res):
            for j in range(ncol):
                a_ = ax[ia, j]
                if j < len(labs):
                    g, jj, title = RT[an][dim], j, labs[j]
                elif j == len(labs):
                    g, jj, title = RT[an]["window"], 0, "whole window"
                else:
                    g, jj, title = RT[an]["well"], 0, "3–13 M_J, 1–8 yr"
                for it in sel:
                    name, da, db = TILTS[it]
                    s = da if kind == "mass" else db
                    c = "0.1" if it == 0 else cmap(0.5 + s / 1.4)
                    lw = 2.4 if it == 0 else 1.4
                    a_.plot(F * 100, g["plA"][it, :, jj], "-", color=c, lw=lw, label=f"{name}: power law")
                    if np.isfinite(g["binA"][it, :, jj]).any():
                        a_.plot(F * 100, g["binA"][it, :, jj], "s", ms=4, mfc="white", mec=c,
                                label=f"{name}: binned")
                a_.axhline(1, color="0.5", lw=0.8, ls=":")
                a_.set_yscale("log")
                a_.set_ylim(0.25, 2)
                if ia == 0:
                    a_.set_title(title, fontsize=12)
                if j == 0:
                    a_.set_ylabel(f"{S3.LABEL[an]}\nreal systems ÷ truth")
                if ia == len(res) - 1:
                    a_.set_xlabel("hosts with a second giant (%)")
        h, l = ax[0, -1].get_legend_handles_labels()
        fig.legend(h, l, loc="outside right center", fontsize=9.5, frameon=False)
        what = "mass slope (dalpha; + = more massive giants)" if kind == "mass" else \
            "period slope (dbeta; + = more long-period giants)"
        fig.suptitle(f"Population-mix test: inferred ÷ true rate when the CLS {what} is tilted (black = CLS)",
                     fontsize=13)
        fig.savefig(out / f"mix_test_{kind}.png", dpi=140)
        plt.close(fig)
    print(f"saved {out}/mix_test_mass.png and mix_test_period.png")


if __name__ == "__main__":
    main()
