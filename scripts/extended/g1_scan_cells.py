"""Occurrence-rate bias against multiplicity, resolved by planet mass and by orbital period.

Uses only the step 3 results (no new fits):

    python scripts/extended/g1_scan_cells.py [--root results/step3_allsky]

The step 3 multiplicity scan reports one number per multiplicity f (the share of giant-planet
hosts with two or more giants): the rate over the whole window (0.3-13 M_J, 0.1-8 yr) or over the
well-measured cells. This script splits the same scan by mass (0.3-1, 1-3, 3-13 M_J) and by
period (0.1-0.3, 0.3-1, 1-2, 2-4, 4-8 yr), with the same host reweighting as step 3 (nominal
weights, true occurrence rate fixed at every f):

  power law    the forward-modelled power law fitted over the whole window (as in step 3),
               integrated over each mass or period range;
  binned       detections / completeness in the cells where the completeness is >= 5% (the only
               cells that can be measured without assuming a shape), summed over the range;
  truth        the true planets in the same range.

Each is computed for the real systems (A) and for the same planets fitted alone (B), for the
one-orbit fit and the BIC choice. A / truth is what a DR4 user would get; B / truth is the error
of the method with no multiplicity at all; A / B isolates the effect of multiplicity.
Output: results/extended/scan_cells/ (summary.txt, scan_cells.csv, scan_cells_lines.png,
scan_cells_heatmap.png).
"""
from __future__ import annotations

import argparse
import csv

import numpy as np

import common as C

M_LAB = ["0.3–1 M_J", "1–3 M_J", "3–13 M_J"]
P_LAB = ["0.1–0.3 yr", "0.3–1 yr", "1–2 yr", "2–4 yr", "4–8 yr"]


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
    out = C.out_dir("scan_cells")
    n_star = meta["n_stars"]
    m_win, P_win = (S3.M_FINE[0], S3.M_FINE[-1]), (S3.P_FINE[0], S3.P_FINE[-1])
    ME, PE = S3.M_EDGES, S3.P_EDGES
    nm, nP = len(ME) - 1, len(PE) - 1
    F = S3.SCAN_F
    # res[an][quantity] arrays (len(F), nm, nP) per star; quantities: plA, plB, binA, binB, true
    res = {}
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
        t_win = np.bincount(pool.t_host, weights=MU.in_window(pool.t_m, pool.t_P, m_win, P_win), minlength=n_hosts)
        eta = r["eta"]
        region = eta >= S3.MIN_ETA
        inj = r["inj"]
        Rm = R.response_matrix(inj["m_all"], inj["P_all"], inj["t_m"], inj["t_P"], inj["f_m"], inj["f_P"],
                               S3.M_FINE, S3.P_TRUE_EDGES, S3.M_FINE, S3.P_FINE)
        q = {k: np.full((len(F), nm, nP), np.nan) for k in ("plA", "plB", "binA", "binB", "true")}
        for i, f in enumerate(F):
            v0, _, _, _ = MU.host_weights_by_mass(pool, f, host_bin)
            ew = np.bincount(host_bin, weights=v0 * t_win, minlength=S3.N_MBIN)
            with np.errstate(invalid="ignore", divide="ignore"):
                Hj = np.where(ew > 0, lam_bin / ew, 0.0) * frac_bin
            v = v0 * Hj[host_bin]                                    # per-star weights of each host
            q["true"][i] = R.counts_2d(pool.t_m, pool.t_P, ME, PE, weights=v[pool.t_host])
            for tag, mm, PP, hh in (("A", pool.a_m, pool.a_P, pool.a_host), ("B", pool.b_m, pool.b_P, pool.b_host)):
                c = R.counts_2d(mm, PP, ME, PE, weights=v[hh])
                with np.errstate(invalid="ignore", divide="ignore"):
                    q[f"bin{tag}"][i] = np.where(region, c / eta, np.nan)
                try:
                    fit = R.fit_power_law_forward(mm, PP, n_star, Rm, S3.M_FINE, S3.P_TRUE_EDGES, S3.M_FINE,
                                                  S3.P_FINE, weights=v[hh] * n_star, errors=False)
                    Cc = np.exp(fit.lnC)
                    for jm in range(nm):
                        for jp in range(nP):
                            q[f"pl{tag}"][i, jm, jp] = Cc * _int_pow(fit.alpha, ME[jm], ME[jm + 1]) * \
                                _int_pow(fit.beta, PE[jp], PE[jp + 1])
                except Exception as e:                               # noqa: BLE001
                    print(f"  power law failed at f={f} ({tag}): {e}", flush=True)
        q["region"] = region
        res[an] = q

    # marginals: by mass (sum over periods) and by period (sum over masses); binned uses measurable cells only
    def marg(q, axis):
        reg = q["region"]
        o = {}
        for k in ("plA", "plB", "true"):
            o[k] = np.nansum(q[k], axis=axis)
        for tag in ("A", "B"):
            num = np.nansum(np.where(reg, q[f"bin{tag}"], 0.0), axis=axis)
            den = np.nansum(np.where(reg[None], q["true"], 0.0), axis=axis)
            with np.errstate(invalid="ignore", divide="ignore"):
                o[f"bin{tag}_ratio"] = np.where(den > 0, num / den, np.nan)
        return o

    rows, L = [], ["=== OCCURRENCE-RATE BIAS VS MULTIPLICITY, BY MASS AND BY PERIOD ===",
                   f"step 3 results: {a.root} ({n_star:,} stars). Ratios to the true rate in the same range; "
                   "A = real systems, B = planets alone; nominal weights (no bootstrap).",
                   "binned = detections / completeness, only cells with completeness >= 5% ('-' = none in that range).", ""]
    M = {}
    for an, q in res.items():
        M[an] = {"mass": marg(q, 2), "period": marg(q, 1)}
        for dim, labs in (("mass", M_LAB), ("period", P_LAB)):
            g = M[an][dim]
            L.append(f"--- {S3.LABEL[an]}, by {dim} ---")
            L.append("   range          f:   " + "  ".join(f"{f:4.1f}" for f in F))
            for j, lab in enumerate(labs):
                with np.errstate(invalid="ignore", divide="ignore"):
                    plA, plB = g["plA"][:, j] / g["true"][:, j], g["plB"][:, j] / g["true"][:, j]
                bA = g["binA_ratio"][:, j]
                fmt = lambda x: "  ".join("   -" if not np.isfinite(v) else f"{v:4.2f}" for v in x)  # noqa: E731
                L += [f"   {lab:<12} PL A/true  {fmt(plA)}",
                      f"   {'':<12} PL B/true  {fmt(plB)}",
                      f"   {'':<12} PL A/B     {fmt(plA / plB)}",
                      f"   {'':<12} bin A/true {fmt(bA)}"]
                for i, f in enumerate(F):
                    rows.append(dict(analysis=an, by=dim, range=lab, f=f, true_per_star=g["true"][i, j],
                                     pl_A_over_true=plA[i], pl_B_over_true=plB[i], pl_A_over_B=plA[i] / plB[i],
                                     bin_A_over_true=bA[i], bin_B_over_true=g["binB_ratio"][i, j]))
            L.append("")
    with open(out / "scan_cells.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    C.write_summary(out / "summary.txt", L)

    plt, ps = C.setup_plots()
    from matplotlib.colors import LogNorm
    # ---- figure 1: lines, one panel per mass range (top) and per period range (bottom)
    fig, ax = plt.subplots(2, 5, figsize=(22, 8.6), sharex=True)
    col = {"one": "#2a6fdb", "best": "#1a9e77"}
    for row, (dim, labs) in enumerate((("mass", M_LAB), ("period", P_LAB))):
        for j in range(5):
            a_ = ax[row, j]
            if j >= len(labs):
                a_.axis("off")
                continue
            for an in res:
                g = M[an][dim]
                with np.errstate(invalid="ignore", divide="ignore"):
                    a_.plot(F * 100, g["plA"][:, j] / g["true"][:, j], "-o", ms=4, color=col[an],
                            label=f"{S3.LABEL[an]}: real systems, power law")
                    a_.plot(F * 100, g["plB"][:, j] / g["true"][:, j], "--", color=col[an], alpha=0.8,
                            label=f"{S3.LABEL[an]}: planets alone, power law")
                bA = g["binA_ratio"][:, j]
                if np.isfinite(bA).any():
                    a_.plot(F * 100, bA, "s", ms=6, mfc="white", mec=col[an],
                            label=f"{S3.LABEL[an]}: real systems, binned (measured cells)")
            a_.axhline(1, color="0.2", lw=1)
            a_.set_yscale("log")
            a_.set_ylim(0.2, 5)
            a_.set_title(labs[j], fontsize=12)
            if j == 0:
                a_.set_ylabel(f"inferred ÷ true rate\n(by {dim})")
            if row == 1:
                a_.set_xlabel("giant-planet hosts with a second giant (%)")
    h, l = ax[0, 2].get_legend_handles_labels()
    ax[0, 4].legend(h, l, loc="center", fontsize=10.5, frameon=False)
    fig.suptitle("Occurrence-rate bias against multiplicity, split by planet mass (top) and orbital period (bottom); "
                 "1 = right", fontsize=13)
    fig.tight_layout()
    fig.savefig(out / "scan_cells_lines.png", dpi=150)
    plt.close(fig)

    # ---- figure 2: heatmaps, multiplicity on x, period (top) or mass (bottom) on y
    fig, ax = plt.subplots(2 * len(res), 3, figsize=(19, 6.2 * len(res) + 1), constrained_layout=True)
    norm = LogNorm(0.25, 4)
    cmap = plt.get_cmap("RdBu")
    titles = ["real systems ÷ truth\n(what DR4 would give)", "planets alone ÷ truth\n(error of the method alone)",
              "real systems ÷ planets alone\n(effect of multiplicity)"]
    for ia, (an, q) in enumerate(res.items()):
        for row, (dim, labs) in enumerate((("period", P_LAB), ("mass", M_LAB))):
            g = M[an][dim]
            with np.errstate(invalid="ignore", divide="ignore"):
                mats = [g["plA"] / g["true"], g["plB"] / g["true"], g["plA"] / g["plB"]]
            for k in range(3):
                a_ = ax[2 * ia + row, k]
                img = mats[k].T                                      # rows = ranges, columns = f
                im = a_.imshow(img, origin="lower", aspect="auto", cmap=cmap, norm=norm,
                               extent=[-5, 85, -0.5, len(labs) - 0.5])
                for jj in range(len(labs)):
                    for ii, f in enumerate(F):
                        val = img[jj, ii]
                        if np.isfinite(val):
                            a_.text(f * 100, jj, f"{val:.2f}", ha="center", va="center", fontsize=8.5,
                                    color="white" if abs(np.log(val)) > 0.8 else "0.1")
                a_.set_yticks(range(len(labs)), labs)
                a_.set_xticks(F * 100)
                if 2 * ia + row == 2 * len(res) - 1:
                    a_.set_xlabel("giant-planet hosts with a second giant (%)")
                if 2 * ia + row == 0:
                    a_.set_title(titles[k], fontsize=13)
                if k == 0:
                    a_.set_ylabel(f"{S3.LABEL[an]}\nby {dim}", fontsize=12)
    cb = fig.colorbar(im, ax=ax, shrink=0.5, pad=0.01)
    cb.set_label("ratio (1 = equal; blue above, red below)")
    fig.suptitle("Power-law occurrence rate split by orbital period and by planet mass, against multiplicity",
                 fontsize=15)
    fig.savefig(out / "scan_cells_heatmap.png", dpi=140)
    plt.close(fig)
    print(f"saved {out}/scan_cells_lines.png and scan_cells_heatmap.png")


if __name__ == "__main__":
    main()
