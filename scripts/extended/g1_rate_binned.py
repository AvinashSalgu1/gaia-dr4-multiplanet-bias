"""Occurrence rates without assuming a power law, and the injection response matrix.

Uses only the step 3 results (no new fits, a few minutes):

    python scripts/extended/g1_rate_binned.py [--root results/step3_gaia] [--boot 100]

Step 3 infers the rate over the whole giant-planet range by forward-modelling a power
law through the injection response. Here the same detections are fitted with a BINNED
rate model instead: one free occurrence rate per (mass, period) cell, no shape assumed.
The true cells match the fitted cells, except that the shortest- and longest-period
cells are widened (0.05-0.3 yr and 4-16 yr) so that planets scattering into the
window from outside are counted in the cell they enter; their rates are compared
with the truth over the same widened ranges. The expected number of reported planets in fitted cell k is

    N_k = N_star * sum_j R[k, j] r_j ,

with R[k, j] the share of injected planets from true cell j that are reported in fitted
cell k (measured from the single-planet injections), and the r_j are found by maximising
the Poisson likelihood (r_j >= 0). Cells to which DR4 has almost no sensitivity
(sum_k R[k, j] < MIN_RESP) cannot be measured without a shape assumption; they are
reported as unconstrained instead of being extrapolated.

Outputs: the response matrix figure (true cell -> fitted cell, and the share not
detected), and inferred / true per cell for real systems (A) and planets alone (B),
with bootstrap ranges (resampling detections and injections). The ratio A / B in the
constrained cells is the multiplicity effect free of any power-law extrapolation.
"""
from __future__ import annotations

import argparse

import numpy as np

import common as C

M_T = np.array([0.3, 1.0, 3.0, 13.0])
P_T = np.array([0.05, 0.3, 1.0, 2.0, 4.0, 16.0])   # one true cell per fitted cell; the edge cells also take in
                                                   # planets that drift into the window from outside
MIN_RESP = 0.01


def fit_binned(n_k, R, n_star, free):
    """Poisson ML rates (planets per star) for the true cells in `free`; others fixed at 0."""
    from scipy.optimize import minimize
    Rf = R[:, free]

    def nll(x):
        r = np.exp(x)
        mu = n_star * Rf @ r
        return -(np.sum(n_k * np.log(np.maximum(mu, 1e-300))) - mu.sum())

    def grad(x):
        r = np.exp(x)
        mu = np.maximum(n_star * Rf @ r, 1e-300)
        return -((n_k / mu - 1) @ (n_star * Rf)) * r

    x0 = np.full(free.sum(), np.log(max(n_k.sum(), 1) / n_star / max(Rf.sum(), 1e-9)))
    res = minimize(nll, x0, jac=grad, method="L-BFGS-B", bounds=[(-30, 5)] * free.sum())
    r = np.zeros(R.shape[1])
    r[free] = np.exp(res.x)
    return r


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--root", type=C.Path, default=C.STEP3)
    p.add_argument("--boot", type=int, default=100)
    a = p.parse_args()
    import step3 as S3
    from gedi import rates as R_

    z, meta, jobs, d = C.load_step3(a.root)
    n_star = meta["n_stars"]
    out = C.out_dir("rate_binned")
    nMf, nPf = len(S3.M_EDGES) - 1, len(S3.P_EDGES) - 1
    nMt, nPt = len(M_T) - 1, len(P_T) - 1
    rng = np.random.default_rng(5)
    results = {}
    L = ["=== OCCURRENCE RATE WITHOUT A POWER-LAW SHAPE (binned forward model) ===",
         f"step 3 folder: {a.root}; {n_star:,} stars; fitted cells {S3.M_EDGES.tolist()} M_J x "
         f"{S3.P_EDGES.tolist()} yr; true cells extend to {P_T[0]}-{P_T[-1]} yr.",
         f"A true cell is 'unconstrained' if fewer than {100 * MIN_RESP:g}% of its injected planets are reported "
         "anywhere; such cells are not extrapolated.",
         "Ranges: 16-84% bootstrap over detections (Poisson weights) and injections.", ""]
    for an in ("one", "best"):
        r = S3.analyse_one(z, meta, jobs, d, an)
        inj = r["inj"]
        Rm = R_.response_matrix(inj["m_all"], inj["P_all"], inj["t_m"], inj["t_P"], inj["f_m"], inj["f_P"],
                                M_T, P_T, S3.M_EDGES, S3.P_EDGES)
        resp = Rm.sum(axis=0)
        free = resp >= MIN_RESP
        results[an] = dict(R=Rm, free=free, scen={})
        L.append(f"--- {S3.LABEL[an]} ---")
        for scen in S3.SCENARIOS:
            h = r["scen"][scen]["hosts"]
            truth = R_.true_rate(h["pm"], h["pP"], n_star, M_T, P_T).ravel()
            est = {}
            for tag, mm, PP in (("A", h["A_m"], h["A_P"]), ("B", h["B_m"], h["B_P"])):
                nk = R_.counts_2d(mm, PP, S3.M_EDGES, S3.P_EDGES).ravel()
                r0 = fit_binned(nk, Rm, n_star, free)
                bs = []
                for _ in range(a.boot):
                    wi = rng.poisson(1.0, len(inj["m_all"])).astype(float)
                    Rb = R_.response_matrix(inj["m_all"], inj["P_all"], inj["t_m"], inj["t_P"], inj["f_m"],
                                            inj["f_P"], M_T, P_T, S3.M_EDGES, S3.P_EDGES, w_all=wi,
                                            w_det=wi[inj["idx"]])
                    wd = rng.poisson(1.0, len(mm)).astype(float)
                    nb = R_.counts_2d(mm, PP, S3.M_EDGES, S3.P_EDGES, weights=wd).ravel()
                    bs.append(fit_binned(nb, Rb, n_star, free))
                bs = np.array(bs)
                est[tag] = dict(r=r0, lo=np.percentile(bs, 16, axis=0), hi=np.percentile(bs, 84, axis=0), boot=bs)
            results[an]["scen"][scen] = dict(truth=truth, **est)
            L.append(f" [{scen}] inferred / true per true cell, A (real systems) | B (planets alone); "
                     "'-' = unconstrained")
            L.append("      M_J \\ yr  " + "".join(f"{P_T[j]:g}-{P_T[j + 1]:g}".center(19) for j in range(nPt)))
            for i in range(nMt):
                cells = []
                for j in range(nPt):
                    k = i * nPt + j
                    if not free[k] or truth[k] <= 0:
                        cells.append("-".center(19))
                        continue
                    qa, qb = est["A"]["r"][k] / truth[k], est["B"]["r"][k] / truth[k]
                    cells.append(f"{qa:.2f}|{qb:.2f}".center(19))
                L.append(f"    {M_T[i]:>4g}-{M_T[i + 1]:<4g}  " + "".join(cells))
            # well-constrained region: cells whose response is at least 5% (like step 3's eta >= 5%)
            good = (resp >= 0.05) & (truth > 0)
            tA, tB, tT = est["A"]["r"][good].sum(), est["B"]["r"][good].sum(), truth[good].sum()
            bA = est["A"]["boot"][:, good].sum(1) / tT
            bB = est["B"]["boot"][:, good].sum(1) / tT
            ratio = est["A"]["boot"][:, good].sum(1) / est["B"]["boot"][:, good].sum(1)
            L.append(f"   cells with response >= 5% ({int(good.sum())}): inferred / true A {tA / tT:.3f} "
                     f"[{np.percentile(bA, 16):.3f}, {np.percentile(bA, 84):.3f}], B {tB / tT:.3f} "
                     f"[{np.percentile(bB, 16):.3f}, {np.percentile(bB, 84):.3f}]; A / B {tA / tB:.3f} "
                     f"[{np.percentile(ratio, 16):.3f}, {np.percentile(ratio, 84):.3f}]")
            fA = r["scen"][scen]["fA"]
            fT = r["scen"][scen]["fT"]
            if fA and fT:
                L.append(f"   for comparison, power law over the whole window (step 3): inferred / true "
                         f"{fA.rate / fT.rate:.3f}")
            unc = (~free) & (truth > 0)
            L.append(f"   unconstrained cells hold {100 * truth[unc].sum() / truth.sum():.0f}% of the true planets "
                     "in the true-cell grid; no rate can be measured there without assuming a shape")
        L.append("")
    C.write_summary(out / "summary.txt", L)

    # ---------------------------------------------------------------- figures
    plt, ps = C.setup_plots()
    Rm, free = results["best"]["R"], results["best"]["free"]
    fig, ax = plt.subplots(1, 2, figsize=(15, 6.2), constrained_layout=True, gridspec_kw=dict(width_ratios=[1.3, 1]))
    show = np.vstack([Rm, 1 - Rm.sum(axis=0, keepdims=True)])
    im = ax[0].imshow(np.where(show > 0, show, np.nan), cmap="viridis", aspect="auto", vmin=0, vmax=1)
    tl = [f"{M_T[i]:g}-{M_T[i + 1]:g} M_J\n{P_T[j]:g}-{P_T[j + 1]:g} yr" for i in range(nMt) for j in range(nPt)]
    fl = [f"{S3.M_EDGES[i]:g}-{S3.M_EDGES[i + 1]:g} M_J, {S3.P_EDGES[j]:g}-{S3.P_EDGES[j + 1]:g} yr"
          for i in range(nMf) for j in range(nPf)] + ["not detected"]
    ax[0].set_xticks(range(len(tl)), tl, rotation=90, fontsize=6.5)
    ax[0].set_yticks(range(len(fl)), fl, fontsize=7)
    ax[0].set_xlabel("true cell of the injected planet")
    ax[0].set_ylabel("cell in which it is reported (fitted mass, period)")
    for j in np.nonzero(~free)[0]:
        ax[0].axvspan(j - 0.5, j + 0.5, color="0.85", alpha=0.6, lw=0)
    fig.colorbar(im, ax=ax[0], label="share of injected planets of the true cell")
    ax[0].set_title("Injection response (best model by BIC); grey columns: unconstrained", fontsize=10.5)
    cols = {"independent": ps.C_TARGET, "clustered": ps.C_WRONG}
    xs = np.arange(nMt * nPt)
    for k_, scen in enumerate(("independent", "clustered")):
        s = results["best"]["scen"][scen]
        ok = free & (s["truth"] > 0)
        for tag, mk, ls in (("A", "o", "-"), ("B", "s", "none")):
            q = s[tag]["r"] / np.where(s["truth"] > 0, s["truth"], np.nan)
            lo = s[tag]["lo"] / np.where(s["truth"] > 0, s["truth"], np.nan)
            hi = s[tag]["hi"] / np.where(s["truth"] > 0, s["truth"], np.nan)
            off = (k_ - 0.5) * 0.3 + (0.08 if tag == "B" else -0.08)
            ax[1].errorbar(xs[ok] + off, q[ok], [q[ok] - lo[ok], hi[ok] - q[ok]], fmt=mk, color=cols[scen],
                           mfc=cols[scen] if tag == "A" else "white", capsize=2, ms=4,
                           label=f"{scen}: " + ("real systems (A)" if tag == "A" else "planets alone (B)"))
    ax[1].axhline(1, color=ps.INK, lw=1)
    ax[1].set_xticks(xs[free], [tl[k] for k in xs[free]], rotation=90, fontsize=6.5)
    ax[1].set_yscale("log")
    ax[1].set_ylim(0.1, 10)
    ax[1].set_ylabel("inferred / true occurrence rate")
    ax[1].legend(fontsize=8)
    ax[1].set_title("Binned rate (one free rate per cell), best model by BIC", fontsize=10.5)
    fig.savefig(out / "rate_binned.png", dpi=150, bbox_inches="tight")
    print(f"saved {out / 'rate_binned.png'}")


if __name__ == "__main__":
    main()
