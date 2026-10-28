"""Why does the power-law rate of planets fitted ALONE (B) drift upward with multiplicity?

Uses only the step 3 results (no new fits):

    python scripts/extended/g1_pl_drift.py [--root results/step3_allsky]

In the whole-sky step 3 run, the forward-modelled power-law rate over the whole window
(0.3-13 M_J, 0.1-8 yr) for each planet fitted alone (B) rises from 1.01 to about 1.6 of the
true rate as the fraction f of giant hosts with a second giant goes from 0 to 0.8, while the
binned rate in the well-measured cells stays near 1. B has no multiplicity in its fits, so the
drift must come from the population or from the method. For each f of the multiplicity scan
(same host reweighting as step 3, nominal weights) this script compares five power-law rates,
each divided by the true (CLS) in-window rate:

  B          fit to the B detections (as in step 3);
  oracle     fit to the counts the injection response PREDICTS for the true planets
             (expected fitted counts N_star * R @ true counts): what B should give if the
             response measured from single injections is right;
  truth      fit directly to the true planets in the window (no detection at all);
  measured   fit to the true planets, but only in the fine cells inside the well-measured
             coarse cells (completeness >= 5%), extrapolated to the whole window: the pure
             effect of the population shape on the extrapolation;
  B/oracle   ratio of B detections to predicted counts in the well-measured region.

Reading: if `measured` follows B, the drift is the shape of the population (the power law
fitted where Gaia measures does not describe the rest of the window); if `oracle` follows B
but `measured` does not, the response (migration from outside the window) is the cause; if B
exceeds `oracle`, the detections are not what the single-planet injections predict.
Also reported: the share of true in-window planets that lie in the well-measured region, and
the slopes (alpha, beta) of the `truth` and `measured` fits.
Output: results/extended/pl_drift/ (summary.txt, pl_drift.png, pl_drift.csv).
"""
from __future__ import annotations

import argparse
import csv

import numpy as np

import common as C

SNR_EDGES = np.array([0.8, 1.2, 2.0, 3.5, 6.0, 12.0, 30.0, 1e6])


def _wmedian(x, w):
    o = np.argsort(x)
    c = np.cumsum(w[o])
    return float(x[o][np.searchsorted(c, 0.5 * c[-1])])


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--root", default=str(C.REPO / "results" / "step3_allsky"))
    a = p.parse_args()
    import step3 as S3
    from gaiabias import multiplicity as MU
    from gaiabias import rates as R
    from pathlib import Path

    z, meta, jobs, d = C.load_step3(Path(a.root))
    out = C.out_dir("pl_drift")
    n_star = meta["n_stars"]
    m_win, P_win = (S3.M_FINE[0], S3.M_FINE[-1]), (S3.P_FINE[0], S3.P_FINE[-1])
    nmf, npf = len(S3.M_FINE) - 1, len(S3.P_FINE) - 1
    nPt = len(S3.P_TRUE_EDGES) - 1
    fm_c = np.sqrt(S3.M_FINE[:-1] * S3.M_FINE[1:])
    fp_c = np.sqrt(S3.P_FINE[:-1] * S3.P_FINE[1:])
    cm, cp = np.meshgrid(fm_c, fp_c, indexing="ij")              # fitted-cell centres (pseudo-points)
    cm, cp = cm.ravel(), cp.ravel()

    def fit_counts(counts, Rm, true_edges_P):
        """Forward power law fitted to expected counts in the fine fitted cells (no Poisson noise)."""
        f_ = R.fit_power_law_forward(cm, cp, n_star, Rm, S3.M_FINE, true_edges_P, S3.M_FINE, S3.P_FINE,
                                     m_win=m_win, P_win=P_win, weights=counts, errors=False)
        return f_

    rows, L = [], ["=== WHY THE POWER-LAW RATE OF PLANETS ALONE DRIFTS WITH MULTIPLICITY ===",
                   f"step 3 results: {a.root} ({n_star:,} stars). Rates are divided by the true (CLS) in-window rate.", ""]
    for an in S3.SCAN_METHODS:
        print(f"analysis {an} ...", flush=True)
        r = S3.analyse_one(z, meta, jobs, d, an)
        hosts = {sc: r["scen"][sc]["hosts"] for sc in S3.SCENARIOS}
        pool = MU.build_pool(hosts)
        # true period [yr] of the planet behind each B detection, in the order build_pool uses
        b_trueP = np.concatenate([hosts[sc]["pP"][hosts[sc]["B_pl"]] for sc in hosts])
        long_ = b_trueP > S3.INJ_P_YR[1]                      # beyond the longest injected period
        t_snr = np.concatenate([z[f"{sc}_snr"] for sc in hosts])     # per-transit S/N of every true planet
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
        r_cls = float(win0.sum()) / n_star
        t_in = MU.in_window(pool.t_m, pool.t_P, m_win, P_win)
        t_win = np.bincount(pool.t_host, weights=t_in, minlength=n_hosts)
        region = r["eta"] >= S3.MIN_ETA
        # fine fitted cells inside the well-measured coarse cells
        im_c = np.clip(np.digitize(cm, S3.M_EDGES) - 1, 0, len(S3.M_EDGES) - 2)
        ip_c = np.clip(np.digitize(cp, S3.P_EDGES) - 1, 0, len(S3.P_EDGES) - 2)
        meas_fine = region[im_c, ip_c]
        inj = r["inj"]
        Rm = R.response_matrix(inj["m_all"], inj["P_all"], inj["t_m"], inj["t_P"], inj["f_m"], inj["f_P"],
                               S3.M_FINE, S3.P_TRUE_EDGES, S3.M_FINE, S3.P_FINE)
        # response conditioned on S/N as well: R_snr[k, j, b] = reported in fitted cell k per injection in
        # true cell j and S/N bin b (S/N below the pre-fit cut is never fitted, so never detected)
        inj_tab = z["inj"]
        row_of = {int(k): i for i, k in enumerate(z["inj_idx"])}
        nT = (len(S3.M_FINE) - 1) * nPt
        def tcell(m, P):
            jm, jp = np.digitize(m, S3.M_FINE) - 1, np.digitize(P, S3.P_TRUE_EDGES) - 1
            ok = (jm >= 0) & (jm < len(S3.M_FINE) - 1) & (jp >= 0) & (jp < nPt)
            return np.where(ok, jm * nPt + jp, -1)
        def fcell(m, P):
            km, kp = np.digitize(m, S3.M_FINE) - 1, np.digitize(P, S3.P_FINE) - 1
            ok = (km >= 0) & (km < nmf) & (kp >= 0) & (kp < npf)
            return np.where(ok, km * npf + kp, -1)
        snr_b = np.digitize(inj_tab[:, 8], SNR_EDGES) - 1
        j_inj = tcell(inj_tab[:, 1], inj_tab[:, 2] / 365.25)
        den = np.zeros((nT, len(SNR_EDGES) - 1))
        okd = (j_inj >= 0) & (snr_b >= 0) & (snr_b < len(SNR_EDGES) - 1)
        np.add.at(den, (j_inj[okd], snr_b[okd]), 1.0)
        rows_det = np.array([row_of[int(k)] for k in inj["idx"]], int)
        k_det = fcell(inj["f_m"], inj["f_P"])
        okn = (k_det >= 0) & okd[rows_det]
        Rs = np.zeros((nmf * npf, nT, len(SNR_EDGES) - 1))
        np.add.at(Rs, (k_det[okn], j_inj[rows_det[okn]], snr_b[rows_det[okn]]), 1.0)
        with np.errstate(invalid="ignore", divide="ignore"):
            Rs = np.where(den > 0, Rs / den, 0.0)
        j_t = tcell(pool.t_m, pool.t_P)
        b_t = np.digitize(t_snr, SNR_EDGES) - 1
        okt = (j_t >= 0) & (b_t >= 0) & (b_t < len(SNR_EDGES) - 1)
        I_fine = np.eye(nmf * npf)                                  # truth: identity response on the window
        M_fine = np.diag(meas_fine.astype(float))                   # truth seen only where Gaia measures
        L += [f"--- {S3.LABEL[an]} ---",
              "     f    B   B(P<16yr) oracle oracle(S/N) truth measured  B/oracle meas./other  B/oracle(S/N) meas./other"
              "  B from P>16yr  median S/N true low cells"]
        for f in S3.SCAN_F:
            v0, mu, lost, f_real = MU.host_weights_by_mass(pool, f, host_bin)
            ew = np.bincount(host_bin, weights=v0 * t_win, minlength=S3.N_MBIN)
            with np.errstate(invalid="ignore", divide="ignore"):
                Hj = np.where(ew > 0, lam_bin / ew, 0.0) * frac_bin
            v = v0 * Hj[host_bin]
            wb, wt = v[pool.b_host] * n_star, v[pool.t_host] * n_star
            res = {}
            try:
                res["B"] = R.fit_power_law_forward(pool.b_m, pool.b_P, n_star, Rm, S3.M_FINE, S3.P_TRUE_EDGES,
                                                   S3.M_FINE, S3.P_FINE, weights=wb, errors=False)
            except Exception:
                res["B"] = None
            try:
                res["B_inj"] = R.fit_power_law_forward(pool.b_m[~long_], pool.b_P[~long_], n_star, Rm, S3.M_FINE,
                                                       S3.P_TRUE_EDGES, S3.M_FINE, S3.P_FINE,
                                                       weights=wb[~long_], errors=False)
            except Exception:
                res["B_inj"] = None
            c_true = R.counts_2d(pool.t_m, pool.t_P, S3.M_FINE, S3.P_TRUE_EDGES, weights=wt).ravel()
            n_or = Rm @ c_true                                          # predicted fitted counts
            c_win = R.counts_2d(pool.t_m, pool.t_P, S3.M_FINE, S3.P_FINE, weights=wt).ravel()
            res["oracle"] = fit_counts(n_or, Rm, S3.P_TRUE_EDGES)
            c_ts = np.zeros((nT, len(SNR_EDGES) - 1))
            np.add.at(c_ts, (j_t[okt], b_t[okt]), wt[okt])
            n_or2 = np.einsum("kjb,jb->k", Rs, c_ts)                     # S/N-aware prediction
            res["oracle_snr"] = fit_counts(n_or2, Rm, S3.P_TRUE_EDGES)
            res["truth"] = fit_counts(c_win, I_fine, S3.P_FINE)
            res["measured"] = fit_counts(c_win * meas_fine, M_fine, S3.P_FINE)
            nB = R.counts_2d(pool.b_m, pool.b_P, S3.M_FINE, S3.P_FINE, weights=wb).ravel()
            with np.errstate(invalid="ignore", divide="ignore"):
                b_over_or = float(nB[meas_fine].sum() / n_or[meas_fine].sum())
            share = float(c_win[meas_fine].sum() / c_win.sum())
            nB_long = R.counts_2d(pool.b_m[long_], pool.b_P[long_], S3.M_FINE, S3.P_FINE, weights=wb[long_]).ravel()
            with np.errstate(invalid="ignore", divide="ignore"):
                b_over_or_un = float(nB[~meas_fine].sum() / n_or[~meas_fine].sum())
                long_share_un = float(nB_long[~meas_fine].sum() / nB[~meas_fine].sum())
                long_share = float(nB_long.sum() / nB.sum())
                b_over_or2_un = float(nB[~meas_fine].sum() / n_or2[~meas_fine].sum())
                b_over_or2_m = float(nB[meas_fine].sum() / n_or2[meas_fine].sum())
            # median S/N of true planets vs injections, in true cells with low completeness (fitted-cell centre
            # not in the well-measured region), counting only planets above the pre-fit cut
            jp_t = j_t % nPt
            in_win = (jp_t >= 1) & (jp_t <= npf)                  # true P cells 1..npf are the window's
            fine_t = np.clip(j_t // nPt * npf + jp_t - 1, 0, nmf * npf - 1)
            low_t = ~in_win | ~meas_fine[fine_t]
            sel_t = okt & low_t & (t_snr >= S3.S_CUT)
            snr_true_low = _wmedian(t_snr[sel_t], wt[sel_t]) if sel_t.any() else np.nan
            q = {k: (x.rate / r_cls if x is not None else np.nan) for k, x in res.items()}
            row = dict(method=an, f=f, f_real=f_real, **{f"rate_{k}": v_ for k, v_ in q.items()},
                       B_over_oracle_meas=b_over_or, B_over_oracle_unmeas=b_over_or_un,
                       B_from_long_P=long_share, B_from_long_P_unmeas=long_share_un, share_true_in_measured=share,
                       B_over_oracle_snr_meas=b_over_or2_m,
                       B_over_oracle_snr_unmeas=b_over_or2_un, median_snr_true_lowcells=snr_true_low,
                       truth_alpha=res["truth"].alpha, truth_beta=res["truth"].beta,
                       meas_alpha=res["measured"].alpha, meas_beta=res["measured"].beta)
            rows.append(row)
            L.append(f"  {f:4.2f}  {q['B']:5.2f}   {q['B_inj']:5.2f}   {q['oracle']:5.2f}   {q['oracle_snr']:5.2f}    "
                     f"{q['truth']:5.2f}  {q['measured']:5.2f}      {b_over_or:5.2f} / {b_over_or_un:5.2f}        "
                     f"{b_over_or2_m:5.2f} / {b_over_or2_un:5.2f}         {100 * long_share:4.1f}%          "
                     f"{snr_true_low:5.2f}")
        L.append("")
    with open(out / "pl_drift.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    L += [f"oracle(S/N): the same prediction with the response measured separately in S/N bins {list(SNR_EDGES)}",
          "(per-transit S/N of each planet on its own star). If B/oracle(S/N) is ~1 everywhere, the drift comes from",
          "the planets of multi-giant hosts having different S/N within a cell than the injections, i.e. the",
          "completeness of a cell depends on which stars host its planets.",
          "B(P<16yr): B fitted without the detections of planets whose true period exceeds the longest injected",
          "period (16 yr), which the injection response cannot describe. 'other' = fitted cells outside the",
          "well-measured region.",
          "Reading: 'measured' following B = the population shape makes the power-law extrapolation wrong;",
          "'oracle' following B but not 'measured' = the response (migration into the window) is the cause;",
          "B/oracle far from 1 = the detections differ from what single-planet injections predict."]
    C.write_summary(out / "summary.txt", L)

    plt, ps = C.setup_plots()
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True, sharey=True)
    for ax, an in zip(axs, S3.SCAN_METHODS):
        sel = [rw for rw in rows if rw["method"] == an]
        f = np.array([rw["f"] for rw in sel]) * 100
        for k, st, lab in (("B", dict(lw=2.2), "planets fitted alone (B)"),
                           ("B_inj", dict(lw=1.2, alpha=0.7), "B without planets beyond 16 yr"),
                           ("oracle", dict(ls="--"), "predicted from the injections"),
                           ("oracle_snr", dict(ls="--", lw=2.0), "predicted, response per S/N bin"),
                           ("measured", dict(ls="-."), "truth fitted in well-measured cells only"),
                           ("truth", dict(ls=":"), "truth")):
            ax.plot(f, [rw[f"rate_{k}"] for rw in sel], marker="o", ms=3.5, label=lab, **st)
        ax.axhline(1, color="0.3", lw=0.8)
        ax.set_title(S3.LABEL[an])
        ax.set_xlabel("giant-planet hosts with a second giant (%)")
    axs[0].set_ylabel("power-law rate in the window / true rate")
    axs[0].legend(frameon=False, fontsize=8)
    fig.savefig(out / "pl_drift.png", dpi=150)


if __name__ == "__main__":
    main()
