"""Step 4: which DR4 planets to distrust, and what radial-velocity follow-up would reveal.

Uses ONLY the step 3 results (results/step3_gaia/: population.npz, jobs.csv, meta.json,
fits.csv); no new Gaia fits, so it runs in about a minute:

    python scripts/step4.py                       # summary + figures in results/step4/
    python scripts/step4.py --root other/step3/folder --baseline 2

Part 1 - warning flags. For every planet that DR4 would report (best model by BIC, and
the one-planet fit for comparison) on the simulated real Gaia stars, it looks only at
quantities a DR4 user can see (adopted model, significance, period, eccentricity, signal
strength, period precision, RUWE) and measures how much each flag raises the chance that
the planet is not accurate (imprecise or spurious) or that the star hides another giant.

Part 2 - radial-velocity follow-up. For the same planets, using the true systems:
RV semi-amplitudes of the reported planet and of any hidden giant, the V magnitude, and
the number of HARPS or FEROS epochs (and telescope hours) needed to confirm the Gaia
orbit, to reject a spurious one, and to detect the hidden giant (as an orbit, or as a
trend over the programme). Instrument assumptions are in src/gaiabias/rv.py.

Bigger samples: everything scales with the step 3 sample (python scripts/step3.py --stars,
--sky-fraction, --inj-fits); this script only needs the resulting folder.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

ROOT3 = Path(__file__).resolve().parents[1] / "results" / "step3_gaia"
OUT = Path(__file__).resolve().parents[1] / "results" / "step4"
T_YR = 1977.0 / 365.25
ANALYSES = {"best": "best model by BIC (DR4-like)", "one": "one planet (L&W)"}
BIC_NAMES = ["no planet", "acceleration only", "jerk only", "one planet", "one planet + acceleration",
             "one planet + jerk", "two planets"]
FLAGS = {  # name -> (label, function of the planet table returning a boolean array)
    "model": ("adopted model is not a plain orbit\n(acceleration, jerk or two planets)", lambda t: t["bic"] != 3),
    "sig": ("significance < 5\n(DR4: OrbitalPoorlyConstrained)", lambda t: t["sig"] < 5),
    "longP": ("period > 0.7 x time span\n(> 3.8 yr)", lambda t: t["P"] > 0.7 * T_YR),
    "ecc": ("eccentricity > 0.5", lambda t: t["e1"] > 0.5),
    "weak": ("near threshold\n(Δχ² < 100)", lambda t: t["dchi2"] < 100),
    "sigP": ("period error > 10%", lambda t: t["sigP"] > 0.1),
    "ruwe": ("RUWE > 1.2", lambda t: t["ruwe"] > 1.2),
}
SUSPECT = ("model", "sig", "longP", "ecc")          # the combined rule recommended for follow-up
V_REACH = {"HARPS": 12.0, "FEROS": 10.0}            # practical faint limits for precise RV with each instrument


def load(root):
    import step3 as S3
    P = S3.paths(root)
    with np.load(P["pop"]) as f:
        z = {k: f[k] for k in f.files}
    meta = json.loads(P["meta"].read_text())
    with open(P["jobs"]) as f:
        jobs = list(csv.DictReader(f))
    return z, meta, jobs, S3.load_fits(root)


def planet_table(z, meta, jobs, d, analysis, scen, baseline_yr):
    """One row per planet reported in A (real systems) for one analysis and scenario."""
    import step3 as S3
    from gaiabias import rv as RV
    from gaiabias.step2.outcomes import bic_choice, reported

    pos = {j: k for k, j in enumerate(d["job"])}
    rep = reported(d, analysis, "w")
    bic = bic_choice(d, "w")
    pm, pP, pst = z[f"{scen}_m"], z[f"{scen}_P"] / 365.25, z[f"{scen}_star"]
    pe, pinc = z[f"{scen}_e"], z[f"{scen}_inc"]
    order = np.argsort(pst, kind="stable")
    starts = np.searchsorted(pst[order], np.arange(meta["n_stars"] + 1))
    G, plx, mst = z["star_G"], z["star_parallax"], z["star_m_star"]
    get = lambda name, k: float(d[name][k]) if name in d else np.nan  # noqa: E731
    rows = []
    for j in jobs:
        if j["kind"] != "A" or j["scen"] != scen or j["job"] not in pos:
            continue
        k, s = pos[j["job"]], int(j["star"])
        idx = order[starts[s]:starts[s + 1]]
        reps = [(P_[k] / 365.25, m_[k], sg[k]) for mask, P_, m_, sg in rep if mask[k]]
        if not reps:
            continue
        matched = set()
        for Py, m, sg in reps:
            cls = S3.classify(Py, m, pP[idx], pm[idx])
            near = idx[np.argmin(np.abs(np.log(Py / pP[idx])))]
            if cls < 2:
                matched.add(int(near))
            rows.append(dict(star=s, k=k, P=Py, m=m, sig=sg, cls=cls, near=int(near)))
        for r in rows[len(rows) - len(reps):]:
            r["hidden"] = [int(i) for i in idx if int(i) not in matched]
            r.update(bic=int(bic[k]), e1=get("one_w_e", k), dchi2=get("one_w_dchi2", k), ruwe=get("ruwe_dr4", k),
                     sigP=get("one_w_sigP", k) / max(get("one_w_P", k), 1e-9), n_tr=get("n_transits", k))
    t = {key: np.array([r[key] for r in rows]) for key in ("star", "P", "m", "sig", "cls", "near", "bic", "e1",
                                                           "dchi2", "ruwe", "sigP", "n_tr")}
    t["star"], t["near"], t["cls"], t["bic"] = (t[c].astype(int) for c in ("star", "near", "cls", "bic"))
    hidden = [r["hidden"] for r in rows]
    t["n_hidden"] = np.array([len(h) for h in hidden])
    # radial velocities
    ms = mst[t["star"]]
    M_G = G[t["star"]] + 5 * np.log10(plx[t["star"]] / 100.0)
    t["V"] = RV.v_mag(G[t["star"]], RV.bp_rp_from_MG(M_G))
    t["m_star"] = ms
    sin_near = np.sin(pinc[t["near"]])
    t["K_pred"] = RV.semi_amplitude(t["m"], t["P"], ms, 0.0, sin_near)          # what the Gaia orbit predicts
    t["K_near"] = RV.semi_amplitude(pm[t["near"]], pP[t["near"]], ms, pe[t["near"]], sin_near)
    for inst in RV.INSTRUMENTS:
        sig = RV.rv_precision(t["V"], ms, inst)
        t[f"sig_{inst}"] = sig
        t[f"N_conf_{inst}"] = RV.epochs_known_period(t["K_pred"], sig)
        dK = np.sqrt((t["K_pred"] ** 2 + t["K_near"] ** 2) / 2)           # spurious: signals do not match
        t[f"N_rej_{inst}"] = RV.epochs_known_period(dK, sig)
        n_c = np.full(len(rows), np.inf)
        for i, h in enumerate(hidden):
            if h:
                h = np.array(h)
                Kc = RV.semi_amplitude(pm[h], pP[h], ms[i], pe[h], np.sin(pinc[h]))
                n_c[i] = float(np.min(RV.epochs_companion(Kc, pP[h], sig[i], baseline_yr)))
        t[f"N_comp_{inst}"] = n_c
    t["K_hidden_max"] = np.array([float(np.max(RV.semi_amplitude(pm[h], pP[h], ms[i], pe[h], np.sin(pinc[h]))))
                                  if h else 0.0 for i, h in enumerate(hidden)])
    t["P_hidden"] = np.array([float(pP[h][np.argmax(pm[h] * pP[h] ** (-1 / 3))]) if h else np.nan
                              for h in hidden])
    return t


def flag_stats(t):
    bad = t["cls"] >= 1
    spur = t["cls"] == 2
    hid = t["n_hidden"] > 0
    out = {}
    masks = {k: f(t) for k, (_, f) in FLAGS.items()}
    masks["suspect"] = np.any([masks[k] for k in SUSPECT], axis=0)
    for k, m in masks.items():
        m = np.nan_to_num(m, nan=0).astype(bool)
        frac = lambda x, sel: float(x[sel].mean()) if sel.any() else np.nan  # noqa: E731
        out[k] = dict(n=int(m.sum()), share=float(m.mean()), bad_in=frac(bad, m), bad_out=frac(bad, ~m),
                      spur_in=frac(spur, m), spur_out=frac(spur, ~m), hid_in=frac(hid, m), hid_out=frac(hid, ~m),
                      caught_bad=float(m[bad].mean()) if bad.any() else np.nan,
                      caught_spur=float(m[spur].mean()) if spur.any() else np.nan,
                      caught_hid=float(m[hid].mean()) if hid.any() else np.nan)
    out["_base"] = dict(n=len(bad), bad=float(bad.mean()), spur=float(spur.mean()), hid=float(hid.mean()))
    out["_masks"] = masks
    return out


def summary(tabs, stats, meta, baseline_yr, scale):
    import step3 as S3
    from gaiabias import rv as RV
    L = ["=== STEP 4 SUMMARY ===",
         f"Source: step 3 fits ({meta['n_stars']:,} real Gaia stars <= 250 pc). Full-sky scaling x{scale:.2f} "
         f"({meta['n_selected']:,} selected stars / sky fraction {S3.sky_fraction(meta):.2f} / stars used).",
         "Accuracy of a reported planet: correct / imprecise / spurious (step 3 definitions); 'not accurate' = "
         "imprecise or spurious. 'Hidden giant' = another giant planet on the star that was not reported.",
         f"RV programme baseline for companions: {baseline_yr:g} yr; detection = 5 sigma (src/gaiabias/rv.py).", ""]
    for an, lab in ANALYSES.items():
        for scen in S3.SCENARIOS:
            t, st = tabs[(an, scen)], stats[(an, scen)]
            b = st["_base"]
            L.append(f"--- {lab}, {scen} multiplicity: {b['n']} reported planets ({b['n'] * scale:,.0f} full sky); "
                     f"not accurate {100 * b['bad']:.0f}%, spurious {100 * b['spur']:.0f}%, "
                     f"star hides another giant {100 * b['hid']:.0f}% ---")
            L.append("   flag                                    share  | not accurate: flagged vs not | spurious: "
                     "flagged vs not | hidden giant: flagged vs not | share of all not-accurate / spurious / hidden caught")
            for k in list(FLAGS) + ["suspect"]:
                s = st[k]
                name = (FLAGS[k][0].replace("\n", " ") if k in FLAGS else
                        "SUSPECT = any of: model, significance, period, eccentricity")
                L.append(f"   {name[:40]:<40} {100 * s['share']:4.0f}%  | {100 * s['bad_in']:3.0f}% vs "
                         f"{100 * s['bad_out']:3.0f}%            | {100 * s['spur_in']:3.0f}% vs {100 * s['spur_out']:3.0f}%"
                         f"          | {100 * s['hid_in']:3.0f}% vs {100 * s['hid_out']:3.0f}%            | "
                         f"{100 * s['caught_bad']:3.0f}% / {100 * s['caught_spur']:3.0f}% / {100 * s['caught_hid']:3.0f}%")
            sus = st["_masks"]["suspect"]
            for inst in RV.INSTRUMENTS:
                nconf = t[f"N_conf_{inst}"]
                ncomp = t[f"N_comp_{inst}"]
                hid = t["n_hidden"] > 0
                spur = t["cls"] == 2
                L.append(f"   {inst}: per-epoch precision median {np.median(t[f'sig_{inst}']):.1f} m/s (V median "
                         f"{np.median(t['V']):.1f}); confirm the Gaia orbit with <= 10 epochs: "
                         f"{100 * np.mean(nconf <= 10):.0f}% of planets (median {np.median(nconf):.0f} epochs); "
                         f"reject a spurious orbit with <= 10 epochs: "
                         f"{100 * np.mean(t[f'N_rej_{inst}'][spur] <= 10) if spur.any() else np.nan:.0f}%")
                if hid.any():
                    L.append(f"      hidden giants detectable with <= 10 / 20 / 40 epochs over {baseline_yr:g} yr: "
                             f"{100 * np.mean(ncomp[hid] <= 10):.0f}% / {100 * np.mean(ncomp[hid] <= 20):.0f}% / "
                             f"{100 * np.mean(ncomp[hid] <= 40):.0f}%  (among SUSPECT planets: "
                             f"{100 * np.mean(ncomp[hid & sus] <= 20) if (hid & sus).any() else np.nan:.0f}% with <= 20)")
                n_s = int(sus.sum())
                cost = RV.hours(np.minimum(np.where(np.isfinite(ncomp), np.maximum(ncomp, nconf), nconf), 40), inst)
                L.append(f"      follow-up of all SUSPECT planets (<= 40 epochs each): {n_s} here, "
                         f"{n_s * scale:,.0f} full sky, {cost[sus].sum() * scale:,.0f} h full sky "
                         f"({np.median(cost[sus]) if n_s else np.nan:.1f} h per star median)")
                reach = sus & (t["V"] <= V_REACH[inst])
                L.append(f"      ... of these with V <= {V_REACH[inst]:g} (practical {inst} limit): "
                         f"{reach.sum() * scale:,.0f} full sky ({100 * reach.sum() / max(n_s, 1):.0f}% of SUSPECT), "
                         f"{cost[reach].sum() * scale:,.0f} h; their hidden giants found with <= 20 epochs: "
                         f"{100 * np.mean(ncomp[reach & hid] <= 20) if (reach & hid).any() else np.nan:.0f}%")
            L.append("")
    return L


def plot(tabs, stats, png_flags, png_rv, baseline_yr):
    import step3 as S3
    from gaiabias import rv as RV
    plt, ps = S3._setup()
    keys = list(FLAGS) + ["suspect"]
    labels = [FLAGS[k][0] if k in FLAGS else "SUSPECT\n(any of the first four)" for k in keys]
    fig, ax = plt.subplots(1, 3, figsize=(19, 6.2), constrained_layout=True)
    y = np.arange(len(keys))[::-1]
    for c, (what, title) in enumerate((("bad", "share of reported planets that are NOT accurate"),
                                       ("hid", "share of reported planets whose star hides\nanother giant planet (0.3-13 M_J, 0.1-30 au)"),
                                       ("caught_bad", "share of all inaccurate planets\nthat the flag catches"))):
        a = ax[c]
        for i, scen in enumerate(S3.SCENARIOS):
            st = stats[("best", scen)]
            off = (i - 0.5) * 0.36
            if what == "caught_bad":
                v = [st[k]["caught_bad"] for k in keys]
                a.barh(y + off, v, 0.34, color=ps.C_WRONG if i == 0 else "#f3b98c", label=f"{scen} multiplicity")
            else:
                vin = [st[k][f"{what}_in"] for k in keys]
                vout = [st[k][f"{what}_out"] for k in keys]
                a.barh(y + off, vin, 0.34, color=ps.C_WRONG if i == 0 else "#f3b98c",
                       label=f"flagged ({scen})")
                a.scatter(vout, y + off, marker="|", s=160, color=ps.INK, zorder=3,
                          label="not flagged" if i == 0 else None)
        a.set_yticks(y, labels if c == 0 else [""] * len(keys), fontsize=8.5)
        ps.percent_axis(a, "x", "")
        a.set_title(title, fontsize=10.5)
    h, l = ax[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, fontsize=9, bbox_to_anchor=(0.5, -0.07))
    fig.suptitle("Step 4: warning signs in what DR4 publishes (best model by BIC; real Gaia stars, simulated planets)",
                 fontsize=12)
    fig.savefig(png_flags, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {png_flags}")

    fig, ax = plt.subplots(1, 3, figsize=(19, 5.6), constrained_layout=True)
    n = np.arange(1, 81)
    sty = {"HARPS": "-", "FEROS": "--"}
    for i, scen in enumerate(S3.SCENARIOS):
        t = tabs[("best", scen)]
        col = ps.C_TARGET if i == 0 else "#8fb8f0"
        sus = stats[("best", scen)]["_masks"]["suspect"]
        hid = t["n_hidden"] > 0
        for inst in RV.INSTRUMENTS:
            ax[0].plot(n, [np.mean(t[f"N_conf_{inst}"] <= k) for k in n], sty[inst], color=col,
                       label=f"{inst}, {scen}")
            if hid.any():
                ax[1].plot(n, [np.mean(t[f"N_comp_{inst}"][hid] <= k) for k in n], sty[inst], color=col,
                           label=f"{inst}, {scen}")
        ax[2].scatter(t["V"][hid], t["K_hidden_max"][hid], s=6, alpha=0.35, color=col, label=f"hidden giants ({scen})")
        ax[2].scatter(t["V"][sus & hid], t["K_hidden_max"][sus & hid], s=6, alpha=0.6, color=ps.C_WRONG,
                      label="... of SUSPECT planets" if i == 0 else None)
    ax[0].set_title("Confirming the Gaia planet with RV\n(period known from Gaia)", fontsize=10.5)
    ax[1].set_title(f"Finding the hidden giant with RV\n(orbit or trend over a {baseline_yr:g}-yr programme)",
                    fontsize=10.5)
    for a in ax[:2]:
        a.set_xlabel("number of RV epochs")
        ps.percent_axis(a, "y", "share of stars")
        a.set_ylim(0, 1)
        a.legend(fontsize=8, loc="lower right")
    ax[2].set_yscale("log")
    ax[2].set_xlabel("V magnitude")
    ax[2].set_ylabel("RV semi-amplitude of the hidden giant (m/s)")
    for inst in RV.INSTRUMENTS:
        V = np.linspace(4, 16, 50)
        ax[2].plot(V, 5 * RV.rv_precision(V, 1.0, inst) / np.sqrt(10), sty[inst], color=ps.INK, lw=1,
                   label=f"{inst}: 5σ with 10 epochs")
    for inst, vv in V_REACH.items():
        ax[2].axvline(vv, color="0.6", lw=1, ls=":")
        ax[2].text(vv - 0.15, 1.2, f"{inst} limit", rotation=90, fontsize=8, color="0.4", ha="right")
    ax[2].set_title("Signal size of the hidden giants", fontsize=10.5)
    ax[2].legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=3)
    fig.suptitle("Step 4: radial-velocity follow-up of DR4 planets (best model by BIC)", fontsize=12)
    fig.savefig(png_rv, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {png_rv}")


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    p.add_argument("--root", type=Path, default=ROOT3, help="step 3 results folder")
    p.add_argument("--out", type=Path, default=OUT)
    p.add_argument("--baseline", type=float, default=3.0, help="RV programme length in years")
    a = p.parse_args()
    import step3 as S3
    z, meta, jobs, d = load(a.root)
    scale = meta["n_selected"] / S3.sky_fraction(meta) / meta["n_stars"]
    tabs, stats = {}, {}
    for an in ANALYSES:
        for scen in S3.SCENARIOS:
            tabs[(an, scen)] = planet_table(z, meta, jobs, d, an, scen, a.baseline)
            stats[(an, scen)] = flag_stats(tabs[(an, scen)])
    a.out.mkdir(parents=True, exist_ok=True)
    L = summary(tabs, stats, meta, a.baseline, scale)
    print("\n".join(L))
    (a.out / "summary.txt").write_text("\n".join(L) + "\n")
    with open(a.out / "planets.csv", "w", newline="") as f:          # one row per reported planet (best, both)
        w = csv.writer(f)
        cols = ["star", "P", "m", "sig", "cls", "bic", "e1", "dchi2", "ruwe", "sigP", "n_hidden", "V", "K_pred",
                "K_hidden_max", "P_hidden", "N_conf_HARPS", "N_comp_HARPS", "N_conf_FEROS", "N_comp_FEROS"]
        w.writerow(["scenario"] + cols)
        for scen in S3.SCENARIOS:
            t = tabs[("best", scen)]
            for i in range(len(t["P"])):
                w.writerow([scen] + [f"{t[c][i]:.4g}" if np.ndim(t[c][i]) == 0 else "" for c in cols])
    plot(tabs, stats, a.out / "step4_flags.png", a.out / "step4_rv.png", a.baseline)


if __name__ == "__main__":
    main()
