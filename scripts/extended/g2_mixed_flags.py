"""Warning flags with brown-dwarf and stellar-binary contaminants.

    python scripts/extended/g2_mixed_flags.py              # simulate and fit (resumable), then analyse
    python scripts/extended/g2_mixed_flags.py --analyse    # analysis only
    python scripts/extended/g2_mixed_flags.py --quick      # 300 + 600 systems, a short test

The step 4 warning flags were tested on planet-only simulations (g1_flags_heldout.py). Real DR4
orbits with planetary masses also come from brown dwarfs and from stellar binaries whose
photocentre moves much less than the primary, because the companion's light pulls the
photocentre towards it. This script adds such systems to the step 3 real-star simulation:

* host stars: drawn at random from the step 3 sample (its real positions, distances, G and masses);
* brown dwarfs: 13-80 M_J (log-uniform), Kipping (2013) eccentricities;
* stellar companions: mass ratio q uniform in [0.1, 1] (Raghavan et al. 2010) with the companion
  above 0.08 M_sun, eccentricity uniform in [0, 0.8];
* periods log-uniform over 0.05-16 yr, the window of the step 3 injections, isotropic orientations;
* photocentre: the reflex orbit of the primary scaled by (q/(1+q) - beta)/(q/(1+q)), with
  beta = L2/(L1+L2) in the G band from the main-sequence mass-M_G relation used for the step 3
  masses. It is simulated as a dark companion with the same angular photocentre orbit, which is
  how DR4 would interpret it (the mass DR4 infers is the "effective" dark mass);
* the same DR4 fits as step 3 (one orbit, + acceleration, + jerk, two orbits, BIC choice), and the
  same selection: the DR3 RUWE of the simulated system must be below the step 3 cut (1.4), since
  the real stars were selected that way. Systems with a summed signal-to-noise below 0.8 are not
  fitted (as in step 3: they would not be detected).

Rates per star (assumptions, varied in the summary):
* brown dwarf with P = 0.05-16 yr: 0.6% (the brown-dwarf desert, e.g. Grether & Lineweaver 2006);
* stellar companion: binary fraction 0.46 for > 0.6 M_sun and 0.27 below (Raghavan et al. 2010;
  Winters et al. 2019) times the share of the log-normal period distribution (mean log10 P/d = 5.03,
  sigma 2.28) that falls in 0.05-16 yr.

Every planet-mass orbit (m < 13 M_J, passing the detection test) reported for a brown dwarf or
a binary counts as "not a real planet". The flags are then scored on held-out stars for the
planet-only and the mixed population (contaminants weighted by their rates), as in
g1_flags_heldout.py. Output: results/extended/mixed_flags/ (systems.csv, summary.txt,
mixed_flags.png).
"""
from __future__ import annotations

import argparse
import csv
import itertools

import numpy as np

import common as C

P_WIN_YR = (0.05, 16.0)
BD_RATE = 0.006
BIN_FRAC = (0.46, 0.27)          # m_star > 0.6 / <= 0.6 M_sun
LOGP_MU, LOGP_SIG = 5.03, 2.28   # log10 P[d], Raghavan et al. (2010)
RATE_SCALES = (0.5, 1.0, 2.0)
FITTED = ["one", "accel", "jerk", "two"]
SEED = 2031


def window_share():
    from scipy.stats import norm
    lo, hi = np.log10(np.array(P_WIN_YR) * 365.25)
    return norm.cdf((hi - LOGP_MU) / LOGP_SIG) - norm.cdf((lo - LOGP_MU) / LOGP_SIG)


def mg_of_mass(m):
    """Main-sequence M_G for a mass in M_sun (the table behind the step 3 masses, extrapolated)."""
    from gaiabias.population import _MS_MASS, _MS_MG
    m = np.atleast_1d(np.asarray(m, float))
    out = np.interp(m, _MS_MASS, _MS_MG)
    slope = (_MS_MG[1] - _MS_MG[0]) / (_MS_MASS[1] - _MS_MASS[0])
    low = m < _MS_MASS[0]
    out[low] = _MS_MG[0] + slope * (m[low] - _MS_MASS[0])
    return out


def photocentre_factor(q, m1):
    """(q/(1+q) - beta) / (q/(1+q)): photocentre orbit relative to the primary's reflex orbit."""
    q = np.asarray(q, float)
    m2 = q * m1
    lum = m2 >= 0.08
    beta = np.zeros_like(q)
    r = 10 ** (-0.4 * (mg_of_mass(m2[lum]) - mg_of_mass(np.full(lum.sum(), m1))))
    beta[lum] = r / (1 + r)
    f = q / (1 + q)
    return (f - beta) / f


def effective_mass_mj(m1, m2_mj, P_days, factor):
    """Dark-companion mass with the same angular reflex orbit as the photocentre."""
    from scipy.optimize import brentq
    from gaiabias.orbits import reflex_a0_mas
    target = abs(factor) * reflex_a0_mas(m1, m2_mj, P_days, 1.0)
    lo, hi = 1e-6, max(m2_mj, 1e-3) * 1.0001
    if target <= reflex_a0_mas(m1, lo, P_days, 1.0):
        return 0.0
    return brentq(lambda m: reflex_a0_mas(m1, m, P_days, 1.0) - target, lo, hi)


def draw_systems(z, n_bd, n_bin, seed=SEED):
    from gaiabias.constants import MJUP_MSUN
    from gaiabias.orbits import draw_ecc_kipping
    rng = np.random.default_rng(seed)
    n_star = len(z["star_G"])
    rows = []
    for kind, n in (("brown dwarf", n_bd), ("stellar", n_bin)):
        for i in range(n):
            s = int(rng.integers(n_star))
            m1 = float(z["star_m_star"][s])
            P = float(np.exp(rng.uniform(*np.log(P_WIN_YR)))) * 365.25
            if kind == "brown dwarf":
                m2 = float(np.exp(rng.uniform(np.log(13), np.log(80))))
                fac, e = 1.0, float(draw_ecc_kipping(rng))
            else:
                qmin = max(0.1, 0.08 / m1)
                q = float(rng.uniform(qmin, 1.0)) if qmin < 1 else 1.0
                m2 = q * m1 / MJUP_MSUN
                fac, e = float(photocentre_factor(np.array([q]), m1)[0]), float(rng.uniform(0, 0.8))
            omega = float(rng.uniform(0, 2 * np.pi)) + (np.pi if fac < 0 else 0.0)
            m_eff = effective_mass_mj(m1, m2, P, fac)
            rows.append(dict(id=f"{kind[0]}{i}", kind=kind, star=s, m_true=m2, P=P, e=e, factor=fac, m_eff=m_eff,
                             inc=float(np.arccos(rng.uniform(-1, 1))), Omega=float(rng.uniform(0, np.pi)),
                             omega=omega % (2 * np.pi), Tp=float(rng.uniform(0, P))))
    return rows


def run_one(r, z, synthetic):
    import step3 as S3
    from gaiabias.orbits import Companion
    from gaiabias.simulate import snr1
    from gaiabias.orbits import reflex_a0_mas
    s = r["star"]
    star = {k: float(z[f"star_{k}"][s]) for k in ("m_star", "G", "parallax", "ra", "dec", "pmra", "pmdec")}
    out = dict(id=r["id"])
    a0 = reflex_a0_mas(star["m_star"], r["m_eff"], r["P"], star["parallax"]) if r["m_eff"] > 0 else 0.0
    out["snr1"] = snr1(a0, star["G"])
    if out["snr1"] < 0.8:
        out["fitted"] = 0
        return out
    comp = Companion(r["m_eff"], r["P"], r["e"], r["inc"], r["Omega"], r["omega"], r["Tp"])
    payload = dict(job=r["id"], kind="M", star=star, comps=[comp], seed=[SEED, 7, s, int(r["id"][1:]),
                                                                       0 if r["kind"] == "stellar" else 1])
    out.update(S3.run_job(payload, FITTED, synthetic))
    out["fitted"] = 1
    return out


def simulate_and_fit(z, meta, n_bd, n_bin, jobs, synthetic, out):
    from joblib import Parallel, delayed
    import step3 as S3
    sys_path, fit_path = out / "systems.csv", out / "fits.csv"
    if not sys_path.exists():
        rows = draw_systems(z, n_bd, n_bin)
        with open(sys_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    with open(sys_path) as f:
        rows = [dict(r, star=int(r["star"]), **{k: float(r[k]) for k in ("m_true", "P", "e", "factor", "m_eff",
                                                                        "inc", "Omega", "omega", "Tp")})
                for r in csv.DictReader(f)]
    cols = ["id", "fitted", "snr1"] + S3.fit_columns(FITTED)[1:]
    done = set()
    if fit_path.exists():
        with open(fit_path) as f:
            done = {r["id"] for r in csv.DictReader(f)}
    else:
        with open(fit_path, "w", newline="") as f:
            csv.DictWriter(f, fieldnames=cols).writeheader()
    todo = [r for r in rows if r["id"] not in done]
    print(f"{len(done)} systems done, {len(todo)} to go", flush=True)
    t = C.Timer()
    chunk = 400
    for i in range(0, len(todo), chunk):
        res = Parallel(n_jobs=jobs)(delayed(run_one)(r, z, synthetic) for r in todo[i:i + chunk])
        with open(fit_path, "a", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
            for r in res:
                w.writerow({k: r.get(k, "") for k in cols})
        print(f"  {min(i + chunk, len(todo))}/{len(todo)} ({t})", flush=True)


def contaminant_table(z, meta, out):
    """One row per planet-mass orbit reported for a brown dwarf or binary that passes the DR3 RUWE cut."""
    from gaiabias.step2 import analyses as A
    from gaiabias.step2.outcomes import bic_choice, reported
    with open(out / "systems.csv") as f:
        sysrows = {r["id"]: r for r in csv.DictReader(f)}
    with open(out / "fits.csv") as f:
        fits = list(csv.DictReader(f))
    n_sim = {"brown dwarf": 0, "stellar": 0}
    for r in sysrows.values():
        n_sim[r["kind"]] += 1
    fitted = [r for r in fits if r["fitted"] == "1"]
    ruwe_max = meta.get("ruwe_max") or 1.4
    keep = [r for r in fitted if r.get("ruwe_dr3", "") not in ("", "nan") and float(r["ruwe_dr3"]) < ruwe_max]
    info = dict(n_sim=n_sim, n_fitted={k: sum(sysrows[r["id"]]["kind"] == k for r in fitted) for k in n_sim},
                n_ruwe_ok={k: sum(sysrows[r["id"]]["kind"] == k for r in keep) for k in n_sim})
    if not keep:
        return {}, info
    d = {}
    for k in keep[0]:
        if k == "id" or k in A.TEXT_FIELDS:
            d[k] = np.array([r[k] for r in keep])
        else:
            d[k] = np.array([float(r[k]) if r[k] not in ("", "nan") else np.nan for r in keep])
    rep = reported(d, "best", "w")
    bic = bic_choice(d, "w")
    t = {k: [] for k in ("P", "m", "sig", "bic", "e1", "kind", "m_true", "factor", "star")}
    for mask, P, m, sg in rep:
        for i in np.nonzero(mask)[0]:
            s = sysrows[d["id"][i]]
            t["P"].append(P[i] / 365.25)
            t["m"].append(m[i])
            t["sig"].append(sg[i])
            t["bic"].append(bic[i])
            t["e1"].append(d["one_w_e"][i] if "one_w_e" in d else np.nan)
            t["kind"].append(s["kind"])
            t["m_true"].append(float(s["m_true"]))
            t["factor"].append(float(s["factor"]))
            t["star"].append(int(s["star"]))
    t = {k: np.array(v) for k, v in t.items()}
    info["n_reported"] = {k: int((t["kind"] == k).sum()) for k in n_sim}
    return t, info


def wscores(flagged, bad, w):
    tp, fp = np.sum(w * (flagged & bad)), np.sum(w * (flagged & ~bad))
    fn, tn = np.sum(w * (~flagged & bad)), np.sum(w * (~flagged & ~bad))
    return dict(recall=tp / max(tp + fn, 1e-12), precision=tp / max(tp + fp, 1e-12),
                fpr=fp / max(fp + tn, 1e-12), share=(tp + fp) / max(np.sum(w), 1e-12), bad=(tp + fn) / max(np.sum(w), 1e-12))


def analyse(root, out):
    import g1_flags_heldout as G1
    import step3 as S3
    import step4 as S4

    z, meta, jobs, d = C.load_step3(root)
    ct, info = contaminant_table(z, meta, out)
    rates = {}
    for kind in ("brown dwarf", "stellar"):
        if kind == "brown dwarf":
            rate = BD_RATE
        else:
            m = z["star_m_star"]
            rate = float(np.mean(np.where(m > 0.6, BIN_FRAC[0], BIN_FRAC[1]))) * window_share()
        rates[kind] = rate
    n_star = meta["n_stars"]
    L = ["=== WARNING FLAGS WITH BROWN-DWARF AND STELLAR-BINARY CONTAMINANTS ===",
         f"step 3 folder: {root}; {n_star:,} real stars; best model by BIC (DR4-like).",
         f"Contaminants simulated on random step 3 stars: {info['n_sim']}; fitted (summed SNR >= 0.8): "
         f"{info['n_fitted']}; passing the DR3 RUWE cut: {info['n_ruwe_ok']}; reported as planets (m < 13 M_J): "
         f"{info.get('n_reported', {})}.",
         f"Rates per star in {P_WIN_YR[0]}-{P_WIN_YR[1]:g} yr: brown dwarf {100 * rates['brown dwarf']:.2f}%, "
         f"stellar companion {100 * rates['stellar']:.1f}% (binary fraction x {window_share():.2f} of the period "
         "distribution).", ""]
    if ct:
        for kind in ("brown dwarf", "stellar"):
            sel = ct["kind"] == kind
            if sel.any():
                L.append(f"  {kind}s reported as planets: reported mass median {np.median(ct['m'][sel]):.1f} M_J "
                         f"(true companion median {np.median(ct['m_true'][sel]):.0f} M_J); adopted model plain orbit "
                         f"{100 * np.mean(ct['bic'][sel] == 3):.0f}%; period median {np.median(ct['P'][sel]):.2f} yr")
        L.append("")
    res = {}
    for scen in S3.SCENARIOS:
        t = S4.planet_table(z, meta, jobs, d, "best", scen, 2.0)
        bad_p = t["cls"] >= 1
        rng = np.random.default_rng([2027, S3.SCENARIOS.index(scen)])
        us = np.unique(t["star"])
        train_stars = set(rng.choice(us, len(us) // 2, replace=False).tolist())
        tr_p = np.array([s in train_stars for s in t["star"]])
        nc = len(ct["P"]) if ct else 0
        rngc = np.random.default_rng([2027, 9, S3.SCENARIOS.index(scen)])
        tr_c = rngc.random(nc) < 0.5
        w_c = np.array([rates[k] * n_star / info["n_sim"][k] for k in ct["kind"]]) if nc else np.zeros(0)
        comb = {k: np.concatenate([np.asarray(t[k], float), np.asarray(ct[k], float) if nc else np.zeros(0)])
                for k in ("P", "sig", "bic", "e1")}
        bad = np.concatenate([bad_p, np.ones(nc, bool)])
        is_c = np.concatenate([np.zeros(len(bad_p), bool), np.ones(nc, bool)])
        tr = np.concatenate([tr_p, tr_c])
        te = ~tr

        def weights(scale):
            return np.concatenate([np.ones(len(bad_p)), scale * w_c])

        w1 = weights(1.0)
        best, best_j = None, -1
        for cut in itertools.product(G1.SIG_GRID, G1.PF_GRID, G1.E_GRID):
            f = G1.flags(comb, *cut).any(axis=0)
            sc = wscores(f[tr], bad[tr], w1[tr])
            if sc["recall"] - sc["fpr"] > best_j:
                best, best_j = cut, sc["recall"] - sc["fpr"]
        rows = {}
        for name, cut in (("default", G1.DEFAULT), ("tuned on mixed training half", best)):
            f = G1.flags(comb, *cut).any(axis=0)
            rows[name] = dict(cut=cut, planet=wscores(f[te & ~is_c], bad[te & ~is_c], w1[te & ~is_c]),
                              **{f"mixed x{s:g}": wscores(f[te], bad[te], weights(s)[te]) for s in RATE_SCALES})
            rows[name]["caught_contaminants"] = float(np.average(f[te & is_c], weights=w1[te & is_c])) \
                if (te & is_c).any() else np.nan
        nflag = G1.flags(comb, *G1.DEFAULT).sum(axis=0)
        cal = []
        for k in range(5):
            sel = (nflag == k) if k < 4 else (nflag >= 4)
            cal.append((k, *(float(np.average(bad[te & sel & m], weights=w1[te & sel & m]))
                             if (te & sel & m).sum() >= 10 else np.nan for m in (~is_c, np.ones_like(is_c)))))
        res[scen] = dict(rows=rows, cal=cal)
        cont_share = np.sum(w1[te & is_c]) / max(np.sum(w1[te]), 1e-12)
        L.append(f"--- {scen} multiplicity: {len(bad_p)} reported planets from planet systems, {nc} from "
                 f"contaminants; contaminants are {100 * cont_share:.1f}% of reported planets (rates x1) ---")
        for name, r in rows.items():
            c = r["cut"]
            L.append(f"  {name}: significance < {c[0]:g}, period > {c[1]:g} x span, eccentricity > {c[2]:g}; "
                     f"flags {100 * r['caught_contaminants']:.0f}% of contaminant 'planets'")
            for pop in ["planet"] + [f"mixed x{s:g}" for s in RATE_SCALES]:
                s = r[pop]
                lab = "planet-only" if pop == "planet" else f"mixed, contaminant rates {pop[7:]}"
                L.append(f"    {lab:<34} recall {100 * s['recall']:5.1f}%  precision {100 * s['precision']:5.1f}%  "
                         f"FPR {100 * s['fpr']:5.1f}%  flagged {100 * s['share']:5.1f}%  inaccurate "
                         f"{100 * s['bad']:5.1f}%")
        L.append("  calibration, test half (default thresholds): P(inaccurate | flags) planet-only / mixed")
        for k, pp, pm in cal:
            L.append(f"    {k if k < 4 else '4+':>2} flags: {pp:.2f} / {pm:.2f}")
        L.append("")
    L.append("Contaminant rates and the binary photocentre model are simple assumptions; the 'mixed' rows show how "
             "much the scores change, not a prediction of the real DR4 contamination.")
    C.write_summary(out / "summary.txt", L)

    plt, ps = C.setup_plots()
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8), constrained_layout=True)
    metrics = ["recall", "precision", "fpr"]
    x = np.arange(len(metrics))
    cols = [ps.C_TARGET, "#8fb8f0", ps.C_WRONG, "#f3b98c"]
    for i, scen in enumerate(S3.SCENARIOS):
        for j, pop in enumerate(["planet", "mixed x1"]):
            r = res[scen]["rows"]["default"][pop]
            ax[0].bar(x + (2 * i + j - 1.5) * 0.19, [r[m] for m in metrics], 0.18, color=cols[2 * i + j],
                      label=f"{scen}, {'planet-only' if pop == 'planet' else 'with BDs and binaries'}")
    ax[0].set_xticks(x, ["recall", "precision", "false-positive rate"])
    ps.percent_axis(ax[0], "y", "test half, default rule")
    ax[0].set_ylim(0, 1)
    ax[0].legend(fontsize=8)
    ax[0].set_title("Suspect rule with and without non-planet companions", fontsize=10.5)
    for i, scen in enumerate(S3.SCENARIOS):
        cal = res[scen]["cal"]
        k = [c[0] for c in cal]
        ax[1].plot(k, [c[1] for c in cal], "o--", mfc="white", color=[ps.C_TARGET, ps.C_WRONG][i],
                   label=f"{scen}: planet-only")
        ax[1].plot(k, [c[2] for c in cal], "o-", color=[ps.C_TARGET, ps.C_WRONG][i], label=f"{scen}: mixed")
    ax[1].set_xticks(range(5), ["0", "1", "2", "3", "4+"])
    ax[1].set_xlabel("number of warning flags raised")
    ps.percent_axis(ax[1], "y", "share of reported planets not accurate")
    ax[1].set_ylim(0, 1)
    ax[1].legend(fontsize=8)
    ax[1].set_title("Calibration", fontsize=10.5)
    fig.savefig(out / "mixed_flags.png", dpi=150, bbox_inches="tight")
    print(f"saved {out / 'mixed_flags.png'}")


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--root", type=C.Path, default=C.STEP3)
    p.add_argument("--n-bd", type=int, default=3000)
    p.add_argument("--n-bin", type=int, default=6000)
    p.add_argument("-j", "--jobs", type=int, default=-1)
    p.add_argument("--quick", action="store_true")
    p.add_argument("--synthetic", action="store_true", help="fake scanning law (code test only)")
    p.add_argument("--analyse", action="store_true")
    a = p.parse_args()
    out = C.out_dir("mixed_flags_quick" if a.quick else "mixed_flags")
    z, meta, _, _ = C.load_step3(a.root)
    if not a.analyse:
        n_bd, n_bin = (300, 600) if a.quick else (a.n_bd, a.n_bin)
        simulate_and_fit(z, meta, n_bd, n_bin, a.jobs, a.synthetic, out)
    analyse(a.root, out)


if __name__ == "__main__":
    main()
