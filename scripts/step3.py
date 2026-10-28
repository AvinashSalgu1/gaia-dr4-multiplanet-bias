"""Step 3: do multi-planet systems bias the giant-planet occurrence rate inferred from Gaia DR4?

ONE command does everything, and every stage resumes where it stopped:

    python scripts/step3.py                 # every selected star over the whole sky -> results/step3_allsky/
    python scripts/step3.py --stars 300000  # a random subset of them
    python scripts/step3.py --analyse       # only redo the summary and figures
    python scripts/step3.py --download-only # only fetch the Gaia stars (light on the CPU)
    python scripts/step3.py --circles --stars 300000   # the 25% sky sample of versions <= 0.2.0
                                                       # -> results/step3_gaia/

Stages:
 1. STARS   Gaia DR3 stars closer than 250 pc, selected as Lammers & Winn
            Sec. 4.1 (main sequence, their distance limit), RUWE <= 1.4.
            Downloaded once with astroquery into data/gaia_dr3_250pc/ (skipped
            if already there): the whole sky by default, in 768 HEALPix tiles.
            See src/gedi/gaia_stars.py.
 2. BUILD   a random subset of those stars; giant planets from the Lammers &
            Winn Eq. 23 occurrence model (CLS-normalised; 0.3-13 M_J, 0.1-30 au)
            under two multiplicity scenarios with the SAME mean number of
            planets per star (independent / clustered as Bryan et al. 2016),
            all multi-planet systems Hill-stable. Jobs to fit:
              A  each star with ALL its planets (what Gaia sees)
              B  each planet of a multi-planet star ALONE (same epochs, same noise)
              V  validation: systems just below the pre-fit cut
              I  single-planet injections -> completeness, as a DR4 user would
              N  stars without planets -> RUWE from noise alone
 3. FIT     each job analysed with the orbit models Gaia DR4 publishes: one
            planet (gaiamock, Lammers & Winn criteria), + acceleration, + jerk,
            two planets; the "best model by BIC" (DR4-like choice) is picked
            from those fits (as in step 2). Every reported planet also gets the
            significance a0/sigma_a0 (DR4 flags < 5 as poorly constrained).
            Also the RUWE the star would have in DR4 and in DR3.
 4. ANALYSE occurrence rates in mass x period bins from A (reality) and B
            (planets alone) with the completeness from I (measured at the
            FITTED parameters, so planets scattering between bins are
            accounted for), a power law forward-modelled through the injection
            response, compared with the true rate; the RUWE study; and a
            MULTIPLICITY SCAN: occurrence bias and DR4 yield versus the fraction
            of giant hosts with a second giant (0-80%) at the same CLS rate,
            built by reweighting the fitted systems (no new fits; see
            src/gedi/multiplicity.py). Summary + figures in results/step3_gaia/.

Outputs: results/step3_gaia/summary.txt and the figures results/step3_gaia/*.png.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

RESULTS = Path(__file__).resolve().parents[1] / "results" / "step3_gaia"           # --circles (<= 0.2.0)
RESULTS_ALLSKY = Path(__file__).resolve().parents[1] / "results" / "step3_allsky"   # default (tiles)
SCENARIOS = ["independent", "clustered"]
FITTED = ["one", "accel", "jerk", "two"]          # DR4's orbit models (step 2 code)
JUDGED = ["one", "accel", "jerk", "two", "best"]  # "best" = lowest BIC among all models (DR4-like choice)
LABEL = {"one": "one planet (L&W)", "accel": "one planet + acceleration", "jerk": "one planet + jerk",
         "two": "two planets (if justified)", "best": "best model by BIC (DR4-like)", "dr3": "DR3-style pipeline"}
S_CUT = 0.8                       # summed SNR_1 below which a star is not fitted (Delta chi2 ~ 22 SNR^2 << 50)
VAL_RANGE = (0.5, S_CUT)          # validation sample just below the cut
N_VAL = 400                       # validation stars per scenario
N_NULL = 1500                     # stars without planets, for the RUWE null distribution
M_EDGES = np.array([0.3, 1.0, 3.0, 13.0])                 # M_J
P_EDGES = np.array([0.1, 0.3, 1.0, 2.0, 4.0, 8.0])        # years
M_FINE = np.geomspace(0.3, 13.0, 9)
P_FINE = np.geomspace(0.1, 8.0, 13)
INJ_P_YR = (0.05, 16.0)           # injections extend past the shown window (planets scatter between bins)
MIN_ETA, MIN_DET = 0.05, 10       # a bin is "reliable" above these; others are drawn grey with their counts
RUWE_CUT = 1.4
SEC_PER_JOB = 4.4                 # one + accel + jerk + two fits, one core (rough, from step 2)


def paths(root: Path):
    return dict(pop=root / "population.npz", jobs=root / "jobs.csv", meta=root / "meta.json",
                fits=root / "fits.csv", summary=root / "summary.txt",
                rates_png=root / "step3_rates.png", overview_png=root / "step3_overview.png",
                ruwe_png=root / "step3_ruwe.png", scan_png=root / "step3_multiplicity.png",
                scan_csv=root / "multiplicity_scan.csv")


# ============================================================================ build
def _snr_sum(stars, pl):
    from gedi.population import planet_snr1
    snr = planet_snr1(stars, pl)
    return snr, np.bincount(pl.star, weights=snr, minlength=len(stars)), np.bincount(pl.star, minlength=len(stars))


def _draw_injections(stars, rng, target_fits, block=100_000):
    """Single planets, log-uniform in m and P, until `target_fits` pass the cut.
    Hosts drawn proportional to M_star (the occurrence model's stellar-mass scaling)."""
    from gedi.orbits import draw_ecc_kipping, reflex_a0_mas
    from gedi.simulate import sigma_al_per_transit

    w = stars.m_star / stars.m_star.sum()
    rows, idx, all_m, all_P = [], [], [], []
    n_fit = offset = 0
    while n_fit < target_fits:
        s = rng.choice(len(stars), block, p=w)
        m = np.exp(rng.uniform(np.log(0.3), np.log(13.0), block))
        P = 365.25 * np.exp(rng.uniform(*np.log(INJ_P_YR), block))
        snr = reflex_a0_mas(stars.m_star[s], m, P, stars.parallax[s]) / sigma_al_per_transit(stars.G[s])
        e = draw_ecc_kipping(rng, block)
        inc = np.arccos(rng.uniform(-1, 1, block))
        Om, om = rng.uniform(0, np.pi, block), rng.uniform(0, 2 * np.pi, block)
        Tp = rng.uniform(0, 1, block) * P
        above = np.nonzero(snr >= S_CUT)[0]
        stop = block
        if len(above) > target_fits - n_fit:
            above = above[:target_fits - n_fit]
            stop = above[-1] + 1
        rows.append(np.column_stack([s, m, P, e, inc, Om, om, Tp, snr])[above])
        idx.append(offset + above)
        all_m.append(m[:stop].astype(np.float32))
        all_P.append(P[:stop].astype(np.float32))
        offset += stop
        n_fit += len(above)
    return np.concatenate(rows), np.concatenate(idx), np.concatenate(all_m), np.concatenate(all_P)


def build(root: Path, n_stars: int, inj_fits: int, seed: int, ruwe_max, stars_override=None,
          sky=("circles", None), variant=(), n_jobs: int = 1):
    """sky = (layout, sky fraction asked for); n_stars <= 0 means every selected star;
    variant: population variants of gedi.population (coplanar, stable, lowecc)."""
    from gedi.population import draw_planets, host_fraction_clustered

    P = paths(root)
    root.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    if stars_override is None:
        from gedi import gaia_stars as GS
        stars, info = GS.load_sample(n_stars, seed, ruwe_max, layout=sky[0])
    else:
        stars, info = stars_override
    save = {f"star_{k}": getattr(stars, k) for k in stars.__dataclass_fields__}
    save["star_ruwe_dr3_real"] = info["ruwe"]
    jobs = []
    meta = dict(n_stars=len(stars), n_stars_asked=n_stars, sky=list(sky), seed=seed, ruwe_max=ruwe_max,
                s_cut=S_CUT, inj_fits=inj_fits, variant=sorted(variant),
                host_fraction=host_fraction_clustered(), selection=info["counts"],
                n_selected=info["n_selected"], frac_flame_mass=info["frac_flame_mass"], scen={})
    for si, scen in enumerate(SCENARIOS):
        vinfo = {}
        pl = draw_planets(stars, np.random.default_rng([seed, 1 + si]), scen, variant=variant,
                          jobs=n_jobs, info=vinfo)
        for k in ("star", "m", "P", "a", "e", "inc", "Omega", "omega", "Tp"):
            save[f"{scen}_{k}"] = getattr(pl, k)
        snr, ssum, npl = _snr_sum(stars, pl)
        save[f"{scen}_snr"] = snr
        fit_stars = np.nonzero(ssum >= S_CUT)[0]
        multi = npl[pl.star] >= 2
        order = np.argsort(pl.star, kind="stable")
        starts = np.searchsorted(pl.star[order], np.arange(len(stars) + 1))

        def planets_of(s):
            return order[starts[s]:starts[s + 1]]

        for s in fit_stars:
            jobs.append(dict(job=f"A-{scen}-{s}", kind="A", scen=scen, star=int(s),
                             planets=";".join(map(str, planets_of(s)))))
        for i in np.nonzero(multi & (snr >= S_CUT))[0]:
            jobs.append(dict(job=f"B-{scen}-{i}", kind="B", scen=scen, star=int(pl.star[i]), planets=str(i)))
        rng_v = np.random.default_rng([seed, 10 + si])
        below = np.nonzero((ssum >= VAL_RANGE[0]) & (ssum < VAL_RANGE[1]))[0]
        for s in rng_v.choice(below, size=min(N_VAL, len(below)), replace=False):
            jobs.append(dict(job=f"V-{scen}-{s}", kind="V", scen=scen, star=int(s),
                             planets=";".join(map(str, planets_of(s)))))
        meta["scen"][scen] = dict(n_planets=len(pl), dropped=pl.n_dropped, redrawn=pl.n_redrawn,
                                  hosts=int((npl > 0).sum()), multi_hosts=int((npl >= 2).sum()),
                                  n_A=int(len(fit_stars)), n_B=int((multi & (snr >= S_CUT)).sum()),
                                  n_V=int(min(N_VAL, len(below))), stability=vinfo)
    rng_n = np.random.default_rng([seed, 30])
    for s in rng_n.choice(len(stars), size=min(N_NULL, len(stars)), replace=False):
        jobs.append(dict(job=f"N-{s}", kind="N", scen="", star=int(s), planets=""))
    inj, inj_idx, all_m, all_P = _draw_injections(stars, np.random.default_rng([seed, 20]), inj_fits)
    save.update(inj=inj, inj_idx=inj_idx, inj_all_m=all_m, inj_all_P=all_P)
    for i in range(len(inj)):
        jobs.append(dict(job=f"I-{i}", kind="I", scen="", star=int(inj[i, 0]), planets=str(i)))
    meta["n_inj"], meta["n_I"] = int(len(all_m)), int(len(inj))
    np.savez_compressed(P["pop"], **save)
    with open(P["jobs"], "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["job", "kind", "scen", "star", "planets"])
        w.writeheader()
        w.writerows(jobs)
    P["meta"].write_text(json.dumps(meta, indent=1))

    d = 1000 / stars.parallax
    print(f"built in {time.perf_counter() - t0:.0f} s: {len(stars)} real stars (of {info['n_selected']} selected; "
          f"distance 5/50/95% {np.percentile(d, 5):.0f}/{np.median(d):.0f}/{np.percentile(d, 95):.0f} pc; "
          f"mass from Gaia FLAME for {100 * info['frac_flame_mass']:.0f}%)")
    for k, v in info["counts"].items():
        print(f"    {k}: {v}")
    for scen, dd in meta["scen"].items():
        print(f"  {scen:<11} planets/star {dd['n_planets'] / len(stars):.3f}  hosts {dd['hosts']}  multi-planet hosts "
              f"{dd['multi_hosts']}  fits: A {dd['n_A']}  B {dd['n_B']}  validation {dd['n_V']}")
    n_fit = sum(j["kind"] != "N" for j in jobs)
    ncpu = os.cpu_count() or 1
    print(f"  injections to fit: {meta['n_I']};  total fits: {n_fit} -> roughly "
          f"{n_fit * SEC_PER_JOB / 3600 / ncpu:.1f} h on the {ncpu} cores of this machine")


def built_matches(root: Path, n_stars, inj_fits, seed, ruwe_max, sky=("circles", None), variant=()) -> bool:
    P = paths(root)
    if not (P["meta"].exists() and P["pop"].exists() and P["jobs"].exists()):
        return False
    m = json.loads(P["meta"].read_text())
    if list(m.get("sky", ["circles", None])) != list(sky) or m.get("variant", []) != sorted(variant):
        return False
    if n_stars <= 0:
        same_n = m.get("n_stars_asked", m.get("n_stars")) == n_stars or m.get("n_stars") == m.get("n_selected")
    else:
        same_n = m.get("n_stars") == min(n_stars, m.get("n_selected", n_stars))
    return (same_n and m.get("inj_fits") == inj_fits and m.get("seed") == seed
            and m.get("ruwe_max", "unset") == ruwe_max)


# ============================================================================ fit
def _payloads(root: Path):
    from gedi.orbits import Companion

    P = paths(root)
    with np.load(P["pop"]) as f:
        z = {k: f[k] for k in f.files if not k.startswith("inj_all")}
    meta = json.loads(P["meta"].read_text())
    with open(P["jobs"]) as f:
        jobs = list(csv.DictReader(f))
    out = []
    for j in jobs:
        s = int(j["star"])
        star = {k: float(z[f"star_{k}"][s]) for k in ("m_star", "G", "parallax", "ra", "dec", "pmra", "pmdec")}
        if j["kind"] == "N":
            comps, seed = [], [meta["seed"], 60, s]
        elif j["kind"] == "I":
            r = z["inj"][int(j["planets"])]
            comps = [Companion(r[1], r[2], r[3], r[4], r[5], r[6], r[7])]
            seed = [meta["seed"], 99, int(j["planets"])]
        else:
            sc = j["scen"]
            comps = [Companion(*(float(z[f"{sc}_{k}"][i]) for k in ("m", "P", "e", "inc", "Omega", "omega", "Tp")))
                     for i in map(int, j["planets"].split(";"))]
            seed = [meta["seed"], 50 + SCENARIOS.index(sc), s]      # A, B and V of one star share epochs and noise
        out.append(dict(job=j["job"], kind=j["kind"], star=star, comps=comps, seed=seed))
    return out


def _data(star, comps, rng, dr, synthetic):
    from gedi import epochs as E
    from gedi.simulate import Star, sigma_al_per_transit, simulate
    ep = E.synthetic(n_transits=80 if dr == "dr4" else 45, rng=rng) if synthetic else \
        E.from_gaiamock(star["ra"], star["dec"], dr, rng=rng)
    st = Star(parallax_mas=star["parallax"], m_star_msun=star["m_star"], G_mag=star["G"],
              pmra=star["pmra"], pmdec=star["pmdec"])
    sim = simulate(ep, st, comps, add_noise=False)
    noise = float(sigma_al_per_transit(star["G"])) * rng.standard_normal(len(ep))
    if star["G"] < 13:
        noise += rng.uniform(0, 0.04) * rng.standard_normal(len(ep))
    return ep, sim.obs + noise, sim.err


def run_job(payload, which, synthetic=False) -> dict:
    from gedi.fitting import fit_5par, ruwe_from_chi2
    from gedi.step2 import analyses as A

    rng = np.random.default_rng(payload["seed"])
    ep, obs, err = _data(payload["star"], payload["comps"], rng, "dr4", synthetic)
    row = dict(job=payload["job"], n_transits=len(ep),
               ruwe_dr4=ruwe_from_chi2(fit_5par(ep, obs, err)[0], len(ep)))
    try:   # the same planets seen by DR3 (34 months), independent noise
        ep3, obs3, err3 = _data(payload["star"], payload["comps"], np.random.default_rng(payload["seed"] + [3]),
                                "dr3", synthetic)
        row["ruwe_dr3"] = ruwe_from_chi2(fit_5par(ep3, obs3, err3)[0], len(ep3))
    except Exception:
        row["ruwe_dr3"] = np.nan
    if payload["kind"] != "N":
        row.update(A.run(ep, obs, err, payload["star"]["m_star"], "w", set(which)))
    return row


def fit_columns(which):
    from gedi.step2 import analyses as A
    return ["job", "n_transits", "ruwe_dr4", "ruwe_dr3"] + [c for a in which for c in A.cols(a, "w")]


def fit(root: Path, which, n_jobs: int, chunk: int = 600, synthetic=False, part=None):
    """Fit all jobs not yet in fits.csv. part=(k, n): only jobs with index % n == k, written to
    fits.part<k>of<n>.csv (for running the n parts on separate machines or cluster nodes; join them
    afterwards with merge_parts)."""
    from joblib import Parallel, delayed

    P = paths(root)
    if part is not None:
        P = dict(P, fits=root / f"fits.part{part[0]}of{part[1]}.csv")
    cols = fit_columns(which)
    done = set()
    if P["fits"].exists():
        with open(P["fits"]) as f:
            rd = csv.DictReader(f)
            if set(cols) - set(rd.fieldnames or []):
                raise SystemExit(f"{P['fits']} was made with different analyses; move it away to refit")
            done = {r["job"] for r in rd}
    else:
        with open(P["fits"], "w", newline="") as f:
            csv.DictWriter(f, fieldnames=cols).writeheader()
    allp = _payloads(root)
    if part is not None:
        allp = [p for i, p in enumerate(allp) if i % part[1] == part[0]]
    todo = [p for p in allp if p["job"] not in done]
    rank = {"A": 0, "V": 1, "B": 2, "I": 3, "N": 4}         # partial runs are still usable; N (no fit) is quick
    todo.sort(key=lambda p: rank[p["kind"]])
    print(f"fitting: {len(done)} jobs already done, {len(todo)} to go", flush=True)
    t0 = time.perf_counter()
    for i in range(0, len(todo), chunk):
        rows = Parallel(n_jobs=n_jobs)(delayed(run_job)(p, which, synthetic) for p in todo[i:i + chunk])
        with open(P["fits"], "a", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
            for r in rows:
                w.writerow({k: r.get(k, "") for k in cols})
        n = min(i + chunk, len(todo))
        el = (time.perf_counter() - t0) / 60
        print(f"  {n}/{len(todo)} ({el:.1f} min, ~{el / n * (len(todo) - n):.0f} min left)", flush=True)


def merge_parts(root: Path):
    """Join fits.part*of*.csv (from --part runs) into fits.csv (rows already in fits.csv are kept)."""
    parts = sorted(root.glob("fits.part*of*.csv"))
    if not parts:
        raise SystemExit(f"no fits.part*of*.csv files in {root}")
    P = paths(root)
    seen, rows, header = set(), [], None
    for f in ([P["fits"]] if P["fits"].exists() else []) + parts:
        with open(f) as fh:
            rd = csv.DictReader(fh)
            header = header or rd.fieldnames
            for r in rd:
                if r["job"] not in seen:
                    seen.add(r["job"])
                    rows.append(r)
    with open(P["fits"], "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header)
        w.writeheader()
        w.writerows(rows)
    print(f"merged {len(parts)} part files -> {P['fits']} ({len(rows)} fitted jobs)")


# ============================================================================ analyse
def load_fits(root: Path):
    from gedi.step2 import analyses as A
    with open(paths(root)["fits"]) as f:
        rows = list(csv.DictReader(f))
    text = A.TEXT_FIELDS
    d = {}
    for k in rows[0].keys() if rows else []:
        if k == "job" or k in text:
            d[k] = np.array([r[k] for r in rows])
        else:
            d[k] = np.array([float(r[k]) if r[k] not in ("", "nan") else np.nan for r in rows])
    return d


def reported_planets(d, analysis):
    """Per job: list of (fitted P [yr], m [M_J], significance) of every planet the analysis reports."""
    from gedi.step2.outcomes import reported
    rep = reported(d, analysis, "w")
    n = len(d["job"])
    out = [[] for _ in range(n)]
    for mask, P, m, sig in rep:
        for i in np.nonzero(mask)[0]:
            out[i].append((P[i] / 365.25, m[i], sig[i]))
    return out


TOL_OK = (np.log(1.1), np.log(1.2))       # "correct": period within 10%, mass within 20%
TOL_ROUGH = (np.log(1.5), np.log(2.0))    # "imprecise": same planet, but period within x1.5 / mass within x2 only
P_TRUE_EDGES = np.r_[INJ_P_YR[0], P_FINE, INJ_P_YR[1]]   # power-law true cells reach past the window


def classify(Py, m, P_true, m_true):
    """0 = correct, 1 = imprecise (a real planet, poorly measured), 2 = spurious (no real planet resembles it)."""
    if len(P_true) == 0:
        return 2
    dP, dm = np.abs(np.log(Py / P_true)), np.abs(np.log(m / m_true))
    if np.any((dP < TOL_OK[0]) & (dm < TOL_OK[1])):
        return 0
    if np.any((dP < TOL_ROUGH[0]) & (dm < TOL_ROUGH[1])):
        return 1
    return 2


def injection_response(z, jobs, rep_of):
    """True and fitted parameters of the single-planet injections (completeness in fitted space)."""
    m_all = z["inj_all_m"].astype(float)
    P_all = z["inj_all_P"].astype(float) / 365.25
    det_t_m, det_t_P, det_f_m, det_f_P, det_idx = [], [], [], [], []
    det_any = np.zeros(len(m_all), bool)
    n_done = 0
    for j in jobs:
        if j["kind"] != "I":
            continue
        r = rep_of(j["job"])
        if r is None:
            continue
        n_done += 1
        i = int(j["planets"])
        k = int(z["inj_idx"][i])
        det_any[k] = len(r) > 0
        for Py, m, _ in r:
            det_t_m.append(m_all[k]), det_t_P.append(P_all[k]), det_f_m.append(m), det_f_P.append(Py)
            det_idx.append(k)
    arr = lambda x: np.array(x, float)  # noqa: E731
    return dict(m_all=m_all, P_all=P_all, det_any=det_any, n_done=n_done, t_m=arr(det_t_m), t_P=arr(det_t_P),
                f_m=arr(det_f_m), f_P=arr(det_f_P), idx=np.array(det_idx, int))


def eta_fitted(inj, w_all=None):
    """Completeness per coarse bin with detections at FITTED parameters (optionally bootstrap-weighted)."""
    from gedi import rates as R
    w_det = None if w_all is None else w_all[inj["idx"]]
    return R.completeness_fitted(inj["m_all"], inj["P_all"], None, inj["f_m"], inj["f_P"], M_EDGES, P_EDGES,
                                 w_inj=None if w_all is None else (w_all, w_det))


def analyse_one(z, meta, jobs, d, analysis):
    """Rates and diagnostics for one analysis, both scenarios; also per-host tables for the multiplicity scan."""
    from gedi import rates as R

    n_star = meta["n_stars"]
    pos = {j: k for k, j in enumerate(d["job"])}
    reps = reported_planets(d, analysis)

    def rep_of(job):
        k = pos.get(job)
        return reps[k] if k is not None else None          # None = not fitted (yet)

    inj = injection_response(z, jobs, rep_of)
    eta_old, _, _ = R.completeness(inj["m_all"], inj["P_all"], inj["det_any"], M_EDGES, P_EDGES)
    eta, eta_err, _ = eta_fitted(inj)
    Rm = R.response_matrix(inj["m_all"], inj["P_all"], inj["t_m"], inj["t_P"], inj["f_m"], inj["f_P"],
                           M_FINE, P_TRUE_EDGES, M_FINE, P_FINE)

    out = dict(I_done=inj["n_done"], scen={}, eta=eta, eta_old=eta_old, inj=inj)
    for scen in SCENARIOS:
        pm, pP, pst = z[f"{scen}_m"], z[f"{scen}_P"] / 365.25, z[f"{scen}_star"]
        npl = np.bincount(pst, minlength=n_star)
        order = np.argsort(pst, kind="stable")
        starts = np.searchsorted(pst[order], np.arange(n_star + 1))
        planets_of = lambda s: order[starts[s]:starts[s + 1]]  # noqa: E731

        A_rows, A_jobs_done = [], 0           # (star, P, m, sig, class)
        for j in jobs:
            if j["kind"] == "A" and j["scen"] == scen:
                r = rep_of(j["job"])
                if r is None:
                    continue
                A_jobs_done += 1
                s_ = int(j["star"])
                idx = planets_of(s_)
                for Py, m, sg in r:
                    A_rows.append((s_, Py, m, sg, classify(Py, m, pP[idx], pm[idx])))
        A = np.array(A_rows, float).reshape(-1, 5)
        A_star, A_P, A_m, A_sig, A_cls = A[:, 0].astype(int), A[:, 1], A[:, 2], A[:, 3], A[:, 4].astype(int)

        # B: each planet alone (single-planet stars: the A fit IS the planet alone)
        B_rows = []                            # (planet, P, m, class)
        single = set(np.nonzero(npl == 1)[0].tolist())
        for k in range(len(A_star)):
            if A_star[k] in single:
                i = int(planets_of(A_star[k])[0])
                B_rows.append((i, A_P[k], A_m[k], A_cls[k]))
        for j in jobs:
            if j["kind"] == "B" and j["scen"] == scen:
                i = int(j["planets"])
                for Py, m, _ in (rep_of(j["job"]) or []):
                    B_rows.append((i, Py, m, classify(Py, m, pP[[i]], pm[[i]])))
        B = np.array(B_rows, float).reshape(-1, 4)
        B_pl, B_P, B_m, B_cls = B[:, 0].astype(int), B[:, 1], B[:, 2], B[:, 3].astype(int)

        truth = R.true_rate(pm, pP, n_star, M_EDGES, P_EDGES)
        rA = R.binned_rate(A_m, A_P, n_star, eta, eta_err, M_EDGES, P_EDGES)
        rB = R.binned_rate(B_m, B_P, n_star, eta, eta_err, M_EDGES, P_EDGES)
        reliable = (eta >= MIN_ETA) & (rB.n_det >= MIN_DET)
        with np.errstate(invalid="ignore"):
            n_wrong_conf = int(np.sum((A_cls == 2) & (A_sig >= 5.0)))
            n_imp_conf = int(np.sum((A_cls == 1) & (A_sig >= 5.0)))
        in_multi = npl[A_star] >= 2 if len(A_star) else np.zeros(0, bool)
        v = []
        for j in jobs:
            if j["kind"] == "V" and j["scen"] == scen:
                r = rep_of(j["job"])
                if r is not None:
                    v.append(len(r) > 0)

        def pl_fwd(m, Pyr):
            try:
                return R.fit_power_law_forward(m, Pyr, n_star, Rm, M_FINE, P_TRUE_EDGES, M_FINE, P_FINE)
            except Exception:
                return None

        try:
            fT = R.fit_power_law(pm, pP, n_star, M_FINE, P_FINE, None)
        except Exception:
            fT = None
        with np.errstate(invalid="ignore", divide="ignore"):
            qA, eA = rA.f / truth, rA.err / truth
            qB, eB = rB.f / truth, rB.err / truth
        window = (eta >= MIN_ETA) & (truth > 0)          # bins where a DR4 user can correct for completeness
        cls_counts = lambda c: np.bincount(c, minlength=3) if len(c) else np.zeros(3, int)  # noqa: E731
        out["scen"][scen] = dict(
            truth=truth, rA=rA, rB=rB, qA=qA, eA=eA, qB=qB, eB=eB, reliable=reliable, window=window,
            nA=len(A_m), nB=len(B_m), clsA=cls_counts(A_cls), clsB=cls_counts(B_cls),
            n_ok_A=int(np.sum(A_cls == 0)), n_ok_B=int(np.sum(B_cls == 0)),
            n_wrong_conf=n_wrong_conf, n_imp_conf=n_imp_conf,
            frac_multi=float(in_multi.mean()) if len(in_multi) else np.nan, A_done=A_jobs_done,
            val=(int(sum(v)), len(v)), fT=fT, fA=pl_fwd(A_m, A_P), fB=pl_fwd(B_m, B_P),
            hosts=dict(npl=npl, pm=pm, pP=pP, pst=pst, m_star=z["star_m_star"], A_star=A_star, A_P=A_P, A_m=A_m, A_cls=A_cls,
                       B_pl=B_pl, B_P=B_P, B_m=B_m, B_cls=B_cls))
    return out


def ruwe_study(z, meta, jobs, d):
    """Who has RUWE > 1.4: stars with no giant, one giant, or several?"""
    pos = {j: k for k, j in enumerate(d["job"])}
    res = {}
    nul = [pos[j["job"]] for j in jobs if j["kind"] == "N" and j["job"] in pos]
    res["null"] = (d["ruwe_dr4"][nul], d["ruwe_dr3"][nul])
    for scen in SCENARIOS:
        npl = np.bincount(z[f"{scen}_star"], minlength=meta["n_stars"])
        snr_sum = np.bincount(z[f"{scen}_star"], weights=z[f"{scen}_snr"], minlength=meta["n_stars"])
        rows = [(pos[j["job"]], int(j["star"])) for j in jobs if j["kind"] == "A" and j["scen"] == scen
                and j["job"] in pos]
        k = np.array([r[0] for r in rows], int)
        s = np.array([r[1] for r in rows], int)
        res[scen] = dict(r4=d["ruwe_dr4"][k], r3=d["ruwe_dr3"][k], n_pl=npl[s], snr=snr_sum[s],
                         hosts_total={1: int((npl == 1).sum()), 2: int((npl >= 2).sum())})
    return res


def summary_lines(meta, res, ruwe, which_done):
    L = ["=== STEP 3 SUMMARY ===",
         f"real Gaia DR3 stars used: {meta['n_stars']} (random subset of {meta['n_selected']} selected); "
         + (f"mass from Gaia FLAME for {100 * meta['frac_flame_mass']:.0f}%, otherwise from absolute G"
            if meta['frac_flame_mass'] > 0 else "masses from absolute G with the main-sequence mass-magnitude "
            "relation (as Lammers & Winn)"),
         "selection: " + "; ".join(f"{k} {v}" for k, v in meta["selection"].items())]
    for scen, dd in meta["scen"].items():
        L.append(f"  {scen:<11} giants per star {dd['n_planets'] / meta['n_stars']:.3f}; stars with a giant "
                 f"{dd['hosts']}, of which with 2 or more {100 * dd['multi_hosts'] / max(dd['hosts'], 1):.0f}%")
    L += ["", "HOW TO READ: 'inferred / true' is the occurrence rate a DR4 user would infer (detections counted at",
          "their FITTED mass and period, corrected with the completeness measured the same way from single-planet",
          "injections, so planets scattering between bins are accounted for) divided by the true rate. 1.00 = right.",
          "A = real systems (all planets on each star). B = the same planets, each ALONE on its star (same data",
          "otherwise). A vs B isolates the effect of planets sharing a star. Reported planets are 'correct' (period",
          "within 10%, mass within 20% of a real planet), 'imprecise' (within x1.5 in period and x2 in mass: a real",
          "planet, poorly measured) or 'spurious' (no real planet resembles it). Power laws are fitted by forward",
          "modelling through the injection response (true -> fitted parameters).", ""]
    if meta.get("variant"):
        L.append(f"POPULATION VARIANT: {', '.join(meta['variant'])} (see src/gedi/population.py)")
        for scen, dd in meta["scen"].items():
            st = dd.get("stability") or {}
            if st:
                L.append(f"  [{scen}] N-body screening ({st['stable_orbits']:.0e} inner orbits): "
                         f"{st['n_unstable_first']} of {st['n_multi']} multi-planet systems unstable when first "
                         f"drawn; {st['n_reduced']} still unstable after redrawing (kept their first planet only)")
        L.append("")
    for a, r in res.items():
        L.append(f"--- {LABEL[a]} ---  (injections fitted: {r['I_done']}/{meta['n_I']})")
        for scen, s in r["scen"].items():
            rel, win = s["reliable"], s["window"]
            tot = lambda x: float(np.nansum(np.where(win, x, np.nan)))  # noqa: E731
            tT, tA, tB = tot(s["truth"]), tot(s["rA"].f), tot(s["rB"].f)
            nwA, nwB = int(np.sum(np.where(win, s["rA"].n_det, 0))), int(np.sum(np.where(win, s["rB"].n_det, 0)))
            L.append(f" [{scen}] planets reported: A {s['nA']}  B {s['nB']}  (A from {s['A_done']} fitted stars); "
                     f"{100 * s['frac_multi']:.0f}% of A detections are on stars with 2+ giants")
            cA, cB = s["clsA"], s["clsB"]
            pc = lambda c, n: f"{c} ({100 * c / max(n, 1):.0f}%)"  # noqa: E731
            L.append(f"   reported planets: correct / imprecise / spurious:  A {pc(cA[0], s['nA'])} / "
                     f"{pc(cA[1], s['nA'])} / {pc(cA[2], s['nA'])};  B {pc(cB[0], s['nB'])} / {pc(cB[1], s['nB'])} / "
                     f"{pc(cB[2], s['nB'])}")
            L.append(f"   real planets found (correct + imprecise): A vs B {(cA[0] + cA[1]) / max(cB[0] + cB[1], 1) - 1:+.1%};"
                     f"  spurious A vs B: {cA[2]} vs {cB[2]}")
            L.append(f"   spurious planets DR4 would publish as a normal orbit (significance >= 5): {s['n_wrong_conf']} "
                     f"of {cA[2]}")
            L.append(f"   pre-fit cut check: {s['val'][0]} of {s['val'][1]} systems just below the cut detected")
            if win.any() and tT > 0 and tB > 0:
                L.append(f"   occurrence summed over the {int(win.sum())} bins with completeness >= {MIN_ETA:g} "
                         f"(planets per star): true {tT:.4f}, inferred A {tA:.4f} ({tA / tT - 1:+.1%}), "
                         f"B {tB:.4f} ({tB / tT - 1:+.1%}); A vs B {tA / tB - 1:+.1%}   "
                         f"[from {nwA} / {nwB} detections{'; FEW - treat as rough' if min(nwA, nwB) < 50 else ''}]")
            else:
                L.append("   no totals yet (no detections in bins with measurable completeness)")
            L.append("   inferred / true per bin, A | B  ('~' = fewer than "
                     f"{MIN_DET} detections, unreliable):")
            L.append("      M_J \\ yr   " + "".join(f"{P_EDGES[k]:g}-{P_EDGES[k + 1]:g}".center(16)
                                                     for k in range(len(P_EDGES) - 1)))
            for i in range(len(M_EDGES) - 1):
                cells = []
                for j in range(len(P_EDGES) - 1):
                    if not np.isfinite(s["qA"][i, j]) and not np.isfinite(s["qB"][i, j]):
                        cells.append("-".center(16))
                        continue
                    mark = "" if rel[i, j] else "~"
                    cells.append(f"{mark}{s['qA'][i, j]:.2f}|{s['qB'][i, j]:.2f}".center(16))
                L.append(f"    {M_EDGES[i]:>4g}-{M_EDGES[i + 1]:<4g}   " + "".join(cells))
            for name, ff in (("true", s["fT"]), ("B", s["fB"]), ("A", s["fA"])):
                L.append(f"   power law {name:<4}: " + (ff.summary() if ff else "fit failed"))
        L.append("")
    if ruwe:
        n4, n3 = ruwe["null"]
        L += [f"--- RUWE: what pushes a star above {RUWE_CUT}? (simulated; real stars, our planets) ---",
              f"  no planet: RUWE > {RUWE_CUT} in {100 * np.mean(n4 > RUWE_CUT):.2f}% (DR4) / "
              f"{100 * np.nanmean(n3 > RUWE_CUT):.2f}% (DR3) of {len(n4)} stars (noise and bright-star jitter only)"]
        for scen in SCENARIOS:
            r = ruwe[scen]
            for lab, sel in (("1 giant", r["n_pl"] == 1), ("2+ giants", r["n_pl"] >= 2)):
                tot_hosts = r["hosts_total"][1 if lab == "1 giant" else 2]
                hi4, hi3 = (sel & (r["r4"] > RUWE_CUT)).sum(), (sel & (r["r3"] > RUWE_CUT)).sum()
                L.append(f"  [{scen}] stars with {lab}: {tot_hosts}; RUWE > {RUWE_CUT}: {hi4} in DR4 "
                         f"({100 * hi4 / max(tot_hosts, 1):.2f}% of such stars), {hi3} in DR3")
            hi = r["r4"] > RUWE_CUT
            if hi.any():
                L.append(f"  [{scen}] of all simulated stars with DR4 RUWE > {RUWE_CUT}: "
                         f"{100 * np.mean(r['n_pl'][hi] >= 2):.0f}% have 2+ giants")
        L.append("  (stellar and brown-dwarf companions are not simulated; in real data they dominate RUWE > 1.4)")
    return L


# ============================================================================ multiplicity scan
SCAN_F = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])   # fraction of giant hosts with 2+ giants
SCAN_METHODS = ["one", "best"]
SCAN_BOOT = 50
N_MBIN = 4                          # stellar-mass bins for the scan
LW_YIELD = (7500, 2100)            # Lammers & Winn (2026): predicted DR4 planet discoveries (all distances, no RUWE cut)
F_MARK = {"independent placement": None, "Bryan et al. (2016)": None}     # filled from the two simulated scenarios


def sky_fraction(meta) -> float:
    import re
    for k in meta.get("selection", {}):
        mm = re.search(r"([0-9.]+)% of the sky", k)
        if mm:
            return float(mm.group(1)) / 100
    return 1.0


def scan_multiplicity(res, meta, f_grid=SCAN_F, n_boot=SCAN_BOOT, seed=11):
    """Occurrence bias and yield vs multiplicity, by reweighting the fitted hosts (no new fits).
    Returns {method: dict of arrays over f (nominal, and bootstrap 16/84 percentiles)}."""
    from gedi import multiplicity as MU
    from gedi import rates as R

    n_star = meta["n_stars"]
    n_full = meta["n_selected"] / sky_fraction(meta)        # stars passing the cuts, whole sky
    m_win, P_win = (M_FINE[0], M_FINE[-1]), (P_FINE[0], P_FINE[-1])
    out = {}
    for an in [a for a in SCAN_METHODS if a in res]:
        r = res[an]
        hosts = {sc: r["scen"][sc]["hosts"] for sc in SCENARIOS}
        pool = MU.build_pool(hosts)
        h0 = hosts["independent"]
        n_hosts = len(pool.k)
        # stellar-mass bins: the CLS in-window rate is imposed per bin, so planets keep their M_star dependence
        ms_all = h0["m_star"]
        edges = np.quantile(ms_all, np.linspace(0, 1, N_MBIN + 1))
        edges[-1] += 1e-9
        star_bin = np.clip(np.digitize(ms_all, edges) - 1, 0, N_MBIN - 1)
        frac_bin = np.bincount(star_bin, minlength=N_MBIN) / n_star
        host_bin = np.clip(np.digitize(pool.m_star, edges) - 1, 0, N_MBIN - 1)
        win0 = MU.in_window(h0["pm"], h0["pP"], m_win, P_win)
        lam_bin = np.bincount(star_bin[h0["pst"][win0]], minlength=N_MBIN) / np.maximum(
            np.bincount(star_bin, minlength=N_MBIN), 1)         # CLS in-window giants per star, per mass bin
        r_cls = float(win0.sum()) / n_star
        t_win = np.bincount(pool.t_host, weights=MU.in_window(pool.t_m, pool.t_P, m_win, P_win), minlength=n_hosts)
        region = (r["eta"] >= MIN_ETA)                           # coarse bins a DR4 user can trust
        im, ip = np.digitize(pool.t_m, M_EDGES) - 1, np.digitize(pool.t_P, P_EDGES) - 1
        okc = (im >= 0) & (im < len(M_EDGES) - 1) & (ip >= 0) & (ip < len(P_EDGES) - 1)
        t_reg = np.zeros(n_hosts)
        sel = okc.copy()
        sel[okc] = region[im[okc], ip[okc]]
        np.add.at(t_reg, pool.t_host[sel], 1.0)
        inj = r["inj"]
        rng = np.random.default_rng(seed)
        reps = []                                                # (host multipliers, eta, R) per replicate
        for b in range(n_boot + 1):
            if b == 0:
                g, wi = np.ones(n_hosts), None
            else:
                g, wi = rng.poisson(1.0, n_hosts).astype(float), rng.poisson(1.0, len(inj["m_all"])).astype(float)
            eta, _, _ = eta_fitted(inj, wi)
            Rm = R.response_matrix(inj["m_all"], inj["P_all"], inj["t_m"], inj["t_P"], inj["f_m"], inj["f_P"],
                                   M_FINE, P_TRUE_EDGES, M_FINE, P_FINE, w_all=wi,
                                   w_det=None if wi is None else wi[inj["idx"]])
            reps.append((g, eta, Rm))
        keys = ["pl_A", "pl_B", "bin_A", "bin_B", "yield_A", "yield_B", "y_ok", "y_imp", "y_spu", "mu", "lost",
                "multi_det", "f_real", "bin_true"]
        res_f = {k: np.full((len(f_grid), n_boot + 1), np.nan) for k in keys}
        for i, f in enumerate(f_grid):
            v0, mu, lost, f_real = MU.host_weights_by_mass(pool, f, host_bin)
            kmax = int(pool.k.max())
            Nkj = np.zeros((N_MBIN, kmax + 1))
            np.add.at(Nkj, (host_bin, pool.k), 1.0)
            for b, (g, eta, Rm) in enumerate(reps):
                # re-normalise within each (mass bin, k) so the multiplicity distribution stays exact
                gkj = np.zeros((N_MBIN, kmax + 1))
                np.add.at(gkj, (host_bin, pool.k), g)
                with np.errstate(invalid="ignore", divide="ignore"):
                    v = v0 * g * np.where(gkj > 0, Nkj / gkj, 0.0)[host_bin, pool.k]
                # hosts per star in each mass bin so that the in-window rate is the CLS one there
                ew = np.bincount(host_bin, weights=v * t_win, minlength=N_MBIN)
                with np.errstate(invalid="ignore", divide="ignore"):
                    Hj = np.where(ew > 0, lam_bin / ew, 0.0) * frac_bin
                v = v * Hj[host_bin]
                Hs = 1.0
                wa, wb = Hs * v[pool.a_host], Hs * v[pool.b_host]
                with np.errstate(invalid="ignore", divide="ignore"):
                    def binned(mm, PP, ww):
                        c = R.counts_2d(mm, PP, M_EDGES, P_EDGES, weights=ww)
                        return float(np.nansum(np.where(region, c / eta, np.nan)))
                    tr = Hs * float(np.sum(v * t_reg))
                    res_f["bin_A"][i, b] = binned(pool.a_m, pool.a_P, wa) / tr
                    res_f["bin_B"][i, b] = binned(pool.b_m, pool.b_P, wb) / tr
                    res_f["bin_true"][i, b] = tr
                for tag, mm, PP, ww in (("A", pool.a_m, pool.a_P, wa), ("B", pool.b_m, pool.b_P, wb)):
                    try:
                        fit_ = R.fit_power_law_forward(mm, PP, n_star, Rm, M_FINE, P_TRUE_EDGES, M_FINE, P_FINE,
                                                       weights=ww * n_star, errors=False)
                        res_f[f"pl_{tag}"][i, b] = fit_.rate / r_cls
                    except Exception:
                        pass
                res_f["yield_A"][i, b] = n_full * wa.sum()
                res_f["yield_B"][i, b] = n_full * wb.sum()
                for c, k in enumerate(("y_ok", "y_imp", "y_spu")):
                    res_f[k][i, b] = n_full * wa[pool.a_cls == c].sum()
                res_f["multi_det"][i, b] = wa[pool.k[pool.a_host] >= 2].sum() / max(wa.sum(), 1e-30)
                res_f["mu"][i, b], res_f["lost"][i, b], res_f["f_real"][i, b] = mu, lost, f_real
        summ = {}
        for k, a in res_f.items():
            summ[k] = a[:, 0]
            summ[k + "_lo"], summ[k + "_hi"] = np.nanpercentile(a[:, 1:], [16, 84], axis=1)
        summ.update(f=np.asarray(f_grid), r_cls=r_cls, n_full=n_full,
                    f_scen={sc: float(np.mean(hosts[sc]["npl"][hosts[sc]["npl"] > 0] >= 2)) for sc in SCENARIOS})
        out[an] = summ
    return out


def write_scan_csv(scan, path):
    cols = ["method", "f_multi", "mu"] + [k for k in next(iter(scan.values())) if isinstance(
        next(iter(scan.values()))[k], np.ndarray) and k not in ("f", "mu", "mu_lo", "mu_hi")]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for an, s in scan.items():
            for i, f in enumerate(s["f"]):
                w.writerow([an, f"{f:.2f}", f"{s['mu'][i]:.4f}"] + [f"{s[k][i]:.5g}" for k in cols[3:]])


def scan_lines(scan):
    L = ["--- MULTIPLICITY SCAN (reweighting the fitted systems; same CLS occurrence rate at every level) ---",
         "f = fraction of giant-planet hosts with 2 or more giants. Occurrence: inferred / true in the window "
         "0.3-13 M_J, 0.1-8 yr",
         "(power law, forward-modelled) and in the bins with completeness >= 5% (binned). A = real systems, "
         "B = planets alone.",
         "Yield: planets reported by Gaia DR4 for stars within 250 pc and RUWE <= 1.4, scaled to the whole sky "
         f"(Lammers & Winn predict {LW_YIELD[0]:,} +/- {LW_YIELD[1]:,} for all distances).",
         "Ranges are bootstrap 16-84% intervals of THIS simulation."]
    for an, s in scan.items():
        fs = ", ".join(f"{k} {100 * v:.0f}%" for k, v in s["f_scen"].items())
        L.append(f" [{LABEL[an]}]  (simulated scenarios: {fs}; CLS in-window rate {s['r_cls']:.4f} per star)")
        L.append("     f    rate A (PL)          rate B (PL)    rate A (binned)  rate B (binned)   "
                 "yield A   correct/imprecise/spurious   on 2+ hosts")
        for i, f in enumerate(s["f"]):
            rng_ = lambda k: f"{s[k][i]:.2f} [{s[k + '_lo'][i]:.2f},{s[k + '_hi'][i]:.2f}]"  # noqa: E731
            L.append(f"   {s['f_real'][i]:4.2f}  {rng_('pl_A'):<19} {s['pl_B'][i]:.2f}            {rng_('bin_A'):<19} "
                     f"{s['bin_B'][i]:.2f}            {s['yield_A'][i]:7,.0f}   "
                     f"{s['y_ok'][i]:6,.0f}/{s['y_imp'][i]:6,.0f}/{s['y_spu'][i]:5,.0f}       "
                     f"{100 * s['multi_det'][i]:3.0f}%" + (f"   (k-tail {100 * s['lost'][i]:.1f}% dropped)"
                                                            if s["lost"][i] > 0.005 else ""))
    return L


def plot_scan(scan, png):
    plt, ps = _setup()
    fig, ax = plt.subplots(1, 3, figsize=(18, 5.4), constrained_layout=True)
    col = {"one": "#2a6fdb", "best": "#1a9e77"}
    for an, s in scan.items():
        f = 100 * s["f_real"]
        c = col[an]
        for a_, key, ttl in ((ax[0], "pl", "whole giant range 0.3-13 M_J, 0.1-8 yr (power-law fit)"),
                             (ax[1], "bin", "well-measured bins only (completeness >= 5%)")):
            true = 100 * (np.full(len(f), s["r_cls"]) if key == "pl" else s["bin_true"])
            if an == "one":
                a_.plot(f, true, color=ps.INK, lw=2, label="true rate (input)")
            a_.fill_between(f, true * s[f"{key}_A_lo"], true * s[f"{key}_A_hi"], color=c, alpha=0.18, lw=0)
            a_.plot(f, true * s[f"{key}_A"], "-o", color=c, ms=4, label=f"{LABEL[an]}: real systems")
            a_.plot(f, true * s[f"{key}_B"], "--", color=c, lw=1.2, label=f"{LABEL[an]}: planets alone")
            a_.set_title(ttl, fontsize=10.5)
        ax[2].fill_between(f, s["yield_A_lo"], s["yield_A_hi"], color=c, alpha=0.18, lw=0)
        ax[2].plot(f, s["yield_A"], "-o", color=c, ms=4, label=f"{LABEL[an]}: all reported")
        ax[2].plot(f, s["y_ok"] + s["y_imp"], ":", color=c, lw=1.5, label=f"{LABEL[an]}: real planets only")
    for a_ in ax[:2]:
        a_.set_ylabel("occurrence rate (planets per 100 stars)")
        a_.set_ylim(0, None)
    ax[2].axhspan(LW_YIELD[0] - LW_YIELD[1], LW_YIELD[0] + LW_YIELD[1], color="0.85", lw=0)
    ax[2].axhline(LW_YIELD[0], color="0.4", lw=1)
    ax[2].text(1, LW_YIELD[0] + LW_YIELD[1] * 0.1, "Lammers & Winn prediction\n(all distances)", fontsize=8, color="0.3")
    ax[2].set_ylabel("planets reported in DR4 (whole sky, ≤ 250 pc)")
    ax[2].set_title("Yield at the same (CLS) occurrence rate", fontsize=10.5)
    ax[2].set_ylim(0, None)
    any_s = next(iter(scan.values()))
    for a_ in ax:
        a_.set_xlabel("giant-planet hosts with a second giant (%)")
        for lab, fv, ls in (("independent", any_s["f_scen"]["independent"], ":"),
                            ("Bryan et al. 2016", any_s["f_scen"]["clustered"], "--")):
            a_.axvline(100 * fv, color="0.5", ls=ls, lw=1)
            a_.text(100 * fv + 0.8, a_.get_ylim()[1] - 0.02 * np.diff(a_.get_ylim())[0], lab, rotation=90,
                    fontsize=8, color="0.4", va="top")
        a_.legend(fontsize=7.5, loc="lower left")
    fig.suptitle("How giant-planet multiplicity changes what Gaia DR4 infers (same true occurrence rate at every "
                 "multiplicity; shaded: 16-84% bootstrap)", fontsize=12)
    fig.savefig(png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {png}")


# ============================================================================ figures
def _setup():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from gedi import plotstyle as ps
    ps.apply()
    return plt, ps


MCOL = ["#8fb8f0", "#2a6fdb", "#123a7a"]
TICK = {"one": "one\nplanet", "accel": "+ accel-\neration", "jerk": "+ jerk", "two": "two\nplanets",
        "best": "best\n(BIC)", "dr3": "DR3-\nstyle"}


def plot_rates(res, png):
    """Occurrence rate per bin (planets per 100 stars): true (black bars), real systems (A), planets alone (B)."""
    plt, ps = _setup()
    an = list(res)
    fig, axes = plt.subplots(len(SCENARIOS), len(an), figsize=(4.6 * len(an) + 1, 8.6), squeeze=False,
                             sharey=True, constrained_layout=True)
    x = np.sqrt(P_EDGES[1:] * P_EDGES[:-1])
    lo_y, hi_y = 0.02, 3.0
    for c, a in enumerate(an):
        for r, scen in enumerate(SCENARIOS):
            ax = axes[r, c]
            s = res[a]["scen"][scen]
            for k in range(len(M_EDGES) - 1):
                off = 1 + 0.07 * (k - 1)
                rel = s["reliable"][k]
                tr = 100 * s["truth"][k]
                for j in range(len(x)):
                    if rel[j] and tr[j] > 0:
                        ax.plot([x[j] * off / 1.12, x[j] * off * 1.12], [tr[j]] * 2, color=ps.INK, lw=2.2,
                                solid_capstyle="butt", label="true rate (input)" if (k, j) == (2, 2) else None)
                for rr, mfc, ls, tag in ((s["rA"], MCOL[k], "-", "A"), (s["rB"], "white", "--", "B")):
                    f, e = 100 * rr.f[k], 100 * rr.err[k]
                    ok = np.isfinite(f) & (f > 0)
                    good = ok & rel
                    if good.any():
                        ax.errorbar(x[good] * off, f[good], e[good], color=MCOL[k], marker="o", mfc=mfc, ls=ls,
                                    lw=1.5 if tag == "A" else 1.0, capsize=2,
                                    label=f"{M_EDGES[k]:g}-{M_EDGES[k + 1]:g} M_J, "
                                          + ("real systems (A)" if tag == "A" else "planets alone (B)"))
                    weak = ok & ~rel
                    if weak.any():
                        ax.scatter(x[weak] * off, np.clip(f[weak], lo_y, hi_y), s=22, facecolor="none",
                                   edgecolor=ps.C_NONE, zorder=3)
                nA = s["rA"].n_det[k]
                for j in range(len(x)):
                    if nA[j] > 0:
                        ax.text(x[j] * off, 0.02 + 0.045 * k, f"{int(nA[j])}", fontsize=6.5, ha="center",
                                color=MCOL[k], transform=ax.get_xaxis_transform())
            ax.set_yscale("log")
            ps.years_axis(ax, "x", lo=0.1, hi=8, label="orbital period")
            ps.mark_dr4(ax, text=(r == 0 and c == 0), y=0.75)
            ax.set_ylim(lo_y, hi_y)
            ax.set_title(f"{LABEL[a]}\n{scen} multiplicity", fontsize=10.5)
            if c == 0:
                ax.set_ylabel("occurrence rate\n(planets per 100 stars, in each bin)")
    h, l = axes[0, 0].get_legend_handles_labels()
    seen = {}
    for hh, ll in zip(h, l):
        seen.setdefault(ll, hh)
    fig.legend(list(seen.values()), list(seen.keys()), loc="lower center", ncol=4, fontsize=8.5,
               bbox_to_anchor=(0.5, -0.08))
    fig.suptitle("Giant-planet occurrence rate inferred from simulated DR4 data on real Gaia stars (≤ 250 pc); "
                 "black bars: true rate.\n"
                 f"Grey open circles: completeness < {100 * MIN_ETA:g}% or fewer than {MIN_DET} detections "
                 "(unreliable). Small numbers: detections in A.", fontsize=11.5)
    fig.savefig(png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {png}")


def plot_overview(res, png):
    plt, ps = _setup()
    an = list(res)
    fig, ax = plt.subplots(1, 3, figsize=(19, 5.6), constrained_layout=True)
    top = 1.0
    w = 0.38
    for k, scen in enumerate(SCENARIOS):
        for i, a in enumerate(an):
            s = res[a]["scen"][scen]
            win = s["window"]
            tT = np.nansum(np.where(win, s["truth"], np.nan))
            ax[k].plot([i - w, i + w], [100 * tT] * 2, color=ps.INK, lw=2.2, solid_capstyle="butt",
                       label="true rate (input)" if i == 0 else None)
            for dx, rr, tag in ((-w / 2, s["rA"], "A"), (w / 2, s["rB"], "B")):
                v = 100 * np.nansum(np.where(win, rr.f, np.nan))
                e = 100 * np.sqrt(np.nansum(np.where(win, rr.err, np.nan) ** 2))
                n = int(np.sum(np.where(win, rr.n_det, 0)))
                ax[k].bar(i + dx, v, w * 0.95, yerr=e, capsize=3, color=ps.C_TARGET if tag == "A" else "#bcd3f5",
                          label=("real systems (A)" if tag == "A" else "planets alone (B)") if i == 0 else None)
                ax[k].text(i + dx, v + e + 0.03, f"{v:.2f}", ha="center", fontsize=8)
                ax[k].text(i + dx, 0.06, f"n={n}", ha="center", va="bottom", rotation=90, fontsize=7,
                           color="white" if tag == "A" else ps.INK)
                top = max(top, v + e + 0.5)
        ax[k].set_xticks(range(len(an)), [TICK[a] for a in an], fontsize=8.5)
        ax[k].set_ylabel("occurrence rate, planets per 100 stars\n(3-13 M_J, 1-8 yr: the bins Gaia measures well)")
        ax[k].set_title(f"{scen} multiplicity")
        ax[k].legend(fontsize=8.5, loc="upper center", ncol=2)
    for k in range(len(SCENARIOS)):
        ax[k].set_ylim(0, top + 0.3)
    a = ax[2]
    for i, an_ in enumerate(an):
        for k, scen in enumerate(SCENARIOS):
            s = res[an_]["scen"][scen]
            n = max(s["nA"], 1)
            imp, spu = s["clsA"][1] / n, s["clsA"][2] / n
            x = i + (k - 0.5) * w
            a.bar(x, spu, w * 0.95, color=ps.C_WRONG, alpha=1.0 if k == 0 else 0.6,
                  label=f"spurious ({scen})" if i == 0 else None)
            a.bar(x, imp, w * 0.95, bottom=spu, color="#f3b98c", alpha=1.0 if k == 0 else 0.6, hatch="" if k == 0
                  else "//", edgecolor="white", label=f"imprecise ({scen})" if i == 0 else None)
    a.set_xticks(range(len(an)), [TICK[x] for x in an], fontsize=8.5)
    ps.percent_axis(a, "y", "share of reported planets (A)")
    a.set_title("Reported planets that are not accurate (A)\nimprecise = real planet, x1.5 in P / x2 in mass")
    a.set_ylim(0, 0.65)
    a.legend(fontsize=8, ncol=2, loc="upper center")
    fig.suptitle("Step 3 overview: does assuming at most one planet per star bias the inferred occurrence rate? "
                 "(black bars: true rate)",
                 fontsize=12)
    fig.savefig(png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {png}")


def plot_ruwe(ruwe, png):
    plt, ps = _setup()
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8), constrained_layout=True)
    bins = np.geomspace(0.8, 30, 50)
    n4, _ = ruwe["null"]
    r = ruwe["independent"]
    for vals, lab, col in ((n4, "no planet", ps.C_NONE), (r["r4"][r["n_pl"] == 1], "1 giant", ps.C_TARGET),
                           (r["r4"][r["n_pl"] >= 2], "2+ giants", ps.C_WRONG)):
        if len(vals):
            ax[0].hist(np.clip(vals, bins[0], bins[-1]), bins=bins, histtype="step", lw=1.8, color=col,
                       density=True, label=f"{lab} (n={len(vals)})")
    ax[0].set_xscale("log")
    ax[0].axvline(RUWE_CUT, color=ps.INK, ls="--", lw=1)
    ax[0].set_xlabel("simulated DR4 RUWE")
    ax[0].set_ylabel("share of stars (normalised)")
    ax[0].set_title("RUWE of stars that could host a detectable giant\n(independent multiplicity)")
    ax[0].legend(fontsize=9)
    sb = np.geomspace(0.8, 50, 12)
    xc = np.sqrt(sb[1:] * sb[:-1])
    for scen, ls in (("independent", "-"), ("clustered", "--")):
        r = ruwe[scen]
        for lab, sel, col in (("1 giant", r["n_pl"] == 1, ps.C_TARGET), ("2+ giants", r["n_pl"] >= 2, ps.C_WRONG)):
            fr = [np.mean(r["r4"][sel & (r["snr"] >= lo) & (r["snr"] < hi)] > RUWE_CUT)
                  if (sel & (r["snr"] >= lo) & (r["snr"] < hi)).sum() >= 10 else np.nan for lo, hi in zip(sb[:-1], sb[1:])]
            ax[1].plot(xc, fr, ls=ls, marker="o", ms=3.5, color=col, label=f"{lab}, {scen}")
    ax[1].set_xscale("log")
    ps.percent_axis(ax[1], "y", f"share with DR4 RUWE > {RUWE_CUT}")
    ax[1].set_xlabel("summed signal strength of the star's giants (SNR per transit)")
    ax[1].set_title(f"When do planets alone push RUWE above {RUWE_CUT}?")
    ax[1].legend(fontsize=8.5)
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
    for a_, ticks in ((ax[0], [1, 1.4, 2, 3, 5, 10, 20]), (ax[1], [1, 2, 3, 5, 10, 20, 50])):
        a_.xaxis.set_major_locator(FixedLocator(ticks))
        a_.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        a_.xaxis.set_minor_formatter(NullFormatter())
    fig.savefig(png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {png}")


def analyse(root: Path, which):
    P = paths(root)
    if not P["fits"].exists():
        raise SystemExit("no fits yet")
    meta = json.loads(P["meta"].read_text())
    with np.load(P["pop"]) as f:
        z = {k: f[k] for k in f.files}
    with open(P["jobs"]) as f:
        jobs = list(csv.DictReader(f))
    d = load_fits(root)
    judged = [a for a in JUDGED if a in which or (a == "best" and {"accel", "jerk", "two"} <= set(which))]
    if "dr3" in which:
        judged.append("dr3")
    res = {a: analyse_one(z, meta, jobs, d, a) for a in judged}
    ruwe = ruwe_study(z, meta, jobs, d) if "ruwe_dr4" in d else None
    L = summary_lines(meta, res, ruwe, which)
    scan = None
    if all(res[a]["I_done"] > 0 for a in SCAN_METHODS if a in res):
        print("multiplicity scan (reweighting, bootstrap) ...", flush=True)
        scan = scan_multiplicity(res, meta)
        L += [""] + scan_lines(scan)
        write_scan_csv(scan, P["scan_csv"])
    print("\n".join(L))
    P["summary"].write_text("\n".join(L) + "\n")
    plot_rates({a: res[a] for a in ("one", "best") if a in res} or res, P["rates_png"])
    plot_rates(res, root / "step3_rates_all.png")
    plot_overview(res, P["overview_png"])
    if ruwe:
        plot_ruwe(ruwe, P["ruwe_png"])
    if scan:
        plot_scan(scan, P["scan_png"])


# ============================================================================ main
def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    p.add_argument("--stars", type=int, default=0,
                   help="random subset of the selected Gaia stars (0 = all of them, the default)")
    p.add_argument("--inj-fits", type=int, default=15_000, help="single-planet injections to fit (completeness)")
    p.add_argument("--ruwe-max", type=float, default=1.4, help="drop real stars above this RUWE (0 = no cut)")
    p.add_argument("--with-dr3-pipeline", action="store_true", help="also run the (slow) DR3-style cascade")
    p.add_argument("--seed", type=int, default=3)
    p.add_argument("-j", "--jobs", type=int, default=-1)
    p.add_argument("--analyse", action="store_true", help="only redo summary and figures")
    p.add_argument("--download-only", action="store_true", help="only fetch the Gaia stars (no fitting)")
    p.add_argument("--build-only", action="store_true",
                   help="download if needed and build stars, planets and the job list, then stop "
                        "(before starting --part jobs on a cluster)")
    p.add_argument("--rebuild", action="store_true", help="redraw stars and planets (deletes old fits)")
    p.add_argument("--source", choices=["auto", "vizier", "esa"], default="auto",
                   help="archive for the Gaia download (auto: VizieR first, ESA if a piece stalls)")
    p.add_argument("--sky-fraction", type=float, default=None,
                   help="fraction of the sky to download (default: 1 = whole sky; 0.25 with --circles)")
    p.add_argument("--circles", action="store_true",
                   help="use the 3-deg sky circles of versions <= 0.2.0 (max 0.52 of the sky) "
                        "and results/step3_gaia/")
    p.add_argument("--variant", default="",
                   help="population variants, comma-separated: coplanar, stable, lowecc "
                        "(results go to results/step3_allsky_<variants> unless --out is given)")
    p.add_argument("--out", type=Path, default=None,
                   help="results folder (default results/step3_allsky/, or results/step3_gaia/ with --circles)")
    p.add_argument("--download-jobs", type=int, default=6, help="sky circles requested at the same time")
    p.add_argument("--part", default=None, help="k/n: fit only part k of n (0-based), e.g. one cluster job each")
    p.add_argument("--merge-parts", action="store_true", help="join the --part results into fits.csv")
    a = p.parse_args()
    layout = "circles" if a.circles else "tiles"
    if a.sky_fraction is None:
        a.sky_fraction = 0.25 if a.circles else 1.0
    sky = (layout, a.sky_fraction)
    variant = tuple(sorted(v.strip() for v in a.variant.split(",") if v.strip()))
    default_root = RESULTS if a.circles else RESULTS_ALLSKY
    if variant:
        default_root = default_root.with_name(default_root.name + "_" + "_".join(variant))
    root = a.out or default_root
    ruwe_max = None if a.ruwe_max <= 0 else a.ruwe_max
    part = tuple(int(x) for x in a.part.split("/")) if a.part else None
    if a.merge_parts:
        merge_parts(root)
        analyse(root, FITTED + (["dr3"] if a.with_dr3_pipeline else []))
        return
    which = FITTED + (["dr3"] if a.with_dr3_pipeline else [])
    if not a.analyse:
        from gedi import gaia_stars as GS
        if not GS.have_stars(layout, a.sky_fraction):
            print(f"1. downloading Gaia DR3 stars (once; {100 * a.sky_fraction:.0f}% of the sky in "
                  f"{layout}, resumable) ...", flush=True)
            GS.download(sky_fraction=a.sky_fraction, source=a.source, jobs=a.download_jobs, layout=layout)
        else:
            print(f"1. Gaia DR3 stars already downloaded ({GS.merged_path(layout)})")
        if a.download_only:
            import numpy as _np
            sel, counts = GS.select(dict(_np.load(GS.merged_path(layout))), ruwe_max)
            for k, v in counts.items():
                print(f"    {k}: {v}")
            return
        if a.rebuild or not built_matches(root, a.stars, a.inj_fits, a.seed, ruwe_max, sky, variant):
            if paths(root)["fits"].exists():
                old = paths(root)["fits"].with_suffix(f".old{int(time.time())}.csv")
                paths(root)["fits"].rename(old)
                print(f"   settings changed: old fits moved to {old.name}")
            print("2. building stars, planets and the job list ...", flush=True)
            build(root, a.stars, a.inj_fits, a.seed, ruwe_max, sky=sky, variant=variant, n_jobs=a.jobs)
        else:
            print("2. stars and planets already built with these settings")
        if a.build_only:
            return
        print("3. fitting ..." + (f" (part {part[0]} of {part[1]})" if part else ""), flush=True)
        fit(root, which, a.jobs, part=part)
        if part:
            print(f"part {part[0]}/{part[1]} done; when all parts are done, copy the fits.part*.csv files into "
                  f"{root} and run: python scripts/step3.py --merge-parts")
            return
    print("4. analysing ...", flush=True)
    analyse(root, which)


if __name__ == "__main__":
    main()
