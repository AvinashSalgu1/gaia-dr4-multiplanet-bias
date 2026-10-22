"""Census of known inner-planet hosts and their DR4 sensitivity to outer companions.

Runs on any machine with gaiamock (for each host's real Gaia scanning pattern) and internet
access (NASA Exoplanet Archive). Resumable; a few hours on a 12-core machine with the defaults:

    python scripts/extended/g2_census.py                  # download + injections + summary
    python scripts/extended/g2_census.py --analyse        # only redo summary and figures
    python scripts/extended/g2_census.py --quick          # 40 hosts, 2 injections per cell (test)

1. HOSTS. All planets in the NASA Exoplanet Archive composite table (pscomppars) are
   downloaded once (data/extended/exoplanet_archive.csv; --refresh to update). Hosts with
   at least one planet with P < 200 d are classed by their inner planets:
   hot Jupiter (P < 10 d, mass > 0.3 M_J) > warm Jupiter (10-200 d, > 0.3 M_J) >
   small planet (mass < 30 M_earth or radius < 4 R_earth) > other. Hosts beyond
   --max-dist pc or without distance, stellar mass or Gaia G are listed but not simulated.
   Hosts that already have a known outer giant (P > 200 d, > 0.3 M_J) are flagged.
2. SENSITIVITY. For each host, outer companions are injected into simulated DR4 epoch
   astrometry with the host's real scanning pattern, brightness, distance and mass, on a
   grid of masses (0.3, 1, 3, 10 M_J) and periods (0.5-12 yr), --n-inj per grid cell with
   random orientation, phase and Kipping eccentricity, and fitted with the one-orbit model
   (Lammers & Winn detection criteria). Cells whose signal-to-noise per transit is below
   0.4 count as non-detections without a fit (none is detected there in step 1).
3. FORECAST. Expected numbers of detected outer giants per class, integrating the
   detection probability over the field occurrence model (CLS-normalised, as in step 3),
   with Poisson 68% and 95% intervals, for the field rate and for 3 times the field rate;
   the smallest rate ratio each class could distinguish from the field at 3 sigma; and the
   number of hosts that pass detection-probability thresholds at 3, 1 and 0.3 M_J (P = 3 yr).

Outputs in results/extended/census/: hosts.csv (one row per host, with its detection
probabilities), summary.txt, census_sensitivity.png, census_forecast.png.
Assumptions: the known inner planet's own astrometric signal is ignored (micro-arcseconds);
stellar companions are not simulated; the field occurrence model is used for every class.
"""
from __future__ import annotations

import argparse
import csv
import io
import urllib.parse
import urllib.request
import zlib

import numpy as np

import common as C

TAP = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
COLS = ("pl_name,hostname,ra,dec,sy_dist,st_mass,sy_gaiamag,sy_vmag,pl_orbper,pl_bmassj,pl_bmasse,pl_rade,"
        "discoverymethod")
M_GRID = np.array([0.3, 1.0, 3.0, 10.0])
P_GRID = np.array([0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0])
SNR_SKIP = 0.4
CLASSES = ["hot Jupiter", "warm Jupiter", "small planet", "other"]
ARCHIVE = C.DATA / "exoplanet_archive.csv"


# ---------------------------------------------------------------- hosts
def download(refresh=False):
    if ARCHIVE.exists() and not refresh:
        return
    q = f"select {COLS} from pscomppars"
    url = TAP + "?" + urllib.parse.urlencode(dict(query=q, format="csv"))
    print("downloading the NASA Exoplanet Archive composite table ...", flush=True)
    with urllib.request.urlopen(url, timeout=300) as r:
        text = r.read().decode()
    rows = list(csv.DictReader(io.StringIO(text)))
    if not rows:
        raise SystemExit("empty answer from the archive")
    C.DATA.mkdir(parents=True, exist_ok=True)
    ARCHIVE.write_text(text)
    print(f"saved {ARCHIVE} ({len(rows)} planets)")


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


def host_table(max_dist):
    with open(ARCHIVE) as f:
        planets = list(csv.DictReader(f))
    by_host = {}
    for p in planets:
        by_host.setdefault(p["hostname"], []).append(p)
    hosts = []
    for name, pl in by_host.items():
        P = np.array([_f(p["pl_orbper"]) for p in pl])
        mj = np.array([_f(p["pl_bmassj"]) for p in pl])
        me = np.array([_f(p["pl_bmasse"]) for p in pl])
        re_ = np.array([_f(p["pl_rade"]) for p in pl])
        inner = P < 200
        if not inner.any():
            continue
        giant = mj > 0.3
        small = (me < 30) | (re_ < 4)
        if np.any(inner & (P < 10) & giant):
            cls = "hot Jupiter"
        elif np.any(inner & (P >= 10) & giant):
            cls = "warm Jupiter"
        elif np.any(inner & small):
            cls = "small planet"
        else:
            cls = "other"
        p0 = pl[0]
        h = dict(host=name, cls=cls, ra=_f(p0["ra"]), dec=_f(p0["dec"]), dist=_f(p0["sy_dist"]),
                 m_star=_f(p0["st_mass"]), G=_f(p0["sy_gaiamag"]), V=_f(p0["sy_vmag"]),
                 n_planets=len(pl), known_outer_giant=int(np.any((P > 200) & giant)),
                 methods=";".join(sorted({p["discoverymethod"] for p, inn in zip(pl, inner) if inn})))
        h["simulate"] = int(np.isfinite([h["ra"], h["dec"], h["dist"], h["m_star"], h["G"]]).all()
                            and h["dist"] <= max_dist and 0.08 < h["m_star"] < 3)
        hosts.append(h)
    return hosts


# ---------------------------------------------------------------- sensitivity
def host_epochs(h, synthetic, rng):
    from gaiabias import epochs as E
    if synthetic:
        return E.synthetic(n_transits=80, rng=rng)
    return E.from_gaiamock(h["ra"], h["dec"], "dr4", reject_frac=0.0, rng=rng)


def sensitivity(h, n_inj, seed, synthetic=False):
    """Detections per (mass, period) cell for one host: array (len(M_GRID), len(P_GRID))."""
    from gaiabias.epochs import Epochs
    from gaiabias.fitting import detected_lw, fit_orbit_trend
    from gaiabias.orbits import Companion, draw_ecc_kipping, reflex_a0_mas
    from gaiabias.simulate import Star, sigma_al_per_transit, simulate

    rng = np.random.default_rng([seed, zlib.crc32(h["host"].encode())])
    ep_all = host_epochs(h, synthetic, rng)
    plx = 1000.0 / h["dist"]
    sig = float(sigma_al_per_transit(h["G"]))
    star = Star(parallax_mas=plx, m_star_msun=h["m_star"], G_mag=h["G"])
    det = np.zeros((len(M_GRID), len(P_GRID)), int)
    for i, m in enumerate(M_GRID):
        for j, Py in enumerate(P_GRID):
            P = 365.25 * Py
            if reflex_a0_mas(h["m_star"], m, P, plx) / sig < SNR_SKIP:
                continue
            for _ in range(n_inj):
                keep = rng.uniform(0, 1, len(ep_all)) > 0.1                 # 10% of transits lost, as gaiamock
                ep = Epochs(ep_all.t_yr[keep], ep_all.psi[keep], ep_all.plx_factor[keep])
                c = Companion(m, P, float(draw_ecc_kipping(rng)), np.arccos(rng.uniform(-1, 1)),
                              rng.uniform(0, np.pi), rng.uniform(0, 2 * np.pi), rng.uniform(0, P))
                sim = simulate(ep, star, [c], add_noise=False)
                noise = sig * rng.standard_normal(len(ep))
                if h["G"] < 13:
                    noise += rng.uniform(0, 0.04) * rng.standard_normal(len(ep))
                try:
                    f = fit_orbit_trend(ep, sim.obs + noise, sim.err, h["m_star"], n_trend=0)
                    det[i, j] += int(detected_lw(f.delta_chi2, f.m_mjup, f.P_days, f.sigma_P))
                except Exception:
                    pass
    return det


def run(hosts, n_inj, jobs, seed, synthetic, out):
    from joblib import Parallel, delayed
    path = out / "sensitivity.csv"
    done = set()
    if path.exists():
        with open(path) as f:
            done = {r["host"] for r in csv.DictReader(f)}
    todo = [h for h in hosts if h["simulate"] and h["host"] not in done]
    print(f"sensitivity: {len(done)} hosts done, {len(todo)} to go ({n_inj} injections per grid cell)", flush=True)
    cols = ["host", "n_inj"] + [f"det_{m:g}MJ_{p:g}yr" for m in M_GRID for p in P_GRID]
    if not path.exists():
        with open(path, "w", newline="") as f:
            csv.writer(f).writerow(cols)
    t = C.Timer()
    chunk = 24
    for k in range(0, len(todo), chunk):
        part = todo[k:k + chunk]
        res = Parallel(n_jobs=jobs)(delayed(sensitivity)(h, n_inj, seed, synthetic) for h in part)
        with open(path, "a", newline="") as f:
            w = csv.writer(f)
            for h, d in zip(part, res):
                w.writerow([h["host"], n_inj] + d.ravel().tolist())
        print(f"  {k + len(part)}/{len(todo)} hosts ({t})", flush=True)


# ---------------------------------------------------------------- forecast
def p_det_table(out):
    tab = {}
    with open(out / "sensitivity.csv") as f:
        for r in csv.DictReader(f):
            n = int(r["n_inj"])
            d = np.array([float(r[f"det_{m:g}MJ_{p:g}yr"]) for m in M_GRID for p in P_GRID]).reshape(len(M_GRID), -1)
            tab[r["host"]] = d / n
    return tab


def expected(pdet, m_star, m_range=(0.3, 13.0), n_m=40, n_p=60):
    """Expected detected outer giants for one host at the field rate (step 3 occurrence model)."""
    from gaiabias import population as POP
    lm = np.linspace(np.log(m_range[0]), np.log(m_range[1]), n_m)
    lp = np.linspace(np.log(0.3), np.log(20.0), n_p)
    dlm, dlp = lm[1] - lm[0], lp[1] - lp[0]
    a = (np.exp(lp)) ** (2 / 3) * m_star ** (1 / 3)
    dens = POP.LW_C * (m_star / POP.M_STAR_NORM) * POP.a_shape(a) * (2 / 3)   # per ln m per ln P
    # bilinear interpolation of the detection probability in (ln m, ln P); zero beyond the grid in P
    gm, gp = np.log(M_GRID), np.log(P_GRID)
    pm = np.array([[np.interp(x, gp, pdet[i], left=pdet[i][0], right=0.0) for x in lp] for i in range(len(gm))])
    full = np.array([[np.interp(y, gm, pm[:, k], left=0.0, right=pm[-1, k]) for k in range(len(lp))] for y in lm])
    return float(np.sum(full * dens[None, :]) * dlm * dlp)


def analyse(hosts, out):
    from scipy.stats import poisson
    tab = p_det_table(out)
    iP3 = list(P_GRID).index(3.0)
    rows = []
    for h in hosts:
        p = tab.get(h["host"])
        if p is None:
            h.update(p3_3MJ=np.nan, p3_1MJ=np.nan, p3_03MJ=np.nan, mu=np.nan, mu_low=np.nan)
        else:
            h.update(p3_3MJ=p[2, iP3], p3_1MJ=p[1, iP3], p3_03MJ=p[0, iP3], mu=expected(p, h["m_star"]),
                     mu_low=expected(p, h["m_star"], m_range=(0.3, 1.0)))
        rows.append(h)
    keys = ["host", "cls", "dist", "m_star", "G", "V", "n_planets", "known_outer_giant", "methods", "simulate",
            "p3_3MJ", "p3_1MJ", "p3_03MJ", "mu", "mu_low"]
    with open(out / "hosts.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(keys)
        for h in rows:
            w.writerow([f"{h[k]:.4g}" if isinstance(h[k], float) else h[k] for k in keys])
    L = ["=== CENSUS OF KNOWN INNER-PLANET HOSTS AND DR4 SENSITIVITY ===",
         f"archive file: {ARCHIVE.name}; hosts with a planet at P < 200 d: {len(rows)}; simulated: "
         f"{sum(h['simulate'] for h in rows)} (distance, mass and G known, within the distance limit).",
         "p(3 yr) = DR4 detection probability of an outer companion with a 3-yr period (one-orbit fit).",
         "Expected detections integrate the detection probability over the field occurrence model "
         "(0.3-13 M_J, periods 0.3-20 yr); ranges are Poisson 68% / 95%.", ""]
    L.append(f"{'class':<14}{'hosts':>7}{'sim.':>6}{'outer giant':>12}{'p(3yr)>0.5 at 3/1/0.3 M_J':>28}"
             f"{'expected (field)':>18}{'68%':>10}{'95%':>12}{'x3 field':>10}{'min ratio 3sig':>16}{'0.3-1 M_J':>11}")
    forecast = {}
    for c in CLASSES:
        hs = [h for h in rows if h["cls"] == c]
        sim = [h for h in hs if h["simulate"] and np.isfinite(h["mu"])]
        n3 = sum(h["p3_3MJ"] >= 0.5 for h in sim)
        n1 = sum(h["p3_1MJ"] >= 0.5 for h in sim)
        n03 = sum(h["p3_03MJ"] >= 0.5 for h in sim)
        mu = float(np.nansum([h["mu"] for h in sim]))
        mul = float(np.nansum([h["mu_low"] for h in sim]))
        i68 = poisson.ppf([0.16, 0.84], mu) if mu > 0 else (0, 0)
        i95 = poisson.ppf([0.025, 0.975], mu) if mu > 0 else (0, 0)
        kmin = 1 + 3 / np.sqrt(mu) if mu > 0 else np.inf
        forecast[c] = dict(mu=mu, i68=i68, i95=i95, n=len(sim))
        L.append(f"{c:<14}{len(hs):>7}{len(sim):>6}{sum(h['known_outer_giant'] for h in hs):>12}"
                 f"{f'{n3} / {n1} / {n03}':>28}{mu:>18.1f}{f'{i68[0]:.0f}-{i68[1]:.0f}':>10}"
                 f"{f'{i95[0]:.0f}-{i95[1]:.0f}':>12}{3 * mu:>10.1f}{kmin:>16.2f}{mul:>11.2f}")
    L += ["", "'min ratio 3sig': the smallest ratio of a class's outer-giant rate to the field rate that its expected",
          "number of detections could distinguish at 3 sigma (Poisson only). '0.3-1 M_J': expected detections of",
          "Saturn-to-Jupiter-mass companions at the field rate (the H3 regime).",
          "Field occurrence model for every class; known outer giants are not removed (they would be among the",
          "detections); stellar companions not simulated."]
    C.write_summary(out / "summary.txt", L)

    plt, ps = C.setup_plots()
    col = {"hot Jupiter": ps.C_WRONG, "warm Jupiter": "#7a5195", "small planet": ps.C_TARGET, "other": ps.C_NONE}
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.8), constrained_layout=True)
    for k, (key, lab) in enumerate((("p3_3MJ", "3 M_J"), ("p3_1MJ", "1 M_J"), ("p3_03MJ", "0.3 M_J"))):
        for c in CLASSES[:3]:
            sim = [h for h in rows if h["cls"] == c and h["simulate"] and np.isfinite(h[key])]
            ax[k].scatter([h["dist"] for h in sim], [h[key] for h in sim], s=9, alpha=0.6, color=col[c],
                          label=f"{c} ({len(sim)})")
        ax[k].set_xscale("log")
        ax[k].axhline(0.5, color="0.6", lw=0.8, ls=":")
        ax[k].set_xlabel("distance (pc)")
        ps.percent_axis(ax[k], "y", f"DR4 detection probability, {lab}, P = 3 yr")
        ax[k].set_ylim(-0.02, 1.02)
        ax[k].set_title(f"Outer companion of {lab} with a 3-yr period", fontsize=10.5)
    ax[0].legend(fontsize=8)
    fig.savefig(out / "census_sensitivity.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 4.4), constrained_layout=True)
    for i, c in enumerate(CLASSES[:3]):
        fc = forecast[c]
        for k, fac in enumerate((1, 3)):
            mu = fac * fc["mu"]
            lo, hi = poisson.ppf([0.025, 0.975], mu) if mu > 0 else (0, 0)
            l6, h6 = poisson.ppf([0.16, 0.84], mu) if mu > 0 else (0, 0)
            x = i + (k - 0.5) * 0.3
            ax.plot([x, x], [lo, hi], color=col[c], lw=1.2)
            ax.plot([x, x], [l6, h6], color=col[c], lw=5, alpha=0.6)
            ax.scatter([x], [mu], color=col[c], marker="o" if fac == 1 else "s", zorder=3,
                       label=("field rate" if fac == 1 else "3 x field rate") if i == 0 else None)
    ax.set_xticks(range(3), [f"{c}\n({forecast[c]['n']} hosts)" for c in CLASSES[:3]])
    ax.set_ylabel("expected outer giants detected by DR4")
    ax.set_yscale("symlog", linthresh=1)
    ax.legend(fontsize=8.5)
    ax.set_title("Forecast detections (68% thick, 95% thin)", fontsize=10.5)
    fig.savefig(out / "census_forecast.png", dpi=150, bbox_inches="tight")
    print(f"saved {out / 'census_sensitivity.png'} and {out / 'census_forecast.png'}")


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    p.add_argument("--n-inj", type=int, default=6, help="injections per (mass, period) grid cell")
    p.add_argument("--max-dist", type=float, default=200.0, help="simulate hosts within this distance (pc)")
    p.add_argument("-j", "--jobs", type=int, default=-1)
    p.add_argument("--seed", type=int, default=21)
    p.add_argument("--refresh", action="store_true", help="download the archive table again")
    p.add_argument("--analyse", action="store_true")
    p.add_argument("--quick", action="store_true", help="40 hosts, 2 injections per cell (test run)")
    p.add_argument("--synthetic", action="store_true", help="synthetic epochs instead of gaiamock (code test only)")
    p.add_argument("--archive", type=C.Path, default=None, help="use this archive CSV instead of downloading")
    a = p.parse_args()
    global ARCHIVE
    out = C.out_dir("census_quick" if a.quick else "census")
    if a.archive:
        ARCHIVE = a.archive
    else:
        download(a.refresh)
    hosts = host_table(a.max_dist)
    if a.quick:
        rng = np.random.default_rng(0)
        sim = [h for h in hosts if h["simulate"]]
        keep = {h["host"] for h in rng.choice(sim, min(40, len(sim)), replace=False)}
        hosts = [h for h in hosts if h["host"] in keep]
        a.n_inj = 2
    if not a.analyse:
        run(hosts, a.n_inj, a.jobs, a.seed, a.synthetic, out)
    analyse(hosts, out)


if __name__ == "__main__":
    main()
