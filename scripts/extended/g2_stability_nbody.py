"""Are the simulated multi-planet systems dynamically stable? (N-body check)

Needs the `rebound` package (pip install rebound); no gaiamock. About 10-40 minutes:

    python scripts/extended/g2_stability_nbody.py        # [--orbits 1e5] [-j -1]

The simulations keep a planet pair only if the gap between the inner planet's apocentre
and the outer planet's pericentre exceeds 2 sqrt(3) mutual Hill radii (the circular-orbit
criterion of Gladman 1993 applied at closest approach), ignoring mutual inclination. This
script integrates 600 pairs drawn from the step 3 population (data/extended/
pairs_for_nbody.csv, written by g1_export_pairs.py): 300 just above that limit, 200 near a
2:1 period ratio and 100 well separated, with their full 3D orbits (the simulated
inclinations and nodes, so mutual inclinations are included).

Each pair is integrated with WHFast for --orbits orbits of the inner planet (time step 1/30
of the inner period) with the MEGNO chaos indicator. A pair is UNSTABLE if the planets come
closer than one mutual Hill radius, one is ejected (beyond 100 x the outer semi-major
axis) or MEGNO exceeds 4 (strongly chaotic); CHAOTIC if MEGNO is between 2.5 and 4.

Output in results/extended/stability_nbody/: pairs.csv (result per pair), summary.txt,
stability_nbody.png.

With --all, every neighbouring pair of the step 3 population (data/extended/pairs_all.csv, about
38,000 pairs, written by g1_export_pairs.py --all) is integrated instead, in resumable chunks
(about 2 hours on 12 cores). Output in results/extended/stability_all/: pairs.csv (scenario,
star, status per pair) and summary.txt. g1_stable_scan.py then repeats the multiplicity scan
with only the systems whose pairs all survive.
"""
from __future__ import annotations

import argparse
import csv

import numpy as np

import common as C

MJUP_MSUN = 9.547919e-4


def integrate(row, n_orbits):
    import rebound
    ms = float(row["m_star"])
    sim = rebound.Simulation()
    sim.units = ("yr", "AU", "Msun")
    sim.add(m=ms)
    for k in ("1", "2"):
        P_yr = float(row[f"P{k}"]) / 365.25
        Tp = float(row[f"Tp{k}"]) / 365.25
        sim.add(m=float(row[f"m{k}"]) * MJUP_MSUN, P=P_yr, e=float(row[f"e{k}"]), inc=float(row[f"inc{k}"]),
                Omega=float(row[f"Omega{k}"]), omega=float(row[f"omega{k}"]) + np.pi,   # planet = star's omega + pi
                M=(-2 * np.pi * Tp / P_yr) % (2 * np.pi))
    sim.move_to_com()
    p1, p2 = sim.particles[1], sim.particles[2]
    r_hill = ((p1.m + p2.m) / (3 * ms)) ** (1 / 3) * 0.5 * (p1.a + p2.a)
    sim.integrator = "whfast"
    sim.dt = p1.P / 30.0
    sim.exit_max_distance = 100 * p2.a
    sim.exit_min_distance = r_hill
    sim.init_megno()
    t_end = n_orbits * p1.P
    status, megno = "stable", np.nan
    try:
        sim.integrate(t_end, exact_finish_time=0)
        megno = sim.megno()
        if megno > 4:
            status = "unstable"
        elif megno > 2.5:
            status = "chaotic"
    except rebound.Escape:
        status = "unstable"
    except rebound.Encounter:
        status = "unstable"
    except Exception:
        status = "failed"
    return status, megno, sim.t / p1.P


def run_all(n_orbits, jobs, quick=False):
    """All neighbouring pairs, resumable; one row per pair with scenario, star and status."""
    from joblib import Parallel, delayed
    from gaiabias.orbits import mutual_inclination
    out = C.out_dir("stability_all_quick" if quick else "stability_all")
    pairs = list(csv.DictReader(open(C.DATA / "pairs_all.csv")))
    for i, r in enumerate(pairs):
        r["pair"] = str(i)
    if quick:
        pairs = pairs[::200]
    path = out / "pairs.csv"
    cols = ["pair", "scenario", "star", "hill_ratio", "period_ratio", "e1", "e2", "mutual_inc_deg", "status", "megno",
            "orbits_done"]
    done = set()
    if path.exists():
        done = {r["pair"] for r in csv.DictReader(open(path))}
    else:
        with open(path, "w", newline="") as f:
            csv.writer(f).writerow(cols)
    todo = [r for r in pairs if r["pair"] not in done]
    print(f"{len(done)} pairs done, {len(todo)} to go", flush=True)
    t = C.Timer()
    chunk = 1200
    for i in range(0, len(todo), chunk):
        part = todo[i:i + chunk]
        res = Parallel(n_jobs=jobs)(delayed(integrate)(r, n_orbits) for r in part)
        with open(path, "a", newline="") as f:
            w = csv.writer(f)
            for r, (st, mg, nd) in zip(part, res):
                mi = np.degrees(mutual_inclination(float(r["inc1"]), float(r["Omega1"]), float(r["inc2"]),
                                                   float(r["Omega2"])))
                w.writerow([r["pair"], r["scenario"], r["star"], r["hill_ratio"], r["period_ratio"], r["e1"], r["e2"],
                            f"{mi:.1f}", st, f"{mg:.3f}", f"{nd:.0f}"])
        print(f"  {min(i + chunk, len(todo))}/{len(todo)} ({t})", flush=True)
    rows = list(csv.DictReader(open(path)))
    L = ["=== N-BODY CHECK OF EVERY SIMULATED GIANT PAIR ===",
         f"{len(rows)} neighbouring pairs of the step 3 population, {n_orbits:g} inner orbits (WHFast, MEGNO).",
         "unstable = close encounter (< 1 mutual Hill radius), ejection or MEGNO > 4; chaotic = MEGNO 2.5-4.", ""]
    for scen in ("independent", "clustered"):
        sel = [r for r in rows if r["scenario"] == scen]
        if not sel:
            continue
        un = np.array([r["status"] == "unstable" for r in sel])
        hosts = {}
        for r, u in zip(sel, un):
            hosts[r["star"]] = hosts.get(r["star"], False) or bool(u)
        mi = np.array([float(r["mutual_inc_deg"]) for r in sel])
        L.append(f"  {scen}: {len(sel)} pairs, unstable {100 * un.mean():.0f}%; multi-giant systems with an unstable "
                 f"pair: {sum(hosts.values())} of {len(hosts)} ({100 * np.mean(list(hosts.values())):.0f}%)")
        for lo, hi in ((0, 30), (30, 60), (60, 180)):
            m = (mi >= lo) & (mi < hi)
            if m.sum():
                L.append(f"      mutual inclination {lo}-{hi} deg: unstable {100 * un[m].mean():.0f}% of {m.sum()}")
    C.write_summary(out / "summary.txt", L)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--all", action="store_true", help="integrate every neighbouring pair (pairs_all.csv)")
    p.add_argument("--orbits", type=float, default=1e5, help="integration length in inner orbits")
    p.add_argument("-j", "--jobs", type=int, default=-1)
    p.add_argument("--quick", action="store_true", help="30 pairs, 1e3 orbits (test)")
    p.add_argument("--plot-only", action="store_true")
    a = p.parse_args()
    if a.all:
        run_all(1e3 if a.quick else a.orbits, a.jobs, a.quick)
        return
    from joblib import Parallel, delayed
    out = C.out_dir("stability_nbody")
    pairs = list(csv.DictReader(open(C.DATA / "pairs_for_nbody.csv")))
    if a.quick:
        pairs, a.orbits = pairs[::20], 1e3
    path = out / "pairs.csv"
    if not a.plot_only:
        t = C.Timer()
        res = Parallel(n_jobs=a.jobs, verbose=2)(delayed(integrate)(r, a.orbits) for r in pairs)
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["group", "scenario", "hill_ratio", "period_ratio", "e1", "e2", "mutual_inc_deg", "status",
                        "megno", "orbits_done"])
            from gaiabias.orbits import mutual_inclination
            for r, (st, mg, nd) in zip(pairs, res):
                mi = np.degrees(mutual_inclination(float(r["inc1"]), float(r["Omega1"]), float(r["inc2"]),
                                                   float(r["Omega2"])))
                w.writerow([r["group"], r["scenario"], r["hill_ratio"], r["period_ratio"], r["e1"], r["e2"],
                            f"{mi:.1f}", st, f"{mg:.3f}", f"{nd:.0f}"])
        print(f"integrated {len(pairs)} pairs in {t}")
    rows = list(csv.DictReader(open(path)))
    L = ["=== N-BODY CHECK OF THE HILL-STABILITY FILTER ===",
         f"{len(rows)} planet pairs from the step 3 population, integrated for {a.orbits:g} inner orbits (WHFast, MEGNO).",
         "unstable = close encounter (< 1 mutual Hill radius), ejection or MEGNO > 4; chaotic = MEGNO 2.5-4.", ""]
    for g in ("borderline", "near 2:1", "control"):
        sel = [r for r in rows if r["group"] == g]
        if not sel:
            continue
        cnt = {s: sum(r["status"] == s for r in sel) for s in ("stable", "chaotic", "unstable", "failed")}
        L.append(f"  {g:<11} n={len(sel):>4}: stable {cnt['stable']}, chaotic {cnt['chaotic']}, unstable "
                 f"{cnt['unstable']} ({100 * cnt['unstable'] / len(sel):.0f}%), failed {cnt['failed']}")
        mi = np.array([float(r["mutual_inc_deg"]) for r in sel])
        un = np.array([r["status"] == "unstable" for r in sel])
        for lo, hi in ((0, 30), (30, 60), (60, 180)):
            m = (mi >= lo) & (mi < hi)
            if m.sum() >= 10:
                L.append(f"      mutual inclination {lo}-{hi} deg: unstable {100 * un[m].mean():.0f}% of {m.sum()}")
    L += ["", "In the full simulation 43% of neighbouring giant pairs lie in the borderline range (1.0-1.6 times",
          "the required gap); the unstable share there sets how many simulated systems are unphysical."]
    C.write_summary(out / "summary.txt", L)
    plt, ps = C.setup_plots()
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6), constrained_layout=True)
    col = {"stable": ps.C_TARGET, "chaotic": ps.C_ACCEL, "unstable": ps.C_WRONG, "failed": ps.C_NONE}
    for s in col:
        sel = [r for r in rows if r["status"] == s]
        ax[0].scatter([float(r["hill_ratio"]) for r in sel], [max(float(r["e1"]), float(r["e2"])) for r in sel],
                      s=10, color=col[s], label=f"{s} ({len(sel)})", alpha=0.7)
        ax[1].scatter([float(r["period_ratio"]) for r in sel], [float(r["mutual_inc_deg"]) for r in sel], s=10,
                      color=col[s], alpha=0.7)
    ax[0].set_xscale("log")
    ax[0].set_xlabel("gap / gap required by the Hill filter")
    ax[0].set_ylabel("larger eccentricity of the pair")
    ax[0].legend(fontsize=8.5)
    ax[1].set_xscale("log")
    ax[1].set_xlabel("period ratio")
    ax[1].set_ylabel("mutual inclination (deg)")
    fig.suptitle(f"N-body check of simulated planet pairs ({a.orbits:g} inner orbits)", fontsize=12)
    fig.savefig(out / "stability_nbody.png", dpi=150, bbox_inches="tight")
    print(f"saved {out / 'stability_nbody.png'}")


if __name__ == "__main__":
    main()
