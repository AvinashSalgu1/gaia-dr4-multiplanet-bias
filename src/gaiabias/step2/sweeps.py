"""Controlled experiments: fix one two-planet system, change ONE thing at a
time, and repeat each setting many times with fresh noise (everything else
identical), so each point is a probability with an error bar.

Reference systems (Sun-like star, G = 10, target 3 M_J at 1.5 yr, SNR1 = 4,
e = 0.2, tilted 60 deg; the same Gaia epochs for every run):
    outer : companion at 15 yr with 3x the target's wobble (about 1.9 M_J)
    inner : companion at 0.3 yr with half the target's wobble (about 4.3 M_J)

1D sweeps: target SNR1 (system moved closer / further), target period,
target eccentricity, wobble ratio (= companion mass), companion period (at
fixed wobble ratio), companion eccentricity, companion position in its orbit
at mid-mission, angle between the orbits, tilt of the whole system.
Settings that give a companion above 13 M_J or an unstable pair are skipped.

Maps: companion period x wobble ratio, in slices of a third parameter
(default companion eccentricity 0 / 0.3 / 0.6).

All four analyses are run on the data with the companion (the Gaia-style
pipeline in sweeps only); the one-planet fit is also run on the same data
without the companion, which is the reference for "found correctly".

--vary-angles: draw the orientations and orbital phases at random for every
repeat (except the swept one) instead of fixing them.
"""
from __future__ import annotations

import csv
import time
import zlib

import numpy as np

from . import analyses as A
from .sim import RESULTS

CSV_PATH = RESULTS / "step2_sweeps.csv"
T_MID_DAYS = 0.5 * (2456891.5 + 2458868.5) - 2457936.875   # DR4 mid-point vs reference epoch
YR = 365.25

COMMON = dict(snr1_t=4.0, P_t_yr=1.5, e_t=0.2, inc_t_deg=60.0, Omega_t_deg=40.0, omega_t_deg=60.0,
              M_t_deg=0.0, mutual_deg=30.0, azimuth_deg=90.0, omega_c_deg=120.0,
              G=10.0, m_star=1.0, m_t_ref=3.0, jitter=0.02)
BASELINES = {
    "outer": dict(COMMON, P_c_yr=15.0, rho=3.0, e_c=0.2, M_c_deg=90.0),
    "inner": dict(COMMON, P_c_yr=0.3, rho=0.5, e_c=0.1, M_c_deg=90.0),
}
SWEEPS = {
    "snr1_t": ("target's signal strength (SNR per transit)\n(same system moved closer / further)",
               np.geomspace(1.5, 10, 9), True),
    "P_t_yr": ("target's period (years)", np.geomspace(0.3, 3.8, 9), True),
    "e_t": ("target's eccentricity", np.linspace(0, 0.8, 9), False),
    "rho": ("companion's wobble ÷ target's wobble\n(= heavier companion)", np.geomspace(0.01, 100, 13), True),
    "P_c_yr": ("companion's period (years)\n(wobble ratio kept fixed)", np.geomspace(0.1, 54, 13), True),
    "e_c": ("companion's eccentricity", np.linspace(0, 0.9, 10), False),
    "M_c_deg": ("companion's place in its orbit at mid-mission\n(degrees; 0 = closest to the star)",
                np.arange(0, 360, 30.0), False),
    "mutual_deg": ("angle between the two orbits (degrees)", np.arange(0, 181, 20.0), False),
    "inc_t_deg": ("tilt of the whole system (0 = face-on,\n90 = edge-on), degrees", np.linspace(0, 90, 7), False),
}
MAP_X = np.geomspace(0.1, 54, 10)        # companion period [yr]
MAP_Y = np.geomspace(0.01, 100, 9)       # wobble ratio
ANGLE_KEYS = ["inc_t_deg", "Omega_t_deg", "omega_t_deg", "M_t_deg", "mutual_deg", "azimuth_deg",
              "omega_c_deg", "M_c_deg"]

FIT_COLS = (A.cols("one", "w") + A.cols("one", "a") + A.cols("accel", "w") + A.cols("jerk", "w")
            + A.cols("two", "w"))
COLUMNS = (["task", "study", "kind", "param", "i", "value", "param2", "j", "value2", "slice", "k",
            "slice_value", "rep", "valid", "n_transits", "m_t", "P_t", "m_c", "P_c", "rho"] + FIT_COLS)


# ----------------------------------------------------------------------------- systems
def reference_parallax(p):
    """Parallax giving the reference target (m_t_ref at P_t with SNR1 = 4)."""
    from ..orbits import reflex_a0_mas
    from ..simulate import sigma_al_per_transit
    a0_per_mas = reflex_a0_mas(p["m_star"], p["m_t_ref"], p["P_t_yr"] * YR, 1.0)
    return COMMON["snr1_t"] * float(sigma_al_per_transit(p["G"])) / a0_per_mas


def build(p, plx_ref):
    """(target, companion, parallax) for parameter set p, or None if not a stable planet pair."""
    from ..orbits import (Companion, Tp_for_mean_anomaly, companion_mass_from_a0, hill_stable,
                          orientation_relative_to)
    from ..simulate import sigma_al_per_transit

    plx = plx_ref * p["snr1_t"] / COMMON["snr1_t"]       # SNR sweep = move the system
    a0_t = p["snr1_t"] * float(sigma_al_per_transit(p["G"]))
    P_t, P_c = p["P_t_yr"] * YR, p["P_c_yr"] * YR
    m_t = companion_mass_from_a0(a0_t, P_t, p["m_star"], plx)
    m_c = companion_mass_from_a0(p["rho"] * a0_t, P_c, p["m_star"], plx)
    if not (0.1 < m_t < 13.0 and 0.03 < m_c < 13.0):
        return None
    d = np.radians
    inc_t, Om_t = d(p["inc_t_deg"]), d(p["Omega_t_deg"])
    inc_c, Om_c = orientation_relative_to(inc_t, Om_t, d(p["mutual_deg"]), d(p["azimuth_deg"]))
    t = Companion(m_t, P_t, p["e_t"], inc_t, Om_t, d(p["omega_t_deg"]),
                  Tp_for_mean_anomaly(d(p["M_t_deg"]), T_MID_DAYS, P_t))
    c = Companion(m_c, P_c, p["e_c"], inc_c, Om_c, d(p["omega_c_deg"]),
                  Tp_for_mean_anomaly(d(p["M_c_deg"]), T_MID_DAYS, P_c))
    if not hill_stable(p["m_star"], t, c):
        return None
    return t, c, plx


def randomise_angles(p, rng, keep):
    q = dict(p)
    for k in ANGLE_KEYS:
        if k == keep:
            continue
        q[k] = np.degrees(np.arccos(rng.uniform(-1, 1))) if k in ("inc_t_deg", "mutual_deg") \
            else rng.uniform(0, 360)
    return q


def reference_epochs():
    """Fixed DR4 epochs at a sky position with a typical number of transits (~78)."""
    from .. import epochs as E
    best = None
    for ra, dec in [(120.0, -30.0), (45.0, 20.0), (200.0, 10.0), (300.0, -50.0), (80.0, 60.0),
                    (160.0, -10.0), (250.0, 35.0), (10.0, -70.0)]:
        ep = E.from_gaiamock(ra, dec, "dr4", rng=np.random.default_rng(2024))
        if best is None or abs(len(ep) - 78) < abs(len(best[2]) - 78):
            best = (ra, dec, ep)
    return best


def simulate_pair(p, plx, t, c, ep, rng):
    """(obs alone, obs with companion, err) with one shared noise draw."""
    from ..simulate import Star, sigma_al_per_transit, simulate
    star = Star(parallax_mas=plx, m_star_msun=p["m_star"], G_mag=p["G"], pmra=20.0, pmdec=-10.0)
    single = simulate(ep, star, [t], add_noise=False)
    both = simulate(ep, star, [t, c], add_noise=False)
    noise = float(sigma_al_per_transit(p["G"])) * rng.standard_normal(len(ep))
    if p["G"] < 13:
        noise += p["jitter"] * rng.standard_normal(len(ep))
    return single.obs + noise, both.obs + noise, both.err


# ----------------------------------------------------------------------------- one run
def run_task(task, ep, plx_ref, seed0, with_gaia=False):
    rng = np.random.default_rng([seed0, zlib.crc32(task["task"].encode())])
    p = dict(task["params"])
    if task["vary"]:
        p = randomise_angles(p, rng, keep=task["param"])
    row = {k: v for k, v in task.items() if k in COLUMNS}
    row.update(valid=0, n_transits=len(ep), m_t=np.nan, P_t=p["P_t_yr"] * YR, m_c=np.nan, P_c=p["P_c_yr"] * YR, rho=p["rho"])
    sysm = build(p, plx_ref)
    if sysm is not None:
        t, c, plx = sysm
        row.update(valid=1, m_t=t.m_mjup, m_c=c.m_mjup)
        obs_a, obs_w, err = simulate_pair(p, plx, t, c, ep, rng)
        which = {"one", "accel", "jerk", "two"}
        row.update(A.run(ep, obs_w, err, p["m_star"], "w", which))
        row.update(A.run(ep, obs_a, err, p["m_star"], "a", {"one"}))
    for k in COLUMNS:
        row.setdefault(k, np.nan)
    return row


def _task(study, kind, params, vary, rep, param="", i=-1, value=np.nan, param2="", j=-1, value2=np.nan,
          slc="", k=-1, slice_value=np.nan):
    key = f"{study}|{kind}|{param}|{i}|{param2}|{j}|{slc}|{k}|{rep}"
    return dict(task=key, study=study, kind=kind, param=param, i=i, value=value, param2=param2, j=j,
                value2=value2, slice=slc, k=k, slice_value=slice_value, rep=rep, params=params, vary=vary)


def sweep_tasks(base_name, vary, reps, quick=False):
    base = BASELINES[base_name]
    study = f"{base_name}-{'varied' if vary else 'fixed'}"
    tasks = [_task(study, "baseline", base, vary, r) for r in range(reps * 2)]
    for key, (_, values, _) in SWEEPS.items():
        for i, v in enumerate(values[::3] if quick else values):
            for r in range(reps):
                tasks.append(_task(study, "sweep", dict(base, **{key: float(v)}), vary, r, key, i, float(v)))
    return tasks


def map_tasks(base_name, vary, reps, slc, slice_values, quick=False):
    base = BASELINES[base_name]
    study = f"{base_name}-{'varied' if vary else 'fixed'}"
    xs, ys = (MAP_X[::3], MAP_Y[::3]) if quick else (MAP_X, MAP_Y)
    tasks = []
    for k, sv in enumerate(slice_values):
        for i, x in enumerate(xs):
            for j, y in enumerate(ys):
                for r in range(reps):
                    p = dict(base, P_c_yr=float(x), rho=float(y), **{slc: float(sv)})
                    tasks.append(_task(study, "map", p, vary, r, "P_c_yr", i, float(x), "rho", j, float(y),
                                       slc, k, float(sv)))
    return tasks


def run(tasks, csv_path=CSV_PATH, n_jobs=-1, seed0=13, with_gaia=False, chunk=800):
    from joblib import Parallel, delayed

    RESULTS.mkdir(exist_ok=True)
    done = set()
    if csv_path.exists():
        with open(csv_path) as f:
            done = {r["task"] for r in csv.DictReader(f)}
    else:
        with open(csv_path, "w", newline="") as f:
            csv.DictWriter(f, fieldnames=COLUMNS).writeheader()
    todo = [t for t in tasks if t["task"] not in done]
    ra, dec, ep = reference_epochs()
    print(f"reference star at RA {ra}, Dec {dec}: {len(ep)} transits; "
          f"{len(tasks) - len(todo)} runs already done, {len(todo)} to go", flush=True)
    plx = {b: reference_parallax(BASELINES[b]) for b in BASELINES}
    t0 = time.perf_counter()
    for s in range(0, len(todo), chunk):
        part = todo[s:s + chunk]
        rows = Parallel(n_jobs=n_jobs)(
            delayed(run_task)(t, ep, plx[t["study"].split("-")[0]], seed0, with_gaia) for t in part)
        with open(csv_path, "a", newline="") as f:
            w = csv.DictWriter(f, fieldnames=COLUMNS)
            for r in rows:
                w.writerow({k: r.get(k, "") for k in COLUMNS})
        n = s + len(part)
        el = (time.perf_counter() - t0) / 60
        print(f"  {n}/{len(todo)}  ({el:.1f} min, ~{el / n * (len(todo) - n):.0f} min left)", flush=True)


# ----------------------------------------------------------------------------- analysis
def load(csv_path=CSV_PATH):
    with open(csv_path) as f:
        rows = list(csv.DictReader(f))
    txt = {"task", "study", "kind", "param", "param2", "slice"}
    return {k: np.array([r[k] for r in rows]) if k in txt else
            np.array([float(r[k]) if r[k] not in ("", "nan") else np.nan for r in rows]) for k in COLUMNS}


def target_correct(d):
    """Outcomes per analysis; reference = the one-planet fit on the target-alone data
    (the other analyses are not run on the target-alone data in the sweeps)."""
    from .outcomes import judge, reported
    alone = reported(d, "one", "a")
    out = {}
    for a in A.ANALYSES:
        out[a] = judge(d, a, alone=alone)
    return out


def _rate(mask, flag):
    n = mask.sum()
    if n == 0:
        return np.nan, np.nan
    p = flag[mask].mean()
    return p, np.sqrt(max(p * (1 - p), 1.0 / n) / n)


STYLE = {"one": dict(color="#9bbcf0", ls="-"), "accel": dict(color="#4f86dc", ls="-"),
         "jerk": dict(color="#7a5195", ls="-"),
         "two": dict(color="#1d3f8f", ls="-"),
         "best": dict(color="#2ca58d", ls="-")}


def summary(d, outs):
    valid = d["valid"] == 1
    lines = ["=== STEP 2 SWEEPS SUMMARY ===",
             "% of repeats in which the TARGET is found correctly: "
             + " / ".join(A.LABEL[a] for a in outs) + "   [one planet, companion removed]"]
    for study in sorted(set(d["study"])):
        S = d["study"] == study
        b = S & (d["kind"] == "baseline") & valid
        if b.any():
            lines.append(f"\n[{study}] reference system (target {np.nanmedian(d['m_t'][b]):.1f} M_J, companion "
                         f"{np.nanmedian(d['m_c'][b]):.1f} M_J, {b.sum()} repeats): "
                         + " / ".join(f"{100 * outs[a]['target'][b].mean():.0f}" for a in outs)
                         + f"   [{100 * outs['one']['alone_detected'][b].mean():.0f}]")
        for key, (label, _, _) in SWEEPS.items():
            m = S & (d["kind"] == "sweep") & (d["param"] == key)
            if not m.any():
                continue
            vals = np.unique(d["value"][m])
            rows = []
            for v in vals:
                mv = m & valid & (d["value"] == v)
                if mv.any():
                    rows.append((v, {a: outs[a]["target"][mv].mean() for a in outs}))
            if not rows:
                continue
            skipped = [f"{v:g}" for v in vals if not (m & valid & (d["value"] == v)).any()]
            txt = f"  {label.split(chr(10))[0][:46]:<46}"
            for a in outs:
                w = min(rows, key=lambda r: r[1][a])
                bst = max(rows, key=lambda r: r[1][a])
                txt += f" | {a}: {100 * w[1][a]:3.0f}% at {w[0]:.3g} .. {100 * bst[1][a]:3.0f}% at {bst[0]:.3g}"
            if skipped:
                txt += f"  (not a stable planet at: {', '.join(skipped)})"
            lines.append(txt)
        mp = S & (d["kind"] == "map") & valid
        for sv in np.unique(d["slice_value"][mp]):
            sl = mp & (d["slice_value"] == sv)
            for a in outs:
                cells = []
                for i, j in set(zip(d["i"][sl].astype(int), d["j"][sl].astype(int))):
                    c = sl & (d["i"] == i) & (d["j"] == j)
                    cells.append((outs[a]["target"][c].mean(), d["value"][c][0], d["value2"][c][0]))
                w, bb = min(cells), max(cells)
                lines.append(f"  map {d['slice'][sl][0]} = {sv:g}, {A.LABEL[a]}: worst {100 * w[0]:.0f}% "
                             f"(companion {w[1]:.2g} yr, wobble ratio {w[2]:.2g}), best {100 * bb[0]:.0f}% "
                             f"(companion {bb[1]:.2g} yr, ratio {bb[2]:.2g})")
    return lines


def plot_sweeps(d, outs, study, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
    from .. import plotstyle as ps

    S = (d["study"] == study) & (d["kind"] == "sweep")
    if not S.any():
        return
    ps.apply()
    base = BASELINES[study.split("-")[0]]
    keys = [k for k in SWEEPS if (S & (d["param"] == k)).any()]
    ncol = 3
    nrow = int(np.ceil(len(keys) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(15, 3.9 * nrow), constrained_layout=True, sharey=True,
                            squeeze=False)
    axs = axs.ravel()
    valid = d["valid"] == 1
    for a, key in zip(axs, keys):
        label, _, logx = SWEEPS[key]
        m = S & (d["param"] == key)
        vals = np.unique(d["value"][m])
        good = np.array([(m & valid & (d["value"] == v)).any() for v in vals])
        for an in outs:
            r = np.array([_rate(m & valid & (d["value"] == v), outs[an]["target"]) for v in vals])
            a.plot(vals[good], 100 * r[good, 0], marker="o", ms=3.5, lw=1.8, label=A.LABEL[an], **STYLE[an])
            if an == "two":
                a.fill_between(vals[good], 100 * (r[good, 0] - r[good, 1]), 100 * (r[good, 0] + r[good, 1]),
                               color=STYLE[an]["color"], alpha=0.12, lw=0)
        ra = np.array([_rate(m & valid & (d["value"] == v), outs["one"]["alone_detected"]) for v in vals])
        a.plot(vals[good], 100 * ra[good, 0], color=ps.C_NONE, lw=1.3, ls="--", label="one planet, companion removed")
        rw = np.array([_rate(m & valid & (d["value"] == v), outs["one"]["wrong"]) for v in vals])
        a.plot(vals[good], 100 * rw[good, 0], color=ps.C_WRONG, lw=1.5, marker="s", ms=3,
               label="WRONG planet reported (one-planet fit)")
        step = np.diff(vals).min() if len(vals) > 1 else 1.0
        for v, g in zip(vals, good):
            if not g:
                a.axvspan(*((v / 1.25, v * 1.25) if logx else (v - 0.5 * step, v + 0.5 * step)),
                          color=ps.GRID, lw=0, zorder=0)
        if logx:
            a.set_xscale("log")
            ticks = [t for t in (0.01, 0.03, 0.1, 0.3, 1, 2, 3, 5, 10, 30, 100) if vals.min() <= t <= vals.max()]
            a.xaxis.set_major_locator(FixedLocator(ticks))
            a.xaxis.set_major_formatter(FuncFormatter(lambda v, _: {0.01: "1/100", 0.1: "1/10"}.get(v, f"{v:g}")))
            a.xaxis.set_minor_formatter(NullFormatter())
        if key == "P_c_yr":
            ps.mark_dr4(a, text=False)
        a.axvline(base[key], color=ps.MUTED, lw=1, ls="-.", zorder=1)
        a.set_xlabel(label, fontsize=9.5)
        a.set_ylim(-2, 102)
        a.grid(axis="y", color=ps.GRID, lw=0.8)
    for a in axs[len(keys):]:
        a.axis("off")
    for a in axs[::ncol]:
        a.set_ylabel("% of repeats: target\nfound correctly")
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.07), fontsize=9.5)
    fig.suptitle(f"Change one thing at a time ({study} reference: companion at {base['P_c_yr']:g} yr with "
                 f"{base['rho']:g}× the target's wobble; dash-dot = reference value; grey band = not a stable "
                 f"planet)", fontsize=11.5)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")


def plot_maps(d, outs, study, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from .. import plotstyle as ps

    M = (d["study"] == study) & (d["kind"] == "map")
    if not M.any():
        return
    ps.apply()
    valid = d["valid"] == 1
    svals = np.unique(d["slice_value"][M])
    slc = d["slice"][M][0]
    xv, yv = np.unique(d["value"][M]), np.unique(d["value2"][M])

    def edges(v):
        if len(v) == 1:
            return np.array([v[0] / 1.5, v[0] * 1.5])
        mid = np.sqrt(v[1:] * v[:-1])
        return np.concatenate([[v[0] ** 2 / mid[0]], mid, [v[-1] ** 2 / mid[-1]]])

    xe, ye = edges(xv), edges(yv)
    pipes = [a for a in ("one", "accel", "jerk", "two", "best") if a in outs]
    fig, axs = plt.subplots(len(pipes), len(svals), figsize=(4.6 * len(svals) + 1, 4.0 * len(pipes)),
                            constrained_layout=True, squeeze=False)
    im = None
    for r, an in enumerate(pipes):
        for c, sv in enumerate(svals):
            a = axs[r, c]
            g = np.full((len(yv), len(xv)), np.nan)
            for i, x in enumerate(xv):
                for j, y in enumerate(yv):
                    m = M & valid & (d["slice_value"] == sv) & (d["value"] == x) & (d["value2"] == y)
                    if m.any():
                        g[j, i] = outs[an]["target"][m].mean()
            a.set_facecolor(ps.GRID)
            im = a.pcolormesh(xe, ye, 100 * g, cmap="Blues", vmin=0, vmax=100)
            for i, x in enumerate(xv):
                for j, y in enumerate(yv):
                    if np.isfinite(g[j, i]):
                        a.text(x, y, f"{100 * g[j, i]:.0f}", ha="center", va="center", fontsize=7,
                               color="white" if g[j, i] > 0.6 else ps.INK)
            ps.years_axis(a, "x", lo=xe[0], hi=xe[-1])
            ps.ratio_axis(a, "y", lo=ye[0], hi=ye[-1])
            ps.mark_dr4(a, text=False)
            a.set_title(f"{A.LABEL[an]}\n{slc.replace('_', ' ')} = {sv:g}", fontsize=10.5)
    fig.colorbar(im, ax=axs, label="% of repeats: target found correctly", shrink=0.8)
    fig.suptitle(f"Best and worst cases: companion period × wobble ratio ({study}); "
                 "grey = not a stable planet (> 13 M_J or unstable)", fontsize=11.5)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")


def analyse(csv_path=CSV_PATH, make_plots=True):
    d = load(csv_path)
    outs = target_correct(d)
    lines = summary(d, outs)
    print("\n".join(lines))
    stem = csv_path.stem
    (RESULTS / f"{stem}_summary.txt").write_text("\n".join(lines) + "\n")
    if make_plots:
        for study in sorted(set(d["study"])):
            plot_sweeps(d, outs, study, RESULTS / f"{stem}_{study}.png")
            plot_maps(d, outs, study, RESULTS / f"{stem}_maps_{study}.png")
