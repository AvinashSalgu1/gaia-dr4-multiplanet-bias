"""What else decides the outcome, beyond the companion's period and wobble ratio?

Companion period P_c and wobble ratio rho dominate. Each other property is
judged AFTER removing their effect: every system's outcome minus the average
outcome of the OTHER systems in its (P_c, rho) cell (12 x 10 log cells as in
the step-2 map). Binning that residual against a property shows whether the
property matters (flat = no). The spread expected from noise alone is measured
by shuffling the property (permutation test).
"""
from __future__ import annotations

import numpy as np

from .sim import T_DR4_DAYS as T

REGIMES = [("companion < 3.8 yr", lambda xc: xc < 0.7, "#3b7d6e"),
           ("companion 3.8-11 yr", lambda xc: (xc >= 0.7) & (xc < 2), "#7a5195"),
           ("companion > 11 yr", lambda xc: xc >= 2, "#d1495b")]
LOG_KEYS = {"snr1_t", "ratio", "seen_ratio", "P_t"}


def cell_residual(y, P_c, rho):
    xe = np.geomspace(30, 10 * T, 13)
    ye = np.geomspace(0.01, 100, 11)
    cell = np.clip(np.digitize(P_c, xe) - 1, 0, 11) * 10 + np.clip(np.digitize(rho, ye) - 1, 0, 9)
    r = np.full(len(y), np.nan)
    for k in np.unique(cell):
        m = cell == k
        n = m.sum()
        if n >= 10:
            r[m] = y[m] - (y[m].sum() - y[m]) / (n - 1)
    return r


def factor_table(d):
    """(key, plain-language label, values, bin edges or None for quantile bins)."""
    deg = 180 / np.pi
    left = d["dchi2_left_c_accel"] / np.maximum(d["dchi2_seen_c"], 1e-12)
    return [
        ("e_c", "companion's eccentricity", d["e_c"], None),
        ("M_mid_c", "companion's distance from periastron\nat mid-mission (degrees of orbit)",
         d["M_mid_c"] * deg, np.linspace(0, 180, 7)),
        ("cosi_c", "companion's orbit: edge-on (0)\nto face-on (1)", np.abs(np.cos(d["inc_c"])),
         np.linspace(0, 1, 6)),
        ("mutual", "angle between the two orbits (degrees)", d["mutual_inc"] * deg, np.linspace(0, 180, 7)),
        ("e_t", "target's eccentricity", d["e_t"], None),
        ("cosi_t", "target's orbit: edge-on (0)\nto face-on (1)", np.abs(np.cos(d["inc_t"])),
         np.linspace(0, 1, 6)),
        ("snr1_t", "target's signal strength (SNR per transit)", d["snr1_t"], None),
        ("P_t", "target's period (years)", d["P_t"] / 365.25, None),
        ("ratio", "companion period ÷ target period", d["P_c"] / d["P_t"], None),
        ("n_tr", "number of Gaia transits", d["n_transits"], None),
        ("G", "star brightness (G magnitude)", d["G"], np.array([8, 10, 12, 13, 14, 15.0])),
        ("seen_ratio", "companion's VISIBLE wobble ÷ target's\n(after removing straight-line motion)",
         np.sqrt(d["dchi2_seen_c"] / np.maximum(d["dchi2_seen_t"], 1e-12)), None),
        ("left_accel", "share of companion's visible wobble an\nacceleration term can NOT absorb",
         np.clip(left, 0, 1), None),
    ]


def binned(x, r, edges, valid, min_n=30):
    if not valid.any():
        return np.array([]), np.array([]), np.array([]), edges
    if edges is None:
        edges = np.unique(np.nanquantile(x[valid], np.linspace(0, 1, 7)))
    if len(edges) < 2:
        return np.array([]), np.array([]), np.array([]), edges
    ix = np.clip(np.digitize(x, edges) - 1, 0, len(edges) - 2)
    xs, ms, ses = [], [], []
    for k in range(len(edges) - 1):
        m = valid & (ix == k)
        if m.sum() < min_n:
            continue
        xs.append(np.median(x[m]))
        ms.append(r[m].mean())
        ses.append(r[m].std() / np.sqrt(m.sum()))
    return np.array(xs), np.array(ms), np.array(ses), edges


def spread(x, r, edges, valid):
    _, ms, _, _ = binned(x, r, edges, valid)
    return (ms.max() - ms.min()) if len(ms) > 1 else np.nan


def importance(d, outcome, n_perm=20, seed=0):
    """{(factor, regime): (spread, chance spread)} for one 0/1 outcome array."""
    rng = np.random.default_rng(seed)
    xc = d["P_c"] / T
    r = cell_residual(outcome.astype(float), d["P_c"], d["rho"])
    table = {}
    for key, _, x, edges in factor_table(d):
        for name, sel, _ in REGIMES:
            v = sel(xc) & np.isfinite(r) & np.isfinite(x)
            st = spread(x, r, edges, v)
            perm = []
            for _ in range(n_perm):
                xp = x.copy()
                xp[v] = rng.permutation(x[v])
                perm.append(spread(xp, r, edges, v))
            table[(key, name)] = (st, np.nanmean(perm) if np.isfinite(perm).any() else np.nan)
    return r, table


def ranking(table, keys):
    score = {}
    for k in keys:
        vals = [table[(k, n)][0] - table[(k, n)][1] for n, _, _ in REGIMES]
        vals = [v for v in vals if np.isfinite(v)]
        if vals:
            score[k] = max(vals)
    return sorted(score, key=lambda k: -score[k]), score


def plot_grid(d, r, table, order, path, title_outcome="the TARGET planet is found correctly"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
    from .. import plotstyle as ps

    if not order:
        print("  (too few systems for the 'what else matters' figure)")
        return
    ps.apply()
    xc = d["P_c"] / T
    byk = {f[0]: f for f in factor_table(d)}
    n, ncol = len(order), 4
    nrow = int(np.ceil(n / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(4.1 * ncol, 3.3 * nrow), constrained_layout=True,
                            sharey=True, squeeze=False)
    axs = axs.ravel()
    for a, key in zip(axs, order):
        _, label, x, edges = byk[key]
        for name, sel, col in REGIMES:
            v = sel(xc) & np.isfinite(r) & np.isfinite(x)
            xs, ms, ses, _ = binned(x, r, edges, v)
            if len(xs) < 2:
                continue
            a.fill_between(xs, 100 * (ms - ses), 100 * (ms + ses), color=col, alpha=0.18, lw=0)
            a.plot(xs, 100 * ms, "o-", color=col, ms=3.5, lw=1.6, label=name)
        a.axhline(0, color=ps.MUTED, lw=1)
        if key in LOG_KEYS:
            a.set_xscale("log")
            lo, hi = np.nanmin(x[x > 0]), np.nanmax(x)
            ticks = [t for t in (0.01, 0.03, 0.1, 0.3, 1, 2, 3, 5, 10, 30, 100) if lo <= t <= hi]
            a.xaxis.set_major_locator(FixedLocator(ticks))
            a.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
            a.xaxis.set_minor_formatter(NullFormatter())
        a.set_xlabel(label, fontsize=9.5)
        best = np.nanmax([table[(key, nm)][0] for nm, _, _ in REGIMES])
        chance = np.nanmean([table[(key, nm)][1] for nm, _, _ in REGIMES])
        a.set_title(f"spread up to {100 * best:.0f} points (chance alone: ~{100 * chance:.0f})",
                    fontsize=10, fontweight="normal", color=ps.INK)
        a.grid(axis="y", color=ps.GRID, lw=0.8)
    for a in axs[n:]:
        a.axis("off")
    for a in axs[::ncol]:
        a.set_ylabel("better (+) or worse (−) than\ntypical, percentage points", fontsize=9.5)
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.04), fontsize=10)
    fig.suptitle(f"What else decides whether {title_outcome}? (one-planet fit)\n"
                 "Each line compares systems with the same companion period and wobble ratio; "
                 "flat = does not matter.  Most important first.", fontsize=12)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {path}")
