"""Figure: N-body stability of the simulated giant pairs and its effect on the bias.

Reads results/extended/stability_all/pairs.csv (g2_stability_nbody.py --all) and
results/extended/stable_scan_nbody/summary.txt (g1_stable_scan.py --mask nbody); seconds:

    python scripts/extended/fig_stability.py

Output: results/extended/fig_stability.png
"""
from __future__ import annotations

import csv
import re

import numpy as np

import common as C


def scan_ratios(path):
    """{(method, sample): [(f, A/B)]} from the stable-scan summary (binned columns)."""
    out, method = {}, None
    for line in path.read_text().splitlines():
        if line.startswith("--- "):
            method = "one orbit" if "one planet" in line else "best by BIC"
        m = re.match(r"\s+([\d.]+)\s+(all|kept)\s+([\d.]+) \[[\d.,]+\]\s+([\d.]+)", line)
        if m and method:
            out.setdefault((method, m.group(2)), []).append((float(m.group(1)), float(m.group(3)) / float(m.group(4))))
    return out


def main():
    rows = list(csv.DictReader(open(C.OUT / "stability_all" / "pairs.csv")))
    mi = np.array([float(r["mutual_inc_deg"]) for r in rows])
    hr = np.array([float(r["hill_ratio"]) for r in rows])
    un = np.array([r["status"] == "unstable" for r in rows])
    sc = np.array([r["scenario"] for r in rows])
    plt, ps = C.setup_plots()
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.8), constrained_layout=True)
    edges = np.arange(0, 181, 15)
    for s, col in (("independent", ps.C_TARGET), ("clustered", ps.C_WRONG)):
        sel = sc == s
        frac = [un[sel & (mi >= lo) & (mi < hi)].mean() if (sel & (mi >= lo) & (mi < hi)).sum() > 20 else np.nan
                for lo, hi in zip(edges[:-1], edges[1:])]
        ax[0].plot(0.5 * (edges[:-1] + edges[1:]), frac, "o-", color=col, label=f"{s} ({sel.sum():,} pairs)")
    ps.percent_axis(ax[0], "y", "share of pairs that go unstable")
    ax[0].set_ylim(0, 1)
    ax[0].set_xlabel("mutual inclination of the pair (deg)")
    ax[0].set_title("(a) tilted pairs are unstable", fontsize=10.5)
    ax[0].legend(fontsize=8.5)
    hedges = np.geomspace(1, 30, 12)
    for s, col in (("independent", ps.C_TARGET), ("clustered", ps.C_WRONG)):
        for lo_i, hi_i, ls, lab in ((0, 30, "-", "inclined < 30 deg"), (30, 181, "--", "inclined > 30 deg")):
            sel = (sc == s) & (mi >= lo_i) & (mi < hi_i)
            frac = [un[sel & (hr >= lo) & (hr < hi)].mean() if (sel & (hr >= lo) & (hr < hi)).sum() > 20 else np.nan
                    for lo, hi in zip(hedges[:-1], hedges[1:])]
            ax[1].plot(np.sqrt(hedges[:-1] * hedges[1:]), frac, ls, color=col, marker="o", ms=3,
                       label=f"{s}, {lab}")
    ax[1].set_xscale("log")
    from matplotlib.ticker import FixedLocator, NullLocator
    ax[1].xaxis.set_major_locator(FixedLocator([1, 1.5, 2, 3, 4, 6]))
    ax[1].xaxis.set_minor_locator(NullLocator())
    ax[1].set_xticklabels(["1", "1.5", "2", "3", "4", "6"])
    ps.percent_axis(ax[1], "y", "share of pairs that go unstable")
    ax[1].set_ylim(0, 1)
    ax[1].set_xlabel("gap / gap required by the Hill filter")
    ax[1].set_title("(b) passing the Hill filter is not enough", fontsize=10.5)
    ax[1].legend(fontsize=7.5)
    r = scan_ratios(C.OUT / "stable_scan_nbody" / "summary.txt")
    for (method, sample), col, ls in ((("one orbit", "all"), ps.C_TARGET, "-"), (("one orbit", "kept"), ps.C_TARGET, "--"),
                                      (("best by BIC", "all"), ps.C_WRONG, "-"), (("best by BIC", "kept"), ps.C_WRONG, "--")):
        v = r.get((method, sample), [])
        if v:
            ax[2].plot([f for f, _ in v], [x for _, x in v], ls, marker="o", color=col,
                       label=f"{method}, {'all systems' if sample == 'all' else 'stable systems only'}")
    ax[2].axhline(1, color=ps.MUTED, lw=0.8, ls=":")
    ax[2].set_ylim(0.85, 1.05)
    ax[2].set_xlabel("fraction of giant hosts with a second giant")
    ax[2].set_ylabel("rate with companions / rate of same planets alone")
    ax[2].set_title("(c) the bias in the well-measured region is unchanged", fontsize=10.5)
    ax[2].legend(fontsize=8)
    out = C.OUT / "fig_stability.png"
    fig.savefig(out, dpi=170, bbox_inches="tight")
    print(f"saved {out}")


if __name__ == "__main__":
    main()
