"""Figure: DR4 sensitivity of known inner-planet hosts and forecast outer-giant detections.

Reads results/extended/census/hosts.csv (made by g2_census.py); seconds:

    python scripts/extended/fig_census.py

Output: results/extended/fig_census.png
"""
from __future__ import annotations

import csv

import numpy as np

import common as C

CLASSES = ["hot Jupiter", "warm Jupiter", "small planet"]


def main():
    from scipy.stats import poisson
    with open(C.OUT / "census" / "hosts.csv") as f:
        rows = [r for r in csv.DictReader(f) if r["simulate"] in ("1", "True") and r["mu"] not in ("", "nan")]
    plt, ps = C.setup_plots()
    col = {"hot Jupiter": ps.C_WRONG, "warm Jupiter": "#7a5195", "small planet": ps.C_TARGET}
    fig, ax = plt.subplots(1, 4, figsize=(18, 4.4), gridspec_kw=dict(width_ratios=[1, 1, 1, 1.1]),
                           constrained_layout=True)
    rng = np.random.default_rng(1)
    for k, (key, lab) in enumerate([("p3_3MJ", "3 M_J"), ("p3_1MJ", "1 M_J"), ("p3_03MJ", "0.3 M_J")]):
        a = ax[k]
        for c in ["small planet", "warm Jupiter", "hot Jupiter"]:
            sel = [r for r in rows if r["cls"] == c]
            d = np.array([float(r["dist"]) for r in sel])
            p = np.array([float(r[key]) for r in sel]) + rng.uniform(-0.015, 0.015, len(sel))
            a.scatter(d, p, s=9, alpha=0.55, color=col[c], label=f"{c} ({len(sel)})" if k == 0 else None)
        a.set_xscale("log")
        a.axhline(0.5, color=ps.MUTED, lw=0.8, ls=":")
        a.set_xlabel("distance (pc)")
        ps.percent_axis(a, "y", "DR4 detection probability" if k == 0 else None)
        a.set_ylim(-0.05, 1.05)
        a.set_title(f"({'abc'[k]}) outer companion of {lab}, P = 3 yr", fontsize=10.5)
    ax[0].legend(fontsize=8, loc="lower left", markerscale=1.5)
    a = ax[3]
    for i, c in enumerate(CLASSES):
        mu = float(np.sum([float(r["mu"]) for r in rows if r["cls"] == c]))
        n = sum(r["cls"] == c for r in rows)
        for j, (scale, mk) in enumerate([(1, "o"), (3, "s")]):
            m = scale * mu
            lo68, hi68 = poisson.ppf([0.16, 0.84], m)
            lo95, hi95 = poisson.ppf([0.025, 0.975], m)
            x = i + (j - 0.5) * 0.28
            a.plot([x, x], [max(lo95, 0.5), hi95], color=col[c], lw=1.2)
            a.plot([x, x], [max(lo68, 0.5), hi68], color=col[c], lw=6, alpha=0.6, solid_capstyle="butt")
            a.plot(x, m, mk, color=col[c], ms=7, label=("field rate" if scale == 1 else "3 x field rate") if i == 0 else None)
    a.set_yscale("log")
    a.set_ylim(0.5, 200)
    a.set_xticks(range(len(CLASSES)), [f"{c}\n({sum(r['cls'] == c for r in rows)} hosts)" for c in CLASSES])
    a.set_ylabel("expected outer giants detected by DR4")
    a.legend(fontsize=8, loc="upper left")
    a.set_title("(d) forecast detections (68% thick, 95% thin)", fontsize=10.5)
    out = C.OUT / "fig_census.png"
    fig.savefig(out, dpi=170, bbox_inches="tight")
    print(f"saved {out}")


if __name__ == "__main__":
    main()
