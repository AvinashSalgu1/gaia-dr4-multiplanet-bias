"""Planet pairs for the N-body stability check (g2_stability_nbody.py).

Reads the step 3 population (results/step3_gaia/population.npz) and writes a small table
data/extended/pairs_for_nbody.csv, so that the N-body check can run on any machine without
the step 3 results:

    python scripts/extended/g1_export_pairs.py [--root results/step3_gaia]

Pairs of neighbouring planets (sorted by period) in multi-planet systems, in three groups:
* borderline: spacing just above the Hill criterion used in the simulation
  (gap / required gap between 1.0 and 1.6);
* near 2:1: period ratio 1.9-2.15 (where adding Keplerian orbits is least accurate);
* control: well separated (gap / required gap above 3).

With --all, every neighbouring pair of both scenarios is written to data/extended/pairs_all.csv
instead, for the full N-body run (g2_stability_nbody.py --all), which marks each simulated
multi-giant system as stable or not.
"""
from __future__ import annotations

import argparse
import csv

import numpy as np

import common as C

N_PER_GROUP = {"borderline": 300, "near 2:1": 200, "control": 100}
FIELDS = ("m", "P", "e", "inc", "Omega", "omega", "Tp")


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--root", type=C.Path, default=C.STEP3)
    p.add_argument("--all", action="store_true",
                   help="write every neighbouring pair to data/extended/pairs_all.csv (for the full N-body run)")
    a = p.parse_args()
    from gedi.constants import MJUP_MSUN
    from gedi.orbits import relative_semimajor_axis_au

    with np.load(a.root / "population.npz") as f:
        z = {k: f[k] for k in f.files}
    rng = np.random.default_rng(8)
    rows = []
    for scen in ("independent", "clustered"):
        st = z[f"{scen}_star"]
        order = np.lexsort((z[f"{scen}_P"], st))
        st_o = st[order]
        same = st_o[1:] == st_o[:-1]
        for k in np.nonzero(same)[0]:
            i, j = order[k], order[k + 1]
            ms = float(z["star_m_star"][st[i]])
            m1, m2 = z[f"{scen}_m"][i] * MJUP_MSUN, z[f"{scen}_m"][j] * MJUP_MSUN
            a1 = relative_semimajor_axis_au(z[f"{scen}_P"][i], ms + m1)
            a2 = relative_semimajor_axis_au(z[f"{scen}_P"][j], ms + m1 + m2)
            r_h = ((m1 + m2) / (3 * ms)) ** (1 / 3) * 0.5 * (a1 + a2)
            gap = a2 * (1 - z[f"{scen}_e"][j]) - a1 * (1 + z[f"{scen}_e"][i])
            ratio = gap / (2 * np.sqrt(3) * r_h)
            pr = z[f"{scen}_P"][j] / z[f"{scen}_P"][i]
            row = dict(scenario=scen, star=int(st[i]), m_star=ms, hill_ratio=float(ratio), period_ratio=float(pr))
            for n, idx in (("1", i), ("2", j)):
                for fld in FIELDS:
                    row[f"{fld}{n}"] = float(z[f"{scen}_{fld}"][idx])
            rows.append(row)
    if a.all:
        C.DATA.mkdir(parents=True, exist_ok=True)
        path = C.DATA / "pairs_all.csv"
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["group"] + list(rows[0]))
            w.writeheader()
            for r in rows:
                w.writerow({"group": "all", **{k: (f"{v:.8g}" if isinstance(v, float) else v) for k, v in r.items()}})
        print(f"saved {path} ({len(rows)} neighbouring pairs)")
        return
    hr = np.array([r["hill_ratio"] for r in rows])
    pr = np.array([r["period_ratio"] for r in rows])
    groups = {"borderline": np.nonzero((hr >= 1) & (hr < 1.6))[0],
              "near 2:1": np.nonzero((pr >= 1.9) & (pr <= 2.15))[0],
              "control": np.nonzero(hr > 3)[0]}
    out_rows = []
    for g, idx in groups.items():
        pick = rng.choice(idx, min(N_PER_GROUP[g], len(idx)), replace=False)
        for k in pick:
            out_rows.append(dict(group=g, **rows[k]))
        print(f"{g}: {len(idx)} pairs available, {len(pick)} exported")
    C.DATA.mkdir(parents=True, exist_ok=True)
    path = C.DATA / "pairs_for_nbody.csv"
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0]))
        w.writeheader()
        for r in out_rows:
            w.writerow({k: (f"{v:.8g}" if isinstance(v, float) else v) for k, v in r.items()})
    print(f"saved {path} ({len(out_rows)} pairs; {len(rows)} neighbouring pairs in total)")


if __name__ == "__main__":
    main()
