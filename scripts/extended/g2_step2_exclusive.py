"""Step 2 outcomes as mutually exclusive categories.

Reads results/step2_population.csv (made by scripts/step2.py run); no new fits, a minute:

    python scripts/extended/g2_step2_exclusive.py

The step 2 summary gives, per analysis, the share of systems in which the target is
correct, the companion is correct, or a wrong planet is reported; a two-orbit solution can
do several of these at once, so those shares can add to more than 100%. This prints the
same results split into categories that add up to exactly 100% of the systems:
both planets correct / target only / companion only (target missed) / a correct planet
plus a wrong one / wrong planet only / acceleration or jerk only / nothing reported,
for all systems and for each kind of companion, plus the false-alarm rate.
Output: results/extended/step2_exclusive/summary.txt.
"""
from __future__ import annotations

import numpy as np

import common as C


def main():
    from gaiabias.step2 import population as POP
    from gaiabias.step2.analyses import ANALYSES, LABEL
    from gaiabias.step2.outcomes import exclusive_categories, judge_all
    from gaiabias.step2.plots import regimes

    if not POP.OUT_CSV.exists():
        raise SystemExit(f"{POP.OUT_CSV} not found: this needs the machine where step 2 was run")
    d = POP.load()
    outs = judge_all(d)
    out = C.out_dir("step2_exclusive")
    n = len(d["seed"])
    cats = list(exclusive_categories(outs["one"]).keys())
    short = ["both ok", "target", "comp. only", "ok + wrong", "wrong only", "acc/jerk", "nothing"]
    L = ["=== STEP 2 OUTCOMES AS MUTUALLY EXCLUSIVE CATEGORIES (% of systems; rows add to 100) ===",
         f"{n} systems. Columns: " + "; ".join(f"{s} = {c}" for s, c in zip(short, cats)) + ".",
         "'any planet' = first five columns; [false alarm] = a non-target or second planet reported when the "
         "companion is removed.", ""]
    for name, mask in [("all systems", np.ones(n, bool))] + [(l.replace("\n", " "), m) for l, m in regimes(d)]:
        if not mask.any():
            continue
        L.append(f"--- {name} (n={int(mask.sum())}) ---")
        L.append(f"   {'analysis':<30}" + "".join(f"{s:>11}" for s in short) + f"{'any planet':>12}{'[false al.]':>12}")
        for a in ANALYSES:
            ex = exclusive_categories(outs[a])
            v = [100 * ex[c][mask].mean() for c in cats]
            L.append(f"   {LABEL[a]:<30}" + "".join(f"{x:>11.1f}" for x in v) + f"{sum(v[:5]):>12.1f}"
                     f"{100 * outs[a]['false_alarm'][mask].mean():>12.1f}")
        L.append("")
    C.write_summary(out / "summary.txt", L)


if __name__ == "__main__":
    main()
