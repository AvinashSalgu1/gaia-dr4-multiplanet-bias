"""Which choice is right, the BIC or the Bayesian evidence?

Reads results/step2_population.csv and the evidence runs of g2_evidence.py
(results/extended/evidence*/evidence.csv: the default priors and, if present, the runs with
--prior-scale 3 and --prior-scale 0.3). No new fits, seconds:

    python scripts/extended/g2_evidence_truth.py

For each system the planets reported by the model that the BIC adopts and by the model that the
evidence adopts are judged against the true system with the step 2 rules (target correct,
companion correct, wrong planet; the reference for the target is the BIC choice on the
target-alone data, as in step 2). The summary gives, per prior scale and per group, the share of
systems in each outcome for both choices, and the same for the systems where the two choices
differ. Output: results/extended/evidence_truth/summary.txt.
"""
from __future__ import annotations

import csv

import numpy as np

import common as C

MODELS = ["no planet", "acceleration only", "jerk only", "one orbit", "orbit + acceleration", "orbit + jerk",
          "two orbits"]
GROUPS = ["weak", "eccentric", "near 2:1", "random"]


def judge_with_choice(d, choice):
    """Step 2 outcome flags when the adopted model (index into MODELS) is `choice` on the data with companion."""
    from gaiabias.step2 import outcomes as O
    original = O.bic_choice

    def forced(dd, v):
        return np.asarray(choice) if v == "w" else original(dd, v)

    O.bic_choice = forced
    try:
        return O.judge(d, "best")
    finally:
        O.bic_choice = original


def outcome_shares(o, sel):
    t, c, w = o["target"][sel], o["companion"][sel], o["wrong"][sel]
    n = max(int(sel.sum()), 1)
    return dict(n=int(sel.sum()), target=100 * t.sum() / n, companion=100 * c.sum() / n,
                both=100 * (t & c & ~w).sum() / n, wrong=100 * w.sum() / n,
                none=100 * (~t & ~c & ~w).sum() / n)


def main():
    from gaiabias.step2 import population as POP
    if not POP.OUT_CSV.exists():
        raise SystemExit(f"{POP.OUT_CSV} not found: this needs the machine where step 2 was run")
    full = POP.load()
    index = {int(s): i for i, s in enumerate(full["seed"])}
    runs = sorted(p for p in C.OUT.glob("evidence*/evidence.csv") if "quick" not in p.parent.name)
    if not runs:
        raise SystemExit("no evidence runs found: run g2_evidence.py first")
    out = C.out_dir("evidence_truth")
    L = ["=== BIC CHOICE VS EVIDENCE CHOICE, JUDGED AGAINST THE TRUE SYSTEMS ===",
         "Percent of systems: target correct / companion correct / both planets correct / a wrong planet reported / "
         "nothing correct or wrong (no planet reported).",
         "Every system has both planets; 'target' and 'companion' as in step 2.", ""]
    for path in runs:
        with open(path) as f:
            ev = [r for r in csv.DictReader(f) if int(r["seed"]) in index]
        idx = np.array([index[int(r["seed"])] for r in ev])
        d = {k: (np.asarray(v)[idx] if hasattr(v, "__len__") and len(v) == len(full["seed"]) else v)
             for k, v in full.items()}
        ch_bic = np.array([MODELS.index(r["bic_choice"]) for r in ev])
        ch_ev = np.array([MODELS.index(r["evidence_choice"]) for r in ev])
        o_bic, o_ev = judge_with_choice(d, ch_bic), judge_with_choice(d, ch_ev)
        groups = np.array([r["group"] for r in ev])
        differ = ch_bic != ch_ev
        L.append(f"--- {path.parent.name} ({len(ev)} systems; choices differ in {int(differ.sum())}) ---")
        L.append(f"   {'subset':<34}{'choice':<10}{'target':>8}{'comp.':>8}{'both':>8}{'wrong':>8}{'nothing':>9}")
        subsets = [("all", np.ones(len(ev), bool)), ("choices differ", differ)] + \
                  [(g, groups == g) for g in GROUPS if (groups == g).any()] + \
                  [("BIC trend model, evidence two orbits", np.isin(ch_bic, [4, 5]) & (ch_ev == 6)),
                   ("BIC two orbits, evidence one orbit", (ch_bic == 6) & (ch_ev == 3))]
        for name, sel in subsets:
            if not sel.any():
                continue
            for lab, o in (("BIC", o_bic), ("evidence", o_ev)):
                s = outcome_shares(o, sel)
                L.append(f"   {(name + f' (n={s['n']})') if lab == 'BIC' else '':<34}{lab:<10}{s['target']:>8.1f}"
                         f"{s['companion']:>8.1f}{s['both']:>8.1f}{s['wrong']:>8.1f}{s['none']:>9.1f}")
        L.append("")
    C.write_summary(out / "summary.txt", L)


if __name__ == "__main__":
    main()
