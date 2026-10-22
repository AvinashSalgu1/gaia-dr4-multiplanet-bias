"""Figure: robustness of the step 2 results to the noise model, and BIC vs Bayesian evidence.

Reads results/extended/noise_robustness/summary.txt and results/extended/evidence/evidence.csv
(made by g2_noise_robustness.py and g2_evidence.py); no fits, seconds:

    python scripts/extended/fig_noise_evidence.py

Output: results/extended/fig_noise_evidence.png
"""
from __future__ import annotations

import csv
import re

import numpy as np

import common as C

LABELS = {"baseline": "reference", "ccd_corr": "correlated\nCCDs", "transit_corr": "correlated\ntransits",
          "excess": "excess\nnoise", "underest": "errors 30%\nunderestimated"}
MODELS = ["no planet", "acceleration only", "jerk only", "one orbit", "orbit + acceleration", "orbit + jerk",
          "two orbits"]
SHORT = ["none", "acc.", "jerk", "orbit", "orbit\n+ acc.", "orbit\n+ jerk", "two\norbits"]


def read_noise(path):
    rows, cur = {}, None
    for line in path.read_text().splitlines():
        m = re.match(r"\s+(\w+)\s+\(n=\d+\): detected ([\d.]+)% -> ([\d.]+)%", line)
        if m:
            cur = m.group(1)
            rows[cur] = dict(alone=float(m.group(2)), with_c=float(m.group(3)))
            continue
        m = re.match(r"\s+(one orbit|best by BIC)\s+([\d.]+) /\s+([\d.]+) /\s+([\d.]+) /\s+([\d.]+)\s+\[([\d.]+)\]", line)
        if m and cur:
            key = "one" if m.group(1) == "one orbit" else "bic"
            rows[cur][key] = dict(target=float(m.group(2)), comp=float(m.group(3)), wrong=float(m.group(4)),
                                  nothing=float(m.group(5)), fa=float(m.group(6)))
    return rows


def main():
    noise = read_noise(C.OUT / "noise_robustness" / "summary.txt")
    with open(C.OUT / "evidence" / "evidence.csv") as f:
        ev = list(csv.DictReader(f))
    conf = np.zeros((7, 7), int)
    for r in ev:
        conf[MODELS.index(r["bic_choice"]), MODELS.index(r["evidence_choice"])] += 1

    plt, ps = C.setup_plots()
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.9), gridspec_kw=dict(width_ratios=[1.35, 1.0, 1.05]),
                           constrained_layout=True)
    keys = [k for k in LABELS if k in noise]
    x = np.arange(len(keys))
    w = 0.2
    for j, (an, met, col, lab) in enumerate([("one", "target", "#8fb8f0", "target correct, one orbit"),
                                             ("bic", "target", ps.C_TARGET, "target correct, best by BIC"),
                                             ("one", "wrong", "#f3b98c", "wrong planet, one orbit"),
                                             ("bic", "wrong", ps.C_WRONG, "wrong planet, best by BIC")]):
        v = [noise[k][an][met] / 100 for k in keys]
        ax[0].bar(x + (j - 1.5) * w, v, w * 0.95, color=col, label=lab)
    ax[0].set_xticks(x, [LABELS[k] for k in keys], fontsize=9)
    ps.percent_axis(ax[0], "y", "share of all systems")
    ax[0].set_ylim(0, 0.8)
    ax[0].legend(fontsize=8, ncol=2, loc="upper right")
    ax[0].set_title("(a) Step 2 outcomes under five noise models", fontsize=10.5)

    cost = [noise[k]["alone"] - noise[k]["with_c"] for k in keys]
    gain = [noise[k]["bic"]["target"] - noise[k]["one"]["target"] for k in keys]
    ax[1].plot(x, cost, "o-", color=ps.INK, label="detections lost to the companion\n(one orbit, alone - with companion)")
    ax[1].plot(x, gain, "s-", color=ps.C_TARGET, label="targets gained by the BIC choice\n(BIC - one orbit)")
    ax[1].set_xticks(x, [LABELS[k] for k in keys], fontsize=9)
    ax[1].set_ylabel("percentage points of all systems")
    ax[1].set_ylim(0, 20)
    ax[1].legend(fontsize=8, loc="lower left")
    ax[1].set_title("(b) Cost of the companion and gain from the BIC choice", fontsize=10.5)

    frac = conf / np.maximum(conf.sum(axis=1, keepdims=True), 1)
    ax[2].imshow(frac, cmap="Blues", vmin=0, vmax=1)
    for i in range(7):
        for j in range(7):
            if conf[i, j]:
                ax[2].text(j, i, str(conf[i, j]), ha="center", va="center", fontsize=8.5,
                           color="white" if frac[i, j] > 0.55 else ps.INK)
    ax[2].set_xticks(range(7), SHORT, fontsize=8)
    ax[2].set_yticks(range(7), [s.replace("\n", " ") for s in SHORT], fontsize=8)
    ax[2].set_xlabel("model with the highest Bayesian evidence")
    ax[2].set_ylabel("model with the lowest BIC")
    ax[2].set_title(f"(c) BIC vs evidence, {len(ev)} systems", fontsize=10.5)
    out = C.OUT / "fig_noise_evidence.png"
    fig.savefig(out, dpi=170, bbox_inches="tight")
    print(f"saved {out}")


if __name__ == "__main__":
    main()
