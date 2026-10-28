"""Figure: how stellar binaries can look like planets in Gaia astrometry.

Rebuilds the stellar-binary systems of g2_mixed_flags.py (same seed, same step 3 stars; no fits)
and shows why some of them have planet-like photocentre orbits:

    python scripts/extended/fig_binaries.py [--root results/step3_gaia]

The fitted outcome (how many were reported as planets) is read from
results/extended/mixed_flags/summary.txt. Output: results/extended/fig_binaries.png
"""
from __future__ import annotations

import argparse
import re

import numpy as np

import common as C


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--root", type=C.Path, default=C.STEP3)
    a = p.parse_args()
    import g2_mixed_flags as M
    from gedi.orbits import reflex_a0_mas
    from gedi.simulate import snr1
    with np.load(a.root / "population.npz") as f:
        z = {k: f[k] for k in f.files if k.startswith("star_")}
    rows = [r for r in M.draw_systems(z, 3000, 6000) if r["kind"] == "stellar"]
    m_true = np.array([r["m_true"] for r in rows])
    m_eff = np.array([r["m_eff"] for r in rows])
    m1 = np.array([z["star_m_star"][r["star"]] for r in rows])
    q = m_true * 9.547919e-4 / m1
    s = np.array([snr1(reflex_a0_mas(z["star_m_star"][r["star"]], r["m_eff"], r["P"], z["star_parallax"][r["star"]])
                       if r["m_eff"] > 0 else 0.0, z["star_G"][r["star"]]) for r in rows])
    txt = (C.OUT / "mixed_flags" / "summary.txt").read_text()
    mm = re.search(r"passing the DR3 RUWE cut: \{[^}]*'stellar': (\d+)\}; reported as planets \(m < 13 M_J\): \{[^}]*'stellar': (\d+)\}", txt)
    n_ruwe, n_rep = (int(mm.group(1)), int(mm.group(2))) if mm else (None, None)

    plt, ps = C.setup_plots()
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.8), gridspec_kw=dict(width_ratios=[1, 1.15, 0.9]),
                           constrained_layout=True)
    for m, ls in ((1.0, "-"), (0.5, "--")):
        qq = np.linspace(max(0.1, 0.081 / m), 1, 200)
        ax[0].plot(qq, M.photocentre_factor(qq, m), ls, color=ps.INK, label=f"primary {m:g} M_sun")
    ax[0].axhline(0, color=ps.MUTED, lw=0.8)
    ax[0].set_xlabel("mass ratio q = M2 / M1")
    ax[0].set_ylabel("photocentre wobble / wobble of the primary")
    ax[0].set_title("(a) the companion's light cancels the wobble", fontsize=10.5)
    ax[0].legend(fontsize=8.5)
    planet = m_eff < 13
    ax[1].scatter(q[~planet], m_eff[~planet], s=4, alpha=0.3, color=ps.C_NONE, label="looks like a star or brown dwarf")
    ax[1].scatter(q[planet], m_eff[planet], s=6, alpha=0.6, color=ps.C_WRONG,
                  label=f"looks like a planet (< 13 M_J): {planet.sum()} of {len(rows)}")
    ax[1].axhline(13, color=ps.C_WRONG, lw=1, ls=":")
    ax[1].set_yscale("log")
    ax[1].set_ylim(0.05, 2000)
    ax[1].set_xlabel("mass ratio q = M2 / M1")
    ax[1].set_ylabel("mass Gaia would infer (M_J, dark companion)")
    ax[1].set_title("(b) apparent mass of 6,000 simulated stellar companions", fontsize=10.5)
    ax[1].legend(fontsize=8.5, loc="lower left")
    det = planet & (s >= 0.8)
    stages = [("simulated\nbinaries", len(rows)), ("planet-like\nwobble", int(planet.sum())),
              ("planet-like and\nlarge enough\nto detect", int(det.sum()))]
    if n_rep is not None:
        stages.append(("reported as a\nplanet after DR4\nfits and RUWE cut", n_rep))
    x = np.arange(len(stages))
    ax[2].bar(x, [v for _, v in stages], color=[ps.C_NONE, ps.C_WRONG, ps.C_WRONG, ps.C_TARGET][:len(stages)])
    for i, (_, v) in enumerate(stages):
        ax[2].text(i, v * 1.15, f"{v:,}", ha="center", fontsize=9)
    ax[2].set_yscale("log")
    ax[2].set_ylim(5, 20000)
    ax[2].set_xticks(x, [l for l, _ in stages], fontsize=8.5)
    ax[2].set_ylabel("number of systems")
    ax[2].set_title("(c) from binaries to false planets", fontsize=10.5)
    out = C.OUT / "fig_binaries.png"
    fig.savefig(out, dpi=170, bbox_inches="tight")
    print(f"saved {out}; planet-like {planet.sum()}, detectable {det.sum()}, reported {n_rep}")


if __name__ == "__main__":
    main()
