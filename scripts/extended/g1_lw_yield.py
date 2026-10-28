"""Why our DR4 yield differs from Lammers & Winn (2026): the same planets counted their way.

Uses only the step 3 population (no fits; seconds):

    python scripts/extended/g1_lw_yield.py [--root results/step3_allsky]

L&W count a planet as a DR4 discovery when SNR_1 > 1.5 and P < 4 yr (their semi-analytic
yield), with sigma_fov = 54 micro-arcsec x max[1, 10^(0.2 (G - 14))]. Step 3 instead counts
detections from full orbit fits (all periods up to 8 yr, imprecise and spurious planets
included) with the gaiamock noise table. This script applies, to the real stars and simulated
planets of step 3 in each multiplicity scenario, (1) L&W's counting rule with the gaiamock
noise, and (2) L&W's counting rule with L&W's noise, and gives (2) by distance. All numbers are
for the stars within 250 pc with RUWE <= 1.4, scaled to the whole sky (sky fraction from the run).
Output: results/extended/lw_yield/summary.txt.
"""
from __future__ import annotations

import argparse

import numpy as np

import common as C


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--root", type=C.Path, default=C.STEP3)
    a = p.parse_args()
    import step3 as S3
    from gaiabias.orbits import reflex_a0_mas
    from gaiabias.simulate import sigma_al_per_transit

    z, meta, _, _ = C.load_step3(a.root)
    scale = 1.0 / S3.sky_fraction(meta)
    G, plx, ms = z["star_G"], z["star_parallax"], z["star_m_star"]
    out = C.out_dir("lw_yield")
    L = ["=== DR4 YIELD: STEP 3 PLANETS COUNTED WITH THE LAMMERS & WINN RULE ===",
         f"step 3 results: {a.root} ({meta['n_stars']:,} stars, sky fraction {1 / scale:.3f}); whole-sky numbers.",
         "Rule: SNR_1 > 1.5 and P < 4 yr. Full-fit yields are in the step 3 summary.", ""]
    edges = [0, 50, 100, 150, 200, 250]
    for scen in S3.SCENARIOS:
        s, m, P = z[f"{scen}_star"], z[f"{scen}_m"], z[f"{scen}_P"]
        a0 = reflex_a0_mas(ms[s], m, P, plx[s])
        snr_gm = a0 / sigma_al_per_transit(G[s])
        snr_lw = a0 / (0.054 * np.maximum(1.0, 10 ** (0.2 * (G[s] - 14.0))))
        short = P / 365.25 < 4.0
        n_gm = scale * np.sum((snr_gm > 1.5) & short)
        sel = (snr_lw > 1.5) & short
        n_lw = scale * np.sum(sel)
        d = 1000.0 / plx[s][sel]
        by_d = [scale * np.sum((d >= lo) & (d < hi)) for lo, hi in zip(edges[:-1], edges[1:])]
        med = float(np.median(1000 * sigma_al_per_transit(G[(G < 14)])))
        L += [f"[{scen}]",
              f"  L&W rule, gaiamock noise:  {n_gm:8,.0f}",
              f"  L&W rule, L&W noise:       {n_lw:8,.0f}   (gaiamock per-transit noise for G < 14: median {med:.0f} "
              "micro-arcsec, against L&W's 54)",
              "  by distance (L&W rule and noise): " + ", ".join(f"{lo}-{hi} pc {n:,.0f}" for (lo, hi), n in
                                                              zip(zip(edges[:-1], edges[1:]), by_d)), ""]
    L += [f"Lammers & Winn predict {S3.LW_YIELD[0]:,} +/- {S3.LW_YIELD[1]:,} for all distances; the remaining difference",
          "comes from the stellar sample (real Gaia DR3 stars here, a Galaxy model in L&W)."]
    C.write_summary(out / "summary.txt", L)


if __name__ == "__main__":
    main()
