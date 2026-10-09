"""Smoke test and timing, on any machine with gaiamock installed.

Usage (from the repository folder, with the package installed):
    python scripts/check_setup.py

Prints: whether gaiamock and the scan files are found, one single-planet
injection/recovery, one two-planet system fitted with a one-planet model,
and the time per fit (needed to size the full injection-recovery runs).
"""
import os
import time

import numpy as np

from gaiabias import epochs as E
from gaiabias import gaiamock_bridge as gb
from gaiabias.fitting import fit_one_companion, lammers_winn_detected
from gaiabias.orbits import Companion
from gaiabias.simulate import Star, simulate, snr1


def main():
    print(f"gaiamock dir : {gb.gaiamock_dir()}")
    gb.gaiamock()
    gb.c_funcs()
    n = gb.n_scan_files()
    print(f"scan files   : {n} / 49152")
    if n < 49152:
        raise SystemExit("healpix_scans/ incomplete - unzip healpix_scans.zip inside gaiamock/")

    rng = np.random.default_rng(42)
    ra, dec = 120.0, -35.0
    star = Star(parallax_mas=40.0, m_star_msun=0.7, G_mag=10.5, pmra=50.0, pmdec=-80.0)

    print("\n--- 1 planet injected, 1-planet fit ---")
    ep = E.from_gaiamock(ra, dec, "dr4", rng=rng)
    c = Companion(m_mjup=5.0, P_days=900.0, ecc=0.2, inc=1.0, Omega=0.5, omega=1.5, Tp_days=200.0)
    sim = simulate(ep, star, [c], rng=rng)
    t0 = time.perf_counter()
    fit = fit_one_companion(ep, sim.obs, sim.err, star.m_star_msun)
    dt = time.perf_counter() - t0
    print(f"N transits {len(ep)}, SNR1 = {snr1(sim.a0_mas[0], star.G_mag):.1f}")
    print(f"P   true {c.P_days:8.1f}  fit {fit.P_days:8.1f} +/- {fit.sigma_P:.1f} d")
    print(f"m   true {c.m_mjup:8.2f}  fit {fit.m_mjup:8.2f} M_J")
    print(f"dchi2 {fit.delta_chi2:.0f}  detected (L&W): {lammers_winn_detected(fit)}  [{dt:.2f} s]")

    print("\n--- 2 planets injected, 1-planet fit ---")
    outer = Companion(m_mjup=3.0, P_days=3000.0, ecc=0.1, inc=0.6, Omega=2.0, omega=0.3, Tp_days=900.0)
    sim2 = simulate(ep, star, [c, outer], rng=rng)
    fit2 = fit_one_companion(ep, sim2.obs, sim2.err, star.m_star_msun)
    print(f"inner P {c.P_days:.0f} d, m {c.m_mjup} M_J | outer P {outer.P_days:.0f} d, m {outer.m_mjup} M_J")
    print(f"fit   P {fit2.P_days:.1f} d, m {fit2.m_mjup:.2f} M_J, e {fit2.ecc:.2f}, "
          f"dchi2 {fit2.delta_chi2:.0f}, detected: {lammers_winn_detected(fit2)}")
    print(f"reduced chi2 of 1-planet fit: {fit2.chi2 / (len(ep) - 12):.2f}  (1.0 = adequate)")

    print("\n--- timing: 20 fits ---")
    t0 = time.perf_counter()
    for _ in range(20):
        sim = simulate(ep, star, [c], rng=rng)
        fit_one_companion(ep, sim.obs, sim.err, star.m_star_msun, uncertainties=True)
    per = (time.perf_counter() - t0) / 20
    ncpu = os.cpu_count() or 1
    print(f"{per:.2f} s per fit on one core -> ~{1e5 * per / 3600 / ncpu:.1f} h for 1e5 fits on the "
          f"{ncpu} cores of this machine")


if __name__ == "__main__":
    main()
