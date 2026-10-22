"""How many RV epochs confirm a Gaia orbit, for realistic targets, instruments and seasons.

No Gaia fits and no step 3 results needed (a few minutes):

    python scripts/extended/g1_rv_grid.py [--n-mc 200] [--years 2] [--edge-on]

Grid: planet mass (1 and 3 M_J), semi-major axis 1-5 au around a 1 M_sun star, V = 9-12,
stellar jitter 2, 5 and 10 m/s, instruments FEROS, CORALIE and HARPS, a programme of
`--years` years from La Silla. For every grid point, `--n-mc` simulated campaigns:

* the planet's orbit is drawn (eccentricity from Kipping 2013, random phase and argument
  of periastron) with a fixed inclination, because Gaia measures it and targets are
  chosen with it: 60 deg by default (the median of random orientations), or 90 deg with
  --edge-on; the RV semi-amplitude follows from the true mass and inclination;
* Gaia provides the orbit shape with realistic errors (period 1%, periastron time 2% of
  the period, eccentricity 0.03, argument of periastron 5 deg);
* epochs fall on nights when the target (declination -30 deg, random right ascension) is
  above 30 deg altitude within 4 h of midnight at La Silla, with 20% of nights lost to
  weather, either at random among those nights or placed near the predicted RV extremes;
* epochs are scheduled with the Gaia ephemeris propagated over the ~10 years between the
  DR4 reference epoch (2017.5) and the programme (from late 2027), so ephemeris errors
  move the planned epochs away from the true extremes;
* the RVs (truth + noise of the instrument and jitter) are fitted at the Gaia period with a
  free phase and amplitude (sine and cosine terms, plus the second harmonic when Gaia gives
  e > 0.3) and an offset. All terms are linear, so the posterior is Gaussian; the orbit
  counts as confirmed when the improvement over a constant RV corresponds to 5 sigma
  (false-alarm probability 5.7e-7 for the number of added terms).

The result is the share of campaigns that confirm the orbit with 6, 10 and 20 epochs. When the
period is much longer than the programme, part of the orbit is never covered and some
phases cannot be confirmed with any number of new epochs; archival RVs or a longer
baseline are then needed.
Instrument precision (per epoch, before jitter) scales with V from: HARPS 1.0 m/s at V = 8
(floor 0.5), FEROS 5.0 m/s (floor 4.0), CORALIE 6.0 m/s (floor 5.0). These are planning
assumptions to be replaced by the measured performance of each instrument.
"""
from __future__ import annotations

import argparse

import numpy as np

import common as C

INST = {"HARPS": (1.0, 0.5), "FEROS": (5.0, 4.0), "CORALIE": (6.0, 5.0)}
A_AU = [1.0, 2.0, 3.0, 4.0, 5.0]
V_MAG = [9.0, 10.0, 11.0, 12.0]
JITTER = [2.0, 5.0, 10.0]
MASSES = [1.0, 3.0]
N_GRID = np.array([4, 6, 8, 10, 12, 15, 20, 25, 30, 40])
LAT = np.radians(-29.26)
DEC = np.radians(-30.0)
K_CONST = 28.4329


def sigma_inst(V, inst):
    s8, floor = INST[inst]
    return np.hypot(s8 * 10 ** (0.2 * (V - 8)), floor)


def observable_nights(ra_deg, years, rng, weather=0.2):
    """Night numbers (days from programme start) when the target can be observed."""
    days = np.arange(int(365.25 * years))
    ra_sun = 360.0 * ((days + 280) % 365.25) / 365.25          # programme starts ~1 October
    ha_mid = (ra_sun + 180.0 - ra_deg + 540.0) % 360.0 - 180.0  # target hour angle at local midnight [deg]
    alt_min = np.radians(30.0)
    cos_h = (np.sin(alt_min) - np.sin(LAT) * np.sin(DEC)) / (np.cos(LAT) * np.cos(DEC))
    ha_lim = np.degrees(np.arccos(np.clip(cos_h, -1, 1)))
    ok = np.abs(ha_mid) <= ha_lim + 60.0
    ok &= rng.uniform(0, 1, len(days)) > weather
    return days[ok].astype(float)


def kepler_rv(t, P, Tp, e, om):
    from gaiabias.orbits import solve_kepler
    M = 2 * np.pi * (t - Tp) / P
    E = solve_kepler(M, e)
    nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
    return np.cos(nu + om) + e * np.cos(om)


def campaign(rng, m, a_au, V, jit, inst, years, n_list, strategy, inc):
    """Fraction-confirmed flags for each N in n_list, for one simulated system."""
    from gaiabias.orbits import draw_ecc_kipping
    P = 365.25 * a_au ** 1.5
    e = float(draw_ecc_kipping(rng))
    om, Tp = rng.uniform(0, 2 * np.pi), rng.uniform(0, P)
    K = K_CONST * m * np.sin(inc) * (P / 365.25) ** (-1 / 3) / np.sqrt(1 - e ** 2)
    # Gaia's version of the orbit shape
    Pg = P * (1 + 0.01 * rng.standard_normal())
    Tg = Tp + 0.02 * P * rng.standard_normal()
    eg = float(np.clip(e + 0.03 * rng.standard_normal(), 0, 0.95))
    og = om + np.radians(5) * rng.standard_normal()
    nights = observable_nights(rng.uniform(0, 360), years, rng)
    sig = np.hypot(sigma_inst(V, inst), jit)
    out = np.zeros(len(n_list), bool)
    if len(nights) < 2:
        return out
    t_off = 10.3 * 365.25                                       # DR4 reference epoch (2017.5) to programme start
    shape_pred = kepler_rv(nights + t_off, Pg, Tg, eg, og)
    n_harm = 2 if eg > 0.3 else 1
    from scipy.stats import chi2 as chi2_dist
    thr = chi2_dist.isf(5.733e-7, 2 * n_harm)                  # 5 sigma for 2 x n_harm added terms
    for i, n in enumerate(n_list):
        n_eff = min(n, len(nights))
        if n_eff < 2 * n_harm + 2:
            continue
        if strategy == "random":
            pick = np.sort(rng.choice(len(nights), n_eff, replace=False))
        else:   # half near the predicted maximum, half near the minimum
            order_hi, order_lo = np.argsort(-shape_pred), np.argsort(shape_pred)
            k_hi = n_eff // 2
            cand_hi = order_hi[:max(3 * k_hi, 1)]
            cand_lo = order_lo[:max(3 * (n_eff - k_hi), 1)]
            pick = np.unique(np.concatenate([rng.choice(cand_hi, min(k_hi, len(cand_hi)), replace=False),
                                             rng.choice(cand_lo, min(n_eff - k_hi, len(cand_lo)), replace=False)]))
        t = nights[pick] + t_off
        y = K * kepler_rv(t, P, Tp, e, om) + sig * rng.standard_normal(len(t))
        cols = [np.ones(len(t))]
        for h in range(1, n_harm + 1):
            cols += [np.cos(2 * np.pi * h * t / Pg), np.sin(2 * np.pi * h * t / Pg)]
        X = np.column_stack(cols)
        beta, *_ = np.linalg.lstsq(X / sig, y / sig, rcond=None)
        chi2_full = np.sum(((y - X @ beta) / sig) ** 2)
        chi2_const = np.sum(((y - y.mean()) / sig) ** 2)
        out[i] = chi2_const - chi2_full > thr
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--n-mc", type=int, default=200)
    p.add_argument("--edge-on", action="store_true", help="inclination 90 deg instead of 60 deg")
    p.add_argument("--years", type=float, default=2.0)
    a = p.parse_args()
    inc = np.radians(90.0 if a.edge_on else 60.0)
    out = C.out_dir("rv_grid_edge_on" if a.edge_on else "rv_grid")
    rng = np.random.default_rng(42)
    rows = []
    for strategy in ("extremes", "random"):
        for inst in INST:
            for jit in JITTER:
                for m in MASSES:
                    for a_au in A_AU:
                        for V in V_MAG:
                            ok = np.array([campaign(rng, m, a_au, V, jit, inst, a.years, N_GRID, strategy, inc)
                                           for _ in range(a.n_mc)])
                            frac = ok.mean(axis=0)
                            n90 = N_GRID[np.argmax(frac >= 0.9)] if (frac >= 0.9).any() else np.inf
                            n50 = N_GRID[np.argmax(frac >= 0.5)] if (frac >= 0.5).any() else np.inf
                            rows.append(dict(strategy=strategy, inst=inst, jitter=jit, m=m, a=a_au, V=V,
                                             n90=n90, n50=n50, frac10=frac[list(N_GRID).index(10)],
                                             frac6=frac[list(N_GRID).index(6)],
                                             frac20=frac[list(N_GRID).index(20)], curve=frac))
    import csv
    with open(out / "rv_grid.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["strategy", "instrument", "jitter_ms", "mass_MJ", "a_au", "V", "epochs_90pct", "epochs_50pct",
                    "confirmed_with_6", "confirmed_with_10", "confirmed_with_20"])
        for r in rows:
            w.writerow([r["strategy"], r["inst"], r["jitter"], r["m"], r["a"], r["V"], r["n90"], r["n50"],
                        f"{r['frac6']:.2f}", f"{r['frac10']:.2f}", f"{r['frac20']:.2f}"])
    L = ["=== RV CAMPAIGNS TO CONFIRM A GAIA ORBIT (5 sigma on the amplitude) ===",
         f"{a.n_mc} simulated campaigns per point; {a.years:g}-yr programme from La Silla; target at dec -30; "
         f"inclination {np.degrees(inc):.0f} deg; Gaia orbit shape with realistic errors.",
         "Entries: share of campaigns confirming the orbit with 6 / 10 / 20 epochs.", ""]
    for strategy in ("extremes", "random"):
        for m in MASSES:
            L.append(f"--- {m:g} M_J planet, 1 M_sun star, epochs "
                     f"{'near the predicted RV extremes' if strategy == 'extremes' else 'at random in the season'} ---")
            L.append("   instrument jitter  a[au] " + "".join(f"{'V=' + format(V, 'g'):>16}" for V in V_MAG))
            for inst in INST:
                for jit in JITTER:
                    for a_au in A_AU:
                        sel = [r for r in rows if r["strategy"] == strategy and r["inst"] == inst and r["jitter"] == jit
                               and r["m"] == m and r["a"] == a_au]
                        L.append(f"   {inst:<9} {jit:4.0f} m/s  {a_au:4.0f}  " + "".join(
                            f"{r['frac6']:>6.2f}/{r['frac10']:.2f}/{r['frac20']:.2f}" for r in sel))
            L.append("")
    sel = [r for r in rows if r["strategy"] == "extremes" and r["m"] == 1.0 and r["inst"] in ("FEROS", "CORALIE")]
    worst = [r for r in sel if r["frac10"] < 0.9]
    L.append(f"1 M_J with FEROS or CORALIE, epochs near the extremes: 10 epochs confirm fewer than 90% of orbits in "
             f"{len(worst)} of {len(sel)} grid points.")
    C.write_summary(out / "summary.txt", L)

    plt, ps = C.setup_plots()
    fig, axes = plt.subplots(2, 3, figsize=(16, 8.6), sharex=True, sharey=True, constrained_layout=True)
    cmap = plt.get_cmap("viridis")
    for r_, m in enumerate(MASSES):
        for c_, inst in enumerate(INST):
            ax = axes[r_, c_]
            for k, V in enumerate(V_MAG):
                for jit, ls in ((5.0, "-"), (10.0, "--")):
                    y = [next(r["frac10"] for r in rows if r["strategy"] == "extremes" and r["inst"] == inst and
                              r["jitter"] == jit and r["m"] == m and r["a"] == a_au and r["V"] == V) for a_au in A_AU]
                    ax.plot(A_AU, y, ls, marker="o", ms=3.5, color=cmap(k / (len(V_MAG) - 0.5)),
                            label=f"V = {V:g}, jitter {jit:g} m/s" if (r_, c_) == (0, 0) else None)
            ax.axhline(0.9, color="0.6", lw=0.8, ls=":")
            ax.set_title(f"{inst}, {m:g} M_J planet", fontsize=10.5)
            ax.set_ylim(0, 1.02)
            if r_ == 1:
                ax.set_xlabel("semi-major axis (au), 1 M_sun star")
            if c_ == 0:
                ps.percent_axis(ax, "y", "orbits confirmed at 5 sigma\nwith 10 RV epochs")
    axes[0, 0].legend(fontsize=7.5, ncol=2, loc="lower left")
    fig.suptitle(f"Confirming Gaia orbits from La Silla with 10 epochs ({a.years:g}-yr programme, epochs near the "
                 f"predicted extremes, inclination {np.degrees(inc):.0f} deg)", fontsize=11.5)
    fig.savefig(out / "rv_grid.png", dpi=150, bbox_inches="tight")
    print(f"saved {out / 'rv_grid.png'}")


if __name__ == "__main__":
    main()
