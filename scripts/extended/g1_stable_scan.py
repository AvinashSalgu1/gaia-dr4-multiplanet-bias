"""Does dynamical instability of the simulated systems change the multiplicity bias?

Repeats the step 3 multiplicity scan (reweighting the existing fits, no new fits) with only the
multi-giant systems that are dynamically plausible, and compares it with the scan over all systems:

    python scripts/extended/g1_stable_scan.py --mask nbody      # needs results/extended/stability_all/pairs.csv
    python scripts/extended/g1_stable_scan.py --mask coplanar   # proxy: every neighbouring pair inclined < 30 deg

nbody: a system is dropped if any of its neighbouring pairs was unstable in the full N-body run
(g2_stability_nbody.py --all: close encounter, ejection or MEGNO > 4 within 1e5 inner orbits).
coplanar: a system is kept only if every neighbouring pair has a mutual inclination below 30 deg
(26-37% of such pairs were unstable in the 600-pair check, against 72-98% above 60 deg); this needs
no N-body run and gives an early answer. Single-giant hosts are always kept. The multiplicity
distribution at each f is then built from the kept systems only.

Output: results/extended/stable_scan_<mask>/summary.txt
"""
from __future__ import annotations

import argparse
import csv

import numpy as np

import common as C


def bad_stars(mask):
    """{scenario: set of star indices whose system is dropped}."""
    out = {"independent": set(), "clustered": set()}
    if mask == "nbody":
        path = C.OUT / "stability_all" / "pairs.csv"
        if not path.exists():
            raise SystemExit(f"{path} not found: run g2_stability_nbody.py --all first")
        for r in csv.DictReader(open(path)):
            if r["status"] == "unstable":
                out[r["scenario"]].add(int(r["star"]))
    else:
        from gaiabias.orbits import mutual_inclination
        for r in csv.DictReader(open(C.DATA / "pairs_all.csv")):
            mi = np.degrees(mutual_inclination(float(r["inc1"]), float(r["Omega1"]), float(r["inc2"]),
                                               float(r["Omega2"])))
            if mi >= 30:
                out[r["scenario"]].add(int(r["star"]))
    return out


def patched_weights(keep):
    """host_weights_by_mass with the multiplicity distribution built from kept hosts only."""
    from gaiabias import multiplicity as MU

    def hw(pool, f, mbin_host, k_min_pool=10):
        mu = MU.mu_for(f)
        kmax = int(pool.k.max())
        q_all = MU.ztp_pmf(np.arange(kmax + 1), mu)
        w = np.zeros(len(pool.k))
        lost, f_real, n_tot = 0.0, 0.0, 0
        for j in np.unique(mbin_host):
            sel = (mbin_host == j) & keep
            N_k = np.bincount(pool.k[sel], minlength=kmax + 1).astype(float)
            usable = N_k >= k_min_pool
            lost = max(lost, float(q_all[~usable].sum() + max(0.0, 1 - q_all.sum())))
            q = np.where(usable, q_all, 0.0)
            q /= q.sum()
            f_real += float(q[2:].sum()) * sel.sum()
            n_tot += int(sel.sum())
            with np.errstate(invalid="ignore", divide="ignore"):
                per = np.where(N_k > 0, q / N_k, 0.0)
            w[sel] = per[pool.k[sel]]
        return w, mu, lost, f_real / max(n_tot, 1)
    return hw


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--root", type=C.Path, default=C.STEP3)
    p.add_argument("--mask", choices=["nbody", "coplanar"], default="nbody")
    p.add_argument("--boot", type=int, default=20)
    a = p.parse_args()
    import step3 as S3
    from gaiabias import multiplicity as MU

    z, meta, jobs, d = C.load_step3(a.root)
    res = {an: S3.analyse_one(z, meta, jobs, d, an) for an in S3.SCAN_METHODS}
    bad = bad_stars(a.mask)
    keep = []
    n_multi, n_drop = {}, {}
    for sc in S3.SCENARIOS:
        npl = res[S3.SCAN_METHODS[0]]["scen"][sc]["hosts"]["npl"]
        hosts = np.nonzero(npl > 0)[0]
        drop = np.array([(npl[s] >= 2) and (int(s) in bad[sc]) for s in hosts])
        keep.append(~drop)
        n_multi[sc], n_drop[sc] = int((npl[hosts] >= 2).sum()), int(drop.sum())
    keep = np.concatenate(keep)
    grid = [0.1, 0.3, 0.6, 0.79]
    scan_all = S3.scan_multiplicity(res, meta, f_grid=grid, n_boot=a.boot)
    original = MU.host_weights_by_mass
    MU.host_weights_by_mass = patched_weights(keep)
    try:
        scan_keep = S3.scan_multiplicity(res, meta, f_grid=grid, n_boot=a.boot)
    finally:
        MU.host_weights_by_mass = original
    out = C.out_dir(f"stable_scan_{a.mask}")
    what = {"nbody": "systems with an unstable pair in the N-body run (1e5 inner orbits) removed",
            "coplanar": "only multi-giant systems whose neighbouring pairs are all inclined < 30 deg kept"}[a.mask]
    L = ["=== MULTIPLICITY SCAN WITH DYNAMICALLY PLAUSIBLE SYSTEMS ONLY ===",
         f"step 3 folder: {a.root}; {what}.",
         "Dropped multi-giant systems: " + "; ".join(f"{sc} {n_drop[sc]} of {n_multi[sc]}" for sc in S3.SCENARIOS) + ".",
         "Inferred / true occurrence rate: well-measured bins (binned) and full giant range (power law); "
         "A = real systems, B = same planets fitted alone. Ranges: bootstrap 16-84%.", ""]
    for an in S3.SCAN_METHODS:
        L.append(f"--- {S3.LABEL[an]} ---")
        L.append(f"   {'f':>5}  {'sample':<10}{'binned A':>22}{'binned B':>10}{'power law A':>22}{'power law B':>13}")
        for i, f in enumerate(grid):
            for lab, s in (("all", scan_all[an]), ("kept", scan_keep[an])):
                r = lambda k: f"{s[k][i]:.2f} [{s[k + '_lo'][i]:.2f},{s[k + '_hi'][i]:.2f}]"  # noqa: E731
                L.append(f"   {f:>5.2f}  {lab:<10}{r('bin_A'):>22}{s['bin_B'][i]:>10.2f}{r('pl_A'):>22}"
                         f"{s['pl_B'][i]:>13.2f}" + (f"   (realised f {s['f_real'][i]:.2f})" if lab == "kept" else ""))
        L.append("")
    if a.mask == "coplanar":
        L.append("The coplanar mask is a proxy (pairs below 30 deg can still be unstable); the N-body mask is the "
                 "real test.")
    C.write_summary(out / "summary.txt", L)


if __name__ == "__main__":
    main()
