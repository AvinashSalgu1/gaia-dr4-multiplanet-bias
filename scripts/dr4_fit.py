"""Fit the DR4 orbit models of this project to REAL Gaia DR4 epoch astrometry.

    python scripts/dr4_fit.py FILE [FILE ...] [--masses masses.csv] [--part k/n] [-j N] [--out DIR]

FILE: epoch-astrometry files with the DR4 data-model columns (the pre-release VOTable now; the
bulk DR4 files once released, see src/gaiabias/dr4.py). Every source is reduced to per-transit
along-scan positions and fitted with the seven models of step 2 (no planet, acceleration, jerk,
one orbit, orbit + acceleration, orbit + jerk, two orbits); the lowest BIC is adopted and its
planets are reported with the project's detection test, as in the simulations.

Host masses (needed only to turn orbit sizes into companion masses): --masses CSV with columns
source_id, m_star; otherwise the mass of the source in the whole-sky DR3 star file (FLAME or the
main-sequence relation, data/gaia_dr3_250pc/stars_tiles.npz); otherwise the main-sequence
relation applied to G and the parallax of the 5-parameter fit.

Resumable: results are appended to OUT/fits.csv (one row per source) as sources finish.
--part k/n fits only sources with index k, k+n, ... (one cluster job per part); the part
files (fits.part<k>of<n>.csv) can simply be concatenated.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def host_masses(sids, masses_csv=None):
    out = {}
    if masses_csv:
        for r in csv.DictReader(open(masses_csv)):
            out[int(r["source_id"])] = float(r["m_star"])
    from gaiabias import gaia_stars as GS
    p = GS.merged_path("tiles")
    if p.exists():
        sel, _ = GS.select(dict(np.load(p)), ruwe_max=None)
        idx = {int(s): i for i, s in enumerate(sel["source_id"])}
        for s in sids:
            if s not in out and int(s) in idx:
                out[int(s)] = float(sel["mass"][idx[int(s)]])
    return out


def one(args):
    sid, ep, y, err, ncc, G, m_star = args
    from gaiabias import dr4
    from gaiabias.fitting import fit_5par
    from gaiabias.gaia_stars import mass_from_MG
    chi2_5, mu = fit_5par(ep, y, err)
    plx = mu[4]
    source_of_mass = "given"
    if m_star is None:
        m_star = float(mass_from_MG(np.array([G + 5 * np.log10(plx / 100.0)]))[0]) if plx > 0 else 1.0
        source_of_mass = "M_G relation"
    row = dict(source_id=sid, G=G, parallax_5p=plx, chi2_5p=chi2_5, m_star=m_star, m_star_from=source_of_mass,
               n_ccd_median=float(np.median(ncc)))
    row.update(dr4.fit_source(ep, y, err, m_star))
    return row


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("files", nargs="+")
    p.add_argument("--masses", default=None)
    p.add_argument("--part", default=None)
    p.add_argument("-j", "--jobs", type=int, default=-1)
    p.add_argument("--excess", action="store_true", help="add agis_source_excess_noise to the CCD errors")
    p.add_argument("--out", type=Path, default=ROOT / "results" / "dr4")
    a = p.parse_args()
    from joblib import Parallel, delayed
    from gaiabias import dr4
    a.out.mkdir(parents=True, exist_ok=True)
    k, n = (int(x) for x in a.part.split("/")) if a.part else (0, 1)
    path = a.out / ("fits.csv" if n == 1 else f"fits.part{k}of{n}.csv")
    done = {int(r["source_id"]) for r in csv.DictReader(open(path))} if path.exists() else set()
    tasks = []
    for f in a.files:
        d = dr4.load_ccd_table(f)
        sids = [int(s) for s in dr4.source_ids(d)]
        masses = host_masses(sids, a.masses)
        for i, s in enumerate(sids):
            if i % n != k or s in done:
                continue
            ep, y, err, ncc, G = dr4.source_transits(d, s, excess=a.excess)
            if len(ep) < 20:
                continue
            tasks.append((s, ep, y, err, ncc, G, masses.get(s)))
    print(f"{len(done)} sources done, {len(tasks)} to fit", flush=True)
    w = fh = None
    for i, row in enumerate(Parallel(n_jobs=a.jobs, return_as="generator_unordered")(delayed(one)(t) for t in tasks), 1):
        if w is None:
            new = not path.exists()
            fh = open(path, "a", newline="")
            w = csv.DictWriter(fh, fieldnames=list(row), extrasaction="ignore", restval="")
            if new:
                w.writeheader()
        w.writerow(row)
        fh.flush()
        if i % 100 == 0 or i == len(tasks):
            print(f"  {i}/{len(tasks)}", flush=True)
    if fh:
        fh.close()
    rows = list(csv.DictReader(open(path))) if path.exists() else []
    from collections import Counter
    L = [f"=== DR4 FITS: {len(rows)} sources ===", "adopted model (lowest BIC):"]
    for mname, c in Counter(r["model"] for r in rows).most_common():
        L.append(f"  {mname:<28}{c:>8}")
    L.append(f"sources with a reported planet: {sum(int(r['n_planets'] or 0) > 0 for r in rows)}")
    (a.out / "summary.txt").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
