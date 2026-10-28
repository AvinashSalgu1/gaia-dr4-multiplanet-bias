"""The step-2 population run: every system, all four analyses, with and
without the companion. Resumable: systems already in the output file are
skipped. If an older results/step2_pairs.csv is present, its one-planet fits are
reused, since the same seed gives exactly the same data.

Output: results/step2_population.csv (one row per system).
"""
from __future__ import annotations

import csv
import time

import numpy as np

from . import analyses as A
from .sim import META_COLS, RESULTS, describe, make_data, population_seeds

OUT_CSV = RESULTS / "step2_population.csv"
COLUMNS = META_COLS + A.all_cols()

# old file -> {new column: old column}. Only the one-planet (gaiamock) fits of the
# original step 2 are reused: the other analyses now store more numbers (BIC chi2
# values, significance) than the old files have, so they are refitted.
_OLD = {
    "step2_pairs.csv": {
        **{f"one_w_{k}": f"b_{o}" for k, o in [("det", "det"), ("P", "P"), ("m", "m"), ("e", "ecc"),
                                               ("dchi2", "dchi2"), ("F2", "F2"), ("a0", "a0"), ("sa0", "sig_a0")]},
        **{f"one_a_{k}": f"s_{o}" for k, o in [("det", "det"), ("P", "P"), ("m", "m"), ("e", "ecc"),
                                               ("dchi2", "dchi2"), ("F2", "F2"), ("a0", "a0"), ("sa0", "sig_a0")]},
    },
}


def _val(c, x):
    if c in A.TEXT_FIELDS:
        return x
    return float(x) if x not in ("", None) else np.nan


def load_old_results() -> dict:
    """seed -> {column: value} from the earlier separate step-2 runs."""
    cache: dict[int, dict] = {}
    for name, mapping in _OLD.items():
        path = RESULTS / name
        if not path.exists():
            continue
        with open(path) as f:
            for r in csv.DictReader(f):
                row = cache.setdefault(int(float(r["seed"])), {})
                for new, old in mapping.items():
                    if old in r and new not in row:
                        row[new] = _val(new, r[old])
    return cache


def missing(row: dict) -> list[tuple[str, str]]:
    """(analysis, version) pairs that still need fitting for this system."""
    out = []
    for a in A.FITTED:
        for v in A.VERSIONS:
            # present = this analysis was run (a failed or empty fit is stored as NaN)
            if not all(c in row for c in A.cols(a, v, required_only=True)):
                out.append((a, v))
    return out


def run_one(seed: int, cached: dict | None = None) -> dict:
    s, ep, single, both, noise = make_data(seed)
    row = describe(seed, s, ep, single, both)
    row.update(cached or {})
    todo = missing(row)
    for v, obs in (("w", both.obs + noise), ("a", single.obs + noise)):
        which = {a for a, vv in todo if vv == v}
        if which:
            row.update(A.run(ep, obs, both.err, s["m_star"], v, which))
    for c in COLUMNS:
        row.setdefault(c, np.nan)
    return row


def _read_existing() -> tuple[list[dict], dict]:
    """Rows already in the output file: (complete rows, {seed: values} of incomplete ones).

    A row is incomplete when a newer version of the code stores more numbers
    than it has (e.g. the chi2 values needed for the BIC choice); those
    systems are re-run, reusing everything that is already there.
    """
    if not OUT_CSV.exists():
        return [], {}
    complete, partial = [], {}
    with open(OUT_CSV) as f:
        for r in csv.DictReader(f):
            vals = {k: _val(k, v) for k, v in r.items() if k in COLUMNS and v not in ("", None)}
            if missing(vals):
                partial[int(float(r["seed"]))] = vals
            else:
                complete.append(r)
    return complete, partial


def run(n: int, n_jobs: int = -1, chunk: int = 400):
    from joblib import Parallel, delayed

    RESULTS.mkdir(exist_ok=True)
    seeds = population_seeds(n)
    complete, partial = _read_existing()
    # rewrite the file with complete rows only (current column set); incomplete ones are redone below
    tmp = OUT_CSV.with_suffix(".tmp")
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        w.writeheader()
        for r in complete:
            w.writerow({k: r.get(k, "") for k in COLUMNS})
    tmp.replace(OUT_CSV)
    done = {int(float(r["seed"])) for r in complete}
    todo = [s for s in seeds if s not in done]
    old = load_old_results()
    for s_, vals in partial.items():
        old.setdefault(s_, {}).update(vals)
    reuse = sum(1 for s in todo if s in old)
    print(f"step 2 population: {len(seeds)} systems requested, {len(seeds) - len(todo)} already done, "
          f"{len(todo)} to go ({reuse} of them can reuse earlier fits)", flush=True)
    t0 = time.perf_counter()
    for i in range(0, len(todo), chunk):
        part = todo[i:i + chunk]
        rows = Parallel(n_jobs=n_jobs)(delayed(run_one)(s, old.get(s)) for s in part)
        with open(OUT_CSV, "a", newline="") as f:
            w = csv.DictWriter(f, fieldnames=COLUMNS)
            for r in rows:
                w.writerow({k: r[k] for k in COLUMNS})
        d = i + len(part)
        el = (time.perf_counter() - t0) / 60
        print(f"  {d}/{len(todo)}  ({el:.1f} min, ~{el / d * (len(todo) - d):.0f} min left)", flush=True)
    print(f"saved {OUT_CSV}")


def load(n: int | None = None) -> dict:
    """Population results as numpy arrays (optionally only the first n systems)."""
    with open(OUT_CSV) as f:
        rows = list(csv.DictReader(f))
    if n is not None:
        keep = set(population_seeds(n))
        rows = [r for r in rows if int(float(r["seed"])) in keep]
    return {k: (np.array([r[k] for r in rows]) if k in A.TEXT_FIELDS
                else np.array([_val(k, r[k]) for r in rows], float)) for k in COLUMNS}
