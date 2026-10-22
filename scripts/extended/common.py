"""Shared helpers for the extended analyses in scripts/extended/.

Every script writes into results/extended/<name>/ (summary.txt, figures, CSV).
g1_* scripts only read existing step 2 or step 3 results (no gaiamock needed);
g2_* scripts make new Gaia fits and need gaiamock (see README.md).
"""
from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

STEP3 = REPO / "results" / "step3_gaia"
OUT = REPO / "results" / "extended"
DATA = REPO / "data" / "extended"


def out_dir(name: str) -> Path:
    d = OUT / name
    d.mkdir(parents=True, exist_ok=True)
    return d


def load_step3(root: Path = STEP3):
    """(z, meta, jobs, d) from a step 3 results folder."""
    import step3 as S3
    P = S3.paths(root)
    if not P["fits"].exists():
        raise SystemExit(f"no step 3 results in {root} (run scripts/step3.py first, or pass --root)")
    with np.load(P["pop"]) as f:
        z = {k: f[k] for k in f.files}
    meta = json.loads(P["meta"].read_text())
    with open(P["jobs"]) as f:
        jobs = list(csv.DictReader(f))
    return z, meta, jobs, S3.load_fits(root)


def setup_plots():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from gaiabias import plotstyle as ps
    ps.apply()
    return plt, ps


def write_summary(path: Path, lines):
    text = "\n".join(lines) + "\n"
    print(text)
    path.write_text(text)


class Timer:
    def __init__(self):
        self.t0 = time.perf_counter()

    def __str__(self):
        return f"{(time.perf_counter() - self.t0) / 60:.1f} min"


def poisson_interval(mu, q=(0.16, 0.84)):
    """Central interval of a Poisson count with mean mu."""
    from scipy.stats import poisson
    return tuple(float(poisson.ppf(x, mu)) for x in q)
