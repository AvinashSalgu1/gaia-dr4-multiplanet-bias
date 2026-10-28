"""Step 2: what an unmodelled second planet does to Gaia DR4 planet fits.

Simulated two-planet systems (target planet in Gaia's detectable range + a
companion anywhere from 30 days to 54 years), each analysed with the orbit
models Gaia DR4 publishes - one planet, one planet + acceleration, one planet
+ jerk, two planets - and the DR4-like choice between them (lowest BIC),
with and without the companion (identical noise). See src/gedi/step2/.

Run from the repository folder. Costs are rough, in CPU-core hours (divide by the
number of cores for the wall-clock time):

    python scripts/step2.py run -n 3000        # first 3000 systems (~10 core-hours)
    python scripts/step2.py run                # all 20,000 (~70 core-hours; resumes, keeps what is done)
    python scripts/step2.py plot               # summary + figures from what is done so far
    python scripts/step2.py sweeps --baseline outer     # controlled sweeps (~10 core-hours each)
    python scripts/step2.py sweeps --baseline inner
    python scripts/step2.py maps               # best/worst-case maps (~8 core-hours)
    python scripts/step2.py examples           # MCMC corner plots for 3 example stars (~1 core-hour each)

Add --quick to sweeps / maps / examples for a 1-3 minute smoke test.
Everything is written to results/ (step2_*.csv, step2_*.png, step2_summary.txt).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    sub = p.add_subparsers(dest="what", required=True)
    r = sub.add_parser("run", help="population run (resumable), then summary + figures")
    r.add_argument("-n", type=int, default=20000)
    r.add_argument("-j", "--jobs", type=int, default=-1)
    pl = sub.add_parser("plot", help="summary + figures from the population results")
    pl.add_argument("-n", type=int, default=None, help="use only the first n systems")
    for name in ("sweeps", "maps"):
        s = sub.add_parser(name)
        s.add_argument("--baseline", choices=["outer", "inner"], default="outer")
        s.add_argument("--vary-angles", action="store_true")
        s.add_argument("-r", "--repeats", type=int, default=None)
        s.add_argument("--quick", action="store_true")
        s.add_argument("-j", "--jobs", type=int, default=-1)
        s.add_argument("--plot-only", action="store_true")
        if name == "maps":
            s.add_argument("--slice", default="e_c")
            s.add_argument("--slice-values", type=float, nargs="+", default=[0.0, 0.3, 0.6])
    e = sub.add_parser("examples")
    e.add_argument("--quick", action="store_true")
    e.add_argument("--plot-only", action="store_true")
    a = p.parse_args()

    from gedi.step2 import population as POP
    from gedi.step2 import sweeps as SW
    from gedi.step2.plots import analyse_population
    from gedi.step2.sim import RESULTS

    if a.what == "run":
        POP.run(a.n, a.jobs)
        analyse_population(POP.load(a.n))
    elif a.what == "plot":
        if not POP.OUT_CSV.exists():
            raise SystemExit(f"{POP.OUT_CSV} not found - run 'python scripts/step2.py run -n 3000' first")
        analyse_population(POP.load(a.n))
    elif a.what in ("sweeps", "maps"):
        csv_path = RESULTS / ("step2_sweeps_quick.csv" if a.quick else "step2_sweeps.csv")
        if not a.plot_only:
            if a.what == "sweeps":
                reps = a.repeats or (4 if a.quick else 60)
                SW.run(SW.sweep_tasks(a.baseline, a.vary_angles, reps, a.quick), csv_path, a.jobs)
            else:
                if a.slice not in SW.BASELINES["outer"]:
                    raise SystemExit(f"unknown --slice {a.slice}; choose from {sorted(SW.BASELINES['outer'])}")
                reps = a.repeats or (3 if a.quick else 40)
                SW.run(SW.map_tasks(a.baseline, a.vary_angles, reps, a.slice, a.slice_values, a.quick),
                       csv_path, a.jobs)
        SW.analyse(csv_path)
    elif a.what == "examples":
        from gedi.step2 import examples as EX
        if a.plot_only:
            EX.plot()
        else:
            EX.run(quick=a.quick)


if __name__ == "__main__":
    main()
