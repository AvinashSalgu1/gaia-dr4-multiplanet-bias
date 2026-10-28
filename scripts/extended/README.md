# Extended analyses

Robustness tests, comparisons and forecasts that build on steps 2 and 3: model selection, noise,
rate inference, the drift of the power law of planets fitted alone, comparisons with Lammers &
Winn (2026) and Greenfield & Owen (2026), stability of the simulated systems, contamination by stellar companions, a census of
known planet hosts, what DR4 can measure, and radial-velocity follow-up. Results go to
`results/extended/<name>/` (`summary.txt`, figures and small tables). The full description and
the results are in [docs/METHODS_AND_RESULTS.md](../../docs/METHODS_AND_RESULTS.md),
section "Extended analyses".

Scripts named `g1_*` reuse existing step 2 or step 3 fits and run in a minute or two; `g2_*`
scripts make new Gaia fits (gaiamock) or integrations (rebound, dynesty). `fig_*` scripts make
the combined figures from the results.

| Script | What it answers | Needs | Time (12 cores) |
|---|---|---|---|
| `g2_step2_exclusive.py` | Step 2 outcomes as categories that add up to 100% | step 2 results | ~1 min |
| `g2_noise_robustness.py` | Do the step 2 results survive correlated, excess or underestimated noise? | gaiamock | 1.5–3 h |
| `g2_evidence.py` | Does the Bayesian evidence choose the same model as the BIC? (`--prior-scale` for other priors) | gaiamock, dynesty | ~1 h per prior scale |
| `g2_evidence_truth.py` | Which choice is right, judged against the true systems? | step 2 results, evidence runs | seconds |
| `g1_rate_binned.py` | Occurrence rate without assuming a power law; injection response matrix | step 3 results | ~1 min |
| `g1_pl_drift.py` | Why does the power-law rate of planets fitted alone drift with multiplicity? | step 3 results | ~1 min |
| `g1_lw_yield.py` | The step 3 planets counted with the Lammers & Winn rule and noise | step 3 results | seconds |
| `g2_go_checks.py` | Comparison with Greenfield & Owen: `resoecc`, `single`, `second` (`--prior-scale`) | gaiamock, dynesty | 1–6 h per part (64 cores) |
| `g1_multiplicity_signal.py` | Does the mix of DR4 solution types constrain how many hosts have a second giant? | step 3 results | ~1 min |
| `g1_export_pairs.py` | Writes the planet pairs for the N-body check (`--all` for every pair) | step 3 results | seconds |
| `g2_stability_nbody.py` | Are the simulated pairs dynamically stable? (`--all`: every pair; superseded by `step3.py --variant stable`) | rebound | 2 min (600 pairs), ~2 h (all) |
| `g1_stable_scan.py` | Multiplicity scan with only the stable systems | step 3 results, N-body run | ~1 min |
| `g1_flags_heldout.py` | Warning flags scored on stars not used to choose them | step 3 results | ~1 min |
| `g2_mixed_flags.py` | Warning flags with brown dwarfs and stellar binaries added | step 3 results, gaiamock | ~1 h |
| `g2_census.py` | Known inner-planet hosts, their DR4 sensitivity to outer giants, forecast detections | gaiamock, internet | ~3 h |
| `g2_observability.py` | Which companions DR4 measures as full orbits and which only as accelerations | gaiamock | ~10 min |
| `g1_rv_grid.py` | RV epochs needed to confirm a Gaia orbit (HARPS, FEROS, CORALIE; La Silla seasons) | nothing | ~10 min |
| `fig_noise_evidence.py`, `fig_census.py`, `fig_stability.py`, `fig_binaries.py` | Combined figures | the results above | seconds |

Every `g2_*` script resumes where it stopped and has `--quick` for a short test run; most also
have `--synthetic` (fake scanning law, code test only). `--root /path/to/step3_results` points the
`g1_*` scripts at a step 3 results folder elsewhere.

## Installation

On top of the main installation (top-level README):

```bash
pip install -r requirements-extended.txt     # rebound, dynesty
pytest tests/test_extended.py
```

## Running

`run_extended.sh` runs any list of steps in order, logs to `results/extended/logs/`, and can
be started again after an interruption:

```bash
bash scripts/extended/run_extended.sh all                 # everything, in a sensible order
bash scripts/extended/run_extended.sh census noise        # chosen steps
```

Steps: `table`, `noise`, `evidence`, `evidence3`, `evidence03`, `truth`, `census`, `stability`,
`stabilityall`, `observability`, `mixed`. The steps are independent apart from `truth` (after
the evidence runs) and the stable-system scan (after `stabilityall`), so they can be split
between machines by giving each machine a different list. Use `tmux` or `nohup` for long runs.

The `g1_*` scripts and the figures are run directly:

```bash
python scripts/extended/g1_rate_binned.py
python scripts/extended/g1_multiplicity_signal.py
python scripts/extended/g1_flags_heldout.py
python scripts/extended/g1_rv_grid.py            # --edge-on for inclination 90 deg
python scripts/extended/g1_export_pairs.py --all
python scripts/extended/g1_stable_scan.py --mask nbody
python scripts/extended/fig_noise_evidence.py
python scripts/extended/fig_census.py
python scripts/extended/fig_stability.py
python scripts/extended/fig_binaries.py
```
