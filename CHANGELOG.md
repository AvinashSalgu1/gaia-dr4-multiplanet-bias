# Changelog

## 0.3.0 (27 October 2026)

Whole-sky analyses, alternative planet populations, comparisons with published work, the input
pipeline for real DR4 data, and templates for running at scale. The fitting code is unchanged;
the default stellar sample and output folders changed (see below), so 0.2.0 numbers are
reproduced only with `--circles`.

### Changed
- The project is renamed GEDI (Gaia Exoplanet Detection and Inference): the Python package `gaiabias` is
  now `gedi` (`import gedi`; `pip install -e .` again after updating), and the repository is now `gedi`.
- Step 3 uses the whole sky within 250 pc by default: the Gaia DR3 stars are downloaded in 768
  HEALPix level-3 tiles, each queried as a range of Gaia source_ids, and every selected star is
  fitted (`--stars` default 0 = all; 2,311,408 stars). Results go to `results/step3_allsky/`.
  Whole-sky numbers are counted directly instead of scaled up from a 25% sample. `--circles`
  reproduces the earlier 25% sample of 3° circles in `results/step3_gaia/`.
- Extended scripts read the whole-sky step 3 run by default (`--root` for another folder).
- `g2_evidence.py` keeps all cores busy (tasks streamed as they finish).

### Added
- Population variants in step 3 (`--variant`): `coplanar` (mutual inclinations from a Rayleigh
  distribution with σ = 2°), `stable` (every multi-planet system integrated with rebound/WHFast
  for 10^4 inner orbits and redrawn if MEGNO shows chaos), `lowecc` (eccentricities of planets in
  multi-planet systems from a Rayleigh distribution with σ = 0.1). Default draws are unchanged.
- `step3.py --build-only`: writes stars, planets and the job list, for clusters.
- `scripts/extended/g1_pl_drift.py`: why the power-law rate of planets fitted alone drifts with
  multiplicity (population-averaged against S/N-aware completeness).
- `scripts/extended/g1_lw_yield.py`: the step 3 planets counted with the Lammers & Winn rule.
- `scripts/extended/g2_go_checks.py`: comparison with Greenfield & Owen (2026) — eccentricity
  and mass of one-orbit fits to circular 2:1 pairs, the single-planet detection threshold, and
  the threshold for a second planet near and away from 2:1.
- `src/gaiabias/dr4.py` and `scripts/dr4_fit.py`: real DR4 epoch astrometry (pre-release
  VOTable; ECSV, CSV or FITS through astropy) reduced to per-transit positions and fitted with the
  seven models; resumable and splittable. `data/dr4_prerelease/` holds the 12-source pre-release.
- `tools/slurm/`: SLURM job arrays for step 3 and the DR4 fits.
- `tools/remote_runs/`: template for running jobs on a shared server without a queue through a
  GitHub Actions runner and a private repository.
- Tests for the variants, the N-body screen and the DR4 pipeline.

### Removed
- The pairwise stability figure from the documentation; the `stable` variant replaces the
  post-hoc check (the scripts remain).

## 0.2.0 (22 October 2026)

Robustness tests of the step 2 and step 3 results, contamination by stellar companions, a census of
known planet hosts, and forecasts for follow-up. The simulation, fitting and step 1–4 code in
`src/gaiabias/` and `scripts/` is unchanged, so earlier results reproduce exactly.

### Added
- `scripts/extended/`, with its own README and a runner (`run_extended.sh`) that resumes after
  interruptions and can split steps between machines:
  - `g2_step2_exclusive.py`: step 2 outcomes as mutually exclusive categories (both planets
    correct, target only, companion only, correct plus wrong, wrong only, acceleration or jerk
    only, nothing reported), for all systems and per kind of companion.
  - `g2_noise_robustness.py`: step 2 refitted under correlated CCD errors, errors shared between
    transits, magnitude-dependent excess noise and underestimated uncertainties.
  - `g2_evidence.py`: Bayesian evidence of the seven models by nested sampling (dynesty), with
    the linear parameters integrated analytically; `--prior-scale` for prior sensitivity.
    `g2_evidence_truth.py` scores the BIC and evidence choices against the true systems.
  - `g1_rate_binned.py`: occurrence rate with one free rate per mass-period cell and the
    injection response matrix.
  - `g1_multiplicity_signal.py`: how the mix of DR4 solution types depends on the fraction of
    hosts with a second giant, and how well that fraction can be measured.
  - `g1_export_pairs.py` and `g2_stability_nbody.py`: N-body integrations (rebound, WHFast,
    MEGNO) of simulated planet pairs, for a sample or (`--all`) every neighbouring pair.
    `g1_stable_scan.py` repeats the multiplicity scan with only the stable systems.
  - `g1_flags_heldout.py`: warning flags chosen on half of the stars and scored on the other
    half (recall, precision, false-positive rate, calibration by number of flags).
  - `g2_mixed_flags.py`: brown dwarfs and stellar binaries added to the real-star simulation,
    with the photocentre reduced by the companion's light and the DR3 RUWE selection applied;
    the flags scored for planet-only and mixed populations.
  - `g2_census.py`: census of known inner-planet hosts (NASA Exoplanet Archive) with each host's
    DR4 sensitivity to outer companions and forecast detections per class.
  - `g2_observability.py`: which companions DR4 measures as full orbits and which only as
    accelerations, for Sun-like and M-dwarf hosts at 25, 50 and 100 pc.
  - `g1_rv_grid.py`: simulated RV campaigns from La Silla (HARPS, FEROS, CORALIE) with seasonal
    visibility, stellar jitter and a free-phase confirmation test.
  - `fig_noise_evidence.py`, `fig_census.py`, `fig_stability.py`, `fig_binaries.py`: combined
    figures.
- `requirements-extended.txt` (rebound, dynesty) and `tests/test_extended.py`: the linear
  evidence integral against brute force, the evidence model choice on a clear signal, RV
  campaigns, census forecasts, the binned rate, observable nights and the binary photocentre.
- Documentation: a section "Extended analyses (version 0.2.0)" in
  `docs/METHODS_AND_RESULTS.md` with the results and nine new figures in `docs/figures/ext_*.png`.

### Changed
- Step 2 results table in `docs/METHODS_AND_RESULTS.md`: now in mutually exclusive categories, so
  each row adds up to 100%. Version 0.1.0 listed target correct, companion correct and wrong
  planet as overlapping shares, so the two-orbit and BIC rows added up to more than 100%. The
  numbers themselves are unchanged.
- Caveats and possible extensions updated with the new results: most simulated pairs are
  dynamically unstable but the bias in the well-measured region is unchanged; near-twin stellar
  binaries can be reported as planets; the BIC adopts the two-orbit model too easily compared
  with the Bayesian evidence. New extensions: a stable and coplanar reference population, binaries
  from Gaia DR3 non-single-star solutions, and planets in binaries.
- README: the extended analyses with their expected plots, and a link to this changelog.
- Version 0.1.0 → 0.2.0.

### Main new results
- The detections lost to an unmodelled companion (about 13 percentage points) are the same under
  all five noise models; the benefit of DR4's alternative models drops from 17 to 8 points when
  the uncertainties are underestimated by 30%.
- The BIC and the Bayesian evidence agree for single orbits and do equally well overall against
  the truth, but the BIC adopts two orbits too often.
- Without an assumed power law, multiplicity changes the rate in the measurable cells (3–13 M_J)
  by 1–6%; 62–68% of the simulated giants lie in cells where no rate can be measured.
- The share of non-plain-orbit DR4 solutions rises from 5% to 25% with multiplicity, enough to
  constrain it to about ±0.02–0.03 with 3,000 planets (statistical only).
- 72% of simulated neighbouring pairs are unstable, mostly inclined ones; the bias in the
  well-measured region is unchanged when only stable systems are kept.
- On held-out stars the flags keep a recall of 72–75% and a precision of 48%. With stellar
  binaries added, near-twin binaries make up 8–10% of reported planets and 40% of them pass every
  flag.
- About 5, 4 and 33 outer giants are expected around known hot-Jupiter, warm-Jupiter and
  small-planet hosts within 200 pc at the field rate; Saturn-mass companions are out of reach
  for DR4 alone.
- Beyond about 5 au DR4 alone measures only accelerations, not full orbits.

## 0.1.0

First public version: simulation of Gaia DR4 epoch astrometry with any number of companions,
DR4-like orbit fits and model choice, steps 1–4 (single-planet validation, controlled two-planet
experiment, real Gaia stars with occurrence rates, yield and RUWE, warning flags and RV
follow-up) and the DR4 pre-release check.
