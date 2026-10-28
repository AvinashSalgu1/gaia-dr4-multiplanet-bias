"""Step 2: two-planet systems analysed as if they had fewer planets.

Modules
-------
sim        exact system draws and data (same random-number order as the
           original step-2 script, so earlier results stay valid)
analyses   the four analyses (one planet, one planet + acceleration,
           two planets, Gaia-style pipeline), each run on the data with the
           companion and on the same data without it
outcomes   which planet each analysis reports and whether it is right
population the resumable 20,000-system run
factors    what else matters, beyond companion period and wobble ratio
sweeps     controlled one-parameter-at-a-time experiments and slice maps
examples   posterior samples (MCMC) and corner plots for three example stars
plots      overview and detail figures, text summary

Entry point: scripts/step2.py
"""
