# Running on a SLURM cluster

The long runs split into independent parts, so they map onto SLURM job arrays: one array task per
part, each on the cores of one node, and a short merge job at the end. Nothing needs MPI.

| Script | What | Typical use |
|---|---|---|
| `step3_parts.sbatch` + `step3_merge.sbatch` | step 3 (real stars, occurrence rates) split with `--part k/n` | whole sky within 250 pc (~2.3 M stars): ~1,100 core-hours |
| `dr4_fit_parts.sbatch` | seven-model fits of real DR4 epoch astrometry (`scripts/dr4_fit.py`) | all DR4 sources within 250 pc |

Set `CONDA_DIR` and `CONDA_ENV` (or edit the scripts) and adapt `--cpus-per-task` and `--time` to the
cluster's limits. Extra step 3 options go through `STEP3_ARGS`, for example
`--export=ALL,STEP3_ARGS="--variant coplanar,stable --inj-fits 50000"`. Create `logs/` first.
