#!/bin/bash
# Seven-model fits of real DR4 epoch astrometry, one array task per part:
#   sbatch --array=0-99 --export=ALL,FILES="/path/to/dr4/epoch/*.ecsv" tools/slurm/dr4_fit_parts.sbatch
# Each part writes results/dr4/fits.part<k>of<N>.csv; concatenate them when all are done.
#SBATCH --job-name=dr4fit
#SBATCH --cpus-per-task=32
#SBATCH --time=72:00:00
#SBATCH --output=logs/dr4fit_%A_%a.out
set -euo pipefail
N_PARTS=${N_PARTS:-100}
source "${CONDA_DIR:-$HOME/miniforge3}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV:-gaia}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 LOKY_MAX_CPU_COUNT=$SLURM_CPUS_PER_TASK
python scripts/dr4_fit.py $FILES --part "$SLURM_ARRAY_TASK_ID/$N_PARTS" -j "$SLURM_CPUS_PER_TASK"
