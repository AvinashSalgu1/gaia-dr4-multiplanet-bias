#!/bin/bash
# Step 3 split over a SLURM job array: one array task per part, then a merge job.
#   sbatch --array=0-19 tools/slurm/step3_parts.sbatch            # 20 parts
#   sbatch --dependency=afterok:<array job id> tools/slurm/step3_merge.sbatch
# Run from the repository folder, after building the population once with the same options:
#   python scripts/step3.py --build-only [STEP3_ARGS]
#SBATCH --job-name=step3
#SBATCH --cpus-per-task=32          # cores per part
#SBATCH --time=72:00:00             # within the cluster's limit
#SBATCH --output=logs/step3_%A_%a.out
set -euo pipefail
N_PARTS=${N_PARTS:-20}
module purge 2>/dev/null || true
source "${CONDA_DIR:-$HOME/miniforge3}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV:-gaia}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 LOKY_MAX_CPU_COUNT=$SLURM_CPUS_PER_TASK
python scripts/step3.py --part "$SLURM_ARRAY_TASK_ID/$N_PARTS" --jobs "$SLURM_CPUS_PER_TASK" ${STEP3_ARGS:-}
