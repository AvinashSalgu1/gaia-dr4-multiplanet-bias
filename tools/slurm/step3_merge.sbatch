#!/bin/bash
# Joins the parts written by step3_parts.sbatch and makes the summary and figures.
#SBATCH --job-name=step3-merge
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --output=logs/step3_merge_%j.out
set -euo pipefail
source "${CONDA_DIR:-$HOME/miniforge3}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV:-gaia}"
python scripts/step3.py --merge-parts ${STEP3_ARGS:-}
