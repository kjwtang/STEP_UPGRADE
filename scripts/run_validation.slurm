#!/bin/bash
# Adjust account, partition, time, memory, input path, and tuning parameters
# for the local HPC site before submission with: sbatch scripts/run_validation.slurm
#SBATCH --job-name=step-smoke
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --output=step-smoke-%j.out

set -euo pipefail

source .venv/bin/activate

python scripts/run_npy_validation.py /path/to/precip.npy \
  --output-dir results/smoke-${SLURM_JOB_ID} \
  --start 0 --stop 24 \
  --threshold 0.6 --bridge-radius 1 --min-size 4 \
  --workers "${SLURM_CPUS_PER_TASK}" --tau 0.35 --max-displacement 20
