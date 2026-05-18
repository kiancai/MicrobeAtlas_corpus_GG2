#!/bin/bash
#SBATCH -J sd100_bc
#SBATCH -p corexd192
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH -t 06:00:00
#SBATCH --mem 450G
#SBATCH --output=/hpcdisk1/limk_group/caiqy/project/260428_greengene2/logs/sd100_bc_%j.log
#SBATCH --error=/hpcdisk1/limk_group/caiqy/project/260428_greengene2/logs/sd100_bc_%j.log
# Submit with explicit resources:
#   sbatch -c 64 --mem=450G -t 06:00:00 scripts/SampleDistance100k/03_compute_bc_100k.sbatch

set -euo pipefail

PROJECT_DIR="/hpcdisk1/limk_group/caiqy/project/260428_greengene2"
mkdir -p "$PROJECT_DIR/logs"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate baseBio

N="${SLURM_CPUS_PER_TASK:-64}"
export N_JOBS="$N"
export BATCH="${BATCH:-250}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export PYTHONUNBUFFERED=1

echo "Python: $(python --version)"
echo "N_JOBS: $N_JOBS"
echo "BATCH: $BATCH"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-N/A}"
echo "Host: $(hostname)"

cd "$PROJECT_DIR"
python -u scripts/SampleDistance100k/03_compute_bc_100k.py
