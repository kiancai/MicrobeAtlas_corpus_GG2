#!/bin/bash
#SBATCH -J sd100_subset
#SBATCH -p corexd192
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH -t 02:00:00
#SBATCH --mem 128G
#SBATCH --output=/hpcdisk1/limk_group/caiqy/project/260428_greengene2/logs/sd100_subset_%j.log
#SBATCH --error=/hpcdisk1/limk_group/caiqy/project/260428_greengene2/logs/sd100_subset_%j.log
# Submit with explicit resources:
#   sbatch -c 8 --mem=128G -t 02:00:00 scripts/SampleDistance100k/02_build_subset_anndata_100k.sbatch

set -euo pipefail

PROJECT_DIR="/hpcdisk1/limk_group/caiqy/project/260428_greengene2"
mkdir -p "$PROJECT_DIR/logs"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate baseBio

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export OPENBLAS_NUM_THREADS="$OMP_NUM_THREADS"
export MKL_NUM_THREADS="$OMP_NUM_THREADS"
export PYTHONUNBUFFERED=1

cd "$PROJECT_DIR"
python -u scripts/SampleDistance100k/02_build_subset_anndata_100k.py
