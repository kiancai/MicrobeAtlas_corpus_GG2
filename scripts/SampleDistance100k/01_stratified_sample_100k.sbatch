#!/bin/bash
#SBATCH -J sd100_sample
#SBATCH -p corexd192
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH -t 01:00:00
#SBATCH --mem 64G
#SBATCH --output=/hpcdisk1/limk_group/caiqy/project/260428_greengene2/logs/sd100_sample_%j.log
#SBATCH --error=/hpcdisk1/limk_group/caiqy/project/260428_greengene2/logs/sd100_sample_%j.log
# Submit with explicit resources:
#   sbatch -c 4 --mem=64G -t 01:00:00 scripts/SampleDistance100k/01_stratified_sample_100k.sbatch

set -euo pipefail

PROJECT_DIR="/hpcdisk1/limk_group/caiqy/project/260428_greengene2"
mkdir -p "$PROJECT_DIR/logs"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate baseBio

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export PYTHONUNBUFFERED=1

cd "$PROJECT_DIR"
python -u scripts/SampleDistance100k/01_stratified_sample_100k.py
