#!/bin/bash
#SBATCH -J sd100_pcoa
#SBATCH -p corexd192
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH -t 06:00:00
#SBATCH --mem 450G
#SBATCH --output=/hpcdisk1/limk_group/caiqy/project/260428_greengene2/logs/sd100_pcoa_%j.log
#SBATCH --error=/hpcdisk1/limk_group/caiqy/project/260428_greengene2/logs/sd100_pcoa_%j.log
# Submit with explicit resources:
#   sbatch -c 64 --mem=450G -t 06:00:00 scripts/SampleDistance100k/05_pcoa_100k.sbatch

set -euo pipefail

PROJECT_DIR="/hpcdisk1/limk_group/caiqy/project/260428_greengene2"
mkdir -p "$PROJECT_DIR/logs"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate baseBio

N="${SLURM_CPUS_PER_TASK:-64}"
export OMP_NUM_THREADS="$N"
export OPENBLAS_NUM_THREADS="$N"
export MKL_NUM_THREADS="$N"
export PCOA_K="${PCOA_K:-10}"
export PCOA_EXTRA="${PCOA_EXTRA:-5}"
export SYMM_BATCH="${SYMM_BATCH:-2000}"
export AUDIT_N="${AUDIT_N:-3000}"
export PYTHONUNBUFFERED=1

echo "Python: $(python --version)"
echo "OMP_NUM_THREADS: $OMP_NUM_THREADS"
echo "PCOA_K: $PCOA_K"
echo "SYMM_BATCH: $SYMM_BATCH"
echo "AUDIT_N: $AUDIT_N"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-N/A}"
echo "Host: $(hostname)"

cd "$PROJECT_DIR"
python -u scripts/SampleDistance100k/05_pcoa_100k.py
