#!/bin/bash
#SBATCH -J sd_pcoa
#SBATCH -p corexd192
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH -t 02:00:00
#SBATCH --mem 512G
#SBATCH --output=/hpcdisk1/limk_group/caiqy/project/260428_greengene2/logs/sd_pcoa_%j.log
#SBATCH --error=/hpcdisk1/limk_group/caiqy/project/260428_greengene2/logs/sd_pcoa_%j.log
# 提交方式（必须显式带资源 flag，即使 #SBATCH 已写）：
#   sbatch -c 64 --mem=512G -t 02:00:00 scripts/SampleDistance/06_pcoa.sbatch

set -euo pipefail

PROJECT_DIR="/hpcdisk1/limk_group/caiqy/project/260428_greengene2"
mkdir -p "$PROJECT_DIR/logs"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate baseBio

# BLAS 线程数：用 SLURM 分配的 CPU
N="${SLURM_CPUS_PER_TASK:-64}"
export OMP_NUM_THREADS="$N"
export OPENBLAS_NUM_THREADS="$N"
export MKL_NUM_THREADS="$N"
export PYTHONUNBUFFERED=1

echo "=== 环境信息 ==="
echo "Python:           $(python --version)"
echo "OMP_NUM_THREADS:  $OMP_NUM_THREADS"
echo "SLURM_JOB_ID:     ${SLURM_JOB_ID:-N/A}"
echo "Host:             $(hostname)"
echo ""

cd "$PROJECT_DIR"
python scripts/SampleDistance/06_pcoa.py
