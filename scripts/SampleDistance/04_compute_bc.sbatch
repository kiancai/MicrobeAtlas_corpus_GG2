#!/bin/bash
#SBATCH -J sd_bc
#SBATCH -p corexd192
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH -t 02:00:00
#SBATCH --mem 512G
#SBATCH --output=/hpcdisk1/limk_group/caiqy/project/260428_greengene2/logs/sd_bc_%j.log
#SBATCH --error=/hpcdisk1/limk_group/caiqy/project/260428_greengene2/logs/sd_bc_%j.log
# 提交方式（必须显式带资源 flag，即使 #SBATCH 已写）：
#   sbatch -c 64 --mem=512G -t 02:00:00 scripts/SampleDistance/04_compute_bc.sbatch

set -euo pipefail

PROJECT_DIR="/hpcdisk1/limk_group/caiqy/project/260428_greengene2"
mkdir -p "$PROJECT_DIR/logs"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate baseBio

# multiprocessing.Pool 用 sched_getaffinity 自动检测；显式 N_JOBS 保险
N="${SLURM_CPUS_PER_TASK:-64}"
export N_JOBS="$N"
# 避免 scipy/numpy BLAS 内部多线程跟 multiprocessing 抢 CPU
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# 让 Python print 不缓冲，方便看进度
export PYTHONUNBUFFERED=1

echo "=== 环境信息 ==="
echo "Python:           $(python --version)"
echo "N_JOBS:           $N"
echo "SLURM_JOB_ID:     ${SLURM_JOB_ID:-N/A}"
echo "Host:             $(hostname)"
echo ""

cd "$PROJECT_DIR"
python -u scripts/SampleDistance/04_compute_bc.py
