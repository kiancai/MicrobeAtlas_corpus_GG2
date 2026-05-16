#!/bin/bash
#SBATCH -J sd_wunifrac
#SBATCH -p corexd192
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH -t 02:00:00
#SBATCH --mem 512G
#SBATCH --output=/hpcdisk1/limk_group/caiqy/project/260428_greengene2/logs/sd_wunifrac_%j.log
#SBATCH --error=/hpcdisk1/limk_group/caiqy/project/260428_greengene2/logs/sd_wunifrac_%j.log
# 提交方式（必须显式带资源 flag，即使 #SBATCH 已写）：
#   sbatch -c 64 --mem=512G -t 02:00:00 scripts/SampleDistance/05_compute_wunifrac.sbatch

set -euo pipefail

PROJECT_DIR="/hpcdisk1/limk_group/caiqy/project/260428_greengene2"
mkdir -p "$PROJECT_DIR/logs"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate baseBio

# 把 SLURM 分到的 CPU 数交给 OpenMP / unifrac
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-64}"
export N_JOBS="$OMP_NUM_THREADS"
export PYTHONUNBUFFERED=1

echo "=== 环境信息 ==="
echo "Python:           $(python --version)"
echo "OMP_NUM_THREADS:  $OMP_NUM_THREADS"
echo "SLURM_JOB_ID:     ${SLURM_JOB_ID:-N/A}"
echo "Host:             $(hostname)"
python -c "import unifrac; print('unifrac:', unifrac.__version__ if hasattr(unifrac, '__version__') else 1.5)"
echo ""

cd "$PROJECT_DIR"
python scripts/SampleDistance/05_compute_wunifrac.py
