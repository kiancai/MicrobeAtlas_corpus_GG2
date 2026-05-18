#!/bin/bash
#SBATCH -J sd100_wuni
#SBATCH -p corexd192
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH -t 04:00:00
#SBATCH --mem 450G
#SBATCH --output=/hpcdisk1/limk_group/caiqy/project/260428_greengene2/logs/sd100_wunifrac_%j.log
#SBATCH --error=/hpcdisk1/limk_group/caiqy/project/260428_greengene2/logs/sd100_wunifrac_%j.log
# Submit with explicit resources:
#   sbatch -c 64 --mem=450G -t 04:00:00 scripts/SampleDistance100k/04_compute_wunifrac_100k.sbatch

set -euo pipefail

PROJECT_DIR="/hpcdisk1/limk_group/caiqy/project/260428_greengene2"
mkdir -p "$PROJECT_DIR/logs"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate baseBio

N="${SLURM_CPUS_PER_TASK:-64}"
export OMP_NUM_THREADS="$N"
export WRITE_BATCH="${WRITE_BATCH:-500}"
export PYTHONUNBUFFERED=1

echo "Python: $(python --version)"
echo "OMP_NUM_THREADS: $OMP_NUM_THREADS"
echo "WRITE_BATCH: $WRITE_BATCH"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-N/A}"
echo "Host: $(hostname)"
python -c "import unifrac; print('unifrac:', getattr(unifrac, '__version__', 'unknown'))"

cd "$PROJECT_DIR"
python -u scripts/SampleDistance100k/04_compute_wunifrac_100k.py
