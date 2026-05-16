#!/bin/bash
#SBATCH -J mcf_final
#SBATCH -p corexd192
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH -t 02:00:00
#SBATCH --mem=64G
#SBATCH --output=/hpcdisk1/limk_group/caiqy/project/260428_greengene2/logs/mcf_final_%j.log
#SBATCH --error=/hpcdisk1/limk_group/caiqy/project/260428_greengene2/logs/mcf_final_%j.log

set -euo pipefail

ROOT=/hpcdisk1/limk_group/caiqy/project/260428_greengene2
cd "$ROOT"
mkdir -p logs

PYBIO=/hpcdisk1/limk_group/caiqy/miniforge3/envs/baseBio/bin/python

echo "=== 11_finalize_corpus  $(date) ==="
echo "node : $(hostname)"
echo "cpus : $SLURM_CPUS_PER_TASK"

"$PYBIO" scripts/Merge/11_finalize_corpus.py

echo "=== done  $(date) ==="
