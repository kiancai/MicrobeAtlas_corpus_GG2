#!/bin/bash
#SBATCH --job-name=mcf_meta
#SBATCH --output=scripts/merge_v2/logs/02_meta_%j.out
#SBATCH --error=scripts/merge_v2/logs/02_meta_%j.err
#SBATCH -p corexd192
#SBATCH -c 4
#SBATCH --mem=180G
#SBATCH -t 03:00:00

set -euo pipefail
cd /hpcdisk1/limk_group/caiqy/project/260428_greengene2
source /hpcdisk1/limk_group/caiqy/miniforge3/etc/profile.d/conda.sh
conda activate baseBio
echo "开始 $(date)  节点 $(hostname)"
python3 scripts/merge_v2/02_merge_metadata.py
echo "结束 $(date)"
