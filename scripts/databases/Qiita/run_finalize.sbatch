#!/bin/bash
#SBATCH -J qiita_finalize
#SBATCH -p corexd192
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem 256G
#SBATCH -t 12:00:00
#SBATCH --output=/hpcdisk1/limk_group/caiqy/project/260428_greengene2/logs/qiita_finalize_%j.log
#SBATCH --error=/hpcdisk1/limk_group/caiqy/project/260428_greengene2/logs/qiita_finalize_%j.log
# 分类 array 全部成功后跑: 04 合并 taxonomy → 05 per-study collapse → 06 拼表
# 提交: sbatch --dependency=afterok:<ARRAY_JID> -c 8 --mem=256G -t 12:00:00 scripts/Qiita/run_finalize.sbatch
set -euo pipefail
PROJECT_DIR="/hpcdisk1/limk_group/caiqy/project/260428_greengene2"
cd "$PROJECT_DIR"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate baseBio
PY=python
echo "=== $(date '+%F %T') 04 merge taxonomy ==="
$PY scripts/Qiita/04_merge_taxonomy.py
echo "=== $(date '+%F %T') 05 per-study genus collapse ==="
$PY scripts/Qiita/05_build_genus.py
echo "=== $(date '+%F %T') 06 concat → qiita.gg2.genus.h5ad ==="
$PY scripts/Qiita/06_concat_genus.py
echo "=== $(date '+%F %T') DONE ==="
