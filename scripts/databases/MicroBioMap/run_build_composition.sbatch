#!/bin/bash
#SBATCH -J mbm_build
#SBATCH -p corexd192
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem 256G
#SBATCH -t 08:00:00
#SBATCH --output=/hpcdisk1/limk_group/caiqy/project/260428_greengene2/logs/mbm_build_%j.log
#SBATCH --error=/hpcdisk1/limk_group/caiqy/project/260428_greengene2/logs/mbm_build_%j.log
# 驱动 03→04→05，一键产出 GG2 重注释完整物种组成表
# 提交: sbatch -c 8 --mem=256G -t 08:00:00 scripts/databases/MicroBioMap/run_build_composition.sbatch
# 前置: 02 分类 array 全部完成（work/tax 下 31 片 taxonomy.tsv 齐）
set -euo pipefail

PROJECT_DIR="/hpcdisk1/limk_group/caiqy/project/260428_greengene2"
PY="$PROJECT_DIR/scripts/databases/MicroBioMap"
TAX_DIR="$PROJECT_DIR/rawdata/MicroBioMap/work/tax"

echo "=== host=$(hostname)  start=$(date '+%F %T') ==="

# 前置检查：31 片 taxonomy 齐
n_tax=$(ls "$TAX_DIR"/shard_*.taxonomy.tsv 2>/dev/null | wc -l)
echo "已完成 shard taxonomy: $n_tax/31"
[[ "$n_tax" -eq 31 ]] || { echo "ERROR: 分类未全部完成（$n_tax/31），先等 02 array 跑完"; exit 1; }

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate baseBio

FT="$PROJECT_DIR/results/feature_table"
GLOBAL_TAX="$PROJECT_DIR/rawdata/MicroBioMap/work/taxonomy_gg2.global.tsv"

echo "=== 03: 合并 shard taxonomy → 全局 taxonomy ==="
if [[ -f "$GLOBAL_TAX" ]]; then echo "已存在 $GLOBAL_TAX，跳过 03"; else python "$PY/03_merge_taxonomy.py"; fi

echo "=== 04: 项目表 + taxonomy → sample × ASV ==="
if [[ -f "$FT/microbiomap.gg2.asv.h5ad" ]]; then echo "已存在 asv.h5ad，跳过 04"; else python "$PY/04_build_asv_anndata.py"; fi

echo "=== 05: ASV → genus collapse（完整组成表，全样本）==="
if [[ -f "$FT/microbiomap.gg2.genus.h5ad" ]]; then echo "已存在 genus.h5ad，跳过 05"; else python "$PY/05_build_genus_anndata.py"; fi

echo "=== 完成: results/feature_table/microbiomap.gg2.genus.h5ad ==="
echo "end=$(date '+%F %T')"
