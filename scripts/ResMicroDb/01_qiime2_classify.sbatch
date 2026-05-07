#!/bin/bash
#SBATCH -J resmicrodb_gg2
#SBATCH -p corexd192
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem 128G
#SBATCH -t 04:00:00
#SBATCH --output=/hpcdisk1/limk_group/caiqy/project/260428_greengene2/logs/resmicrodb_gg2_%j.log
#SBATCH --error=/hpcdisk1/limk_group/caiqy/project/260428_greengene2/logs/resmicrodb_gg2_%j.log
# 提交方式: bash scripts/ResMicroDb/01_run_loop.sh
# 单独提交: sbatch scripts/ResMicroDb/01_qiime2_classify.sbatch <PROJECT>

set -euo pipefail

PROJECT="${1:?用法: sbatch 01_qiime2_classify.sbatch <PROJECT>}"

PROJECT_DIR="/hpcdisk1/limk_group/caiqy/project/260428_greengene2"
DATA_DIR="$PROJECT_DIR/rawdata/ResMicroDb/16S"
CLASSIFIER="/hpcdisk1/limk_group/caiqy/tools/qiime2_database/greengenes2/2024.09.backbone.full-length.nb.sklearn-1.4.2.qza"

PROJ_DIR="$DATA_DIR/$PROJECT"
RES_DIR="$PROJ_DIR/results"
LOG_DIR="$PROJ_DIR/log_gg2"
ASV_FA="$RES_DIR/asv.fa"
TMP_REPSEQS="$LOG_DIR/rep_seqs.gg2tmp.qza"
OUT_QZA="$RES_DIR/taxonomy_gg2.qza"
OUT_TXT="$RES_DIR/taxonomy_gg2.txt"
EXPORT_DIR="$LOG_DIR/export.gg2tmp"

mkdir -p "$LOG_DIR"

echo "=== job=$SLURM_JOB_ID  project=$PROJECT  host=$(hostname) ==="
echo "Started: $(date '+%F %T')"

[[ -f "$ASV_FA" ]] || { echo "ERROR: 缺少 $ASV_FA"; exit 1; }

# 断点续投：已完成则跳过
if [[ -f "$OUT_QZA" && -f "$OUT_TXT" ]]; then
    n_asv=$(grep -c "^>" "$ASV_FA")
    n_tax=$(($(wc -l < "$OUT_TXT") - 1))
    if [[ "$n_asv" -eq "$n_tax" ]]; then
        echo "已存在且行数匹配 (asv=$n_asv tax=$n_tax)，跳过"
        exit 0
    fi
    echo "WARN: 已有输出但行数不匹配 (asv=$n_asv tax=$n_tax)，重做"
fi

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate qiime2-2026.1

echo "QIIME2: $(qiime info 2>/dev/null | grep -i 'QIIME 2 release' || echo unknown)"

echo "=== Step 1: import asv.fa ==="
rm -f "$TMP_REPSEQS"
qiime tools import \
    --type 'FeatureData[Sequence]' \
    --input-path "$ASV_FA" \
    --output-path "$TMP_REPSEQS"

echo "=== Step 2: classify-sklearn (GG2 NB, $SLURM_CPUS_PER_TASK jobs) ==="
qiime feature-classifier classify-sklearn \
    --i-classifier "$CLASSIFIER" \
    --i-reads "$TMP_REPSEQS" \
    --o-classification "$OUT_QZA" \
    --p-n-jobs "$SLURM_CPUS_PER_TASK" \
    --p-read-orientation both \
    --verbose

echo "=== Step 3: export to TSV ==="
rm -rf "$EXPORT_DIR"
qiime tools export --input-path "$OUT_QZA" --output-path "$EXPORT_DIR"
mv "$EXPORT_DIR/taxonomy.tsv" "$OUT_TXT"
rm -rf "$EXPORT_DIR" "$TMP_REPSEQS"

echo "=== Step 4: verify ==="
n_asv=$(grep -c "^>" "$ASV_FA")
n_tax=$(($(wc -l < "$OUT_TXT") - 1))
echo "  ASV 数: $n_asv   taxonomy 行数(去表头): $n_tax"
[[ "$n_asv" -eq "$n_tax" ]] || { echo "ERROR: 行数不一致"; exit 2; }

head -3 "$OUT_TXT"

echo "Finished: $(date '+%F %T')"
