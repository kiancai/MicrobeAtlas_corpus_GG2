#!/bin/bash
#SBATCH -J mbm_gg2
#SBATCH -p corexd192
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem 128G
#SBATCH -t 04:00:00
#SBATCH --output=/hpcdisk1/limk_group/caiqy/project/260428_greengene2/logs/mbm_gg2_%A_%a.log
#SBATCH --error=/hpcdisk1/limk_group/caiqy/project/260428_greengene2/logs/mbm_gg2_%A_%a.log
# 提交: sbatch -c 8 --mem=128G -t 04:00:00 --array=0-30 scripts/databases/MicroBioMap/02_classify_shard.sbatch
# 复用 ResMicroDb 的 GG2 24.09 NB(both) 配方，分片版（全局 ASV）
set -euo pipefail

PROJECT_DIR="/hpcdisk1/limk_group/caiqy/project/260428_greengene2"
WORK="$PROJECT_DIR/rawdata/MicroBioMap/work"
CLASSIFIER="/hpcdisk1/limk_group/caiqy/tools/qiime2_database/greengenes2/2024.09.backbone.full-length.nb.sklearn-1.4.2.qza"

SHARD=$(printf "%02d" "${SLURM_ARRAY_TASK_ID:?需要 --array}")
FA="$WORK/shards/shard_${SHARD}.fa"
TAX_DIR="$WORK/tax"
TMP_DIR="$WORK/tmp_import"
OUT_QZA="$TAX_DIR/shard_${SHARD}.taxonomy.qza"
OUT_TXT="$TAX_DIR/shard_${SHARD}.taxonomy.tsv"
TMP_REPSEQS="$TMP_DIR/shard_${SHARD}.repseqs.qza"
EXPORT_DIR="$TMP_DIR/export_${SHARD}"

mkdir -p "$TAX_DIR" "$TMP_DIR"
echo "=== job=$SLURM_JOB_ID array=$SLURM_ARRAY_TASK_ID shard=$SHARD host=$(hostname) ==="
echo "Started: $(date '+%F %T')"

[[ -f "$FA" ]] || { echo "ERROR: 缺少 $FA"; exit 1; }
n_in=$(grep -c '^>' "$FA")

# 断点续投：已完成且行数匹配则跳过
if [[ -f "$OUT_TXT" ]]; then
    n_out=$(($(wc -l < "$OUT_TXT") - 1))
    if [[ "$n_in" -eq "$n_out" ]]; then
        echo "shard $SHARD 已完成 ($n_out)，跳过"; exit 0
    fi
    echo "WARN: 已有输出但行数不匹配 ($n_in vs $n_out)，重做"
fi

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate qiime2-2026.1
echo "QIIME2: $(qiime info 2>/dev/null | grep -i 'QIIME 2 release' || echo unknown)"

echo "=== Step 1: import shard_${SHARD}.fa ($n_in seq) ==="
rm -f "$TMP_REPSEQS"
qiime tools import \
    --type 'FeatureData[Sequence]' \
    --input-path "$FA" \
    --output-path "$TMP_REPSEQS"

echo "=== Step 2: classify-sklearn (GG2 24.09 NB both, $SLURM_CPUS_PER_TASK jobs) ==="
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
n_out=$(($(wc -l < "$OUT_TXT") - 1))
echo "  in=$n_in  out=$n_out"
[[ "$n_in" -eq "$n_out" ]] || { echo "ERROR: 行数不一致"; exit 2; }
head -3 "$OUT_TXT"
echo "Finished: $(date '+%F %T')"
