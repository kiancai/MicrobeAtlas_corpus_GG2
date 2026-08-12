#!/bin/bash
#SBATCH -J mbm_orient
#SBATCH -p corexd192
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem 128G
#SBATCH -t 04:00:00
#SBATCH --output=/hpcdisk1/limk_group/caiqy/project/260428_greengene2/results/orientation_check_mbm/logs/orient_%A_%a.log
#SBATCH --error=/hpcdisk1/limk_group/caiqy/project/260428_greengene2/results/orientation_check_mbm/logs/orient_%A_%a.log
# MBM auto vs both 方向对照（每个 array 任务对 1 个输入跑两遍：auto 和 both）
# 提交: sbatch -c 8 --mem=128G -t 04:00:00 --array=0-15 scripts/databases/MicroBioMap/orientation_check.sbatch
set -euo pipefail

PROJECT_DIR="/hpcdisk1/limk_group/caiqy/project/260428_greengene2"
CLASSIFIER="/hpcdisk1/limk_group/caiqy/tools/qiime2_database/greengenes2/2024.09.backbone.full-length.nb.sklearn-1.4.2.qza"
OUT="$PROJECT_DIR/results/orientation_check_mbm"
MANIFEST="$OUT/manifest.tsv"

TID="${SLURM_ARRAY_TASK_ID:?需要 --array}"
LINE=$(awk -F'\t' -v t="$TID" '$1==t{print; exit}' "$MANIFEST")
[[ -n "$LINE" ]] || { echo "ERROR: manifest 无 task_id=$TID"; exit 1; }
LABEL=$(echo "$LINE" | cut -f2)
FA=$(echo "$LINE" | cut -f3)
PREFIX=$(echo "$LINE" | cut -f4)

TMP="$OUT/tmp_${PREFIX}"; TAX="$OUT/tax"
mkdir -p "$TMP" "$TAX"
echo "=== job=$SLURM_JOB_ID task=$TID label=$LABEL host=$(hostname) $(date '+%F %T') ==="
[[ -f "$FA" ]] || { echo "ERROR: 缺 $FA"; exit 1; }
n_in=$(grep -c '^>' "$FA"); echo "输入序列数: $n_in"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate qiime2-2026.1

REPSEQS="$TMP/repseqs.qza"
rm -f "$REPSEQS"
qiime tools import --type 'FeatureData[Sequence]' --input-path "$FA" --output-path "$REPSEQS"

run_mode () {
    local MODE="$1"
    local out_txt="$TAX/${PREFIX}.${MODE}.tsv"
    if [[ -f "$out_txt" ]]; then
        local n_out=$(($(wc -l < "$out_txt") - 1))
        [[ "$n_out" -eq "$n_in" ]] && { echo "[$MODE] 已完成,跳过"; return; }
    fi
    local qza="$TMP/tax_${MODE}.qza"; local exp="$TMP/export_${MODE}"
    echo "=== classify ($MODE) $(date '+%T') ==="
    qiime feature-classifier classify-sklearn \
        --i-classifier "$CLASSIFIER" --i-reads "$REPSEQS" \
        --o-classification "$qza" --p-n-jobs "$SLURM_CPUS_PER_TASK" \
        --p-read-orientation "$MODE" --verbose
    rm -rf "$exp"; qiime tools export --input-path "$qza" --output-path "$exp"
    mv "$exp/taxonomy.tsv" "$out_txt"; rm -rf "$exp" "$qza"
    local n_out=$(($(wc -l < "$out_txt") - 1))
    echo "[$MODE] 完成: $out_txt ($n_out)"
    [[ "$n_out" -eq "$n_in" ]] || { echo "ERROR: $MODE 行数不一致"; exit 2; }
}

run_mode auto
run_mode both
rm -rf "$TMP"
echo "=== task $TID 完成 $(date '+%F %T') ==="
