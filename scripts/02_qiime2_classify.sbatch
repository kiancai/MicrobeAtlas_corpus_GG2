#!/bin/bash
#SBATCH -J gg2_classify
#SBATCH -p corexd192
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH -t 96:00:00
#SBATCH --mem 128G
#SBATCH --output=/hpcdisk1/limk_group/caiqy/project/260428_greengene2/logs/gg2_classify_%j.log
#SBATCH --error=/hpcdisk1/limk_group/caiqy/project/260428_greengene2/logs/gg2_classify_%j.log
# 提交方式: sbatch scripts/02_qiime2_classify.sbatch

set -euo pipefail

PROJECT_DIR="/hpcdisk1/limk_group/caiqy/project/260428_greengene2"
FASTA="$PROJECT_DIR/results/rep_seqs/otus97_rep.fasta"
CLASSIFIER="/hpcdisk1/limk_group/caiqy/tools/qiime2_database/greengenes2/2024.09.backbone.full-length.nb.sklearn-1.4.2.qza"

mkdir -p "$PROJECT_DIR/results/taxonomy" "$PROJECT_DIR/logs"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate qiime2-2026.1

echo "=== 环境信息 ==="
echo "Python: $(python --version)"
python -c "import sklearn; print('sklearn:', sklearn.__version__)"
echo "QIIME2: $(qiime info 2>/dev/null | grep 'QIIME 2 release' || echo 'unknown')"
echo ""

echo "=== Step 1: 导入序列 ==="
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path "$FASTA" \
  --output-path "$PROJECT_DIR/results/rep_seqs/otus97_rep.qza"

echo "=== Step 2: 分类注释（使用 $SLURM_CPUS_PER_TASK 核）==="
qiime feature-classifier classify-sklearn \
  --i-classifier "$CLASSIFIER" \
  --i-reads "$PROJECT_DIR/results/rep_seqs/otus97_rep.qza" \
  --o-classification "$PROJECT_DIR/results/taxonomy/taxonomy.qza" \
  --p-n-jobs "$SLURM_CPUS_PER_TASK" \
  --verbose

echo "=== Step 3: 导出 TSV ==="
qiime tools export \
  --input-path "$PROJECT_DIR/results/taxonomy/taxonomy.qza" \
  --output-path "$PROJECT_DIR/results/taxonomy/"

LINES=$(wc -l < "$PROJECT_DIR/results/taxonomy/taxonomy.tsv")
echo "完成. taxonomy.tsv 行数: $LINES  (期望: 111871，含表头)"
head -3 "$PROJECT_DIR/results/taxonomy/taxonomy.tsv"
