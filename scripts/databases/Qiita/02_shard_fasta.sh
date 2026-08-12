#!/bin/bash
# 全局 ASV fasta 切片(供 GG2 NB array 分类)。照搬 MBM 01_shard_fasta.sh, 改路径。
# 用法: bash scripts/Qiita/02_shard_fasta.sh [SHARD_SEQS]
set -euo pipefail
PROJECT_DIR="/hpcdisk1/limk_group/caiqy/project/260428_greengene2"
WORK="$PROJECT_DIR/rawdata/Qiita/work"
FA="$WORK/asv_global.fa"
SHARD_DIR="$WORK/shards"
SHARD_SEQS="${1:-120000}"
[[ -f "$FA" ]] || { echo "缺少 $FA"; exit 1; }
n_seq=$(grep -c '^>' "$FA")
echo "全局 ASV: $n_seq; 每片 $SHARD_SEQS"
rm -rf "$SHARD_DIR"; mkdir -p "$SHARD_DIR"
split -l "$((SHARD_SEQS * 2))" -d -a 2 --additional-suffix=.fa "$FA" "$SHARD_DIR/shard_"
n_shard=$(ls "$SHARD_DIR"/shard_*.fa | wc -l)
echo "生成 $n_shard 片"
tot=$(cat "$SHARD_DIR"/shard_*.fa | grep -c '^>')
[[ "$tot" -eq "$n_seq" ]] || { echo "ERROR: 分片总数 $tot != $n_seq"; exit 2; }
echo "总数守恒 ✓  提交: sbatch -c 8 --mem=128G -t 04:00:00 --array=0-$((n_shard-1)) scripts/Qiita/03_classify_shard.sbatch"
