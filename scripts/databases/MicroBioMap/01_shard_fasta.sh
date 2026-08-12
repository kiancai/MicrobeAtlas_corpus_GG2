#!/bin/bash
# 把全局 ASV fasta 切成固定 seq 数的分片，供 GG2 NB array 分类
# 用法: bash scripts/databases/MicroBioMap/01_shard_fasta.sh [SHARD_SEQS]
set -euo pipefail

PROJECT_DIR="/hpcdisk1/limk_group/caiqy/project/260428_greengene2"
WORK="$PROJECT_DIR/rawdata/MicroBioMap/work"
FA="$WORK/asv_global.fa"
SHARD_DIR="$WORK/shards"
SHARD_SEQS="${1:-120000}"          # 每片序列数（默认 12 万）

[[ -f "$FA" ]] || { echo "缺少 $FA"; exit 1; }

n_seq=$(grep -c '^>' "$FA")
echo "全局 ASV fasta: $n_seq 序列；每片 $SHARD_SEQS 条"

rm -rf "$SHARD_DIR"; mkdir -p "$SHARD_DIR"
# fasta 每序列 2 行（header + 单行 seq），按 2*SHARD_SEQS 行切
split -l "$((SHARD_SEQS * 2))" -d -a 2 --additional-suffix=.fa "$FA" "$SHARD_DIR/shard_"

n_shard=$(ls "$SHARD_DIR"/shard_*.fa | wc -l)
echo "生成 $n_shard 片:"
for f in "$SHARD_DIR"/shard_*.fa; do
    printf "  %s  %s seq\n" "$(basename "$f")" "$(grep -c '^>' "$f")"
done
# 校验总数守恒
tot=$(cat "$SHARD_DIR"/shard_*.fa | grep -c '^>')
echo "分片总序列数: $tot (应 = $n_seq)"
[[ "$tot" -eq "$n_seq" ]] || { echo "ERROR: 分片总数不守恒"; exit 2; }
echo "ARRAY_MAX=$((n_shard - 1))  # 提交用 --array=0-$((n_shard - 1))"
