#!/usr/bin/env bash
# 从 otus.97.allinfo 提取 97% OTU 代表全长 16S 序列
# 序列 ID 使用 97_XXXXX 格式，与 BIOM 矩阵 OTU ID 对齐
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
ALLINFO="$PROJECT_DIR/rawdata/MicrobeAtlas/OTU_count/otus.97.allinfo"
OUT="$PROJECT_DIR/results/rep_seqs/otus97_rep.fasta"

mkdir -p "$PROJECT_DIR/results/rep_seqs"

echo "[1/2] 提取序列..."
awk -F'\t' '{
    n = split($1, parts, ";")
    if ($7 != "") printf ">%s\n%s\n", parts[n], $7
}' "$ALLINFO" > "$OUT"

echo "[2/2] 验证输出..."
SEQ_COUNT=$(grep -c "^>" "$OUT")
echo "  提取序列数: $SEQ_COUNT  (期望: 111870)"

echo "  --- ID 样例 ---"
grep "^>" "$OUT" | head -3

echo "  --- 序列长度样例 (bp) ---"
awk '/^>/{if(seq) print length(seq); seq=""} !/^>/{seq=seq$0} END{if(seq) print length(seq)}' \
  "$OUT" | head -5 || true

echo "完成: $OUT"
