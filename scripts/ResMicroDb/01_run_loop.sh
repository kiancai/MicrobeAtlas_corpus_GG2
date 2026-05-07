#!/usr/bin/env bash
# 简单轮询提交器：维护用户队列在 CAP 上限，每 INTERVAL 秒补位一次。
#
# 用法:
#   bash scripts/ResMicroDb/01_run_loop.sh                # 默认 CAP=90, INTERVAL=30s
#   CAP=80 INTERVAL=60 bash scripts/ResMicroDb/01_run_loop.sh
#   bash scripts/ResMicroDb/01_run_loop.sh --dry-run      # 只打印待提交列表
#
# 后台运行:
#   nohup bash scripts/ResMicroDb/01_run_loop.sh > logs/run_loop.log 2>&1 &
#
# 中断: Ctrl-C 或 kill <PID>。已提交到 SLURM 的任务不受影响；下次再跑会跳过已完成项目。

set -uo pipefail

PROJECT_DIR="/hpcdisk1/limk_group/caiqy/project/260428_greengene2"
DATA_DIR="$PROJECT_DIR/rawdata/ResMicroDb/16S"
SCRIPTS_DIR="$PROJECT_DIR/scripts/ResMicroDb"
SBATCH_FILE="$SCRIPTS_DIR/01_qiime2_classify.sbatch"

CAP="${CAP:-90}"             # 用户队列上限（QOS 限制大致 100，留 buffer）
INTERVAL="${INTERVAL:-30}"   # 轮询间隔（秒）
JOB_NAME="resmicrodb_gg2"    # 与 sbatch -J 保持一致，用于精准过滤队列计数

# 集群本地 sbatch wrapper 要求命令行显式带资源 flag (即使 #SBATCH 里已写)
SBATCH_RES=( -c 8 --mem=128G -t 04:00:00 )

USER_NAME=$(whoami)
DRY_RUN=0
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=1

mkdir -p "$PROJECT_DIR/logs"

# ---------- 1. 收集 pending 项目 ----------
echo "=== 扫描 pending 项目 ==="
pending=()
for d in "$DATA_DIR"/*/; do
    proj=$(basename "$d")
    asv_fa="$d/results/asv.fa"
    [[ -f "$asv_fa" ]] || continue
    tax_qza="$d/results/taxonomy_gg2.qza"
    tax_txt="$d/results/taxonomy_gg2.txt"
    if [[ -f "$tax_qza" && -f "$tax_txt" ]]; then
        n_asv=$(grep -c "^>" "$asv_fa")
        n_tax=$(($(wc -l < "$tax_txt") - 1))
        [[ "$n_asv" -eq "$n_tax" ]] && continue
    fi
    pending+=("$proj")
done
TOTAL=${#pending[@]}
echo "  待处理: $TOTAL"

if [[ "$TOTAL" -eq 0 ]]; then
    echo "全部已完成，退出"
    exit 0
fi

if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "--dry-run: 前 5 个 = ${pending[@]:0:5}"
    exit 0
fi

# ---------- 2. 轮询提交 ----------
echo "=== 开始轮询 (CAP=$CAP, INTERVAL=${INTERVAL}s) ==="
echo "时间起点: $(date '+%F %T')"

submitted=0
idx=0
while [[ "$idx" -lt "$TOTAL" ]]; do
    # 当前队列里属于本任务的数量（已提交未完成）
    running=$(squeue -u "$USER_NAME" -h -o "%j" 2>/dev/null | grep -c "^$JOB_NAME$" || true)
    slots=$(( CAP - running ))

    if [[ "$slots" -le 0 ]]; then
        printf "[%s] 队列满 (%d/%d)，等 %ds\n" "$(date '+%H:%M:%S')" "$running" "$CAP" "$INTERVAL"
        sleep "$INTERVAL"
        continue
    fi

    batch_n=$(( TOTAL - idx ))
    [[ "$batch_n" -gt "$slots" ]] && batch_n=$slots

    printf "[%s] 队列 %d/%d，提交 %d 个 (剩 %d/%d)\n" \
        "$(date '+%H:%M:%S')" "$running" "$CAP" "$batch_n" "$((TOTAL-idx))" "$TOTAL"

    end_idx=$(( idx + batch_n ))
    while [[ "$idx" -lt "$end_idx" ]]; do
        proj="${pending[$idx]}"
        # 提交失败时不增 idx，下轮重试；连续失败也只是慢，不会卡死
        if jid=$(sbatch --parsable "${SBATCH_RES[@]}" "$SBATCH_FILE" "$proj" 2>&1); then
            if [[ "$jid" =~ ^[0-9]+$ ]]; then
                submitted=$(( submitted + 1 ))
                idx=$(( idx + 1 ))
            else
                echo "  WARN: 提交 $proj 返回非数字: $jid，1s 后重试"
                sleep 1
                break
            fi
        else
            echo "  WARN: 提交 $proj 失败: $jid，5s 后重试"
            sleep 5
            break
        fi
    done

    # 即使没填满也短暂歇一下，避免对调度器太凶
    [[ "$idx" -lt "$TOTAL" ]] && sleep 2
done

echo
echo "=== 全部提交完毕 ==="
echo "已提交: $submitted / $TOTAL"
echo "时间结束: $(date '+%F %T')"
echo
echo "查看队列: squeue -u $USER_NAME -n $JOB_NAME"
echo "查看日志: ls $PROJECT_DIR/logs/resmicrodb_gg2_*.log"
