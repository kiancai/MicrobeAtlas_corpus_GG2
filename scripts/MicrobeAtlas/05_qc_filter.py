# ---
# jupyter:
#   jupytext:
#     formats: py:percent,ipynb
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: baseBio
#     language: python
#     name: python3
# ---

# %% [markdown]
# # QC 过滤（v1：删浅层注释）
#
# **过滤步骤**
# 1. **结构性**：删所有"未注释到 genus"的 var（Genus 字段长度 ≤ 3，即 `''` 或 `'g__'` 占位符）
# 2. **数据驱动（迭代到不动点）**：
#    - 删 prevalence == 0 的 var
#    - 删 `row_sum < MIN_READS` 或 `row_nnz < MIN_FEATURES` 的 sample
#    - 重新统计 → 重新过滤 → 直到 shape 不再变化
#
# **为什么必须迭代**：删 var 改变 sample 的 row_sum / row_nnz；删 sample 又改变 var 的
# prevalence。一轮过完未必收敛——某个 var 在被删的几个样本里恰好是它仅有的出现，
# 第二轮就会变成 prevalence==0。不迭代会留下"漏网之鱼"。
#
# **输入**：`results/feature_table/gg2.{full,minfilter}.nonzero.h5ad`（来自 04）
# **输出**：`results/feature_table/gg2.{full,minfilter}.qc.h5ad`

# %%
from pathlib import Path

import numpy as np
import anndata as ad

PROJECT_DIR = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
OUT_DIR     = PROJECT_DIR / "results/feature_table"

MIN_READS    = 1000
MIN_FEATURES = 10

INPUTS = [
    ("full",      OUT_DIR / "gg2.full.nonzero.h5ad",      OUT_DIR / "gg2.full.qc.h5ad"),
    ("minfilter", OUT_DIR / "gg2.minfilter.nonzero.h5ad", OUT_DIR / "gg2.minfilter.qc.h5ad"),
]

print(f"阈值: MIN_READS = {MIN_READS}, MIN_FEATURES = {MIN_FEATURES}, 删除 shallow = True")


# %%
def compute_row_stats(X):
    """对 CSR 矩阵流式算 row_nnz 和 row_sum，避免一次性 dense 化。"""
    row_nnz = np.diff(X.indptr)
    cs = np.concatenate(([0], np.cumsum(X.data, dtype=np.int64)))
    row_sum = np.diff(cs[X.indptr])
    return row_sum, row_nnz


def compute_var_prevalence(X, n_var):
    """每个 var 在多少个样本中非零（CSR 下 indices 列号的频次）。"""
    return np.bincount(X.indices, minlength=n_var)


# %%
for tag, in_path, out_path in INPUTS:
    print(f"\n{'='*60}\n  {tag}\n{'='*60}")
    if not in_path.exists():
        print(f"  跳过：{in_path} 不存在，请先跑 04_drop_empty_samples.py")
        continue

    print(f"  读取 {in_path.name}")
    adata = ad.read_h5ad(in_path)
    print(f"  初始: {adata.shape}")

    # ---- Step 1: 结构性过滤——删 shallow var（未到 genus）----
    # 约定（来自 03 的 fillna('')）：Genus 列要么是 'g__<name>'（长度 > 3），
    # 要么是 'g__' 占位符（长度 == 3），要么是 ''（长度 == 0，原 taxon 没 g__ 段）
    is_real_genus = adata.var['Genus'].str.len() > 3
    n_shallow = int((~is_real_genus).sum())
    print(f"\n  Step 1: 删 shallow var (无真实 genus)")
    print(f"    删除 {n_shallow:,} / {adata.n_vars:,} 个 var")
    adata = adata[:, is_real_genus.values].copy()
    print(f"    →  {adata.shape}")

    # ---- Step 2: 迭代直到收敛 ----
    print(f"\n  Step 2: 迭代过滤（prevalence + min_reads + min_features）")
    prev_shape = None
    n_iter = 0
    while adata.shape != prev_shape:
        n_iter += 1
        prev_shape = adata.shape

        # 2a. 删 prevalence == 0 的 var（在当前 sample 集上）
        prev = compute_var_prevalence(adata.X, adata.n_vars)
        var_keep = prev > 0
        n_drop_var = int((~var_keep).sum())
        if n_drop_var > 0:
            adata = adata[:, var_keep].copy()

        # 2b. 删失败 sample（在当前 var 集上）
        row_sum, row_nnz = compute_row_stats(adata.X)
        fail_reads = row_sum < MIN_READS
        fail_feats = row_nnz < MIN_FEATURES
        sample_keep = ~(fail_reads | fail_feats)
        n_drop_sample = int((~sample_keep).sum())
        if n_drop_sample > 0:
            adata = adata[sample_keep].copy()

        print(f"    iter {n_iter}: 删 {n_drop_var:>4} var + {n_drop_sample:>10,} sample → {adata.shape}")

    print(f"\n  收敛：{n_iter} 轮")
    print(f"  最终形状: {adata.shape}")

    # ---- Step 3: 终态校验 ----
    prev = compute_var_prevalence(adata.X, adata.n_vars)
    row_sum, row_nnz = compute_row_stats(adata.X)
    assert (prev > 0).all(),                  "存在 prevalence==0 的 var，迭代未收敛"
    assert (row_sum >= MIN_READS).all(),      "存在 row_sum < MIN_READS 的 sample"
    assert (row_nnz >= MIN_FEATURES).all(),   "存在 row_nnz < MIN_FEATURES 的 sample"
    assert (adata.var['Genus'].str.len() > 3).all(), "存在 shallow var"
    print(f"  ✓ 校验通过：所有阈值都严格满足")

    # ---- 写出 ----
    print(f"  写出 {out_path.name} ...")
    adata.write_h5ad(out_path, compression='gzip')
    del adata

print("\n完成。")
