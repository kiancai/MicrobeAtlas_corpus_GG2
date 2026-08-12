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
# # ResMicroDb 05: QC 过滤（v1：删浅层注释）
#
# 镜像 MicrobeAtlas 05_qc_filter.py。差异：
# - 输入路径换成 ResMicroDb 04 的输出
# - 阈值保持一致；ResMicroDb 项目间深度差异较大，未来若需要一份「保留浅深度
#   样本」的版本，再仿造 MicrobeAtlas 05_v2 增加同名脚本即可
#
# **过滤步骤**
# 1. **结构性**：删所有"未注释到 genus"的 var（Genus 字段长度 ≤ 3）
# 2. **数据驱动（迭代到不动点）**：
#    - 删 prevalence == 0 的 var
#    - 删 row_sum < MIN_READS 或 row_nnz < MIN_FEATURES 的 sample
#    - 重复直到 shape 不变
#
# **输入**: `results/feature_table/resmicrodb.gg2.genus.nonzero.h5ad`
# **输出**: `results/feature_table/resmicrodb.gg2.genus.qc.h5ad`

# %%
from pathlib import Path

import numpy as np
import anndata as ad

PROJECT_DIR = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
IN_PATH     = PROJECT_DIR / "results/feature_table/resmicrodb.gg2.genus.nonzero.h5ad"
OUT_PATH    = PROJECT_DIR / "results/feature_table/resmicrodb.gg2.genus.qc.h5ad"

MIN_READS    = 1000
MIN_FEATURES = 5

print(f"阈值: MIN_READS = {MIN_READS}, MIN_FEATURES = {MIN_FEATURES}, 删除 shallow = True")


# %%
def compute_row_stats(X):
    """对 CSR 矩阵流式算 row_nnz 和 row_sum。"""
    row_nnz = np.diff(X.indptr)
    cs = np.concatenate(([0], np.cumsum(X.data, dtype=np.int64)))
    row_sum = np.diff(cs[X.indptr])
    return row_sum, row_nnz


def compute_var_prevalence(X, n_var):
    return np.bincount(X.indices, minlength=n_var)


# %%
print(f"\n读取 {IN_PATH.name}")
adata = ad.read_h5ad(IN_PATH)
print(f"初始: {adata.shape}")

# ---- Step 1: 删 shallow var ----
is_real_genus = adata.var['Genus'].str.len() > 3
n_shallow = int((~is_real_genus).sum())
print(f"\nStep 1: 删 shallow var (无真实 genus)")
print(f"  删除 {n_shallow:,} / {adata.n_vars:,} 个 var")
adata = adata[:, is_real_genus.values].copy()
print(f"  →  {adata.shape}")

# ---- Step 2: 迭代直到收敛 ----
print(f"\nStep 2: 迭代过滤")
prev_shape = None
n_iter = 0
while adata.shape != prev_shape:
    n_iter += 1
    prev_shape = adata.shape

    prev = compute_var_prevalence(adata.X, adata.n_vars)
    var_keep = prev > 0
    n_drop_var = int((~var_keep).sum())
    if n_drop_var > 0:
        adata = adata[:, var_keep].copy()

    row_sum, row_nnz = compute_row_stats(adata.X)
    fail_reads = row_sum < MIN_READS
    fail_feats = row_nnz < MIN_FEATURES
    sample_keep = ~(fail_reads | fail_feats)
    n_drop_sample = int((~sample_keep).sum())
    if n_drop_sample > 0:
        adata = adata[sample_keep].copy()

    print(f"  iter {n_iter}: 删 {n_drop_var:>4} var + {n_drop_sample:>8,} sample → {adata.shape}")

print(f"\n收敛: {n_iter} 轮，最终形状 {adata.shape}")

# ---- 校验 ----
prev = compute_var_prevalence(adata.X, adata.n_vars)
row_sum, row_nnz = compute_row_stats(adata.X)
assert (prev > 0).all(),                "存在 prevalence==0 的 var"
assert (row_sum >= MIN_READS).all(),    "存在 row_sum < MIN_READS"
assert (row_nnz >= MIN_FEATURES).all(), "存在 row_nnz < MIN_FEATURES"
assert (adata.var['Genus'].str.len() > 3).all(), "存在 shallow var"
print("✓ 校验通过")

print(f"\n写出 {OUT_PATH.name} ...")
adata.write_h5ad(OUT_PATH, compression='gzip')
print("完成。")
