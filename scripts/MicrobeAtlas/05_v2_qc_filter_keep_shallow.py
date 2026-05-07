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
# # QC 过滤（v2：保留浅层注释）
#
# 与 `05_qc_filter.py` 唯一区别：**不删未到 genus 的 var**。其它逻辑（迭代、阈值、校验）完全一致。
#
# **过滤步骤**
# 1. 删 prevalence == 0 的 var（在当前 sample 集上）
# 2. 删 `row_sum < MIN_READS` 或 `row_nnz < MIN_FEATURES` 的 sample
# 3. 重复 1-2 直到 shape 不再变化
#
# **输入**：`results/feature_table/gg2.{full,minfilter}.nonzero.h5ad`
# **输出**：`results/feature_table/gg2.{full,minfilter}.qc_v2.h5ad`

# %%
from pathlib import Path

import numpy as np
import anndata as ad

PROJECT_DIR = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
OUT_DIR     = PROJECT_DIR / "results/feature_table"

MIN_READS    = 1000
MIN_FEATURES = 5

INPUTS = [
    ("full",      OUT_DIR / "gg2.full.nonzero.h5ad",      OUT_DIR / "gg2.full.qc_v2.h5ad"),
    ("minfilter", OUT_DIR / "gg2.minfilter.nonzero.h5ad", OUT_DIR / "gg2.minfilter.qc_v2.h5ad"),
]

print(f"阈值: MIN_READS = {MIN_READS}, MIN_FEATURES = {MIN_FEATURES}, 删除 shallow = False")


# %%
def compute_row_stats(X):
    row_nnz = np.diff(X.indptr)
    cs = np.concatenate(([0], np.cumsum(X.data, dtype=np.int64)))
    row_sum = np.diff(cs[X.indptr])
    return row_sum, row_nnz


def compute_var_prevalence(X, n_var):
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

    # ---- 迭代直到收敛 ----
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

        print(f"  iter {n_iter}: 删 {n_drop_var:>4} var + {n_drop_sample:>10,} sample → {adata.shape}")

    print(f"\n  收敛：{n_iter} 轮，最终形状 {adata.shape}")

    # ---- 终态校验 ----
    prev = compute_var_prevalence(adata.X, adata.n_vars)
    row_sum, row_nnz = compute_row_stats(adata.X)
    assert (prev > 0).all(),                "存在 prevalence==0 的 var，迭代未收敛"
    assert (row_sum >= MIN_READS).all(),    "存在 row_sum < MIN_READS 的 sample"
    assert (row_nnz >= MIN_FEATURES).all(), "存在 row_nnz < MIN_FEATURES 的 sample"
    print(f"  ✓ 校验通过：所有阈值都严格满足")

    print(f"  写出 {out_path.name} ...")
    adata.write_h5ad(out_path, compression='gzip')
    del adata

print("\n完成。")
