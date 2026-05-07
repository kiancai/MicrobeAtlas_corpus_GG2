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
# # ResMicroDb 04: 剔除"零 taxon"样本
#
# **目的**：03 步过滤 mito/chloro/非 B-A 后，某些样本的 reads 全部落在被剔除
# 的 ASV 上，聚合后 nnz=0。这些样本对下游分析无意义，必须删掉。
#
# **本脚本只做这一件事**：丢 `nnz == 0` 的行；其它 QC 阈值在 05 处理。
#
# **输入**: `results/feature_table/resmicrodb.gg2.genus.h5ad`
# **输出**: `results/feature_table/resmicrodb.gg2.genus.nonzero.h5ad`
#
# 不覆盖原文件，便于回退。
#
# 与 MicrobeAtlas 04 的差异：
# - 同时报告"哪些项目最容易出零样本"（提示宿主 16S 污染严重的项目，便于
#   下游针对性 QC）

# %%
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad

PROJECT_DIR = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
IN_PATH     = PROJECT_DIR / "results/feature_table/resmicrodb.gg2.genus.h5ad"
OUT_PATH    = PROJECT_DIR / "results/feature_table/resmicrodb.gg2.genus.nonzero.h5ad"
ASV_PATH    = PROJECT_DIR / "results/feature_table/resmicrodb.gg2.asv.h5ad"  # 借 var.project 反查

# %%
print(f"读取 {IN_PATH.name}")
adata = ad.read_h5ad(IN_PATH)
print(adata)

row_nnz = np.diff(adata.X.indptr)
keep    = row_nnz > 0
n_total = adata.n_obs
n_keep  = int(keep.sum())
n_drop  = n_total - n_keep

print(f"\n样本数         : {n_total:,}")
print(f"零 taxon 丢弃   : {n_drop:,}  ({n_drop/n_total*100:.2f}%)")
print(f"保留           : {n_keep:,}")

# %% [markdown]
# ## 项目维度的零样本分布
#
# 反查每个零样本属于哪个项目（通过 02 输出的 sample × ASV AnnData，因为
# 03 输出的 obs 已经丢掉了 project 信息）。

# %%
if n_drop > 0:
    drop_samples = adata.obs.index[~keep]
    print(f"\n=== 零样本所属项目 (top 15) ===")
    asv_ad = ad.read_h5ad(ASV_PATH, backed='r')
    sample_to_proj = {}
    # 任取每个样本第一个非零 ASV 的 project；零样本不应在 ASV 层也是空，
    # 但若空就标 'UNKNOWN'
    asv_proj_arr = asv_ad.var['project'].values
    # backed 模式下 adata.X 是 _CSRDataset 代理，没有 indptr/indices；
    # 直接访问其底层 h5 group。indptr 整体很小，一次性载入。
    x_group = asv_ad.X.group
    indptr = x_group['indptr'][:]
    indices_ds = x_group['indices']
    sample_idx_map = {s: i for i, s in enumerate(asv_ad.obs.index)}
    for s in drop_samples:
        i = sample_idx_map.get(s)
        if i is None:
            sample_to_proj[s] = 'NOT_IN_ASV'
            continue
        row_start = int(indptr[i])
        row_end   = int(indptr[i + 1])
        if row_end == row_start:
            # ASV 层也空，从 sample id 前缀已不可得（项目前缀只在 var）
            # 退化方案：跨 obs 没法直接回查；标 NO_NONZERO
            sample_to_proj[s] = 'NO_NONZERO_ASV'
        else:
            asv_pos = int(indices_ds[row_start])
            sample_to_proj[s] = asv_proj_arr[asv_pos]
    asv_ad.file.close()

    proj_counts = pd.Series(sample_to_proj).value_counts()
    for proj, n in proj_counts.head(15).items():
        print(f"  {proj:<18}  {n:>5d}")
else:
    print("无零样本，跳过项目分布报告")

# %%
adata_kept = adata[keep].copy()
print(f"\n写出 {OUT_PATH.name} ...")
adata_kept.write_h5ad(OUT_PATH, compression='gzip')
print(f"形状: {adata_kept.shape}")

print("\n完成。")
