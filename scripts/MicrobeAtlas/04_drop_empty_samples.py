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
# # 剔除"零 taxon"样本
#
# **目的**：03 步在源头过滤了 mitochondria / chloroplast / non-B-A，并剔除了
# BIOM 里的 `Unmapped` 桶。某些样本的 reads 全部落在被剔除的条目上，剩 0 个
# var——这些样本对下游任何分析都无意义，必须删掉。
#
# **本脚本只做这一件事**：丢 `nnz == 0` 的行。其它阈值（min reads / 低流行率
# 菌 / 低注释率样本等）请用 `05_qc_filter.py` 或后续脚本处理。
#
# **输入** (`results/feature_table/`)
# - `gg2.full.h5ad`
# - `gg2.minfilter.h5ad`
#
# **输出** (`results/feature_table/`)
# - `gg2.full.nonzero.h5ad`
# - `gg2.minfilter.nonzero.h5ad`
#
# 不覆盖原文件，便于回退。

# %%
from pathlib import Path

import numpy as np
import anndata as ad

PROJECT_DIR = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
OUT_DIR     = PROJECT_DIR / "results/feature_table"

INPUTS = [
    ("full",      OUT_DIR / "gg2.full.h5ad",      OUT_DIR / "gg2.full.nonzero.h5ad"),
    ("minfilter", OUT_DIR / "gg2.minfilter.h5ad", OUT_DIR / "gg2.minfilter.nonzero.h5ad"),
]

# %%
for tag, in_path, out_path in INPUTS:
    print(f"\n=== {tag} ===")
    print(f"  读取 {in_path.name}")
    adata = ad.read_h5ad(in_path)

    row_nnz  = np.diff(adata.X.indptr)
    keep     = row_nnz > 0
    n_total  = adata.n_obs
    n_keep   = int(keep.sum())
    n_drop   = n_total - n_keep

    print(f"  样本数: {n_total:,}")
    print(f"  零 taxon 丢弃: {n_drop:,}  ({n_drop/n_total*100:.2f}%)")
    print(f"  保留: {n_keep:,}")

    if n_drop == 0:
        print(f"  无可丢样本，直接复制即可（仍写出以保持文件命名一致）")

    adata_kept = adata[keep].copy()
    print(f"  写出 {out_path.name} ...")
    adata_kept.write_h5ad(out_path, compression='gzip')
    print(f"  形状: {adata_kept.shape}")
    del adata, adata_kept

print("\n完成。")
