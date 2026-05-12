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
# # MicrobeAtlas 07: metadata-based 样本筛选
#
# 输入 `gg2.full.qc.with_meta.h5ad` (1,762,635 × 6,306, obs 26 列)，按 metadata 规则筛
# → `gg2.full.qc.with_meta.filtered.h5ad`。
#
# **筛选规则**（语料库立场——只删"不可能进入任何训练分析"的样本，其他都保留）：
#
# | # | 规则 | 估计影响 |
# |---|---|---|
# | 1 | drop `Sequencing_Type ∈ {RNAseq, NaN}` | -29,934 (RNAseq 9,415 + NaN 20,519) |
# | — | 保留 `AMPLICON`（1,587,494）+ `WGS`（145,207），WGS 经 MicrobeAtlas 流水线提取 16S read 后与 AMPLICON 在 OTU 表语义上对齐 |
# | — | 阴性对照不删（保留全部部位，含可能的 control），下游训练时按 `MA_SampleSite` 等再筛 |
# | — | 异常 lat/lng 不处理 |
# | — | 宿主限定（人/动物/环境）不在此步做 |
#
# **预期产物**：约 1,732,701 × 6,306（var 不变）。
# 注意：本步可能产生新的全零 var（极小概率），但不在此步做 var QC——若需要可后续再走一遍 `04_drop_empty_samples` / QC。

# %%
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd

ad.settings.allow_write_nullable_strings = True

ROOT = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
IN_H5AD = ROOT / "results/feature_table/gg2.full.qc.with_meta.h5ad"
OUT_H5AD = ROOT / "results/feature_table/gg2.full.qc.with_meta.filtered.h5ad"

# %% [markdown]
# ## §1 读入

# %%
adata = ad.read_h5ad(IN_H5AD)
print(f"in: {adata.shape}  obs cols: {len(adata.obs.columns)}")

# %% [markdown]
# ## §2 当前 Sequencing_Type 分布

# %%
vc = adata.obs["Sequencing_Type"].value_counts(dropna=False)
print(vc.to_string())

# %% [markdown]
# ## §3 应用筛选

# %%
seq = adata.obs["Sequencing_Type"]
# 保留 AMPLICON, WGS；删 RNAseq, NaN
keep = seq.isin(["AMPLICON", "WGS"])
n_before = adata.n_obs
n_drop = (~keep).sum()
print(f"删除 {n_drop:,} 个样本（RNAseq + NaN）")
print(f"  其中 RNAseq:  {(seq == 'RNAseq').sum():,}")
print(f"  其中 NaN:     {seq.isna().sum():,}")

adata_f = adata[keep].copy()
print(f"after: {adata_f.shape}  obs cols: {len(adata_f.obs.columns)}")
print(f"保留率: {adata_f.n_obs / n_before * 100:.2f}%")

# %% [markdown]
# ## §4 保留后的 Sequencing_Type 分布

# %%
print(adata_f.obs["Sequencing_Type"].value_counts(dropna=False).to_string())

# %% [markdown]
# ## §5 写出 + 回读校验

# %%
adata_f.write_h5ad(OUT_H5AD, compression="gzip")
print(f"written: {OUT_H5AD}")
print(f"size: {OUT_H5AD.stat().st_size / 1024**2:.1f} MB")

b = ad.read_h5ad(OUT_H5AD, backed="r")
print(f"回读: {b.shape}  obs cols: {len(b.obs.columns)}")
assert b.shape == adata_f.shape
assert "RNAseq" not in set(b.obs["Sequencing_Type"].dropna().unique())
assert b.obs["Sequencing_Type"].isna().sum() == 0
print("OK")
