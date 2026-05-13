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
# # ResMicroDb 07: attach metadata 到 anndata.obs
#
# 把 06b 修正后的 `metadata_all.standardized.fixed.parquet` (135,746 × 36) 按 `Run` 左 join 进
# `resmicrodb.gg2.genus.qc.h5ad` (93,425 × 4,952, obs=0 列)，得到
# `resmicrodb.gg2.genus.qc.with_meta.h5ad`。
#
# **设计要点**：
# - parquet 保留了 dtype（category/nullable string/boolean/float），直接 read 不需重转
# - anndata.obs_names 就是 Run，metadata 的 Run 列也是 100% 覆盖（前一步实测：
#   交集 = anndata 全集 93,425；无样本缺 metadata）
# - 列顺序：`Database` (常量 'ResMicroDb') 放第一列，之后保留 metadata 自身 36 列顺序
# - 不做任何 filter，纯粹 attach；filter 走 08 步

# %%
from pathlib import Path
import pandas as pd
import anndata as ad

# anndata < 0.11 默认禁写 pd.StringArray (string[python])，显式打开
ad.settings.allow_write_nullable_strings = True

ROOT = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
IN_H5AD = ROOT / "results/feature_table/resmicrodb.gg2.genus.qc.h5ad"
IN_META = ROOT / "results/feature_table/metadata_all.standardized.fixed.parquet"
OUT_H5AD = ROOT / "results/feature_table/resmicrodb.gg2.genus.qc.with_meta.h5ad"

# %% [markdown]
# ## §1 读入 anndata + metadata

# %%
adata = ad.read_h5ad(IN_H5AD)
print(f"anndata: {adata.shape}  obs cols: {len(adata.obs.columns)}")

md = pd.read_parquet(IN_META)
print(f"metadata: {md.shape}")
print(f"metadata dtypes 摘要:")
print(md.dtypes.value_counts())

# %% [markdown]
# ## §2 按 Run 左 join

# %%
runs = adata.obs_names.to_series()
md_idx = md.set_index("Run", drop=False)  # drop=False 让 Run 既作 index 又作 obs 列

missing = runs[~runs.isin(md_idx.index)]
assert len(missing) == 0, f"{len(missing)} anndata 样本在 metadata 里找不到"

obs_new = md_idx.loc[runs.values].copy()
obs_new.index = adata.obs_names
obs_new.insert(0, "Database", pd.Categorical(["ResMicroDb"] * len(obs_new)))

print(f"obs_new: {obs_new.shape}")
print(f"列顺序: {list(obs_new.columns)}")

# %% [markdown]
# ## §3 校验：dtype 与覆盖率

# %%
print("=== dtype 保留检查 ===")
for c in obs_new.columns:
    print(f"  {c:<30} {str(obs_new[c].dtype)}")

print()
print("=== 非空率（top 10 缺失最严重的列） ===")
nonnull_pct = (obs_new.notna().sum() / len(obs_new) * 100).sort_values()
print(nonnull_pct.head(10).round(1).to_string())

# %% [markdown]
# ## §4 写出 anndata（obs 37 列：Database + 36 metadata；Run 既作 obs_names 又作列）

# %%
adata.obs = obs_new
print(f"final adata: {adata.shape}  obs cols: {len(adata.obs.columns)}")

adata.write_h5ad(OUT_H5AD, compression="gzip")
print(f"written: {OUT_H5AD}")
print(f"size: {OUT_H5AD.stat().st_size / 1024**2:.1f} MB")

# %% [markdown]
# ## §5 回读 sanity check

# %%
b = ad.read_h5ad(OUT_H5AD, backed="r")
print(f"回读: {b.shape}  obs cols: {len(b.obs.columns)}")
assert b.shape == adata.shape
assert list(b.obs.columns) == list(obs_new.columns)
print("OK")
