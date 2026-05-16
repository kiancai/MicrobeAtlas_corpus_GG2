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
# # 02: 按 index 切 50k 子集 anndata
#
# 输入：
# - `results/feature_table/merged.gg2.with_phylo.h5ad`  (1,826,126 × 8,114)
# - `results/sample_distance/subset_50k_index.tsv`     50,000 行（01 步产物）
#
# 输出：
# - `results/sample_distance/subset_50k.h5ad`  (50,000 × 8,114)：
#   - `X`：从主表 copy（CSR int32 稀疏）
#   - `obs`：54 列 + 新增 `stratum_id` + `sub_stratum` 共 56 列
#   - `var`：8,114 行（保留 6 级 taxonomy + `observed` bool）
#   - `varp`：`taxo_dist` + `phylo_dist`（从主表搬过来，**不变**）
#
# **不做距离计算**，距离矩阵由 04/05 步追加到 `obsp`。
# **不写 obsm**，PCoA 坐标由 06 步追加。

# %%
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad

ad.settings.allow_write_nullable_strings = True

ROOT = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
ANN_IN  = ROOT / "results/feature_table/merged.gg2.with_phylo.h5ad"
IDX_IN  = ROOT / "results/sample_distance/subset_50k_index.tsv"
ANN_OUT = ROOT / "results/sample_distance/subset_50k.h5ad"

# %% [markdown]
# ## §1 读 index + 主表（backed）

# %%
idx = pd.read_csv(IDX_IN, sep="\t", dtype={"obs_name": str})
print(f"index: {len(idx):,} 行")
print(idx.head())

print(f"\n读 {ANN_IN.name} (backed) ...")
adata = ad.read_h5ad(ANN_IN, backed="r")
print(f"  shape: {adata.shape}")
print(f"  varp: {list(adata.varp.keys())}")

# %% [markdown]
# ## §2 找子集行号并按 index 顺序切片

# %%
all_names = adata.obs_names.astype(str).values
name_to_pos = pd.Series(np.arange(len(all_names)), index=all_names)

# 校验 index 里所有 obs_name 都在主表里
missing = idx["obs_name"].values[~idx["obs_name"].isin(name_to_pos.index)]
assert len(missing) == 0, f"index 里 {len(missing)} 个 obs_name 不在主表中: {missing[:5]}"

row_pos = name_to_pos.loc[idx["obs_name"].values].values.astype(np.int64)
print(f"行号范围: min={row_pos.min()}  max={row_pos.max()}")
assert len(np.unique(row_pos)) == len(row_pos), "row_pos 有重复"

# %% [markdown]
# ## §3 切 X / obs / var / varp

# %%
print("切 X ...")
# backed 模式下 adata.X 是 SparseDataset，需要先按行索引取出
X_sub = adata.X[row_pos, :]
# 转 CSR（如果原本不是）
import scipy.sparse as sp
if not sp.isspmatrix_csr(X_sub):
    X_sub = sp.csr_matrix(X_sub)
print(f"  X_sub shape: {X_sub.shape}  nnz: {X_sub.nnz:,}  dtype: {X_sub.dtype}")

# %%
print("切 obs ...")
obs_sub = adata.obs.iloc[row_pos].copy()
obs_sub = obs_sub.reset_index(drop=True)
obs_sub.index = idx["obs_name"].values  # 用 obs_name 作 index
obs_sub.index.name = None
# 追加抽样信息
obs_sub["stratum_id"]  = pd.Categorical(idx["stratum_id"].values)
obs_sub["sub_stratum"] = pd.array(idx["sub_stratum"].astype(str).values, dtype="string")
print(f"  obs_sub shape: {obs_sub.shape}")
print(f"  新列 stratum_id 取值数: {obs_sub['stratum_id'].nunique()}")

# %%
print("切 var ...")
var_sub = adata.var.copy()
print(f"  var_sub shape: {var_sub.shape}")

# %%
print("拷 varp ...")
# varp 是 var × var 矩阵，不随 obs 切片变化，原样保留
varp_dict = {}
for k in adata.varp.keys():
    arr = np.asarray(adata.varp[k])
    varp_dict[k] = arr
    print(f"  varp['{k}']: shape={arr.shape}  dtype={arr.dtype}  mem={arr.nbytes/1024**2:.1f} MB")

# %% [markdown]
# ## §4 构造子集 anndata 并写盘

# %%
sub = ad.AnnData(
    X=X_sub,
    obs=obs_sub,
    var=var_sub,
)
for k, v in varp_dict.items():
    sub.varp[k] = v

print(f"\nsub anndata: {sub.shape}")
print(f"  obs 列: {len(sub.obs.columns)}")
print(f"  var 列: {len(sub.var.columns)}")
print(f"  varp:   {list(sub.varp.keys())}")

# %%
print(f"\n写出 {ANN_OUT.name} (compression=gzip) ...")
sub.write_h5ad(ANN_OUT, compression="gzip")
print(f"  大小: {ANN_OUT.stat().st_size / 1024**2:.1f} MB")

# %% [markdown]
# ## §5 读回 sanity

# %%
b = ad.read_h5ad(ANN_OUT, backed="r")
print(f"\n读回 shape: {b.shape}")
print(f"  obs[:3] index: {list(b.obs_names[:3])}")
print(f"  stratum_id 唯一数: {b.obs['stratum_id'].nunique()}")
print(f"  varp: {list(b.varp.keys())}")
# 距离矩阵随机抽对
print(f"  varp['phylo_dist'][0,1]: {float(np.asarray(b.varp['phylo_dist'][0,1])):.4f}")
print(f"  varp['taxo_dist'][0,1]:  {int(np.asarray(b.varp['taxo_dist'][0,1]))}")
