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
# # 05: 50k × 50k Weighted Normalized UniFrac 距离矩阵
#
# 输入：
# - `results/sample_distance/subset_50k.h5ad`        (50,000 × 8,114)
# - `results/sample_distance/genus_tree.nwk`         8,114-tip folded tree
#
# 输出：
# - 把 50k × 50k weighted UniFrac 写回到同一个 anndata 的 `obsp['distance_wunifrac']`（float16）
#
# ## 实现
#
# 用 **Striped Fast UniFrac**（McDonald 2018 / Sfiligoi 2022）via `unifrac` python 包
# `unifrac.weighted_normalized(table, tree)`：
# - 输入 1：BIOM-Format v2.1 文件（features=8,114 g__）
# - 输入 2：Newick 树文件（8,114 个 g__ tip）
# - 用 `OMP_NUM_THREADS` 控制并行
# - 返回 `skbio.DistanceMatrix`，shape (50000, 50000) float64
#
# **var_id → g__ 映射**：anndata var_names 是完整 6 级路径 (`d__...;g__Foo`)，
# 但树的 tip 是裸 g__ token。所以建 BIOM 时用 `vocab["Genus"]` (= g__ token) 作 obs_id。

# %%
from pathlib import Path
import os
import sys
import time
import tempfile
import numpy as np
import pandas as pd
import anndata as ad

ad.settings.allow_write_nullable_strings = True

ROOT = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
ANN_PATH  = ROOT / "results/sample_distance/subset_50k.h5ad"
TREE_PATH = ROOT / "results/sample_distance/genus_tree.nwk"
VOCAB_IN  = ROOT / "results/phylogeny/genus_vocab.tsv"
TMP_DIR   = ROOT / "results/sample_distance" / "_tmp_unifrac"
TMP_DIR.mkdir(exist_ok=True)
BIOM_TMP  = TMP_DIR / "subset_50k.biom"

# 并行：用尽可能多的 CPU
N_THREADS = os.environ.get("OMP_NUM_THREADS", str(os.cpu_count()))
os.environ["OMP_NUM_THREADS"] = N_THREADS
print(f"OMP_NUM_THREADS = {N_THREADS}")

# %% [markdown]
# ## §1 读 anndata + vocab，建 var_id → g__ 映射

# %%
print(f"读 {ANN_PATH.name} ...")
adata = ad.read_h5ad(ANN_PATH)
print(f"  shape: {adata.shape}")

vocab = pd.read_csv(VOCAB_IN, sep="\t", index_col="var_id")
# vocab.index 是完整 var_id，vocab["Genus"] 是 g__ token
# anndata var_names 应该就是 vocab.index 的子集（理论上完全一致）
var_ids = adata.var_names.astype(str).values
miss_in_vocab = set(var_ids) - set(vocab.index.astype(str))
assert len(miss_in_vocab) == 0, f"anndata 有 {len(miss_in_vocab)} 个 var 不在 vocab"

# 取每个 var_id 对应的 g__ token，作为 BIOM 表的 observation id
genus_for_var = vocab.loc[var_ids, "Genus"].astype(str).values
assert len(set(genus_for_var)) == len(genus_for_var), "var → genus 出现重复"
print(f"var_id → g__ 映射建立 ({len(genus_for_var)} 个一一对应)")

# %% [markdown]
# ## §2 转 relative abundance + 写 BIOM 临时文件

# %%
import scipy.sparse as sp
import biom

print("\n转 relative abundance ...")
t0 = time.time()
X = adata.X.astype(np.float64)
row_sum = np.asarray(X.sum(axis=1)).ravel()
assert (row_sum > 0).all(), "存在零行"
inv = sp.diags(1.0 / row_sum)
P = inv @ X  # CSR sparse, float64, row-normalized
print(f"  P shape: {P.shape}  nnz: {P.nnz:,}  耗时 {time.time() - t0:.1f}s")

# %%
print("\n构造 BIOM Table 并写盘 ...")
t0 = time.time()
# biom.Table 接受 (features, samples) 形状的数据，所以要转置
data_T = P.T.tocsr()
sample_ids = adata.obs_names.astype(str).tolist()
table = biom.Table(
    data=data_T,
    observation_ids=genus_for_var.tolist(),
    sample_ids=sample_ids,
)
# 写 HDF5 BIOM v2.1
import h5py
with h5py.File(BIOM_TMP, "w") as h:
    table.to_hdf5(h, generated_by="SampleDistance/05_compute_wunifrac")
print(f"  BIOM 写盘耗时 {time.time() - t0:.1f}s  大小 {BIOM_TMP.stat().st_size/1024**2:.0f} MB")

# %% [markdown]
# ## §3 调 unifrac.weighted_normalized

# %%
import unifrac

print(f"\n开始计算 weighted normalized UniFrac (threads={N_THREADS}) ...")
t0 = time.time()
dm = unifrac.weighted_normalized(
    table=str(BIOM_TMP),
    phylogeny=str(TREE_PATH),
    threads=int(N_THREADS),
    variance_adjusted=False,
    bypass_tips=False,
)
elapsed = time.time() - t0
print(f"  耗时 {elapsed/60:.1f} min")
print(f"  DistanceMatrix shape: {dm.shape}")
print(f"  ids[:3]: {dm.ids[:3]}")

# %% [markdown]
# ## §4 把 distance 矩阵按 anndata obs_names 顺序对齐

# %%
print("\n对齐 obs_names 顺序 ...")
D = dm.data
ids = list(dm.ids)
pos = {s: i for i, s in enumerate(ids)}
order = np.array([pos[s] for s in sample_ids])
assert len(order) == adata.n_obs
D_aligned = D[np.ix_(order, order)].astype(np.float32)

print(f"  D_aligned shape: {D_aligned.shape}  dtype: {D_aligned.dtype}")
print(f"  对角 max: {np.diag(D_aligned).max():.6e}")
print(f"  对称偏差: {np.abs(D_aligned - D_aligned.T).max():.6e}")
print(f"  非对角 min: {D_aligned[np.triu_indices_from(D_aligned, k=1)].min():.6f}")
print(f"  非对角 max: {D_aligned[np.triu_indices_from(D_aligned, k=1)].max():.6f}")
print(f"  非对角 median: {np.median(D_aligned[np.triu_indices_from(D_aligned, k=1)]):.4f}")

# 强制对称 + 对角清零（消除 fp 舍入）
D_aligned = (D_aligned + D_aligned.T) / 2
np.fill_diagonal(D_aligned, 0.0)

# %% [markdown]
# ## §5 写回 obsp（float16）

# %%
print(f"\n转 float16 并写回 obsp['distance_wunifrac'] ...")
D16 = D_aligned.astype(np.float16)
print(f"  float16 mem: {D16.nbytes/1024**3:.2f} GB")
adata.obsp["distance_wunifrac"] = D16

# %%
print(f"\n写回 {ANN_PATH.name} (compression=gzip) ...")
adata.write_h5ad(ANN_PATH, compression="gzip")
print(f"  新文件大小: {ANN_PATH.stat().st_size/1024**3:.2f} GB")

# %% [markdown]
# ## §6 清理临时文件

# %%
BIOM_TMP.unlink()
TMP_DIR.rmdir()
print(f"已清理 {BIOM_TMP}")

# %% [markdown]
# ## §7 读回验证

# %%
b = ad.read_h5ad(ANN_PATH, backed="r")
print(f"\nobsp keys: {list(b.obsp.keys())}")
arr = np.asarray(b.obsp["distance_wunifrac"][:5, :5])
print(f"distance_wunifrac[:5,:5] (dtype={arr.dtype}):")
print(arr)
