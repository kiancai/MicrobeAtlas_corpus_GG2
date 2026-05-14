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
# # 10: anndata 扩到 GG2 24.09 全 genus 空间 + 挂距离矩阵
#
# 输入：
# - `results/feature_table/merged.gg2.h5ad`               (1,826,126 × 6,857)
# - `results/phylogeny/genus_vocab.tsv`                   GG2 24.09 全 ~8,114 genus
# - `results/phylogeny/genus_taxo_dist.npz`               8114² int8
# - `results/phylogeny/genus_phylo_dist.npz`              8114² float32
#
# 输出：
# - `results/feature_table/merged.gg2.with_phylo.h5ad`    (1,826,126 × 8,114)
#
# **核心动作**：把 anndata var 从 6,857 扩到 GG2 24.09 的全 8,114 个 genus，
# 补的 1,257 个 var 全部是 0 count（X 这部分是 CSR 全零，几乎不占空间）。
# 这样 var × var 的距离矩阵就可以原生挂 `varp`，全部信息留在 anndata 里。
#
# **var 排序约定**：
# - 原 6,857 个 var 位置不动（保留 09_merge 的原始顺序）
# - 新增 1,257 个排到末尾，按 var_id (taxonomy path) 字母序
# - 新增 `observed` bool 列：原 6,857 True，新增 1,257 False
#
# **为什么补回 1,257 不冲突 05 步的过滤**：
# - 05 步 feature 端只删了 `prevalence == 0` 和 `shallow`（无真实 g__ 标签）的 var
# - 没有"feature 至少在 5 个样本出现"这一条
# - 所以这 1,257 个补回来的 genus 真的在所有存活样本里都是 0
# - 补回去等于扩 var namespace 到 GG2 24.09 完整空间，无信号损失

# %%
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix, hstack

ad.settings.allow_write_nullable_strings = True

ROOT = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
ANN_IN   = ROOT / "results/feature_table/merged.gg2.h5ad"
VOCAB_IN = ROOT / "results/phylogeny/genus_vocab.tsv"
TAXO_IN  = ROOT / "results/phylogeny/genus_taxo_dist.npz"
PHYLO_IN = ROOT / "results/phylogeny/genus_phylo_dist.npz"
ANN_OUT  = ROOT / "results/feature_table/merged.gg2.with_phylo.h5ad"

RANK_COLS = ["Domain", "Phylum", "Class", "Order", "Family", "Genus"]

# %% [markdown]
# ## §1 读入 anndata、vocab 和两个距离矩阵

# %%
adata = ad.read_h5ad(ANN_IN)
print(f"anndata: {adata.shape}")
print(f"  X dtype: {adata.X.dtype}, type: {type(adata.X).__name__}")
print(f"  obs cols: {len(adata.obs.columns)}  var cols: {len(adata.var.columns)}")
print(f"  obsm: {list(adata.obsm.keys())}")
print(f"  obsp: {list(adata.obsp.keys())}")
print(f"  uns: {list(adata.uns.keys())}")
print(f"  layers: {list(adata.layers.keys())}")

# %%
vocab = pd.read_csv(VOCAB_IN, sep="\t", index_col="var_id")
print(f"vocab: {vocab.shape}")

taxo_npz  = np.load(TAXO_IN)
phylo_npz = np.load(PHYLO_IN)
taxo_mat  = taxo_npz["dist"]
phylo_mat = phylo_npz["dist"]
taxo_ids  = taxo_npz["var_id"].astype(str)
phylo_ids = phylo_npz["var_id"].astype(str)
print(f"taxo_dist  shape={taxo_mat.shape}  dtype={taxo_mat.dtype}")
print(f"phylo_dist shape={phylo_mat.shape}  dtype={phylo_mat.dtype}")

# %% [markdown]
# ### sanity：vocab.index 与两个 npz 的 var_id 完全一致

# %%
assert np.array_equal(vocab.index.values.astype(str), taxo_ids), "vocab vs taxo_dist 顺序不一致"
assert np.array_equal(vocab.index.values.astype(str), phylo_ids), "vocab vs phylo_dist 顺序不一致"
print("✅ vocab / taxo / phylo 的 var_id 顺序一致")

# %% [markdown]
# ## §2 找出 anndata 缺的 genus（应当 = 1,257 个）

# %%
original_var_ids = adata.var_names.values.astype(str)
original_set = set(original_var_ids)
all_set = set(vocab.index.astype(str))

# 缺失：vocab 有 anndata 没
missing_var_ids = sorted(all_set - original_set)
# 反向应该为空（anndata 都该在 vocab 里）
extra_in_adata = original_set - all_set

print(f"anndata var:        {len(original_set):,}")
print(f"vocab (GG2 全):     {len(all_set):,}")
print(f"anndata 缺的 genus: {len(missing_var_ids):,}")
print(f"anndata 多的 genus: {len(extra_in_adata):,}  (应为 0)")
assert len(extra_in_adata) == 0, f"⚠️  anndata 有 {len(extra_in_adata)} 个 GG2 没的 genus"

# %% [markdown]
# ## §3 构造新的 var 行（1,257 个 missing genus）

# %%
missing_var = vocab.loc[missing_var_ids, RANK_COLS].copy()
print(f"missing_var shape: {missing_var.shape}")
print(f"\n前 3 个 missing:")
print(missing_var.head(3))

# %% [markdown]
# ## §4 构造新的完整 var DataFrame
#
# - 原 6,857 行保持原顺序、加 `observed=True`
# - 新 1,257 行按字母序追加、`observed=False`
# - 总计 8,114 行，index 全局唯一

# %%
orig_var = adata.var[RANK_COLS].copy()
orig_var["observed"] = True
missing_var["observed"] = False

new_var = pd.concat([orig_var, missing_var], axis=0)
new_var.index.name = "var_id"

print(f"new_var: {new_var.shape}")
print(f"  observed=True:  {new_var['observed'].sum():,}")
print(f"  observed=False: {(~new_var['observed']).sum():,}")
assert new_var.index.is_unique, "var_id 不唯一"
assert len(new_var) == len(vocab), "总数 != vocab"

# %% [markdown]
# ## §5 扩展 X：右拼 1,257 个零列
#
# CSR 全零块基本不占空间。

# %%
n_obs = adata.n_obs
n_missing = len(missing_var_ids)
print(f"准备右拼 {n_missing} 个零列，dtype={adata.X.dtype}")

zero_block = csr_matrix((n_obs, n_missing), dtype=adata.X.dtype)
new_X = hstack([adata.X, zero_block], format="csr")
print(f"new_X shape: {new_X.shape}  dtype: {new_X.dtype}  nnz: {new_X.nnz:,}")
assert new_X.nnz == adata.X.nnz, "扩展后 nnz 变了，零列不该贡献 nnz"
assert new_X.shape == (n_obs, len(new_var))

# %% [markdown]
# ## §6 把两个距离矩阵按 new_var 顺序重排

# %%
# vocab/npz 的顺序是 var_id 字母序；new_var 的顺序是 [原 6857 在原顺序 + 1257 字母序]
target_order = list(new_var.index)
pos_in_npz = {v: i for i, v in enumerate(taxo_ids)}
order_idx = np.array([pos_in_npz[v] for v in target_order])

taxo_aligned  = taxo_mat[np.ix_(order_idx, order_idx)]
phylo_aligned = phylo_mat[np.ix_(order_idx, order_idx)]

print(f"taxo_aligned  shape={taxo_aligned.shape}  dtype={taxo_aligned.dtype}  "
      f"mem={taxo_aligned.nbytes / 1024**2:.1f} MB")
print(f"phylo_aligned shape={phylo_aligned.shape} dtype={phylo_aligned.dtype} "
      f"mem={phylo_aligned.nbytes / 1024**2:.1f} MB")

# sanity
assert taxo_aligned.shape  == (len(new_var), len(new_var))
assert phylo_aligned.shape == (len(new_var), len(new_var))
assert (np.diag(taxo_aligned)  == 0).all(), "taxo 对角线非 0"
assert (np.diag(phylo_aligned) == 0).all(), "phylo 对角线非 0"
print("✅ 对齐矩阵 sanity 通过")

# %% [markdown]
# ## §7 构造新的 AnnData，迁移 obs/obsm/obsp/uns/layers

# %%
new_adata = ad.AnnData(
    X=new_X,
    obs=adata.obs.copy(),
    var=new_var,
)

# 迁移其他槽位
for k in adata.obsm.keys():
    new_adata.obsm[k] = adata.obsm[k]
for k in adata.obsp.keys():
    new_adata.obsp[k] = adata.obsp[k]
for k in adata.uns.keys():
    new_adata.uns[k] = adata.uns[k]
# 如果有 layers，每个都要同样 hstack 零块
for k in adata.layers.keys():
    L = adata.layers[k]
    if hasattr(L, "tocsr"):  # 稀疏
        L_csr = L.tocsr() if not isinstance(L, csr_matrix) else L
        new_adata.layers[k] = hstack([L_csr, csr_matrix((n_obs, n_missing), dtype=L_csr.dtype)],
                                     format="csr")
    else:
        new_adata.layers[k] = np.hstack([L, np.zeros((n_obs, n_missing), dtype=L.dtype)])

# 挂上距离矩阵
new_adata.varp["taxo_dist"]  = taxo_aligned
new_adata.varp["phylo_dist"] = phylo_aligned

print(f"new_adata: {new_adata.shape}")
print(f"  obs: {new_adata.obs.shape}  var: {new_adata.var.shape}")
print(f"  obsm: {list(new_adata.obsm.keys())}")
print(f"  obsp: {list(new_adata.obsp.keys())}")
print(f"  varp: {list(new_adata.varp.keys())}")
print(f"  uns:  {list(new_adata.uns.keys())}")
print(f"  layers: {list(new_adata.layers.keys())}")

# %% [markdown]
# ## §8 写盘

# %%
print(f"写出 {ANN_OUT} ...")
new_adata.write_h5ad(ANN_OUT, compression="gzip")
size_mb = ANN_OUT.stat().st_size / 1024**2
orig_size = ANN_IN.stat().st_size / 1024**2
print(f"  原 {ANN_IN.name}:        {orig_size:.0f} MB")
print(f"  新 {ANN_OUT.name}:  {size_mb:.0f} MB")

# %% [markdown]
# ## §9 读回验证

# %%
b = ad.read_h5ad(ANN_OUT, backed="r")
print(f"shape: {b.shape}")
print(f"varp keys: {list(b.varp.keys())}")
print(f"var cols:  {list(b.var.columns)}")
print(f"observed=True:  {b.var['observed'].sum():,}")
print(f"observed=False: {(~b.var['observed']).sum():,}")
print(f"\ntaxo_dist[:3,:3]:\n{np.asarray(b.varp['taxo_dist'][:3, :3])}")
print(f"\nphylo_dist[:3,:3]:\n{np.asarray(b.varp['phylo_dist'][:3, :3])}")
