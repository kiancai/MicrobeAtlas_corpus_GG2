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
# # ResMicroDb 03: ASV → Genus AnnData (镜像 MicrobeAtlas 03)
#
# **输入** `results/feature_table/resmicrodb.gg2.asv.h5ad`  (02 输出)
# - sample × ASV，var 含 `gg2_Domain..gg2_Species` + `gg2_Confidence` + `silva_*` + `project`
#
# **输出** `results/feature_table/resmicrodb.gg2.genus.h5ad`
# - sample × genus_var
# - var 列与 MicrobeAtlas `gg2.full.h5ad.var` 一致: `[Domain, Phylum, Class, Order, Family, Genus]`
#   (字符串值，缺位用 `'d__'`/`'p__'` 等空占位符)
# - 不保留 SILVA 注释（聚合后 GG2 genus 与 SILVA genus 多对一，无意义）
# - 也不保留 project 列（一个 var 可能跨多项目）
#
# **过滤逻辑（与 MicrobeAtlas 03 完全一致）**
# 1. `gg2_Domain ∈ {d__Bacteria, d__Archaea}` —— 丢 Unassigned / Eukaryota
# 2. 丢 mitochondria / chloroplast (Taxon 含相应字符串)
#
# 浅层注释（任意深度未到 genus）的 ASV **不丢弃**，通过空占位符保留 counts。

# %%
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

PROJECT_DIR = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
IN_PATH     = PROJECT_DIR / "results/feature_table/resmicrodb.gg2.asv.h5ad"
OUT_PATH    = PROJECT_DIR / "results/feature_table/resmicrodb.gg2.genus.h5ad"

# %% [markdown]
# ## Step A: 读 02 输出 + 重建 6 级 var_id
#
# 与 MicrobeAtlas 03 保持一致的常量与函数。

# %%
RANK_COLS     = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
PREFIXES      = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
RANK_TO_GENUS = RANK_COLS[:6]
PFX_TO_GENUS  = PREFIXES[:6]


def deepest_rank(row, rank_cols=RANK_COLS, prefixes=PREFIXES):
    for col, pfx in zip(reversed(rank_cols), reversed(prefixes)):
        v = row[col]
        if v is not None and v != pfx and not (isinstance(v, float) and np.isnan(v)):
            return col
    return 'None'


def make_var_id(row, rank_cols=RANK_TO_GENUS, prefixes=PFX_TO_GENUS):
    parts = []
    for col, pfx in zip(rank_cols, prefixes):
        v = row[col]
        if v is None or (isinstance(v, float) and np.isnan(v)):
            parts.append(pfx)
        else:
            parts.append(v)
    return ';'.join(parts)


# %%
adata_in = ad.read_h5ad(IN_PATH)
print(adata_in)

# 把 gg2_* 列重命名为无前缀，方便复用 MicrobeAtlas 同款函数
gg2_cols = [f'gg2_{c}' for c in RANK_COLS]
tax = adata_in.var[gg2_cols + [c for c in adata_in.var.columns if c.startswith('silva_')] + ['project']].copy()
tax.columns = RANK_COLS + [c for c in tax.columns if c.startswith('silva_')] + ['project']
tax['Taxon_full'] = adata_in.var[gg2_cols].apply(
    lambda r: ';'.join('' if pd.isna(v) else v for v in r), axis=1
)

print(f"原始 ASV 数: {len(tax):,}")

# %% [markdown]
# ## Step B: 应用过滤 (Domain∈{B,A}、去 mito/chloro)

# %%
mask_dom    = tax['Domain'].isin(['d__Bacteria', 'd__Archaea'])
mask_mito   = ~tax['Taxon_full'].str.contains('mitochondri', case=False, na=False)
mask_chloro = ~tax['Taxon_full'].str.contains('chloroplast',  case=False, na=False)

print(f"过滤前 ASV 数: {len(tax):,}")
print(f"  丢 Domain 非 B/A : {(~mask_dom).sum():,}")
print(f"  丢 mitochondria  : {(~mask_mito).sum():,}")
print(f"  丢 chloroplast   : {(~mask_chloro).sum():,}")

keep_mask = (mask_dom & mask_mito & mask_chloro).values
tax_kept  = tax[keep_mask].copy()
print(f"过滤后 ASV 数: {len(tax_kept):,}")

# %%
tax_kept['deepest_rank'] = tax_kept.apply(deepest_rank, axis=1)
tax_kept['var_id']       = tax_kept.apply(make_var_id, axis=1)

RANK_ORDER = ['Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Domain', 'None']
print("\n=== ASV 按最深可用层级分布 ===")
depth_dist = tax_kept['deepest_rank'].value_counts().reindex(RANK_ORDER, fill_value=0)
for r in RANK_ORDER:
    n = depth_dist[r]
    print(f"  {r:<10} {n:>10,}  ({n/len(tax_kept)*100:5.2f}%)")

n_var = tax_kept['var_id'].nunique()
print(f"\nvar_id 总数（聚合后特征数）: {n_var:,}")

# %% [markdown]
# ## Step C: 构 var DataFrame
#
# 每个 var_id 一行，列与 MicrobeAtlas 一致。同一 var_id 下所有 ASV 的前 6 级
# 完全相同（构造方式保证），取首条即可。

# %%
var_list = sorted(tax_kept['var_id'].unique())
var_to_idx = {v: i for i, v in enumerate(var_list)}

# 直接从 var_id 拆 6 列，避免 h5ad Categorical 列 fillna 报错；
# make_var_id 已在缺位处填好层级前缀（d__/p__/...），拆回来即所需占位符。
var_df = pd.DataFrame(
    [v.split(';') for v in var_list],
    columns=RANK_TO_GENUS,
    index=pd.Index(var_list, name='var_id'),
).astype(str)
print(f"var_df shape: {var_df.shape}")
var_df.head(3)

# %% [markdown]
# ## Step D: 子集 X + 构聚合矩阵 → 输出 AnnData

# %%
asv_keep_idx = np.where(keep_mask)[0]
print(f"保留 ASV 索引数: {len(asv_keep_idx):,}")

X_kept = adata_in.X[:, asv_keep_idx].tocsc().tocsr()
print(f"X_kept: {X_kept.shape}, nnz={X_kept.nnz:,}")

# 聚合矩阵 M: (n_var, n_kept_asv)
var_idx_per_asv = np.array(
    [var_to_idx[v] for v in tax_kept['var_id']],
    dtype=np.int64
)
n_kept = len(asv_keep_idx)
M = sp.csr_matrix(
    (np.ones(n_kept, dtype=np.int32),
     (var_idx_per_asv, np.arange(n_kept))),
    shape=(n_var, n_kept),
)
print(f"聚合矩阵 M: {M.shape}, nnz={M.nnz:,}")

print("聚合中 ...")
X_var = (X_kept @ M.T).tocsr()
max_val = int(X_var.data.max()) if X_var.nnz else 0
print(f"聚合后单格最大值: {max_val:,}  (int32 上限 2,147,483,647)")
assert max_val < 2_147_483_647, "聚合值超 int32"
X_var.data = X_var.data.astype(np.int32, copy=False)
print(f"sample × var: {X_var.shape}, nnz={X_var.nnz:,}")

# %% [markdown]
# ### 丢弃 reads 统计 (含按项目分组)

# %%
sample_total_full = np.asarray(adata_in.X.sum(axis=1)).ravel()
sample_total_kept = np.asarray(X_var.sum(axis=1)).ravel()
total_full = int(sample_total_full.sum())
total_kept = int(sample_total_kept.sum())
print(f"丢弃 reads: {total_full - total_kept:,} / {total_full:,} "
      f"({(total_full - total_kept) / total_full * 100:.2f}%)")

nz = sample_total_full > 0
loss_frac = np.zeros_like(sample_total_full, dtype=float)
loss_frac[nz] = 1 - sample_total_kept[nz] / sample_total_full[nz]
print("\n=== per-sample 丢失比例分布 ===")
for q in [0.5, 0.75, 0.9, 0.95, 0.99]:
    print(f"  {int(q*100):>3d}% 分位: {np.quantile(loss_frac, q)*100:7.4f}%")
print(f"  最大       : {loss_frac.max()*100:7.4f}%")
print(f"  >50% 的样本: {(loss_frac > 0.5).sum():,} / {len(loss_frac):,}")

# 项目维度（loss 主要来自 mito/chloro，按项目集中度可揭示宿主污染）
proj_per_sample = adata_in.var['project'].iloc[adata_in.X.indices] if False else None  # 占位，下面用更高效的方式
# 用 ASV-level project + count 加权
asv_project = adata_in.var['project'].values
sample_total_kept_by_proj = pd.Series(
    np.asarray((adata_in.X[:, asv_keep_idx]).sum(axis=0)).ravel(),
    index=tax_kept['project'].values,
).groupby(level=0).sum()
sample_total_full_by_proj = pd.Series(
    np.asarray(adata_in.X.sum(axis=0)).ravel(),
    index=asv_project,
).groupby(level=0).sum()
proj_loss = (1 - sample_total_kept_by_proj.reindex(sample_total_full_by_proj.index)
                  / sample_total_full_by_proj).sort_values(ascending=False)
print("\n=== 丢失比例最高的 10 个项目 ===")
for proj, loss in proj_loss.head(10).items():
    print(f"  {proj:<14}  {loss*100:6.2f}%")

# %% [markdown]
# ## Step E: 写出

# %%
adata_out = ad.AnnData(
    X   = X_var,
    obs = pd.DataFrame(index=adata_in.obs.index.copy()),
    var = var_df.copy(),
)
print(adata_out)

OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
adata_out.write_h5ad(OUT_PATH, compression='gzip')
print(f"\n已保存: {OUT_PATH}")

# %% [markdown]
# ## 完成
#
# 下一步: `04_drop_empty_samples.py` 剔除聚合后零 taxon 的样本。
