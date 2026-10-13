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
# # MicroBioMap 05: ASV → Genus AnnData  (= ResMicroDb 03 逐字逻辑)
#
# **这是交付物**：GG2 24.09 重新注释的**完整物种组成表**（genus 级）。
# 逻辑与 ResMicroDb 03 / MicrobeAtlas 03 **完全一致**，只改输入/输出路径。
#
# **输入** `results/feature_table/microbiomap.gg2.asv.h5ad`（04 输出）
# - sample × ASV，var 含 `gg2_Domain..gg2_Species` + `gg2_Confidence` + `project`
#
# **输出** `results/feature_table/microbiomap.gg2.genus.h5ad`
# - sample(SRR) × genus_var（168,464 行全保留）
# - var 列与 MA/RM 一致：`[Domain, Phylum, Class, Order, Family, Genus]`（字符串，缺位空占位符）
#
# **过滤逻辑（与 RM/MA 03 完全一致，仅 taxon 级清理，不动样本）**
# 1. `gg2_Domain ∈ {d__Bacteria, d__Archaea}` —— 丢 Unassigned / Eukaryota
# 2. 丢 mitochondria / chloroplast（Taxon 含相应字符串）
# 浅层注释（未到 genus）的 ASV **不丢弃**，通过空占位符 `g__` 保留 counts。
#
# **到此为止**：sample 级筛选（剔零 / QC）与合并属"待商议"，本步不做。
# 注：MBM ASV 为全局共享 → `project` 列恒为 'MicroBioMap'，故"按项目丢失率"
# 诊断只会显示一组（无意义但无害）；若需按真实项目可由 sample_metadata 派生。

# %%
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

PROJECT_DIR = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
IN_PATH     = PROJECT_DIR / "results/feature_table/microbiomap.gg2.asv.h5ad"
OUT_PATH    = PROJECT_DIR / "results/feature_table/microbiomap.gg2.genus.h5ad"

# %% [markdown]
# ## Step A: 读 04 输出 + 重建 6 级 var_id（与 RM/MA 03 同款常量与函数）

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
mask_mito   = ~tax['Taxon_full'].str.contains(r'(?:^|;)\s*[dkpcofgs]__Mitochondria\s*(?=;|$)', case=False, na=False)
mask_chloro = ~tax['Taxon_full'].str.contains(r'(?:^|;)\s*[dkpcofgs]__Chloroplast\s*(?=;|$)', case=False, na=False)

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

# %%
var_list = sorted(tax_kept['var_id'].unique())
var_to_idx = {v: i for i, v in enumerate(var_list)}

var_df = pd.DataFrame(
    [v.split(';') for v in var_list],
    columns=RANK_TO_GENUS,
    index=pd.Index(var_list, name='var_id'),
).astype(str)
print(f"var_df shape: {var_df.shape}")

# %% [markdown]
# ## Step D: 子集 X + 构聚合矩阵 → 输出 AnnData

# %%
asv_keep_idx = np.where(keep_mask)[0]
print(f"保留 ASV 索引数: {len(asv_keep_idx):,}")

X_kept = adata_in.X[:, asv_keep_idx].tocsc().tocsr()
print(f"X_kept: {X_kept.shape}, nnz={X_kept.nnz:,}")

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
# ### 丢弃 reads 统计

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
print(f"  零 taxon 样本(聚合后全 0): {(sample_total_kept == 0).sum():,} / {len(loss_frac):,}  (本步不剔，留待筛选讨论)")

# %% [markdown]
# ## Step E: 写出（全样本保留）

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
print("\n=== 这是交付物：GG2 24.09 重注释完整物种组成表（genus 级，全样本）===")
print("后续 sample 级筛选 / metadata / 合并 —— 待与用户商议，本流程不做。")
