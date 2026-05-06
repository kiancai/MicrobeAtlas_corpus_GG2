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
# # OTU → Genus AnnData 构建（信息无损方案）
#
# **输入**
# - `results/taxonomy/taxonomy.tsv` — GG2 NB 注释结果（111,870 OTU）
# - `rawdata/.../samples-otus.97.mapped.biom.gz` — 全量 BIOM (269 万样本)
# - `rawdata/.../samples-otus.97.mapped.metag.minfilter.refilt.biom.gz` — 过滤 BIOM (188 万)
#
# **输出** (`results/feature_table/`)
# - `otu_taxonomy_full.tsv` — 全部 OTU 的注释 + var_id 映射（含 Confidence，下游可自行 filter）
# - `var_summary.tsv` — 每个 var_id 包含的 OTU 数 + counts 占比 + 最深可用层级
# - `gg2.full.h5ad` — sample × var AnnData (full BIOM)
# - `gg2.minfilter.h5ad` — sample × var AnnData (minfilter BIOM)
#
# **聚合 key（var_id）= QIIME2 风格的"前 6 级完整路径"**
# - 完整到 genus: `d__Bacteria;p__Firmicutes;...;f__Streptococcaceae;g__Streptococcus`
# - 仅到 family : `d__Bacteria;p__Firmicutes;...;f__Streptococcaceae;g__`
# - 仅到 phylum : `d__Bacteria;p__Firmicutes;c__;o__;f__;g__`
# - 完全 Unassigned: `Unassigned`
#
# 这样**同名 genus 跨 family**会被分开（GG2 没有这种情况，但更稳）；不同 family 下的 `g__` 也不会错合并。
#
# **AnnData 内容**
# - `X`: scipy.sparse.csr_matrix, dtype=int32, raw counts（**不归一化**）
# - `obs`: 空 DataFrame，index = sample IDs（不读 metadata）
# - `var`: DataFrame, index = var_id, 列 = [Domain, Phylum, Class, Order, Family, Genus, deepest_rank, is_resolved_genus, n_otu]
#
# **不做的事**
# - Confidence 过滤（保留全部 111,870 OTU，Confidence 信息存在 OTU 映射表里）
# - Domain 过滤（保留全部，含 Unassigned）
# - 归一化、`layers['counts']`、metadata 合并（后续单独处理）
#
# **环境**: `baseBio` (h5py, pandas, scipy, anndata, jupytext)

# %%
import os
import subprocess
import tempfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

PROJECT_DIR  = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
TAX_TSV      = PROJECT_DIR / "results/taxonomy/taxonomy.tsv"
BIOM_FULL_GZ = PROJECT_DIR / "rawdata/MicrobeAtlas/OTU_count/samples-otus.97.mapped.biom.gz"
BIOM_MIN_GZ  = PROJECT_DIR / "rawdata/MicrobeAtlas/OTU_count/samples-otus.97.mapped.metag.minfilter.refilt.biom.gz"
OUT_DIR      = PROJECT_DIR / "results/feature_table"
OUT_DIR.mkdir(parents=True, exist_ok=True)

print(f"Project: {PROJECT_DIR}")
print(f"Output : {OUT_DIR}")

# %% [markdown]
# ## Step A: 读 taxonomy + 拆 7 级 + 计算 var_id

# %%
RANK_COLS = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
PREFIXES  = ['d__',    'p__',    'c__',   'o__',   'f__',    'g__',   's__']
RANK_TO_GENUS = RANK_COLS[:6]      # Domain..Genus，构造 var_id 用
PFX_TO_GENUS  = PREFIXES[:6]

def parse_gg2(taxon_str):
    """按 7 级前缀拆解 GG2 taxonomy 字符串。无对应前缀返回 None。"""
    parts = [p.strip() for p in str(taxon_str).split(';')]
    return {col: next((p for p in parts if p.startswith(pfx)), None)
            for col, pfx in zip(RANK_COLS, PREFIXES)}

def deepest_rank(row):
    """返回最深可用层级名，如 'Genus' / 'Family' / 'None'。"""
    for col, pfx in zip(reversed(RANK_COLS), reversed(PREFIXES)):
        v = row[col]
        if v is not None and v != pfx:
            return col
    return 'None'

def make_var_id(row):
    """QIIME2 风格的前 6 级完整路径作为聚合 key。
    完全 Unassigned (无任何层级) → 'Unassigned'。"""
    if all(row[c] is None for c in RANK_TO_GENUS):
        return 'Unassigned'
    parts = []
    for col, pfx in zip(RANK_TO_GENUS, PFX_TO_GENUS):
        v = row[col]
        parts.append(v if v is not None else pfx)
    return ';'.join(parts)

tax = pd.read_csv(TAX_TSV, sep='\t', index_col=0)
print(f"原始 OTU 数: {len(tax):,}")
print(f"taxonomy.tsv 列: {list(tax.columns)}")

ranks = tax['Taxon'].apply(parse_gg2).apply(pd.Series)
tax = pd.concat([tax, ranks], axis=1)
tax['deepest_rank']      = tax.apply(deepest_rank, axis=1)
tax['var_id']            = tax.apply(make_var_id, axis=1)
tax['is_resolved_genus'] = tax['deepest_rank'].isin(['Genus', 'Species'])

tax.head(3)

# %% [markdown]
# ### 注释深度 + var_id 多样性

# %%
RANK_ORDER = ['Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Domain', 'None']
print("=== OTU 按最深可用层级分布 ===")
depth_dist = tax['deepest_rank'].value_counts().reindex(RANK_ORDER, fill_value=0)
for r in RANK_ORDER:
    n = depth_dist[r]
    print(f"  {r:<10} {n:>10,}  ({n/len(tax)*100:5.2f}%)")

print()
print(f"=== var_id 总数（聚合后的特征数）===")
print(f"  unique var_id: {tax['var_id'].nunique():,}")
print(f"  其中 is_resolved_genus=True 的 var_id: {tax[tax['is_resolved_genus']]['var_id'].nunique():,}")
print(f"  其中 is_resolved_genus=False 的 var_id: {tax[~tax['is_resolved_genus']]['var_id'].nunique():,}")

# %% [markdown]
# ## Step B: 构建 var DataFrame
#
# 每个 var_id 一行。同一 var_id 下所有 OTU 必然有相同的前 6 级注释（因为 var_id 就是这 6 级），所以 6 级列取首条即可。

# %%
var_list = sorted(tax['var_id'].unique())
var_to_idx = {v: i for i, v in enumerate(var_list)}
n_var = len(var_list)
print(f"var 总数: {n_var:,}")

# n_otu 列：每个 var_id 包含的 OTU 数（聚合度量）
n_otu_per_var = tax.groupby('var_id').size().rename('n_otu')

var_df = (
    tax.drop_duplicates(subset='var_id', keep='first')
       .set_index('var_id')
       [RANK_TO_GENUS + ['deepest_rank', 'is_resolved_genus']]
       .join(n_otu_per_var)
       .loc[var_list]
)
print(f"var_df shape: {var_df.shape}")
var_df.head(3)

# %% [markdown]
# ### var_id 聚合度诊断（看哪些 var_id 是"大杂烩"）

# %%
print("Top 10 包含最多 OTU 的 var_id:")
top10 = var_df.sort_values('n_otu', ascending=False).head(10)
for vid, row in top10.iterrows():
    print(f"  n_otu={row['n_otu']:>5}  deepest={row['deepest_rank']:<8}  {vid[:100]}")

print()
print("=== var_id 包含 OTU 数的分布 ===")
print(var_df['n_otu'].describe().to_string())

# %% [markdown]
# ## Step C: 保存 OTU 级映射表

# %%
out_cols = ['Taxon', 'Confidence'] + RANK_COLS + ['deepest_rank', 'var_id', 'is_resolved_genus']
map_path = OUT_DIR / "otu_taxonomy_full.tsv"
tax[out_cols].to_csv(map_path, sep='\t')
print(f"已保存 OTU 级映射: {map_path}  ({len(tax):,} 行)")

summary_path = OUT_DIR / "var_summary.tsv"
var_df.to_csv(summary_path, sep='\t')
print(f"已保存 var 汇总: {summary_path}  ({len(var_df):,} 行)")

# %% [markdown]
# ## Step D: BIOM 读取助手
#
# 解压 `.biom.gz` → 读取 HDF5 为 (OTU × sample) CSR 稀疏矩阵。
#
# > 注：解压后可能 10–30 GB，确保 `$TMPDIR` / `/tmp` 足够。

# %%
def biom_gz_to_otu_csr(biom_gz_path):
    """解压 .biom.gz 到临时文件，读取为 OTU × sample 的 CSR 矩阵。"""
    print(f"  解压 {Path(biom_gz_path).name} ...")
    with tempfile.NamedTemporaryFile(suffix='.biom', delete=False) as tmp:
        tmp_path = tmp.name
    with open(tmp_path, 'wb') as fout:
        subprocess.run(['zcat', str(biom_gz_path)], stdout=fout, check=True)

    print(f"  读取 HDF5 ...")
    with h5py.File(tmp_path, 'r') as f:
        obs_ids    = np.array([x.decode() for x in f['observation/ids'][:]])
        sample_ids = np.array([x.decode() for x in f['sample/ids'][:]])
        data       = f['observation/matrix/data'][:]
        indices    = f['observation/matrix/indices'][:]
        indptr     = f['observation/matrix/indptr'][:]
    os.unlink(tmp_path)

    mat = sp.csr_matrix((data, indices, indptr),
                        shape=(len(obs_ids), len(sample_ids)))
    print(f"  OTU × sample = {mat.shape}, nnz = {mat.nnz:,}")
    return mat, obs_ids, sample_ids

# %% [markdown]
# ## Step E-1: 处理 full BIOM
#
# `samples-otus.97.mapped.biom.gz`（269 万样本）

# %%
mat_full, otu_ids_full, sample_ids_full = biom_gz_to_otu_csr(BIOM_FULL_GZ)

# %%
# 在 BIOM OTU 中找出与 taxonomy 的交集
otu_pos_full = pd.Series(np.arange(len(otu_ids_full)), index=otu_ids_full)
kept_in_biom = tax.index.intersection(otu_ids_full)
missing      = tax.index.difference(otu_ids_full)

print(f"taxonomy OTU 总数: {len(tax):,}")
print(f"BIOM OTU 总数    : {len(otu_ids_full):,}")
print(f"交集             : {len(kept_in_biom):,}")
if len(missing) > 0:
    print(f"NOTE: {len(missing):,} 个 taxonomy OTU 不在 BIOM 中（参考库 OTU 但未被检出）")

idx_subset = otu_pos_full.loc[kept_in_biom].values
mat_kept   = mat_full[idx_subset]
print(f"\n子集后矩阵: {mat_kept.shape}, nnz = {mat_kept.nnz:,}")

del mat_full

# %%
# 构建聚合矩阵 M: (n_var, n_kept_otu)，M[v, o] = 1 当且仅当 OTU o 属于 var_id v
var_idx_per_otu = np.array(
    [var_to_idx[v] for v in tax.loc[kept_in_biom, 'var_id']],
    dtype=np.int64
)
n_kept = len(kept_in_biom)

M = sp.csr_matrix(
    (np.ones(n_kept, dtype=np.int32),
     (var_idx_per_otu, np.arange(n_kept))),
    shape=(n_var, n_kept)
)
print(f"聚合矩阵 M: {M.shape}, nnz = {M.nnz:,}")

print("聚合中 ...")
mat_var_full      = (M @ mat_kept).T.tocsr()
mat_var_full.data = mat_var_full.data.astype(np.int32)
print(f"sample × var: {mat_var_full.shape}, nnz = {mat_var_full.nnz:,}")

del mat_kept, M

# %%
adata_full = ad.AnnData(
    X   = mat_var_full,
    obs = pd.DataFrame(index=pd.Index(sample_ids_full, name='sample_id')),
    var = var_df.copy()
)
print(adata_full)

out_path = OUT_DIR / "gg2.full.h5ad"
adata_full.write_h5ad(out_path, compression='gzip')
print(f"\n已保存: {out_path}")

del mat_var_full, adata_full

# %% [markdown]
# ## Step E-2: 处理 minfilter BIOM
#
# `samples-otus.97.mapped.metag.minfilter.refilt.biom.gz`（188 万样本）。
#
# 逻辑与 E-1 完全一致，重新展开一遍以便看到中间结果。

# %%
mat_min, otu_ids_min, sample_ids_min = biom_gz_to_otu_csr(BIOM_MIN_GZ)

# %%
otu_pos_min    = pd.Series(np.arange(len(otu_ids_min)), index=otu_ids_min)
kept_in_biom_m = tax.index.intersection(otu_ids_min)
missing_m      = tax.index.difference(otu_ids_min)

print(f"taxonomy OTU 总数: {len(tax):,}")
print(f"BIOM OTU 总数    : {len(otu_ids_min):,}")
print(f"交集             : {len(kept_in_biom_m):,}")
if len(missing_m) > 0:
    print(f"NOTE: {len(missing_m):,} 个 taxonomy OTU 不在 BIOM 中")

idx_subset_m = otu_pos_min.loc[kept_in_biom_m].values
mat_kept_m   = mat_min[idx_subset_m]
print(f"\n子集后矩阵: {mat_kept_m.shape}, nnz = {mat_kept_m.nnz:,}")

del mat_min

# %%
var_idx_per_otu_m = np.array(
    [var_to_idx[v] for v in tax.loc[kept_in_biom_m, 'var_id']],
    dtype=np.int64
)
n_kept_m = len(kept_in_biom_m)

M_m = sp.csr_matrix(
    (np.ones(n_kept_m, dtype=np.int32),
     (var_idx_per_otu_m, np.arange(n_kept_m))),
    shape=(n_var, n_kept_m)
)
print(f"聚合矩阵 M: {M_m.shape}, nnz = {M_m.nnz:,}")

print("聚合中 ...")
mat_var_min      = (M_m @ mat_kept_m).T.tocsr()
mat_var_min.data = mat_var_min.data.astype(np.int32)
print(f"sample × var: {mat_var_min.shape}, nnz = {mat_var_min.nnz:,}")

del mat_kept_m, M_m

# %%
adata_min = ad.AnnData(
    X   = mat_var_min,
    obs = pd.DataFrame(index=pd.Index(sample_ids_min, name='sample_id')),
    var = var_df.copy()
)
print(adata_min)

out_path = OUT_DIR / "gg2.minfilter.h5ad"
adata_min.write_h5ad(out_path, compression='gzip')
print(f"\n已保存: {out_path}")

# %% [markdown]
# ## 完成
#
# **输出文件**
# - `results/feature_table/otu_taxonomy_full.tsv` — OTU 级映射，含 Confidence、deepest_rank、var_id
# - `results/feature_table/var_summary.tsv` — var 汇总
# - `results/feature_table/gg2.full.h5ad` — sample × var AnnData (full)
# - `results/feature_table/gg2.minfilter.h5ad` — sample × var AnnData (minfilter)
#
# **下游用法（示例）**
#
# ```python
# import anndata as ad
# adata = ad.read_h5ad("results/feature_table/gg2.full.h5ad")
#
# # 只要"真 genus"特征
# adata_genus = adata[:, adata.var['is_resolved_genus']].copy()
#
# # 只要 confidence ≥ 0.7 的 OTU 聚合（需重跑或在 OTU 映射表里 filter 后重聚合）
# # —— 这一步 confidence 是 OTU 级的，已在 otu_taxonomy_full.tsv 中
# ```
