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
# # OTU → Genus AnnData 构建
#
# **输入**
# - `results/taxonomy/taxonomy.tsv` — GG2 NB 注释结果（111,870 OTU）
# - `rawdata/.../samples-otus.97.mapped.biom.gz` — 全量 BIOM (269 万样本)
# - `rawdata/.../samples-otus.97.mapped.metag.minfilter.refilt.biom.gz` — 过滤 BIOM (188 万)
#
# **输出** (`results/feature_table/`)
# - `gg2.full.h5ad` — sample × var AnnData (full BIOM)
# - `gg2.minfilter.h5ad` — sample × var AnnData (minfilter BIOM)
#
# **过滤标准（参考 easyAmp.r / Yong-Xin Liu 16S 流程）**
# 1. Domain ∈ {Bacteria, Archaea} —— 丢 Unassigned 与 Eukaryota
# 2. 丢 mitochondria（taxonomy 字符串含 'mitochondri'，不区分大小写）
# 3. 丢 chloroplast（taxonomy 字符串含 'chloroplast'，不区分大小写）
#
# 之所以源头过滤 mito/chloro：GG2 沿 GTDB 把它们嵌在 `d__Bacteria` 下面，
# 共 ~1,750 个 var_id，看上去跟普通细菌无异，下游做 phylum/class 级丰度统计
# 时若忘记过滤会污染 Cyanobacteriota / Alphaproteobacteria 的结果。
#
# 不用 Confidence 阈值：实测它命中的 OTU 全部已被 Domain 过滤覆盖。
#
# **聚合 key（var_id）= QIIME2 风格的"前 6 级完整路径"**
# - 完整到 genus: `d__Bacteria;p__Firmicutes;...;f__Streptococcaceae;g__Streptococcus`
# - 仅到 family : `d__Bacteria;p__Firmicutes;...;f__Streptococcaceae;g__`
# - 仅到 phylum : `d__Bacteria;p__Firmicutes;c__;o__;f__;g__`
#
# 浅层注释（任意深度未到 genus）的 OTU **不丢弃**，通过空占位符保留其 counts。
#
# **AnnData 内容**
# - `X`: scipy.sparse.csr_matrix, dtype=int32, raw counts（不归一化）
# - `obs`: 空 DataFrame，index = sample IDs（不读 metadata）
# - `var`: DataFrame, index = var_id, 列 = [Domain, Phylum, Class, Order, Family, Genus]
#
# **不做的事**
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
    过滤 Domain∈{B,A} 后不会出现完全 Unassigned 的行。"""
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
tax['deepest_rank'] = tax.apply(deepest_rank, axis=1)
tax['var_id']       = tax.apply(make_var_id, axis=1)

tax.head(3)

# %% [markdown]
# ### 应用过滤：Domain ∈ {B,A}、去 mito/chloro

# %%
print(f"过滤前 OTU 数: {len(tax):,}")

mask_dom    = tax['Domain'].isin(['d__Bacteria', 'd__Archaea'])
mask_mito   = ~tax['Taxon'].str.contains('mitochondri', case=False, na=False)
mask_chloro = ~tax['Taxon'].str.contains('chloroplast',  case=False, na=False)

print(f"  丢 Domain 非 B/A : {(~mask_dom).sum():,}")
print(f"  丢 mitochondria  : {(~mask_mito).sum():,}")
print(f"  丢 chloroplast   : {(~mask_chloro).sum():,}")

tax = tax[mask_dom & mask_mito & mask_chloro].copy()
print(f"过滤后 OTU 数: {len(tax):,}")

RANK_ORDER = ['Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Domain', 'None']
print("\n=== OTU 按最深可用层级分布 ===")
depth_dist = tax['deepest_rank'].value_counts().reindex(RANK_ORDER, fill_value=0)
for r in RANK_ORDER:
    n = depth_dist[r]
    print(f"  {r:<10} {n:>10,}  ({n/len(tax)*100:5.2f}%)")

print(f"\nvar_id 总数（聚合后特征数）: {tax['var_id'].nunique():,}")

# %% [markdown]
# ## Step B: 构建 var DataFrame
#
# 每个 var_id 一行，只保留 6 级注释。同一 var_id 下所有 OTU 的前 6 级完全一致，取首条即可。

# %%
var_list = sorted(tax['var_id'].unique())
var_to_idx = {v: i for i, v in enumerate(var_list)}
n_var = len(var_list)
print(f"var 总数: {n_var:,}")

var_df = (
    tax.drop_duplicates(subset='var_id', keep='first')
       .set_index('var_id')
       [RANK_TO_GENUS]
       .loc[var_list]
       .fillna('')
)
print(f"var_df shape: {var_df.shape}")
var_df.head(3)

# %% [markdown]
# ## Step C: BIOM 读取助手
#
# 解压 `.biom.gz` → 读取 HDF5 为 (OTU × sample) CSR 稀疏矩阵。
#
# > 注：解压后可能 10–30 GB，确保 `$TMPDIR` / `/tmp` 足够。

# %%
def biom_gz_to_otu_csr(biom_gz_path):
    """解压 .biom.gz 到临时文件，读取为 OTU × sample 的 CSR 矩阵。"""
    print(f"  解压 {Path(biom_gz_path).name} ...")
    with tempfile.NamedTemporaryFile(suffix='.biom', delete=False,
                                      dir=PROJECT_DIR) as tmp:
        tmp_path = tmp.name
    try:
        with open(tmp_path, 'wb') as fout:
            subprocess.run(['zcat', str(biom_gz_path)], stdout=fout, check=True)

        print(f"  读取 HDF5 ...")
        with h5py.File(tmp_path, 'r') as f:
            obs_ids    = np.array([x.decode() for x in f['observation/ids'][:]])
            sample_ids = np.array([x.decode() for x in f['sample/ids'][:]])
            data       = f['observation/matrix/data'][:]
            indices    = f['observation/matrix/indices'][:]
            indptr     = f['observation/matrix/indptr'][:]
    finally:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)

    mat = sp.csr_matrix((data, indices, indptr),
                        shape=(len(obs_ids), len(sample_ids)))
    print(f"  OTU × sample = {mat.shape}, nnz = {mat.nnz:,}")
    return mat, obs_ids, sample_ids

# %% [markdown]
# ## Step D-1: 处理 full BIOM
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

# 总体丢失统计
sample_total_full = np.asarray(mat_full.sum(axis=0)).ravel()
sample_total_kept = np.asarray(mat_kept.sum(axis=0)).ravel()
total_reads       = int(sample_total_full.sum())
kept_reads        = int(sample_total_kept.sum())
dropped_reads     = total_reads - kept_reads
print(f"丢弃 reads: {dropped_reads:,} / {total_reads:,} "
      f"({dropped_reads/total_reads*100:.2f}%)")

# per-sample 丢失比例分布（更能区分"少数样本拖累池子" vs "普遍丢失"）
nz = sample_total_full > 0
loss_frac = np.zeros_like(sample_total_full, dtype=float)
loss_frac[nz] = 1 - sample_total_kept[nz] / sample_total_full[nz]
print("\n=== per-sample 丢失比例分布 ===")
for q in [0.5, 0.75, 0.9, 0.95, 0.99]:
    print(f"  {int(q*100):>3d}% 分位: {np.quantile(loss_frac, q)*100:7.4f}%")
print(f"  最大       : {loss_frac.max()*100:7.4f}%")
print(f"  >50% 的样本: {(loss_frac > 0.5).sum():,} / {len(loss_frac):,}")

del mat_full, sample_total_full, sample_total_kept, loss_frac

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
max_val = int(mat_var_full.data.max())
print(f"聚合后单格最大值: {max_val:,}  (int32 上限 2,147,483,647)")
assert max_val < 2_147_483_647, "聚合值超 int32，需要改 dtype 或重新审视数据"
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
# ## Step D-2: 处理 minfilter BIOM
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

sample_total_full_m = np.asarray(mat_min.sum(axis=0)).ravel()
sample_total_kept_m = np.asarray(mat_kept_m.sum(axis=0)).ravel()
total_reads_m       = int(sample_total_full_m.sum())
kept_reads_m        = int(sample_total_kept_m.sum())
dropped_reads_m     = total_reads_m - kept_reads_m
print(f"丢弃 reads: {dropped_reads_m:,} / {total_reads_m:,} "
      f"({dropped_reads_m/total_reads_m*100:.2f}%)")

nz_m = sample_total_full_m > 0
loss_frac_m = np.zeros_like(sample_total_full_m, dtype=float)
loss_frac_m[nz_m] = 1 - sample_total_kept_m[nz_m] / sample_total_full_m[nz_m]
print("\n=== per-sample 丢失比例分布 ===")
for q in [0.5, 0.75, 0.9, 0.95, 0.99]:
    print(f"  {int(q*100):>3d}% 分位: {np.quantile(loss_frac_m, q)*100:7.4f}%")
print(f"  最大       : {loss_frac_m.max()*100:7.4f}%")
print(f"  >50% 的样本: {(loss_frac_m > 0.5).sum():,} / {len(loss_frac_m):,}")

del mat_min, sample_total_full_m, sample_total_kept_m, loss_frac_m

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
max_val_m = int(mat_var_min.data.max())
print(f"聚合后单格最大值: {max_val_m:,}  (int32 上限 2,147,483,647)")
assert max_val_m < 2_147_483_647, "聚合值超 int32，需要改 dtype 或重新审视数据"
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
# - `results/feature_table/gg2.full.h5ad` — sample × var AnnData (full)
# - `results/feature_table/gg2.minfilter.h5ad` — sample × var AnnData (minfilter)
#
# **下游用法（示例）**
#
# ```python
# import anndata as ad
# adata = ad.read_h5ad("results/feature_table/gg2.full.h5ad")
#
# # 已在源头过滤 Unassigned / Eukaryota / mito / chloro，可直接使用
#
# # 只保留"真 genus"特征（Genus 列非空占位符）
# genus_mask = adata.var['Genus'].str.len() > 3  # 'g__' 是 3 字符
# adata_genus = adata[:, genus_mask].copy()
# ```
