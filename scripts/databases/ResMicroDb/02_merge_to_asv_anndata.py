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
# # ResMicroDb 02: 合并 398 项目 → sample × ASV AnnData (双注释)
#
# **输入** (每项目 `rawdata/ResMicroDb/16S/<PROJECT>/results/`)
# - `asv.fa`             — ASV 序列 (项目内 ID `ASV_N`)
# - `otutab.txt`         — ASV × sample dense TSV，第一列 `#OTU ID`
# - `taxonomy_gg2.txt`   — 01 步生成: 3 列 `Feature ID / Taxon / Confidence`
# - `taxonomy_silva.txt` — 原有 8 列: `OTUID / Kingdom..Species`
#
# **输出** `results/feature_table/resmicrodb.gg2.asv.h5ad`
# - `X`: sample × ASV  (CSR, int32)
# - `obs`: index = sample_id (空 DataFrame)
# - `var`: index = `<PROJECT>__<ASV>`，列 = `[project, gg2_*, gg2_Confidence, silva_*]`
#
# **不做的事**
# - 不丢任何 ASV / 不过滤 mito/chloro / 不归一化 / 不聚合 → 这些在 03 做

# %%
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

PROJECT_DIR = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
DATA_DIR    = PROJECT_DIR / "rawdata/ResMicroDb/16S"
OUT_DIR     = PROJECT_DIR / "results/feature_table"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_PATH = OUT_DIR / "resmicrodb.gg2.asv.h5ad"

# %% [markdown]
# ## Step A: 列项目并核对输入完整性

# %%
projects = sorted(p.name for p in DATA_DIR.iterdir()
                  if (p / "results/asv.fa").is_file())
print(f"有 asv.fa 的项目: {len(projects)}")

REQUIRED = ["asv.fa", "otutab.txt", "taxonomy_gg2.txt", "taxonomy_silva.txt"]
missing = {f: [] for f in REQUIRED}
for proj in projects:
    for f in REQUIRED:
        if not (DATA_DIR / proj / "results" / f).is_file():
            missing[f].append(proj)

for f, lst in missing.items():
    print(f"  缺 {f:<22}: {len(lst)}")
    if lst:
        print(f"    示例: {lst[:5]}")

if missing["taxonomy_gg2.txt"]:
    raise RuntimeError(
        f"{len(missing['taxonomy_gg2.txt'])} 个项目还没跑完 01 (缺 taxonomy_gg2.txt)，先重投这些任务"
    )

# %% [markdown]
# ## Step B: GG2 Taxon 解析助手
#
# 与 MicrobeAtlas 03 中的 `parse_gg2` 一致；这里返回 dict 直接展开成列。

# %%
RANK_COLS = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
PREFIXES  = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']


def parse_gg2(taxon_str):
    """按 7 级前缀拆 GG2 taxonomy 字符串，缺位返回 None。"""
    if pd.isna(taxon_str):
        return {c: None for c in RANK_COLS}
    parts = [p.strip() for p in str(taxon_str).split(';')]
    return {c: next((p for p in parts if p.startswith(pfx)), None)
            for c, pfx in zip(RANK_COLS, PREFIXES)}


# %% [markdown]
# ## Step C: 逐项目处理 → 累积 COO 三元组
#
# 用 (row=sample_global_idx, col=asv_global_idx, val=count) 在每项目内
# 计算偏移后 append。最后一次性 `coo_matrix → csr`。
#
# 同时累积每项目的 var DataFrame，最后 concat。

# %%
sample_id_list   = []   # 全局 sample order
asv_id_list      = []   # 全局 ASV order (已 namespace)
var_dfs          = []   # 每项目一份 var DataFrame
coo_rows_list    = []
coo_cols_list    = []
coo_data_list    = []

sample_offset = 0
asv_offset    = 0
total_reads   = 0
total_nnz     = 0

SILVA_RANKS = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

for i, proj in enumerate(projects, 1):
    rdir = DATA_DIR / proj / "results"

    otu = pd.read_csv(rdir / "otutab.txt", sep='\t', index_col=0)
    otu.index.name = 'asv_id'
    otu = otu.astype(np.int32, copy=False)

    asv_local    = otu.index.tolist()
    sample_local = otu.columns.tolist()
    n_asv, n_sample = otu.shape

    # 防御 1: QIIME2 某些版本 export 的 taxonomy.tsv 第二行是 #q2:types 元数据，跳过之
    gg2 = pd.read_csv(rdir / "taxonomy_gg2.txt", sep='\t')
    gg2 = gg2.rename(columns={gg2.columns[0]: 'asv_id'}).set_index('asv_id')
    if str(gg2.index[0]).startswith('#'):
        gg2 = gg2.iloc[1:]
    gg2['Confidence'] = pd.to_numeric(gg2['Confidence'], errors='coerce')
    gg2_ranks = gg2['Taxon'].apply(parse_gg2).apply(pd.Series)
    gg2_ranks.columns = [f'gg2_{c}' for c in RANK_COLS]
    gg2_ranks['gg2_Confidence'] = gg2['Confidence'].astype(np.float32)

    silva = pd.read_csv(rdir / "taxonomy_silva.txt", sep='\t')
    silva = silva.rename(columns={silva.columns[0]: 'asv_id'}).set_index('asv_id')
    if str(silva.index[0]).startswith('#'):
        silva = silva.iloc[1:]
    silva = silva.rename(columns={c: f'silva_{c}' for c in silva.columns})

    # 防御 2: otutab.txt / taxonomy_gg2.txt / taxonomy_silva.txt 的 ASV ID 集合一致性
    set_otu, set_gg2, set_silva = set(asv_local), set(gg2.index), set(silva.index)
    missing_in_gg2   = set_otu - set_gg2
    missing_in_silva = set_otu - set_silva
    extra_in_gg2     = set_gg2 - set_otu
    extra_in_silva   = set_silva - set_otu
    if missing_in_gg2 or missing_in_silva or extra_in_gg2 or extra_in_silva:
        raise RuntimeError(
            f"项目 {proj} ASV ID 集合不一致:\n"
            f"  otu - gg2  : {len(missing_in_gg2):>6} (例 {list(missing_in_gg2)[:3]})\n"
            f"  otu - silva: {len(missing_in_silva):>6} (例 {list(missing_in_silva)[:3]})\n"
            f"  gg2 - otu  : {len(extra_in_gg2):>6} (例 {list(extra_in_gg2)[:3]})\n"
            f"  silva - otu: {len(extra_in_silva):>6} (例 {list(extra_in_silva)[:3]})"
        )

    var_proj = pd.DataFrame({'project': proj}, index=asv_local)
    var_proj = var_proj.join(gg2_ranks.reindex(asv_local))
    var_proj = var_proj.join(silva.reindex(asv_local))
    var_proj.index = [f"{proj}__{a}" for a in asv_local]

    mat = sp.csr_matrix(otu.values)             # ASV × sample
    coo = mat.T.tocoo()                          # sample × ASV
    coo_rows_list.append(coo.row.astype(np.int64) + sample_offset)
    coo_cols_list.append(coo.col.astype(np.int64) + asv_offset)
    coo_data_list.append(coo.data.astype(np.int32, copy=False))

    sample_id_list.extend(sample_local)
    asv_id_list.extend(var_proj.index)
    var_dfs.append(var_proj)

    sample_offset += n_sample
    asv_offset    += n_asv
    total_reads   += int(otu.values.sum())
    total_nnz     += coo.nnz

    if i % 50 == 0 or i == len(projects):
        print(f"  [{i:>3d}/{len(projects)}] {proj:<14}  "
              f"ASV={n_asv:>6d}  sample={n_sample:>5d}  "
              f"累计 sample={sample_offset:,}  ASV={asv_offset:,}  nnz={total_nnz:,}")

# %% [markdown]
# ## Step D: 校验全局 sample/ASV ID 唯一性

# %%
print(f"\n全局 sample 数: {len(sample_id_list):,}")
print(f"全局 ASV    数: {len(asv_id_list):,}")
print(f"全局 reads  数: {total_reads:,}")
print(f"全局 nnz    数: {total_nnz:,}")

dup_samples = [k for k, v in Counter(sample_id_list).items() if v > 1]
dup_asvs    = [k for k, v in Counter(asv_id_list).items()    if v > 1]
print(f"重复 sample ID: {len(dup_samples)}  (示例: {dup_samples[:5]})")
print(f"重复 ASV    ID: {len(dup_asvs)}     (示例: {dup_asvs[:5]})")

assert not dup_asvs, "ASV ID 命名空间化后仍有重复，逻辑出错"
if dup_samples:
    print("WARN: 跨项目存在同名 sample，下游使用时需注意 (例如 reindex 或加项目前缀)")

# %% [markdown]
# ## Step E: 构 sparse CSR + 写 AnnData

# %%
print("拼接 COO 三元组 ...")
rows = np.concatenate(coo_rows_list); del coo_rows_list
cols = np.concatenate(coo_cols_list); del coo_cols_list
data = np.concatenate(coo_data_list); del coo_data_list
print(f"  rows.dtype={rows.dtype}  cols.dtype={cols.dtype}  data.dtype={data.dtype}")

n_samples = len(sample_id_list)
n_asvs    = len(asv_id_list)

print(f"构建 COO ({n_samples:,} × {n_asvs:,}) ...")
X_coo = sp.coo_matrix((data, (rows, cols)), shape=(n_samples, n_asvs))
print("→ CSR ...")
X = X_coo.tocsr()
del X_coo, rows, cols, data
print(f"X: shape={X.shape}, nnz={X.nnz:,}, dtype={X.dtype}")
assert X.nnz == total_nnz, f"nnz 不一致: {X.nnz} vs {total_nnz}"
assert int(X.sum()) == total_reads, f"reads 总数不一致: {int(X.sum())} vs {total_reads}"

# %%
var_df = pd.concat(var_dfs, axis=0, copy=False)
var_df.index.name = 'asv_id'
assert len(var_df) == n_asvs
assert (var_df.index == asv_id_list).all(), "var_df 顺序与 asv_id_list 不一致"

print("var_df.shape =", var_df.shape)
print("var_df.columns =", list(var_df.columns))
var_df.head(3)

# %%
adata = ad.AnnData(
    X   = X,
    obs = pd.DataFrame(index=pd.Index(sample_id_list, name='sample_id')),
    var = var_df,
)
print(adata)

print(f"\n写出: {OUT_PATH}")
adata.write_h5ad(OUT_PATH, compression='gzip')
print("完成")

# %% [markdown]
# ## Step F: 抽样校验
#
# 抽 5 个项目，把 02 输出的 sample × ASV 子矩阵列和与原 `otutab.txt` 列和对比。

# %%
import random
random.seed(0)
check_projs = random.sample(projects, 5)
print(f"抽样校验 5 个项目: {check_projs}")

for proj in check_projs:
    otu = pd.read_csv(DATA_DIR / proj / "results/otutab.txt", sep='\t', index_col=0)
    samples_proj = otu.columns.tolist()
    asvs_proj    = [f"{proj}__{a}" for a in otu.index.tolist()]
    sub = adata[adata.obs.index.isin(samples_proj),
                adata.var.index.isin(asvs_proj)]
    sub_total  = int(np.asarray(sub.X.sum()).item())
    orig_total = int(otu.values.sum())
    status = "OK" if sub_total == orig_total else "MISMATCH"
    print(f"  {proj:<14}  AnnData={sub_total:,}  otutab={orig_total:,}  {status}")

# %% [markdown]
# ## 完成
#
# 后续：03 步用 `parse_gg2` + `make_var_id` (来自 MicrobeAtlas 03) 把
# `gg2_*` 6 级折叠成 var_id，过滤 mito/chloro/非 BA，输出 sample × genus AnnData。
