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
# # MicroBioMap 04: 项目表 + 全局 taxonomy → sample × ASV AnnData  (≡ ResMicroDb 02)
#
# 产出与 RM `resmicrodb.gg2.asv.h5ad` **同一套 var schema**，使下游 05/06/07 能
# 直接复用 RM 03/04/05 原逻辑。**`parse_gg2` 与 RM 02 逐字相同。**
#
# **输入**
# - `rawdata/MicroBioMap/15122187.zip::project_asv_tables.tar.gz` — 483 个
#   `asv/PRJ*.txt`（行=全局 ASV id，列=SRS，值=float count）
# - `rawdata/MicroBioMap/work/taxonomy_gg2.global.tsv` — 03 输出（全局 asv→Taxon/Confidence）
# - `rawdata/MicroBioMap/work/sample_metadata.tsv` — srs↔srr↔project 映射
#
# **输出** `results/feature_table/microbiomap.gg2.asv.h5ad`
# - `X`: sample(SRR) × ASV(全局) CSR int32
# - `obs`: index = SRR（run accession，= corpus 的 Run 键）
# - `var`: index = 全局 ASV id；列 = `[project, gg2_Domain..gg2_Species, gg2_Confidence]`
#
# **与 RM 02 的差异（仅"读输入"层，不改变分类/聚合语义）**
# 1. ASV 为**全局共享 id**（同序列跨项目同号、只分类一次）→ 不做 `PROJECT__ASV`
#    命名空间化；`var['project']` 设常量 `'MicroBioMap'`（下游仅用于诊断，不进丰度）。
# 2. project 表 count 是 **float 文本**（整数值）→ 读入断言整性后 `astype(int32)`。
# 3. 列名是 **SRS**，按 sample_metadata 映射成 **SRR**（corpus 主键）。
# 4. 不带 `silva_*`（MBM 自带 SILVA 来自 DADA2 不同分类器，下游不用；RM 03 容忍 0 个 silva 列）。
#
# **流式策略（磁盘/内存友好）**：zip→tar 单趟流式迭代成员，每个项目表 `pd.read_csv`
# 分块（chunksize=ASV 行）→ 取非零 → 累积 COO。全程不解压到磁盘，峰值内存 ≈ 一个 chunk。

# %%
import zipfile
import tarfile
import io
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

PROJECT_DIR = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
MBM_DIR     = PROJECT_DIR / "rawdata/MicroBioMap"
WORK        = MBM_DIR / "work"
ZIP_PATH    = MBM_DIR / "15122187.zip"
TAX_PATH    = WORK / "taxonomy_gg2.global.tsv"
SAMPLE_META = WORK / "sample_metadata.tsv"
OUT_DIR     = PROJECT_DIR / "results/feature_table"
OUT_PATH    = OUT_DIR / "microbiomap.gg2.asv.h5ad"
OUT_DIR.mkdir(parents=True, exist_ok=True)

TAR_MEMBER = "project_asv_tables.tar.gz"
CHUNK      = 20000          # 每次读入的 ASV 行数（内存上界 ≈ CHUNK × n_sample_in_proj）

# %% [markdown]
# ## Step A: GG2 Taxon 解析助手（与 RM 02 逐字一致）

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
# ## Step B: 读全局 taxonomy → 固定 ASV 列空间

# %%
tax = pd.read_csv(TAX_PATH, sep='\t')
tax = tax.rename(columns={tax.columns[0]: 'asv_id'}).set_index('asv_id')
if str(tax.index[0]).startswith('#'):
    tax = tax.iloc[1:]
tax['Confidence'] = pd.to_numeric(tax['Confidence'], errors='coerce')

asv_ids   = tax.index.to_numpy()
n_asv     = len(asv_ids)
asv_to_col = {a: i for i, a in enumerate(asv_ids)}
print(f"全局 ASV（列空间）: {n_asv:,}")
assert len(asv_to_col) == n_asv, "taxonomy 存在重复 ASV id"

# %% [markdown]
# ## Step C: SRS→SRR / SRS→project 映射 + 全局 sample(SRR) 行空间

# %%
sm = pd.read_csv(SAMPLE_META, sep='\t', usecols=['srs', 'srr', 'project'], dtype=str)
srs2srr  = dict(zip(sm['srs'], sm['srr']))
srs2proj = dict(zip(sm['srs'], sm['project']))

srr_all = sm['srr'].tolist()
assert len(set(srr_all)) == len(srr_all), "sample_metadata 的 srr 不唯一"
srr_to_row = {s: i for i, s in enumerate(srr_all)}
n_sample = len(srr_all)
print(f"全局 sample（行空间，SRR）: {n_sample:,}")

# %% [markdown]
# ## Step D: 单趟流式读 483 项目表 → 累积 COO

# %%
rows_list, cols_list, data_list = [], [], []
total_reads = 0
n_proj = 0
seen_srr_cols = 0          # 所有项目表的样本列总数（应 == n_sample）
proj_names = []

zf  = zipfile.ZipFile(ZIP_PATH)
gz  = zf.open(TAR_MEMBER)
tar = tarfile.open(fileobj=gz, mode='r|gz')

for member in tar:
    if not (member.isfile() and member.name.endswith('.txt')):
        continue
    if Path(member.name).name.startswith('._'):    # 跳过 macOS AppleDouble 垃圾文件
        continue
    proj = Path(member.name).stem
    proj_names.append(proj)
    # 流式 tar 的 extractfile 不支持 seekable() → 先读进 BytesIO（可 seek）再喂 pandas
    buf = io.BytesIO(tar.extractfile(member).read())

    proj_cols_done = None
    proj_reads = 0
    for chunk in pd.read_csv(buf, sep='\t', index_col=0, chunksize=CHUNK):
        if proj_cols_done is None:
            # 列：SRS → SRR → 全局 row idx（每个项目只算一次）
            cols_srs = chunk.columns.to_numpy()
            cols_srr = np.array([srs2srr.get(s, None) for s in cols_srs], dtype=object)
            if any(v is None for v in cols_srr):
                bad = [s for s, v in zip(cols_srs, cols_srr) if v is None][:5]
                raise RuntimeError(f"项目 {proj} 有 SRS 无法映射到 SRR，例: {bad}")
            proj_col_row = np.array([srr_to_row[s] for s in cols_srr], dtype=np.int32)
            proj_cols_done = len(cols_srs)
            seen_srr_cols += len(cols_srs)

        # ASV 行 → 全局列 idx（每 chunk ≤ CHUNK 次查表）
        asv_rows = chunk.index.to_numpy()
        try:
            chunk_asv_col = np.array([asv_to_col[a] for a in asv_rows], dtype=np.int32)
        except KeyError as e:
            raise RuntimeError(f"项目 {proj} 出现不在全局 taxonomy 的 ASV: {e}")

        vals = chunk.to_numpy()                       # float (ASV × sample)
        rr, cc = np.nonzero(vals)
        if rr.size == 0:
            continue
        v = vals[rr, cc]
        # 整性断言（read count 必为整数）
        if not np.all(v == np.rint(v)):
            raise RuntimeError(f"项目 {proj} 出现非整数 count")
        rows_list.append(proj_col_row[cc])            # sample 行
        cols_list.append(chunk_asv_col[rr])           # ASV 列
        data_list.append(v.astype(np.int32))
        proj_reads += int(v.sum())

    total_reads += proj_reads
    n_proj += 1
    if n_proj % 50 == 0 or n_proj <= 3:
        print(f"  [{n_proj:>3d}] {proj:<14} samples={proj_cols_done:<5} "
              f"reads={proj_reads:>12,} 累计 reads={total_reads:,}")

tar.close(); gz.close(); zf.close()
print(f"\n项目数: {n_proj}  |  样本列总数: {seen_srr_cols:,} (sample_metadata: {n_sample:,})")
if seen_srr_cols != n_sample:
    print(f"  WARN: 差 {n_sample - seen_srr_cols:,}（>0=有样本无项目表→全0行；<0=有跨表重复样本，已 sum_duplicates）")

# %% [markdown]
# ## Step E: 拼 COO → CSR + 守恒校验

# %%
rows = np.concatenate(rows_list); del rows_list
cols = np.concatenate(cols_list); del cols_list
data = np.concatenate(data_list); del data_list
print(f"非零三元组: {len(data):,}  (rows {rows.dtype} / cols {cols.dtype} / data {data.dtype})")

X = sp.coo_matrix((data, (rows, cols)), shape=(n_sample, n_asv)).tocsr()
X.sum_duplicates()      # 理论上无 (sample,asv) 重复（每样本仅属一个项目），保险
del rows, cols, data
print(f"X: {X.shape}, nnz={X.nnz:,}, dtype={X.dtype}")
assert int(X.sum()) == total_reads, f"reads 守恒失败: {int(X.sum())} vs {total_reads}"
print(f"✓ reads 守恒: {total_reads:,}")

# 列覆盖：有多少全局 ASV 从未在任何项目表出现（期望 0，因 obs_md=项目表里的所有 ASV）
col_nnz = np.asarray((X != 0).sum(axis=0)).ravel()
n_empty_col = int((col_nnz == 0).sum())
print(f"从未出现的全局 ASV 列: {n_empty_col:,} / {n_asv:,}")

# %% [markdown]
# ## Step F: 构 var（gg2_* + Confidence + project）→ 写 AnnData

# %%
# 与 RM 02 同用 parse_gg2 解析；展开用 DataFrame(list-of-dict) 替代 .apply(pd.Series)
# （结果完全一致，仅为在 366 万行规模上避开 .apply(pd.Series) 的极慢路径）
gg2_ranks = pd.DataFrame([parse_gg2(t) for t in tax['Taxon']], index=tax.index)[RANK_COLS]
gg2_ranks.columns = [f'gg2_{c}' for c in RANK_COLS]
gg2_ranks['gg2_Confidence'] = tax['Confidence'].astype(np.float32)

var_df = pd.DataFrame({'project': 'MicroBioMap'},
                      index=pd.Index(asv_ids, name='asv_id'))
var_df = var_df.join(gg2_ranks.reindex(asv_ids))
assert (var_df.index.to_numpy() == asv_ids).all(), "var 顺序与列空间不一致"
print("var_df.shape =", var_df.shape, "| 列 =", list(var_df.columns))

adata = ad.AnnData(
    X   = X,
    obs = pd.DataFrame(index=pd.Index(srr_all, name='sample_id')),
    var = var_df,
)
print(adata)

print(f"\n写出: {OUT_PATH}")
adata.write_h5ad(OUT_PATH, compression='gzip')
print("完成。下一步: 05_build_genus_anndata.py（RM 03 逐字逻辑）")

# %% [markdown]
# ## Step G: 抽样校验（5 个项目，子矩阵列和 vs 原表总和）

# %%
import io
zf  = zipfile.ZipFile(ZIP_PATH)
gz  = zf.open(TAR_MEMBER)
tar = tarfile.open(fileobj=gz, mode='r|gz')
check = set(np.random.RandomState(0).choice(proj_names, size=5, replace=False))
done = 0
for member in tar:
    if not (member.isfile() and member.name.endswith('.txt')):
        continue
    if Path(member.name).name.startswith('._'):
        continue
    proj = Path(member.name).stem
    if proj not in check:
        continue
    df = pd.read_csv(io.BytesIO(tar.extractfile(member).read()), sep='\t', index_col=0)
    orig_total = int(np.rint(df.to_numpy().sum()))
    srrs = [srs2srr[s] for s in df.columns]
    rows_idx = [srr_to_row[s] for s in srrs]
    sub_total = int(X[rows_idx, :].sum())
    status = "OK" if sub_total == orig_total else "MISMATCH"
    print(f"  {proj:<14} AnnData={sub_total:,}  原表={orig_total:,}  {status}")
    done += 1
    if done == len(check):
        break
tar.close(); gz.close(); zf.close()
