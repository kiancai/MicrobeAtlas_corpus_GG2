#!/usr/bin/env python3
"""PRIME Layer1：gg_genus_absolute.csv → prime.gg2.genus.h5ad（符合组成表契约）。

PRIME 已是 **GG2 24.09** genus counts 表（实测 5,308 深层真属全部 ⊆ 8,114 词表，零版本漂移），
故 Tier1 零重处理，只需 reshape + taxon 清理 + 对齐契约 var/obs schema（对标 MBM/RM 表）。

契约（见 .claude/merge_architecture.md）：
- counts（int32）、**全样本保留**（含零 taxon）、不归一、不筛选
- taxon 清理：去 Mitochondria / Chloroplast / 非 Bacteria&Archaea（PRIME 仅 3 列）
- var = var_id(6级路径) 索引 + [Domain,Phylum,Class,Order,Family,Genus]；浅层 `g__` 桶保留
- obs 主键 = Run（PRIME 的 Sample 列即 run，如 DRR396973）
"""
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

ROOT = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
IN_CSV = ROOT / "rawdata/PRIME/gg_genus_absolute.csv"
OUT    = ROOT / "results/feature_table/prime.gg2.genus.h5ad"
CHUNK  = 4000

LEVELS = ["Domain", "Phylum", "Class", "Order", "Family", "Genus"]

# ---- 1. 读列头 → var_id + 解析六级 + 决定清理掩码 ----
with open(IN_CSV) as f:
    header = f.readline().rstrip("\n").split(",")
var_ids = header[1:]
n_var_raw = len(var_ids)
print(f"原始 var 列 = {n_var_raw:,}")

parsed = []
for p in var_ids:
    segs = p.split(";")
    segs = (segs + ["", "", "", "", "", ""])[:6]   # 补齐到 6 级
    parsed.append(segs)
var_df = pd.DataFrame(parsed, columns=LEVELS)
var_df["var_id"] = var_ids

# taxon 清理掩码：保留 Domain∈{Bacteria,Archaea} 且 不含 Mitochondria/Chloroplast
dom_ok = var_df["Domain"].isin(["d__Bacteria", "d__Archaea"])
has_mito = var_df[LEVELS].apply(lambda c: c.str.contains("Mitochondria", case=False, na=False)).any(axis=1)
has_chlo = var_df[LEVELS].apply(lambda c: c.str.contains("Chloroplast",  case=False, na=False)).any(axis=1)
keep_var = dom_ok & ~has_mito & ~has_chlo
drop_idx = np.where(~keep_var.values)[0]
print(f"清理掉 {len(drop_idx)} 列: {[var_ids[i] for i in drop_idx]}")
keep_pos = np.where(keep_var.values)[0]            # 在原始列中的位置
print(f"保留 var = {len(keep_pos):,}")

# ---- 2. 分块读 counts → 保留列 → 稀疏 ----
print(f"分块读 {IN_CSV.name} (chunk={CHUNK}) ...", flush=True)
obs_names = []
blocks = []
reader = pd.read_csv(IN_CSV, chunksize=CHUNK, dtype={0: str})
for i, chunk in enumerate(reader):
    obs_names.extend(chunk.iloc[:, 0].astype(str).tolist())
    mat = chunk.iloc[:, 1:].to_numpy(dtype=np.int32, copy=False)   # (rows, n_var_raw)
    mat = mat[:, keep_pos]                                         # 列清理
    blocks.append(sp.csr_matrix(mat))
    if (i + 1) % 5 == 0:
        print(f"  已读 {len(obs_names):,} 行", flush=True)
X = sp.vstack(blocks, format="csr")
del blocks
print(f"X shape = {X.shape}  dtype={X.dtype}  nnz={X.nnz:,}")
assert X.shape[0] == len(obs_names)
assert X.shape[1] == len(keep_pos)

# ---- 3. 组装 AnnData（契约 schema）----
var_keep = var_df.iloc[keep_pos].copy().set_index("var_id")
for c in LEVELS:
    var_keep[c] = var_keep[c].astype("category")
obs = pd.DataFrame(index=pd.Index(obs_names, name="sample_id"))

A = ad.AnnData(X=X, obs=obs, var=var_keep)
A.var_names = var_keep.index
print(f"\nAnnData: {A.shape}")

# ---- 4. 契约校验 ----
assert A.X.dtype == np.int32, "X 必须 int32 counts"
assert A.obs_names.is_unique, "obs(Run) 必须唯一"
assert not A.var_names.duplicated().any(), "var_id 必须唯一"
real_genus = (A.var["Genus"].astype(str).str.len() > 3).sum()
shallow = A.n_vars - real_genus
zero_rows = int((np.asarray(X.sum(axis=1)).ravel() == 0).sum())
print(f"  真 genus = {real_genus:,}  浅层桶 = {shallow:,}  零 taxon 样本(保留) = {zero_rows:,}")
# domain 纯度
doms = set(A.var["Domain"].astype(str).unique())
assert doms <= {"d__Bacteria", "d__Archaea"}, f"残留非 B&A domain: {doms}"
print(f"  domain = {doms}")

# ---- 5. 写出 ----
OUT.parent.mkdir(parents=True, exist_ok=True)
A.write_h5ad(OUT, compression="gzip")
print(f"\n✓ 写出 {OUT}  ({A.shape[0]:,} × {A.shape[1]:,})")
