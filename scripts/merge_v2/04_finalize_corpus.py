#!/usr/bin/env python3
"""Layer2 步4：finalize → 最终语料 MCFCorpusV3。

- X = 相对丰度（默认 genus/deep_sum，已解析属内部组成）；layers['counts'] = 整数 counts
- obs = merged.obs.parquet（旗标 + 公共/前缀 metadata + GMrepo/mBodyMap 富集）
- var = 8,114（Domain..Genus；删旧 MCFCorpusV2 统计列，重算 prevalence/n_samples_observed）
- varp = taxo_dist + phylo_dist；varm = position_encoding —— **从 MCFCorpusV2 直接拷贝**（var 顺序一致）
- uns = provenance + 相对丰度口径说明（unresolved_reads 支持其它分母）
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, h5py
ad.settings.allow_write_nullable_strings = True
ROOT = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
FT = ROOT / "results/feature_table"
MERGED = FT / "merged.counts.h5ad"
OBSP   = FT / "metadata/merged.obs.parquet"
V2     = FT / "MCFCorpusV2.gg2.h5ad"
OUT    = FT / "MCFCorpusV3.gg2.h5ad"
LEVELS = ["Domain","Phylum","Class","Order","Family","Genus"]

print("读 merged.counts ...", flush=True)
A = ad.read_h5ad(MERGED)
counts = A.X.tocsr().astype(np.int32)
n_obs, n_var = A.shape
print(f"  {A.shape}  nnz={counts.nnz:,}")

# ---- 相对丰度（按 deep 行和归一；零行保持全零）----
row_sum = np.asarray(counts.sum(axis=1)).ravel().astype(np.float64)
inv = np.divide(1.0, row_sum, out=np.zeros_like(row_sum), where=row_sum>0)
X_rel = (sp.diags(inv) @ counts).astype(np.float32)
print(f"  相对丰度 OK；零行(deep_sum=0)={int((row_sum==0).sum()):,}")

# ---- obs ----
print("读 merged.obs.parquet ...", flush=True)
obs = pd.read_parquet(OBSP).reindex(A.obs_names)
assert len(obs)==n_obs and (obs.index==A.obs_names).all()
for c in obs.columns:               # object→string（含 pd.NA），数值列保持
    if obs[c].dtype=="object":
        obs[c]=obs[c].astype("string")
print(f"  obs {obs.shape}")

# ---- var：删旧统计列，保 Domain..Genus + 重算 prevalence ----
var = A.var[[c for c in LEVELS if c in A.var.columns]].copy()
nnz_per_var = counts.getnnz(axis=0)
var["n_samples_observed"] = nnz_per_var.astype(np.int64)
var["prevalence"] = (nnz_per_var / n_obs).astype(np.float32)
print(f"  var {var.shape}  非零属(prevalence>0)={int((nnz_per_var>0).sum()):,}/{n_var}")

# ---- 组装 ----
Af = ad.AnnData(X=X_rel, obs=obs, var=var)
Af.var_names = A.var_names
Af.layers["counts"] = counts

# ---- varp / varm 从 V2 拷贝 ----
print("拷贝 varp/varm（V2）...", flush=True)
with h5py.File(V2,"r") as f:
    Af.varp["taxo_dist"]  = f["varp/taxo_dist"][:]
    Af.varp["phylo_dist"] = f["varp/phylo_dist"][:]
    Af.varm["position_encoding"] = f["varm/position_encoding"][:]
print(f"  varp: taxo_dist{Af.varp['taxo_dist'].shape} phylo_dist{Af.varp['phylo_dist'].shape}  varm: pe{Af.varm['position_encoding'].shape}")

# ---- uns ----
Af.uns["build"] = "merge_v2 (PRIME L1 → 01_concat → 02_merge_metadata → 04_finalize)"
Af.uns["databases"] = "MA+RM+MBM+PRIME+Qiita（全留不删，QC=旗标）"
Af.uns["X_definition"] = ("X = 相对丰度 = counts / deep_genus_sum（已解析属内部组成）。"
    "layers['counts']=整数 counts；obs['unresolved_reads']=浅层桶 reads。"
    "下游可选分母：(a) genus/genus.sum [=X] (b) genus/(genus.sum+unresolved_reads) (c) [genus,unresolved]/total")
Af.uns["qc_def"] = "obs['qc_pass'] = n_deep_reads>=1000 & n_deep_features>=5（仅旗标，一行未删）"

print(f"\n最终: {Af.shape}")
print(obs["Database"].value_counts().to_string())
print(f"写出 {OUT} ...", flush=True)
Af.write_h5ad(OUT, compression="gzip")
print(f"✓ 完成 {OUT}  ({Af.n_obs:,} × {Af.n_vars:,})")
