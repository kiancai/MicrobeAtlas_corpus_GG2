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
# # 11: Finalize → MCFCorpus.gg2.h5ad
#
# 输入：
# - `results/feature_table/merged.gg2.with_phylo.h5ad`  (1,826,126 × 8,114)
#
# 输出：
# - `results/feature_table/MCFCorpus.gg2.h5ad`          (1,826,126 × 8,114)
#
# **核心动作**：在 10 步产物上做 corpus finalize，让下游 MiCoFormer 拿来即用：
# 1. `X` 改成 **relative abundance**（l1-normalize per sample, float32 CSR），原 int32 counts 备份到 `layers['counts']`
# 2. obs 加 3 列：`total_reads`(int64) / `n_taxa`(int32) / `Run_paired_id`(str, NaN if not in cross-DB pair)
# 3. var 加 4 列：`n_samples_observed` / `prevalence` / `mass_fraction` / `mean_rel_abundance_when_present`
# 4. uns 加 `provenance` dict（版本/构建日期/git commit/源文件/语义说明）
#
# **不变**：obs 原 54 列 / var 原 7 列 (6-rank + observed) / varp (taxo_dist, phylo_dist) / shape
#
# **重命名约定**：`MCFCorpus.gg2.h5ad` 是 MiCoFormer 用的最终语料库，phylo 是默认且唯一保留版本。

# %%
from pathlib import Path
from datetime import datetime
import subprocess
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix, diags

ad.settings.allow_write_nullable_strings = True

ROOT    = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
ANN_IN  = ROOT / "results/feature_table/merged.gg2.with_phylo.h5ad"
ANN_OUT = ROOT / "results/feature_table/MCFCorpus.gg2.h5ad"

# %% [markdown]
# ## §1 读入

# %%
adata = ad.read_h5ad(ANN_IN)
print(f"input: {adata.shape}  X dtype={adata.X.dtype}  nnz={adata.X.nnz:,}")
print(f"  obs cols: {len(adata.obs.columns)}  var cols: {len(adata.var.columns)}")
print(f"  varp: {list(adata.varp.keys())}  layers: {list(adata.layers.keys())}")

assert isinstance(adata.X, csr_matrix), f"X 应为 CSR，实际 {type(adata.X).__name__}"
assert np.issubdtype(adata.X.dtype, np.integer), f"X 应为整数 counts，实际 {adata.X.dtype}"

counts = adata.X.copy()

# %% [markdown]
# ## §2 obs 统计：total_reads / n_taxa

# %%
total_reads = np.asarray(counts.sum(axis=1)).flatten().astype(np.int64)
n_taxa      = counts.getnnz(axis=1).astype(np.int32)

print(f"total_reads: min={total_reads.min():,}  median={int(np.median(total_reads)):,}  max={total_reads.max():,}")
print(f"n_taxa     : min={n_taxa.min()}     median={int(np.median(n_taxa))}     max={n_taxa.max()}")

assert (total_reads > 0).all(), f"{(total_reads == 0).sum()} 个样本 total_reads==0（QC 后不应存在）"
assert (n_taxa > 0).all(),      f"{(n_taxa == 0).sum()} 个样本 n_taxa==0"

adata.obs["total_reads"] = total_reads
adata.obs["n_taxa"]      = n_taxa

# %% [markdown]
# ## §3 var 统计：n_samples_observed / prevalence / mass_fraction / mean_when_present
#
# - `n_samples_observed[j]` = genus j 出现的样本数 = `(X[:,j] > 0).sum()` = `getnnz(axis=0)[j]`
# - `prevalence[j]` = `n_samples_observed[j] / n_obs`
# - `mass_fraction[j]` = `X[:,j].sum() / X.sum()`（全语料 reads 占比）
# - `mean_rel_abundance_when_present[j]` = `X_rel[:,j].sum() / n_samples_observed[j]`（出现时均值）

# %%
n_samples_observed = counts.getnnz(axis=0).astype(np.int32)
prevalence         = (n_samples_observed / adata.n_obs).astype(np.float32)

genus_total_reads  = np.asarray(counts.sum(axis=0)).flatten().astype(np.int64)
corpus_total_reads = int(genus_total_reads.sum())
mass_fraction      = (genus_total_reads / max(corpus_total_reads, 1)).astype(np.float32)

print(f"var: {(n_samples_observed > 0).sum():,} / {adata.n_vars:,} genus 至少在 1 个样本出现")
print(f"     未观测 (n=0): {(n_samples_observed == 0).sum():,}（应 ≈ 1,257，对应 observed=False）")
print(f"     prevalence:    min={prevalence.min():.6f}  max={prevalence.max():.4f}  median={np.median(prevalence):.6f}")
print(f"     mass_fraction: top-1 = {mass_fraction.max():.4f}  sum = {mass_fraction.sum():.6f}")
assert abs(mass_fraction.sum() - 1.0) < 1e-3, "mass_fraction 总和不约等于 1"

# %% [markdown]
# ## §4 X → relative abundance（l1-normalize per row），原 counts 备份到 layers

# %%
print("Normalizing X to relative abundance (l1, axis=1) ...")
# 手写 l1-normalize: X_rel = diag(1/row_sums) @ counts (float32)
row_sums_int = np.asarray(counts.sum(axis=1)).flatten().astype(np.float64)
assert (row_sums_int > 0).all(), f"{(row_sums_int == 0).sum()} 个样本 row_sum==0"
inv_row = (1.0 / row_sums_int).astype(np.float32)
X_rel = (diags(inv_row) @ counts.astype(np.float32)).tocsr()
print(f"  X_rel: dtype={X_rel.dtype}  nnz={X_rel.nnz:,}  shape={X_rel.shape}")

# sanity: 每行 sum ≈ 1（QC 已保证 row_sum > 0）
row_sums = np.asarray(X_rel.sum(axis=1)).flatten()
print(f"  row_sum: min={row_sums.min():.6f}  max={row_sums.max():.6f}")
assert ((row_sums > 0.999) & (row_sums < 1.001)).all(), "存在 row_sum 偏离 1 的样本"

# %%
# mean_rel_abundance_when_present 依赖 X_rel
genus_rel_sum = np.asarray(X_rel.sum(axis=0)).flatten()
mean_when_present = np.zeros(adata.n_vars, dtype=np.float32)
nz = n_samples_observed > 0
mean_when_present[nz] = (genus_rel_sum[nz] / n_samples_observed[nz]).astype(np.float32)

adata.var["n_samples_observed"]              = n_samples_observed
adata.var["prevalence"]                      = prevalence
adata.var["mass_fraction"]                   = mass_fraction
adata.var["mean_rel_abundance_when_present"] = mean_when_present

print(f"mean_when_present: min={mean_when_present[nz].min():.2e}  max={mean_when_present[nz].max():.4f}")

# %% [markdown]
# ## §5 obs.Run_paired_id：跨库 Run 重复对的 group id

# %%
run_db_unique = adata.obs.groupby("Run")["Database"].nunique()
paired_runs   = set(run_db_unique[run_db_unique >= 2].index.astype(str))
print(f"cross-DB paired Runs: {len(paired_runs):,}  (期望 32,698)")

run_str = adata.obs["Run"].astype(str)
paired_id = run_str.where(run_str.isin(paired_runs), other=np.nan)
adata.obs["Run_paired_id"] = paired_id.astype("string")

n_paired_rows = adata.obs["Run_paired_id"].notna().sum()
print(f"  paired rows: {n_paired_rows:,}  (期望 ~65,396 = 32,698 × 2)")
assert n_paired_rows >= 2 * len(paired_runs), "paired rows < 2 × paired Runs"

# %% [markdown]
# ## §6 把 X 换成 rel_abundance，原 counts 移到 layers

# %%
adata.X = X_rel
adata.layers["counts"] = counts
print(f"X        : {adata.X.dtype}  nnz={adata.X.nnz:,}")
print(f"layers['counts']: {adata.layers['counts'].dtype}  nnz={adata.layers['counts'].nnz:,}")

# %% [markdown]
# ## §7 uns.provenance

# %%
try:
    git_commit = subprocess.check_output(
        ["git", "-C", str(ROOT), "rev-parse", "HEAD"],
        text=True,
    ).strip()
except Exception as e:
    print(f"⚠️  无法获取 git commit: {e}")
    git_commit = "unknown"

adata.uns["provenance"] = {
    "corpus_name":     "MCFCorpus.gg2",
    "version":         "v1",
    "build_date":      datetime.now().strftime("%Y-%m-%d"),
    "pipeline_commit": git_commit,
    "gg2_version":     "24.09",
    "n_obs":           int(adata.n_obs),
    "n_vars":          int(adata.n_vars),
    "sources": [
        "results/feature_table/gg2.full.qc.with_meta.filtered.h5ad",
        "results/feature_table/resmicrodb.gg2.genus.qc.with_meta.filtered.h5ad",
    ],
    "X_semantics":      "relative_abundance (l1-normalized per sample, float32 CSR, row sums = 1)",
    "counts_semantics": "layers['counts'] = raw int32 CSR (original reads)",
    "varp_keys":        list(adata.varp.keys()),
    "obs_added_cols":   ["total_reads", "n_taxa", "Run_paired_id"],
    "var_added_cols":   ["n_samples_observed", "prevalence", "mass_fraction", "mean_rel_abundance_when_present"],
}
print("provenance:")
for k, v in adata.uns["provenance"].items():
    print(f"  {k}: {v}")

# %% [markdown]
# ## §8 写盘

# %%
print(f"\nWriting {ANN_OUT} ...")
adata.write_h5ad(ANN_OUT, compression="gzip")
size_mb = ANN_OUT.stat().st_size / 1024**2
in_mb   = ANN_IN.stat().st_size / 1024**2
print(f"  input  ({ANN_IN.name}):  {in_mb:.0f} MB")
print(f"  output ({ANN_OUT.name}): {size_mb:.0f} MB")

# %% [markdown]
# ## §9 读回验证

# %%
b = ad.read_h5ad(ANN_OUT, backed="r")
print(f"shape: {b.shape}")
print(f"  X dtype={b.X.dtype}  layers={list(b.layers.keys())}  varp={list(b.varp.keys())}")
print(f"  obs cols ({len(b.obs.columns)}): ... 末三列 = {list(b.obs.columns[-3:])}")
print(f"  var cols ({len(b.var.columns)}): {list(b.var.columns)}")

row0_rel = b.X[0].toarray().flatten() if hasattr(b.X[0], "toarray") else np.asarray(b.X[0]).flatten()
row0_cnt = b.layers["counts"][0].toarray().flatten() if hasattr(b.layers["counts"][0], "toarray") else np.asarray(b.layers["counts"][0]).flatten()
print(f"\n样本 0: X[0] 前 5 个非零 rel_abun = {np.sort(row0_rel[row0_rel > 0])[-5:][::-1]}")
print(f"        counts[0].sum = {int(row0_cnt.sum()):,}  obs.total_reads[0] = {int(b.obs['total_reads'].iloc[0]):,}")
print(f"\n✅ all checks done")
