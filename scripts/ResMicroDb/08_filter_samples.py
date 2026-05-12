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
# # ResMicroDb 08: metadata-based 样本筛选（passthrough）
#
# 输入 `resmicrodb.gg2.genus.qc.with_meta.h5ad` (93,425 × 4,952, obs 36 列)
# → `resmicrodb.gg2.genus.qc.with_meta.filtered.h5ad`。
#
# **当前 filter 规则**：**无 row-level 删除**。
#
# - jxt 上游 ASV 流水线已限定 16S，anndata 内 `Sequencing_Type` 100% = `16S`
#   → 与 MA 07 步对齐的"去 RNAseq/NaN"在 RM 这里**已自然满足**
# - 阴性对照（`Sample_Site == 'Negative Control'`, 1,508 个）按"语料库立场"保留，
#   下游训练时再筛
# - 异常坐标不处理（RM 这边坐标本来就没异常值）
#
# 本步保留脚本骨架（与 MA 07 编号对仗，便于 09 merge 直接消费 `.filtered.` 后缀），
# 同时做 sanity checks（Sequencing_Type、obs 完整性）。

# %%
from pathlib import Path
import anndata as ad

ad.settings.allow_write_nullable_strings = True

ROOT = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
IN_H5AD = ROOT / "results/feature_table/resmicrodb.gg2.genus.qc.with_meta.h5ad"
OUT_H5AD = ROOT / "results/feature_table/resmicrodb.gg2.genus.qc.with_meta.filtered.h5ad"

# %% [markdown]
# ## §1 读入 + sanity

# %%
adata = ad.read_h5ad(IN_H5AD)
print(f"in: {adata.shape}  obs cols: {len(adata.obs.columns)}")

# Assert: anndata 全部是 16S（jxt 流水线保证）
seq_unique = set(adata.obs["Sequencing_Type"].dropna().unique())
print(f"Sequencing_Type unique: {seq_unique}")
assert seq_unique == {"16S"}, f"意料外的测序类型: {seq_unique}"

n_seq_nan = adata.obs["Sequencing_Type"].isna().sum()
assert n_seq_nan == 0, f"{n_seq_nan} 行 Sequencing_Type 为空"

# %% [markdown]
# ## §2 应用筛选（当前 = no-op）

# %%
keep = slice(None)  # 全保留
n_drop = 0
print(f"row filter: 删除 {n_drop} 样本（passthrough）")
adata_f = adata[keep].copy()
print(f"after: {adata_f.shape}")

# %% [markdown]
# ## §3 信息性诊断（不 filter，只报告）

# %%
print("=== Sample_Site top 10 ===")
print(adata_f.obs["Sample_Site"].value_counts(dropna=False).head(10).to_string())
print()
n_nc = (adata_f.obs["Sample_Site"] == "Negative Control").sum()
print(f"⚠ Negative Control: {n_nc:,} 个（语料库阶段保留；下游训练时按需筛）")

# %% [markdown]
# ## §4 写出 + 回读校验

# %%
adata_f.write_h5ad(OUT_H5AD, compression="gzip")
print(f"written: {OUT_H5AD}")
print(f"size: {OUT_H5AD.stat().st_size / 1024**2:.1f} MB")

b = ad.read_h5ad(OUT_H5AD, backed="r")
print(f"回读: {b.shape}  obs cols: {len(b.obs.columns)}")
assert b.shape == adata_f.shape
print("OK")
