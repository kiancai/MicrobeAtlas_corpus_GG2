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
# # 09: 合并 MicrobeAtlas + ResMicroDb
#
# 输入：
# - `gg2.full.qc.with_meta.filtered.h5ad`               MA (1,732,701 × 6,306, obs 26 列)
# - `resmicrodb.gg2.genus.qc.with_meta.filtered.h5ad`   RM (   93,425 × 4,952, obs 36 列)
#
# 输出：`merged.gg2.h5ad` (1,826,126 × 6,857, obs 54 列)
#
# **合并设计**：
#
# 1. **obs_names**：全局唯一加前缀编号（`MA_0000000..MA_1732700`, `RM_0000000..RM_0093424`）。
#    Run 字段独立保留，允许重复（实测两库 Run 交集 32,698 个，按"语料库立场"两份都留，
#    视为同样本不同流水线视角）。
# 2. **obs 列对齐**：
#    - 公共组 9 列不前缀：`Database, Run, BioSample, Project_ID, Sequencing_Type, Sex,
#      Smoking, Latitude, Longitude`
#    - MA 独有 17 列保持已有 `MA_*` 前缀
#    - RM 独有 28 列新加 `RM_*` 前缀（产物 07 步无前缀；前缀只在 09 合并步引入）
# 3. **标准化动作**（合并前应用）：
#    - `Sex`: MA `female/male` → `Female/Male`（统一首字母大写，按 RM 标准）
#    - `Smoking`: 取值集 ⊇ {Smoker, Non-smoker, Ex-smoker}（MA 只出现前两个，不动）
#    - `Sequencing_Type`: 不映射，取值集 = {AMPLICON, WGS, 16S}
#    - `BioSample`: SRS（MA） / SAMN（RM）不同 ID 系统并存
# 4. **var 合并**：`ad.concat(join='outer')` → 6,857 = 4,401 共享 + 1,905 MA 独有 + 551 RM 独有

# %%
from pathlib import Path
import anndata as ad
import pandas as pd
import numpy as np

ad.settings.allow_write_nullable_strings = True

ROOT = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
MA_IN = ROOT / "results/feature_table/gg2.full.qc.with_meta.filtered.h5ad"
RM_IN = ROOT / "results/feature_table/resmicrodb.gg2.genus.qc.with_meta.filtered.h5ad"
OUT = ROOT / "results/feature_table/merged.gg2.h5ad"

COMMON_COLS = [
    "Database", "Run", "BioSample", "Project_ID", "Sequencing_Type",
    "Sex", "Smoking", "Latitude", "Longitude",
]

# %% [markdown]
# ## §1 读入两库

# %%
ma = ad.read_h5ad(MA_IN)
rm = ad.read_h5ad(RM_IN)
print(f"MA: {ma.shape}  obs cols: {len(ma.obs.columns)}")
print(f"RM: {rm.shape}  obs cols: {len(rm.obs.columns)}")

print()
print(f"MA var ∩ RM var: {len(set(ma.var_names) & set(rm.var_names)):,}")
print(f"MA var ∪ RM var: {len(set(ma.var_names) | set(rm.var_names)):,}")

# %% [markdown]
# ## §2 obs schema 一致性检查

# %%
ma_cols = set(ma.obs.columns)
rm_cols = set(rm.obs.columns)
common_present = ma_cols & rm_cols
print(f"两库实际同名列: {sorted(common_present)}")
expected_common = set(COMMON_COLS)
missing_in_either = expected_common - common_present
assert not missing_in_either, f"预期公共列缺失: {missing_in_either}"
print(f"预期 9 个公共列全部对齐 ✓")

print(f"\nMA 独有列 ({len(ma_cols - rm_cols)}): {sorted(ma_cols - rm_cols)}")
print(f"\nRM 独有列 ({len(rm_cols - ma_cols)}): {sorted(rm_cols - ma_cols)}")

# %% [markdown]
# ## §3 Sex 标准化（MA: female/male → Female/Male）

# %%
print("MA Sex 标准化前:")
print(ma.obs["Sex"].value_counts(dropna=False).to_string())

# MA Sex 是 category，重命名 categories（不改 dtype，不需重建）
ma.obs["Sex"] = ma.obs["Sex"].cat.rename_categories({"female": "Female", "male": "Male"})

print("\nMA Sex 标准化后:")
print(ma.obs["Sex"].value_counts(dropna=False).to_string())

# RM Sex 已经是 Female/Male
print("\nRM Sex (对照):")
print(rm.obs["Sex"].value_counts(dropna=False).to_string())

# %% [markdown]
# ## §4 RM 独有列加 RM_ 前缀

# %%
rm_only = [c for c in rm.obs.columns if c not in COMMON_COLS]
print(f"RM 待加前缀的列 ({len(rm_only)}): {rm_only}")
rm.obs = rm.obs.rename(columns={c: f"RM_{c}" for c in rm_only})
print(f"\n重命名后 RM obs cols: {list(rm.obs.columns)}")

# %% [markdown]
# ## §5 obs_names 全局重编号

# %%
ma.obs_names = pd.Index([f"MA_{i:07d}" for i in range(ma.n_obs)])
rm.obs_names = pd.Index([f"RM_{i:07d}" for i in range(rm.n_obs)])

print(f"MA obs_names: {ma.obs_names[0]} ... {ma.obs_names[-1]}")
print(f"RM obs_names: {rm.obs_names[0]} ... {rm.obs_names[-1]}")

# %% [markdown]
# ## §6 ad.concat（outer-join var；obs 列各自对齐）

# %%
print("合并中...")
merged = ad.concat(
    [ma, rm],
    axis=0,
    join="outer",          # var 取并集
    merge="unique",        # var attrs 用唯一值合并
    fill_value=0,          # X 缺失补 0（合理：未观测 = 0 计数）
)
print(f"merged: {merged.shape}  obs cols: {len(merged.obs.columns)}")
print(f"\nobs cols 完整列表:")
for c in merged.obs.columns:
    print(f"  {c}")

# %% [markdown]
# ## §7 合并后 sanity

# %%
# Database 分布
print("Database 分布:")
print(merged.obs["Database"].value_counts(dropna=False).to_string())

# obs_names 唯一性
assert merged.obs_names.is_unique, "obs_names 不唯一！"
print("\nobs_names 全局唯一 ✓")

# Sex 取值集
sex_vals = set(merged.obs["Sex"].dropna().unique())
print(f"\nSex 取值集: {sex_vals}")
assert sex_vals <= {"Female", "Male"}, f"Sex 含非标准值: {sex_vals}"

# Run 重复（预期 32,698）
n_run_dup = merged.obs["Run"].duplicated(keep=False).sum()
n_run_uniq = merged.obs["Run"].nunique()
print(f"\nRun 行数: {merged.n_obs:,}  唯一 Run: {n_run_uniq:,}  重复行: {n_run_dup:,}")

# var 完整性
print(f"\nvar 来源:")
ma_vars_in_merged = set(ma.var_names) & set(merged.var_names)
rm_vars_in_merged = set(rm.var_names) & set(merged.var_names)
print(f"  MA var → merged: {len(ma_vars_in_merged):,} (MA 原 {ma.n_vars:,})")
print(f"  RM var → merged: {len(rm_vars_in_merged):,} (RM 原 {rm.n_vars:,})")
print(f"  仅 MA: {len(ma_vars_in_merged - rm_vars_in_merged):,}")
print(f"  仅 RM: {len(rm_vars_in_merged - ma_vars_in_merged):,}")
print(f"  两者皆有: {len(ma_vars_in_merged & rm_vars_in_merged):,}")

# X 类型
print(f"\nX dtype: {merged.X.dtype}  format: {type(merged.X).__name__}")
print(f"X nnz: {merged.X.nnz:,}  density: {merged.X.nnz / (merged.n_obs * merged.n_vars) * 100:.3f}%")

# %% [markdown]
# ## §8 reads 守恒检查

# %%
ma_reads = ma.X.sum()
rm_reads = rm.X.sum()
merged_reads = merged.X.sum()
print(f"MA reads:     {ma_reads:>15,}")
print(f"RM reads:     {rm_reads:>15,}")
print(f"sum:          {ma_reads + rm_reads:>15,}")
print(f"merged reads: {merged_reads:>15,}")
assert int(merged_reads) == int(ma_reads + rm_reads), "reads 不守恒！"
print("reads 守恒 ✓")

# %% [markdown]
# ## §9 写出 + 回读校验

# %%
print(f"写出 {OUT}...")
merged.write_h5ad(OUT, compression="gzip")
print(f"size: {OUT.stat().st_size / 1024**2:.1f} MB")

b = ad.read_h5ad(OUT, backed="r")
print(f"回读: {b.shape}  obs cols: {len(b.obs.columns)}  var cols: {len(b.var.columns)}")
assert b.shape == merged.shape
assert list(b.obs.columns) == list(merged.obs.columns)
print("OK")
