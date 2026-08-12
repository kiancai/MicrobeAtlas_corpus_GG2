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
# # MicroBioMap 03: 合并 31 片 shard taxonomy → 全局 asv→taxonomy
#
# **输入** `rawdata/MicroBioMap/work/tax/shard_*.taxonomy.tsv`（02 array 分类输出）
# - 每片 3 列：`Feature ID / Taxon / Confidence`（QIIME2 export 格式）
#
# **输出** `rawdata/MicroBioMap/work/taxonomy_gg2.global.tsv`
# - 366 万行，3 列，**与 ResMicroDb 各项目的 `taxonomy_gg2.txt` 同格式**
# - 供 04 步用同一个 `parse_gg2` 解析（与 RM 02 一致）
#
# **校验**：合并行数 == 3,663,370；Feature ID 唯一；== 全局 fasta 的 ASV 集合。
#
# 说明：NB 分类是逐序列独立的，分片只是把同一份全局 ASV 拆开并行，合并后
# 与"一次性整体分类"结果**完全等价**（与 RM 逐项目分类同理，sharding-invariant）。

# %%
from pathlib import Path
import numpy as np
import pandas as pd

PROJECT_DIR = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
WORK        = PROJECT_DIR / "rawdata/MicroBioMap/work"
TAX_DIR     = WORK / "tax"
FASTA       = WORK / "asv_global.fa"
OUT_PATH    = WORK / "taxonomy_gg2.global.tsv"

N_EXPECT = 3_663_370

# %% [markdown]
# ## Step A: 读全部 shard taxonomy 并 concat

# %%
shard_files = sorted(TAX_DIR.glob("shard_*.taxonomy.tsv"))
print(f"找到 {len(shard_files)} 片 taxonomy")
assert len(shard_files) == 31, f"期望 31 片，实际 {len(shard_files)}"

parts = []
for f in shard_files:
    df = pd.read_csv(f, sep='\t', dtype=str)
    # 防御：QIIME 某些版本 export 第二行是 #q2:types 元数据
    if len(df) and str(df.iloc[0, 0]).startswith('#'):
        df = df.iloc[1:]
    parts.append(df)
    print(f"  {f.name:<28} {len(df):>8,} 行  列={list(df.columns)}")

tax = pd.concat(parts, axis=0, ignore_index=True)
tax = tax.rename(columns={tax.columns[0]: 'Feature ID'})
print(f"\n合并后总行数: {len(tax):,}")

# %% [markdown]
# ## Step B: 校验（行数 / 唯一性 / 与 fasta ASV 集一致）

# %%
assert len(tax) == N_EXPECT, f"行数 {len(tax):,} != {N_EXPECT:,}"

dup = tax['Feature ID'].duplicated().sum()
assert dup == 0, f"存在 {dup} 个重复 Feature ID"
print("✓ 行数与唯一性通过")

# 与 fasta ASV 集合一致性
fa_ids = set()
with open(FASTA) as fh:
    for line in fh:
        if line.startswith('>'):
            fa_ids.add(line[1:].strip())
tax_ids = set(tax['Feature ID'])
print(f"fasta ASV: {len(fa_ids):,} | taxonomy ASV: {len(tax_ids):,}")
assert fa_ids == tax_ids, (
    f"ASV 集合不一致: fasta-only={len(fa_ids - tax_ids)}, tax-only={len(tax_ids - fa_ids)}"
)
print("✓ 与 fasta ASV 集合完全一致")

# %% [markdown]
# ## Step C: Confidence 数值化分布 + 写出

# %%
conf = pd.to_numeric(tax['Confidence'], errors='coerce')
print("Confidence 分布:")
for q in [0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1.0]:
    print(f"  {int(q*100):>3d}% 分位: {np.nanquantile(conf, q):.4f}")
print(f"  NaN: {conf.isna().sum():,}")

# genus 级注释率（g__ 后跟真实名；注意 GG2 输出分隔是 "; " 带空格）
has_genus = tax['Taxon'].str.contains(r'g__[A-Za-z0-9]', na=False)
print(f"\n注释到 genus 的 ASV: {has_genus.sum():,} ({has_genus.mean()*100:.2f}%)")

# %%
tax[['Feature ID', 'Taxon', 'Confidence']].to_csv(OUT_PATH, sep='\t', index=False)
print(f"\n已写出: {OUT_PATH}  ({len(tax):,} 行)")
print("下一步: 04_build_asv_anndata.py（项目表 + 此 taxonomy → sample × ASV）")
