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
# # 02: 分类法层级距离矩阵（GG2 24.09 全量）
#
# 输入：
# - `results/phylogeny/genus_vocab.tsv`     GG2 24.09 全 ~8,114 个 genus 的 6 级 rank 路径
#
# 输出：
# - `results/phylogeny/genus_taxo_dist.npz` int8 矩阵 (~8114×8114) + var_id 索引
#
# **距离定义**：两个 genus 在 GG2 7 级分类法层级上"差几个 rank 级别"
#
# | 取值 | 含义 |
# |---|---|
# | 0 | 同 genus（自己 vs 自己） |
# | 1 | 同 family 不同 genus |
# | 2 | 同 order 不同 family |
# | 3 | 同 class 不同 order |
# | 4 | 同 phylum 不同 class |
# | 5 | 同 domain 不同 phylum |
# | 6 | 跨 domain（Archaea ↔ Bacteria） |
#
# 算法：对每对 (i, j)，从最深 rank (Genus) 向上找第一个相等的 rank。
# 向量化实现：从 Domain 到 Genus 依次比较，每一级 same 都覆盖前一次结果，最后一次写入即最深匹配。
#
# 矩阵尺寸：8114² × 1 byte (int8) ≈ 65 MB 内存，gzip 落盘后估计 ~10–15 MB。

# %%
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
VOCAB_IN = ROOT / "results/phylogeny/genus_vocab.tsv"
OUT = ROOT / "results/phylogeny/genus_taxo_dist.npz"

RANK_COLS = ["Domain", "Phylum", "Class", "Order", "Family", "Genus"]

# %% [markdown]
# ## §1 读 vocab

# %%
vocab = pd.read_csv(VOCAB_IN, sep="\t", index_col="var_id")
print(f"vocab shape: {vocab.shape}")
N = len(vocab)
print(f"genus 数 N = {N}")
vocab.head(3)

# %% [markdown]
# ## §2 距离矩阵向量化计算
#
# 思路：初始化距离 = 6（跨 domain，最远）；
# 从 Domain (r=0) 到 Genus (r=5) 依次比较，
# 凡是在 rank r 上相等的 pair，把距离覆盖为 `5 - r`。
# 由于 7 级分类法是严格嵌套（家庭层相等必然分类层也相等），
# 越深的 rank 写入越晚，最终保留的是最深的匹配。

# %%
ranks = vocab[RANK_COLS].values  # (N, 6) object array

# 各 rank 编码为整数（更快的等值比较）
codes = np.zeros((N, 6), dtype=np.int32)
for r, col in enumerate(RANK_COLS):
    # pd.factorize 给出 (codes, uniques)
    codes[:, r], _ = pd.factorize(vocab[col])
print(f"各 rank 唯一值数: {[codes[:, r].max() + 1 for r in range(6)]}")

# %%
dist = np.full((N, N), 6, dtype=np.int8)

for r in range(6):
    same = codes[:, r][:, None] == codes[:, r][None, :]   # (N, N) bool
    dist[same] = 5 - r
    del same

print(f"dist 矩阵 dtype={dist.dtype} shape={dist.shape}  "
      f"占用 {dist.nbytes / 1024**2:.1f} MB")

# %% [markdown]
# ## §3 sanity check

# %%
# 对角线应该全 0
assert (np.diag(dist) == 0).all(), "对角线非 0"

# 对称性
assert (dist == dist.T).all(), "矩阵非对称"

# 取值集合
unique_vals, counts = np.unique(dist, return_counts=True)
print("距离取值分布:")
for v, c in zip(unique_vals, counts):
    pct = 100 * c / dist.size
    print(f"  d = {v}: {c:>13,}  ({pct:5.2f}%)")

# %% [markdown]
# ### 抽样验证：几个手挑 pair

# %%
def lookup(genus_name):
    idx = vocab.index[vocab["Genus"] == genus_name]
    if len(idx) == 0:
        return None
    return vocab.index.get_loc(idx[0])

# 跨 domain → 应为 6
g1 = vocab.index[vocab["Domain"] == "d__Archaea"][0]
g2 = vocab.index[vocab["Domain"] == "d__Bacteria"][0]
i1, i2 = vocab.index.get_loc(g1), vocab.index.get_loc(g2)
print(f"跨 domain 示例: {vocab.loc[g1, 'Genus']} vs {vocab.loc[g2, 'Genus']}  d = {dist[i1, i2]}")

# 同 family 不同 genus → 应为 1
fams = vocab.groupby("Family").size()
big_fam = fams[fams >= 2].index[0]
two = vocab[vocab["Family"] == big_fam].index[:2]
i1, i2 = vocab.index.get_loc(two[0]), vocab.index.get_loc(two[1])
print(f"同 family ({big_fam}) 示例: {vocab.loc[two[0], 'Genus']} vs {vocab.loc[two[1], 'Genus']}  d = {dist[i1, i2]}")

# %% [markdown]
# ## §4 落盘
#
# 用 `np.savez_compressed` 存矩阵 + var_id 索引（保证后面挂回 anndata 时 var_names 对齐）。

# %%
np.savez_compressed(
    OUT,
    dist=dist,
    var_id=vocab.index.values.astype(str),
)
size_mb = OUT.stat().st_size / 1024**2
print(f"已写出: {OUT}")
print(f"  文件大小: {size_mb:.1f} MB（gzip 压缩后）")
