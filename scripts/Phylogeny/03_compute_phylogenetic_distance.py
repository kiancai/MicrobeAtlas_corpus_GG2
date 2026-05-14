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
# # 03: 系统发育（patristic）距离矩阵（GG2 24.09 全量）
#
# 输入：
# - `results/phylogeny/genus_vocab.tsv`                                  8,114 个 genus
# - `tools/qiime2_database/greengenes2/tree/2024.09.phylogeny.id.nwk`   GG2 系统发育树
#
# 输出：
# - `results/phylogeny/genus_phylo_dist.npz`   float32 矩阵 (8114×8114) + var_id 索引
#
# ## 为什么重写：旧版的失败
#
# 上一版 03 用 skbio 的 `tree.shear(target_genera)` + `pruned.tip_tip_distances()`。
# 在 GG2 这种"23M tip + ~30k 内部节点 + 31.7% 0 枝长 + 巨型多分叉"的病态结构下，
# skbio 的 unifurcation collapse 时丢失枝长，**50% 的非对角 pair 算出 0**。
#
# 比如 `g__MCBC01` vs `g__Altiarchaeum`（两个不同 order 的 Archaea genus）：
# - skbio shear+tip_tip_distances: 0.0  ✗
# - skbio node.distance() 直接算:   30.854 ✓
# - 手写解析 + 链遍历:               30.854 ✓
#
# ## 新版方案
#
# **完全不用 skbio**。
#
# 1. 手写一个轻量 newick 解析器，**只构建内部节点**（跳过所有 tip）。GG2 树 23M tip 但只
#    ~100 万 internal node（包括根附近未标 label 的结构节点），实际内存占用小很多。
# 2. 把内部节点的 parent-child + 枝长构造成 scipy 稀疏邻接矩阵
# 3. 用 `scipy.sparse.csgraph.shortest_path` (Dijkstra) 算 8,114 个 g__ 节点到所有
#    内部节点的距离，再切片得到 8114×8114 对称矩阵
# 4. **独立交叉验证**：随机抽 10 个 pair 用手写"祖先链遍历"算法重算，必须完全吻合
#    才允许落盘。否则 raise。

# %%
from pathlib import Path
import re
import time
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path

ROOT = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
VOCAB_IN = ROOT / "results/phylogeny/genus_vocab.tsv"
TREE_IN = Path("/hpcdisk1/limk_group/caiqy/tools/qiime2_database/greengenes2/tree/2024.09.phylogeny.id.nwk")
OUT = ROOT / "results/phylogeny/genus_phylo_dist.npz"

# %% [markdown]
# ## §1 读 vocab

# %%
vocab = pd.read_csv(VOCAB_IN, sep="\t", index_col="var_id")
target_genera = list(vocab["Genus"])
print(f"待处理 genus 数: {len(target_genera)}")
assert len(target_genera) == len(set(target_genera)), "vocab 里 Genus 列重复"

# %% [markdown]
# ## §2 手写 newick 解析器（只保留内部节点）
#
# Newick 语法状态机：
# - `(`  → 开始新内部节点（push 到 stack，等 `)` 时再填 label/length）
# - `)`  → 关闭当前内部节点，紧跟读 label（可能是 `'...'` 引号包裹）和 `:length`
# - `,`  → 兄弟分隔
# - `;`  → 结束
# - 其它：tip 名字，整个 tip 块（label + 可选 `:length`）**跳过不存**
#
# 内部节点结构：`{'label': str, 'parent_id': int, 'length': float}`。
# `parent_id == -1` 表示根。`length` = 该节点到 parent 的枝长。

# %%
print(f"读取 {TREE_IN.name} ({TREE_IN.stat().st_size / 1024**2:.0f} MB) ...")
t0 = time.time()
with open(TREE_IN) as f:
    text = f.read()
print(f"  载入字符数: {len(text):,}  耗时 {time.time() - t0:.1f}s")

# %%
print("解析 newick...")
t0 = time.time()

nodes = []     # 内部节点列表
stack = []     # 当前 open 的内部节点 id 栈
pos = 0
n = len(text)

while pos < n:
    c = text[pos]
    if c == "(":
        new_id = len(nodes)
        parent_id = stack[-1] if stack else -1
        nodes.append({"label": "", "parent_id": parent_id, "length": 0.0})
        stack.append(new_id)
        pos += 1
    elif c == ")":
        pos += 1
        # 读 label
        if pos < n and text[pos] == "'":
            end = text.find("'", pos + 1)
            label = text[pos + 1:end]
            pos = end + 1
        else:
            end = pos
            while end < n and text[end] not in ",():;":
                end += 1
            label = text[pos:end].strip()
            pos = end
        # 读 length
        length = 0.0
        if pos < n and text[pos] == ":":
            pos += 1
            end = pos
            while end < n and text[end] not in ",():;":
                end += 1
            try:
                length = float(text[pos:end])
            except ValueError:
                pass
            pos = end
        # 弹栈，回填
        cid = stack.pop()
        nodes[cid]["label"] = label
        nodes[cid]["length"] = length
    elif c == ",":
        pos += 1
    elif c == ";":
        break
    else:
        # tip 块：跳过 label
        if c == "'":
            end = text.find("'", pos + 1)
            pos = end + 1
        else:
            end = pos
            while end < n and text[end] not in ",():;":
                end += 1
            pos = end
        # 跳过 tip 的 :length
        if pos < n and text[pos] == ":":
            pos += 1
            end = pos
            while end < n and text[end] not in ",():;":
                end += 1
            pos = end

# 释放大字符串内存
text = None

print(f"  完成: {len(nodes):,} 个内部节点  耗时 {time.time() - t0:.1f}s")

# 基本 sanity：负枝长 clamp 到 0（GG2 有 4 个 -1e-6 的浮点噪音，clamp 不影响下游）
neg_count = sum(1 for n in nodes if n["length"] < 0)
if neg_count > 0:
    neg_lens = [n["length"] for n in nodes if n["length"] < 0]
    print(f"  ⚠️  {neg_count} 个内部节点出现负枝长（min={min(neg_lens):.2e}），clamp 到 0")
    for n in nodes:
        if n["length"] < 0:
            n["length"] = 0.0
parent_ids = np.array([n["parent_id"] for n in nodes])
assert (parent_ids[1:] < np.arange(1, len(nodes))).all(), "parent_id 不满足 parent_id < child_id（解析栈顺序应保证这一点）"
print(f"  parent_id < child_id 检查通过")
print(f"  根节点数 (parent_id == -1): {(parent_ids == -1).sum()}")

# %% [markdown]
# ## §3 建 `g__Foo → 节点 id` 映射 + vocab cross-check

# %%
g_re = re.compile(r"g__[A-Za-z0-9_\-\.]+")
g_to_id = {}
dup = []
for i, nd in enumerate(nodes):
    if "g__" not in nd["label"]:
        continue
    m = g_re.search(nd["label"])
    if m:
        g = m.group(0)
        if g in g_to_id:
            dup.append(g)
        g_to_id[g] = i

print(f"找到 g__ 节点: {len(g_to_id):,}")
assert len(dup) == 0, f"⚠️  发现重复 g__ token（GG2 单系应保证唯一）: {dup[:5]}"

missing = [g for g in target_genera if g not in g_to_id]
extra = [g for g in g_to_id if g not in set(target_genera)]
print(f"vocab 缺（必须为 0）: {len(missing)}")
print(f"树里多 vocab 没的:    {len(extra)}")
assert len(missing) == 0, f"vocab 有但树里没有: {missing[:5]}"

# %% [markdown]
# ## §4 构造 sparse 邻接矩阵
#
# 每条 parent-child 边记两次（无向图）：`(child, parent, length)` 和 `(parent, child, length)`。
# 后面 scipy.csgraph 当无向图处理。

# %%
print("构造邻接矩阵...")
t0 = time.time()
M = len(nodes)
# 预估边数 ≈ M-1 条树边（少一个根没有 parent），每条记两次
edge_cnt = 2 * sum(1 for n in nodes if n["parent_id"] != -1)
row = np.empty(edge_cnt, dtype=np.int32)
col = np.empty(edge_cnt, dtype=np.int32)
data = np.empty(edge_cnt, dtype=np.float64)
k = 0
for i, nd in enumerate(nodes):
    p = nd["parent_id"]
    if p == -1:
        continue
    L = nd["length"]
    row[k] = i; col[k] = p; data[k] = L; k += 1
    row[k] = p; col[k] = i; data[k] = L; k += 1
assert k == edge_cnt
adj = csr_matrix((data, (row, col)), shape=(M, M))
print(f"  完成  shape={adj.shape}  nnz={adj.nnz:,}  耗时 {time.time() - t0:.1f}s")
del row, col, data

# %% [markdown]
# ## §5 批量 scipy.shortest_path
#
# 一次性算 `(8114, M)` 的距离会吃 ~65 GB（M ≈ 1M）。分批：BATCH=500，每批 4 GB 中间量。
# 用 Dijkstra（method='D'），无向图。

# %%
target_indices = np.array([g_to_id[g] for g in target_genera], dtype=np.int32)
N = len(target_indices)
dist = np.zeros((N, N), dtype=np.float32)

BATCH = 500
print(f"批量计算 shortest_path: N={N}, BATCH={BATCH}, 总批次={ (N + BATCH - 1) // BATCH }")
t_total = time.time()
for start in range(0, N, BATCH):
    t0 = time.time()
    end = min(start + BATCH, N)
    batch = target_indices[start:end]
    bdist = shortest_path(adj, directed=False, indices=batch, method="D")
    dist[start:end, :] = bdist[:, target_indices].astype(np.float32)
    print(f"  batch [{start:>5}:{end:>5}]  {time.time() - t0:.1f}s")
print(f"总耗时 {time.time() - t_total:.1f}s")

# %% [markdown]
# ## §6 自动 sanity

# %%
print("Sanity checks:")
print(f"  shape:       {dist.shape}")
print(f"  dtype:       {dist.dtype}")
print(f"  对角 max:    {np.diag(dist).max()}  (应为 0)")
print(f"  非负:        {(dist >= 0).all()}")
print(f"  NaN/Inf:     {np.isnan(dist).sum()} / {np.isinf(dist).sum()}")
diff = dist - dist.T
print(f"  对称偏差:    {np.abs(diff).max():.6e}")
assert (np.diag(dist) == 0).all()
assert (dist >= 0).all()
assert np.isnan(dist).sum() == 0 and np.isinf(dist).sum() == 0
assert np.allclose(dist, dist.T, atol=1e-4)

upper = dist[np.triu_indices_from(dist, k=1)]
print(f"  非自身距离 min={upper.min():.6f}  median={np.median(upper):.4f}  max={upper.max():.4f}")
print(f"  0 距离非对角 pair: {(upper == 0).sum():,}  (应该很少，<<{upper.size:,})")

# %% [markdown]
# ## §7 独立交叉验证（chain walking 重算 10 个随机 pair）
#
# **必须全部吻合才落盘**。任一不一致就 raise。

# %%
def chain(i):
    """从节点 i 走到根，返回 (ancestor_id, cumulative_distance_from_i) 列表。"""
    out = [(i, 0.0)]
    cum = 0.0
    while nodes[i]["parent_id"] != -1:
        cum += nodes[i]["length"]
        i = nodes[i]["parent_id"]
        out.append((i, cum))
    return out

def dist_via_chain(gA, gB):
    iA, iB = g_to_id[gA], g_to_id[gB]
    cA = chain(iA)
    cB = chain(iB)
    a_ancestors = {ai: d for ai, d in cA}
    for bi, db in cB:
        if bi in a_ancestors:
            return a_ancestors[bi] + db
    return None

np.random.seed(20260514)
sample_idx = np.random.choice(N, size=(20, 2), replace=True)
print(f"\n{'pair':<18} {'scipy':>14} {'chain':>14} {'diff':>14} {'ok':>4}")
print("-" * 70)
all_ok = True
for i, j in sample_idx:
    if i == j: continue
    gA, gB = target_genera[i], target_genera[j]
    d_scipy = float(dist[i, j])
    d_chain = float(dist_via_chain(gA, gB))
    diff = abs(d_scipy - d_chain)
    ok = diff < 1e-3
    if not ok: all_ok = False
    print(f"({i:>5}, {j:>5})   {d_scipy:>14.6f} {d_chain:>14.6f} {diff:>14.2e}  {'✓' if ok else '✗':>4}")

if not all_ok:
    raise RuntimeError("scipy 与 chain walking 出现不一致，距离矩阵不可信，拒绝落盘")
print("\n✅ 所有抽样 pair 两种方法完全吻合")

# %% [markdown]
# ## §8 phylo dist 跟 taxo hop 的语义一致性
#
# 期望：分类法距离越大（同 family 1 → 跨 domain 6），phylo 距离平均也越大。

# %%
taxo = np.load(ROOT / "results/phylogeny/genus_taxo_dist.npz")["dist"]
print(f"按 taxo hop 分桶（应单调递增）:")
print(f"  {'taxo_hop':>10} {'pair count':>14} {'phylo mean':>12} {'phylo median':>14}")
prev_mean = -1
for hop in [1, 2, 3, 4, 5, 6]:
    mask = np.triu(taxo == hop, k=1)
    if mask.sum() == 0:
        continue
    vals = dist[mask]
    m = vals.mean()
    print(f"  {hop:>10} {int(mask.sum()):>14,} {m:>12.4f} {np.median(vals):>14.4f}")
    assert m > prev_mean, f"分桶均值非单调（hop {hop} mean={m} <= prev {prev_mean}），疑似仍有 bug"
    prev_mean = m
print("✅ 单调性 OK")

# %% [markdown]
# ## §9 落盘

# %%
np.savez_compressed(
    OUT,
    dist=dist,
    var_id=vocab.index.values.astype(str),
)
print(f"已写出: {OUT}")
print(f"  文件大小: {OUT.stat().st_size / 1024**2:.1f} MB")
