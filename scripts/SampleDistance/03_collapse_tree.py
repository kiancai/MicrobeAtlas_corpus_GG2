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
# # 03: 把 GG2 23M-tip 树折叠到 8,114 g__ 节点为叶子的 genus tree
#
# 输入：
# - `tools/qiime2_database/greengenes2/tree/2024.09.phylogeny.id.nwk`  GG2 系统发育树
# - `results/phylogeny/genus_vocab.tsv`                                  8,114 g__ vocabulary
# - `results/phylogeny/genus_phylo_dist.npz`                             cross-check 用 patristic
#
# 输出：
# - `results/sample_distance/genus_tree.nwk`  叶子 = 8,114 个 g__ 节点（标签 = g__Xxx）
#
# ## 算法
#
# 1. 复用 `Phylogeny/03_compute_phylogenetic_distance.py` 的手写 newick 解析器：
#    只构建内部节点（~100 万），跳过所有 tip（~23M）。
# 2. 找出全部 g__ 内部节点（应为 8,114 个）。
# 3. 自底向上计算每个内部节点 `has_g_desc`（自身或子代里至少有一个 g__ 节点）。
# 4. 从根递归输出 newick：
#    - g__ 节点 → 作为叶子 emit `g__Xxx:length`
#    - 非 g__ 内部节点：只递归 `has_g_desc=True` 的孩子；若孩子全无 g__ desc 则丢弃
# 5. 交叉验证：从新树重算 20 个随机 pair 的 patristic，与 `genus_phylo_dist.npz` 对比，
#    误差必须 < 1e-3 才允许落盘。

# %%
from pathlib import Path
import re
import sys
import time
import numpy as np
import pandas as pd

ROOT = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
TREE_IN  = Path("/hpcdisk1/limk_group/caiqy/tools/qiime2_database/greengenes2/tree/2024.09.phylogeny.id.nwk")
VOCAB_IN = ROOT / "results/phylogeny/genus_vocab.tsv"
PHYLO_IN = ROOT / "results/phylogeny/genus_phylo_dist.npz"
OUT      = ROOT / "results/sample_distance/genus_tree.nwk"

sys.setrecursionlimit(2_000_000)

# %% [markdown]
# ## §1 读 vocab 和 phylo_dist（用于 cross-check）

# %%
vocab = pd.read_csv(VOCAB_IN, sep="\t", index_col="var_id")
target_genera = list(vocab["Genus"].astype(str))
target_set = set(target_genera)
print(f"vocab 8,114 g__: {len(target_genera):,}")

phylo_npz = np.load(PHYLO_IN)
phylo_ref = phylo_npz["dist"]           # 8114² float32
phylo_ids = phylo_npz["var_id"].astype(str)
print(f"phylo_dist ref shape: {phylo_ref.shape}")

# %% [markdown]
# ## §2 手写 newick 解析（只内部节点）—— 同 Phylogeny/03

# %%
print(f"\n读取 {TREE_IN.name} ({TREE_IN.stat().st_size / 1024**2:.0f} MB) ...")
t0 = time.time()
with open(TREE_IN) as f:
    text = f.read()
print(f"  载入字符数: {len(text):,}  耗时 {time.time() - t0:.1f}s")

# %%
print("解析 newick ...")
t0 = time.time()
nodes = []     # 内部节点列表
stack = []
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
        cid = stack.pop()
        nodes[cid]["label"] = label
        nodes[cid]["length"] = length
    elif c == ",":
        pos += 1
    elif c == ";":
        break
    else:
        # tip 块：跳过 label 和可选 :length
        if c == "'":
            end = text.find("'", pos + 1)
            pos = end + 1
        else:
            end = pos
            while end < n and text[end] not in ",():;":
                end += 1
            pos = end
        if pos < n and text[pos] == ":":
            pos += 1
            end = pos
            while end < n and text[end] not in ",():;":
                end += 1
            pos = end

text = None  # free
print(f"  完成: {len(nodes):,} 个内部节点  耗时 {time.time() - t0:.1f}s")

# 负枝长 clamp（同 Phylogeny/03）
neg_count = sum(1 for nd in nodes if nd["length"] < 0)
if neg_count > 0:
    for nd in nodes:
        if nd["length"] < 0:
            nd["length"] = 0.0
    print(f"  ⚠️  {neg_count} 个负枝长已 clamp 到 0")

M = len(nodes)

# %% [markdown]
# ## §3 找 g__ 节点 + 构造 children list

# %%
print("\n找 g__ 节点 + 建子节点表 ...")
t0 = time.time()
g_re = re.compile(r"g__[A-Za-z0-9_\-\.]+")
is_genus = np.zeros(M, dtype=bool)
g_label = [""] * M
for i, nd in enumerate(nodes):
    if "g__" in nd["label"]:
        m = g_re.search(nd["label"])
        if m:
            is_genus[i] = True
            g_label[i] = m.group(0)

found_set = {g for g in g_label if g}
missing = target_set - found_set
extra = found_set - target_set
print(f"  树里 g__ 节点: {is_genus.sum():,}")
print(f"  vocab 缺（必须 0）: {len(missing)}")
print(f"  树里多 vocab 没的: {len(extra)}")
assert len(missing) == 0, f"vocab 有但树里没有: {list(missing)[:5]}"

# children list
parent_ids = np.array([nd["parent_id"] for nd in nodes], dtype=np.int64)
children = [[] for _ in range(M)]
roots = []
for i in range(M):
    p = parent_ids[i]
    if p == -1:
        roots.append(i)
    else:
        children[p].append(i)
assert len(roots) == 1, f"不止一个根：{len(roots)} 个"
ROOT_ID = roots[0]
print(f"  根节点 id: {ROOT_ID}")
print(f"  耗时 {time.time() - t0:.1f}s")

# %% [markdown]
# ## §4 自底向上计算 has_g_desc

# %%
# parent_id < child_id 已保证（栈顺序），所以倒序遍历即可
has_g_desc = is_genus.copy()
for i in range(M - 1, -1, -1):
    if has_g_desc[i] and parent_ids[i] != -1:
        has_g_desc[parent_ids[i]] = True
print(f"\nhas_g_desc=True 节点数: {int(has_g_desc.sum()):,}  (包含 g__ + 它们的祖先)")
print(f"完全可丢弃的内部节点: {int((~has_g_desc).sum()):,}")

# %% [markdown]
# ## §5 迭代式 emit newick（避免深递归）
#
# 用显式栈做后序遍历输出。每个节点入栈两次：第一次 push children，第二次 emit 自己的合并字符串。

# %%
print("\n生成折叠 newick ...")
t0 = time.time()

# 原树内部节点 label 可能含 `;`、空格、复合 rank (`o__Foo; f__Bar`) → 必须引号包起来
# g__ tip label 只含 alnum/_/-/.，不需要引号
_SAFE = re.compile(r"^[A-Za-z0-9_\-\.]*$")
def quote_label(s: str) -> str:
    if s == "" or _SAFE.fullmatch(s):
        return s
    # 内部单引号转义为两个
    return "'" + s.replace("'", "''") + "'"

# 每个节点最终输出的字符串（占内存：8114 × 平均 30 字符 + 内部 join）
# 为了节省内存，我们在子节点都处理完后立刻把它们的字符串拼成 parent 的字符串，并释放子的字符串
out_str = [None] * M

# 处理嵌套 g__：如果一个 g__ 节点还有 g__ 后代，把它拆成
#   synthetic_internal(枝长 = 原 g__ 枝长)
#     ├── g__Xxx tip (length=0)              ← 保留 self
#     └── recurse(原 g__ 节点的 g__ 后代...)   ← 保留嵌套
# 这样 8,114 个 g__ 全部当 tip，patristic 距离不变。

def kids_with_g(i):
    return [c for c in children[i] if has_g_desc[c]]

# 用 "事件" 栈：(node_id, phase)。phase=0 = pre-visit, phase=1 = post-visit
WORK = [(ROOT_ID, 0)]
while WORK:
    i, phase = WORK.pop()
    L = nodes[i]["length"]
    if phase == 0:
        kw = kids_with_g(i)
        if is_genus[i] and len(kw) == 0:
            # 真叶 g__：直接 emit `g__Xxx:length`
            if parent_ids[i] == -1:
                out_str[i] = g_label[i]
            else:
                out_str[i] = f"{g_label[i]}:{L:g}"
        else:
            # 非 g__ 内部 OR 嵌套 g__：先 push 自己 phase=1，再 push 子节点
            WORK.append((i, 1))
            for c in kw:
                WORK.append((c, 0))
    else:
        kw = kids_with_g(i)
        kids_str = [out_str[c] for c in kw]
        for c in kw:
            out_str[c] = None
        if is_genus[i]:
            # 嵌套 g__：在子节点开头插一个 0 枝长 self-tip
            kids_str = [f"{g_label[i]}:0"] + kids_str
            inner = ",".join(kids_str)
            # 合成内部节点不写 label（g__ label 已被 self-tip 拿走）
            if parent_ids[i] == -1:
                out_str[i] = f"({inner})"
            else:
                out_str[i] = f"({inner}):{L:g}"
        else:
            inner = ",".join(kids_str)
            label = quote_label(nodes[i]["label"])
            if parent_ids[i] == -1:
                out_str[i] = f"({inner}){label}"
            else:
                out_str[i] = f"({inner}){label}:{L:g}"

nwk = out_str[ROOT_ID] + ";"
out_str = None  # free
print(f"  newick 字符数: {len(nwk):,}  耗时 {time.time() - t0:.1f}s")

# %% [markdown]
# ## §6 写盘

# %%
with open(OUT, "w") as f:
    f.write(nwk)
print(f"\n已写出 {OUT}")
print(f"  大小: {OUT.stat().st_size / 1024**2:.2f} MB")

# %% [markdown]
# ## §7 Cross-check: 从新树重算 20 对 patristic 距离，对比 phylo_dist
#
# 用 skbio 读新树（8,114 叶子，~几万节点，skbio 完全 handle 得了），
# `node.distance()` 两两算后跟 ref 比对。

# %%
print("\n用 skbio 重读新树做 cross-check ...")
t0 = time.time()
from skbio import TreeNode
tree = TreeNode.read(str(OUT), convert_underscores=False)
# 建 label → node 的字典
tip_map = {tip.name: tip for tip in tree.tips()}
print(f"  叶子数: {len(tip_map):,}  耗时 {time.time() - t0:.1f}s")
assert len(tip_map) == len(target_genera), \
    f"叶子数 {len(tip_map)} != target {len(target_genera)}"
miss = target_set - set(tip_map.keys())
assert len(miss) == 0, f"新树缺 {len(miss)} 个 target g__: {list(miss)[:5]}"

# %%
# phylo_dist.npz 的索引是 var_id (完整 6 级路径)，需要 g__Genus → index 的映射
# vocab 的行序就是 phylo_dist 的行/列序
assert np.array_equal(vocab.index.values.astype(str), phylo_ids), "vocab vs phylo_dist var_id 顺序不一致"
g_to_idx = {g: i for i, g in enumerate(vocab["Genus"].astype(str).values)}
np.random.seed(20260515)
samples = []
while len(samples) < 20:
    i, j = np.random.choice(len(target_genera), 2, replace=False)
    if i == j:
        continue
    samples.append((i, j))

print(f"\n{'pair':<25} {'tree dist':>14} {'phylo ref':>14} {'diff':>14} {'ok':>4}")
print("-" * 75)
all_ok = True
for i, j in samples:
    gi, gj = target_genera[i], target_genera[j]
    d_tree = float(tip_map[gi].distance(tip_map[gj]))
    d_ref  = float(phylo_ref[g_to_idx[gi], g_to_idx[gj]])
    diff = abs(d_tree - d_ref)
    ok = diff < 1e-3
    if not ok: all_ok = False
    print(f"({gi[:12]:<12}, {gj[:12]:<12})  {d_tree:>14.6f} {d_ref:>14.6f} {diff:>14.2e}  {'✓' if ok else '✗':>4}")

if not all_ok:
    raise RuntimeError("新树 patristic 与 phylo_dist 不一致，折叠有 bug，拒绝接受")
print("\n✅ 全部 20 对吻合，folded tree 通过 cross-check")
