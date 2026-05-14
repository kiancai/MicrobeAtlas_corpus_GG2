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
# # 01: 从 GG2 24.09 系统发育树提取全 genus 词表
#
# 输入：
# - `tools/qiime2_database/greengenes2/tree/2024.09.phylogeny.id.nwk`   GG2 系统发育树
#
# 输出：
# - `results/phylogeny/genus_vocab.tsv`     GG2 24.09 全 ~8,114 个 genus 的 6 级 rank 路径
#
# **为什么不用 `taxonomy.id.tsv.gz`**：
# TSV 是 NB 在 tip 级别给出的分类结果，~9M tip 因置信度不够只到 domain/phylum/...，
# 没分到 genus 级别，TSV 里只能"看到" 6,757 个 genus。
# 但 GG2 24.09 实际定义的完整 genus 空间是 8,114 个——其余 1,357 个只出现在
# 树的内部节点 label 里（monophyletic clade 客观存在但 NB 在 tip 级别没把任何 tip
# 分到那里）。
# **NB classifier 是可以输出全部 8,114 个 genus 标签的**，所以下游 anndata 的
# var.Genus 可能包含 TSV 缺的那 1,357 个里的某些——必须按树来取词表。
#
# **核心动作**：
#
# 1. 用 skbio 载入树（~2 min，~10 GB RAM）
# 2. 遍历所有内部节点，从 label 抽 `g__Foo` token
# 3. 对每个 g__Foo 的节点，**walk up 到根**，从祖先 label 里收集 d/p/c/o/f rank token
#    拼成完整 6 级 path（GG2 monophyly decoration 保证每个 rank 只在唯一一条祖先链上出现）
# 4. **单系性 sanity check**：每个 g__ token 应该只在树里出现 1 次（已实测过，但还是查一遍）
# 5. 落盘 vocab.tsv，index = `var_id` = `;`-连接的 6 级 path 字符串
#
# 本步骤**不依赖 anndata**，输出是"GG2 24.09 全量 genus 副产物"，可独立发布。

# %%
from pathlib import Path
import re
import time
import pandas as pd
from skbio import TreeNode

ROOT = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
TREE_IN = Path("/hpcdisk1/limk_group/caiqy/tools/qiime2_database/greengenes2/tree/2024.09.phylogeny.id.nwk")
OUT = ROOT / "results/phylogeny/genus_vocab.tsv"
OUT.parent.mkdir(parents=True, exist_ok=True)

RANK_PREFIXES = ["d__", "p__", "c__", "o__", "f__", "g__"]
RANK_COLS = ["Domain", "Phylum", "Class", "Order", "Family", "Genus"]

# %% [markdown]
# ## §1 载入树

# %%
t0 = time.time()
print(f"载入 {TREE_IN.name} ({TREE_IN.stat().st_size / 1024**2:.0f} MB) ...")
tree = TreeNode.read(str(TREE_IN), format="newick")
print(f"  耗时 {time.time() - t0:.1f}s")
print(f"  tip 数:    {tree.count(tips=True):,}")
print(f"  节点总数:  {tree.count():,}")

# %% [markdown]
# ## §2 抽 rank token 的辅助
#
# 内部节点 label 形如 `'o__SURF-38; f__SURF-38; g__SURF-38; s__SURF-38 sp003599335'`，
# 多 rank 用 `; ` 分隔。同一节点可能同时是 order/family/genus/species 的代表（多个 rank
# 在同一节点合并的情况是 GG2 monophyly decoration 的常见现象——当一个 clade 在多个 rank
# 上都只有 1 个子簇时，这些 rank 都标在同一个内部节点上）。

# %%
TOKEN_RE = {p: re.compile(rf"{re.escape(p)}[A-Za-z0-9_\-\.]*") for p in RANK_PREFIXES}

def extract_ranks(label):
    """从一个 label 字符串里抽出 d/p/c/o/f/g 的 token（缺的返回 ''）。"""
    out = {}
    if not label:
        return out
    for prefix, rex in TOKEN_RE.items():
        m = rex.search(label)
        if m:
            out[prefix] = m.group(0)
    return out


# 简单单元自测
_t = extract_ranks("'o__SURF-38; f__SURF-38; g__SURF-38; s__SURF-38 sp003599335'")
print("test extract_ranks:", _t)
assert _t.get("g__") == "g__SURF-38"
assert _t.get("o__") == "o__SURF-38"

# %% [markdown]
# ## §3 一次遍历，找所有 g__ 节点 + walk up 拼完整 path

# %%
t0 = time.time()

genus_to_node = {}      # g__Foo -> TreeNode
genus_paths = {}        # g__Foo -> (Domain, Phylum, Class, Order, Family, Genus)
dup_warning = []

# 第一遍：找 g__ 节点
for node in tree.non_tips():
    name = node.name
    if not name or "g__" not in name:
        continue
    g_match = TOKEN_RE["g__"].search(name)
    if not g_match:
        continue
    g = g_match.group(0)
    if g in genus_to_node:
        dup_warning.append(g)
    else:
        genus_to_node[g] = node

print(f"第一遍：找到 {len(genus_to_node):,} 个 g__ 节点  耗时 {time.time() - t0:.1f}s")
if dup_warning:
    print(f"  ⚠️  {len(dup_warning)} 个 g__ 在多个节点出现（前 5）: {dup_warning[:5]}")
else:
    print(f"  ✅ 所有 g__ token 都唯一出现（GG2 单系性符合预期）")

# %%
# 第二遍：对每个 g__ 节点 walk up 拼 d/p/c/o/f
t0 = time.time()
for g, node in genus_to_node.items():
    ranks = {"g__": g}
    # 包括节点本身的 label（可能 carry f__ / o__ 等）
    cur = node
    while cur is not None and not all(p in ranks for p in RANK_PREFIXES):
        more = extract_ranks(cur.name)
        for p, v in more.items():
            if p not in ranks:
                ranks[p] = v
        cur = cur.parent
    genus_paths[g] = tuple(ranks.get(p, "") for p in RANK_PREFIXES)

print(f"第二遍：拼 6 级 path  耗时 {time.time() - t0:.1f}s")

# %% [markdown]
# ## §4 缺级 sanity check
#
# 每个 g__Foo 都应该有完整的 d/p/c/o/f 链。如果某 rank 缺失（空字符串），
# 说明 walk up 没找到对应的 ancestor label——这通常是 GG2 树根附近的特殊节点
# （或我们的 regex/数据假设有问题），需要看一下。

# %%
n_incomplete = 0
incomplete_examples = []
for g, path in genus_paths.items():
    if any(v == "" for v in path):
        n_incomplete += 1
        if len(incomplete_examples) < 5:
            incomplete_examples.append((g, path))

print(f"path 完整 (6 级齐全): {len(genus_paths) - n_incomplete:,}")
print(f"path 不完整:          {n_incomplete:,}")
if incomplete_examples:
    print("\n不完整示例:")
    for g, p in incomplete_examples:
        print(f"  {g}: {p}")

# %% [markdown]
# ## §5 构造 vocab DataFrame

# %%
records = []
for gn, path in genus_paths.items():
    d, p, c, o, fam, gx = path
    var_id = ";".join(path)
    records.append({
        "var_id": var_id,
        "Domain": d, "Phylum": p, "Class": c,
        "Order": o, "Family": fam, "Genus": gn,
    })

vocab = pd.DataFrame.from_records(records).set_index("var_id")
vocab = vocab.sort_index()    # 字母序，保证可重现
print(f"vocab shape: {vocab.shape}")

print(f"\nDomain 分布:\n{vocab['Domain'].value_counts()}")
print(f"\n前 5 行:")
print(vocab.head())
print(f"\n后 5 行:")
print(vocab.tail())

# 带 _A/_B/_C 后缀（GTDB monophyly-split）数量
n_suffix = vocab["Genus"].str.contains(r"_[A-Z]$", regex=True).sum()
print(f"\n带 _A/_B/_C 后缀: {n_suffix} / {len(vocab)}")

# %% [markdown]
# ## §6 落盘

# %%
vocab.to_csv(OUT, sep="\t")
print(f"已写出: {OUT}")
print(f"  shape: {vocab.shape}")
print(f"  文件大小: {OUT.stat().st_size / 1024**2:.2f} MB")
