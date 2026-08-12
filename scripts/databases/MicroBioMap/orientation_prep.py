#!/usr/bin/env python
"""MBM auto-vs-both spot-check 数据准备：
- 选 ~15 个跨 V 区/规模的 project，重建 per-project asv.fa + per-ASV reads（从 asv.h5ad 切）
- 全局随机 60k ASV 样本 + 全局 per-ASV reads
orientation 是 per-project 属性 → per-project auto 才是 both 的公平对手。
"""
import os, random
import numpy as np
import h5py
from scipy.sparse import csr_matrix

ROOT = "/hpcdisk1/limk_group/caiqy/project/260428_greengene2"
H5 = f"{ROOT}/results/feature_table/microbiomap.gg2.asv.h5ad"
META = f"{ROOT}/rawdata/MicroBioMap/work/sample_metadata.tsv"
PROJCSV = f"{ROOT}/rawdata/MicroBioMap/work/projects.csv"
ASVFA = f"{ROOT}/rawdata/MicroBioMap/work/asv_global.fa"
OUT = f"{ROOT}/results/orientation_check_mbm"
os.makedirs(f"{OUT}/samples", exist_ok=True)
os.makedirs(f"{OUT}/tax", exist_ok=True)
os.makedirs(f"{OUT}/logs", exist_ok=True)

# ---------- 1. sample(srr) -> project ----------
srr2proj = {}
with open(META) as f:
    hdr = f.readline().rstrip("\n").split("\t")
    i_srr, i_proj = hdr.index("srr"), hdr.index("project")
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) > max(i_srr, i_proj):
            srr2proj[p[i_srr]] = p[i_proj]
print(f"sample_metadata: {len(srr2proj):,} srr→project")

# project -> region
proj_region = {}
with open(PROJCSV, encoding="utf-8-sig") as f:
    hdr = f.readline().rstrip("\n").split(",")
    i_p, i_a = hdr.index("project"), hdr.index("amplicon")
    for line in f:
        p = line.rstrip("\n").split(",")
        if len(p) > max(i_p, i_a):
            proj_region[p[i_p]] = p[i_a]

# ---------- 2. 读 asv.h5ad 的 X / obs / var（用 h5py，省内存） ----------
print("读 asv.h5ad ...")
with h5py.File(H5, "r") as h:
    g = h["X"]
    shape = tuple(int(x) for x in g.attrs["shape"])
    X = csr_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=shape)
    oidx = h["obs"].attrs.get("_index", "_index")
    obs_names = np.array([x.decode() if isinstance(x, bytes) else x for x in h["obs"][oidx][:]])
    vidx = h["var"].attrs.get("_index", "_index")
    var_names = np.array([x.decode() if isinstance(x, bytes) else x for x in h["var"][vidx][:]])
print(f"X {X.shape}  nnz={X.nnz:,}")

# obs row -> project
row_proj = np.array([srr2proj.get(s, "NA") for s in obs_names])
from collections import Counter
proj_rows = {}
for i, p in enumerate(row_proj):
    proj_rows.setdefault(p, []).append(i)
proj_sizes = {p: len(r) for p, r in proj_rows.items() if p != "NA"}
print(f"project 数: {len(proj_sizes)}")

# ---------- 3. 选 ~15 个 project：按 region 分层 + 规模适中 ----------
# 先算每个 project 的 ASV 数（count>0），挑 ASV 数在 [800, 60000] 的
cand = []
for p, rows in proj_rows.items():
    if p == "NA" or len(rows) < 5:
        continue
    sub = X[rows]
    nz = (np.asarray(sub.sum(axis=0)).ravel() > 0).sum()
    cand.append((p, len(rows), int(nz), proj_region.get(p, "?")))
cand = [c for c in cand if 800 <= c[2] <= 60000]
# 按 region 分组，每组挑规模分散的几个，总共 ~15
from collections import defaultdict
byreg = defaultdict(list)
for c in cand:
    byreg[c[3]].append(c)
sel = []
regions = sorted(byreg, key=lambda r: -len(byreg[r]))
# 轮流从各 region 取（按 ASV 数排序后等距），直到 15 个
per_reg = max(1, 15 // max(1, len(regions)))
for r in regions:
    lst = sorted(byreg[r], key=lambda x: x[2])
    idx = np.linspace(0, len(lst) - 1, min(per_reg + 1, len(lst))).astype(int)
    for j in sorted(set(idx)):
        sel.append(lst[j])
    if len(sel) >= 18:
        break
# 去重 + 截到 15
seen = set(); final = []
for c in sel:
    if c[0] not in seen:
        seen.add(c[0]); final.append(c)
final = final[:15]
print(f"\n选中 {len(final)} 个 project (project, n_samp, n_asv, region):")
for p, ns, na, reg in final:
    print(f"  {p:>14}  samp={ns:>5}  asv={na:>6}  region={reg}")

# ---------- 4. 重建 per-project counts + 收集需要的 ASV id ----------
needed = defaultdict(list)   # asv_id -> [dest_label,...]
proj_counts = {}             # label -> dict(asv_id->reads)
for p, ns, na, reg in final:
    rows = proj_rows[p]
    col = np.asarray(X[rows].sum(axis=0)).ravel()
    nz = np.where(col > 0)[0]
    ids = var_names[nz]; reads = col[nz]
    label = f"proj_{p}"
    proj_counts[label] = dict(zip(ids, reads.astype(np.int64)))
    for a in ids:
        needed[a].append(label)

# ---------- 5. 全局随机 60k ASV + 全局 reads ----------
gcol = np.asarray(X.sum(axis=0)).ravel()
rng = np.random.default_rng(42)
gidx = rng.choice(len(var_names), size=60000, replace=False)
gids = var_names[gidx]
greads = gcol[gidx]
gcounts = dict(zip(gids, greads.astype(np.int64)))
for a in gids:
    needed[a].append("global60k")
del X, gcol

# ---------- 6. 流式扫 asv_global.fa，写各 fasta ----------
print(f"\n需要 {len(needed):,} 条唯一 ASV 序列，流式扫 asv_global.fa ...")
handles = {}
def H(label):
    if label not in handles:
        handles[label] = open(f"{OUT}/samples/{label}.fa", "w")
    return handles[label]
cur_id = None
written = 0
with open(ASVFA) as f:
    for line in f:
        if line.startswith(">"):
            cur_id = line[1:].strip().split()[0]
        else:
            if cur_id in needed:
                for label in needed[cur_id]:
                    h = H(label)
                    h.write(f">{cur_id}\n{line if line.endswith(chr(10)) else line+chr(10)}")
                written += 1
for h in handles.values():
    h.close()
print(f"写出序列条目 {written:,}（含跨 project 重复）")

# ---------- 7. 写 counts tsv + manifest ----------
for label, d in list(proj_counts.items()) + [("global60k", gcounts)]:
    with open(f"{OUT}/samples/{label}.counts.tsv", "w") as f:
        f.write("asv_id\treads\n")
        for a, r in d.items():
            f.write(f"{a}\t{r}\n")

labels = ["global60k"] + [f"proj_{p}" for p, _, _, _ in final]
with open(f"{OUT}/manifest.tsv", "w") as f:
    for i, lab in enumerate(labels):
        f.write(f"{i}\t{lab}\t{OUT}/samples/{lab}.fa\t{lab}\n")
print(f"\nmanifest: {len(labels)} 行 → {OUT}/manifest.tsv")
print("准备完成。")
