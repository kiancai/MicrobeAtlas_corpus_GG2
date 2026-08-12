#!/usr/bin/env python
"""Qiita Layer1 步0a(v2): 每个 16S 样本选一个 Deblur 表(90nt reference-hit 优先),
**只在 16S 表里选**(混合 study 用 artifact↔prep 的 Jaccard 锁定 gene=16S,纯16S study 全部算16S)。
输出 sample→artifact 指派 + 覆盖报告。"""
import os, glob, h5py
from collections import Counter, defaultdict

ROOT = "/hpcdisk1/limk_group/caiqy/project/260428_greengene2"
BIOMDIR = f"{ROOT}/rawdata/Qiita/BIOM"
TPL = f"{ROOT}/rawdata/Qiita/templates"
NETNEW = f"{ROOT}/rawdata/Qiita/qiita_netnew_16s.tsv"
OUT = f"{ROOT}/rawdata/Qiita/work/sample_artifact_map.tsv"
os.makedirs(os.path.dirname(OUT), exist_ok=True)

def norm(g):
    g = (g or "").lower().replace(" ", "").replace("_", "")
    if g.startswith("16s"): return "16S"
    if g.startswith("18s"): return "18S"
    if "its" in g: return "ITS"
    if "lsu" in g or "23s" in g or "28s" in g: return "LSU"
    if g in ("", "notapplicable", "notprovided", "xxebixx", "na"): return "NA"
    return "other"

# 16S 样本全集
sn16 = set()
with open(NETNEW) as f:
    h = f.readline().rstrip("\n").split("\t"); si = h.index("sample_name")
    for line in f:
        sn16.add(line.split("\t")[si])
print(f"16S 样本全集: {len(sn16):,}", flush=True)

# prep info: study -> [(gene, frozenset(samples))]; study gene 集合
study_preps = defaultdict(list)
study_genes = defaultdict(set)
for pf in glob.glob(f"{TPL}/*_prep_*.txt"):
    study = os.path.basename(pf).split("_prep_")[0]
    with open(pf, encoding="utf-8", errors="replace") as fh:
        hd = fh.readline().rstrip("\n").split("\t"); idx = {c: i for i, c in enumerate(hd)}
        si2, gi = idx.get("sample_name"), idx.get("target_gene")
        if si2 is None: continue
        samps = set(); genes = set()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) <= si2: continue
            samps.add(p[si2])
            genes.add(norm(p[gi]) if (gi is not None and len(p) > gi) else "NA")
    g = next((x for x in genes if x != "NA"), "NA")
    study_preps[study].append((g, frozenset(samps)))
    study_genes[study].add(g)
mixed = {s for s, gs in study_genes.items() if "16S" in gs and (gs - {"16S", "NA"})}
pure16 = {s for s, gs in study_genes.items() if "16S" in gs and not (gs - {"16S", "NA"})}
print(f"study: 纯16S {len(pure16)} / 混合(含16S+其它) {len(mixed)} / 无16S {len(study_genes)-len(pure16)-len(mixed)}", flush=True)

def art_gene(study, S):
    if study in pure16: return "16S"
    if study not in mixed: return None
    best_g, best_j = None, -1.0
    for g, ps in study_preps[study]:
        u = len(S | ps)
        j = (len(S & ps) / u) if u else 0.0
        if j > best_j: best_j, best_g = j, g
    return best_g

TRIMPRI = {90: 0, 100: 1, 150: 2, 125: 3, 250: 4, 200: 5, 300: 6}
def artkey(trim, bt): return (0 if bt == "reference-hit.biom" else 1, TRIMPRI.get(trim, 9))

best = {}
skipped_non16 = 0
done = 0
for a in os.listdir(BIOMDIR):
    for bt in ("reference-hit.biom", "all.biom"):
        p = f"{BIOMDIR}/{a}/{bt}"
        if not os.path.exists(p): continue
        try:
            with h5py.File(p, "r") as hf:
                oids = hf["observation/ids"]
                if oids.shape[0] == 0: continue
                first = oids[0].decode()
                if len(first) < 50 or first[:20].upper().strip("ACGTN") != "": continue
                trim = len(first)
                S = {s.decode() for s in hf["sample/ids"][:]}
        except Exception:
            continue
        study = next(iter(S)).split(".")[0] if S else ""
        g = art_gene(study, S)
        if g != "16S":
            skipped_non16 += 1
            continue
        k = artkey(trim, bt)
        for ss in S:
            if ss in sn16 and (ss not in best or k < best[ss][0]):
                best[ss] = (k, a, bt, trim)
    done += 1
    if done % 4000 == 0: print(f"  ...{done} artifact", flush=True)

with open(OUT, "w") as fo:
    fo.write("sample_name\tartifact\tbiom_type\ttrim\n")
    for s, (k, a, bt, trim) in best.items():
        fo.write(f"{s}\t{a}\t{bt}\t{trim}\n")

print("\n" + "=" * 55)
print(f"覆盖 16S 样本: {len(best):,} / {len(sn16):,}")
print(f"跳过的非16S表(混合study里的18S/ITS): {skipped_non16:,}")
print("选中表 trim(按样本):", dict(Counter(trim for k, a, bt, trim in best.values())))
print("选中表类型(按样本):", dict(Counter(bt for k, a, bt, trim in best.values())))
print(f"被选中 artifact 数: {len(set((a,bt) for k,a,bt,trim in best.values())):,}")
print(f"指派表写: {OUT}")
