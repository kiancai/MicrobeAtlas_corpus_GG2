#!/usr/bin/env python
"""Qiita 数据全貌 profiling: study×测序方法混合 / 16S study 数 / metadata 覆盖 / 能注释的16S总量。
不管去重——注释阶段全注释, 重复是下游筛选的事。"""
import os, glob, h5py
from collections import Counter, defaultdict

ROOT = "/hpcdisk1/limk_group/caiqy/project/260428_greengene2"
TPL = f"{ROOT}/rawdata/Qiita/templates"
BIOMDIR = f"{ROOT}/rawdata/Qiita/BIOM"

def rows(fn):
    with open(fn, encoding="utf-8", errors="replace") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(hdr)}
        for line in fh:
            yield idx, line.rstrip("\n").split("\t")

def norm_gene(g):
    g = (g or "").lower().replace(" ", "").replace("_", "")
    if g.startswith("16s"): return "16S"
    if g.startswith("18s"): return "18S"
    if "its" in g: return "ITS"
    if "lsu" in g or "largesubunit" in g or "23s" in g or "28s" in g: return "LSU/23S"
    if g in ("", "notapplicable", "notprovided", "xxebixx", "na"): return "NA"
    if "metagenom" in g or "shotgun" in g: return "WGS/meta"
    return "other"

prep_files = glob.glob(f"{TPL}/*_prep_*.txt")
samp_files = [f for f in glob.glob(f"{TPL}/*.txt") if "_prep_" not in os.path.basename(f)]

# ---- A. study × 测序方法 ----
study_genes = defaultdict(Counter)     # study -> gene -> #样本(按 prep 行)
sn_16s = set()
sn_anyprep = set()
for pf in prep_files:
    study = os.path.basename(pf).split("_prep_")[0]
    for idx, r in rows(pf):
        si, gi = idx.get("sample_name"), idx.get("target_gene")
        if si is None or len(r) <= si: continue
        sn = r[si]; sn_anyprep.add(sn)
        g = norm_gene(r[gi]) if (gi is not None and len(r) > gi) else "NA"
        study_genes[study][g] += 1
        if g == "16S": sn_16s.add(sn)

studies = set(study_genes)
def has(s, g): return study_genes[s].get(g, 0) > 0
s_16s = [s for s in studies if has(s, "16S")]
others = lambda s: {g for g in study_genes[s] if g not in ("16S", "NA")}
s_16only = [s for s in s_16s if not others(s)]
s_mixed = [s for s in s_16s if others(s)]
s_no16 = [s for s in studies if not has(s, "16S")]

print("=" * 60)
print("【A. study × 测序方法】")
print(f"有 prep 的 study 总数: {len(studies)}")
print(f"  含 16S 的 study: {len(s_16s)}")
print(f"    纯 16S(只有16S): {len(s_16only)}")
print(f"    16S + 其它方法混合: {len(s_mixed)}")
print(f"  完全没有 16S 的 study: {len(s_no16)}")
gene_in_studies = Counter()
for s in studies:
    for g in study_genes[s]:
        if g != "NA": gene_in_studies[g] += 1
print("  各方法出现在多少 study:")
for g, n in gene_in_studies.most_common():
    print(f"    {g}: {n} study")
# 混合 study 举例
print("  混合 study 举例(study: 方法分布):")
for s in s_mixed[:6]:
    print(f"    {s}: {dict(study_genes[s])}")

# ---- B. metadata 覆盖(对 16S 样本) ----
KEYCOLS = ["empo_1", "empo_3", "qiita_empo_3", "host_taxid", "sample_type",
           "scientific_name", "sample_scientific_name", "country", "geo_loc_name",
           "collection_timestamp", "host_subject_id", "latitude", "longitude",
           "host_age", "sex", "env_material", "host_body_site"]
cover = Counter()
env1 = Counter(); env3 = Counter()
sn_haveinfo = set()
BAD = {"", "not applicable", "not provided", "not collected", "missing",
       "unknown", "Not applicable", "NA", "na", "none", "None", "no_data"}
for sf in samp_files:
    for idx, r in rows(sf):
        si = idx.get("sample_name")
        if si is None or len(r) <= si: continue
        sn = r[si]
        if sn not in sn_16s: continue
        sn_haveinfo.add(sn)
        for c in KEYCOLS:
            i = idx.get(c)
            if i is not None and len(r) > i and r[i] not in BAD:
                cover[c] += 1
        def gv(c):
            i = idx.get(c); return r[i] if (i is not None and len(r) > i and r[i] not in BAD) else None
        e1 = gv("empo_1") or gv("qiita_empo_1"); e3 = gv("empo_3") or gv("qiita_empo_3")
        if e1: env1[e1] += 1
        if e3: env3[e3] += 1

n16 = len(sn_16s)
print("\n" + "=" * 60)
print(f"【B. metadata 覆盖】(分母 = 16S 样本 {n16:,}; 有样本信息的 {len(sn_haveinfo):,})")
for c in KEYCOLS:
    print(f"  {c:24s}: {cover[c]:>8,}  ({100*cover[c]/n16:4.1f}%)")

print("\n【环境大类 empo_1 (全部 16S 样本)】")
for k, v in env1.most_common(12): print(f"  {v:>8,}  {k}")
print("\n【细分 empo_3 (全部 16S 样本, top 20)】")
for k, v in env3.most_common(20): print(f"  {v:>8,}  {k}")

# ---- C. 能注释的 16S 总量(有 Deblur ASV 表的 16S 样本) ----
seen_asv = set()
nbiom = 0
for a in os.listdir(BIOMDIR):
    d = f"{BIOMDIR}/{a}"
    bt = "reference-hit.biom" if os.path.exists(f"{d}/reference-hit.biom") else ("all.biom" if os.path.exists(f"{d}/all.biom") else None)
    if bt is None: continue
    try:
        with h5py.File(f"{d}/{bt}", "r") as hf:
            oids = hf["observation/ids"]
            if oids.shape[0] == 0: continue
            if oids[0].decode()[:20].upper().strip("ACGTN") != "" or len(oids[0].decode()) < 50:
                continue  # 非序列表
            for s in hf["sample/ids"][:]:
                ss = s.decode()
                if ss in sn_16s: seen_asv.add(ss)
        nbiom += 1
    except Exception:
        pass

print("\n" + "=" * 60)
print("【C. 能注释的 16S 总量】(注释全注释, 不管重复)")
print(f"  16S 样本总数(有16S prep): {n16:,}")
print(f"  其中有 Deblur ASV 表 = 能注释: {len(seen_asv):,}")
print(f"  只有 OTU 表/无序列 = 注释不了: {n16-len(seen_asv):,}")
print(f"  (扫描的 Deblur 序列表 artifact: {nbiom:,})")
