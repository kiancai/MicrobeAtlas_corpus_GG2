#!/usr/bin/env python
"""Qiita 全样本净新增评估（不碰 BIOM，仅用 templates 元数据 + ENA 桥表）。

输入:
  rawdata/Qiita/templates/                 # tar 里抽出的 sample-info + prep-info 模板
  rawdata/Qiita/qiita_ena_bridge.tsv       # ENA UCSDMI: sample_alias(=qiita 样本名) -> run_accession + first_public
  /tmp/corpus_runs.sorted.txt              # 当前语料 obs['Run'] (1,793,428)
  results/feature_table/microbiomap.gg2.genus.h5ad  # MicroBioMap 重注释产物, obs_names=run

逻辑: Qiita 16S 样本(任一 prep 为 16S) -> 经桥表拿 run -> 与 (corpus ∪ MicroBioMap) run 取差集
      -> 按环境(empo)/年份/人体拆分; 无桥 run 的算"净新增(无ENA桥, 未验证)"。
"""
import os, glob, csv, sys
from collections import Counter, defaultdict

ROOT = "/hpcdisk1/limk_group/caiqy/project/260428_greengene2"
TPL = f"{ROOT}/rawdata/Qiita/templates"
BRIDGE = f"{ROOT}/rawdata/Qiita/qiita_ena_bridge.tsv"
CORPUS = "/tmp/corpus_runs.sorted.txt"
MBM = f"{ROOT}/results/feature_table/microbiomap.gg2.genus.h5ad"
OUT = f"{ROOT}/rawdata/Qiita/qiita_netnew_16s.tsv"

def rows(fn):
    with open(fn, encoding="utf-8", errors="replace") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(hdr)}
        for line in fh:
            yield idx, line.rstrip("\n").split("\t")

def is16s(g):
    g = (g or "").lower().replace(" ", "")
    return g.startswith("16s") or g == "16srrna"

prep_files = glob.glob(f"{TPL}/*_prep_*.txt")
samp_files = [f for f in glob.glob(f"{TPL}/*.txt") if "_prep_" not in os.path.basename(f)]
print(f"prep 文件 {len(prep_files)} / 样本信息文件 {len(samp_files)}", flush=True)

# 1) sample_name -> 是否 16S (任一 prep 为 16S)
sn_16s = set()
sn_anygene = set()
for pf in prep_files:
    first = True
    for idx, r in rows(pf):
        si = idx.get("sample_name"); gi = idx.get("target_gene")
        if si is None or len(r) <= si:
            continue
        sn = r[si]; sn_anygene.add(sn)
        if gi is not None and len(r) > gi and is16s(r[gi]):
            sn_16s.add(sn)
print(f"出现在 prep 的样本 {len(sn_anygene):,}; 其中 16S 样本 {len(sn_16s):,}", flush=True)

# 2) sample_name -> 环境/host/study (来自样本信息文件)
sn_env1 = {}; sn_env3 = {}; sn_type = {}; sn_htaxid = {}; sn_study = {}
for sf in samp_files:
    for idx, r in rows(sf):
        si = idx.get("sample_name")
        if si is None or len(r) <= si:
            continue
        sn = r[si]
        def g(col):
            i = idx.get(col)
            return r[i] if (i is not None and len(r) > i and r[i] not in ("", "not applicable", "not provided", "Not applicable")) else None
        e1 = g("empo_1") or g("qiita_empo_1")
        e3 = g("empo_3") or g("qiita_empo_3")
        if e1: sn_env1[sn] = e1
        if e3: sn_env3[sn] = e3
        t = g("sample_type");
        if t: sn_type[sn] = t
        h = g("host_taxid");
        if h: sn_htaxid[sn] = h
        st = g("qiita_study_id")
        if st: sn_study[sn] = st

# 3) 桥表: sample_name -> runs, year
sn_runs = defaultdict(list); run_year = {}
with open(BRIDGE, encoding="utf-8", errors="replace") as fh:
    h = fh.readline().rstrip("\n").split("\t")
    ai = h.index("sample_alias"); ri = h.index("run_accession"); fi = h.index("first_public")
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) <= fi: continue
        alias = p[ai]
        sn = alias.split(":", 1)[1] if (":" in alias and alias.startswith("qiita_sid_")) else alias
        sn_runs[sn].append(p[ri]); run_year[p[ri]] = (p[fi] or "")[:4]
print(f"桥表覆盖样本 {len(sn_runs):,}", flush=True)

# 4) 去重基线
corpus = set(open(CORPUS).read().split())
import anndata as ad
mbm = set(map(str, ad.read_h5ad(MBM, backed="r").obs_names))
print(f"corpus run {len(corpus):,} / MicroBioMap run {len(mbm):,}", flush=True)

# 5) 逐 16S 样本判定
status = Counter()
nn_env1 = Counter(); nn_env3 = Counter(); nn_year = Counter(); nn_human = Counter()
fout = open(OUT, "w")
fout.write("sample_name\tqiita_study_id\tstatus\trun\tfirst_public\tempo_1\tempo_3\tsample_type\thost_taxid\n")
for sn in sn_16s:
    runs = sn_runs.get(sn, [])
    incorp = [r for r in runs if r in corpus]
    inmbm = [r for r in runs if r in mbm]
    if incorp:
        st = "in_corpus"; run = incorp[0]
    elif inmbm:
        st = "in_microbiomap"; run = inmbm[0]
    elif runs:
        st = "netnew_bridged"; run = runs[0]
    else:
        st = "netnew_nobridge"; run = ""
    status[st] += 1
    e1 = sn_env1.get(sn, "NA"); e3 = sn_env3.get(sn, "NA")
    yr = run_year.get(run, "NA") if run else "NA"
    if st.startswith("netnew"):
        nn_env1[e1] += 1; nn_env3[e3] += 1; nn_year[yr] += 1
        nn_human["human" if sn_htaxid.get(sn) == "9606" else "non-human/unknown"] += 1
    fout.write(f"{sn}\t{sn_study.get(sn,'')}\t{st}\t{run}\t{yr}\t{e1}\t{e3}\t{sn_type.get(sn,'')}\t{sn_htaxid.get(sn,'')}\n")
fout.close()

def show(title, c, n=15):
    print(f"\n{title}")
    tot = sum(c.values())
    for k, v in c.most_common(n):
        print(f"  {v:>8,}  {k}")
    if len(c) > n: print(f"  ... (+{len(c)-n} more)")

print("\n" + "=" * 60)
print(f"Qiita 16S 样本总数: {len(sn_16s):,}")
show("按去重状态:", status)
netnew = status["netnew_bridged"] + status["netnew_nobridge"]
print(f"\n>>> 净新增 16S 样本(bridged + nobridge) = {netnew:,}")
print(f"    其中有 ENA 桥可定位 run = {status['netnew_bridged']:,}")
print(f"    无桥(Qiita-only / 非UCSDMI缴存, 未验证) = {status['netnew_nobridge']:,}")
show("净新增 · 按 empo_1(宿主/环境大类):", nn_env1)
show("净新增 · 按 empo_3(细分部位/环境):", nn_env3, 20)
show("净新增 · 人体 vs 非人体(host_taxid):", nn_human)
show("净新增(bridged) · 按 first_public 年份:", nn_year, 20)
print(f"\n清单已写: {OUT}")
