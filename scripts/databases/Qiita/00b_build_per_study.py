#!/usr/bin/env python
"""Qiita Layer1 步0b: 按 00a 的 sample→artifact 指派, **按 study 分组**, 逐 study 产出
RM 样式输入: work/studies/<study>/{asv.fa, asv.h5ad}(sample×ASV counts)。内存友好(逐 study)。"""
import os, numpy as np, scipy.sparse as sp, pandas as pd, anndata as ad, biom
from collections import defaultdict

ROOT = "/hpcdisk1/limk_group/caiqy/project/260428_greengene2"
BIOMDIR = f"{ROOT}/rawdata/Qiita/BIOM"
MAP = f"{ROOT}/rawdata/Qiita/work/sample_artifact_map.tsv"
OUTROOT = f"{ROOT}/rawdata/Qiita/work/studies"
os.makedirs(OUTROOT, exist_ok=True)

# study -> {(artifact,biom_type) -> [samples]}
study_arts = defaultdict(lambda: defaultdict(list))
with open(MAP) as f:
    next(f)
    for line in f:
        s, a, bt, trim = line.rstrip("\n").split("\t")
        study_arts[s.split(".")[0]][(a, bt)].append(s)
print(f"study 数: {len(study_arts):,}; 样本 {sum(len(v) for d in study_arts.values() for v in d.values()):,}", flush=True)

done = 0; tot_samp = 0; tot_asv = 0
for study, art_samps in sorted(study_arts.items()):
    d = f"{OUTROOT}/{study}"
    os.makedirs(d, exist_ok=True)
    if os.path.exists(f"{d}/asv.h5ad") and os.path.exists(f"{d}/asv.fa"):
        done += 1; continue
    seq2id = {}; sample_list = []; sidx = {}; blocks = []
    for (a, bt), samps in art_samps.items():
        t = biom.load_table(f"{BIOMDIR}/{a}/{bt}")
        obs = t.ids("observation")
        oidx = np.fromiter((seq2id.setdefault(o, len(seq2id)) for o in obs), dtype=np.int64, count=len(obs))
        spos = {s: i for i, s in enumerate(t.ids("sample"))}
        these = [s for s in samps if s in spos]
        if not these: continue
        cols = [spos[s] for s in these]
        l2g = np.empty(len(these), dtype=np.int64)
        for li, s in enumerate(these):
            if s not in sidx: sidx[s] = len(sample_list); sample_list.append(s)
            l2g[li] = sidx[s]
        sub = t.matrix_data.tocsc()[:, cols].tocoo()
        blocks.append((l2g[sub.col], oidx[sub.row], sub.data))
    if not blocks: continue
    rows = np.concatenate([b[0] for b in blocks])
    cols_ = np.concatenate([b[1] for b in blocks])
    data = np.concatenate([b[2] for b in blocks])
    assert np.allclose(data, np.round(data)), f"{study}: counts 非整数"
    X = sp.csr_matrix((data.astype(np.int32), (rows, cols_)), shape=(len(sample_list), len(seq2id)), dtype=np.int32)
    X.sum_duplicates()
    id2seq = [None] * len(seq2id)
    for seq, i in seq2id.items(): id2seq[i] = seq
    with open(f"{d}/asv.fa", "w") as fo:
        for i, seq in enumerate(id2seq):
            fo.write(f">asv_{i}\n{seq}\n")
    ad.AnnData(X=X,
               obs=pd.DataFrame(index=pd.Index(sample_list, name="sample_name")),
               var=pd.DataFrame(index=pd.Index([f"asv_{i}" for i in range(len(seq2id))], name="asv_id"))
               ).write(f"{d}/asv.h5ad")
    tot_samp += len(sample_list); tot_asv += len(seq2id)
    done += 1
    if done % 100 == 0: print(f"  ...{done}/{len(study_arts)} study", flush=True)

print(f"\n完成 {done} study; 累计样本 {tot_samp:,}; 累计(study内唯一)ASV {tot_asv:,}")
print(f"输出: {OUTROOT}/<study>/{{asv.fa, asv.h5ad}}")
