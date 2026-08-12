#!/usr/bin/env python
"""Qiita Layer1 步05(B): per-study ASV → genus collapse,用全局 taxonomy。
每 study: asv.h5ad(var=asv_i) + asv.fa(asv_i→seq) → md5(seq) → 全局 taxonomy(md5→Taxon)
→ parse_gg2 → 过滤(B/A, 去mito/chloro) → 6级 var_id 聚合 → genus.h5ad。
collapse 逻辑 = MBM05/RM03 逐字。一次性加载全局 taxonomy, 循环所有 study(可断点续)。"""
import os, glob, hashlib, numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad

ROOT = "/hpcdisk1/limk_group/caiqy/project/260428_greengene2"
STUDIES = f"{ROOT}/rawdata/Qiita/work/studies"
GLOBAL_TAX = f"{ROOT}/rawdata/Qiita/work/taxonomy_gg2.global.tsv"

RANK_COLS = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
PREFIXES  = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
RANK_TO_GENUS, PFX_TO_GENUS = RANK_COLS[:6], PREFIXES[:6]

def parse_gg2(t):
    if t is None or (isinstance(t, float) and pd.isna(t)):
        return {c: None for c in RANK_COLS}
    parts = [p.strip() for p in str(t).split(';')]
    return {c: next((p for p in parts if p.startswith(pfx)), None) for c, pfx in zip(RANK_COLS, PREFIXES)}

def make_var_id(row):
    return ';'.join((row[c] if (row[c] is not None and not (isinstance(row[c], float) and np.isnan(row[c]))) else pfx)
                    for c, pfx in zip(RANK_TO_GENUS, PFX_TO_GENUS))

print("加载全局 taxonomy ...", flush=True)
md5_tax = {}
with open(GLOBAL_TAX) as f:
    next(f)
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            md5_tax[p[0]] = p[1]
print(f"全局 taxonomy: {len(md5_tax):,}", flush=True)

studies = sorted(os.path.basename(os.path.dirname(p)) for p in glob.glob(f"{STUDIES}/*/asv.h5ad"))
print(f"study 数: {len(studies)}", flush=True)
n_missing_tot = 0
for i, study in enumerate(studies):
    d = f"{STUDIES}/{study}"
    out = f"{d}/genus.h5ad"
    if os.path.exists(out):
        continue
    adata = ad.read_h5ad(f"{d}/asv.h5ad")
    # asv_i -> seq
    asv2seq = {}
    with open(f"{d}/asv.fa") as fh:
        hid = None
        for line in fh:
            if line[0] == ">": hid = line[1:].strip()
            else: asv2seq[hid] = line.strip()
    # var(asv_i) -> Taxon(全局)
    taxon = []
    miss = 0
    for v in adata.var_names:
        seq = asv2seq.get(v)
        t = md5_tax.get(hashlib.md5(seq.encode()).hexdigest()) if seq else None
        if t is None: miss += 1
        taxon.append(t)
    n_missing_tot += miss
    parsed = pd.DataFrame([parse_gg2(t) for t in taxon], index=list(adata.var_names))
    parsed['Taxon_full'] = parsed[RANK_COLS].apply(lambda r: ';'.join('' if pd.isna(v) else v for v in r), axis=1)
    keep = (parsed['Domain'].isin(['d__Bacteria', 'd__Archaea'])
            & ~parsed['Taxon_full'].str.contains('mitochondri', case=False, na=False)
            & ~parsed['Taxon_full'].str.contains('chloroplast', case=False, na=False)).values
    tk = parsed[keep].copy()
    tk['var_id'] = tk.apply(make_var_id, axis=1)
    var_list = sorted(tk['var_id'].unique())
    v2i = {v: j for j, v in enumerate(var_list)}
    keep_idx = np.where(keep)[0]
    Xk = adata.X[:, keep_idx].tocsr()
    vidx = np.array([v2i[v] for v in tk['var_id']], dtype=np.int64)
    M = sp.csr_matrix((np.ones(len(keep_idx), np.int32), (vidx, np.arange(len(keep_idx)))), shape=(len(var_list), len(keep_idx)))
    Xv = (Xk @ M.T).tocsr(); Xv.data = Xv.data.astype(np.int32, copy=False)
    var_df = pd.DataFrame([v.split(';') for v in var_list], columns=RANK_TO_GENUS,
                          index=pd.Index(var_list, name='var_id')).astype(str)
    ad.AnnData(X=Xv, obs=pd.DataFrame(index=adata.obs.index.copy()), var=var_df).write(out, compression='gzip')
    if (i + 1) % 100 == 0:
        print(f"  ...{i+1}/{len(studies)} (累计未匹配ASV {n_missing_tot})", flush=True)
print(f"完成。未在全局taxonomy匹配到的ASV总数: {n_missing_tot}(应≈0)")
