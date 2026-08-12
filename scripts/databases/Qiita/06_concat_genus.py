#!/usr/bin/env python
"""Qiita Layer1 步06(B): 拼所有 per-study genus 表 → qiita.gg2.genus.h5ad。
outer-join var_id(6级路径), vstack 样本; counts/全样本/不筛不归一(按 merge 契约)。"""
import os, glob, numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad

ROOT = "/hpcdisk1/limk_group/caiqy/project/260428_greengene2"
STUDIES = f"{ROOT}/rawdata/Qiita/work/studies"
OUT = f"{ROOT}/results/feature_table/qiita.gg2.genus.h5ad"

paths = sorted(glob.glob(f"{STUDIES}/*/genus.h5ad"))
print(f"per-study genus 表: {len(paths)}", flush=True)
allvar = set()
for p in paths:
    allvar |= set(ad.read_h5ad(p, backed='r').var_names)
var_list = sorted(allvar); v2i = {v: i for i, v in enumerate(var_list)}
print(f"union genus var: {len(var_list):,}", flush=True)

obs_index = []; blocks = []
for k, p in enumerate(paths):
    a = ad.read_h5ad(p)
    cols = np.array([v2i[v] for v in a.var_names], dtype=np.int64)
    X = a.X.tocoo()
    blocks.append(sp.csr_matrix((X.data, (X.row, cols[X.col])), shape=(a.n_obs, len(var_list)), dtype=np.int32))
    obs_index.extend(list(a.obs_names))
    if (k + 1) % 100 == 0: print(f"  ...{k+1}/{len(paths)} 累计样本 {len(obs_index):,}", flush=True)

X_all = sp.vstack(blocks, format='csr')
var_df = pd.DataFrame([v.split(';') for v in var_list],
                      columns=['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus'],
                      index=pd.Index(var_list, name='var_id')).astype(str)
adata = ad.AnnData(X=X_all, obs=pd.DataFrame(index=pd.Index(obs_index, name='sample_name')), var=var_df)
dup = pd.Index(obs_index).duplicated().sum()
print(adata); print(f"重复 obs 名: {dup} (应0); 总 counts {int(X_all.sum()):,}")
os.makedirs(os.path.dirname(OUT), exist_ok=True)
adata.write_h5ad(OUT, compression='gzip')
print(f"已保存: {OUT}")
