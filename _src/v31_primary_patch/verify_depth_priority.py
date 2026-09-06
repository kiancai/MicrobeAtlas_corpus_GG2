"""Independent full-row selection and summary check for the second-round trial."""
from pathlib import Path
import hashlib
import json
import h5py
from anndata.io import read_elem
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[4]
EXP = ROOT/'tmp/20260905_mcfcorpus_v31_primary_run_patch'
OUT = EXP/'intermediate/depth_priority'
AUD = EXP/'audit/depth_priority'

def sha(p):
    h = hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(8*1024**2), b''): h.update(b)
    return h.hexdigest()

def main():
    for p,h in json.loads((AUD/'manifest.json').read_text()).items(): assert sha(ROOT/p)==h
    with h5py.File(ROOT/'data/MCFCorpusV3.1/MCFCorpusV3.1.gg2.h5ad','r') as f:
        df = pd.DataFrame({c:read_elem(f['obs'][c]) for c in ['Run','Database','qc_pass','n_deep_reads']})
    for c in ['Run','Database']: df[c]=df[c].astype(object).where(df[c].notna(),'').astype(str)
    df.Run=df.Run.replace({'nan':'','<NA>':'','None':''})
    df['row']=np.arange(len(df))
    df['rm']=df.Database.eq('RM')
    qc=df[df.qc_pass & df.Run.ne('')]
    # Independent reconstruction of the original priority from the corpus.
    baseline=qc.sort_values(['Run','rm','n_deep_reads','row'], ascending=[True,False,False,True]).drop_duplicates('Run')
    bi=baseline.row.to_numpy()
    assert np.array_equal(bi,np.load(OUT/'baseline_indices.npy'))
    assert np.array_equal(np.flatnonzero(df.qc_pass & df.Run.eq('')),np.load(OUT/'qc_no_run_indices.npy'))
    identity=pd.read_csv(EXP/'intermediate/run_selection_comparison.tsv.gz',sep='\t',usecols=['Run','study_key'],keep_default_na=False).set_index('Run').study_key
    baseline=baseline.copy(); baseline['study']=baseline.Run.map(identity)
    selection_runs=set(baseline.loc[~baseline.rm & baseline.study.ne(''),'Run'])
    pool=qc[qc.Run.isin(selection_runs)].copy();pool['study']=pool.Run.map(identity)
    multi=pool.groupby('study').Database.nunique().loc[lambda s:s.gt(1)].index
    pool=pool[pool.study.isin(multi)]
    maxdepth=pool.groupby('Run').n_deep_reads.max()
    dup=set(pool.groupby('Run').size().loc[lambda s:s.gt(1)].index)
    primary={}
    # Independently calculate one score per source in a study, including empty shared sets.
    for study, rows in pool.groupby('study'):
        ranking=[]
        for db, source in rows.groupby('Database'):
            shared=source[source.Run.isin(dup)]
            ratios=shared.n_deep_reads.to_numpy()/shared.Run.map(maxdepth).to_numpy()
            ranking.append((-float(np.median(ratios)) if len(ratios) else 0., db))
        primary[study]=min(ranking)[1]
    best={}
    for row in pool.itertuples():
        if row.Database!=primary[row.study]: continue
        current=best.get(row.Run)
        rank=(-row.n_deep_reads,row.row)
        if current is None or rank<current: best[row.Run]=rank
    expected=bi.copy()
    for i,run in enumerate(baseline.Run):
        if run in best: expected[i]=best[run][1]
    assert np.array_equal(expected,np.load(OUT/'study_strict_indices.npy'))
    dep=df.n_deep_reads.to_numpy(); db=df.Database.to_numpy()
    bd=dep[bi]; scoped=baseline.study.isin(multi).to_numpy() & ~baseline.rm.to_numpy()
    summary=pd.read_csv(OUT/'summary.tsv',sep='\t').set_index('variant')
    for variant, s in summary.iterrows():
        idx=np.load(OUT/f'{variant}_indices.npy')
        target=bi if variant=='baseline' else expected.copy()
        if variant.startswith('study_loss_'):
            loss=int(variant.split('_')[2][:-3])/100
            target=np.where(dep[expected]/bd>=1-loss,expected,bi)
        assert np.array_equal(idx,target)
        assert int((idx!=bi).sum())==s.changed_runs
        assert int((dep[idx]<bd).sum())==s.shallower_runs
        counts=pd.crosstab(baseline.study.to_numpy()[scoped],db[idx][scoped])
        assert int((counts.gt(0).sum(axis=1)>1).sum())==s.mixed_studies
        assert int((counts.sum(axis=1)-counts.max(axis=1)).sum())==s.minority_runs
    result=dict(status='PASS', full_known_runs=len(bi), qc_no_run=int((df.qc_pass & df.Run.eq('')).sum()), study_groups=len(multi), scoped_runs=int(scoped.sum()), variants=len(summary), verifier_sha256=sha(Path(__file__)))
    (AUD/'independent_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)

if __name__=='__main__': main()
