#!/usr/bin/env python3
"""Independent relational checks of the exploratory row selections."""
from pathlib import Path
import hashlib
import json
import h5py
from anndata.io import read_elem
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[4]
EXP=ROOT/'tmp/20260905_mcfcorpus_v31_primary_run_patch'
CORPUS=ROOT/'data/MCFCorpusV3.1/MCFCorpusV3.1.gg2.h5ad'


def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(8*1024**2),b''): h.update(chunk)
    return h.hexdigest()


def main():
    s=json.loads((EXP/'audit/summary.json').read_text())
    assert digest(CORPUS)==s['input_sha256']
    m=json.loads((EXP/'audit/output_manifest.json').read_text())
    for p,h in m['files_sha256'].items(): assert digest(EXP/p)==h,p
    for p,h in s['supporting_inputs_sha256'].items(): assert digest(ROOT/p)==h,p
    with h5py.File(CORPUS) as f:
        d=pd.DataFrame({k:read_elem(f['obs'][k]) for k in ['Database','Run','qc_pass','n_deep_reads','biome_depth']})
    for k in ['Database','Run']: d[k]=d[k].astype(object).where(pd.notna(d[k]),'').astype(str)
    d['row']=np.arange(len(d))
    q=d[d.qc_pass].copy();q.Run=q.Run.replace({'nan':'','None':'','<NA>':''})
    w=q[q.Run.ne('')].copy();n=q[q.Run.eq('')]
    c=pd.read_csv(EXP/'intermediate/run_selection_comparison.tsv.gz',sep='\t',keep_default_na=False)
    assert c.Run.is_unique and c.obs_index.is_unique and c.baseline_obs_index.is_unique
    assert set(c.Run)==set(w.Run)
    chosen=d.loc[c.obs_index];old=d.loc[c.baseline_obs_index]
    assert chosen.qc_pass.all() and old.qc_pass.all()
    assert np.array_equal(chosen.Run.to_numpy(),c.Run.to_numpy())
    assert np.array_equal(old.Run.to_numpy(),c.Run.to_numpy())
    assert np.array_equal(chosen.Database.to_numpy(),c.Database.to_numpy())
    assert np.array_equal(chosen.n_deep_reads.to_numpy(),c.n_deep_reads.to_numpy())
    rm=q[q.Database.eq('RM')]
    assert set(rm.row)<=set(c.obs_index) and len(rm)==s['rm_priority_runs']
    rm_runs=set(rm.Run)
    wr=w[~w.Run.isin(rm_runs)].copy()
    maxima=wr.groupby('Run').n_deep_reads.max()
    oldby=old.set_index('Run')
    assert np.array_equal(oldby.loc[maxima.index,'n_deep_reads'],maxima)
    # Canonical group assignment is explicitly exported; these checks do not assert a biological cohort identity.
    lookup=c.set_index('Run')
    wr['study']=wr.Run.map(lookup.study_key)
    usable=wr[wr.study.ne('')].copy()
    cv=pd.read_csv(EXP/'intermediate/study_source_coverage.tsv',sep='\t')
    expected=usable.groupby(['study','Database']).Run.nunique()
    got=cv.set_index(['study_key','Database']).covered_runs.sort_index()
    assert np.array_equal(expected.sort_index().values,got.values)
    assert list(expected.sort_index().index)==list(got.index)
    leaders=pd.read_csv(EXP/'intermediate/study_policy_comparison.tsv',sep='\t').set_index('study_key')
    for k, g in cv.groupby('study_key'):
        lead=leaders.loc[k]
        assert lead.covered_runs==g.covered_runs.max()
        tied=g[g.covered_runs.eq(g.covered_runs.max())]
        best=tied[tied.depth_tie_score.eq(tied.depth_tie_score.max())].Database.min()
        assert lead.Database==best
    usable['primary']=usable.study.map(leaders.Database)
    available=usable[usable.Database.eq(usable.primary)]
    primary_runs=set(available.Run)
    p=c[c.selection_reason.eq('study_primary')]
    assert set(p.Run)==primary_runs
    assert (p.Database==p.study_key.map(leaders.Database)).all()
    fallback=c[c.selection_reason.isin(['primary_missing_run_depth_fallback','unresolved_study_depth_fallback'])]
    assert np.array_equal(fallback.obs_index,fallback.baseline_obs_index)
    assert not set(fallback.Run)&primary_runs
    delta=(c.obs_index!=c.baseline_obs_index)
    changes=pd.read_csv(EXP/'intermediate/changed_run_choices.tsv.gz',sep='\t')
    assert set(changes.Run)==set(c.loc[delta,'Run']) and len(changes)==s['changed_runs']
    assert len(chosen)+len(n)==s['policies']['study_aware_trial']['selected_rows']
    assert int((chosen.biome_depth>0).sum()+(n.biome_depth>0).sum())==s['policies']['study_aware_trial']['any_biome']
    checks=['full corpus/input/output hashes','exact Run unique/full coverage','selected obs identity and QC',
            'all RM QC rows retained','baseline max depth independently checked',
            'study coverage reconstructed','primary source ranking checked',
            'primary availability and baseline fallback checked','changed Run set','no-Run accounting and biome totals']
    result={'status':'PASS','scope':'first-round exploratory selection; not a released patch',
            'checks':checks,'verifier_sha256':digest(Path(__file__)),'runs_checked':len(c),'new_qc_or_biome_changes':False}
    (EXP/'audit/independent_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
