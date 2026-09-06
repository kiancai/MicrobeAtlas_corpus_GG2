"""Reproducible same-version primary-row patch: build, independently verify, publish."""
from pathlib import Path
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
import shutil
import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem, write_elem

ROOT = Path(__file__).resolve().parents[4]
EXP = ROOT/'tmp/20260905_mcfcorpus_v31_primary_run_patch'
RELEASE = ROOT/'data/MCFCorpusV3.1'
NAME = 'MCFCorpusV3.1.gg2.h5ad'
OLD_HASH = 'b0612b46abd0d734ed501c4c312ba76c18f8e8c61034d33cae66af8ee7ba97e3'
PATCH = 'analysis_primary_depth_rm_20260905'
HISTORY = RELEASE/'history/pre_analysis_primary_20260905'
STAGE = EXP/'staging/primary_patch'
POLICY = 'QC first; exact nonempty Run; eligible RM first; maximum n_deep_reads; original row order for ties; retain every QC-pass no-Run observation. No Study, coverage or biome ranking.'

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''): h.update(b)
    return h.hexdigest()

def dump(p,obj):
    with p.open('x') as f: json.dump(obj,f,indent=2,ensure_ascii=False); f.write('\n')

def obs(path):
    with h5py.File(path,'r') as f:
        d=pd.DataFrame({c:read_elem(f['obs'][c]) for c in ['Database','Run','qc_pass','n_deep_reads','biome_depth']})
    for c in ['Run','Database']: d[c]=d[c].astype(object).where(d[c].notna(),'').astype(str)
    d.Run=d.Run.replace({'nan':'','None':'','<NA>':''})
    return d

def select(d):
    flag=np.zeros(len(d),dtype=bool)
    valid=d[d.qc_pass].copy(); valid['row']=np.flatnonzero(d.qc_pass)
    assert np.isfinite(valid.n_deep_reads).all() and valid.n_deep_reads.ge(0).all()
    flag[valid.loc[valid.Run.eq(''),'row']]=True
    valid=valid[valid.Run.ne('')];valid['rm']=valid.Database.eq('RM')
    win=valid.sort_values(['Run','rm','n_deep_reads','row'],ascending=[True,False,False,True],kind='stable').drop_duplicates('Run')
    flag[win.row]=True
    return flag

def independent_select(d):
    """Streaming tuple minimization, independent of pandas sorting implementation."""
    flag=np.zeros(len(d),dtype=bool); winners={}
    for i,r in enumerate(d.itertuples(index=False)):
        if not r.qc_pass: continue
        if not r.Run: flag[i]=True;continue
        rank=(0 if r.Database=='RM' else 1,-float(r.n_deep_reads),i)
        if r.Run not in winners or rank<winners[r.Run]: winners[r.Run]=rank
    for rank in winners.values(): flag[rank[2]]=True
    return flag

def eq(a,b):
    a,b=np.asarray(a),np.asarray(b)
    return np.array_equal(a,b,equal_nan=True) if a.dtype.kind in 'fc' and b.dtype.kind in 'fc' else np.array_equal(a,b)

def preserve(a,b,records):
    assert type(a)==type(b),(a.name,'type')
    assert set(a.attrs)==set(b.attrs),(a.name,'attrs')
    for k,v in a.attrs.items():
        if a.name=='/obs' and k=='column-order': assert list(b.attrs[k])==list(v)+['analysis_primary']
        else: assert eq(v,b.attrs[k]),(a.name,k)
    if isinstance(a,h5py.Group):
        extra={'analysis_primary'} if a.name=='/obs' else {'analysis_primary_policy'} if a.name=='/uns' else set()
        assert set(b)==set(a)|extra,(a.name,'children')
        for k in a: preserve(a[k],b[k],records)
    else:
        assert a.shape==b.shape and a.dtype==b.dtype,a.name
        if not a.shape: assert eq(a[()],b[()]),a.name
        else:
            step=max(1,min(100000,2000000//max(1,int(np.prod(a.shape[1:])))))
            for i in range(0,a.shape[0],step): assert eq(a[i:i+step],b[i:i+step]),(a.name,i)
        records.append(dict(path=a.name,shape=list(a.shape),dtype=str(a.dtype)))

def build():
    assert sha(RELEASE/NAME)==OLD_HASH
    assert not HISTORY.exists() and not STAGE.exists()
    HISTORY.mkdir(parents=True);STAGE.mkdir(parents=True)
    prior=json.loads((RELEASE/'release_manifest.json').read_text())
    for name,record in prior['files'].items():
        assert sha(RELEASE/name)==record['sha256']
        shutil.copyfile(RELEASE/name,HISTORY/name)
        assert sha(HISTORY/name)==record['sha256']
    shutil.copyfile(RELEASE/'release_manifest.json',HISTORY/'release_manifest.json')
    docdir=EXP/'archive/pre_patch_release_docs';docdir.mkdir()
    for name in ['README.md','DATA_DICTIONARY_MCFCorpusV3.1.md','DATA_DICTIONARY_MCFCorpusV3.1.zh.md']:
        shutil.copyfile(RELEASE/name,docdir/name)
    shutil.copyfile(HISTORY/NAME,STAGE/NAME)
    d=obs(HISTORY/NAME);flag=select(d)
    assert flag.sum()==2033749
    provenance=dict(patch_id=PATCH,created_at_utc=datetime.now(timezone.utc).isoformat(),rule=POLICY,input_sha256=OLD_HASH,producer_sha256=sha(Path(__file__)),version='MCFCorpusV3.1',missing_run_policy='QC-pass no-Run observations retained individually; no inferred identity')
    with h5py.File(STAGE/NAME,'r+') as f:
        assert 'analysis_primary' not in f['obs'] and 'analysis_primary_policy' not in f['uns']
        write_elem(f['obs'],'analysis_primary',flag,dataset_kwargs={'compression':'gzip'})
        f['obs'].attrs['column-order']=list(f['obs'].attrs['column-order'])+['analysis_primary']
        write_elem(f['uns'],'analysis_primary_policy',json.dumps(provenance,sort_keys=True))
    dump(STAGE/'patch_provenance.json',provenance)
    print('BUILD PASS: backup verified, primary rows',int(flag.sum()),flush=True)

def verify():
    assert sha(HISTORY/NAME)==OLD_HASH
    records=[]
    with h5py.File(HISTORY/NAME,'r') as a,h5py.File(STAGE/NAME,'r') as b:
        preserve(a,b,records)
        flag=read_elem(b['obs/analysis_primary'])
        assert flag.dtype==bool
        policy=json.loads(read_elem(b['uns/analysis_primary_policy']))
        assert policy['input_sha256']==OLD_HASH and policy['producer_sha256']==sha(Path(__file__))
    print('All original datasets preserved:',len(records),flush=True)
    d=obs(HISTORY/NAME); expected=independent_select(d)
    assert np.array_equal(flag,expected) and not (flag & ~d.qc_pass).any()
    chosen=d[flag];known=chosen.Run.ne('')
    assert known.sum()==chosen.loc[known,'Run'].nunique()==1971188
    assert (~known).sum()==62561 and flag.sum()==2033749
    assert np.array_equal(np.flatnonzero(flag & d.Run.ne('')),np.sort(np.load(EXP/'intermediate/depth_priority/baseline_indices.npy')))
    import anndata as ad
    a=ad.read_h5ad(STAGE/NAME,backed='r')
    assert a.shape==(3334541,8114) and len(a.obs.columns)==203
    assert a.uns['corpus_version']=='MCFCorpusV3.1' and np.array_equal(a.obs.analysis_primary.to_numpy(),flag)
    a.file.close()
    stats=dict(primary_rows=int(flag.sum()),unique_nonempty_runs=int(known.sum()),qc_no_run_rows=int((~known).sum()),qc_pass_nonprimary=int((d.qc_pass & ~flag).sum()),by_database={k:int(v) for k,v in chosen.Database.value_counts().items()},primary_with_biome=int(chosen.biome_depth.gt(0).sum()),primary_at_least_l3=int(chosen.biome_depth.ge(3).sum()))
    report=dict(status='PASS',patch_id=PATCH,input_sha256=OLD_HASH,output_sha256=sha(STAGE/NAME),output_bytes=(STAGE/NAME).stat().st_size,unchanged_datasets=records,unchanged_dataset_count=len(records),full_selection_independently_checked=True,anndata_readback='PASS',counts=stats,original_release_verification=dict(path=str((HISTORY/'verification.json').relative_to(ROOT)),sha256=sha(HISTORY/'verification.json')),producer_sha256=sha(Path(__file__)))
    dump(STAGE/'verification.json',report)
    print(json.dumps(stats),flush=True)

def publish():
    report=json.loads((STAGE/'verification.json').read_text())
    assert report['status']=='PASS' and report['producer_sha256']==sha(Path(__file__))
    assert sha(STAGE/NAME)==report['output_sha256']
    assert sha(RELEASE/NAME)==sha(HISTORY/NAME)==OLD_HASH
    assert sha(RELEASE/'release_manifest.json')==sha(HISTORY/'release_manifest.json')
    old=json.loads((HISTORY/'release_manifest.json').read_text()); manifest=json.loads(json.dumps(old))
    manifest['patch_id']=PATCH;manifest['patched_at_utc']=datetime.now(timezone.utc).isoformat()
    manifest['patch_provenance']=json.loads((STAGE/'patch_provenance.json').read_text())
    manifest['previous_release']=dict(path=str(HISTORY.relative_to(ROOT)),h5ad_sha256=OLD_HASH,manifest_sha256=sha(HISTORY/'release_manifest.json'))
    manifest['analysis_primary']=report['counts']
    manifest['scope']='V3.1 biome/QC release plus default representative-row patch; all existing data unchanged by patch.'
    manifest['patch_experiment']=str(EXP.relative_to(ROOT))
    for name in [NAME,'verification.json','patch_provenance.json']:
        p=STAGE/name
        manifest['files'][name]=dict(bytes=p.stat().st_size,sha256=sha(p))
    manifest['verification_report_sha256']=sha(STAGE/'verification.json')
    for name in ['patch_release.py','test_patch_release.py']:
        p=Path(__file__).with_name(name);manifest['source_files_sha256'][str(p.relative_to(ROOT))]=sha(p)
    dump(STAGE/'release_manifest.json',manifest)
    # Stage all replacement files on the target filesystem; the manifest commits last.
    names=[NAME,'verification.json','patch_provenance.json','release_manifest.json']
    for name in names:
        p=RELEASE/(name+'.primary-patch-pending')
        with p.open('xb') as dst,(STAGE/name).open('rb') as src: shutil.copyfileobj(src,dst,8*1024**2)
        assert sha(p)==sha(STAGE/name)
    assert sha(RELEASE/NAME)==OLD_HASH
    for name in names: os.replace(RELEASE/(name+'.primary-patch-pending'),RELEASE/name)
    for name,rec in manifest['files'].items(): assert sha(RELEASE/name)==rec['sha256']
    dump(EXP/'audit/published_patch.json',dict(status='RELEASED_AND_VERIFIED',patch_id=PATCH,sha256=sha(RELEASE/NAME),manifest_sha256=sha(RELEASE/'release_manifest.json'),history_sha256=sha(HISTORY/NAME),counts=report['counts']))
    print('PUBLISHED',PATCH,report['output_sha256'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['build','verify','publish']);args=p.parse_args()
    globals()[args.action]()
