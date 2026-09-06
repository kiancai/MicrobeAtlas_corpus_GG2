#!/usr/bin/env python3
"""Explain unresolved project relationships and large depth changes, without new choices."""
from pathlib import Path
from collections import Counter
import json
import pandas as pd

ROOT=Path(__file__).resolve().parents[4]
EXP=ROOT/'tmp/20260905_mcfcorpus_v31_primary_run_patch'
u=pd.read_csv(EXP/'intermediate/unresolved_run_project_evidence.tsv.gz',sep='\t',dtype=str,keep_default_na=False)
pairs=u[['Run','project_key']].drop_duplicates()
types=Counter(); signatures=Counter()
for run,g in pairs.groupby('Run',sort=False):
    keys=sorted(k for k in g.project_key if k)
    prj=[k for k in keys if k.startswith('PRJ')]
    sec=[k for k in keys if not k.startswith('PRJ')]
    if len(prj)==1 and len(sec)>=1: kind='one_primary_plus_unmapped_secondary'
    elif len(prj)>1: kind='multiple_primary_projects'
    elif len(sec)>1: kind='multiple_secondary_projects'
    else: kind='conflict_in_other_nonQC_copy_or_incomplete_key'
    types[kind]+=1
    signatures[(kind,' | '.join(keys))]+=1
sig=pd.DataFrame([{'type':k[0],'projects':k[1],'runs':v} for k,v in signatures.items()]).sort_values('runs',ascending=False)
sig.to_csv(EXP/'intermediate/unresolved_project_signatures.tsv',sep='\t',index=False)
c=pd.read_csv(EXP/'intermediate/changed_run_choices.tsv.gz',sep='\t')
studies=pd.read_csv(EXP/'intermediate/study_policy_comparison.tsv',sep='\t')
g=c.groupby('study_key').agg(changed_runs=('Run','size'),median_depth_ratio=('depth_ratio','median'),minimum_depth_ratio=('depth_ratio','min'))
g['below_half']=c[c.depth_ratio.lt(.5)].groupby('study_key').size().reindex(g.index,fill_value=0)
g['below_tenth']=c[c.depth_ratio.lt(.1)].groupby('study_key').size().reindex(g.index,fill_value=0)
g=g.reset_index().merge(studies.drop(columns=['changed_runs']),on='study_key').sort_values('changed_runs',ascending=False)
g.to_csv(EXP/'intermediate/depth_change_by_study.tsv',sep='\t',index=False)
result={'unresolved_qc_evidence_types':dict(types),'top_unresolved_signatures':sig.head(10).to_dict('records'),
        'top_changed_studies':g.head(10).to_dict('records'),
        'extreme_loss_studies':g[g.below_tenth.gt(0)].sort_values('below_tenth',ascending=False).head(10).to_dict('records')}
(EXP/'audit/diagnostics.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(result,ensure_ascii=False,indent=2))
