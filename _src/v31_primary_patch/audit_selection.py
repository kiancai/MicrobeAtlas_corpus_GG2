#!/usr/bin/env python3
"""Read-only, full-universe V3.1 representative-row policy comparison.

Writes exploratory outputs only to the new experiment. Never patches a corpus.
"""
from pathlib import Path
from collections import defaultdict
import csv
import gzip
import hashlib
import json
import re
import time

import h5py
from anndata.io import read_elem
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[4]
EXP = ROOT / 'tmp/20260905_mcfcorpus_v31_primary_run_patch'
CORPUS = ROOT / 'data/MCFCorpusV3.1/MCFCorpusV3.1.gg2.h5ad'
ENA = ROOT / 'tmp/20260827_v3_mgnify_gold_standardization/intermediate/ena_read_run_mgnify_study_gold_catalog.tsv.gz'
EXPECTED = 'b0612b46abd0d734ed501c4c312ba76c18f8e8c61034d33cae66af8ee7ba97e3'
TOKEN = re.compile(r'(?<![A-Za-z0-9])(?:PRJ[A-Z]{2}\d+|[SED]RP\d+)(?![A-Za-z0-9])')


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(8*1024**2), b''): h.update(b)
    return h.hexdigest()


def table(name, df):
    p = EXP / 'intermediate' / name
    p.parent.mkdir(parents=True, exist_ok=True)
    assert not p.exists(), p
    df.to_csv(p, sep='\t', index=False)


def quantiles(s):
    s = pd.Series(s).dropna()
    return {str(q): float(s.quantile(q)) for q in [0, .01, .05, .25, .5, .75, .95, 1]} if len(s) else {}


def main():
    started = time.monotonic()
    assert sha(CORPUS) == EXPECTED
    hashes = {str(p.relative_to(ROOT)): sha(p) for p in [ENA, Path(__file__)]}
    aliases = defaultdict(set)
    with gzip.open(ENA, 'rt') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            primary = row['ena_primary_study_accession']
            secondary = row['ena_secondary_study_accession']
            if re.fullmatch(r'PRJ[A-Z]{2}\d+', primary) and re.fullmatch(r'[SED]RP\d+', secondary):
                aliases[secondary].add(primary)
    alias = {k: next(iter(v)) for k,v in aliases.items() if len(v) == 1}
    table('ena_project_aliases.tsv', pd.DataFrame([
        {'secondary': k, 'primary_candidates': json.dumps(sorted(v)), 'unambiguous': len(v)==1}
        for k,v in sorted(aliases.items())]))
    print('Official direct aliases', len(alias), flush=True)
    columns = ['Database','Run','Project_ID','QIITA_Study_ID','qc_pass','n_deep_reads','biome_depth','gold_l1','gold_l2','gold_l3','gold_l4']
    with h5py.File(CORPUS, 'r') as f:
        data = {k: read_elem(f['obs'][k]) for k in columns}
        data['obs_name'] = read_elem(f['obs'][f['obs'].attrs['_index']])
    df = pd.DataFrame(data)
    df['obs_index'] = np.arange(len(df), dtype=np.int64)
    for k in ['Database','Run','Project_ID','QIITA_Study_ID','gold_l1','gold_l2','gold_l3','gold_l4']:
        df[k] = df[k].astype(object).where(pd.notna(df[k]), '').astype(str)
    df.Run = df.Run.replace({'nan':'','<NA>':'','None':''})
    assert len(df) == 3334541 and int(df.qc_pass.sum()) == 2323396
    project_map, project_status = {}, {}
    for raw in df.Project_ID.unique():
        tokens = TOKEN.findall(raw)
        rest = TOKEN.sub('', raw).strip(' \t\r\n,;|/[]{}()\'"')
        if not raw or raw in ['nan','None','<NA>']:
            value, state = '', 'missing'
        elif not tokens or rest:
            value, state = '', 'unparsed'
        elif any(t in aliases and len(aliases[t]) > 1 for t in tokens):
            value, state = '', 'official_alias_ambiguous'
        else:
            mapped = sorted(set(alias.get(t, t) for t in tokens))
            if len(mapped) == 1:
                value = mapped[0]
                state = 'primary_or_official_alias' if value.startswith('PRJ') else 'literal_secondary_unmapped'
            else:
                value, state = '', 'multiple_unreconciled_projects'
        project_map[raw], project_status[raw] = value, state
    df['project_key'] = df.Project_ID.map(project_map)
    df['project_state'] = df.Project_ID.map(project_status)
    pstat = df.groupby(['Project_ID','project_state'], dropna=False).agg(rows=('Run','size'),qc_rows=('qc_pass','sum')).reset_index()
    pstat['project_key'] = pstat.Project_ID.map(project_map)
    table('project_token_status.tsv', pstat)
    # Run identity uses all original observations, not just the chosen profile or QC survivors.
    edges = df.loc[df.Run.ne('') & df.project_key.ne(''), ['Run','project_key']].drop_duplicates()
    nkeys = edges.groupby('Run').size()
    single = edges[edges.Run.map(nkeys).eq(1)].set_index('Run').project_key
    bad = set(df.loc[df.Run.ne('') & df.project_state.isin(['unparsed','official_alias_ambiguous','multiple_unreconciled_projects']), 'Run'])
    valid = single[~single.index.isin(bad)]
    qc = df[df.qc_pass].copy()
    nonempty = qc[qc.Run.ne('')].copy()
    norun = qc[qc.Run.eq('')].copy()
    nonempty['is_rm'] = nonempty.Database.eq('RM').astype('int8')
    nonempty = nonempty.sort_values(['Run','is_rm','n_deep_reads','obs_index'], ascending=[True,False,False,True], kind='stable')
    baseline = nonempty.drop_duplicates('Run').set_index('Run', drop=False).copy()
    baseline['study_key'] = baseline.Run.map(valid).fillna('')
    baseline['identity_status'] = np.where(baseline.study_key.ne(''), 'resolved_explicit_project', 'no_usable_project')
    baseline.loc[baseline.Run.map(nkeys).gt(1), 'identity_status'] = 'multiple_unreconciled_projects'
    baseline.loc[baseline.Run.isin(bad), 'identity_status'] = 'unparsed_or_multi_project_record'
    # If every copy lacks a project, a unique explicitly recorded Qiita study is only a source-scoped fallback.
    qe = df.loc[df.Run.ne('') & df.QIITA_Study_ID.ne(''), ['Run','QIITA_Study_ID']].drop_duplicates()
    qn = qe.groupby('Run').size()
    qsingle = qe[qe.Run.map(qn).eq(1)].set_index('Run').QIITA_Study_ID
    qfill = baseline.identity_status.eq('no_usable_project') & baseline.Run.isin(qsingle.index) & ~baseline.Run.isin(nkeys.index)
    baseline.loc[qfill, 'study_key'] = 'QIITA:' + baseline.loc[qfill, 'Run'].map(qsingle)
    baseline.loc[qfill, 'identity_status'] = 'source_scoped_qiita_study'
    rm_runs = set(baseline.loc[baseline.Database.eq('RM'),'Run'])
    nonempty['study_key'] = nonempty.Run.map(baseline.study_key)
    eligible = nonempty[~nonempty.Run.isin(rm_runs) & nonempty.study_key.ne('')].copy()
    # Match depth comparisons within Run; all unique-run profiles have ratio 1 but are excluded from depth tie ranking.
    copies = eligible.groupby('Run').size()
    eligible['relative_depth'] = eligible.n_deep_reads / eligible.groupby('Run').n_deep_reads.transform('max')
    matched = eligible[eligible.Run.map(copies).gt(1)]
    coverage = eligible.groupby(['study_key','Database']).agg(covered_runs=('Run','nunique')).reset_index()
    depth = matched.groupby(['study_key','Database']).relative_depth.median().rename('matched_depth_median').reset_index()
    coverage = coverage.merge(depth, how='left')
    union = eligible.groupby('study_key').Run.nunique()
    coverage['study_runs'] = coverage.study_key.map(union)
    coverage['coverage_fraction'] = coverage.covered_runs / coverage.study_runs
    coverage['depth_tie_score'] = coverage.matched_depth_median.fillna(0.)
    ranking = coverage.sort_values(['study_key','covered_runs','depth_tie_score','Database'],ascending=[True,False,False,True],kind='stable')
    leaders = ranking.drop_duplicates('study_key').set_index('study_key')
    primary = leaders.Database
    candidate = baseline.copy()
    candidate['selection_reason'] = np.where(candidate.Database.eq('RM'),'rm_priority','unresolved_study_depth_fallback')
    # Start with baseline for coverage preservation; replace only where the study-primary source has an eligible copy.
    from_primary = eligible[eligible.Database.eq(eligible.study_key.map(primary))].drop_duplicates('Run').set_index('Run')
    known_nonrm = ~candidate.Run.isin(rm_runs) & candidate.study_key.ne('')
    candidate.loc[known_nonrm,'selection_reason'] = 'primary_missing_run_depth_fallback'
    for name in ['obs_index','obs_name','Database','n_deep_reads','biome_depth','gold_l1','gold_l2','gold_l3','gold_l4']:
        candidate.loc[from_primary.index, name] = from_primary[name]
    candidate.loc[from_primary.index,'selection_reason'] = 'study_primary'
    candidate['primary_database'] = candidate.study_key.map(primary).fillna('')
    changed = candidate.obs_index.ne(baseline.obs_index)
    changes = candidate.loc[changed,['Run','study_key','obs_index','obs_name','Database','n_deep_reads','selection_reason']].copy()
    changes['baseline_obs_index'] = baseline.loc[changes.index,'obs_index']
    changes['baseline_database'] = baseline.loc[changes.index,'Database']
    changes['baseline_depth'] = baseline.loc[changes.index,'n_deep_reads']
    changes['depth_ratio'] = changes.n_deep_reads / changes.baseline_depth
    table('changed_run_choices.tsv.gz', changes.reset_index(drop=True))
    table('study_source_coverage.tsv',coverage)
    table('run_selection_comparison.tsv.gz',candidate[['Run','study_key','identity_status','obs_index','Database','n_deep_reads','selection_reason','primary_database']].assign(baseline_obs_index=baseline.obs_index).reset_index(drop=True))
    unresolved = nonempty[nonempty.Run.isin(candidate.loc[candidate.study_key.eq(''),'Run'])]
    table('unresolved_run_project_evidence.tsv.gz',unresolved[['Run','Database','Project_ID','project_state','project_key']])
    support = leaders[['Database','covered_runs','study_runs','coverage_fraction']].copy()
    support['n_sources'] = coverage.groupby('study_key').size()
    support['coverage_tie_sources'] = coverage[coverage.covered_runs.eq(coverage.study_key.map(leaders.covered_runs))].groupby('study_key').size()
    support['baseline_sources'] = baseline.loc[baseline.Run.isin(eligible.Run)].groupby('study_key').Database.nunique()
    support['candidate_sources'] = candidate.loc[candidate.Run.isin(eligible.Run)].groupby('study_key').Database.nunique()
    support['supplement_runs'] = support.study_runs - support.covered_runs
    support['changed_runs'] = changes.groupby('study_key').size().reindex(support.index,fill_value=0)
    table('study_policy_comparison.tsv',support.reset_index())
    # Include RM exceptions in the all-known-study source-mixing count separately.
    allstudy = baseline[baseline.study_key.ne('')].groupby('study_key').agg(runs=('Run','size'),baseline_sources=('Database','nunique'))
    allstudy['candidate_sources'] = candidate[candidate.study_key.ne('')].groupby('study_key').Database.nunique()
    allstudy['rm_runs'] = baseline[baseline.Database.eq('RM')].groupby('study_key').size().reindex(allstudy.index,fill_value=0)
    table('all_study_mixing_with_rm.tsv',allstudy.reset_index())
    # Independent baseline check by numeric RM priority + max depth, not by rerunning the sort.
    selected = df.set_index('obs_index').loc[candidate.obs_index]
    assert selected.qc_pass.all() and selected.Run.is_unique
    assert set(candidate.Run) == set(nonempty.Run)
    assert len(baseline)==len(candidate)
    rm_observations = set(qc.loc[qc.Database.eq('RM'),'obs_index'])
    assert rm_observations <= set(candidate.obs_index)
    no_rm = nonempty[~nonempty.Run.isin(rm_runs)]
    maxima = no_rm.groupby('Run').n_deep_reads.max()
    assert baseline.loc[maxima.index,'n_deep_reads'].equals(maxima)
    assert (candidate.n_deep_reads <= baseline.n_deep_reads).all()
    assert len(norun)==int(qc.Run.eq('').sum())
    assert int(support.supplement_runs.sum()) == int(candidate.selection_reason.eq('primary_missing_run_depth_fallback').sum())
    assert len(from_primary) == int(support.covered_runs.sum())
    stats = {}
    for name,chosen in [('run_depth_baseline',baseline),('study_aware_trial',candidate)]:
        idx = np.concatenate([chosen.obs_index.to_numpy(),norun.obs_index.to_numpy()])
        rows = df.set_index('obs_index').loc[idx]
        assert len(idx)==len(set(idx)) and rows.qc_pass.all()
        stats[name] = {'selected_rows':len(idx),'selected_unique_nonempty_runs':len(chosen),'no_run_rows':len(norun),
                       'by_database':{k:int(v) for k,v in rows.Database.value_counts().items()},
                       'any_biome':int(rows.biome_depth.gt(0).sum()),'at_least_L3':int(rows.biome_depth.ge(3).sum())}
        table(f'{name}_biome_counts.tsv',rows.groupby(['gold_l1','gold_l2','gold_l3','gold_l4'],dropna=False).size().rename('rows').reset_index())
    multiple = support.n_sources.gt(1)
    duplicate_sizes = nonempty.groupby('Run').size()
    summary = {
        'status':'VERIFIED_EXPLORATORY_NOT_PATCHED','input_sha256':EXPECTED,'supporting_inputs_sha256':hashes,
        'alias_pairs_unambiguous':len(alias),'alias_secondary_ambiguous':sum(len(v)>1 for v in aliases.values()),
        'qc_pass_rows':len(qc),'qc_rows_with_run':len(nonempty),'qc_rows_without_run':len(norun),
        'qc_unique_nonempty_runs':len(baseline),'qc_duplicate_run_groups':int(duplicate_sizes.gt(1).sum()),
        'qc_redundant_rows':len(nonempty)-len(baseline),'rm_priority_runs':len(rm_runs),
        'run_identity_status':{k:int(v) for k,v in baseline.identity_status.value_counts().items()},
        'nonrm_runs_study_resolved':len(eligible.Run.unique()),'nonrm_runs_study_unresolved':int((~candidate.Run.isin(rm_runs)&candidate.study_key.eq('')).sum()),
        'nonrm_study_groups':len(support),'nonrm_multisource_study_groups':int(multiple.sum()),
        'multisource_study_runs':int(support.loc[multiple,'study_runs'].sum()),
        'multisource_studies_single_source_full_coverage':int((multiple & support.coverage_fraction.eq(1)).sum()),
        'multisource_studies_needing_supplement':int((multiple & support.coverage_fraction.lt(1)).sum()),
        'supplement_runs':int(support.supplement_runs.sum()),
        'coverage_tie_studies':int(support.coverage_tie_sources.gt(1).sum()),
        'mixed_nonrm_studies_baseline':int(support.baseline_sources.gt(1).sum()),
        'mixed_nonrm_studies_trial':int(support.candidate_sources.gt(1).sum()),
        'mixed_studies_with_rm_baseline':int(allstudy.baseline_sources.gt(1).sum()),
        'mixed_studies_with_rm_trial':int(allstudy.candidate_sources.gt(1).sum()),
        'changed_runs':int(changed.sum()),'changed_depth_ratio':quantiles(changes.depth_ratio),
        'changed_depth_below_half':int(changes.depth_ratio.lt(.5).sum()),'changed_depth_below_tenth':int(changes.depth_ratio.lt(.1).sum()),
        'selection_reasons':{k:int(v) for k,v in candidate.selection_reason.value_counts().items()},
        'policies':stats,'elapsed_seconds':time.monotonic()-started,
        'limitations':['local ENA direct aliases only; unclosed project relationships are fallback, not claimed resolved',
                       'strict coverage-first trial; not final production rule',
                       'no-Run retained individually, not proven independent samples',
                       'primary-study grouping is not a clinical cohort or split component'],
    }
    assert sha(CORPUS)==EXPECTED
    assert all(sha(ROOT/k)==v for k,v in hashes.items())
    out = EXP/'audit';out.mkdir(exist_ok=True)
    (out/'summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n')
    payload={str(p.relative_to(EXP)):sha(p) for p in (EXP/'intermediate').iterdir()}
    (out/'output_manifest.json').write_text(json.dumps({'status':'PASS','files_sha256':payload,'summary_sha256':sha(out/'summary.json')},indent=2)+'\n')
    print(json.dumps(summary,indent=2,ensure_ascii=False),flush=True)


if __name__=='__main__': main()
