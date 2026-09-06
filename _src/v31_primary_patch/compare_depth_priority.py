"""Second-round read-only sensitivity comparison; no release mutation."""
from pathlib import Path
import hashlib
import json
import time
import h5py
from anndata.io import read_elem
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[4]
EXP = ROOT / 'tmp/20260905_mcfcorpus_v31_primary_run_patch'
OUT = EXP / 'intermediate/depth_priority'
AUD = EXP / 'audit/depth_priority'
CORPUS = ROOT / 'data/MCFCorpusV3.1/MCFCorpusV3.1.gg2.h5ad'

def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(8*1024**2), b''): h.update(b)
    return h.hexdigest()

def main():
    start = time.monotonic()
    manifest = json.loads((EXP / 'audit/phase1_manifest.json').read_text())
    assert sha(CORPUS) == manifest['input_corpus_sha256']
    for name in ['run_selection_comparison.tsv.gz', 'study_source_coverage.tsv', 'study_policy_comparison.tsv']:
        p = EXP / 'intermediate' / name
        assert sha(p) == manifest['files_sha256'][str(p.relative_to(ROOT))]
    assert not OUT.exists() and not AUD.exists(), 'Do not overwrite completed trial'
    OUT.mkdir(parents=True); AUD.mkdir(parents=True)
    with h5py.File(CORPUS, 'r') as f:
        df = pd.DataFrame({k: read_elem(f['obs'][k]) for k in ['Database','Run','qc_pass','n_deep_reads']})
    for c in ['Run','Database']:
        df[c] = df[c].astype(object).where(df[c].notna(), '').astype(str)
    df.Run = df.Run.replace({'nan':'','<NA>':'','None':''})
    df['obs_index'] = np.arange(len(df))
    old = pd.read_csv(EXP / 'intermediate/run_selection_comparison.tsv.gz', sep='\t', usecols=['Run','study_key','baseline_obs_index'], keep_default_na=False)
    base_idx = old.baseline_obs_index.to_numpy()
    assert len(old) == old.Run.nunique() == 1971188
    assert (df.loc[base_idx, 'Run'].to_numpy() == old.Run.to_numpy()).all()
    base_db = df.loc[base_idx, 'Database'].to_numpy()
    base_depth = df.loc[base_idx, 'n_deep_reads'].to_numpy(dtype=float)
    assert (base_depth > 0).all()
    study_map = old.set_index('Run').study_key
    nonrm = base_db != 'RM'
    eligible = df[df.qc_pass & df.Run.isin(old.loc[nonrm & old.study_key.ne(''), 'Run'])].copy()
    eligible['study_key'] = eligible.Run.map(study_map)
    # Recompute matched depth ranking from actual QC copies, not the old chosen rows.
    eligible['relative_depth'] = eligible.n_deep_reads / eligible.groupby('Run').n_deep_reads.transform('max')
    matched = eligible[eligible.groupby('Run').Run.transform('size').gt(1)]
    scores = matched.groupby(['study_key','Database']).relative_depth.median().rename('score')
    coverage = eligible.groupby(['study_key','Database']).Run.nunique().rename('covered_runs').to_frame().join(scores).fillna({'score':0}).reset_index()
    first = pd.read_csv(EXP / 'intermediate/study_source_coverage.tsv', sep='\t')
    check = coverage.merge(first, on=['study_key','Database'], validate='one_to_one')
    assert len(check) == len(first) == len(coverage)
    assert np.allclose(check.score, check.depth_tie_score)
    ranking = coverage.sort_values(['study_key','score','Database'], ascending=[True,False,True])
    leaders = ranking.drop_duplicates('study_key').set_index('study_key').Database
    ranking['primary'] = ranking.Database.eq(ranking.study_key.map(leaders))
    ranking.to_csv(OUT / 'depth_source_ranking.tsv', sep='\t', index=False)
    primary = eligible[eligible.Database.eq(eligible.study_key.map(leaders))].sort_values(['n_deep_reads','obs_index'], ascending=[False,True]).drop_duplicates('Run').set_index('Run').obs_index
    strict_idx = old.Run.map(primary).fillna(old.baseline_obs_index).to_numpy(dtype=np.int64)
    strict_depth = df.loc[strict_idx,'n_deep_reads'].to_numpy(dtype=float)
    strict_ratio = strict_depth / base_depth
    assert (strict_ratio <= 1).all()
    subset = coverage.groupby('study_key').size().loc[lambda s:s.gt(1)].index
    scoped = old.study_key.isin(subset).to_numpy() & nonrm
    assert len(subset) == 155 and scoped.sum() == 95835
    variants = {'baseline':base_idx, 'study_strict':strict_idx}
    for loss in [0,5,10,20]:
        variants[f'study_loss_{loss}pct'] = np.where(strict_ratio >= 1-loss/100, strict_idx, base_idx)
    np.save(OUT / 'qc_no_run_indices.npy', df.index[df.qc_pass & df.Run.eq('')].to_numpy())
    baseline_mix = None
    details, summary = [], []
    for name, idx in variants.items():
        np.save(OUT / f'{name}_indices.npy', idx)
        chosen = df.loc[idx]
        assert chosen.qc_pass.all() and chosen.Run.nunique() == len(idx)
        assert (chosen.Run.to_numpy() == old.Run.to_numpy()).all()
        assert (idx[~nonrm] == base_idx[~nonrm]).all()
        changed = idx != base_idx
        assert not (changed & ~scoped).any()
        depth = chosen.n_deep_reads.to_numpy(dtype=float)
        ratio = depth / base_depth
        if name.startswith('study_loss_'):
            limit = int(name.split('_')[2][:-3])
            assert (ratio >= 1-limit/100).all()
        view = pd.DataFrame({'study_key':old.study_key.to_numpy()[scoped], 'Database':chosen.Database.to_numpy()[scoped], 'changed':changed[scoped], 'depth_ratio':ratio[scoped]})
        per_source = view.groupby(['study_key','Database']).size()
        group = view.groupby('study_key').agg(runs=('Database','size'), sources=('Database','nunique'), changed_runs=('changed','sum'), mean_depth_ratio=('depth_ratio','mean'), min_depth_ratio=('depth_ratio','min'))
        group['minority_runs'] = group.runs - per_source.groupby(level=0).max()
        if baseline_mix is None: baseline_mix = group.copy()
        group['sources_change'] = group.sources - baseline_mix.sources
        group['minority_change'] = group.minority_runs - baseline_mix.minority_runs
        group['variant'] = name
        details.append(group.reset_index())
        r = ratio[changed]
        row = dict(variant=name, changed_runs=int(changed.sum()), shallower_runs=int((ratio<1).sum()), equal_depth_changes=int((changed & (ratio==1)).sum()), mixed_studies=int(group.sources.gt(1).sum()), minority_runs=int(group.minority_runs.sum()), studies_less_sources=int(group.sources_change.lt(0).sum()), studies_more_sources=int(group.sources_change.gt(0).sum()), below_half=int((ratio<.5).sum()), max_depth_loss=float(1-ratio.min()), changed_ratio_median=float(np.median(r)) if len(r) else None, changed_ratio_mean=float(np.mean(r)) if len(r) else None, mean_ratio_scoped=float(ratio[scoped].mean()), changed_studies=int(group.changed_runs.gt(0).sum()))
        summary.append(row)
        if changed.any():
            changes = old.loc[changed,['Run','study_key','baseline_obs_index']].copy()
            changes['candidate_obs_index'] = idx[changed]
            changes['baseline_database'] = base_db[changed]
            changes['candidate_database'] = chosen.Database.to_numpy()[changed]
            changes['baseline_depth'] = base_depth[changed]
            changes['candidate_depth'] = depth[changed]
            changes['depth_ratio'] = ratio[changed]
            changes.to_csv(OUT / f'{name}_changes.tsv.gz', sep='\t', index=False)
        print(row, flush=True)
    pd.DataFrame(summary).to_csv(OUT / 'summary.tsv', sep='\t', index=False)
    pd.concat(details).to_csv(OUT / 'per_study.tsv', sep='\t', index=False)
    assert sha(CORPUS) == manifest['input_corpus_sha256']
    result = dict(status='COMPUTED_NOT_RELEASED', corpus_sha256=sha(CORPUS), known_runs=len(old), qc_no_run=62561, units=2033749, compared_studies=len(subset), compared_runs=int(scoped.sum()), unclosed_project_runs=int(old.study_key.eq('').sum()), variants=summary, elapsed_seconds=time.monotonic()-start)
    (AUD / 'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    files = [Path(__file__), EXP/'PLAN.md', *OUT.iterdir(), AUD/'summary.json']
    (AUD / 'manifest.json').write_text(json.dumps({str(p.relative_to(ROOT)):sha(p) for p in files}, indent=2)+'\n')

if __name__ == '__main__': main()
