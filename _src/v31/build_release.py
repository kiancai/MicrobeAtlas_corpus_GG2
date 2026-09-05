#!/usr/bin/env python3
"""Hash-bound, metadata-only V3.1 build and exhaustive read-only verification.

Run build then verify in separate processes. No publication or input mutation here.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
from anndata.io import read_elem, write_elem
import h5py
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[4]
PRE = ROOT / 'tmp/20260827_v3_mgnify_gold_standardization'
EXP = ROOT / 'tmp/20260905_mcfcorpus_v31_build'
V3 = ROOT / 'data/MCFCorpusV3/MCFCorpusV3.gg2.h5ad'
CAND = PRE / 'intermediate/p5b/v31_integration/v3_observation_biome_qc.tsv.gz'
MANIFEST = PRE / 'audit/p5b/v31_integration/final_manifest.json'
VOCAB = PRE / 'intermediate/gold_5level_path_vocabulary_20260901.tsv'
EXPECTED_HASHES = {
    V3: 'f10ffe14d2e4384a738e84a9d26694a8e46055520cc30fd5e9fb1d8415d6d78c',
    CAND: '57905c0aa7d5a261fd48e980818fcec823391aed59a0f3beef71cedf107d5fbc',
    MANIFEST: '679fa2dec2e460e1dbd1a59b6e84768cd215ddd8b35689845ac91b259f9dc3ac',
    VOCAB: '378a911bfc8e71423411718ef2bc8b1dbc3d93df5e89cdcd0be8f2aacb6cd81c',
}
BOOL_FIELDS = {'cross_source_disagreement', 'v3_qc_pass', 'candidate_qc_pass', 'new_qc_false'}
RENAMES = {'v3_qc_pass': 'qc_pass_v3', 'candidate_qc_pass': 'qc_pass'}
IDENTITY = ['v3_obs_index', 'v3_obs_name', 'Database', 'Run']
QC_DEF = ('V3.1: qc_pass_v3 AND NOT audited exclusion; qc_pass_v3 retains V3 '
          'technical thresholds. Intrinsic technical controls/defined communities/isolates '
          'propagate by exact Run; RM Oral/Cough swab/missing site are RM-only; '
          'without Run apply only to the evidence observation. Missing biome is not QC failure.')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(8 * 1024**2), b''):
            h.update(b)
    return h.hexdigest()


def save(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
        f.write('\n')


def inputs():
    for p, expected in EXPECTED_HASHES.items():
        assert sha(p) == expected, f'Input hash changed: {p}'
    m = json.loads(MANIFEST.read_text())
    assert m['status'] == 'PASS'
    for relative, record in m['files'].items():
        p = PRE / relative
        assert p.stat().st_size == record['bytes'], str(p)
        assert sha(p) == record['sha256'], str(p)
    print('All fixed inputs and preparation package hashes PASS', flush=True)
    return m


def candidate():
    c = pd.read_csv(CAND, sep='\t', dtype=str, keep_default_na=False)
    assert len(c) == 3334541
    assert np.array_equal(c.v3_obs_index.astype('int64'), np.arange(len(c)))
    for name in BOOL_FIELDS:
        assert set(c[name]) <= {'True', 'False'}, name
    return c


def provenance():
    return {
        'version': 'MCFCorpusV3.1',
        'created_at_utc': datetime.now(timezone.utc).isoformat(),
        'inputs': {str(p.relative_to(ROOT)): h for p, h in EXPECTED_HASHES.items()},
        'producer': str(Path(__file__).relative_to(ROOT)),
        'producer_sha256': sha(__file__),
        'producer_repo_base_commit': '2eeaf01bd19c28c14c541bbec9310bde9b012a98',
        'producer_revision': 'new uncommitted source, bound by SHA256',
        'environment': {p: importlib.metadata.version(p) for p in ['anndata', 'h5py', 'numpy', 'pandas']},
        'GOLD_snapshot': '2026-09-01; official release number not asserted',
        'contract': 'Preserve all V3 rows, order, matrices and old metadata except qc_pass; human max L4, other max L3, L5 empty.',
    }


def arr(f, name):
    value = read_elem(f['obs'][name])
    return np.asarray(value)


def build(output):
    assert not output.exists(), 'Never overwrite a build'
    inputs()
    c = candidate()
    new = [RENAMES.get(k, k) for k in c.columns if k not in IDENTITY and k != 'candidate_qc_pass']
    with h5py.File(V3, 'r') as f:
        assert not set(new) & set(f['obs']), 'Field collision'
        assert not {'corpus_version', 'v31_provenance', 'qc_def_v3'} & set(f['uns'])
        assert np.array_equal(arr(f, 'sample_id'), c.v3_obs_name.to_numpy())
        for k in ['Database', 'Run']:
            assert np.array_equal(pd.Series(arr(f, k)).fillna('').to_numpy(), c[k].to_numpy()), k
        assert np.array_equal(arr(f, 'qc_pass'), c.v3_qc_pass.eq('True').to_numpy())
    output.parent.mkdir(parents=True, exist_ok=True)
    # Exclusive reservation prevents accidentally replacing an existing output.
    with output.open('xb') as dst, V3.open('rb') as src:
        shutil.copyfileobj(src, dst, length=8 * 1024**2)
    with h5py.File(output, 'r+') as f:
        for name in c.columns:
            if name in IDENTITY:
                continue
            if name in BOOL_FIELDS:
                value = c[name].eq('True').to_numpy()
            elif name == 'biome_depth':
                value = c[name].astype('int8').to_numpy()
            else:
                value = pd.Categorical(c[name])
            write_elem(f['obs'], RENAMES.get(name, name), value, dataset_kwargs={'compression': 'gzip'})
        f['obs'].attrs['column-order'] = list(f['obs'].attrs['column-order']) + new
        write_elem(f['uns'], 'qc_def_v3', read_elem(f['uns/qc_def']))
        write_elem(f['uns'], 'qc_def', QC_DEF)
        write_elem(f['uns'], 'corpus_version', 'MCFCorpusV3.1')
        write_elem(f['uns'], 'v31_provenance', json.dumps(provenance(), ensure_ascii=False, sort_keys=True))
        f.flush()
    print(f'BUILD COMPLETE {output}', flush=True)


def equal(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if a.dtype.kind in 'fc' and b.dtype.kind in 'fc':
        return np.array_equal(a, b, equal_nan=True)
    return np.array_equal(a, b)


def compare_original(a, b, records):
    """Compare original HDF5 hierarchy, attrs, dtypes and every dataset value.

    Intentional changes are narrowly enumerated, not a broad obs/uns skip.
    """
    for k, v in a.attrs.items():
        if a.name == '/obs' and k == 'column-order':
            assert list(b.attrs[k])[:len(v)] == list(v)
        else:
            assert k in b.attrs and equal(v, b.attrs[k]), (a.name, k)
    assert set(a.attrs) == set(b.attrs), a.name
    if isinstance(a, h5py.Group):
        if a.name not in ['/obs', '/uns']:
            assert set(a) == set(b), a.name
        for k in a:
            assert k in b, (a.name, k)
            if a[k].name in ['/obs/qc_pass', '/uns/qc_def']:
                continue
            compare_original(a[k], b[k], records)
    else:
        assert a.shape == b.shape and a.dtype == b.dtype, a.name
        if not a.shape:
            assert equal(a[()], b[()]), a.name
        else:
            width = max(1, int(np.prod(a.shape[1:])))
            step = max(1, min(100000, 2000000 // width))
            for start in range(0, a.shape[0], step):
                assert equal(a[start:start+step], b[start:start+step]), (a.name, start)
        records.append({'path': a.name, 'shape': list(a.shape), 'dtype': str(a.dtype)})


def verify(output, report):
    m = inputs()
    c = candidate()
    records = []
    with h5py.File(V3, 'r') as a, h5py.File(output, 'r') as b:
        compare_original(a, b, records)
        new = {RENAMES.get(k, k) for k in c if k not in IDENTITY} - {'qc_pass'}
        assert set(b['obs']) == set(a['obs']) | new
        assert set(b['uns']) == set(a['uns']) | {'corpus_version', 'v31_provenance', 'qc_def_v3'}
        assert equal(arr(b, 'sample_id'), c.v3_obs_name)
        for k in ['Database', 'Run']:
            assert equal(pd.Series(arr(b, k)).fillna(''), c[k]), k
        assert equal(arr(b, 'qc_pass_v3'), arr(a, 'qc_pass'))
        for name in c:
            if name in IDENTITY:
                continue
            expected = c[name]
            if name in BOOL_FIELDS:
                expected = expected.eq('True')
            elif name == 'biome_depth':
                expected = expected.astype('int8')
            assert equal(arr(b, RENAMES.get(name, name)), expected), name
        assert read_elem(b['uns/corpus_version']) == 'MCFCorpusV3.1'
        assert read_elem(b['uns/qc_def']) == QC_DEF
        assert equal(read_elem(b['uns/qc_def_v3']), read_elem(a['uns/qc_def']))
        bound = json.loads(read_elem(b['uns/v31_provenance']))
        assert bound['producer_sha256'] == sha(__file__)
        assert bound['inputs'] == {str(p.relative_to(ROOT)): h for p, h in EXPECTED_HASHES.items()}
    print(f'All {len(records)} unchanged datasets and all label/QC rows PASS', flush=True)
    # Independent vocabulary/prefix and QC invariants, not just writer round-trip.
    vocab = pd.read_csv(VOCAB, sep='\t', dtype=str, keep_default_na=False)
    fields = ['ECOSYSTEM', 'ECOSYSTEM CATEGORY', 'ECOSYSTEM TYPE', 'ECOSYSTEM SUBTYPE', 'SPECIFIC ECOSYSTEM']
    valid = {tuple(row[:d]) for row in vocab[fields].itertuples(index=False, name=None) for d in range(1, 5)}
    gl = [f'gold_l{i}' for i in range(1, 6)]
    for row in c[gl+['biome_depth']].drop_duplicates().itertuples(index=False, name=None):
        d = int(row[5]); path = tuple(row[:d])
        assert 0 <= d <= 4 and all(not x for x in row[d:5]), row
        assert all(x and x != 'Unclassified' for x in path), row
        assert not d or path in valid, row
        assert d <= (4 if row[1] == 'Mammals: Human' else 3), row
    old, qc, delta = (c[k].eq('True') for k in ['v3_qc_pass', 'candidate_qc_pass', 'new_qc_false'])
    assert not (qc & ~old).any()
    assert equal(delta, old & ~qc)
    label = c.biome_depth.astype(int).gt(0)
    actual = {}
    for db in ['ALL', 'MA', 'RM', 'MBM', 'PRIME', 'Qiita']:
        s = pd.Series(True, index=c.index) if db == 'ALL' else c.Database.eq(db)
        count = {'observations': int(s.sum()), 'label_bearing': int((s & label).sum()),
                 'old_qc_false': int((s & ~old).sum()), 'new_qc_false': int((s & delta).sum()),
                 'candidate_qc_false': int((s & ~qc).sum()), 'candidate_qc_pass': int((s & qc).sum()),
                 'label_bearing_candidate_qc_pass': int((s & qc & label).sum()),
                 'missing_Run': int((s & c.Run.eq('')).sum()),
                 'new_qc_false_no_Run': int((s & delta & c.Run.eq('')).sum())}
        assert count == m['counts']['by_database'][db], (db, count)
        actual[db] = count
    del c
    obj = ad.read_h5ad(output, backed='r')
    assert obj.shape == (3334541, 8114)
    assert len(obj.obs.columns) == len(set(obj.obs.columns))
    assert obj.obs.qc_pass.dtype == bool and obj.obs.qc_pass_v3.dtype == bool
    obs_columns = list(obj.obs.columns)
    obj.file.close()
    out_hash = sha(output)
    assert sha(V3) == EXPECTED_HASHES[V3], 'Original V3 changed'
    save(report, {'status': 'PASS', 'verified_at_utc': datetime.now(timezone.utc).isoformat(),
                  'output': str(output.relative_to(ROOT)), 'output_sha256': out_hash,
                  'output_bytes': output.stat().st_size, 'shape': [3334541, 8114],
                  'counts_by_database': actual, 'original_datasets_verified': records,
                  'obs_columns': obs_columns, 'provenance': bound,
                  'checks': ['all unchanged HDF5 values/dtypes/attributes', 'all obs identities and candidate fields',
                             'GOLD legal prefixes and depth', 'original QC retained; exact QC delta',
                             'all database counts', 'AnnData backed read', 'V3 unchanged full SHA256']})
    print(f'VERIFICATION PASS {out_hash}', flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('mode', choices=['build', 'verify'])
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--report', type=Path)
    args = p.parse_args()
    if args.mode == 'build':
        build(args.output.resolve())
    else:
        if args.report is None:
            p.error('--report required for verify')
        verify(args.output.resolve(), args.report.resolve())
