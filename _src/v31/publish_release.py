#!/usr/bin/env python3
"""Publish an already verified V3.1 without replacing any existing release."""
import json
import shutil
from datetime import datetime, timezone

from build_release import EXP, ROOT, PRE, V3, EXPECTED_HASHES, VOCAB, MANIFEST, sha, save


def copy_exclusive(source, dest):
    with source.open('rb') as src, dest.open('xb') as dst:
        shutil.copyfileobj(src, dst, length=8 * 1024**2)
    assert sha(source) == sha(dest), str(dest)


if __name__ == '__main__':
    report_path = EXP / 'audit/verification.json'
    report = json.loads(report_path.read_text())
    assert report['status'] == 'PASS'
    source = ROOT / report['output']
    assert sha(source) == report['output_sha256']
    assert sha(V3) == EXPECTED_HASHES[V3]
    assert sha(ROOT / report['provenance']['producer']) == report['provenance']['producer_sha256']
    target_dir = ROOT / 'data/MCFCorpusV3.1'
    target_dir.mkdir(exist_ok=False)
    output = target_dir / 'MCFCorpusV3.1.gg2.h5ad'
    copy_exclusive(source, output)
    files = {
        'verification.json': report_path,
        'gold_biosample_vocabulary_20260901.tsv': VOCAB,
        'biome_preparation_manifest.json': MANIFEST,
    }
    for name, path in files.items():
        copy_exclusive(path, target_dir / name)
    payload = {
        'MCFCorpusV3.1.gg2.h5ad': {'bytes': output.stat().st_size, 'sha256': sha(output)},
        **{name: {'bytes': p.stat().st_size, 'sha256': sha(p)} for name, p in files.items()},
    }
    code_dir = ROOT / 'data/MCFCorpus_incoming/_src/v31'
    code = {str(p.relative_to(ROOT)): sha(p) for p in sorted(code_dir.glob('*.py'))}
    result = {
        'status': 'RELEASED', 'version': 'MCFCorpusV3.1',
        'released_at_utc': datetime.now(timezone.utc).isoformat(),
        'canonical_path': str(output.relative_to(ROOT)),
        'shape': report['shape'], 'counts_by_database': report['counts_by_database'],
        'files': payload, 'source_files_sha256': code,
        'provenance': report['provenance'],
        'experiment': 'tmp/20260905_mcfcorpus_v31_build/',
        'verification_report_sha256': sha(report_path),
        'scope': 'biome/QC metadata-only update, no row removal; V3 untouched',
        'split_task_status': 'No V3.1 split/task generated; V3 indices are historical, not V3.1 eligibility.',
        'evidence_location': str(PRE.relative_to(ROOT)),
    }
    save(target_dir / 'release_manifest.json', result)
    save(EXP / 'audit/publication.json', result)
    print('RELEASED', output, result['files']['MCFCorpusV3.1.gg2.h5ad'], flush=True)
