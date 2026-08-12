#!/usr/bin/env python
"""Qiita Layer1 步04(B): 合并各 shard taxonomy → 全局 md5→taxonomy。校验 == asv_global.fa。"""
import glob, numpy as np, pandas as pd

ROOT = "/hpcdisk1/limk_group/caiqy/project/260428_greengene2"
WORK = f"{ROOT}/rawdata/Qiita/work"
TAX_DIR = f"{WORK}/tax"
FASTA = f"{WORK}/asv_global.fa"
OUT = f"{WORK}/taxonomy_gg2.global.tsv"

shard_files = sorted(glob.glob(f"{TAX_DIR}/shard_*.taxonomy.tsv"))
print(f"shard taxonomy 片: {len(shard_files)}")
parts = []
for f in shard_files:
    df = pd.read_csv(f, sep='\t', dtype=str)
    if len(df) and str(df.iloc[0, 0]).startswith('#'):
        df = df.iloc[1:]
    parts.append(df)
tax = pd.concat(parts, axis=0, ignore_index=True).rename(columns={parts[0].columns[0]: 'Feature ID'})
print(f"合并行数: {len(tax):,}")

fa_ids = set()
with open(FASTA) as fh:
    for line in fh:
        if line.startswith('>'):
            fa_ids.add(line[1:].strip())
tax_ids = set(tax['Feature ID'])
dup = tax['Feature ID'].duplicated().sum()
print(f"fasta ASV {len(fa_ids):,} | taxonomy ASV {len(tax_ids):,} | 重复 ID {dup}")
assert dup == 0, "重复 Feature ID"
assert fa_ids == tax_ids, f"集合不一致 fasta-only={len(fa_ids-tax_ids)} tax-only={len(tax_ids-fa_ids)}"
print("✓ 与 fasta 完全一致")

conf = pd.to_numeric(tax['Confidence'], errors='coerce')
has_g = tax['Taxon'].str.contains(r'g__[A-Za-z0-9]', na=False)
print(f"注释到 genus: {has_g.sum():,} ({has_g.mean()*100:.2f}%); Confidence 中位 {np.nanmedian(conf):.3f}")
tax[['Feature ID', 'Taxon', 'Confidence']].to_csv(OUT, sep='\t', index=False)
print(f"写出: {OUT}  ({len(tax):,} 行)")
