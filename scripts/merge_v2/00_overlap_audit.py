#!/usr/bin/env python3
"""五库样本量 + 两两重叠 + 去重净规模 + ∩MGnify + GMrepo/mBodyMap 覆盖。按 Run 口径。"""
import h5py, csv, gzip, sys
from itertools import combinations
csv.field_size_limit(1<<24)
FT='results/feature_table/'

def h5_index(p, key):
    with h5py.File(p,'r') as f:
        d=f['obs'][key]
        return [x.decode() if isinstance(x,bytes) else str(x) for x in d[:]]

def up(s): return s.strip().upper()

print('== 读取各库 Run 集 ==', flush=True)
# MA: Run.Sample → 取点前缀
ma = set(up(x.split('.',1)[0]) for x in h5_index(FT+'gg2.full.h5ad','sample_id'))
print(f'MA   读完 unique_run={len(ma)}', flush=True)
rm = set(up(x) for x in h5_index(FT+'resmicrodb.gg2.genus.h5ad','sample_id'))
print(f'RM   读完 unique_run={len(rm)}', flush=True)
mbm = set(up(x) for x in h5_index(FT+'microbiomap.gg2.genus.h5ad','sample_id'))
print(f'MBM  读完 unique_run={len(mbm)}', flush=True)

# PRIME: samples.csv Run
prime=set()
with open('rawdata/PRIME/samples.csv',encoding='utf-8-sig') as f:
    r=csv.DictReader(f)
    for row in r:
        v=row.get('Run','').strip()
        if v: prime.add(up(v))
print(f'PRIME 读完 unique_run={len(prime)}', flush=True)

# Qiita: obs=study.sample(321,502)；经 ENA 桥映射 run
qiita_samples=set(h5_index(FT+'qiita.gg2.genus.h5ad','sample_name'))  # study.sample
print(f'Qiita obs 样本数={len(qiita_samples)}', flush=True)
qiita_runs=set(); qiita_mapped_samples=set()
with open('rawdata/Qiita/qiita_ena_bridge.tsv') as f:
    r=csv.DictReader(f,delimiter='\t')
    for row in r:
        alias=row['sample_alias']
        ss=alias.split(':',1)[1] if ':' in alias else alias  # study.sample
        if ss in qiita_samples:
            run=row['run_accession'].strip()
            if run:
                qiita_runs.add(up(run)); qiita_mapped_samples.add(ss)
qiita_norun = len(qiita_samples)-len(qiita_mapped_samples)
print(f'Qiita 有桥(有run)样本={len(qiita_mapped_samples)}  无run样本={qiita_norun}  → run集={len(qiita_runs)}', flush=True)

DBS={'MA':ma,'RM':rm,'MBM':mbm,'PRIME':prime,'Qiita':qiita_runs}

print('\n== 1. 各库 Run 集规模 ==')
for k,v in DBS.items(): print(f'  {k:6s}: {len(v):>9,}')
print(f'  (Qiita 另有 {qiita_norun:,} 样本无 ENA 桥/无 run，必为净新增)')

print('\n== 2. 两两重叠 (∩, run 口径) ==')
for a,b in combinations(DBS,2):
    inter=len(DBS[a]&DBS[b])
    print(f'  {a:6s} ∩ {b:6s} = {inter:>8,}')

print('\n== 3. 并集 + 按纳入顺序的净新增 (MA→RM→MBM→PRIME→Qiita) ==')
order=['MA','RM','MBM','PRIME','Qiita']
seen=set();
for k in order:
    new=len(DBS[k]-seen)
    seen|=DBS[k]
    print(f'  +{k:6s}: 净新增 run={new:>9,}  累计并集={len(seen):>9,}')
union=seen
total_with_qiita_norun = len(union)+qiita_norun
print(f'  五库 Run 并集(去重) = {len(union):,}')
print(f'  + Qiita 无桥样本 {qiita_norun:,} = 五库唯一样本估计 {total_with_qiita_norun:,}')

print('\n== 4. ∩ MGnify (536,566 run；我们已覆盖多少) ==')
mg=set()
with open('rawdata/MGnify/mgnify_run_accessions.sorted.txt') as f:
    for line in f:
        s=line.strip()
        if s: mg.add(up(s))
print(f'  MGnify run 总数 = {len(mg):,}')
inter_mg=union&mg
print(f'  我们五库并集 ∩ MGnify = {len(inter_mg):,}  (占 MGnify {len(inter_mg)/len(mg)*100:.1f}%)')
for k in order:
    print(f'    其中 {k:6s} ∩ MGnify = {len(DBS[k]&mg):>8,}')
print(f'  MGnify 中我们还没有的 = {len(mg-union):,}')

print('\n== 5. GMrepo / mBodyMap 对五库并集的覆盖 (纯 metadata 富集) ==')
def load_runs(path, run_col, sep, gz=False):
    op = gzip.open(path,'rt') if gz else open(path,encoding='utf-8-sig')
    runs=set()
    with op as f:
        r=csv.DictReader(f,delimiter=sep)
        for row in r:
            v=(row.get(run_col) or '').strip().strip('"')
            if v: runs.add(up(v))
    return runs
gmrepo=load_runs('rawdata/GMrepo/sample_to_run_info.txt.gz','run_id','\t',gz=True)
mbody =load_runs('rawdata/mBodyMap/mbodymap_samples.csv','run_id',',')
print(f'  GMrepo run={len(gmrepo):,}   mBodyMap run={len(mbody):,}')
print(f'  GMrepo ∩ mBodyMap = {len(gmrepo&mbody):,} (doc 称 0，验证)')
for nm,st in [('GMrepo',gmrepo),('mBodyMap',mbody)]:
    cov=st&union
    print(f'\n  --- {nm} ---  总 {len(st):,}')
    print(f'    ∩ 五库并集 (可贴 metadata) = {len(cov):,}  ({len(cov)/len(st)*100:.1f}%)')
    print(f'    不在五库 (贴不上,无序列无法加行) = {len(st-union):,}')
    print(f'    其落在各库:')
    for k in order:
        print(f'      ∩ {k:6s} = {len(st&DBS[k]):>7,}')

print('\n== 6. metadata 重复估计 ==')
gm_cov=gmrepo&union; mb_cov=mbody&union
both_meta = (gm_cov|mb_cov)
print(f'  可被 GMrepo 或 mBodyMap 贴到的语料样本(去重) = {len(both_meta):,}')
print(f'  GMrepo+mBodyMap 同时能贴的同一样本 = {len(gm_cov&mb_cov):,} (两库 run 不重叠→预期 0)')
print('  注：富集"重复"主要指这些样本在母库(MA/RM/MBM/PRIME/Qiita)已自带 metadata，')
print('      GMrepo/mBodyMap 的部位/疾病/性别/地理是补充或冗余，按列级合并原则处理。')
