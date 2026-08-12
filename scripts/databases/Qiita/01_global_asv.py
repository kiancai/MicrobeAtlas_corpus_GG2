#!/usr/bin/env python
"""Qiita Layer1 步01(B): 把 685 个 per-study asv.fa 的序列**全局去重** → 一个全局 fasta。
header = md5(seq)(确定性, collapse 时按 md5 查全局 taxonomy, 不需映射文件)。
长度不齐无妨: 按精确序列去重, 90/100/150nt 各算独立 ASV; full 分类器逐序列注释。"""
import glob, os, hashlib

ROOT = "/hpcdisk1/limk_group/caiqy/project/260428_greengene2"
STUDIES = f"{ROOT}/rawdata/Qiita/work/studies"
OUT = f"{ROOT}/rawdata/Qiita/work/asv_global.fa"

seen = {}            # md5 -> seq
n_in = 0
fas = glob.glob(f"{STUDIES}/*/asv.fa")
print(f"per-study asv.fa: {len(fas)}", flush=True)
for k, fa in enumerate(fas):
    with open(fa) as f:
        for line in f:
            if line[0] == ">":
                continue
            seq = line.strip()
            if not seq:
                continue
            n_in += 1
            m = hashlib.md5(seq.encode()).hexdigest()
            if m not in seen:
                seen[m] = seq
    if (k + 1) % 100 == 0:
        print(f"  ...{k+1}/{len(fas)} study, 累计输入 {n_in:,}, 唯一 {len(seen):,}", flush=True)

with open(OUT, "w") as o:
    for m, seq in seen.items():
        o.write(f">{m}\n{seq}\n")

print(f"\n输入序列(含跨study重复): {n_in:,}")
print(f"全局唯一 ASV: {len(seen):,}  (去重率 {100*(1-len(seen)/n_in):.1f}%)")
# 长度分布
from collections import Counter
lc = Counter(len(s) for s in seen.values())
print("唯一 ASV 长度分布:", dict(sorted(lc.items())))
print(f"写出: {OUT}")
