#!/usr/bin/env python
"""盘点 Qiita Deblur ASV BIOM: 把每个表的样本对到净新增清单, 统计 artifact/样本/ASV/trim。
只读 HDF5 元数据(/sample/ids, /observation/ids), 不加载矩阵。优先 reference-hit, 否则 all.biom。"""
import os, glob, h5py
from collections import Counter, defaultdict

ROOT = "/hpcdisk1/limk_group/caiqy/project/260428_greengene2"
BIOMDIR = f"{ROOT}/rawdata/Qiita/BIOM"
NETNEW = f"{ROOT}/rawdata/Qiita/qiita_netnew_16s.tsv"
OUT = f"{ROOT}/rawdata/Qiita/qiita_biom_inventory.tsv"

# 载入净新增清单: sample_name -> status (只含 16S 样本)
st = {}
with open(NETNEW) as f:
    h = f.readline().rstrip("\n").split("\t")
    si, ti = h.index("sample_name"), h.index("status")
    for line in f:
        p = line.rstrip("\n").split("\t")
        st[p[si]] = p[ti]
print(f"16S 样本(含状态) {len(st):,}", flush=True)
NETNEW_ST = {"netnew_bridged", "netnew_nobridge"}

arts = sorted(os.listdir(BIOMDIR))
print(f"artifact 目录 {len(arts):,}", flush=True)

fout = open(OUT, "w")
fout.write("artifact\tbiom_type\tstudy\tn_samples\tn_16s\tn_netnew\tn_asv\ttrim_len\n")
agg = Counter()
trim_netnew = Counter()
study_netnew = defaultdict(int)
seen_netnew = set()   # 去重: 同一净新增样本可能在多个 artifact
nfail = 0
done = 0
for a in arts:
    d = f"{BIOMDIR}/{a}"
    bt = None
    for cand in ("reference-hit.biom", "all.biom"):
        if os.path.exists(f"{d}/{cand}"):
            bt = cand; break
    if bt is None:
        continue
    try:
        with h5py.File(f"{d}/{bt}", "r") as hf:
            sids = hf["sample/ids"]
            n_samp = sids.shape[0]
            oids = hf["observation/ids"]
            n_asv = oids.shape[0]
            trim = len(oids[0].decode() if n_asv else "")
            # 判 ASV(序列)还是 OTU: 首 obs 是否纯 ACGTN
            first = oids[0].decode() if n_asv else ""
            is_seq = first[:20].upper().strip("ACGTN") == "" and len(first) >= 50
            samp = [s.decode() for s in sids[:]]
    except Exception:
        nfail += 1; continue
    if not is_seq:
        continue  # 不是 Deblur 序列表
    study = samp[0].split(".")[0] if samp else ""
    n16 = sum(1 for s in samp if s in st)
    nnn = [s for s in samp if st.get(s) in NETNEW_ST]
    fout.write(f"{a}\t{bt}\t{study}\t{n_samp}\t{n16}\t{len(nnn)}\t{n_asv}\t{trim}\n")
    agg["artifacts_deblur"] += 1
    if n16: agg["artifacts_with_16s"] += 1
    if nnn:
        agg["artifacts_with_netnew"] += 1
        trim_netnew[trim] += len(nnn)
        study_netnew[study] += len(nnn)
        seen_netnew.update(nnn)
    done += 1
    if done % 2000 == 0: print(f"  ...{done} 表已读", flush=True)
fout.close()

print("\n" + "=" * 55)
print(f"Deblur 序列表 artifact: {agg['artifacts_deblur']:,}  (含16S {agg['artifacts_with_16s']:,}, 含净新增 {agg['artifacts_with_netnew']:,})")
print(f"读取失败/跳过: {nfail}")
print(f"\n唯一净新增 16S 样本(在 Deblur 表里实际能拿到的): {len(seen_netnew):,}")
print(f"覆盖 study 数: {len([s for s,n in study_netnew.items() if n])}")
print("\n净新增样本 · 按 ASV 长度(trim context)分布:")
for t, n in sorted(trim_netnew.items()):
    print(f"  {t}nt: {n:,} 样本-表次")
print(f"\n含净新增的 study(top 15, 按净新增样本数):")
for s, n in sorted(study_netnew.items(), key=lambda x: -x[1])[:15]:
    print(f"  study {s}: {n:,}")
print(f"\n明细写: {OUT}")
