#!/usr/bin/env python3
"""MGnify 全库 run × 本项目 run 列表 求重合。
输入：
  /tmp/mgnify_runs.tsv          (accession, secondary_accession, experiment_type, ena_study, is_private)
  /tmp/ma_runs.txt              MA 全集 2,690,735 run
  /tmp/corpus_runs.txt          最终语料 1,793,428 unique run
输出：终端汇总（供写入 RESULTS.md）。
匹配：用 MGnify 的 primary accession，且对 SRS/ERS 样本号补用 secondary_accession 兜底匹配。
"""
import csv, collections

def load_set(path):
    with open(path) as f:
        return set(l.strip() for l in f if l.strip())

ma = load_set("/tmp/ma_runs.txt")
corpus = load_set("/tmp/corpus_runs.txt")
print(f"MA full set      : {len(ma):,}")
print(f"corpus set       : {len(corpus):,}")

# 读 MGnify
rows = []
with open("/tmp/mgnify_runs.tsv") as f:
    r = csv.DictReader(f, delimiter="\t")
    for d in r:
        rows.append((d["accession"], d.get("secondary_accession", ""), d["experiment_type"] or "(empty)"))
print(f"MGnify runs      : {len(rows):,}")

def hit(acc, sec, refset):
    return (acc in refset) or (sec and sec in refset)

# 总体
in_ma   = sum(1 for a, s, e in rows if hit(a, s, ma))
in_corp = sum(1 for a, s, e in rows if hit(a, s, corpus))
N = len(rows)
print("\n=== 总体重合（全部 experiment_type）===")
print(f"MGnify ∩ MA全集   : {in_ma:,}  ({in_ma/N*100:.1f}%)   net-new vs MA   : {N-in_ma:,}")
print(f"MGnify ∩ corpus   : {in_corp:,}  ({in_corp/N*100:.1f}%)   net-new vs corpus: {N-in_corp:,}")

# 按 experiment_type 分解
print("\n=== 按 experiment_type 分解（∩MA / ∩corpus / 该类总数）===")
bytype = collections.defaultdict(lambda: [0, 0, 0])  # type -> [total, in_ma, in_corp]
for a, s, e in rows:
    bytype[e][0] += 1
    if hit(a, s, ma): bytype[e][1] += 1
    if hit(a, s, corpus): bytype[e][2] += 1
print(f"{'type':<20}{'total':>10}{'∩MA':>10}{'∩MA%':>8}{'∩corpus':>10}{'∩corp%':>8}")
for e, (tot, im, ic) in sorted(bytype.items(), key=lambda x: -x[1][0]):
    print(f"{e:<20}{tot:>10,}{im:>10,}{im/tot*100:>7.1f}%{ic:>10,}{ic/tot*100:>7.1f}%")

# amplicon 类（amplicon + metabarcoding，最可比我们的 16S 语料）
amp = [(a, s, e) for a, s, e in rows if e in ("amplicon", "metabarcoding")]
amp_ma = sum(1 for a, s, e in amp if hit(a, s, ma))
amp_cp = sum(1 for a, s, e in amp if hit(a, s, corpus))
print("\n=== amplicon+metabarcoding（最可比 16S 语料）===")
print(f"总数             : {len(amp):,}")
print(f"∩ MA全集         : {amp_ma:,}  ({amp_ma/len(amp)*100:.1f}%)   net-new vs MA   : {len(amp)-amp_ma:,}")
print(f"∩ corpus         : {amp_cp:,}  ({amp_cp/len(amp)*100:.1f}%)   net-new vs corpus: {len(amp)-amp_cp:,}")

# net-new 的前缀分布（不在 corpus 的那些）
print("\n=== net-new vs corpus 的 accession 前缀分布（全部类型）===")
pre = collections.Counter(a[:3] for a, s, e in rows if not hit(a, s, corpus))
for k, v in pre.most_common():
    print(f"  {k}: {v:,}")
