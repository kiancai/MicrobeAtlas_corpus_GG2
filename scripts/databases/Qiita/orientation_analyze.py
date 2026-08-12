#!/usr/bin/env python
"""auto vs both 方向对照分析：
- 全局 60k 样本：ASV 级 Unassigned / genus 解析率 / 平均 Confidence
- 每 study：ASV 级 + reads 加权（用 asv.h5ad counts）reads→genus 对比
判读：both 大幅压低 Unassigned → 混合方向，both 该用；
      both 与 auto 持平但 genus 解析率反降 → MA 式退化，应改 auto。
"""
import os, glob
import numpy as np
import pandas as pd

ROOT = "/hpcdisk1/limk_group/caiqy/project/260428_greengene2"
TAX = f"{ROOT}/results/orientation_check/tax"
STUDY_BASE = f"{ROOT}/rawdata/Qiita/work/studies"


def parse_tax(path):
    """return DataFrame indexed by Feature ID with cols: assigned(domain), genus_ok, conf"""
    df = pd.read_csv(path, sep="\t")
    df.columns = [c.strip() for c in df.columns]
    fid = df["Feature ID"].astype(str)
    tax = df["Taxon"].astype(str)
    conf = pd.to_numeric(df.get("Confidence", np.nan), errors="coerce")

    def genus_ok(t):
        if t.strip() == "Unassigned":
            return False, False
        domain = False
        g = False
        for tok in t.split(";"):
            tok = tok.strip()
            if tok.startswith("d__") and len(tok) > 3:
                domain = True
            if tok.startswith("g__") and len(tok) > 3:
                g = True
        return domain, g

    dg = tax.map(genus_ok)
    out = pd.DataFrame({
        "domain": [x[0] for x in dg],
        "genus": [x[1] for x in dg],
        "conf": conf.values,
    }, index=fid.values)
    return out


def asv_reads(study):
    """per-ASV total reads from asv.h5ad, index=Feature ID(md5)"""
    import anndata as ad
    a = ad.read_h5ad(f"{STUDY_BASE}/{study}/asv.h5ad")
    X = a.X
    tot = np.asarray(X.sum(axis=0)).ravel()  # sum over samples -> per ASV
    return pd.Series(tot, index=a.var_names.astype(str))


def summ(df, w=None):
    n = len(df)
    if w is None:
        return dict(n=n,
                    unassigned=100*(~df["domain"]).mean(),
                    genus=100*df["genus"].mean(),
                    conf=df["conf"].mean())
    w = w.reindex(df.index).fillna(0).values
    tw = w.sum()
    return dict(n=n, reads=int(tw),
                unassigned=100*(w[(~df["domain"]).values].sum())/tw,
                genus=100*(w[df["genus"].values].sum())/tw,
                conf=df["conf"].mean())


# ---------- 全局 60k 样本（ASV 级）----------
print("="*92)
print("【全局随机 60k ASV 样本 — ASV 级】")
print("="*92)
a = parse_tax(f"{TAX}/global60k.auto.tsv")
b = parse_tax(f"{TAX}/global60k.both.tsv")
sa, sb = summ(a), summ(b)
print(f"{'指标':<22}{'auto':>12}{'both':>12}{'Δ(both-auto)':>16}")
for k, lab in [("unassigned","Unassigned %"),("genus","解析到 genus %"),("conf","平均 Confidence")]:
    print(f"{lab:<22}{sa[k]:>12.3f}{sb[k]:>12.3f}{sb[k]-sa[k]:>16.3f}")
# 逐 ASV 一致性：both 相对 auto 救回多少 / 覆盖掉多少
common = a.index.intersection(b.index)
a2, b2 = a.loc[common], b.loc[common]
rescued = ((~a2["genus"]) & (b2["genus"])).sum()   # auto 没到 genus，both 到了
lost = ((a2["genus"]) & (~b2["genus"])).sum()       # auto 到 genus，both 反而没到
print(f"\n逐 ASV: both 救回(auto无→both有 genus) = {rescued:,}   both 丢失(auto有→both无) = {lost:,}   净 = {rescued-lost:+,}")

# ---------- 每 study（reads 加权）----------
print("\n" + "="*92)
print("【15 个 study — reads 加权 reads→genus（核心口径，与 MA/RM 决策一致）】")
print("="*92)
studies = []
for p in sorted(glob.glob(f"{TAX}/study_*.auto.tsv")):
    s = os.path.basename(p).replace("study_","").replace(".auto.tsv","")
    studies.append(s)

rows = []
print(f"{'study':>8}{'n_asv':>8}{'reads':>12}{'una_auto':>10}{'una_both':>10}"
      f"{'gen_auto':>10}{'gen_both':>10}{'Δgenus':>9}")
for s in studies:
    try:
        a = parse_tax(f"{TAX}/study_{s}.auto.tsv")
        b = parse_tax(f"{TAX}/study_{s}.both.tsv")
        w = asv_reads(s)
    except Exception as e:
        print(f"{s:>8}  ERROR {e}")
        continue
    sa, sb = summ(a, w), summ(b, w)
    rows.append((s, sa, sb))
    print(f"{s:>8}{sa['n']:>8}{sa['reads']:>12,}{sa['unassigned']:>10.2f}{sb['unassigned']:>10.2f}"
          f"{sa['genus']:>10.2f}{sb['genus']:>10.2f}{sb['genus']-sa['genus']:>+9.2f}")

# 汇总（reads 加权跨 study 合并）
if rows:
    tot_reads = sum(r[1]['reads'] for r in rows)
    una_a = sum(r[1]['unassigned']*r[1]['reads'] for r in rows)/tot_reads
    una_b = sum(r[2]['unassigned']*r[2]['reads'] for r in rows)/tot_reads
    gen_a = sum(r[1]['genus']*r[1]['reads'] for r in rows)/tot_reads
    gen_b = sum(r[2]['genus']*r[2]['reads'] for r in rows)/tot_reads
    print("-"*92)
    print(f"{'合并(reads加权)':>8}{'':>8}{tot_reads:>12,}{una_a:>10.2f}{una_b:>10.2f}{gen_a:>10.2f}{gen_b:>10.2f}{gen_b-gen_a:>+9.2f}")
    print("\n判读：")
    print(f"  reads→genus  auto={gen_a:.2f}%  both={gen_b:.2f}%  Δ={gen_b-gen_a:+.2f}pp")
    if gen_b - gen_a > 1:
        print("  → both 显著占优：Qiita 确实混合方向，现有 both 决策正确。")
    elif gen_b - gen_a < -1:
        print("  → both 反而掉点（MA 式退化）：Qiita 方向已统一，应改 auto。")
    else:
        print("  → 两者基本持平：both 无明显收益，考虑改 auto（省一半算力，且避免 RC 覆盖风险）。")
