#!/usr/bin/env python
"""MBM auto vs both 方向对照分析（counts 预存在 samples/<label>.counts.tsv）。"""
import os, glob
import numpy as np
import pandas as pd

ROOT = "/hpcdisk1/limk_group/caiqy/project/260428_greengene2"
BASE = f"{ROOT}/results/orientation_check_mbm"
TAX = f"{BASE}/tax"
SAMP = f"{BASE}/samples"


def parse_tax(path):
    df = pd.read_csv(path, sep="\t")
    df.columns = [c.strip() for c in df.columns]
    fid = df["Feature ID"].astype(str)
    tax = df["Taxon"].astype(str)
    conf = pd.to_numeric(df.get("Confidence", np.nan), errors="coerce")

    def dg(t):
        if t.strip() == "Unassigned":
            return False, False
        d = g = False
        for tok in t.split(";"):
            tok = tok.strip()
            if tok.startswith("d__") and len(tok) > 3: d = True
            if tok.startswith("g__") and len(tok) > 3: g = True
        return d, g
    z = tax.map(dg)
    return pd.DataFrame({"domain": [x[0] for x in z], "genus": [x[1] for x in z],
                         "conf": conf.values}, index=fid.values)


def reads(label):
    p = f"{SAMP}/{label}.counts.tsv"
    s = pd.read_csv(p, sep="\t", index_col="asv_id")["reads"]
    s.index = s.index.astype(str)
    return s


def summ(df, w):
    w = w.reindex(df.index).fillna(0).values
    tw = w.sum()
    return dict(n=len(df), reads=int(tw),
                unassigned=100*w[(~df["domain"]).values].sum()/tw,
                genus=100*w[df["genus"].values].sum()/tw,
                conf=df["conf"].mean())


print("="*94)
print("【全局随机 60k ASV — reads 加权 + ASV 级】")
print("="*94)
a = parse_tax(f"{TAX}/global60k.auto.tsv")
b = parse_tax(f"{TAX}/global60k.both.tsv")
w = reads("global60k")
sa, sb = summ(a, w), summ(b, w)
print(f"{'口径':<20}{'auto':>12}{'both':>12}{'Δ':>12}")
for k, lab in [("unassigned","Unassigned %"),("genus","reads→genus %"),("conf","平均 Conf")]:
    print(f"{lab:<20}{sa[k]:>12.3f}{sb[k]:>12.3f}{sb[k]-sa[k]:>12.3f}")
# ASV 级
ag = 100*a["genus"].mean(); bg = 100*b["genus"].mean()
common = a.index.intersection(b.index)
resc = ((~a.loc[common,"genus"]) & (b.loc[common,"genus"])).sum()
lost = ((a.loc[common,"genus"]) & (~b.loc[common,"genus"])).sum()
print(f"ASV 级 genus%: auto={ag:.2f} both={bg:.2f} Δ={bg-ag:+.2f} | both 救回={resc:,} 丢失={lost:,} 净={resc-lost:+,}")

print("\n" + "="*94)
print("【15 个 project — reads 加权 reads→genus】")
print("="*94)
print(f"{'project':>16}{'region':>9}{'reads':>13}{'una_a':>8}{'una_b':>8}{'gen_a':>8}{'gen_b':>8}{'Δgen':>8}")
rows = []
# region 来自 prep 选择（重新读 manifest label 即 proj_<P>，region 需查 projects.csv）
reg = {}
with open(f"{ROOT}/rawdata/MicroBioMap/work/projects.csv", encoding="utf-8-sig") as f:
    h = f.readline().rstrip("\n").split(","); ip, ia = h.index("project"), h.index("amplicon")
    for line in f:
        p = line.rstrip("\n").split(",")
        if len(p) > max(ip, ia): reg[p[ip]] = p[ia] or "?"
for p in sorted(glob.glob(f"{TAX}/proj_*.auto.tsv")):
    lab = os.path.basename(p).replace(".auto.tsv", "")
    proj = lab.replace("proj_", "")
    try:
        a = parse_tax(f"{TAX}/{lab}.auto.tsv"); b = parse_tax(f"{TAX}/{lab}.both.tsv"); w = reads(lab)
    except Exception as e:
        print(f"{proj:>16}  ERROR {e}"); continue
    sa, sb = summ(a, w), summ(b, w)
    rows.append((proj, sa, sb))
    print(f"{proj:>16}{reg.get(proj,'?'):>9}{sa['reads']:>13,}{sa['unassigned']:>8.2f}{sb['unassigned']:>8.2f}"
          f"{sa['genus']:>8.2f}{sb['genus']:>8.2f}{sb['genus']-sa['genus']:>+8.2f}")

if rows:
    tr = sum(r[1]['reads'] for r in rows)
    una_a = sum(r[1]['unassigned']*r[1]['reads'] for r in rows)/tr
    una_b = sum(r[2]['unassigned']*r[2]['reads'] for r in rows)/tr
    gen_a = sum(r[1]['genus']*r[1]['reads'] for r in rows)/tr
    gen_b = sum(r[2]['genus']*r[2]['reads'] for r in rows)/tr
    print("-"*94)
    print(f"{'合并(reads加权)':>16}{'':>9}{tr:>13,}{una_a:>8.2f}{una_b:>8.2f}{gen_a:>8.2f}{gen_b:>8.2f}{gen_b-gen_a:>+8.2f}")
    print(f"\n判读: reads→genus auto={gen_a:.2f}% both={gen_b:.2f}% Δ={gen_b-gen_a:+.2f}pp")
    worse = sum(1 for r in rows if r[2]['genus'] < r[1]['genus'] - 0.5)
    better = sum(1 for r in rows if r[2]['genus'] > r[1]['genus'] + 0.5)
    print(f"  both 更差的 project: {worse}/{len(rows)}   both 更好的: {better}/{len(rows)}")
    if gen_b - gen_a > 1:
        print("  → both 显著占优：MBM 含混合方向，both 决策正确。")
    elif gen_b - gen_a < -1:
        print("  → both 反而掉点：应改 auto（或 per-project auto）。")
    else:
        print("  → 基本持平：both 无明显收益，可改 auto（省算力 + 规避 RC 覆盖）。")
