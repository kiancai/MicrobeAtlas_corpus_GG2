#!/usr/bin/env python3
"""PRIME metadata 标准化 → 公共 schema。源 = rawdata/PRIME/samples.csv（84 列）。

映射见 .claude/data/metadata_merge_design.md「PRIME」段。两陷阱：
  Sequencing_Type=读长配置→PRIME_Read_Config（公共 Sequencing_Type=Amplicon 常量）；
  Host_Age 多分桶→仅纯数值进 Age（分桶/文本→PRIME_Host_Age_Raw）。
其余非公共列一律 PRIME_<原名>；丢 AvgSpotLen/Bases/Release_Date/Create_Date。
"""
import sys, re
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
sys.path.insert(0, str(ROOT / "scripts/merge_v2"))
import common_schema as cs

SRC = ROOT / "rawdata/PRIME/samples.csv"
OUT = ROOT / cs.OUT_DIR / "prime.std.parquet"
DB  = "PRIME"

COUNTRY_FIX = {
    "America (U.S.)": "United States", "Britain (U.K.)": "United Kingdom",
    "Côte d'Ivoire (Ivory Coast)": "Côte d'Ivoire",
}
def fix_country(v):
    v = (v or "").strip()
    if not v: return pd.NA
    return COUNTRY_FIX.get(v, v)

def num_or_na(v):
    try:
        v = str(v).strip()
        if v == "": return pd.NA
        return float(v)
    except (ValueError, TypeError):
        return pd.NA

SMOKE = {  # Smoke_Status → 公共 Smoking（用户拍板：Rarely/Occasionally→Smoker）
    "no": "Non-smoker", "former": "Ex-smoker",
    "yes": "Smoker", "daily": "Smoker", "current": "Smoker",
    "regularly (3-5 times/week)": "Smoker", "occasionally (1-2 times/week)": "Smoker",
    "rarely (a few times/month)": "Smoker",
}
def map_smoke(v):
    s = (v or "").strip().lower()
    return SMOKE.get(s, pd.NA)

def clean(v):  # 空串 → NA
    v = (v or "").strip()
    return v if v else pd.NA

df = pd.read_csv(SRC, dtype=str, keep_default_na=False)
df.columns = [c.lstrip("﻿") for c in df.columns]
print(f"源 {df.shape}")
idx = ("PRIME:" + df["Run"].astype(str)).values
out = pd.DataFrame(index=idx)

# ---- 公共列 ----
out["Run"]            = df["Run"].values
out["Project_ID"]     = df["BioProject"].map(clean).values
out["BioSample"]      = df["BioSample"].map(clean).values
out["Instrument"]     = df["Instrument"].map(clean).values
out["Platform"]       = df["Instrument"].map(lambda v: "ILLUMINA" if "illumina" in (v or "").lower()
                                             else ("LS454" if "454" in (v or "") else pd.NA)).values
out["Region_16S"]     = df["Variable_Region"].map(lambda v: pd.NA if (v or "").strip() in ("","Unknown") else v.strip()).values
out["Sex"]            = df["Host_Sex"].map(lambda v: v.strip() if (v or "").strip() in ("Male","Female") else pd.NA).values
out["BMI"]            = df["Host_BMI"].map(num_or_na).values
out["Country"]        = df["Country"].map(fix_country).values
out["Continent"]      = df["Continent"].map(clean).values
out["Library_Layout"] = df["Library_Layout"].map(cs.norm_layout).values
out["Subject_ID"]     = df["Participant_Id"].map(lambda v: cs.namespace_id(DB, v)).values
out["Age"]            = df["Host_Age"].map(num_or_na).values          # 仅纯数值
out["Age_Group"]      = cs.age_group_from_age(out["Age"]).values
out["Smoking"]        = df["Smoke_Status"].map(map_smoke).values
out["Collection_Date"]= df["Collection_Date"].map(clean).values
out["Height"]         = df["Host_Height"].map(num_or_na).values
out["Weight"]         = df["Host_Weight"].map(num_or_na).values
out["Sequencing_Type"]= "Amplicon"                                    # PRIME 全 16S 扩增子
# PMID 留 NA（DOI→PMID 解析延后）；Doi 存 PRIME_DOI

# ---- PRIME_ 前缀（其余所有列；陷阱列特殊命名）----
consumed = {"Run","BioProject","BioSample","Instrument","Variable_Region","Host_Sex","Host_BMI",
            "Country","Continent","Library_Layout","Participant_Id","Host_Age","Smoke_Status",
            "Collection_Date","Host_Height","Host_Weight"}
drop = {"AvgSpotLen","Bases","Release_Date","Create_Date"}
prefix_cols = []
out["PRIME_Host_Age_Raw"] = df["Host_Age"].map(clean).values; prefix_cols.append("PRIME_Host_Age_Raw")
for col in df.columns:
    if col in consumed or col in drop: continue
    if col == "Sequencing_Type": dst = "PRIME_Read_Config"
    elif col == "Doi":           dst = "PRIME_DOI"
    else:                        dst = f"PRIME_{col}"
    out[dst] = df[col].map(clean).values
    prefix_cols.append(dst)

out = cs.finalize(out, DB, prefix_cols)

# ---- 校验 ----
assert out.index.is_unique, "PRIME sample_id 不唯一"
for c in ["Sex","Sequencing_Type","Smoking","Library_Layout"]:
    bad = set(out[c].dropna().unique()) - cs.VOCAB[c]
    assert not bad, f"{c} 非法值 {bad}"
print(f"列数 = {out.shape[1]}（公共 {len(cs.COMMON_COLS)} + PRIME_ {len(prefix_cols)}）")
print("非空率(公共关键):")
for c in ["BioSample","Sex","Age","BMI","Country","Region_16S","Subject_ID","Smoking","Collection_Date","Height"]:
    print(f"  {c:16s}: {out[c].notna().mean()*100:5.1f}%")
print("Age 数值占比(vs 分桶丢弃):", f"{out['Age'].notna().mean()*100:.1f}%",
      " Age_Raw 非空:", f"{out['PRIME_Host_Age_Raw'].notna().mean()*100:.1f}%")
OUT.parent.mkdir(parents=True, exist_ok=True)
out.to_parquet(OUT)
print(f"\n✓ 写出 {OUT}  {out.shape}")
