#!/usr/bin/env python3
"""MBM metadata 标准化 → 公共 schema。源 = sample_metadata.tsv + projects.csv + tags.tsv.gz(EAV)。

映射见 .claude/data/metadata_merge_design.md「MicroBioMap」段：
- sample_metadata：srr→Run、project→Project_ID、library_strategy=AMPLICON→Sequencing_Type=Amplicon、
  instrument→Instrument、iso→Country+Continent、srs→MBM_SRS、geo_loc_name→MBM_Geo_Raw、region→MBM_Geo_Region
- projects（按 project 广播）：amplicon→Region_16S、sample_type→MBM_Sample_Type、condition→MBM_Condition、
  bead_beating/kit/avg_length/subjects/notes→MBM_*、link→MBM_DOI
- tags（EAV 抽取）：External Id→BioSample(大写)、host_subject_id→Subject_ID(MBM:)、host_sex/sex→Sex、
  host_body_mass_index/body_mass_index→BMI、host_age/age 数值→Age、collection_date→Collection_Date、host_height→Height、host→MBM_Host
- 丢：pubdate/total_bases/library_source
"""
import sys, gzip, csv
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
sys.path.insert(0, str(ROOT / "scripts/merge_v2"))
import common_schema as cs
import iso_country as iso

WORK = ROOT / "rawdata/MicroBioMap/work"
OUT  = ROOT / cs.OUT_DIR / "mbm.std.parquet"
DB   = "MBM"

_MISSING = {"missing","not applicable","not collected","not provided","na","none","nan",""}
def clean(v):
    if v is None or (isinstance(v,float) and pd.isna(v)): return pd.NA
    v = str(v).strip()
    return v if v and v.lower() not in _MISSING else pd.NA
def num_or_na(v):
    try:
        v=str(v).strip()
        return float(v) if v else pd.NA
    except (ValueError,TypeError): return pd.NA
def clean_sex(v):
    s=("" if v is None else str(v)).strip().lower()
    if s in ("female","f"): return "Female"
    if s in ("male","m"):   return "Male"
    return pd.NA

# ---- 1. sample_metadata（主表，每样本）----
sm = pd.read_csv(WORK/"sample_metadata.tsv", sep="\t", dtype=str, keep_default_na=False)
print(f"sample_metadata {sm.shape}")
out = pd.DataFrame(index=("MBM:"+sm["srr"].astype(str)).values)
srr = sm["srr"].astype(str).values
out["Run"]=srr
out["Project_ID"]=sm["project"].map(clean).values
out["Sequencing_Type"]=[ "Amplicon" if str(s).strip().upper()=="AMPLICON" else pd.NA for s in sm["library_strategy"] ]
out["Instrument"]=sm["instrument"].map(lambda v: pd.NA if str(v).strip().lower() in ("unspecified","") else clean(v)).values
out["Country"]=sm["iso"].map(iso.country_of).values
out["Continent"]=sm["iso"].map(iso.continent_of).values
out["MBM_SRS"]=sm["srs"].map(clean).values
out["MBM_Geo_Raw"]=sm["geo_loc_name"].map(clean).values
out["MBM_Geo_Region"]=sm["region"].map(lambda v: pd.NA if str(v).strip().lower()=="unknown" else clean(v)).values
prefix=["MBM_SRS","MBM_Geo_Raw","MBM_Geo_Region"]

# ---- 2. projects（项目级，按 project 广播）----
pj = pd.read_csv(WORK/"projects.csv", dtype=str, keep_default_na=False)
pj.columns=[c.lstrip("﻿") for c in pj.columns]
pj=pj.set_index("project")
proj_of = sm["project"].astype(str).values
def bcast(col):
    m=pj[col].to_dict() if col in pj.columns else {}
    return pd.Series([clean(m.get(p,"")) for p in proj_of]).values
out["Region_16S"]=bcast("amplicon")
PJ_MAP={"sample_type":"MBM_Sample_Type","condition":"MBM_Condition","bead_beating":"MBM_Bead_Beating",
        "kit":"MBM_Kit","avg_length":"MBM_Avg_Length","subjects":"MBM_Subjects","notes":"MBM_Notes","link":"MBM_DOI"}
for src,dst in PJ_MAP.items():
    out[dst]=bcast(src); prefix.append(dst)

# ---- 3. tags（EAV → 按 srr 抽取选定键）----
WANT={"External Id","host_subject_id","host_sex","sex","host_body_mass_index","body_mass_index",
      "host_age","age","collection_date","host_height","host"}
buckets={k:{} for k in WANT}
with gzip.open(WORK/"tags.tsv.gz","rt") as f:
    rd=csv.reader(f,delimiter="\t"); next(rd)
    for row in rd:
        if len(row)<5: continue
        _,run,_,tag,val = row[0],row[1],row[2],row[3],row[4]
        if tag in WANT and run not in buckets[tag]:
            buckets[tag][run]=val
def col(tag): return pd.Series(srr).map(buckets[tag]).values
# BioSample（External Id 小写 samea → 大写）
ext=[clean(v) for v in col("External Id")]
out["BioSample"]=[ (x.upper() if isinstance(x,str) else pd.NA) for x in ext ]
# Subject_ID
out["Subject_ID"]=[ cs.namespace_id(DB,v) for v in col("host_subject_id") ]
# Sex：host_sex 优先，sex 兜底
hs=col("host_sex"); sx=col("sex")
out["Sex"]=[ clean_sex(a) if clean_sex(a) is not pd.NA else clean_sex(b) for a,b in zip(hs,sx) ]
# BMI：host_body_mass_index 优先
hb=col("host_body_mass_index"); bb=col("body_mass_index")
out["BMI"]=[ (num_or_na(a) if num_or_na(a) is not pd.NA else num_or_na(b)) for a,b in zip(hb,bb) ]
# Age：host_age 优先，纯数值
ha=col("host_age"); ag=col("age")
out["Age"]=[ (num_or_na(a) if num_or_na(a) is not pd.NA else num_or_na(b)) for a,b in zip(ha,ag) ]
out["Age_Group"]=cs.age_group_from_age(out["Age"]).values
out["Collection_Date"]=[ clean(v) for v in col("collection_date") ]
out["Height"]=[ num_or_na(v) for v in col("host_height") ]
out["MBM_Host"]=[ clean(v) for v in col("host") ]; prefix.append("MBM_Host")

out=cs.finalize(out, DB, prefix)
assert out.index.is_unique, "MBM sample_id 不唯一"
for c in ["Sex","Sequencing_Type"]:
    bad=set(pd.Series(out[c]).dropna().unique())-cs.VOCAB[c]
    assert not bad, f"{c} 非法值 {bad}"
print(f"列数={out.shape[1]}（公共{len(cs.COMMON_COLS)}+MBM_{len(prefix)}）")
for c in ["Run","Project_ID","BioSample","Country","Continent","Region_16S","Sex","BMI","Age","Subject_ID","Collection_Date","Height","MBM_Sample_Type","MBM_Condition"]:
    print(f"  {c:18s}: 非空 {pd.Series(out[c]).notna().mean()*100:5.1f}%")
OUT.parent.mkdir(parents=True, exist_ok=True)
out.to_parquet(OUT)
print(f"\n✓ 写出 {OUT}  {out.shape}")
