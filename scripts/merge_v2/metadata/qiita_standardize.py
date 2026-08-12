#!/usr/bin/env python3
"""Qiita metadata 标准化 → 公共 schema。源 = ENA 桥 + 873 sample 模板 + 2044 prep 模板。

映射见 .claude/data/metadata_merge_design.md「Qiita」段：
- ENA 桥：run_accession→Run、sample_accession→BioSample、study_accession→Project_ID、instrument_platform→Platform
- sample 模板：lat/lon、host_subject_id→Subject_ID(QIITA:)、sex→Sex、country→Country(剥 GAZ:)、
  host_age(+units→years)→Age、host_body_mass_index→BMI、collection_timestamp→Collection_Date、host_height/weight→Height/Weight
- prep 模板：target_subfragment→Region_16S、target_gene=16S→Sequencing_Type=Amplicon、instrument_model→Instrument
- QIITA_ 前缀（标准化核心列）：sample_type/host_body_site/empo_1/2/3/env_material/scientific_name/host_taxid/geo_loc_name/qiita_study_id
- obs 主键 = study.sample（原生）；run+BioSample 仅 ENA 缴存子集有
"""
import sys, glob, re
from pathlib import Path
import numpy as np, pandas as pd, h5py
ROOT = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
sys.path.insert(0, str(ROOT / "scripts/merge_v2"))
import common_schema as cs

QR = ROOT / "rawdata/Qiita"
TPL = QR / "templates"
COMP = ROOT / "results/feature_table/qiita.gg2.genus.h5ad"
OUT  = ROOT / cs.OUT_DIR / "qiita.std.parquet"
DB   = "Qiita"

_MISSING={"missing","not applicable","not collected","not provided","na","none","nan","unknown","",
          "control sample","not provided ","restricted access"}
def clean(v):
    if v is None or (isinstance(v,float) and pd.isna(v)): return pd.NA
    s=str(v).strip()
    return s if s and s.lower() not in _MISSING else pd.NA
def num_or_na(v):
    s=clean(v)
    if s is pd.NA: return pd.NA
    try: return float(s)
    except ValueError: return pd.NA
def clean_sex(v):
    s=("" if v is None else str(v)).strip().lower()
    if s in ("female","f"): return "Female"
    if s in ("male","m"):   return "Male"
    return pd.NA

UNIT_Y={"years":1,"year":1,"yr":1,"yrs":1,"y":1,"months":1/12,"month":1/12,"mo":1/12,
        "weeks":1/52.1775,"week":1/52.1775,"wk":1/52.1775,"days":1/365.25,"day":1/365.25}
def age_years(val, unit):
    n=num_or_na(val)
    if n is pd.NA: return pd.NA
    u=("" if unit is None else str(unit)).strip().lower()
    return float(n)*UNIT_Y.get(u, 1.0)            # 缺单位默认 years（占比最高）

CFIX={"usa":"United States","united states of america":"United States","us":"United States",
      "united states":"United States","uk":"United Kingdom","united kingdom":"United Kingdom",
      "britain":"United Kingdom","great britain":"United Kingdom"}
def norm_country(v):
    s=clean(v)
    if s is pd.NA: return pd.NA
    s=re.sub(r"^GAZ:","",str(s)).strip()
    return CFIX.get(s.lower(), s)

# ---- 行集 = qiita 组成表 obs（study.sample）----
with h5py.File(COMP,"r") as f:
    g=f["obs"]; k=g.attrs.get("_index","_index"); k=k.decode() if isinstance(k,bytes) else k
    SS=[x.decode() if isinstance(x,bytes) else str(x) for x in g[k][:]]
ss_set=set(SS)
print(f"Qiita 组成表样本 = {len(SS):,}")
out=pd.DataFrame(index=pd.Index(["Qiita:"+s for s in SS]))
base=pd.Series(SS)   # study.sample 顺序

# ---- 1. ENA 桥 → Run/BioSample/Project_ID/Platform ----
br=pd.read_csv(QR/"qiita_ena_bridge.tsv",sep="\t",dtype=str,keep_default_na=False)
br["ss"]=br["sample_alias"].str.split(":",n=1).str[-1]
br=br[br["ss"].isin(ss_set)].drop_duplicates("ss",keep="first").set_index("ss")
out["Run"]=base.map(br["run_accession"]).map(clean).values
out["BioSample"]=base.map(br["sample_accession"]).map(clean).values
out["Project_ID"]=base.map(br["study_accession"]).map(clean).values
out["Platform"]=base.map(br["instrument_platform"]).map(clean).values
print(f"  桥命中 {br.shape[0]:,} / {len(SS):,}")

# ---- 2. sample 模板 → 公共 + QIITA_ ----
SAMPLE_WANT=["latitude","longitude","host_subject_id","sex","country","host_age","host_age_units",
             "host_body_mass_index","collection_timestamp","collection_date","host_height","host_weight",
             "sample_type","host_body_site","empo_1","empo_2","empo_3","env_material","scientific_name",
             "host_taxid","geo_loc_name","qiita_study_id"]
frames=[]
for fp in glob.glob(str(TPL/"*.txt")):
    if "_prep_" in fp: continue
    try:
        t=pd.read_csv(fp,sep="\t",dtype=str,keep_default_na=False,low_memory=False)
    except Exception: continue
    t.columns=[c.strip().lower() for c in t.columns]
    if "sample_name" not in t.columns: continue
    keep=["sample_name"]+[c for c in SAMPLE_WANT if c in t.columns]
    frames.append(t[keep])
st=pd.concat(frames,ignore_index=True).drop_duplicates("sample_name",keep="first").set_index("sample_name")
print(f"  sample 模板合并 {st.shape}  ∩组成表={len(set(st.index)&ss_set):,}")
def scol(c): return base.map(st[c]) if c in st.columns else pd.Series([pd.NA]*len(SS))

out["Latitude"]=scol("latitude").map(num_or_na).values
out["Longitude"]=scol("longitude").map(num_or_na).values
out["Subject_ID"]=[cs.namespace_id(DB,clean(v)) for v in scol("host_subject_id")]
out["Sex"]=scol("sex").map(clean_sex).values
out["Country"]=scol("country").map(norm_country).values
out["Age"]=[age_years(a,u) for a,u in zip(scol("host_age"),scol("host_age_units"))]
out["Age_Group"]=cs.age_group_from_age(out["Age"]).values
out["BMI"]=scol("host_body_mass_index").map(num_or_na).values
ct=scol("collection_timestamp"); cd=scol("collection_date")
out["Collection_Date"]=[ (clean(a) if clean(a) is not pd.NA else clean(b)) for a,b in zip(ct,cd) ]
out["Height"]=scol("host_height").map(num_or_na).values
out["Weight"]=scol("host_weight").map(num_or_na).values

prefix=[]
QMAP={"sample_type":"QIITA_Sample_Type","host_body_site":"QIITA_Host_Body_Site","empo_1":"QIITA_empo_1",
      "empo_2":"QIITA_empo_2","empo_3":"QIITA_empo_3","env_material":"QIITA_Env_Material",
      "scientific_name":"QIITA_Scientific_Name","host_taxid":"QIITA_Host_Taxid",
      "geo_loc_name":"QIITA_Geo_Raw","qiita_study_id":"QIITA_Study_ID"}
for src,dst in QMAP.items():
    out[dst]=scol(src).map(clean).values; prefix.append(dst)

# ---- 3. prep 模板 → Region_16S/Sequencing_Type/Instrument ----
PREP_WANT=["sample_name","target_subfragment","target_gene","instrument_model"]
pframes=[]
for fp in glob.glob(str(TPL/"*_prep_*.txt")):
    try:
        t=pd.read_csv(fp,sep="\t",dtype=str,keep_default_na=False,low_memory=False)
    except Exception: continue
    t.columns=[c.strip().lower() for c in t.columns]
    if "sample_name" not in t.columns: continue
    keep=[c for c in PREP_WANT if c in t.columns]
    if len(keep)>1: pframes.append(t[keep])
pp=pd.concat(pframes,ignore_index=True).drop_duplicates("sample_name",keep="first").set_index("sample_name")
def pcol(c): return base.map(pp[c]) if c in pp.columns else pd.Series([pd.NA]*len(SS))
def norm_region(v):
    s=clean(v)
    if s is pd.NA: return pd.NA
    s=str(s).upper().replace(" REGION","").replace("V3V4","V3-V4").replace("V3-4","V3-V4").replace("V34","V3-V4")
    return s if s.startswith("V") else pd.NA
out["Region_16S"]=pcol("target_subfragment").map(norm_region).values
out["Instrument"]=pcol("instrument_model").map(clean).values
tg=pcol("target_gene").map(lambda v: "Amplicon" if "16s" in str(v).lower() else pd.NA)
out["Sequencing_Type"]=tg.values

out=cs.finalize(out, DB, prefix)
assert out.index.is_unique, "Qiita sample_id 不唯一"
for c in ["Sex","Sequencing_Type"]:
    bad=set(pd.Series(out[c]).dropna().unique())-cs.VOCAB[c]
    assert not bad, f"{c} 非法值 {bad}"
print(f"列数={out.shape[1]}（公共{len(cs.COMMON_COLS)}+QIITA_{len(prefix)}）")
for c in ["Run","BioSample","Project_ID","Platform","Latitude","Subject_ID","Sex","Country","Age","BMI","Collection_Date","Region_16S","Instrument","Sequencing_Type","QIITA_Host_Body_Site","QIITA_empo_3"]:
    print(f"  {c:20s}: 非空 {pd.Series(out[c]).notna().mean()*100:5.1f}%")
OUT.parent.mkdir(parents=True, exist_ok=True)
out.to_parquet(OUT)
print(f"\n✓ 写出 {OUT}  {out.shape}")
