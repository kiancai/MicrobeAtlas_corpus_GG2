#!/usr/bin/env python3
"""MA metadata 标准化 → 公共 schema（全 2.69M）。源 = samples.env.info.tsv（9 列，全样本）。

复用 databases/MicrobeAtlas/06 的派生逻辑（MA metadata 全部来自这 9 列，无外部源），
但对齐到 **全 2.69M**（gg2.full obs，非 QC 子集），并映射到 canonical：
- 公共：Run/BioSample/Project_ID/Latitude/Longitude；Sequencing_Type(AMPLICON→Amplicon/WGS→Shotgun/RNAseq→NA)；
        Sex(female→Female,male→Male)；Smoking(Smoker→Smoker，其余→NA，**绝不填 Non-smoker**)
- MA_ 前缀：Sample_ID/Env_*(8)/IsHuman/SampleSite(_Raw)/Health/AgeGroup/Keywords/Geo_Raw/Institution/Technology_Raw
- 无：Age/BMI/Country/Continent/PMID/Library_Layout/Platform/Instrument/Region_16S/Subject_ID/Collection_Date/Height/Weight → NA
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, h5py
ROOT = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
sys.path.insert(0, str(ROOT / "scripts/merge_v2"))
import common_schema as cs

TSV  = ROOT / "rawdata/MicrobeAtlas/sample_info/samples.env.info.tsv"
COMP = ROOT / "results/feature_table/gg2.full.h5ad"
OUT  = ROOT / cs.OUT_DIR / "ma.std.parquet"
DB   = "MA"

ENV_MAINS = ["animal","soil","aquatic","plant"]
SITE = {"gut","skin","oral","urogenital","lung","gastric","nose","bone"}
SEX  = {"female","male"}
AGE  = {"infant","baby","toddler","adult","elderly"}
HEALTH = {"healthy","disease","infection","inflammatory bowel disease","dermatitis","cystic fibrosis",
          "cholera","malaria","pneumonia","typhoid fever","tuberculosis","patient"}
SEQ_MAP = {"AMPLICON":"Amplicon","WGS":"Shotgun"}        # RNAseq→NA

# ---- 行集 = gg2.full obs（全 2.69M），对齐 tsv ----
with h5py.File(COMP,"r") as f:
    g=f["obs"]; k=g.attrs.get("_index","_index"); k=k.decode() if isinstance(k,bytes) else k
    obs_ids=[x.decode() if isinstance(x,bytes) else str(x) for x in g[k][:]]
print(f"MA 组成表样本 = {len(obs_ids):,}")

df = pd.read_csv(TSV, sep="\t", dtype=str, keep_default_na=False, na_filter=False, low_memory=False)
df.columns = ["MAP_SID","Environments","_col3","Technology","Keywords","_col6","Project","Institution","_col9"]
df = df.set_index("MAP_SID").reindex(obs_ids)
assert df.notna().any(axis=1).all(), "存在 obs 在 tsv 找不到"
print(f"对齐 {len(df):,} 行")

sid = pd.Series(obs_ids)
split = sid.str.split(".", n=1, expand=True)
run_id, biosample = split[0].values, split[1].values

# ---- col2 Environments → Env_* + IsHuman ----
def parse_env(s):
    if not s: return []
    return [(t.split(";",1)+[""])[:2] if ";" in t else (t,"") for t in s.split("|")]
penv = df["Environments"].fillna("").map(parse_env)
env_cols={}
for m in ENV_MAINS:
    flags=penv.map(lambda lst,mm=m: any(t[0]==mm for t in lst))
    subs =penv.map(lambda lst,mm=m: "|".join(t[1] for t in lst if t[0]==mm))
    subs =subs.where(flags, other=pd.NA)
    env_cols[f"MA_Env_{m.capitalize()}"]=flags.values
    env_cols[f"MA_Env_{m.capitalize()}_Sub"]=subs.values
def human_cls(lst):
    asub=[t[1] for t in lst if t[0]=="animal"]
    if "human" not in asub: return pd.NA
    return "HumanMix" if [s for s in asub if s not in ("human","")] else "Human"
ma_human=penv.map(human_cls).values

# ---- col3 → site/sex/age/smoking/health（白名单）----
def pcol3(s): return [t for t in s.split(";") if t] if s else []
p3=df["_col3"].fillna("").map(pcol3)
def pick(toks,allow):
    for t in toks:
        if t in allow: return t
    return pd.NA
site=p3.map(lambda x:pick(x,SITE)).values
sex =p3.map(lambda x:pick(x,SEX)).values
age =p3.map(lambda x:pick(x,AGE)).values
smk =p3.map(lambda x:"Smoker" if any(t=="smoker" for t in x) else pd.NA).values
health=p3.map(lambda x:";".join(t for t in x if t in HEALTH) or pd.NA).values
site_raw=df["_col3"].where(df["_col3"]!="", other=pd.NA).values

# ---- col9 → lat/lon ----
def pgeo(s):
    p=s.split()
    if len(p)!=2: return (np.nan,np.nan)
    try: return (float(p[0]),float(p[1]))
    except ValueError: return (np.nan,np.nan)
geo=df["_col9"].fillna("").map(pgeo)
lat=geo.map(lambda t:t[0]).astype(float).values
lon=geo.map(lambda t:t[1]).astype(float).values

def na(x): return df[x].where(df[x]!="", other=pd.NA).values

# ---- 组装 ----
out=pd.DataFrame(index=pd.Index(["MA:"+s for s in obs_ids]))
out["Run"]=run_id; out["BioSample"]=biosample
out["Project_ID"]=na("Project")
out["Sequencing_Type"]=[SEQ_MAP.get(t, pd.NA) for t in df["Technology"].values]
out["Sex"]=[ {"female":"Female","male":"Male"}.get(s, pd.NA) for s in sex ]
out["Smoking"]=smk
out["Latitude"]=lat; out["Longitude"]=lon
# MA_ 前缀
out["MA_Sample_ID"]=obs_ids
for kk,vv in env_cols.items(): out[kk]=vv
out["MA_IsHuman"]=ma_human
out["MA_SampleSite_Raw"]=site_raw; out["MA_SampleSite"]=site
out["MA_Health"]=health; out["MA_AgeGroup"]=age
out["MA_Keywords"]=na("Keywords"); out["MA_Geo_Raw"]=na("_col9"); out["MA_Institution"]=na("Institution")
out["MA_Technology_Raw"]=na("Technology")
prefix_cols=["MA_Sample_ID"]+list(env_cols)+["MA_IsHuman","MA_SampleSite_Raw","MA_SampleSite",
             "MA_Health","MA_AgeGroup","MA_Keywords","MA_Geo_Raw","MA_Institution","MA_Technology_Raw"]

out=cs.finalize(out, DB, prefix_cols)
assert out.index.is_unique, "MA sample_id 不唯一"
for c in ["Sex","Sequencing_Type","Smoking"]:
    bad=set(pd.Series(out[c]).dropna().unique())-cs.VOCAB[c]
    assert not bad, f"{c} 非法值 {bad}"
print(f"列数={out.shape[1]}（公共{len(cs.COMMON_COLS)}+MA_{len(prefix_cols)}）")
for c in ["BioSample","Sequencing_Type","Sex","Smoking","Latitude","MA_SampleSite","MA_IsHuman"]:
    print(f"  {c:18s}: 非空 {pd.Series(out[c]).notna().mean()*100:5.1f}%")
print("Smoking 值:", dict(pd.Series(out['Smoking']).value_counts(dropna=False)))
OUT.parent.mkdir(parents=True, exist_ok=True)
out.to_parquet(OUT)
print(f"\n✓ 写出 {OUT}  {out.shape}")
