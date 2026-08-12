#!/usr/bin/env python3
"""Layer2 步2：5 库 std metadata + GMrepo/mBodyMap 富集 → merged.counts 的 obs。

- concat 5 张 *.std.parquet（公共 27 列统一 + 各库前缀列）→ reindex 到 merged.counts.obs_names
- 合并 01 产的旗标列（unresolved_reads/qc_pass/Run_paired_id…），Database/Run 以 merged.counts 为准
- GMrepo/mBodyMap **纯前缀**富集（按 Run join；run-disjoint，不碰公共列）
产物：results/feature_table/metadata/merged.obs.parquet（3,334,541 × 全列）
"""
import sys, gzip
from pathlib import Path
import numpy as np, pandas as pd, h5py
ROOT = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
sys.path.insert(0, str(ROOT / "scripts/merge_v2"))
import common_schema as cs
MD = ROOT / "results/feature_table/metadata"
MERGED = ROOT / "results/feature_table/merged.counts.h5ad"
OUT = MD / "merged.obs.parquet"

# ---- 1. merged.counts 的旗标 obs（顺序基准）----
print("读 merged.counts obs ...", flush=True)
with h5py.File(MERGED,"r") as f:
    g=f["obs"]; k=g.attrs.get("_index","_index"); k=k.decode() if isinstance(k,bytes) else k
    obs_index=[x.decode() if isinstance(x,bytes) else str(x) for x in g[k][:]]
    flags={}
    for c in ["Database","Run","unresolved_reads","total_reads","n_deep_reads","n_deep_features",
              "qc_pass","Run_paired_id","n_run_copies"]:
        d=g[c]
        if isinstance(d,h5py.Group):  # categorical
            cats=[x.decode() if isinstance(x,bytes) else x for x in d["categories"][:]]
            codes=d["codes"][:]; arr=np.array([cats[i] if i>=0 else None for i in codes],dtype=object)
        else:
            arr=d[:]
            if arr.dtype.kind=="S": arr=np.array([x.decode() for x in arr],dtype=object)
        flags[c]=arr
obs=pd.DataFrame(flags, index=pd.Index(obs_index,name="sample_id"))
print(f"  merged obs {obs.shape}")

# ---- 2. concat 5 std → reindex ----
print("concat 5 std metadata ...", flush=True)
parts=[pd.read_parquet(MD/f"{db}.std.parquet") for db in ["ma","rm","mbm","prime","qiita"]]
meta=pd.concat(parts, axis=0)
print(f"  concat {meta.shape}  唯一 index={meta.index.is_unique}")
meta=meta.drop(columns=["Database","Run"], errors="ignore")   # 以 merged.counts 为准
meta=meta.reindex(obs.index)
obs=pd.concat([obs, meta], axis=1)
print(f"  合并旗标+metadata {obs.shape}")

# ---- 公共列 sanity（越界源脏值→NA + Region 大写 + 0,0→NA + Age_Group 重派生）----
obs=cs.apply_common_sanity(obs)
print("  已 sanity：公共数值/地理裁剪 + Region 大写 + 0,0 坐标→NA")

# ---- 3. GMrepo / mBodyMap 富集（纯前缀，按 Run join）----
def up(s): return str(s).strip().upper()
run_key=obs["Run"].map(lambda v: up(v) if pd.notna(v) else None)

def load_enrich(path, sep, gz, pref, cmap):
    op=gzip.open(path,"rt") if gz else open(path,encoding="utf-8-sig")
    df=pd.read_csv(op, sep=sep, dtype=str, keep_default_na=False)
    df.columns=[c.strip().strip('"') for c in df.columns]
    df["__run"]=df["run_id"].map(up)
    df=df.drop_duplicates("__run",keep="first").set_index("__run")
    out={}
    for src,dst in cmap.items():
        if src in df.columns:
            s=df[src].map(lambda v: (v.strip() if isinstance(v,str) and v.strip() and v.strip().lower() not in ("","na","unknown","not provided") else pd.NA))
            out[pref+dst]=run_key.map(s)
    return pd.DataFrame(out, index=obs.index)

GM_MAP={"second_sample_id":"BioSample","experiment_type":"Experiment_Type","instrument_model":"Instrument",
        "disease":"Disease","phenotype":"Phenotype","disease_stage":"Disease_Stage","country":"Country",
        "collection_date":"Collection_Date","sex":"Sex","host_age":"Age","diet":"Diet",
        "longitude":"Longitude","latitude":"Latitude","BMI":"BMI","antibiotics_used":"Antibiotics_Used"}
MB_MAP=dict(GM_MAP); MB_MAP.update({"BodySite":"BodySite","BodySite_info":"BodySite_Info"})

print("富集 GMrepo ...", flush=True)
gm=load_enrich(ROOT/"rawdata/GMrepo/sample_to_run_info.txt.gz","\t",True,"GMREPO_",GM_MAP)
print("富集 mBodyMap ...", flush=True)
mb=load_enrich(ROOT/"rawdata/mBodyMap/mbodymap_samples.csv",",",False,"MBODY_",MB_MAP)
obs=pd.concat([obs, gm, mb], axis=1)
n_gm=gm.notna().any(axis=1).sum(); n_mb=mb.notna().any(axis=1).sum()
print(f"  GMrepo 贴到 {n_gm:,} 行 · mBodyMap 贴到 {n_mb:,} 行")

# ---- 4. 写出 ----
print(f"\n最终 obs {obs.shape}")
print("公共列非空率:")
for c in ["BioSample","Sex","Age","Country","Region_16S","Subject_ID","Collection_Date","BMI","Latitude"]:
    print(f"  {c:16s}: {obs[c].notna().mean()*100:5.1f}%")
obs.to_parquet(OUT)
print(f"✓ 写出 {OUT}  {obs.shape}")
