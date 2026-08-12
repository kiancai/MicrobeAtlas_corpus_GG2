#!/usr/bin/env python3
"""RM metadata 标准化 → 公共 schema。源 = metadata_all.standardized.fixed.parquet（已近乎对齐）。

映射见 .claude/data/metadata_merge_design.md「RM」段：
- 公共：Run/Project_ID/BioSample/PMID/Sex/Seq_Type/Library_Layout/Platform/Lat/Lon/
        Age(+start/end)/BMI(+start/end)/Country/Continent/Region_16S；Model→Instrument；
        Age_Group 中央从 Age 派生；Patient_ID→Subject_ID(RM: 命名空间)
- RM_ 前缀：Phenotype/Phenotype_ID/Disease_Stage/Complication/Intervention/Antibiotics_Used/
            Sample_Site/Sample_Type/Case_Or_Control/Is_Healthy/Recent_Antibiotic_Use/Time_Point/Location
- 无：Collection_Date/Height/Weight → NA
"""
import sys
from pathlib import Path
import pandas as pd
ROOT = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
sys.path.insert(0, str(ROOT / "scripts/merge_v2"))
import common_schema as cs

SRC = ROOT / "results/feature_table/metadata_all.standardized.fixed.parquet"
OUT = ROOT / cs.OUT_DIR / "rm.std.parquet"
DB  = "RM"

rm = pd.read_parquet(SRC)
print(f"源 {rm.shape}")
out = pd.DataFrame(index=("RM:" + rm["Run"].astype(str)).values)

# 直传公共列
for c in ["Run","Project_ID","BioSample","PMID","Platform","Region_16S","Sex",
          "Age","Age_start","Age_end","BMI","BMI_start","BMI_end",
          "Country","Continent","Latitude","Longitude","Smoking"]:
    out[c] = rm[c].values
# 映射/派生
out["Instrument"]      = rm["Model"].values
out["Sequencing_Type"] = [cs.norm_seqtype(v) for v in rm["Sequencing_Type"].values]
out["Library_Layout"]  = [cs.norm_layout(v)  for v in rm["Library_Layout"].values]
out["Age_Group"]       = cs.age_group_from_age(rm["Age"].values).values
out["Subject_ID"]      = [cs.namespace_id(DB, v) for v in rm["Patient_ID"].values]

# RM_ 前缀列
RM_MAP = {
    "Phenotype":"RM_Phenotype", "Phenotype_ID":"RM_Phenotype_ID", "Disease_Stage":"RM_Disease_Stage",
    "Complication":"RM_Complication", "Intervention":"RM_Intervention", "Antibiotics_Used":"RM_Antibiotics_Used",
    "Sample_Site":"RM_Sample_Site", "Sample_Type":"RM_Sample_Type", "Case_Or_Control":"RM_Case_Or_Control",
    "Is_Healthy":"RM_Is_Healthy", "Recent_Antibiotic_Use":"RM_Recent_Antibiotic_Use",
    "Time_Point":"RM_Time_Point", "Location":"RM_Location",
}
for src,dst in RM_MAP.items():
    out[dst] = rm[src].values
prefix_cols = list(RM_MAP.values())

out = cs.finalize(out, DB, prefix_cols)

# ---- 校验 ----
assert out.index.is_unique, "RM sample_id 不唯一"
# Age_Group 中央派生 vs RM 自带（应高度一致）
both = pd.DataFrame({"central": out["Age_Group"].astype("object").values,
                     "rm": rm["Age_Group"].astype("object").values})
mask = both["central"].notna() & both["rm"].notna()
agree = (both.loc[mask,"central"] == both.loc[mask,"rm"]).mean()
print(f"Age_Group 中央派生 vs RM 自带 一致率 = {agree:.4f}（n={int(mask.sum()):,}）")
for c in ["Sex","Sequencing_Type","Smoking","Library_Layout"]:
    bad = set(out[c].dropna().unique()) - cs.VOCAB[c]
    assert not bad, f"{c} 出现非法值 {bad}"
print(f"列数 = {out.shape[1]}（公共 {len(cs.COMMON_COLS)} + RM_ {len(prefix_cols)}）")
print("非空率(公共关键):")
for c in ["BioSample","Sex","Age","Country","Region_16S","Subject_ID","Sequencing_Type"]:
    print(f"  {c:16s}: {out[c].notna().mean()*100:5.1f}%")

OUT.parent.mkdir(parents=True, exist_ok=True)
out.to_parquet(OUT)
print(f"\n✓ 写出 {OUT}  {out.shape}")
