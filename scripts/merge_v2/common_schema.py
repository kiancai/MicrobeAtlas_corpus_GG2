#!/usr/bin/env python3
"""Layer2 metadata canonical schema —— 各库 standardize 脚本共用。

公共列(27)+ 缺失约定 + Age_Group 中央分桶 + 词表。各库脚本：
  产出 per-sample 表，index = 全局 `<DB>:<native>`（对齐 merged.counts.obs_names），
  列 = COMMON_COLS（缺失填 NA）+ 该库 `<DB>_` 前缀列。
落盘 `results/feature_table/metadata/<db>.std.parquet`。
"""
import numpy as np
import pandas as pd

# ---- 公共列（顺序即输出顺序）----
COMMON_COLS = [
    "Database", "Run", "BioSample", "Project_ID", "PMID",
    "Sequencing_Type", "Library_Layout", "Platform", "Instrument", "Region_16S",
    "Sex", "Age", "Age_start", "Age_end", "Age_Group",
    "BMI", "BMI_start", "BMI_end", "Height", "Weight",
    "Smoking", "Country", "Continent", "Latitude", "Longitude", "Collection_Date",
    "Subject_ID",
]

# ---- 受控词表（仅文档/校验用；缺失一律 NA，不编造）----
VOCAB = {
    "Sex": {"Male", "Female"},
    "Sequencing_Type": {"Amplicon", "Shotgun"},
    "Smoking": {"Smoker", "Non-smoker", "Ex-smoker"},          # 缺失=NA（不填 Unknown/Non-smoker）
    "Library_Layout": {"PAIRED", "SINGLE"},
    "Age_Group": {"0-3", "3-18", "18-35", "35-45", "45-60", "60-75", "75+"},
    "Continent": {"Europe", "North America", "Asia", "Africa", "Oceania", "South America", "Antarctica"},
}

# ---- Age_Group 中央分桶（从公共数值 Age(years) 派生，复刻 RM 7 桶）----
_AGE_EDGES  = [0, 3, 18, 35, 45, 60, 75, np.inf]
_AGE_LABELS = ["0-3", "3-18", "18-35", "35-45", "45-60", "60-75", "75+"]

def age_group_from_age(age):
    """数值 Age(years) → 7 桶 Categorical。**右闭 (a,b]**（复刻 RM：边界年龄归低桶，实测一致率 1.0）；NA→NA。"""
    a = pd.to_numeric(pd.Series(age), errors="coerce")
    return pd.cut(a, bins=_AGE_EDGES, labels=_AGE_LABELS, right=True, include_lowest=True)

# ---- Sequencing_Type 归一：amplicon vs shotgun ----
_AMPLICON = {"16s", "full-16s", "amplicon", "its", "18s", "16s rrna", "amplicons"}
_SHOTGUN  = {"metagenomics", "wgts", "wgs", "shotgun", "metatranscriptomics", "virome", "metagenomic"}
def norm_seqtype(v):
    if pd.isna(v): return pd.NA
    s = str(v).strip().lower()
    if s in _AMPLICON: return "Amplicon"
    if s in _SHOTGUN:  return "Shotgun"
    return pd.NA

def norm_layout(v):
    if pd.isna(v): return pd.NA
    s = str(v).strip().upper()
    if s.startswith("PAIR"): return "PAIRED"
    if s.startswith("SING"): return "SINGLE"
    return pd.NA

# ---- 公共数值/地理 sanity（越界→NA；源脏值清洗，2026-06-09）----
SANITY = {   # (lo, hi) 闭区间，外侧或非有限 → NA
    "Age": (0, 122), "Age_start": (0, 122), "Age_end": (0, 122),
    "BMI": (8, 100), "BMI_start": (8, 100), "BMI_end": (8, 100),
    "Height": (30, 250), "Weight": (1, 500),
    "Latitude": (-90, 90), "Longitude": (-180, 180),
}
def apply_common_sanity(df):
    """对公共数值/地理列做 sanity 裁剪 + Region 大写 + (0,0)坐标→NA + Age_Group 重派生。
    源脏值（MA 垃圾坐标 / MBM·Qiita 负年龄 / Weight=999909 哨兵 / inf BMI）越界一律 NA。"""
    df = df.copy()
    for col, (lo, hi) in SANITY.items():
        if col in df.columns:
            v = pd.to_numeric(df[col], errors="coerce")
            df[col] = v.where(np.isfinite(v) & (v >= lo) & (v <= hi))
    if "Latitude" in df.columns and "Longitude" in df.columns:   # 0,0 占位 → NA
        z = (pd.to_numeric(df["Latitude"], errors="coerce") == 0) & \
            (pd.to_numeric(df["Longitude"], errors="coerce") == 0)
        df.loc[z, "Latitude"] = np.nan; df.loc[z, "Longitude"] = np.nan
    if "Region_16S" in df.columns:                               # 大小写统一（MBM 源小写）
        df["Region_16S"] = df["Region_16S"].map(lambda v: str(v).upper() if pd.notna(v) else pd.NA)
    if "Age" in df.columns:                                      # 清洗后重派生 Age_Group
        df["Age_Group"] = age_group_from_age(df["Age"]).astype("object")
    return df

def namespace_id(db, val):
    """id 槽列加 DB 命名空间：`RM:100BAL1`；NA 保持 NA。"""
    if pd.isna(val) or str(val).strip() == "": return pd.NA
    return f"{db}:{str(val).strip()}"

def finalize(df, db, prefix_cols):
    """补全 COMMON_COLS（缺失列→NA）+ 排列：COMMON + 前缀列；index 命名 sample_id。"""
    df = df.copy()
    df["Database"] = db
    for c in COMMON_COLS:
        if c not in df.columns:
            df[c] = pd.NA
    ordered = COMMON_COLS + [c for c in prefix_cols if c in df.columns]
    df = df[ordered]
    df.index.name = "sample_id"
    return df

OUT_DIR = "results/feature_table/metadata"
