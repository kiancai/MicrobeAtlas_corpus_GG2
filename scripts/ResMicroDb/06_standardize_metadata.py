# ---
# jupyter:
#   jupytext:
#     formats: py:percent,ipynb
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: baseBio
#     language: python
#     name: python3
# ---

# %% [markdown]
# # ResMicroDb 06: metadata_all 标准化
#
# 把 `rawdata/ResMicroDb/metadata_all.txt` (UTF-16LE, 135,746 × 34) 整理成
# `results/feature_table/metadata_all.standardized.{tsv,parquet}` (UTF-8, 135,746 × 36)。
#
# 与 anndata 的对接（按 obs_names 左 join）放在 07 步，本步只产 sample 级 metadata 表。
#
# 设计文档（含三/四个输入文件关系、列对照、清洗规则、与 jxt 的差异）：
# `.claude/data/metadata_resmicrodb_standardize.md`
#
# **输入**：
# - `rawdata/ResMicroDb/metadata_all.txt` (UTF-16LE + BOM)
# - `rawdata/ResMicroDb/projectTable_changed_250924.tsv` (study 级 16S_Region 来源)
# - `rawdata/ResMicroDb/sampleTable_changed_250924.tsv` (仅做派生列对照校验)
#
# **输出**：
# - `metadata_all.standardized.tsv` —— UTF-8 人类可读；NA 写成空字符串
# - `metadata_all.standardized.parquet` —— 保留 dtype（含 Categorical 顺序 / nullable string /
#   nullable bool），下游 07 步直接 `pd.read_parquet` 即可，无需重转
#
# **36 列**：见 §schema 表 (md doc §2)。改名 6 处、派生 4 列、丢弃冗余 4 列。

# %%
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_DIR = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
META_ALL    = PROJECT_DIR / "rawdata/ResMicroDb/metadata_all.txt"
PROJ_TABLE  = PROJECT_DIR / "rawdata/ResMicroDb/projectTable_changed_250924.tsv"
SAMPLE_TBL  = PROJECT_DIR / "rawdata/ResMicroDb/sampleTable_changed_250924.tsv"  # 仅对照
OUT_TSV     = PROJECT_DIR / "results/feature_table/metadata_all.standardized.tsv"
OUT_PARQUET = PROJECT_DIR / "results/feature_table/metadata_all.standardized.parquet"

EXPECTED_ROWS = 135746
EXPECTED_COLS = 36

# %% [markdown]
# ## 1. Categorical 顺序常量
#
# 仅对"语义上有顺序 / 取值少且固定"的列手工指定；其它 category 列让 pandas 按出现序自动推。

# %%
CAT_ORDER = {
    "Sequencing_Type":      ["16S", "Metagenomics", "Metatranscriptomics",
                              "ITS", "Virome", "Full-16S", "WGTS", "RNA-Seq",
                              "AMPLICON", "18S", "Nanopore", "MeDIP-Seq",
                              "Tn-Seq", "miWTS"],   # 14 unique（不含空）
    "Library_Layout":       ["PAIRED", "SINGLE"],
    "Platform":             ["ILLUMINA", "LS454", "BGISEQ", "ION_TORRENT",
                              "OXFORD_NANOPORE", "PACBIO_SMRT", "DNBSEQ"],
    "Sex":                  ["Female", "Male"],
    "Smoking":              ["Smoker", "Non-smoker", "Ex-smoker"],
    "Recent_Antibiotic_Use":["Yes", "No"],
    "Sample_Site":          ["Nasopharynx", "Nasal", "Sputum", "Oropharynx",
                              "BALF", "Trachea", "Pharynx", "Throat", "Bronchus",
                              "Lung Tissue", "Negative Control", "Positive Control",
                              "Cough swab", "Oral"],
    "Continent":            ["Europe", "North America", "Asia", "Africa",
                              "Oceania", "South America"],
    "Age_Group":            ["0-3", "3-18", "18-35", "35-45", "45-60", "60-75", "75+"],
    "Case_Or_Control":      ["case", "control"],
    "Region_16S":           ["V4", "V3-V4", "V3", "V1-V2", "V1-V3", "V3-V5",
                              "V4-V5", "V5-V7", "V4-V6", "V5-V6",
                              "V1-V3/V3-V5", "V1-V9", "V1-V2/V3-V4",
                              "V6-V8", "V6-V9", "V6"],
}

# 取值多但仍按 category 编码（节省内存；无强制顺序，按 pandas 默认即可）
CAT_NO_ORDER = ["Project_ID", "Country", "Phenotype", "Phenotype_ID",
                "Disease_Stage", "Complication", "Intervention",
                "Antibiotics_Used", "Sample_Type", "Location", "Model"]

# 自由文本（取值过多 / 无意义穷举）→ pandas nullable string
STRING_COLS = ["Run", "BioSample", "PMID", "Patient_ID", "Time_Point"]

FLOAT_COLS = ["Age", "Age_start", "Age_end", "BMI", "BMI_start", "BMI_end",
              "Latitude", "Longitude"]

# 列序（最终输出 36 列顺序）
FINAL_COLUMNS = [
    "Run", "Project_ID", "BioSample", "PMID",
    "Sequencing_Type", "Library_Layout", "Platform", "Model",
    "Phenotype", "Phenotype_ID", "Disease_Stage", "Complication", "Intervention",
    "Smoking", "Recent_Antibiotic_Use", "Antibiotics_Used",
    "Sample_Site", "Sample_Type", "Sex",
    "Age", "Age_start", "Age_end", "Age_Group",
    "BMI", "BMI_start", "BMI_end",
    "Country", "Continent", "Location", "Latitude", "Longitude",
    "Region_16S", "Patient_ID", "Time_Point",
    "Case_Or_Control", "Is_Healthy",
]
assert len(FINAL_COLUMNS) == EXPECTED_COLS

# %% [markdown]
# ## 2. 读 metadata_all.txt（UTF-16LE）

# %%
df = pd.read_csv(
    META_ALL, sep="\t", encoding="utf-16-le",
    dtype=str, keep_default_na=False, low_memory=False,
)
df.columns = [c.lstrip("﻿") for c in df.columns]   # 第一列名带 BOM

assert df.shape == (EXPECTED_ROWS, 34), f"unexpected shape {df.shape}"
assert df["Run"].is_unique
print(f"loaded metadata_all: {df.shape}")

# %% [markdown]
# ## 3. 全列 strip + 列重命名

# %%
# 3.1 全列 strip 前后空白（实测仅 Disease_Stage 8 行有前导空格，无害普适）
for c in df.columns:
    df[c] = df[c].str.strip()

# 3.2 列重命名（与 sampleTable / 通用约定对齐）
df = df.rename(columns={
    "Smoke":                   "Smoking",
    "Recent_Antibiotics_Use":  "Recent_Antibiotic_Use",
    "Body_Site":               "Sample_Site",
    "Body_Site_Raw":           "Sample_Type",
    "age_start":               "Age_start",
    "age_end":                 "Age_end",
})
print(f"renamed cols: {list(df.columns)}")

# %% [markdown]
# ## 4. 定向值清洗（精确匹配，不污染合法值）

# %%
# 4.1 Sample_Site: Lung → Lung Tissue (1,706 行；与 sampleTable 一致)
n_lung = (df.Sample_Site == "Lung").sum()
df.loc[df.Sample_Site == "Lung", "Sample_Site"] = "Lung Tissue"
print(f"Sample_Site: Lung → Lung Tissue, {n_lung} 行")

# 4.2 Patient_ID: 伪 PID 清洗（"50.2_13.4(mean,sd)"是统计描述误塞进 PID 列；与 jxt 一致）
PSEUDO_PID = "50.2_13.4(mean,sd)"
n_pid = (df.Patient_ID == PSEUDO_PID).sum()
df.loc[df.Patient_ID == PSEUDO_PID, "Patient_ID"] = ""
print(f"Patient_ID: 伪 PID '{PSEUDO_PID}' → NA, {n_pid} 行")

# 4.3 Time_Point: HTML escape 修正（&gt;48 → 48+，与 jxt 一致；仅 2 行）
n_tp = (df.Time_Point == "&gt;48").sum()
df.loc[df.Time_Point == "&gt;48", "Time_Point"] = "48+"
print(f"Time_Point: '&gt;48' → '48+', {n_tp} 行")

# %% [markdown]
# ## 5. dtype 转换（空字符串 → NA / NaN）

# %%
# 5.1 数值列 → float64（'' 自动 NaN）
for c in FLOAT_COLS:
    df[c] = pd.to_numeric(df[c].replace("", np.nan), errors="raise")
print(f"float cols converted: {FLOAT_COLS}")

# 5.2 string 列：'' → pd.NA，dtype=string[python]
for c in STRING_COLS:
    df[c] = df[c].replace("", pd.NA).astype("string")

# 5.3 category 列（有序顺序）
for c, cats in CAT_ORDER.items():
    if c in ("Age_Group", "Case_Or_Control", "Region_16S"):
        continue   # 派生列待 §6 创建后再设 category
    s = df[c].replace("", pd.NA)
    unseen = set(s.dropna().unique()) - set(cats)
    assert not unseen, f"[{c}] 出现未知取值: {unseen}"
    df[c] = pd.Categorical(s, categories=cats)

# 5.4 category 列（无序，pandas 默认顺序）
for c in CAT_NO_ORDER:
    s = df[c].replace("", pd.NA)
    df[c] = s.astype("category")

# %% [markdown]
# ## 6. 派生列（4 个）
#
# 逻辑见 `.claude/data/metadata_resmicrodb_standardize.md` §4-§5。

# %%
# 6.1 Age_Group —— jxt 7 档（1.2_clean_phyloseq.Rmd:343-356）
def _age_group(start, end):
    if pd.isna(start) or pd.isna(end):
        return pd.NA
    if start >= 0  and end <= 3:  return "0-3"
    if start >  3  and end <= 18: return "3-18"
    if start > 18  and end <= 35: return "18-35"
    if start > 35  and end <= 45: return "35-45"
    if start > 45  and end <= 60: return "45-60"
    if start > 60  and end <= 75: return "60-75"
    if start > 75:                return "75+"
    return pd.NA   # 区间跨桶（如 (0,18) / (18,100)）→ NA

age_grp = df.apply(lambda r: _age_group(r.Age_start, r.Age_end), axis=1)
df["Age_Group"] = pd.Categorical(age_grp, categories=CAT_ORDER["Age_Group"])
print(f"Age_Group 分布:\n{df.Age_Group.value_counts(dropna=False)}")

# 6.2 Case_Or_Control —— jxt 1.2_clean_phyloseq.Rmd:365-367
def _case_or_control(p):
    if pd.isna(p):                      return pd.NA
    if p in ("Control", "Health"):      return "control"
    return "case"

coc = df["Phenotype"].map(_case_or_control)
df["Case_Or_Control"] = pd.Categorical(coc, categories=CAT_ORDER["Case_Or_Control"])
print(f"Case_Or_Control 分布:\n{df.Case_Or_Control.value_counts(dropna=False)}")

# 6.3 Is_Healthy —— jxt 1.2_clean_phyloseq.Rmd:368-371（三态 nullable bool）
def _is_healthy(p):
    if p == "Health":              return True
    if pd.isna(p):                 return pd.NA
    if p == "Control":             return pd.NA
    return False

ih = df["Phenotype"].map(_is_healthy)
df["Is_Healthy"] = pd.array(ih.tolist(), dtype="boolean")
print(f"Is_Healthy 分布:\n{df.Is_Healthy.value_counts(dropna=False)}")

# 6.4 Region_16S —— 从 projectTable 按 PID + (ma.Sequencing_Type=='16S') 强制约束
pt = pd.read_csv(PROJ_TABLE, sep="\t", dtype=str, keep_default_na=False, low_memory=False)
pt_16s = pt[pt.Sequencing_Type.str.contains("16S")].copy()  # 含 '16S;Metagenomics' 复合
assert pt_16s.Project_ID.is_unique  # 同 PID 在 pt 内唯一一行
pid2region = pt_16s.set_index("Project_ID")["16S_Region"].to_dict()
print(f"projectTable 中 16S 类项目（含复合）数: {len(pid2region)}")

def _region(row):
    if row.Sequencing_Type != "16S":   # 强制：非 16S 样本一律 NA
        return pd.NA
    region = pid2region.get(str(row.Project_ID), "-")
    if region == "-":
        return pd.NA
    return region

regions = df.apply(_region, axis=1)
unseen = set(regions.dropna().unique()) - set(CAT_ORDER["Region_16S"])
assert not unseen, f"Region_16S 出现未知取值: {unseen}"
df["Region_16S"] = pd.Categorical(regions, categories=CAT_ORDER["Region_16S"])
print(f"Region_16S 分布:\n{df.Region_16S.value_counts(dropna=False)}")

# %% [markdown]
# ## 7. 删冗余列、按最终列序排列

# %%
df = df.drop(columns=["Age_With_Interval", "BMI_With_Interval"])
df = df[FINAL_COLUMNS]   # 36 列最终顺序
assert df.shape == (EXPECTED_ROWS, EXPECTED_COLS), f"final shape {df.shape} != ({EXPECTED_ROWS},{EXPECTED_COLS})"
print(f"final shape: {df.shape}")
print(f"\ndtypes:\n{df.dtypes}")

# %% [markdown]
# ## 8. Sanity asserts（14 条）

# %%
# 1. Run 唯一非空
assert df.Run.is_unique and df.Run.notna().all()

# 2. Age 标量 invariant：Age 非空 → Age == Age_start == Age_end
mask = df.Age.notna()
assert ((df.Age[mask] == df.Age_start[mask]) & (df.Age[mask] == df.Age_end[mask])).all(), \
    "Age 标量 invariant 失败"
print(f"  ✓ Age 标量 invariant ({mask.sum()} 行)")

# 3. Age 区间 invariant：Age 空 + start 非空 → start <= end
mask = df.Age.isna() & df.Age_start.notna()
assert (df.Age_start[mask] <= df.Age_end[mask]).all(), "Age 区间 invariant 失败"
print(f"  ✓ Age 区间 invariant ({mask.sum()} 行)")

# 4-5. BMI 同样两条
mask = df.BMI.notna()
assert ((df.BMI[mask] == df.BMI_start[mask]) & (df.BMI[mask] == df.BMI_end[mask])).all()
print(f"  ✓ BMI 标量 invariant ({mask.sum()} 行)")
mask = df.BMI.isna() & df.BMI_start.notna()
assert (df.BMI_start[mask] <= df.BMI_end[mask]).all()
print(f"  ✓ BMI 区间 invariant ({mask.sum()} 行)")

# 6. Country ↔ Continent 同步空（实测全表 0 不一致）
assert (df.Country.isna() == df.Continent.isna()).all(), \
    "Country/Continent 空非同步"
print(f"  ✓ Country ↔ Continent 同步空")

# 7. Location ⊆ Country
mask = df.Location.notna() & df.Country.isna()
assert mask.sum() == 0, f"{mask.sum()} 行 Location 有值但 Country 空"
print(f"  ✓ Location ⊆ Country")

# 8. Region_16S 仅对 16S 样本非空
mask = df.Region_16S.notna() & (df.Sequencing_Type != "16S")
assert mask.sum() == 0, f"{mask.sum()} 行 Region_16S 非空但非 16S"
print(f"  ✓ Region_16S 仅 16S 样本（{df.Region_16S.notna().sum()} 行有值）")

# 9. Region_16S categories 完整 16 类
assert set(df.Region_16S.cat.categories) == set(CAT_ORDER["Region_16S"])

# 10. Age_Group categories 7 类
assert list(df.Age_Group.cat.categories) == CAT_ORDER["Age_Group"]

# 11. Case_Or_Control categories
assert list(df.Case_Or_Control.cat.categories) == CAT_ORDER["Case_Or_Control"]

# 12. Is_Healthy 三态：dtype 是 nullable boolean
assert df.Is_Healthy.dtype == "boolean"

# 13. Is_Healthy 语义：True 当且仅当 Phenotype == 'Health'
n_health = (df.Phenotype == "Health").sum()
assert (df.Is_Healthy == True).sum() == n_health
print(f"  ✓ Is_Healthy: {n_health} 行 True")

# 14. dtype 整体一致性
assert all(df[c].dtype.name == "float64" for c in FLOAT_COLS), "float 列 dtype 不一致"
assert all(df[c].dtype.name == "string"  for c in STRING_COLS), "string 列 dtype 不一致"

print("\n=== 14 条 sanity assert 全通过 ===")

# %% [markdown]
# ## 9. 跨 sampleTable 验证 —— 派生列与 jxt 输出在交集 Run 上逐行一致

# %%
st = pd.read_csv(SAMPLE_TBL, sep="\t", dtype=str, keep_default_na=False, low_memory=False)
common = set(df.Run) & set(st.Run)
assert len(common) == 106464

# 把两端按相同 Run 顺序拉齐，重置成 RangeIndex 后再做向量比较
ours = (df[df.Run.isin(common)]
        .sort_values("Run")
        .reset_index(drop=True))
jxt  = (st[st.Run.isin(common)]
        .sort_values("Run")
        .reset_index(drop=True))
assert (ours.Run.values == jxt.Run.values).all()

# 9.1 Age_Group 一致性
ours_ag = ours.Age_Group.astype(object).where(ours.Age_Group.notna(), "-").values
jxt_ag  = jxt.Age_Group.values   # jxt 用 '-' 标 NA
diff = (ours_ag != jxt_ag).sum()
assert diff == 0, f"Age_Group 与 jxt 不一致: {diff} 行"
print(f"  ✓ Age_Group 与 jxt sampleTable 在 {len(common)} 行交集上完全一致")

# 9.2 Region_16S 一致性（jxt sampleTable 错填了 1,079 行 metag/metat 的 region；我们应少这部分）
ours_rg = ours.Region_16S.astype(object).where(ours.Region_16S.notna(), "-").values
jxt_rg  = jxt["16S Region"].values
both_filled = (ours_rg != "-") & (jxt_rg != "-")
diff_filled = ((ours_rg != jxt_rg) & both_filled).sum()
assert diff_filled == 0, f"Region_16S 在两边都填值的行上有 {diff_filled} 行不一致"

ours_na_jxt_filled = ((ours_rg == "-") & (jxt_rg != "-")).sum()
print(f"  ✓ Region_16S 取值一致；其中我们 NA / jxt 填值: {ours_na_jxt_filled} 行（jxt 错填给 metag/metat 的样本，符合预期）")

# 9.3 Patient_ID 清洗一致性（伪 PID 112 行）
mask = (jxt.Patient_ID == "-").values
ours_na = ours.Patient_ID.isna().values
assert (ours_na | (~mask)).all(), "jxt 标 '-' 的 Patient_ID 行在我们端也应为 NA"
print(f"  ✓ Patient_ID 清洗与 jxt 一致")

# %% [markdown]
# ## 10. 写盘（TSV + Parquet）

# %%
OUT_TSV.parent.mkdir(parents=True, exist_ok=True)

# 10.1 TSV: NA 写空字符串；UTF-8 + LF
df.to_csv(OUT_TSV, sep="\t", index=False, na_rep="",
          encoding="utf-8", lineterminator="\n")
sz_tsv = OUT_TSV.stat().st_size
print(f"  ✓ wrote {OUT_TSV.name}  {sz_tsv/1024/1024:.1f} MB")

# 10.2 Parquet: 完整保留 dtype（含 Categorical 顺序、nullable string、nullable bool）
df.to_parquet(OUT_PARQUET, engine="pyarrow", compression="snappy", index=False)
sz_pq = OUT_PARQUET.stat().st_size
print(f"  ✓ wrote {OUT_PARQUET.name}  {sz_pq/1024/1024:.1f} MB")

# 10.3 round-trip 校验：从 parquet 重读，与内存中 df 完全一致
df_rt = pd.read_parquet(OUT_PARQUET)
assert df_rt.shape == df.shape
for c in df.columns:
    assert df_rt[c].dtype == df[c].dtype, f"[{c}] dtype 不一致: {df_rt[c].dtype} vs {df[c].dtype}"
    # 用 pd.Series.equals（含 NA 等价）
    assert df_rt[c].equals(df[c]), f"[{c}] 内容不一致"
print(f"  ✓ parquet round-trip 一致")

# %% [markdown]
# ## 11. 摘要报告

# %%
print("\n=== ResMicroDb 06 标准化完成 ===\n")
print(f"输入  metadata_all.txt   {EXPECTED_ROWS:>7,} × 34")
print(f"输出  .standardized.tsv  {len(df):>7,} × {df.shape[1]}    {sz_tsv/1024/1024:.1f} MB")
print(f"输出  .standardized.parquet  同上              {sz_pq/1024/1024:.1f} MB")
print()
print(f"全表非空率（前 12 列）:")
for c in FINAL_COLUMNS[:12]:
    nn = df[c].notna().sum()
    print(f"  {c:25s} {nn:>7,} / {len(df):,}  ({nn/len(df)*100:5.1f}%)")
print(f"  ...")
print()
print(f"派生列分布:")
for c in ["Age_Group", "Case_Or_Control", "Is_Healthy", "Region_16S"]:
    vc = df[c].value_counts(dropna=False)
    print(f"  {c}:")
    for k, v in vc.items():
        print(f"    {str(k):25s} {v:>7,}")
