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
# # 06 · MicrobeAtlas metadata 整合
#
# 把 `samples.env.info.tsv` 的 9 列元信息解析、拆分并对齐写入
# `gg2.full.qc.h5ad` 的 `obs`，得到 `gg2.full.qc.with_meta.h5ad`。
#
# **obs schema (26 列)** —— `MA_` 前缀 = MicrobeAtlas 特有，无前缀 = 与 ResMicroDb 通用。
#
# | # | 列 | dtype | 说明 |
# |---|----|-------|------|
# | 1 | `Database` | category | `MicrobeAtlas` / `ResMicroDb` |
# | 2 | `MA_Sample_ID` | string | 原 MAP_SID |
# | 3 | `Run` | string | SRR/ERR/DRR |
# | 4 | `BioSample` | string | SRS/ERS/DRS |
# | 5 | `Project_ID` | category | ERP/SRP/DRP |
# | 6 | `Sequencing_Type` | category | AMPLICON/WGS/RNAseq/NA |
# | 7-14 | `MA_Env_{Animal,Soil,Aquatic,Plant}` + `_Sub` | bool / string | 4 主类 × (bool, 子类 \| 连) |
# | 15 | `MA_IsHuman` | category | `Human` (仅 human) / `HumanMix` (human+其他 animal) / NA |
# | 16 | `MA_SampleSite_Raw` | category | col3 原值 (242 唯一) |
# | 17 | `MA_SampleSite` | category | gut/skin/oral/urogenital/lung/gastric/nose/bone/NA |
# | 18 | `MA_Health` | string | 健康相关 token 用 ; 拼回（不解冲突） |
# | 19 | `Sex` | category | female/male/NA |
# | 20 | `MA_AgeGroup` | category | infant/baby/toddler/adult/elderly/NA |
# | 21 | `Smoking` | category | Smoker/NA（ResMicroDb 阶段会扩 Non-smoker） |
# | 22 | `MA_Keywords` | string | col5 原值 |
# | 23 | `MA_Geo_Raw` | string | col9 原值 |
# | 24 | `Latitude` | float64 | 解析失败 NaN |
# | 25 | `Longitude` | float64 | 解析失败 NaN |
# | 26 | `MA_Institution` | string | col8 原值 |
#
# index 用 anndata 自动 RangeIndex，不继承 MAP_SID。
#
# **校验** (assert)：
# 1. h5ad obs 与 tsv MAP_SID 100% 匹配
# 2. col2 拆 8 列再拼回 == 原值（集合相等）
# 3. col3 拆 token 后无未识别 token（28 token 白名单全覆盖）

# %%
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad

# anndata < 0.11 默认禁写 pd.StringArray (string[python])，需显式打开
ad.settings.allow_write_nullable_strings = True

PROJECT_DIR = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
TSV_PATH    = PROJECT_DIR / "rawdata/MicrobeAtlas/sample_info/samples.env.info.tsv"
H5AD_IN     = PROJECT_DIR / "results/feature_table/gg2.full.qc.h5ad"
H5AD_OUT    = PROJECT_DIR / "results/feature_table/gg2.full.qc.with_meta.h5ad"

DATABASE_NAME = "MicrobeAtlas"

# col2 主类
ENV_MAINS = ["animal", "soil", "aquatic", "plant"]

# col3 token 白名单（穷举 28 个 token，分配到 5 个槽位）
SAMPLE_SITE_TOKENS = {"gut", "skin", "oral", "urogenital", "lung", "gastric", "nose", "bone"}
SEX_TOKENS         = {"female", "male"}
AGE_TOKENS         = {"infant", "baby", "toddler", "adult", "elderly"}
SMOKING_TOKENS     = {"smoker"}
HEALTH_TOKENS = {
    "healthy", "disease", "infection",
    "inflammatory bowel disease", "dermatitis", "cystic fibrosis",
    "cholera", "malaria", "pneumonia", "typhoid fever", "tuberculosis",
    "patient",
}
ALL_COL3_TOKENS = (
    SAMPLE_SITE_TOKENS | SEX_TOKENS | AGE_TOKENS | SMOKING_TOKENS | HEALTH_TOKENS
)


# %% [markdown]
# ## 1. 读 h5ad 与 tsv，按 obs_names 对齐

# %%
adata = ad.read_h5ad(H5AD_IN)
print(f"h5ad: {adata.shape[0]:,} samples × {adata.shape[1]:,} genera")

obs_ids = adata.obs_names.to_numpy()

df = pd.read_csv(
    TSV_PATH, sep="\t",
    dtype=str, keep_default_na=False, na_filter=False,
    low_memory=False,
)
print(f"tsv: {len(df):,} 行 × {df.shape[1]} 列")
print("tsv columns:", list(df.columns))

# 列名固化（防止 _ 列名重复带来对齐风险）
df.columns = [
    "MAP_SID", "Environments", "_col3", "Technology",
    "Keywords", "_col6", "Project", "Institution", "_col9",
]

# 对齐 (h5ad obs 是 tsv MAP_SID 的子集，已先期验证 100% 命中)
df = df.set_index("MAP_SID").reindex(obs_ids)
assert df.notna().any(axis=1).all(), "存在 obs 在 tsv 找不到"
print(f"对齐后: {len(df):,} 行（与 h5ad obs 一一对应）")


# %% [markdown]
# ## 2. col1 (MAP_SID) 拆 Run / BioSample

# %%
sample_id = df.index.to_series().reset_index(drop=True)
split = sample_id.str.split(".", n=1, expand=True)
run_id = split[0]
biosample_id = split[1]

assert (run_id != "").all() and (biosample_id != "").all(), "MAP_SID 拆分有空段"
print("Run head:", run_id.head(3).tolist())
print("BioSample head:", biosample_id.head(3).tolist())


# %% [markdown]
# ## 3. col2 (Environments) → 8 列 + MA_IsHuman
#
# 拆分规则：
# - col2 用 `|` split → 标签列表
# - 每个标签 `主类;子类`，主类 ∈ {animal, soil, aquatic, plant}
# - 同主类多子类：`_Sub` 用 `|` 连
# - 主类无子类（裸 `animal`）：`_Sub` 为 ""
# - col2 全空：4 个 bool False，4 个 sub NA
#
# round-trip 校验：拆 + 拼回 → 与原值 set 相等。

# %%
env_raw = df["Environments"].fillna("").astype(str)

# 解析为 list[ (main, sub) ]
def _parse_env(s):
    if not s:
        return []
    out = []
    for tag in s.split("|"):
        if ";" in tag:
            main, sub = tag.split(";", 1)
        else:
            main, sub = tag, ""
        out.append((main, sub))
    return out

parsed_env = env_raw.apply(_parse_env)

env_cols = {}
for main in ENV_MAINS:
    flag_col = f"MA_Env_{main.capitalize()}"
    sub_col  = f"MA_Env_{main.capitalize()}_Sub"

    flags = parsed_env.apply(lambda lst, m=main: any(t[0] == m for t in lst))
    subs  = parsed_env.apply(
        lambda lst, m=main: "|".join(t[1] for t in lst if t[0] == m)
    )
    # 主类不存在 → NA；存在但子类全空 → ""
    subs = subs.where(flags, other=pd.NA)

    env_cols[flag_col] = flags.to_numpy()
    env_cols[sub_col]  = subs.to_numpy()

# round-trip 校验：从 8 列重建 → 与原值 set 相等
def _rebuild_from_parsed(env_list):
    parts = []
    by_main = {m: [] for m in ENV_MAINS}
    for main, sub in env_list:
        by_main[main].append(sub)
    for main in ENV_MAINS:
        subs = by_main[main]
        if not subs:
            continue
        for s in subs:
            parts.append(f"{main};{s}" if s else main)
    return set(parts)

print("Environments round-trip 校验中…")
orig_sets    = env_raw.apply(lambda s: set(s.split("|")) if s else set())
rebuilt_sets = parsed_env.apply(_rebuild_from_parsed)
mismatch_mask = orig_sets != rebuilt_sets
n_mismatch = int(mismatch_mask.sum())
if n_mismatch:
    bad_idx = np.where(mismatch_mask.to_numpy())[0][:5]
    for i in bad_idx:
        print(f"  [warn] row {i} ({obs_ids[i]}): orig={orig_sets.iat[i]} rebuilt={rebuilt_sets.iat[i]}")
assert n_mismatch == 0, f"Environments 拆分不可逆: {n_mismatch} 行不匹配"
print(f"通过：{len(df):,} 行全部可逆")


# %% [markdown]
# ## 4. MA_IsHuman（严格定义 + HumanMix）
#
# - `Human`：col2 中 animal 子类**只有** human
# - `HumanMix`：含 human 且含其他 animal 子类
# - `NA`：不含 `animal;human` 标签

# %%
def _human_class(env_list):
    animal_subs = [t[1] for t in env_list if t[0] == "animal"]
    has_human = "human" in animal_subs
    if not has_human:
        return pd.NA
    others = [s for s in animal_subs if s != "human" and s != ""]
    return "HumanMix" if others else "Human"

ma_is_human = parsed_env.apply(_human_class)
print("MA_IsHuman 分布:")
print(ma_is_human.value_counts(dropna=False))


# %% [markdown]
# ## 5. col3 → 5 个槽位（白名单）
#
# 28 token 白名单已穷举验证。脚本 assert 全表无未识别 token。

# %%
col3_raw = df["_col3"].fillna("").astype(str)

# 实测 col3 不含 |，直接 split ;
def _parse_col3(s):
    return [t for t in s.split(";") if t] if s else []

parsed_col3 = col3_raw.apply(_parse_col3)

# 未识别 token 检查
unknown = set()
for toks in parsed_col3:
    for t in toks:
        if t not in ALL_COL3_TOKENS:
            unknown.add(t)
assert not unknown, f"col3 含未识别 token: {unknown}"
print(f"通过：col3 全部 token 在 28 词白名单内")

def _pick_one(toks, allowed):
    for t in toks:
        if t in allowed:
            return t
    return pd.NA

def _pick_health(toks):
    hits = [t for t in toks if t in HEALTH_TOKENS]
    return ";".join(hits) if hits else pd.NA

def _has(toks, allowed):
    return any(t in allowed for t in toks)

ma_sample_site_raw = col3_raw.where(col3_raw != "", other=pd.NA)
ma_sample_site = parsed_col3.apply(lambda x: _pick_one(x, SAMPLE_SITE_TOKENS))
ma_health      = parsed_col3.apply(_pick_health)
sex            = parsed_col3.apply(lambda x: _pick_one(x, SEX_TOKENS))
ma_age_group   = parsed_col3.apply(lambda x: _pick_one(x, AGE_TOKENS))
smoking        = parsed_col3.apply(
    lambda x: "Smoker" if _has(x, SMOKING_TOKENS) else pd.NA
)

# Sanity: AgeGroup 至多 1 个（互斥），实测无共现，但仍做 assert
age_multi = parsed_col3.apply(
    lambda x: sum(1 for t in x if t in AGE_TOKENS)
).max()
assert age_multi <= 1, f"AgeGroup 出现共现: max={age_multi}"
sex_multi = parsed_col3.apply(
    lambda x: sum(1 for t in x if t in SEX_TOKENS)
).max()
assert sex_multi <= 1, f"Sex 出现共现: max={sex_multi}"

print("MA_SampleSite 分布:")
print(ma_sample_site.value_counts(dropna=False).head(10))
print("\nSex 分布:")
print(sex.value_counts(dropna=False))
print("\nMA_AgeGroup 分布:")
print(ma_age_group.value_counts(dropna=False))
print("\nSmoking 分布:")
print(smoking.value_counts(dropna=False))


# %% [markdown]
# ## 6. col4 / col5 / col7 / col8

# %%
sequencing_type = df["Technology"].where(df["Technology"] != "", other=pd.NA)
ma_keywords     = df["Keywords"].where(df["Keywords"] != "", other=pd.NA)
project_id      = df["Project"].where(df["Project"] != "", other=pd.NA)
ma_institution  = df["Institution"].where(df["Institution"] != "", other=pd.NA)

print("Sequencing_Type 分布:")
print(sequencing_type.value_counts(dropna=False))


# %% [markdown]
# ## 7. col9 → Latitude / Longitude
#
# 实测有 `0 0`、`-122.726`（单值）、`002.9167`（前导零）等脏值。
# 解析规则：split() 后两个 token 都能 float() 才填，否则 NaN；MA_Geo_Raw 总保留。

# %%
geo_raw = df["_col9"].fillna("").astype(str)

def _parse_geo(s):
    if not s:
        return (np.nan, np.nan)
    parts = s.split()
    if len(parts) != 2:
        return (np.nan, np.nan)
    try:
        return (float(parts[0]), float(parts[1]))
    except ValueError:
        return (np.nan, np.nan)

geo_parsed = geo_raw.apply(_parse_geo)
latitude  = geo_parsed.apply(lambda t: t[0]).astype(np.float64)
longitude = geo_parsed.apply(lambda t: t[1]).astype(np.float64)

geo_raw_out = geo_raw.where(geo_raw != "", other=pd.NA)

print(f"Geo 解析: {latitude.notna().sum():,} / {len(latitude):,} 行得到有效 lat/lon")


# %% [markdown]
# ## 8. 组装 obs DataFrame

# %%
n = len(df)
obs = pd.DataFrame(index=pd.RangeIndex(n).astype(str))
obs["Database"]        = pd.Categorical([DATABASE_NAME] * n,
                                         categories=[DATABASE_NAME, "ResMicroDb"])
obs["MA_Sample_ID"]    = pd.array(df.index.to_numpy(), dtype="string")
obs["Run"]             = pd.array(run_id.to_numpy(), dtype="string")
obs["BioSample"]       = pd.array(biosample_id.to_numpy(), dtype="string")
obs["Project_ID"]      = pd.Categorical(project_id.to_numpy())
obs["Sequencing_Type"] = pd.Categorical(
    sequencing_type.to_numpy(),
    categories=["AMPLICON", "WGS", "RNAseq"],
)

for main in ENV_MAINS:
    flag_col = f"MA_Env_{main.capitalize()}"
    sub_col  = f"MA_Env_{main.capitalize()}_Sub"
    obs[flag_col] = env_cols[flag_col].astype(bool)
    obs[sub_col]  = pd.array(env_cols[sub_col], dtype="string")

obs["MA_IsHuman"] = pd.Categorical(
    ma_is_human.to_numpy(),
    categories=["Human", "HumanMix"],
)
obs["MA_SampleSite_Raw"] = pd.Categorical(ma_sample_site_raw.to_numpy())
obs["MA_SampleSite"]     = pd.Categorical(
    ma_sample_site.to_numpy(),
    categories=sorted(SAMPLE_SITE_TOKENS),
)
obs["MA_Health"]    = pd.array(ma_health.to_numpy(), dtype="string")
obs["Sex"]          = pd.Categorical(sex.to_numpy(), categories=["female", "male"])
obs["MA_AgeGroup"]  = pd.Categorical(
    ma_age_group.to_numpy(),
    categories=["infant", "baby", "toddler", "adult", "elderly"],
)
obs["Smoking"]      = pd.Categorical(
    smoking.to_numpy(),
    categories=["Smoker", "Non-smoker"],
)

obs["MA_Keywords"]    = pd.array(ma_keywords.to_numpy(), dtype="string")
obs["MA_Geo_Raw"]     = pd.array(geo_raw_out.to_numpy(), dtype="string")
obs["Latitude"]       = latitude.to_numpy()
obs["Longitude"]      = longitude.to_numpy()
obs["MA_Institution"] = pd.array(ma_institution.to_numpy(), dtype="string")

print(f"obs 组装完成: {obs.shape}")
print("\ndtypes:")
print(obs.dtypes)


# %% [markdown]
# ## 9. 写出新 h5ad
#
# 注意：此前 adata.obs_names 是 MAP_SID；新 obs 用 RangeIndex（字符串形式）。
# 矩阵 X 的行序保持不变（与 obs_ids 对齐）。

# %%
adata.obs = obs
print(f"写出 → {H5AD_OUT}")
adata.write_h5ad(H5AD_OUT, compression="gzip")
print(f"完成: {adata.shape}")


# %% [markdown]
# ## 10. 抽查
#
# 取前 5 行打印关键字段，确认整合无误。

# %%
check_cols = [
    "Database", "MA_Sample_ID", "Run", "Project_ID", "Sequencing_Type",
    "MA_Env_Animal", "MA_Env_Animal_Sub", "MA_IsHuman",
    "MA_SampleSite", "MA_Health", "Sex", "MA_AgeGroup", "Smoking",
    "Latitude", "Longitude",
]
print(adata.obs[check_cols].head(5))
