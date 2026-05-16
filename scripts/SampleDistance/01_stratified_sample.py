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
# # 01: 分层抽样选 50k 子集
#
# 输入：
# - `results/feature_table/merged.gg2.with_phylo.h5ad`  (1,826,126 × 8,114)
#
# 输出：
# - `results/sample_distance/subset_50k_index.tsv`  50,000 行：
#   - `obs_name`、`Database`、`stratum_id`、`sub_stratum`
#
# ## 配额方案（与上游讨论敲定）
#
# 总量 50,000：
# - MA 30,000、RM 20,000
#
# MA 30,000 的一级桶（按 animal > soil > aquatic > plant 主类优先级互斥归属）：
# - `Human` (IsHuman == 'Human')                       8,000   user-fixed
# - `Animal_other` (Animal True 且非 pure Human)       8,099   sqrt(N) 分配
# - `Soil`                                              5,184   sqrt(N) 分配
# - `Aquatic`                                           4,969   sqrt(N) 分配
# - `Plant`                                             2,248   sqrt(N) 分配
# - `Unknown` (无任何 env flag)                         1,500   user-fixed
#
# 二级分层：
# - `Human`：按 `MA_SampleSite` 8 部位 + NA，sqrt(N) 分配
# - `Animal_other`：按 `MA_Env_Animal_Sub` top-10 + Other，sqrt(N) 分配
# - `Soil/Aquatic/Plant`：按各自 `_Sub` top-10 + Other，sqrt(N) 分配
# - `Unknown`：桶内直接随机
#
# RM 20,000：
# - 一级 `RM_Sample_Site` 14 类，sqrt(N) 分配
# - 二级 `Project_ID`：每个项目均匀分一份（保证不被单项目主导）
#
# 跨库 Run 重复：按上游意愿不去重，两库独立抽样自然产生重叠。
#
# RANDOM_SEED = 42 全程固定。

# %%
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad

ROOT = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
ANN_IN = ROOT / "results/feature_table/merged.gg2.with_phylo.h5ad"
OUT_DIR = ROOT / "results/sample_distance"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_TSV = OUT_DIR / "subset_50k_index.tsv"

RANDOM_SEED = 42

# 一级配额
MA_QUOTA = {
    "Human":         8_000,
    "Animal_other":  8_099,
    "Soil":          5_184,
    "Aquatic":       4_969,
    "Plant":         2_248,
    "Unknown":       1_500,
}
assert sum(MA_QUOTA.values()) == 30_000, f"MA 配额合计 {sum(MA_QUOTA.values())} != 30000"
RM_TOTAL = 20_000

# %% [markdown]
# ## §1 读 obs（backed 模式 + 只取需要的列）

# %%
print(f"读 {ANN_IN.name} (backed) ...")
adata = ad.read_h5ad(ANN_IN, backed="r")
obs = adata.obs[[
    "Database", "Project_ID",
    "MA_IsHuman", "MA_SampleSite",
    "MA_Env_Animal", "MA_Env_Animal_Sub",
    "MA_Env_Soil", "MA_Env_Soil_Sub",
    "MA_Env_Aquatic", "MA_Env_Aquatic_Sub",
    "MA_Env_Plant", "MA_Env_Plant_Sub",
    "RM_Sample_Site",
]].copy()
obs["obs_name"] = adata.obs_names.astype(str)
print(f"  obs shape: {obs.shape}")
print(f"  Database: {dict(obs['Database'].value_counts())}")

# %% [markdown]
# ## §2 给 MA 行打一级 bucket（互斥优先级）

# %%
ma_mask = (obs["Database"] == "MicrobeAtlas").values
rm_mask = (obs["Database"] == "ResMicroDb").values
ma = obs.loc[ma_mask].copy()
rm = obs.loc[rm_mask].copy()
print(f"MA: {len(ma):,}   RM: {len(rm):,}")

# bool 列 NA → False
animal  = ma["MA_Env_Animal"].fillna(False).astype(bool).values
soil    = ma["MA_Env_Soil"].fillna(False).astype(bool).values
aquatic = ma["MA_Env_Aquatic"].fillna(False).astype(bool).values
plant   = ma["MA_Env_Plant"].fillna(False).astype(bool).values
is_human = (ma["MA_IsHuman"].astype(str) == "Human").values  # NA → "nan" → False

bucket = np.empty(len(ma), dtype=object)
bucket[:] = "Unknown"
# 优先级：Human > Animal_other > Soil > Aquatic > Plant > Unknown
# Human 必须先扣（is_human 已包含 animal flag，但其它 animal_other 也带 animal flag）
bucket[plant]   = "Plant"
bucket[aquatic] = "Aquatic"
bucket[soil]    = "Soil"
bucket[animal]  = "Animal_other"   # 覆盖前面，因为 animal 优先级最高（非 Human 时）
bucket[is_human] = "Human"          # Human 最高优先级
ma["bucket"] = bucket

print("\nMA 一级桶分布：")
print(ma["bucket"].value_counts())

# %% [markdown]
# ## §3 二级分层 + sqrt 配额工具函数

# %%
def sqrt_allocate(group_sizes: pd.Series, total: int) -> pd.Series:
    """sqrt(N) 加权分配 total 名额；自动 round + 调整余项。"""
    if total <= 0:
        return pd.Series(0, index=group_sizes.index, dtype=int)
    w = np.sqrt(group_sizes.astype(float).clip(lower=1))
    raw = w / w.sum() * total
    floored = np.floor(raw).astype(int)
    remainder = total - floored.sum()
    # 余项按 raw 小数部分倒序补
    frac = (raw - floored).sort_values(ascending=False)
    add_idx = frac.head(remainder).index
    floored.loc[add_idx] = floored.loc[add_idx] + 1
    # cap：超过实际样本数的截断
    cap = pd.Series(group_sizes).reindex(floored.index)
    overflow = (floored > cap)
    if overflow.any():
        # 把 overflow 部分裁掉，剩余名额从未被 overflow 的组按 frac 再补
        excess = (floored[overflow] - cap[overflow]).sum()
        floored.loc[overflow] = cap.loc[overflow]
        avail = floored.index[(~overflow) & (cap - floored > 0)]
        while excess > 0 and len(avail) > 0:
            add = min(excess, len(avail))
            top = frac.reindex(avail).sort_values(ascending=False).head(add).index
            floored.loc[top] = floored.loc[top] + 1
            # 再次检查 cap
            ov = floored.loc[top] > cap.loc[top]
            if ov.any():
                floored.loc[top[ov]] = cap.loc[top[ov]]
            excess -= add
            avail = floored.index[(floored < cap)]
    return floored.astype(int)


def collapse_to_top_k(values: pd.Series, k: int = 10, other_label: str = "Other") -> pd.Series:
    """把字符串分类列里 top-k 之外的并到 Other。"""
    s = values.astype(str)
    s = s.where(s != "nan", "NA").where(s != "None", "NA").where(s != "<NA>", "NA")
    top = s.value_counts().head(k).index
    return s.where(s.isin(top), other_label)


def stratified_sample(df: pd.DataFrame, stratum_col: str, quotas: pd.Series,
                       rng: np.random.Generator) -> pd.DataFrame:
    """按 stratum_col 抽样：每个 stratum 抽 quotas[stratum] 行；返回 df 的子视图带 sub_stratum 列。"""
    picked = []
    for sid, q in quotas.items():
        sub = df[df[stratum_col] == sid]
        if q <= 0 or len(sub) == 0:
            continue
        take = min(q, len(sub))
        idx = rng.choice(sub.index.values, size=take, replace=False)
        chunk = sub.loc[idx].copy()
        chunk["sub_stratum"] = str(sid)
        picked.append(chunk)
    if not picked:
        return df.iloc[0:0].copy().assign(sub_stratum=pd.Series(dtype=str))
    return pd.concat(picked, axis=0)


# %% [markdown]
# ## §4 抽样 MA 各 bucket

# %%
rng = np.random.default_rng(RANDOM_SEED)
ma_picked = []

# Human：按 MA_SampleSite 8 部位 + NA
sub_df = ma[ma["bucket"] == "Human"].copy()
sub_df["_l2"] = collapse_to_top_k(sub_df["MA_SampleSite"], k=20)  # 实际只有 8 部位 + NA
sizes = sub_df["_l2"].value_counts()
quotas = sqrt_allocate(sizes, MA_QUOTA["Human"])
picked = stratified_sample(sub_df, "_l2", quotas, rng)
picked["stratum_id"] = "MA::Human"
ma_picked.append(picked)
print(f"\nMA::Human  目标 {MA_QUOTA['Human']:>5} 实抽 {len(picked):>5}")
print(quotas.to_string())

# Animal_other：按 MA_Env_Animal_Sub top-10 + Other
sub_df = ma[ma["bucket"] == "Animal_other"].copy()
sub_df["_l2"] = collapse_to_top_k(sub_df["MA_Env_Animal_Sub"], k=10)
sizes = sub_df["_l2"].value_counts()
quotas = sqrt_allocate(sizes, MA_QUOTA["Animal_other"])
picked = stratified_sample(sub_df, "_l2", quotas, rng)
picked["stratum_id"] = "MA::Animal_other"
ma_picked.append(picked)
print(f"\nMA::Animal_other  目标 {MA_QUOTA['Animal_other']:>5} 实抽 {len(picked):>5}")
print(quotas.to_string())

# Soil
sub_df = ma[ma["bucket"] == "Soil"].copy()
sub_df["_l2"] = collapse_to_top_k(sub_df["MA_Env_Soil_Sub"], k=10)
sizes = sub_df["_l2"].value_counts()
quotas = sqrt_allocate(sizes, MA_QUOTA["Soil"])
picked = stratified_sample(sub_df, "_l2", quotas, rng)
picked["stratum_id"] = "MA::Soil"
ma_picked.append(picked)
print(f"\nMA::Soil  目标 {MA_QUOTA['Soil']:>5} 实抽 {len(picked):>5}")
print(quotas.to_string())

# Aquatic
sub_df = ma[ma["bucket"] == "Aquatic"].copy()
sub_df["_l2"] = collapse_to_top_k(sub_df["MA_Env_Aquatic_Sub"], k=10)
sizes = sub_df["_l2"].value_counts()
quotas = sqrt_allocate(sizes, MA_QUOTA["Aquatic"])
picked = stratified_sample(sub_df, "_l2", quotas, rng)
picked["stratum_id"] = "MA::Aquatic"
ma_picked.append(picked)
print(f"\nMA::Aquatic  目标 {MA_QUOTA['Aquatic']:>5} 实抽 {len(picked):>5}")
print(quotas.to_string())

# Plant
sub_df = ma[ma["bucket"] == "Plant"].copy()
sub_df["_l2"] = collapse_to_top_k(sub_df["MA_Env_Plant_Sub"], k=10)
sizes = sub_df["_l2"].value_counts()
quotas = sqrt_allocate(sizes, MA_QUOTA["Plant"])
picked = stratified_sample(sub_df, "_l2", quotas, rng)
picked["stratum_id"] = "MA::Plant"
ma_picked.append(picked)
print(f"\nMA::Plant  目标 {MA_QUOTA['Plant']:>5} 实抽 {len(picked):>5}")
print(quotas.to_string())

# Unknown：桶内随机（不二级分层）
sub_df = ma[ma["bucket"] == "Unknown"].copy()
take = min(MA_QUOTA["Unknown"], len(sub_df))
idx = rng.choice(sub_df.index.values, size=take, replace=False)
picked = sub_df.loc[idx].copy()
picked["sub_stratum"] = "NA"
picked["stratum_id"] = "MA::Unknown"
ma_picked.append(picked)
print(f"\nMA::Unknown  目标 {MA_QUOTA['Unknown']:>5} 实抽 {len(picked):>5}")

ma_subset = pd.concat(ma_picked, axis=0)
print(f"\nMA 子集合计: {len(ma_subset):,}")

# %% [markdown]
# ## §5 抽样 RM

# %%
sub_df = rm.copy()
sub_df["_l1"] = sub_df["RM_Sample_Site"].astype(str).fillna("NA")
sub_df["_l1"] = sub_df["_l1"].replace({"nan": "NA", "None": "NA", "<NA>": "NA"})
sizes_l1 = sub_df["_l1"].value_counts()
print("RM_Sample_Site 一级 size:")
print(sizes_l1)

l1_quotas = sqrt_allocate(sizes_l1, RM_TOTAL)
print(f"\nRM 一级配额（合计 {l1_quotas.sum()}）:")
print(l1_quotas.to_string())

rm_picked = []
for site, q in l1_quotas.items():
    sub = sub_df[sub_df["_l1"] == site]
    if q <= 0 or len(sub) == 0:
        continue
    # 二级 Project_ID：按项目均匀分（cap 单项目最多 q / sqrt(n_proj) 防止单项目主导）
    proj_sizes = sub["Project_ID"].astype(str).value_counts()
    proj_quotas = sqrt_allocate(proj_sizes, q)
    picked = stratified_sample(sub.assign(_p=sub["Project_ID"].astype(str)),
                                "_p", proj_quotas, rng)
    picked["sub_stratum"] = picked["_p"]
    picked["stratum_id"] = f"RM::{site}"
    rm_picked.append(picked)
    print(f"  RM::{site}  目标 {q:>5} 实抽 {len(picked):>5}  ({len(proj_quotas)} 项目)")

rm_subset = pd.concat(rm_picked, axis=0)
print(f"\nRM 子集合计: {len(rm_subset):,}")

# %% [markdown]
# ## §6 合并 + 写出 index

# %%
subset = pd.concat([ma_subset, rm_subset], axis=0)[
    ["obs_name", "Database", "stratum_id", "sub_stratum"]
].reset_index(drop=True)
print(f"\n子集合计: {len(subset):,}")
print(f"  MA: {(subset['Database'] == 'MicrobeAtlas').sum():,}")
print(f"  RM: {(subset['Database'] == 'ResMicroDb').sum():,}")
print(f"  唯一 obs_name: {subset['obs_name'].nunique():,} (应该 == len)")
assert subset["obs_name"].nunique() == len(subset), "subset 内部出现 obs_name 重复"

# 看 stratum_id 分布
print(f"\nstratum_id 分布:")
print(subset["stratum_id"].value_counts().to_string())

# %%
subset.to_csv(OUT_TSV, sep="\t", index=False)
print(f"\n已写出: {OUT_TSV}")
print(f"  大小: {OUT_TSV.stat().st_size / 1024:.1f} KB")
