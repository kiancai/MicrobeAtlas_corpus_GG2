# ---
# jupyter:
#   jupytext:
#     formats: py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#   kernelspec:
#     display_name: baseBio
#     language: python
#     name: python3
# ---

from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad

ROOT = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
ANN_IN = ROOT / "results/feature_table/merged.gg2.with_phylo.h5ad"
OUT_DIR = ROOT / "results/sample_distance_100k"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_TSV = OUT_DIR / "subset_100k_index.tsv"

RANDOM_SEED = 42
N_PAIRED_RUNS = 10_000
RM_ADDITIONAL = 10_000

MA_ADDITIONAL_QUOTA = {
    "Human": 22_000,
    "Animal_other": 12_000,
    "Soil": 13_000,
    "Aquatic": 12_000,
    "Plant": 6_000,
    "Unknown": 5_000,
}

OUT_COLS = [
    "obs_name",
    "Database",
    "Run",
    "Project_ID",
    "sample_role",
    "overview_bucket",
    "ma_bucket_detail",
    "human_site",
    "rm_site",
    "paired_run_id",
    "stratum_id",
    "sub_stratum",
]


def clean_text(values: pd.Series) -> pd.Series:
    s = values.astype(str)
    return s.replace({"nan": "NA", "None": "NA", "<NA>": "NA"})


def sqrt_allocate(group_sizes: pd.Series, total: int) -> pd.Series:
    """Allocate total counts by sqrt(N), respecting group capacity."""
    group_sizes = group_sizes.astype(int)
    if total <= 0 or len(group_sizes) == 0:
        return pd.Series(0, index=group_sizes.index, dtype=int)
    if total >= int(group_sizes.sum()):
        return group_sizes.astype(int)

    weights = np.sqrt(group_sizes.astype(float).clip(lower=1))
    raw = weights / weights.sum() * total
    out = np.floor(raw).astype(int)
    frac = (raw - out).sort_values(ascending=False)
    remainder = total - int(out.sum())
    if remainder > 0:
        out.loc[frac.head(remainder).index] += 1

    cap = group_sizes.reindex(out.index)
    while True:
        overflow = out > cap
        if not overflow.any():
            break
        excess = int((out[overflow] - cap[overflow]).sum())
        out.loc[overflow] = cap.loc[overflow]
        available = out.index[out < cap]
        if excess <= 0 or len(available) == 0:
            break
        ranked = frac.reindex(available).sort_values(ascending=False).index
        for idx in ranked:
            if excess <= 0:
                break
            add = min(excess, int(cap.loc[idx] - out.loc[idx]))
            out.loc[idx] += add
            excess -= add
    return out.astype(int)


def collapse_to_top_k(values: pd.Series, k: int = 10) -> pd.Series:
    s = clean_text(values)
    top = s.value_counts().head(k).index
    return s.where(s.isin(top), "Other")


def sample_by_quota(df: pd.DataFrame, col: str, quotas: pd.Series,
                    rng: np.random.Generator) -> pd.DataFrame:
    picked = []
    for key, quota in quotas.items():
        if quota <= 0:
            continue
        sub = df[df[col] == key]
        if len(sub) == 0:
            continue
        take = min(int(quota), len(sub))
        idx = rng.choice(sub.index.to_numpy(), size=take, replace=False)
        picked.append(sub.loc[idx].copy())
    if not picked:
        return df.iloc[0:0].copy()
    return pd.concat(picked, axis=0)


def sample_two_level(df: pd.DataFrame, level1: str, level2: str, total: int,
                     rng: np.random.Generator) -> pd.DataFrame:
    l1_sizes = df[level1].value_counts()
    l1_quotas = sqrt_allocate(l1_sizes, total)
    picked = []
    for l1, q1 in l1_quotas.items():
        sub1 = df[df[level1] == l1]
        if q1 <= 0 or len(sub1) == 0:
            continue
        l2_sizes = sub1[level2].value_counts()
        l2_quotas = sqrt_allocate(l2_sizes, int(q1))
        part = sample_by_quota(sub1, level2, l2_quotas, rng)
        picked.append(part)
    out = pd.concat(picked, axis=0) if picked else df.iloc[0:0].copy()
    assert len(out) == min(total, len(df)), f"sample_two_level got {len(out)} != requested {total}"
    return out


def add_ma_fields(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    animal = out["MA_Env_Animal"].fillna(False).astype(bool).to_numpy()
    soil = out["MA_Env_Soil"].fillna(False).astype(bool).to_numpy()
    aquatic = out["MA_Env_Aquatic"].fillna(False).astype(bool).to_numpy()
    plant = out["MA_Env_Plant"].fillna(False).astype(bool).to_numpy()
    is_human = (out["MA_IsHuman"].astype(str) == "Human").to_numpy()

    detail = np.array(["Unknown"] * len(out), dtype=object)
    detail[plant] = "Plant"
    detail[aquatic] = "Aquatic"
    detail[soil] = "Soil"
    detail[animal] = "Animal_other"
    detail[is_human] = "Human"

    overview = detail.copy()
    overview[np.isin(detail, ["Human", "Animal_other"])] = "Animal"

    out["ma_bucket_detail"] = detail
    out["overview_bucket"] = overview
    out["human_site"] = clean_text(out["MA_SampleSite"]).where(detail == "Human", "NA")
    out["rm_site"] = "NA"
    return out


def add_rm_fields(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["ma_bucket_detail"] = "NA"
    out["overview_bucket"] = "RM"
    out["human_site"] = "NA"
    out["rm_site"] = clean_text(out["RM_Sample_Site"])
    return out


def finalize_rows(df: pd.DataFrame, sample_role: str, paired_run_id="NA") -> pd.DataFrame:
    out = df.copy()
    out["sample_role"] = sample_role
    if isinstance(paired_run_id, str):
        out["paired_run_id"] = paired_run_id
    else:
        out["paired_run_id"] = paired_run_id.astype(str)
    for col in OUT_COLS:
        if col not in out.columns:
            out[col] = "NA"
    return out[OUT_COLS]


print(f"Reading obs from {ANN_IN.name} ...")
adata = ad.read_h5ad(ANN_IN, backed="r")
obs = adata.obs[[
    "Database", "Run", "Project_ID",
    "MA_IsHuman", "MA_SampleSite",
    "MA_Env_Animal", "MA_Env_Animal_Sub",
    "MA_Env_Soil", "MA_Env_Soil_Sub",
    "MA_Env_Aquatic", "MA_Env_Aquatic_Sub",
    "MA_Env_Plant", "MA_Env_Plant_Sub",
    "RM_Sample_Site",
]].copy()
obs["obs_name"] = adata.obs_names.astype(str)
obs["Run"] = obs["Run"].astype(str)
obs["Project_ID"] = clean_text(obs["Project_ID"])
obs["Database"] = obs["Database"].astype(str)
print(f"obs: {obs.shape}")
print(obs["Database"].value_counts().to_string())

rng = np.random.default_rng(RANDOM_SEED)
ma = add_ma_fields(obs[obs["Database"] == "MicrobeAtlas"].copy())
rm = add_rm_fields(obs[obs["Database"] == "ResMicroDb"].copy())

# Paired block: choose paired Runs by RM site and Project_ID, then include both rows.
run_db = obs.groupby("Run", observed=True)["Database"].nunique()
paired_runs = set(run_db[run_db >= 2].index.astype(str))
print(f"\nGlobal paired Runs: {len(paired_runs):,}")

ma_pair = ma[ma["Run"].isin(paired_runs)][["Run", "obs_name"]].rename(columns={"obs_name": "ma_obs_name"})
rm_pair = rm[rm["Run"].isin(paired_runs)][["Run", "obs_name", "rm_site", "Project_ID"]].rename(columns={"obs_name": "rm_obs_name"})
pair_candidates = rm_pair.merge(ma_pair, on="Run", how="inner")
assert len(pair_candidates) >= N_PAIRED_RUNS, "Not enough paired Runs for requested paired block"
assert pair_candidates["Run"].is_unique, "Expected one MA and one RM row per paired Run"
pair_candidates["_site"] = pair_candidates["rm_site"]
pair_candidates["_project"] = pair_candidates["Project_ID"]
selected_pairs = sample_two_level(pair_candidates, "_site", "_project", N_PAIRED_RUNS, rng)
selected_pair_runs = set(selected_pairs["Run"].astype(str))
print(f"Selected paired Runs: {len(selected_pair_runs):,}")

paired_ma = ma[ma["Run"].isin(selected_pair_runs)].copy()
paired_rm = rm[rm["Run"].isin(selected_pair_runs)].copy()
paired_ma["stratum_id"] = "paired::MA"
paired_rm["stratum_id"] = "paired::RM"
paired_ma["sub_stratum"] = paired_ma["Run"].astype(str)
paired_rm["sub_stratum"] = paired_rm["Run"].astype(str)
paired_ma_out = finalize_rows(paired_ma, "paired_ma", paired_ma["Run"])
paired_rm_out = finalize_rows(paired_rm, "paired_rm", paired_rm["Run"])
paired_obs_names = set(pd.concat([paired_ma_out["obs_name"], paired_rm_out["obs_name"]]))

# Additional MA block.
ma_remaining = ma[~ma["obs_name"].isin(paired_obs_names)].copy()
ma_parts = []
for bucket, quota in MA_ADDITIONAL_QUOTA.items():
    sub = ma_remaining[ma_remaining["ma_bucket_detail"] == bucket].copy()
    if bucket == "Human":
        sub["_l2"] = clean_text(sub["MA_SampleSite"])
    elif bucket == "Animal_other":
        sub["_l2"] = collapse_to_top_k(sub["MA_Env_Animal_Sub"], k=10)
    elif bucket == "Soil":
        sub["_l2"] = collapse_to_top_k(sub["MA_Env_Soil_Sub"], k=10)
    elif bucket == "Aquatic":
        sub["_l2"] = collapse_to_top_k(sub["MA_Env_Aquatic_Sub"], k=10)
    elif bucket == "Plant":
        sub["_l2"] = collapse_to_top_k(sub["MA_Env_Plant_Sub"], k=10)
    else:
        sub["_l2"] = "NA"
    quotas = sqrt_allocate(sub["_l2"].value_counts(), quota)
    picked = sample_by_quota(sub, "_l2", quotas, rng)
    picked["stratum_id"] = f"MA::{bucket}"
    picked["sub_stratum"] = picked["_l2"].astype(str)
    ma_parts.append(picked)
    print(f"MA::{bucket:<13} target={quota:>6} picked={len(picked):>6}")
ma_add = pd.concat(ma_parts, axis=0)
assert len(ma_add) == sum(MA_ADDITIONAL_QUOTA.values())
ma_add_out = finalize_rows(ma_add, "ma_additional")

# Additional RM block.
rm_remaining = rm[~rm["obs_name"].isin(paired_obs_names)].copy()
rm_remaining["_site"] = rm_remaining["rm_site"]
rm_remaining["_project"] = rm_remaining["Project_ID"]
rm_add = sample_two_level(rm_remaining, "_site", "_project", RM_ADDITIONAL, rng)
rm_add["stratum_id"] = "RM::additional"
rm_add["sub_stratum"] = rm_add["_site"].astype(str) + "|" + rm_add["_project"].astype(str)
rm_add_out = finalize_rows(rm_add, "rm_additional")
print(f"RM additional target={RM_ADDITIONAL:,} picked={len(rm_add_out):,}")

subset = pd.concat([paired_ma_out, paired_rm_out, ma_add_out, rm_add_out], axis=0).reset_index(drop=True)
assert len(subset) == 100_000, f"subset has {len(subset)} rows, expected 100000"
assert subset["obs_name"].is_unique, "subset obs_name is not unique"

print("\nFinal sample_role:")
print(subset["sample_role"].value_counts().to_string())
print("\nDatabase:")
print(subset["Database"].value_counts().to_string())
print("\nMA overview bucket:")
print(subset.loc[subset["Database"] == "MicrobeAtlas", "overview_bucket"].value_counts().to_string())
print("\nHuman sites:")
print(subset.loc[subset["ma_bucket_detail"] == "Human", "human_site"].value_counts().to_string())
print("\nRM sites:")
print(subset.loc[subset["Database"] == "ResMicroDb", "rm_site"].value_counts().to_string())
print("\nPaired rows:")
paired = subset[subset["sample_role"].isin(["paired_ma", "paired_rm"])]
print(paired.groupby(["paired_run_id", "Database"]).size().unstack(fill_value=0).shape)

subset.to_csv(OUT_TSV, sep="\t", index=False)
print(f"\nWrote {OUT_TSV}")
