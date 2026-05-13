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
# # ResMicroDb 06b: 应用 ps.16s_0105_new7 修正到 metadata_all
#
# 把 `ps.16s_0105_new7.rds` 揭示的 metadata 错标修正应用到
# `metadata_all.standardized.parquet`，输出 `.fixed.parquet`。
# 后续 07_attach_metadata 默认消费 fixed 版。
#
# **三个 patch**（详见 `rawdata/ResMicroDb/supplement data/CHANGES_0105_new7.md §6`）：
#
# 1. **Sample_Site**：PRJNA914884 (1,195) + PRJNA1058141 (62) 的 `Nasal` → `Nasopharynx`，共 **1,257** 行
# 2. **Phenotype + Phenotype_ID**：PRJNA801796 (Influenza 细分)，按 schema 标准词 + 真实 ontology ID，共 **255** 行
# 3. **cc + ih + Phenotype + Phenotype_ID**：错标按真实身份修正
#    - PRJNA822681 (113 行) 按 `host.disease`
#    - PRJNA824137 (16 行) 按 `real_group`（`HC*` → Health control / `TBZ*` → TB case）
#
# 总：**1,641** 个 Run 被 patch。其他 134,000+ 行不变。
#
# **依赖**：`06b_export_patches.R`（baseR + phyloseq）先跑出 patch tsv。

# %%
from pathlib import Path

import numpy as np
import pandas as pd

ROOT     = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
IN_META  = ROOT / "results/feature_table/metadata_all.standardized.parquet"
IN_PATCH = ROOT / "results/feature_table/metadata_patches_0105.tsv"
OUT_META = ROOT / "results/feature_table/metadata_all.standardized.fixed.parquet"

# Patch 2 映射: jxt 0105 字符串 → (我们目标 Phenotype, 真实 Phenotype_ID)
P2_MAPPING = {
    "Influenza A":                            ("Influenza A Virus",                     "NCIT_C53454"),
    "Influenza B":                            ("Influenza B Virus",                     "NCIT_C53468"),
    "Rhinovirus Infection":                   ("Rhinovirus Infection",                  "NCIT_C122572"),
    "Respiratory Syncytial Virus Infection":  ("Respiratory Syncytial Virus Infection", "EFO_1001413"),
}

# Patch 3a 映射: PRJNA822681 host.disease → (cc, ih, Phenotype, Phenotype_ID)
P3A_MAPPING = {
    "Healthy":  ("control", True,  "Health",   "EFO_0010130"),
    "COVID-19": ("case",    False, "COVID-19", "MONDO_0100096"),
}
# Patch 3b 映射: PRJNA824137 real_group 前缀 → (cc, ih, Phenotype, Phenotype_ID)
P3B_MAPPING = {
    "HC":  ("control", True,  "Health",       "EFO_0010130"),
    "TBZ": ("case",    False, "Tuberculosis", "MONDO_0018076"),
}

EXPECTED = {"P1": 1257, "P2": 255, "P3A": 113, "P3B": 16, "total": 1641}

# %% [markdown]
# ## 1. 读 metadata + patch tsv

# %%
md = pd.read_parquet(IN_META)
print(f"metadata: {md.shape}  (Run 唯一: {md['Run'].is_unique})")

patches = pd.read_csv(IN_PATCH, sep="\t", dtype=str, keep_default_na=False, na_values=[""])
# 清理 host_disease / real_group 的前导空格（R 端 sample_data 写出时带的）
for c in ["host_disease", "real_group"]:
    patches[c] = patches[c].str.strip()
print(f"patches: {patches.shape}")
assert len(patches) == EXPECTED["total"], f"patch 总行数预期 {EXPECTED['total']}, 实际 {len(patches)}"

# 所有 patch Run 必须能命中 metadata
md_idx = md.set_index("Run", drop=False)
miss = set(patches["Run"]) - set(md_idx.index)
assert not miss, f"{len(miss)} 个 patch Run 在 metadata 里找不到"

# %% [markdown]
# ## 2. 记录原值快照（用于改后对比）

# %%
TOUCHED_COLS = ["Sample_Site", "Phenotype", "Phenotype_ID", "Case_Or_Control", "Is_Healthy"]
runs_all = patches["Run"].values
before = md_idx.loc[runs_all, TOUCHED_COLS].copy()
print(f"待改的 {len(runs_all)} 个 Run 在 5 列上的原值已记录")

# %% [markdown]
# ## 3. Patch 1: Sample_Site `Nasal` → `Nasopharynx` （PRJNA914884 + PRJNA1058141）

# %%
P1 = patches[
    patches["Project_ID"].isin(["PRJNA914884", "PRJNA1058141"]) &
    (patches["Body_Site_old"] == "Nasal") &
    (patches["Body_Site_new"] == "Nasopharynx")
]
assert len(P1) == EXPECTED["P1"], f"Patch 1 预期 {EXPECTED['P1']} 行, 实际 {len(P1)}"

runs = P1["Run"].values
md_idx.loc[runs, "Sample_Site"] = "Nasopharynx"
print(f"Patch 1: Sample_Site Nasal → Nasopharynx, {len(P1)} 行")
print(P1.groupby("Project_ID").size().to_string())

# %% [markdown]
# ## 4. Patch 2: Phenotype + Phenotype_ID （PRJNA801796 Influenza 细分）

# %%
P2 = patches[patches["Project_ID"] == "PRJNA801796"]
assert len(P2) == EXPECTED["P2"], f"Patch 2 预期 {EXPECTED['P2']} 行, 实际 {len(P2)}"
# 检查 Phenotype_new 取值集合与 mapping 一致
got = set(P2["Phenotype_new"].unique())
exp = set(P2_MAPPING.keys())
assert got == exp, f"Patch 2 出现意外 Phenotype_new: {got - exp}; 缺失: {exp - got}"

for jxt_ph, (tgt_ph, tgt_id) in P2_MAPPING.items():
    runs = P2.loc[P2["Phenotype_new"] == jxt_ph, "Run"].values
    md_idx.loc[runs, "Phenotype"]    = tgt_ph
    md_idx.loc[runs, "Phenotype_ID"] = tgt_id
    print(f"Patch 2: {jxt_ph!r} → ({tgt_ph!r}, {tgt_id}), {len(runs)} 行")

# %% [markdown]
# ## 5. Patch 3a: cc + ih + Phenotype + Phenotype_ID （PRJNA822681 按 host.disease）

# %%
P3A = patches[patches["Project_ID"] == "PRJNA822681"]
assert len(P3A) == EXPECTED["P3A"], f"Patch 3a 预期 {EXPECTED['P3A']} 行, 实际 {len(P3A)}"
got = set(P3A["host_disease"].unique())
exp = set(P3A_MAPPING.keys())
assert got == exp, f"Patch 3a 出现意外 host_disease: {got - exp}; 缺失: {exp - got}"

for hd, (cc, ih, ph, ph_id) in P3A_MAPPING.items():
    runs = P3A.loc[P3A["host_disease"] == hd, "Run"].values
    md_idx.loc[runs, "Case_Or_Control"] = cc
    md_idx.loc[runs, "Is_Healthy"]      = ih
    md_idx.loc[runs, "Phenotype"]       = ph
    md_idx.loc[runs, "Phenotype_ID"]    = ph_id
    print(f"Patch 3a: host.disease={hd!r} → cc={cc}, ih={ih}, ph={ph!r}, ph_id={ph_id}, {len(runs)} 行")

# %% [markdown]
# ## 6. Patch 3b: cc + ih + Phenotype + Phenotype_ID （PRJNA824137 按 real_group HC*/TBZ*）

# %%
P3B = patches[patches["Project_ID"] == "PRJNA824137"]
assert len(P3B) == EXPECTED["P3B"], f"Patch 3b 预期 {EXPECTED['P3B']} 行, 实际 {len(P3B)}"

# 按 real_group 前缀分组
def rg_prefix(s: str) -> str:
    s = (s or "").strip()
    for pfx in P3B_MAPPING:
        if s.startswith(pfx):
            return pfx
    return ""

prefixes = P3B["real_group"].apply(rg_prefix)
got = set(prefixes.unique())
exp = set(P3B_MAPPING.keys())
unrecognized = got - exp
assert not unrecognized, f"Patch 3b 出现无法识别的 real_group 前缀: {unrecognized}"

# 验证比例: HC 13 + TBZ 3 = 16
n_hc  = (prefixes == "HC").sum()
n_tbz = (prefixes == "TBZ").sum()
assert n_hc == 13 and n_tbz == 3, f"PRJNA824137 HC/TBZ 分布异常: HC={n_hc}, TBZ={n_tbz}"

for pfx, (cc, ih, ph, ph_id) in P3B_MAPPING.items():
    runs = P3B.loc[prefixes == pfx, "Run"].values
    md_idx.loc[runs, "Case_Or_Control"] = cc
    md_idx.loc[runs, "Is_Healthy"]      = ih
    md_idx.loc[runs, "Phenotype"]       = ph
    md_idx.loc[runs, "Phenotype_ID"]    = ph_id
    print(f"Patch 3b: real_group={pfx}* → cc={cc}, ih={ih}, ph={ph!r}, ph_id={ph_id}, {len(runs)} 行")

# %% [markdown]
# ## 7. 改后 sanity checks

# %%
md_fixed = md_idx.reset_index(drop=True)

# 7.1 shape 不变
assert md_fixed.shape == md.shape, f"shape 变了: {md.shape} → {md_fixed.shape}"
# 7.2 dtype 不变（categorical/boolean 保留）
for c in TOUCHED_COLS:
    assert str(md[c].dtype) == str(md_fixed[c].dtype), f"{c} dtype 变了: {md[c].dtype} → {md_fixed[c].dtype}"
# 7.3 未触及行（即 patch 之外的 Run）在 5 列上完全不变
touched_runs = set(patches["Run"])
untouched_mask = ~md_fixed["Run"].isin(touched_runs)
for c in TOUCHED_COLS:
    a = md.loc[untouched_mask, c]
    b = md_fixed.loc[untouched_mask, c]
    # 处理 NA: 都 NA 等价
    eq = (a.astype(object) == b.astype(object)) | (a.isna() & b.isna())
    assert eq.all(), f"{c} 在未触及行上有 {(~eq).sum()} 处变化"
print(f"未触及 {untouched_mask.sum():,} 行在 5 列上完全不变 ✓")

# 7.4 patch 行上每个项目的目标值
md_idx_fixed = md_fixed.set_index("Run", drop=False)

# Patch 1
v = md_idx_fixed.loc[P1["Run"].values, "Sample_Site"]
assert (v == "Nasopharynx").all(), "Patch 1 应用后存在非 Nasopharynx"
# Patch 2
for jxt_ph, (tgt_ph, tgt_id) in P2_MAPPING.items():
    rs = P2.loc[P2["Phenotype_new"] == jxt_ph, "Run"].values
    assert (md_idx_fixed.loc[rs, "Phenotype"] == tgt_ph).all()
    assert (md_idx_fixed.loc[rs, "Phenotype_ID"] == tgt_id).all()
# Patch 3a
for hd, (cc, ih, ph, ph_id) in P3A_MAPPING.items():
    rs = P3A.loc[P3A["host_disease"] == hd, "Run"].values
    assert (md_idx_fixed.loc[rs, "Case_Or_Control"] == cc).all()
    assert (md_idx_fixed.loc[rs, "Is_Healthy"] == ih).all()
    assert (md_idx_fixed.loc[rs, "Phenotype"] == ph).all()
    assert (md_idx_fixed.loc[rs, "Phenotype_ID"] == ph_id).all()
# Patch 3b
for pfx, (cc, ih, ph, ph_id) in P3B_MAPPING.items():
    rs = P3B.loc[prefixes == pfx, "Run"].values
    assert (md_idx_fixed.loc[rs, "Case_Or_Control"] == cc).all()
    assert (md_idx_fixed.loc[rs, "Is_Healthy"] == ih).all()
    assert (md_idx_fixed.loc[rs, "Phenotype"] == ph).all()
    assert (md_idx_fixed.loc[rs, "Phenotype_ID"] == ph_id).all()
print("所有 patch 行的目标值都已正确写入 ✓")

# 7.5 Is_Healthy 三态语义保持: True ⇔ Phenotype == 'Health'
n_health = (md_fixed["Phenotype"] == "Health").sum()
n_ih_true = (md_fixed["Is_Healthy"] == True).sum()
print(f"Phenotype=='Health': {n_health:,}; Is_Healthy==True: {n_ih_true:,}")
# 我们的 patch 也保持这个不变式：所有 Health phenotype 都对应 Is_Healthy=True
# Patch 1 / Patch 2 不动 Phenotype 与 Is_Healthy 的关系
# Patch 3a/3b 都按"Health ↔ True / 其他 ↔ False" 同步改
assert n_ih_true == n_health, f"Is_Healthy 与 Phenotype=='Health' 不一致: {n_ih_true} vs {n_health}"

# %% [markdown]
# ## 8. 改动 audit 摘要 + 写出

# %%
print("\n=== 改动 audit (按列统计实际变化行数) ===")
changes = []
for c in TOUCHED_COLS:
    a = before[c].reset_index(drop=True).astype(object)
    b = md_idx_fixed.loc[patches["Run"].values, c].reset_index(drop=True).astype(object)
    diff = (a != b) | (a.isna() ^ b.isna())
    n_chg = diff.sum()
    changes.append((c, n_chg))
    print(f"  {c:<20} {n_chg:>6} / {len(patches)} 行被改")

# 总计
total_unique_runs_changed = len(patches)
print(f"\n总命中 Run: {total_unique_runs_changed:,}（不重复，5 项目无交集）")

md_fixed.to_parquet(OUT_META)
sz = OUT_META.stat().st_size / 1024 / 1024
print(f"\n写出: {OUT_META}  ({sz:.1f} MB)")

# %% [markdown]
# ## 9. 回读 sanity

# %%
back = pd.read_parquet(OUT_META)
assert back.shape == md.shape
for c in TOUCHED_COLS:
    assert str(back[c].dtype) == str(md[c].dtype)
print(f"回读 OK: {back.shape}, dtype 保留 ✓")
