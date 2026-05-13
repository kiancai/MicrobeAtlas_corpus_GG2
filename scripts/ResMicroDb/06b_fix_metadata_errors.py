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
# # ResMicroDb 06b: 应用 NCBI 真值 + 0105 phyloseq 提供的修正到 metadata_all
#
# 把 5 个项目的 metadata 错标修正应用到 `metadata_all.standardized.parquet`，
# 输出 `.fixed.parquet`。后续 07_attach_metadata 默认消费 fixed 版。
#
# **四个 patch**（详见 `rawdata/ResMicroDb/supplement data/CHANGES_0105_new7.md §6`）：
#
# 1. **Sample_Site**：PRJNA914884 (1,195) + PRJNA1058141 (62) 的 `Nasal` → `Nasopharynx`
#    — 来源：0105 phyloseq diff；**1,257** 行
# 2. **Phenotype + Phenotype_ID**：PRJNA801796 (Influenza 细分)，按 schema 标准词 + 真实
#    ontology ID — 来源：0105 phyloseq sample_data 的 Run 级映射；**255** 行
# 3. **Phenotype + Phenotype_ID + cc + ih**：PRJNA822681 错标 — 来源：NCBI BioSample
#    `host disease` + SRA RunInfo Run↔BioSample 映射；**152** 行
# 4. **Phenotype + Phenotype_ID + cc + ih**：PRJNA824137 错标 — 来源：NCBI BioSample
#    sample title 前缀（HC/TBZ/TBM/LTBI）；**33** 行
#
# 总：**1,697** 个 Run 被修改（项目不重叠）。其他 134,000+ 行不变。
#
# **依赖**：`06b_export_patches.R` 先跑出 `metadata_patches_0105.tsv`（提供 PRJNA801796/
# 914884/1058141 的 Run 级 patch；PRJNA822681/824137 部分会被忽略，改由 NCBI biosample 接管）。

# %%
import re
from pathlib import Path

import pandas as pd

ROOT     = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
IN_META  = ROOT / "results/feature_table/metadata_all.standardized.parquet"
IN_PATCH = ROOT / "results/feature_table/metadata_patches_0105.tsv"
SUPP     = ROOT / "rawdata/ResMicroDb/supplement data"
NCBI_BS_822681  = SUPP / "PRJNA822681_biosample.txt"
NCBI_RI_822681  = SUPP / "PRJNA822681_runinfo.csv"
NCBI_BS_824137  = SUPP / "PRJNA824137_biosample.txt"
OUT_META = ROOT / "results/feature_table/metadata_all.standardized.fixed.parquet"

# Patch 2 映射: 0105 NEW Phenotype 字符串 → (目标 Phenotype, 目标 Phenotype_ID)
P2_MAPPING = {
    "Influenza A":                            ("Influenza A Virus",                     "NCIT_C53454"),
    "Influenza B":                            ("Influenza B Virus",                     "NCIT_C53468"),
    "Rhinovirus Infection":                   ("Rhinovirus Infection",                  "NCIT_C122572"),
    "Respiratory Syncytial Virus Infection":  ("Respiratory Syncytial Virus Infection", "EFO_1001413"),
}

# Patch 3 映射: PRJNA822681 NCBI host disease → (Phenotype, Phenotype_ID, cc, ih) 或 None
P3_MAPPING = {
    "Healthy":      ("Health",    "EFO_0010130",   "control", True),
    "COVID-19":     ("COVID-19",  "MONDO_0100096", "case",    False),
    "Non COVID-19": None,  # 不动，jxt 已用 Pneumonia 表达
}

# Patch 4 映射: PRJNA824137 NCBI title 前缀 → (Phenotype, Phenotype_ID, cc, ih)
P4_MAPPING = {
    "HC":   ("Health",                        "EFO_0010130",   "control", True),
    "TBZ":  ("Tuberculosis",                  "MONDO_0018076", "case",    False),
    "TBM":  ("Tuberculosis",                  "MONDO_0018076", "case",    False),
    "LTBI": ("Latent Tuberculosis Infection", "MONDO_0040753", "case",    False),
}

EXPECTED = {"P1": 1257, "P2": 255, "P3": 152, "P4": 33, "total": 1697}

# %% [markdown]
# ## 1. 读 metadata + 旧 patch tsv

# %%
md = pd.read_parquet(IN_META)
print(f"metadata: {md.shape}  (Run 唯一: {md['Run'].is_unique})")

patches = pd.read_csv(IN_PATCH, sep="\t", dtype=str, keep_default_na=False, na_values=[""])
for c in ["host_disease", "real_group"]:
    if c in patches.columns:
        patches[c] = patches[c].str.strip()
print(f"0105 patch tsv: {patches.shape}（PRJNA822681/824137 部分会被忽略，由 NCBI 真值接管）")

md_idx = md.set_index("Run", drop=False)

# %% [markdown]
# ## 2. 解析 NCBI biosample/runinfo → Patch 3 / Patch 4 的 Run 级表

# %%
def parse_biosample(fp: Path) -> pd.DataFrame:
    """解析 NCBI biosample 文本 dump（GenBank-style 多行块）"""
    records = []
    text = fp.read_text()
    for block in re.split(r"\n(?=\d+: )", text.strip()):
        rec = {}
        m_title = re.match(r"\d+:\s*(.+?)\n", block)
        if m_title:
            rec["title"] = m_title.group(1).strip()
        m_id = re.search(r"Identifiers:\s*(.+)", block)
        if m_id:
            for part in m_id.group(1).split(";"):
                k, _, v = part.strip().partition(":")
                rec[k.strip()] = v.strip()
        for m in re.finditer(r'/([^=]+)="([^"]*)"', block):
            rec[m.group(1).strip()] = m.group(2).strip()
        records.append(rec)
    return pd.DataFrame(records)

# Patch 3: PRJNA822681 (Run ← runinfo ← biosample.host_disease)
bs1 = parse_biosample(NCBI_BS_822681)[["BioSample", "host disease"]].rename(columns={"host disease": "NCBI_truth"})
ri1 = pd.read_csv(NCBI_RI_822681)[["Run", "BioSample"]]
ncbi1 = ri1.merge(bs1, on="BioSample", how="left")[["Run", "NCBI_truth"]]
assert len(ncbi1) == 221 and ncbi1["NCBI_truth"].notna().all(), "PRJNA822681 NCBI 真值覆盖不全"
print(f"NCBI PRJNA822681: {len(ncbi1)} 行；host disease 分布: {ncbi1['NCBI_truth'].value_counts().to_dict()}")

# Patch 4: PRJNA824137 (Run ← metadata_all 自身 ← biosample.title 前缀)
bs2 = parse_biosample(NCBI_BS_824137)[["BioSample", "title"]]
def title_prefix(t: str) -> str:
    t = (t or "").strip()
    for pfx in ("LTBI", "TBM", "TBZ", "HC"):  # LTBI/TBM 必须在 TBZ/HC 前匹配以免错切
        if t.startswith(pfx):
            return pfx
    return ""
bs2["NCBI_truth"] = bs2["title"].apply(title_prefix)
ri2 = md[md["Project_ID"] == "PRJNA824137"][["Run", "BioSample"]]
ncbi2 = ri2.merge(bs2[["BioSample", "NCBI_truth"]], on="BioSample", how="left")[["Run", "NCBI_truth"]]
assert len(ncbi2) == 67 and ncbi2["NCBI_truth"].notna().all() and (ncbi2["NCBI_truth"] != "").all(), "PRJNA824137 NCBI 真值覆盖不全"
print(f"NCBI PRJNA824137: {len(ncbi2)} 行；title 前缀分布: {ncbi2['NCBI_truth'].value_counts().to_dict()}")

# %% [markdown]
# ## 3. 记录原值快照（用于改后对比）

# %%
TOUCHED_COLS = ["Sample_Site", "Phenotype", "Phenotype_ID", "Case_Or_Control", "Is_Healthy"]
# 拼出所有被 patch 触及的 Run（去重；项目不重叠所以也可直接 concat）
touched_runs_all = pd.Index(
    list(patches.loc[patches["Project_ID"].isin(["PRJNA914884", "PRJNA1058141", "PRJNA801796"]), "Run"]) +
    list(ncbi1["Run"]) +
    list(ncbi2["Run"])
).unique()
before = md_idx.loc[touched_runs_all, TOUCHED_COLS].copy()
print(f"待 patch 涉及 {len(touched_runs_all)} 个 Run（5 列原值已快照）")

# %% [markdown]
# ## 4. Patch 1: Sample_Site `Nasal` → `Nasopharynx`（PRJNA914884 + PRJNA1058141）

# %%
P1 = patches[
    patches["Project_ID"].isin(["PRJNA914884", "PRJNA1058141"]) &
    (patches["Body_Site_old"] == "Nasal") &
    (patches["Body_Site_new"] == "Nasopharynx")
]
assert len(P1) == EXPECTED["P1"], f"Patch 1 预期 {EXPECTED['P1']} 行, 实际 {len(P1)}"

md_idx.loc[P1["Run"].values, "Sample_Site"] = "Nasopharynx"
print(f"Patch 1: Sample_Site Nasal → Nasopharynx, {len(P1)} 行")
print(P1.groupby("Project_ID").size().to_string())

# %% [markdown]
# ## 5. Patch 2: Phenotype + Phenotype_ID（PRJNA801796 Influenza 细分）

# %%
P2 = patches[patches["Project_ID"] == "PRJNA801796"]
assert len(P2) == EXPECTED["P2"], f"Patch 2 预期 {EXPECTED['P2']} 行, 实际 {len(P2)}"
got = set(P2["Phenotype_new"].unique())
exp = set(P2_MAPPING.keys())
assert got == exp, f"Patch 2 出现意外 Phenotype_new: {got - exp}; 缺失: {exp - got}"

for jxt_ph, (tgt_ph, tgt_id) in P2_MAPPING.items():
    runs = P2.loc[P2["Phenotype_new"] == jxt_ph, "Run"].values
    md_idx.loc[runs, "Phenotype"]    = tgt_ph
    md_idx.loc[runs, "Phenotype_ID"] = tgt_id
    print(f"Patch 2: {jxt_ph!r} → ({tgt_ph!r}, {tgt_id}), {len(runs)} 行")

# %% [markdown]
# ## 6. Patch 3: PRJNA822681 按 NCBI host disease（152 行）

# %%
# 152 行 = 76 NCBI Healthy + 76 NCBI COVID-19；NCBI Non COVID-19 (69) 不动
p3_to_change = ncbi1[ncbi1["NCBI_truth"].isin(["Healthy", "COVID-19"])]
assert len(p3_to_change) == EXPECTED["P3"], f"Patch 3 预期 {EXPECTED['P3']} 行, 实际 {len(p3_to_change)}"

for hd, target in P3_MAPPING.items():
    if target is None:
        continue
    ph, ph_id, cc, ih = target
    runs = p3_to_change.loc[p3_to_change["NCBI_truth"] == hd, "Run"].values
    md_idx.loc[runs, "Phenotype"]       = ph
    md_idx.loc[runs, "Phenotype_ID"]    = ph_id
    md_idx.loc[runs, "Case_Or_Control"] = cc
    md_idx.loc[runs, "Is_Healthy"]      = ih
    print(f"Patch 3: NCBI host disease={hd!r} → ph={ph!r}, ph_id={ph_id}, cc={cc}, ih={ih}, {len(runs)} 行")

# %% [markdown]
# ## 7. Patch 4: PRJNA824137 按 NCBI title 前缀（33 行）

# %%
# 33 行 = 13 HC + 3 TBZ (4/22/23) + 17 LTBI；14 TBM 已正确不动；20 TBZ 已正确不动
# 实现策略：对所有 67 行做"目标值赋值"，但因为 26 行已是目标值，set 等同操作（idempotent）
# 实际"被改"统计用 audit (§8) 来算

# 对全部 67 行应用 mapping
for pfx, (ph, ph_id, cc, ih) in P4_MAPPING.items():
    runs = ncbi2.loc[ncbi2["NCBI_truth"] == pfx, "Run"].values
    md_idx.loc[runs, "Phenotype"]       = ph
    md_idx.loc[runs, "Phenotype_ID"]    = ph_id
    md_idx.loc[runs, "Case_Or_Control"] = cc
    md_idx.loc[runs, "Is_Healthy"]      = ih
    print(f"Patch 4: NCBI title={pfx}* → ph={ph!r}, ph_id={ph_id}, cc={cc}, ih={ih}, {len(runs)} 行（含 idempotent）")

# %% [markdown]
# ## 8. 改后 sanity checks

# %%
md_fixed = md_idx.reset_index(drop=True)

# 8.1 shape 不变
assert md_fixed.shape == md.shape, f"shape 变了: {md.shape} → {md_fixed.shape}"

# 8.2 dtype 不变
for c in TOUCHED_COLS:
    assert str(md[c].dtype) == str(md_fixed[c].dtype), f"{c} dtype 变了: {md[c].dtype} → {md_fixed[c].dtype}"

# 8.3 未触及行（patch 之外的 Run）在 5 列上完全不变
untouched_mask = ~md_fixed["Run"].isin(touched_runs_all)
for c in TOUCHED_COLS:
    a = md.loc[untouched_mask, c]
    b = md_fixed.loc[untouched_mask, c]
    eq = (a.astype(object) == b.astype(object)) | (a.isna() & b.isna())
    assert eq.all(), f"{c} 在未触及行上有 {(~eq).sum()} 处变化"
print(f"未触及 {untouched_mask.sum():,} 行在 5 列上完全不变 ✓")

# 8.4 patch 行目标值验证
md_idx_fixed = md_fixed.set_index("Run", drop=False)

# Patch 1
assert (md_idx_fixed.loc[P1["Run"].values, "Sample_Site"] == "Nasopharynx").all()
# Patch 2
for jxt_ph, (tgt_ph, tgt_id) in P2_MAPPING.items():
    rs = P2.loc[P2["Phenotype_new"] == jxt_ph, "Run"].values
    assert (md_idx_fixed.loc[rs, "Phenotype"] == tgt_ph).all()
    assert (md_idx_fixed.loc[rs, "Phenotype_ID"] == tgt_id).all()
# Patch 3
for hd, target in P3_MAPPING.items():
    if target is None: continue
    ph, ph_id, cc, ih = target
    rs = ncbi1.loc[ncbi1["NCBI_truth"] == hd, "Run"].values
    assert (md_idx_fixed.loc[rs, "Phenotype"] == ph).all()
    assert (md_idx_fixed.loc[rs, "Phenotype_ID"] == ph_id).all()
    assert (md_idx_fixed.loc[rs, "Case_Or_Control"] == cc).all()
    assert (md_idx_fixed.loc[rs, "Is_Healthy"] == ih).all()
# Patch 4
for pfx, (ph, ph_id, cc, ih) in P4_MAPPING.items():
    rs = ncbi2.loc[ncbi2["NCBI_truth"] == pfx, "Run"].values
    assert (md_idx_fixed.loc[rs, "Phenotype"] == ph).all()
    assert (md_idx_fixed.loc[rs, "Phenotype_ID"] == ph_id).all()
    assert (md_idx_fixed.loc[rs, "Case_Or_Control"] == cc).all()
    assert (md_idx_fixed.loc[rs, "Is_Healthy"] == ih).all()
print("所有 patch 行的目标值都已正确写入 ✓")

# 8.5 不变式: Is_Healthy == True ⇔ Phenotype == 'Health'
n_health = (md_fixed["Phenotype"] == "Health").sum()
n_ih_true = (md_fixed["Is_Healthy"] == True).sum()
print(f"Phenotype=='Health': {n_health:,}; Is_Healthy==True: {n_ih_true:,}")
assert n_ih_true == n_health, f"Is_Healthy 与 Phenotype=='Health' 不一致: {n_ih_true} vs {n_health}"

# %% [markdown]
# ## 9. 改动 audit 摘要 + 写出

# %%
print("\n=== 改动 audit (按列统计实际变化行数) ===")
after_touched = md_idx_fixed.loc[touched_runs_all, TOUCHED_COLS]
for c in TOUCHED_COLS:
    a = before[c].astype(object).reset_index(drop=True)
    b = after_touched[c].astype(object).reset_index(drop=True)
    diff = (a != b) | (a.isna() ^ b.isna())
    print(f"  {c:<20} {diff.sum():>6} / {len(touched_runs_all)} 行被改")

# 按 patch 统计
print("\n=== 各 patch 实际变化行数 ===")
print(f"  Patch 1 (Sample_Site):         {len(P1)} 行（全改）")
print(f"  Patch 2 (Phenotype 细分):       {len(P2)} 行（全改 Phenotype + Phenotype_ID）")
n_p3 = len(p3_to_change)
print(f"  Patch 3 (PRJNA822681 NCBI):    {n_p3} 行（4 列）")
# Patch 4 细分: 真实"被改"行数（排除 idempotent）
p4_changed = 0
p4_ltbi_partial = 0
for pfx, (ph, ph_id, cc, ih) in P4_MAPPING.items():
    rs = ncbi2.loc[ncbi2["NCBI_truth"] == pfx, "Run"].values
    before_p4 = before.loc[rs] if all(r in before.index for r in rs) else md.set_index("Run").loc[rs, TOUCHED_COLS]
    row_changed = (
        (before_p4["Phenotype"].astype(object) != ph) |
        (before_p4["Phenotype_ID"].astype(object) != ph_id) |
        (before_p4["Case_Or_Control"].astype(object) != cc) |
        (before_p4["Is_Healthy"].astype(object) != ih)
    )
    print(f"    Patch 4-{pfx}: {row_changed.sum()}/{len(rs)} 行实际改变")
    p4_changed += row_changed.sum()
print(f"  Patch 4 合计实际改: {p4_changed} 行")
print(f"\n总实际改 Run: {len(P1) + len(P2) + n_p3 + p4_changed}（预期 {EXPECTED['total']}）")

md_fixed.to_parquet(OUT_META)
sz = OUT_META.stat().st_size / 1024 / 1024
print(f"\n写出: {OUT_META}  ({sz:.1f} MB)")

# %% [markdown]
# ## 10. 回读 sanity

# %%
back = pd.read_parquet(OUT_META)
assert back.shape == md.shape
for c in TOUCHED_COLS:
    assert str(back[c].dtype) == str(md[c].dtype)
print(f"回读 OK: {back.shape}, dtype 保留 ✓")
