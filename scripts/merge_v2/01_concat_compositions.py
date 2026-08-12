#!/usr/bin/env python3
"""Layer2 步1：五库组成表 → 对齐 8,114 GG2 词表 → vstack 大 counts + 旗标。

产物 `results/feature_table/merged.counts.h5ad`：
- X = int32 counts，**纯 8,114 深层属**（varp/varm 待 04 挂）
- obs（全样本，一行不删）：
    sample_id = `<DB>:<native>`（全局唯一）· Database · Run（MA 取点前缀；Qiita 经 ENA 桥）
    unresolved_reads（浅层桶求和）· total_reads · n_deep_reads · n_deep_features
    qc_pass（深层 reads≥1000 且 深层属≥5，**仅旗标不删**）· Run_paired_id（同 Run 跨库副本串联）· n_run_copies
- var = 8,114（Domain..Genus，取自 MCFCorpusV2）

设计依据：QC=打标不删、unresolved 进 obs（见 .claude/merge_architecture.md）。
五库深层属已实测 100% ⊆ 8,114 词表（零损失对齐）。
"""
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import h5py, csv, collections

ROOT = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
FT   = ROOT / "results/feature_table"
VOCAB_SRC = FT / "MCFCorpusV2.gg2.h5ad"          # 8,114 var（含 Domain..Genus schema）
OUT  = FT / "merged.counts.h5ad"
BRIDGE = ROOT / "rawdata/Qiita/qiita_ena_bridge.tsv"
MIN_READS, MIN_FEAT = 1000, 5

DBS = [
    ("MA",    FT / "gg2.full.h5ad"),
    ("RM",    FT / "resmicrodb.gg2.genus.h5ad"),
    ("MBM",   FT / "microbiomap.gg2.genus.h5ad"),
    ("PRIME", FT / "prime.gg2.genus.h5ad"),
    ("Qiita", FT / "qiita.gg2.genus.h5ad"),
]

def up(s): return s.strip().upper()

# ---- 0. 读 8,114 词表 var（保序 + schema）----
print("读 8,114 词表 ...", flush=True)
vocab_ad = ad.read_h5ad(VOCAB_SRC, backed="r")
VOCAB = list(vocab_ad.var_names)
assert len(VOCAB) == 8114, len(VOCAB)
vocab_idx = {p: i for i, p in enumerate(VOCAB)}
VAR_DF = vocab_ad.var.copy()
vocab_ad.file.close()
print(f"  词表 {len(VOCAB)} 属")

# ---- Qiita ENA 桥：study.sample -> 代表 run ----
print("读 Qiita ENA 桥 ...", flush=True)
ss2run = {}
with open(BRIDGE) as f:
    rd = csv.DictReader(f, delimiter="\t")
    for row in rd:
        al = row["sample_alias"]; ss = al.split(":", 1)[1] if ":" in al else al
        r = row["run_accession"].strip()
        if r and ss not in ss2run:        # 取第一个 run 作代表
            ss2run[ss] = up(r)
print(f"  桥覆盖 {len(ss2run):,} 个 Qiita 样本")

def run_of(db, native):
    if db == "MA":    return up(native.split(".", 1)[0])
    if db == "Qiita": return ss2run.get(native)          # 可能 None
    return up(native)                                     # RM/MBM/PRIME

# ---- 1. 逐库对齐 + 旗标 ----
X_blocks = []
obs_rows = []   # dict list
for db, path in DBS:
    print(f"\n=== {db} : {path.name} ===", flush=True)
    A = ad.read_h5ad(path)
    var_ids = list(A.var_names)
    n_obs, n_var = A.shape
    # 投影矩阵 P (n_var × 8114)：深层 j → vocab 列；浅层 j 不进 P
    rows_j, cols_v = [], []
    for j, p in enumerate(var_ids):
        vi = vocab_idx.get(p)
        if vi is not None:
            rows_j.append(j); cols_v.append(vi)
    P = sp.csr_matrix((np.ones(len(rows_j), dtype=np.int32), (rows_j, cols_v)),
                      shape=(n_var, 8114))
    X = A.X.tocsr().astype(np.int32)
    aligned = (X @ P).tocsr()                            # n_obs × 8114（仅深层）
    total      = np.asarray(X.sum(axis=1)).ravel().astype(np.int64)
    deep_reads = np.asarray(aligned.sum(axis=1)).ravel().astype(np.int64)
    unresolved = total - deep_reads
    deep_feats = aligned.getnnz(axis=1).astype(np.int32)
    qc_pass = (deep_reads >= MIN_READS) & (deep_feats >= MIN_FEAT)
    assert (unresolved >= 0).all(), "unresolved 负值"
    print(f"  对齐 {aligned.shape}  深层属投影列={len(rows_j)}  qc_pass={int(qc_pass.sum()):,}/{n_obs:,}")

    natives = list(A.obs_names)
    for i, nat in enumerate(natives):
        obs_rows.append((
            f"{db}:{nat}", db, run_of(db, nat),
            int(unresolved[i]), int(total[i]), int(deep_reads[i]), int(deep_feats[i]),
            bool(qc_pass[i]),
        ))
    X_blocks.append(aligned)
    del A, X, aligned, P

# ---- 2. vstack + obs 组装 ----
print("\nvstack ...", flush=True)
Xall = sp.vstack(X_blocks, format="csr"); del X_blocks
print(f"  合并 X = {Xall.shape}  nnz={Xall.nnz:,}")
obs = pd.DataFrame(obs_rows, columns=[
    "sample_id", "Database", "Run", "unresolved_reads", "total_reads",
    "n_deep_reads", "n_deep_features", "qc_pass"]).set_index("sample_id")
assert obs.index.is_unique, "全局 sample_id 不唯一"
assert Xall.shape[0] == len(obs)

# ---- 3. Run_paired_id（同 Run 跨库副本）----
print("算 Run_paired_id ...", flush=True)
run_counts = obs["Run"].value_counts(dropna=True)
dup_runs = set(run_counts[run_counts >= 2].index)
obs["Run_paired_id"] = np.where(obs["Run"].isin(dup_runs), obs["Run"], pd.NA)
obs["n_run_copies"] = obs["Run"].map(run_counts).fillna(1).astype(int)
n_pair_rows = int(obs["Run_paired_id"].notna().sum())
print(f"  有跨库副本的行 = {n_pair_rows:,}（涉及 {len(dup_runs):,} 个 Run）")

# ---- 4. 组装 AnnData + 写出 ----
A = ad.AnnData(X=Xall, obs=obs, var=VAR_DF)
A.var_names = VOCAB
A.uns["build"] = "merge_v2/01_concat_compositions.py"
A.uns["qc_def"] = f"qc_pass = n_deep_reads>={MIN_READS} & n_deep_features>={MIN_FEAT}（仅旗标，不删行）"

print("\n=== 汇总 ===")
print(A)
print(obs["Database"].value_counts().to_string())
print(f"qc_pass 通过: {int(obs['qc_pass'].sum()):,} / {A.n_obs:,} ({100*obs['qc_pass'].mean():.1f}%)")
print(f"写出 {OUT} ...", flush=True)
A.write_h5ad(OUT, compression="gzip")
print(f"✓ 完成 {OUT}  ({A.n_obs:,} × {A.n_vars:,})")
