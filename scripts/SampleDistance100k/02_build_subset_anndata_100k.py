from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

ad.settings.allow_write_nullable_strings = True

ROOT = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
ANN_IN = ROOT / "results/feature_table/merged.gg2.with_phylo.h5ad"
IDX_IN = ROOT / "results/sample_distance_100k/subset_100k_index.tsv"
ANN_OUT = ROOT / "results/sample_distance_100k/subset_100k.h5ad"

INDEX_META_COLS = [
    "sample_role",
    "overview_bucket",
    "ma_bucket_detail",
    "human_site",
    "rm_site",
    "paired_run_id",
    "stratum_id",
    "sub_stratum",
]

print(f"Reading index: {IDX_IN}")
idx = pd.read_csv(IDX_IN, sep="\t", dtype=str)
assert len(idx) == 100_000, f"Expected 100000 index rows, got {len(idx)}"
assert idx["obs_name"].is_unique, "index obs_name is not unique"

print(f"Reading source anndata backed: {ANN_IN.name}")
adata = ad.read_h5ad(ANN_IN, backed="r")
all_names = adata.obs_names.astype(str).to_numpy()
name_to_pos = pd.Series(np.arange(len(all_names)), index=all_names)
missing = idx.loc[~idx["obs_name"].isin(name_to_pos.index), "obs_name"]
assert len(missing) == 0, f"{len(missing)} obs_names missing from source, first={missing.head().tolist()}"
row_pos = name_to_pos.loc[idx["obs_name"].to_numpy()].to_numpy(dtype=np.int64)
assert len(np.unique(row_pos)) == len(row_pos), "row positions are not unique"

print("Slicing X ...")
X_sub = adata.X[row_pos, :]
if not sp.isspmatrix_csr(X_sub):
    X_sub = sp.csr_matrix(X_sub)
print(f"X_sub: shape={X_sub.shape} nnz={X_sub.nnz:,} dtype={X_sub.dtype}")

print("Slicing obs and adding 100k sample metadata ...")
obs_sub = adata.obs.iloc[row_pos].copy()
obs_sub = obs_sub.reset_index(drop=True)
obs_sub.index = idx["obs_name"].to_numpy()
obs_sub.index.name = None

for col in INDEX_META_COLS:
    values = idx[col].fillna("NA").astype(str).to_numpy()
    if col in {"sample_role", "overview_bucket", "ma_bucket_detail", "stratum_id"}:
        obs_sub[col] = pd.Categorical(values)
    else:
        obs_sub[col] = pd.array(values, dtype="string")

print("Copying var and varp ...")
var_sub = adata.var.copy()
varp_dict = {}
for key in adata.varp.keys():
    arr = np.asarray(adata.varp[key])
    varp_dict[key] = arr
    print(f"varp['{key}']: shape={arr.shape} dtype={arr.dtype} mem={arr.nbytes / 1024**2:.1f} MB")

sub = ad.AnnData(X=X_sub, obs=obs_sub, var=var_sub)
for key, arr in varp_dict.items():
    sub.varp[key] = arr

print(f"Writing {ANN_OUT}")
sub.write_h5ad(ANN_OUT, compression="gzip")
print(f"size={ANN_OUT.stat().st_size / 1024**3:.2f} GB")

print("Read-back sanity ...")
b = ad.read_h5ad(ANN_OUT, backed="r")
print(f"shape={b.shape}")
print(f"obs columns={len(b.obs.columns)}")
print(f"varp={list(b.varp.keys())}")
print(b.obs["sample_role"].value_counts().to_string())
