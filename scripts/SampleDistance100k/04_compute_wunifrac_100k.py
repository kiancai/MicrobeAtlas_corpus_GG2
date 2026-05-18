from pathlib import Path
import os
import shutil
import time
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
import biom
import h5py
import unifrac

ad.settings.allow_write_nullable_strings = True

ROOT = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
OUT_DIR = ROOT / "results/sample_distance_100k"
ANN_PATH = OUT_DIR / "subset_100k.h5ad"
TREE_PATH = OUT_DIR / "genus_tree.nwk"
TREE_FALLBACK = ROOT / "results/sample_distance/genus_tree.nwk"
VOCAB_IN = ROOT / "results/phylogeny/genus_vocab.tsv"
OUT = OUT_DIR / "distance_wunifrac.float16.npy"
TMP_DIR = OUT_DIR / "_tmp_unifrac"
BIOM_TMP = TMP_DIR / "subset_100k.biom"

N_THREADS = int(os.environ.get("OMP_NUM_THREADS", os.cpu_count() or 1))
os.environ["OMP_NUM_THREADS"] = str(N_THREADS)

if not TREE_PATH.exists():
    assert TREE_FALLBACK.exists(), f"Missing tree: {TREE_PATH} and fallback {TREE_FALLBACK}"
    shutil.copy2(TREE_FALLBACK, TREE_PATH)
    print(f"Copied folded genus tree from {TREE_FALLBACK} to {TREE_PATH}")

TMP_DIR.mkdir(parents=True, exist_ok=True)

print(f"OMP_NUM_THREADS={N_THREADS}")
print(f"Reading {ANN_PATH.name} ...")
adata = ad.read_h5ad(ANN_PATH)
assert adata.n_obs == 100_000, f"Expected 100000 samples, got {adata.n_obs}"
print(f"shape={adata.shape}")

vocab = pd.read_csv(VOCAB_IN, sep="\t", index_col="var_id")
var_ids = adata.var_names.astype(str).to_numpy()
missing = set(var_ids) - set(vocab.index.astype(str))
assert len(missing) == 0, f"{len(missing)} vars missing from genus vocab"
genus_for_var = vocab.loc[var_ids, "Genus"].astype(str).to_numpy()
assert len(set(genus_for_var)) == len(genus_for_var), "var to genus mapping is not one-to-one"

print("Building relative-abundance BIOM table ...")
t0 = time.time()
X = adata.X.astype(np.float64)
row_sum = np.asarray(X.sum(axis=1)).ravel()
assert (row_sum > 0).all(), "Zero-count rows found"
P = sp.diags(1.0 / row_sum) @ X
data_T = P.T.tocsr()
sample_ids = adata.obs_names.astype(str).tolist()
table = biom.Table(
    data=data_T,
    observation_ids=genus_for_var.tolist(),
    sample_ids=sample_ids,
)
with h5py.File(BIOM_TMP, "w") as h:
    table.to_hdf5(h, generated_by="SampleDistance100k/04_compute_wunifrac_100k")
print(f"BIOM size={BIOM_TMP.stat().st_size / 1024**2:.1f} MB elapsed={time.time() - t0:.1f}s")
del X, P, data_T, table

print("Computing weighted normalized UniFrac ...")
t0 = time.time()
dm = unifrac.weighted_normalized(
    table=str(BIOM_TMP),
    phylogeny=str(TREE_PATH),
    threads=N_THREADS,
    variance_adjusted=False,
    bypass_tips=False,
)
print(f"UniFrac elapsed={time.time() - t0:.1f}s shape={dm.shape}")

D_source = dm.data
ids = list(dm.ids)
pos = {sid: i for i, sid in enumerate(ids)}
order = np.array([pos[sid] for sid in sample_ids], dtype=np.int64)
assert len(order) == adata.n_obs

N = adata.n_obs
D16 = np.lib.format.open_memmap(OUT, mode="w+", dtype=np.float16, shape=(N, N))
BATCH = int(os.environ.get("WRITE_BATCH", "500"))
print(f"Writing aligned float16 npy: {OUT} WRITE_BATCH={BATCH}")
for start in range(0, N, BATCH):
    end = min(start + BATCH, N)
    rows = order[start:end]
    if np.array_equal(order, np.arange(N)):
        chunk = D_source[start:end, :]
    else:
        chunk = D_source[rows, :][:, order]
    chunk = np.asarray(chunk, dtype=np.float32)
    if start < end:
        diag = np.arange(start, end)
        chunk[np.arange(end - start), diag] = 0.0
    D16[start:end, :] = chunk.astype(np.float16)
    if (start // BATCH + 1) % 20 == 0 or end == N:
        print(f"  wrote rows {end:,}/{N:,}")
D16.flush()
print(f"Wrote {OUT} size={OUT.stat().st_size / 1024**3:.2f} GB")

print("Sanity audit on a 2k random submatrix ...")
rng = np.random.default_rng(20260516)
audit_n = min(2000, N)
idx = np.sort(rng.choice(N, audit_n, replace=False))
S = np.asarray(D16[np.ix_(idx, idx)], dtype=np.float32)
np.fill_diagonal(S, 0.0)
print(f"diag max={np.diag(S).max():.6e}")
print(f"sym max diff={np.abs(S - S.T).max():.6e}")
print(f"min={S.min():.6f} max={S.max():.6f} median_upper={np.median(S[np.triu_indices_from(S, k=1)]):.4f}")
print(f"nan/inf={np.isnan(S).sum()} / {np.isinf(S).sum()}")
assert np.isnan(S).sum() == 0 and np.isinf(S).sum() == 0
assert S.min() >= -1e-3 and S.max() <= 1.001

try:
    BIOM_TMP.unlink()
    TMP_DIR.rmdir()
except OSError as exc:
    print(f"Temporary cleanup skipped: {exc}")
