from pathlib import Path
import os
import time
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np
import anndata as ad
import scipy.sparse as sp
from scipy.spatial.distance import cdist

ROOT = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
ANN_PATH = ROOT / "results/sample_distance_100k/subset_100k.h5ad"
OUT = ROOT / "results/sample_distance_100k/distance_bc.float16.npy"

N_JOBS = int(os.environ.get("N_JOBS", len(os.sched_getaffinity(0))))
BATCH = int(os.environ.get("BATCH", "250"))

print(f"N_JOBS={N_JOBS} BATCH={BATCH}")
print(f"Reading {ANN_PATH.name} ...")
adata = ad.read_h5ad(ANN_PATH)
assert adata.n_obs == 100_000, f"Expected 100000 samples, got {adata.n_obs}"
assert sp.issparse(adata.X), "Expected sparse count matrix"
print(f"shape={adata.shape} X dtype={adata.X.dtype} nnz={adata.X.nnz:,}")

print("Converting X to dense row-normalized relative abundance float32 ...")
t0 = time.time()
X = adata.X.astype(np.float32)
row_sum = np.asarray(X.sum(axis=1)).ravel()
assert (row_sum > 0).all(), "Zero-count rows found"
P = (sp.diags(1.0 / row_sum.astype(np.float64)) @ X).toarray().astype(np.float32)
del X, adata
assert np.allclose(P.sum(axis=1), 1.0, atol=1e-4)
print(f"P: shape={P.shape} mem={P.nbytes / 1024**3:.2f} GB elapsed={time.time() - t0:.1f}s")

# scipy cdist braycurtis internally converts inputs to float64 via C.
# With 64 threads each receiving a float32 P, each thread would create its own
# 6.5 GB float64 copy of P → 416 GB peak → OOM.
# Pre-convert once so threads get a direct float64 pointer with no per-thread copy.
print("Converting P float32 → float64 for cdist (one-time, shared across threads) ...")
P64 = P.astype(np.float64)
del P
print(f"P64: mem={P64.nbytes / 1024**3:.2f} GB")

N = P64.shape[0]
D = np.lib.format.open_memmap(OUT, mode="w+", dtype=np.float16, shape=(N, N))

chunk_args = [(s, min(s + BATCH, N)) for s in range(0, N, BATCH)]
total = len(chunk_args)
print(f"Computing Bray-Curtis: chunks={total} workers={N_JOBS}")

done = 0
lock = threading.Lock()
t0 = time.time()


def compute_chunk(args):
    start, end = args
    chunk = cdist(P64[start:end], P64, metric="braycurtis").astype(np.float32)
    # zero diagonal
    diag_start, diag_end = max(start, 0), min(end, N)
    if diag_start < diag_end:
        local = np.arange(diag_start, diag_end) - start
        chunk[local, np.arange(diag_start, diag_end)] = 0.0
    D[start:end, :] = chunk.astype(np.float16)
    return start


with ThreadPoolExecutor(max_workers=N_JOBS) as pool:
    futures = {pool.submit(compute_chunk, a): a for a in chunk_args}
    for fut in as_completed(futures):
        fut.result()  # re-raises any exception immediately
        with lock:
            done += 1
            if done % 20 == 0 or done == total:
                elapsed = time.time() - t0
                eta = elapsed * (total - done) / max(done, 1)
                print(f"[{done:>4}/{total}] elapsed={elapsed / 60:.1f} min ETA={eta / 60:.1f} min",
                      flush=True)

D.flush()
print(f"Wrote {OUT} size={OUT.stat().st_size / 1024**3:.2f} GB")

print("Sanity audit on a 2k random submatrix ...")
rng = np.random.default_rng(20260516)
audit_n = min(2000, N)
idx = np.sort(rng.choice(N, audit_n, replace=False))
S = np.asarray(D[np.ix_(idx, idx)], dtype=np.float32)
np.fill_diagonal(S, 0.0)
print(f"diag max={np.diag(S).max():.6e}")
print(f"sym max diff={np.abs(S - S.T).max():.6e}")
print(f"min={S.min():.6f} max={S.max():.6f} median_upper={np.median(S[np.triu_indices_from(S, k=1)]):.4f}")
print(f"nan/inf={np.isnan(S).sum()} / {np.isinf(S).sum()}")
assert np.isnan(S).sum() == 0 and np.isinf(S).sum() == 0
assert S.min() >= -1e-3 and S.max() <= 1.001
