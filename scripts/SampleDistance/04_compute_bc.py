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
# # 04: 50k × 50k Bray-Curtis 距离矩阵
#
# 输入：
# - `results/sample_distance/subset_50k.h5ad`  (50,000 × 8,114)
#
# 输出：
# - 把 50k × 50k BC 距离写回到同一个 anndata 的 `obsp['distance_bc']`（float16）
#
# ## 公式 & 实现
#
# Bray-Curtis:  `BC(u, v) = sum|u_i − v_i| / sum(u_i + v_i)`
#
# X 是 raw int32 CSR。先按行归一化为 relative abundance（每行除以 total reads），
# 然后做 dense matrix BC。每行归一后 sum(u+v) = 2，BC = sum|u−v| / 2 = 1 − sum_min(u,v)。
#
# 用 `scipy.spatial.distance.cdist(..., 'braycurtis')` 是 C-level 但单线程；50k full 约
# 3 小时（实测 1k × 50k 222s）。改用 `joblib.Parallel` 按行分块并行：每个 worker 处理
# 一片行，写到一个 numpy memmap 文件，最后合并。8 核理论 ~25 min。

# %%
from pathlib import Path
import os
import time
import multiprocessing as mp
import numpy as np
import anndata as ad
import scipy.sparse as sp
from scipy.spatial.distance import cdist

ad.settings.allow_write_nullable_strings = True

ROOT = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
ANN_PATH = ROOT / "results/sample_distance/subset_50k.h5ad"

# 并行 worker 数：cgroup affinity 限定的 CPU 数
N_JOBS = int(os.environ.get("N_JOBS", len(os.sched_getaffinity(0))))
BATCH = 500  # 每个任务处理 500 行 vs 全 50k；总块数 100，worker 抢占式领任务
print(f"N_JOBS = {N_JOBS}  BATCH = {BATCH}")

# %% [markdown]
# ## §1 读 anndata + 转 relative abundance dense float32

# %%
print(f"\n读 {ANN_PATH.name} ...")
adata = ad.read_h5ad(ANN_PATH)
print(f"  shape: {adata.shape}")
print(f"  X: dtype={adata.X.dtype}  nnz={adata.X.nnz:,}")

# %%
print("\n转 relative abundance dense float32 ...")
t0 = time.time()
X = adata.X.astype(np.float64)
row_sum = np.asarray(X.sum(axis=1)).ravel()
assert (row_sum > 0).all(), "存在零行，无法做 BC（需先剔除）"

inv = sp.diags(1.0 / row_sum)
P_sparse = inv @ X
P = P_sparse.toarray().astype(np.float32)
del X, P_sparse, inv
print(f"  P shape: {P.shape}  mem: {P.nbytes/1024**3:.2f} GB  耗时 {time.time()-t0:.1f}s")
assert np.allclose(P.sum(axis=1), 1.0, atol=1e-4)

# %% [markdown]
# ## §2 并行计算 BC：每个 worker 处理 BATCH 行 vs 全表

# %%
N = P.shape[0]

# 全局共享 P：通过 fork 让子进程共享父进程内存（copy-on-write，只读不会真复制）
_P_GLOBAL = P

def compute_chunk(args):
    """子进程入口：使用 fork 继承的全局 P。"""
    start, end = args
    P_local = _P_GLOBAL
    return start, cdist(P_local[start:end], P_local, metric="braycurtis").astype(np.float32)


chunk_args = [(s, min(s + BATCH, N)) for s in range(0, N, BATCH)]
print(f"\n并行计算 BC: N={N}, N_JOBS={N_JOBS}, BATCH={BATCH}, 总块数={len(chunk_args)}")

D = np.zeros((N, N), dtype=np.float32)
t0 = time.time()
ctx = mp.get_context("fork")
done = 0
with ctx.Pool(N_JOBS) as pool:
    for start, chunk in pool.imap_unordered(compute_chunk, chunk_args, chunksize=1):
        end = min(start + BATCH, N)
        D[start:end] = chunk
        done += 1
        if done % 10 == 0 or done == len(chunk_args):
            elapsed = time.time() - t0
            eta = elapsed * (len(chunk_args) - done) / done
            print(f"  [{done:>3}/{len(chunk_args)}] 已 {elapsed/60:.1f} min  ETA {eta/60:.1f} min")

print(f"\n  并行总耗时 {(time.time()-t0)/60:.1f} min")
print(f"  D shape: {D.shape}  dtype: {D.dtype}  mem: {D.nbytes/1024**3:.2f} GB")
del P, _P_GLOBAL

# %% [markdown]
# ## §3 Sanity check

# %%
print("\nSanity checks:")
print(f"  对角 max: {np.diag(D).max():.6e}  (应 ≈ 0)")
print(f"  非负: {(D >= -1e-6).all()}")
print(f"  对称偏差: {np.abs(D - D.T).max():.6e}")
upper = D[np.triu_indices_from(D, k=1)]
print(f"  非对角 min: {upper.min():.6f}")
print(f"  非对角 max: {upper.max():.6f}")
print(f"  非对角 median: {np.median(upper):.4f}")
print(f"  NaN/Inf: {np.isnan(D).sum()} / {np.isinf(D).sum()}")
assert np.isnan(D).sum() == 0 and np.isinf(D).sum() == 0
assert (D <= 1.0 + 1e-4).all(), "BC 应 ≤ 1"

# 强制对称（消除 fp 舍入）
D = (D + D.T) / 2
np.fill_diagonal(D, 0.0)

# %% [markdown]
# ## §4 写回 obsp（float16 压缩）

# %%
print(f"\n转 float16 + 写回 obsp['distance_bc'] ...")
D16 = D.astype(np.float16)
print(f"  float16 mem: {D16.nbytes/1024**3:.2f} GB")
adata.obsp["distance_bc"] = D16
del D

# %%
print(f"\n写回 {ANN_PATH.name} (compression=gzip) ...")
adata.write_h5ad(ANN_PATH, compression="gzip")
print(f"  新文件大小: {ANN_PATH.stat().st_size/1024**3:.2f} GB")

# %% [markdown]
# ## §5 读回验证

# %%
b = ad.read_h5ad(ANN_PATH, backed="r")
print(f"\nobsp keys: {list(b.obsp.keys())}")
arr = np.asarray(b.obsp["distance_bc"][:5, :5])
print(f"distance_bc[:5,:5] (dtype={arr.dtype}):")
print(arr)
