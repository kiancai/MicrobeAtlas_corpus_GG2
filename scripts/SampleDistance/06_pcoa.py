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
# # 06: 两套距离矩阵的 PCoA 坐标
#
# 输入：
# - `results/sample_distance/subset_50k.h5ad`：含 `obsp['distance_bc']` + `obsp['distance_wunifrac']`
#
# 输出：
# - `obsm['X_pcoa_bc']`        (50000, 10) float32
# - `obsm['X_pcoa_wunifrac']`  (50000, 10) float32
# - `results/sample_distance/pcoa_eigenvalues.tsv`：每个距离的前 10 个特征值 + 解释方差比
#
# ## 算法
#
# 50k × 50k dense matrix 的完整 eigendecomp 不可行（10 GB float64 + 特征分解 O(N^3)）。
# 改用 **truncated PCoA**：scipy.sparse.linalg.eigsh 只算 top-K 特征值/特征向量。
# - K=10
# - 输入：double-centered Gram matrix（PCoA 标准做法）
# - 用 LinearOperator + eigsh，避免显式构造 N×N float64（仍需 dense N×N float32 ≈ 10 GB，可接受）
#
# float16 距离矩阵反取 → float32 计算，避免精度损失。

# %%
from pathlib import Path
import time
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse.linalg import LinearOperator, eigsh

ad.settings.allow_write_nullable_strings = True

ROOT = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
ANN_PATH = ROOT / "results/sample_distance/subset_50k.h5ad"
EIG_OUT  = ROOT / "results/sample_distance/pcoa_eigenvalues.tsv"
K = 10  # 取前 10 维 PCoA

# %% [markdown]
# ## §1 读 anndata（in-memory）

# %%
print(f"读 {ANN_PATH.name} ...")
adata = ad.read_h5ad(ANN_PATH)
print(f"  shape: {adata.shape}")
print(f"  obsp keys: {list(adata.obsp.keys())}")

# %% [markdown]
# ## §2 PCoA via truncated eigensolve
#
# 标准 PCoA：
# 1. D = 距离矩阵；A = -D²/2
# 2. B = A - mean_row - mean_col + mean_all   (double centering)
# 3. eigendecompose B；取 top-K 正特征值 λ 和特征向量 v
# 4. 坐标 = v * sqrt(λ)
#
# 这里 N=50k，B 是 dense float32 = 10 GB。可以承受。eigsh 给 sparse 或 LinearOperator
# 都行；用 LinearOperator 避免重复读 B。但 dense B 已经显式存在，直接给 eigsh 也行。

# %%
def pcoa_truncated(D: np.ndarray, k: int = 10):
    """对距离矩阵做 truncated PCoA，返回 (coords, eigvals, explained_var_ratio)"""
    t0 = time.time()
    N = D.shape[0]
    print(f"  N={N}, k={k}")
    # 取 float32 计算
    A = -0.5 * (D.astype(np.float32) ** 2)
    row_mean = A.mean(axis=1, keepdims=True)
    col_mean = A.mean(axis=0, keepdims=True)
    all_mean = A.mean()
    # In-place 减以省内存
    A -= row_mean
    A -= col_mean
    A += all_mean
    B = A  # 对称
    print(f"    double-center 耗时 {time.time() - t0:.1f}s  mem {B.nbytes/1024**3:.2f} GB")

    # trace(B) = sum 所有特征值 = total variance（PCoA 标准分母）
    trace_B = float(np.trace(B.astype(np.float64)))
    print(f"    trace(B) = {trace_B:.4f}")

    # eigsh 取最大 k 个正特征值（PCoA 需正）
    t1 = time.time()
    eigvals, eigvecs = eigsh(B, k=k + 5, which="LA")
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    print(f"    eigsh 耗时 {time.time() - t1:.1f}s")
    print(f"    top {k + 5} eigvals: {eigvals}")

    pos_mask = eigvals > 0
    if pos_mask.sum() < k:
        print(f"  ⚠️  正特征值只有 {pos_mask.sum()} 个，< k={k}")
    sel = np.arange(min(k, len(eigvals)))
    eigvals_sel = np.where(eigvals[sel] > 0, eigvals[sel], 0.0)
    coords = eigvecs[:, sel] * np.sqrt(eigvals_sel)[None, :]
    coords = coords.astype(np.float32)
    # explained variance ratio = eigval / trace(B)
    explained = eigvals_sel / trace_B if trace_B > 0 else np.zeros_like(eigvals_sel)
    return coords, eigvals_sel, explained


# %% [markdown]
# ## §3 跑 BC 距离的 PCoA

# %%
print("\n=== PCoA: distance_bc ===")
D_bc = np.asarray(adata.obsp["distance_bc"]).astype(np.float32)
print(f"D_bc 取出 dtype: {D_bc.dtype}  shape: {D_bc.shape}")
# 强制对称
D_bc = (D_bc + D_bc.T) / 2
np.fill_diagonal(D_bc, 0.0)
coords_bc, eig_bc, expl_bc = pcoa_truncated(D_bc, k=K)
print(f"\nBC top-{K} eigvals: {eig_bc}")
print(f"BC explained var ratio: {expl_bc}")
adata.obsm["X_pcoa_bc"] = coords_bc
del D_bc

# %% [markdown]
# ## §4 跑 weighted UniFrac 距离的 PCoA

# %%
print("\n=== PCoA: distance_wunifrac ===")
D_wu = np.asarray(adata.obsp["distance_wunifrac"]).astype(np.float32)
print(f"D_wu 取出 dtype: {D_wu.dtype}  shape: {D_wu.shape}")
D_wu = (D_wu + D_wu.T) / 2
np.fill_diagonal(D_wu, 0.0)
coords_wu, eig_wu, expl_wu = pcoa_truncated(D_wu, k=K)
print(f"\nwUniFrac top-{K} eigvals: {eig_wu}")
print(f"wUniFrac explained var ratio: {expl_wu}")
adata.obsm["X_pcoa_wunifrac"] = coords_wu
del D_wu

# %% [markdown]
# ## §5 写回 obsm + 落 eigenvalues.tsv

# %%
print(f"\n写回 {ANN_PATH.name} ...")
adata.write_h5ad(ANN_PATH, compression="gzip")
print(f"  大小: {ANN_PATH.stat().st_size/1024**3:.2f} GB")

# %%
df = pd.DataFrame({
    "axis":     list(range(1, K + 1)),
    "bc_eigval":          eig_bc,
    "bc_explained_var":   expl_bc,
    "wunifrac_eigval":         eig_wu,
    "wunifrac_explained_var":  expl_wu,
})
df.to_csv(EIG_OUT, sep="\t", index=False, float_format="%.6f")
print(f"\n已写出 {EIG_OUT}")
print(df.to_string(index=False))

# %% [markdown]
# ## §6 读回验证

# %%
b = ad.read_h5ad(ANN_PATH, backed="r")
print(f"\nobsm keys: {list(b.obsm.keys())}")
print(f"  X_pcoa_bc shape:        {b.obsm['X_pcoa_bc'].shape}")
print(f"  X_pcoa_wunifrac shape:  {b.obsm['X_pcoa_wunifrac'].shape}")
print(f"\nX_pcoa_bc[:3]:")
print(np.asarray(b.obsm["X_pcoa_bc"][:3]))
