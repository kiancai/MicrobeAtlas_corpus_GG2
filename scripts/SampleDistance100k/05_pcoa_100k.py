from pathlib import Path
import os
import time
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse.linalg import eigsh

ad.settings.allow_write_nullable_strings = True

ROOT = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
OUT_DIR = ROOT / "results/sample_distance_100k"
ANN_PATH = OUT_DIR / "subset_100k.h5ad"
DIST_FILES = {
    "bc": OUT_DIR / "distance_bc.float16.npy",
    "wunifrac": OUT_DIR / "distance_wunifrac.float16.npy",
}
COORD_FILES = {
    "bc": OUT_DIR / "pcoa_coords_bc.npy",
    "wunifrac": OUT_DIR / "pcoa_coords_wunifrac.npy",
}
EIG_OUT = OUT_DIR / "pcoa_eigenvalues.tsv"
AUDIT_OUT = OUT_DIR / "pcoa_negative_eigen_audit.tsv"

K = int(os.environ.get("PCOA_K", "10"))
EXTRA = int(os.environ.get("PCOA_EXTRA", "5"))
SYMM_BATCH = int(os.environ.get("SYMM_BATCH", "2000"))
AUDIT_N = int(os.environ.get("AUDIT_N", "3000"))


def symmetrize_inplace(D: np.ndarray, batch: int) -> None:
    n = D.shape[0]
    for i in range(0, n, batch):
        i_end = min(i + batch, n)
        block = D[i:i_end, i:i_end]
        block_avg = (block + block.T) * 0.5
        D[i:i_end, i:i_end] = block_avg
        for j in range(i_end, n, batch):
            j_end = min(j + batch, n)
            a = D[i:i_end, j:j_end].copy()
            b = D[j:j_end, i:i_end].T.copy()
            avg = (a + b) * 0.5
            D[i:i_end, j:j_end] = avg
            D[j:j_end, i:i_end] = avg.T
        if (i // batch + 1) % 10 == 0 or i_end == n:
            print(f"    symmetrized rows {i_end:,}/{n:,}")
    np.fill_diagonal(D, 0.0)


def centered_gram_from_distance(D: np.ndarray) -> tuple[np.ndarray, float]:
    np.square(D, out=D)
    D *= -0.5
    row_mean = D.mean(axis=1, dtype=np.float64).astype(np.float32)
    all_mean = float(row_mean.mean(dtype=np.float64))
    D -= row_mean[:, None]
    D -= row_mean[None, :]
    D += all_mean
    trace = float(np.diag(D).astype(np.float64).sum())
    return D, trace


def audit_negative_spectrum(metric: str, dist_path: Path) -> dict:
    D16 = np.load(dist_path, mmap_mode="r")
    n = D16.shape[0]
    audit_n = min(AUDIT_N, n)
    rng = np.random.default_rng(20260516)
    idx = np.sort(rng.choice(n, audit_n, replace=False))
    D = np.asarray(D16[np.ix_(idx, idx)], dtype=np.float32)
    symmetrize_inplace(D, batch=min(1000, audit_n))
    B, trace = centered_gram_from_distance(D)
    top = np.sort(eigsh(B, k=min(10, audit_n - 2), which="LA", return_eigenvectors=False))[::-1]
    bottom = np.sort(eigsh(B, k=min(10, audit_n - 2), which="SA", return_eigenvectors=False))
    neg = bottom[bottom < 0]
    return {
        "metric": metric,
        "audit_n": audit_n,
        "trace": trace,
        "top10_sum_trace_ratio": float(top[top > 0].sum() / trace) if trace > 0 else 0.0,
        "bottom_min": float(bottom[0]),
        "bottom10_abs_sum": float(np.abs(neg).sum()),
        "bottom10": ",".join(f"{x:.6g}" for x in bottom),
    }


def pcoa_metric(metric: str, dist_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    print(f"\n=== PCoA {metric}: {dist_path.name} ===")
    t0 = time.time()
    D16 = np.load(dist_path, mmap_mode="r")
    assert D16.shape[0] == D16.shape[1], f"{dist_path} is not square"
    D = np.array(D16, dtype=np.float32, copy=True)
    print(f"  loaded D shape={D.shape} mem={D.nbytes / 1024**3:.2f} GB elapsed={time.time() - t0:.1f}s")

    t1 = time.time()
    symmetrize_inplace(D, SYMM_BATCH)
    print(f"  symmetrize elapsed={time.time() - t1:.1f}s")
    print(f"  range after sym: min={D.min():.6f} max={D.max():.6f}")
    assert np.isnan(D).sum() == 0 and np.isinf(D).sum() == 0

    t2 = time.time()
    B, trace = centered_gram_from_distance(D)
    print(f"  double-center elapsed={time.time() - t2:.1f}s trace={trace:.6f}")

    t3 = time.time()
    eigvals, eigvecs = eigsh(B, k=K + EXTRA, which="LA")
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    print(f"  eigsh elapsed={time.time() - t3:.1f}s")
    print(f"  top {K + EXTRA} eigvals: {eigvals}")

    eigvals_sel = np.where(eigvals[:K] > 0, eigvals[:K], 0.0)
    coords = eigvecs[:, :K] * np.sqrt(eigvals_sel)[None, :]
    coords = coords.astype(np.float32)
    explained = eigvals_sel / trace if trace > 0 else np.zeros_like(eigvals_sel)
    return coords, eigvals_sel, explained, trace


print(f"K={K} EXTRA={EXTRA} SYMM_BATCH={SYMM_BATCH} AUDIT_N={AUDIT_N}")
for path in DIST_FILES.values():
    assert path.exists(), f"Missing distance file: {path}"

audit_rows = []
eig_rows = []
coords_by_metric = {}
for metric, path in DIST_FILES.items():
    print(f"\nAuditing negative spectrum for {metric} ...")
    audit = audit_negative_spectrum(metric, path)
    audit_rows.append(audit)
    print(audit)

    coords, eigvals, explained, trace = pcoa_metric(metric, path)
    np.save(COORD_FILES[metric], coords)
    coords_by_metric[metric] = coords
    for axis in range(K):
        eig_rows.append({
            "metric": metric,
            "axis": axis + 1,
            "eigval": float(eigvals[axis]),
            "explained_trace_ratio": float(explained[axis]),
            "trace": float(trace),
        })
    print(f"  wrote coords: {COORD_FILES[metric]}")

pd.DataFrame(eig_rows).to_csv(EIG_OUT, sep="\t", index=False, float_format="%.8g")
pd.DataFrame(audit_rows).to_csv(AUDIT_OUT, sep="\t", index=False, float_format="%.8g")
print(f"\nWrote {EIG_OUT}")
print(f"Wrote {AUDIT_OUT}")

print("Writing coordinates into subset_100k.h5ad obsm ...")
adata = ad.read_h5ad(ANN_PATH)
adata.obsm["X_pcoa_bc"] = coords_by_metric["bc"]
adata.obsm["X_pcoa_wunifrac"] = coords_by_metric["wunifrac"]
adata.uns["sample_distance_100k"] = {
    "distance_bc": str(DIST_FILES["bc"].relative_to(ROOT)),
    "distance_wunifrac": str(DIST_FILES["wunifrac"].relative_to(ROOT)),
    "pcoa_coords_bc": str(COORD_FILES["bc"].relative_to(ROOT)),
    "pcoa_coords_wunifrac": str(COORD_FILES["wunifrac"].relative_to(ROOT)),
    "explained_variance_note": "explained_trace_ratio is eigval / trace(double-centered Gram); negative eigenvalues are audited separately.",
}
adata.write_h5ad(ANN_PATH, compression="gzip")
print(f"Updated {ANN_PATH}")
