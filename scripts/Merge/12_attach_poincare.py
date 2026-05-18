#!/usr/bin/env python3
"""
把 genus 位置编码（欧氏嵌入坐标）写入 MCFCorpus，生成 MCFCorpusV2。

用法:
    python 12_attach_poincare.py [--dim 32]

输入:
    results/feature_table/MCFCorpus.gg2.h5ad
    results/poincare/euclidean_d{D}.npz

输出:
    results/feature_table/MCFCorpusV2.gg2.h5ad
      - varm['position_encoding']: (8114, D) float32
      - uns['provenance']['position_encoding']: 训练元数据
"""

import argparse
import time
from pathlib import Path

import numpy as np
import anndata as ad

ROOT = Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dim', type=int, default=32)
    args = parser.parse_args()

    npz_path   = ROOT / f"results/poincare/euclidean_d{args.dim}.npz"
    h5ad_in    = ROOT / "results/feature_table/MCFCorpus.gg2.h5ad"
    h5ad_out   = ROOT / "results/feature_table/MCFCorpusV2.gg2.h5ad"

    print(f"=== 12_attach_position_encoding  dim={args.dim} ===\n")

    # ── 加载 npz ──────────────────────────────────────────────────────────────
    print(f"加载欧氏位置编码: {npz_path.name} ...")
    data = np.load(npz_path, allow_pickle=True)
    coords  = data['coords']      # (8114, D) float32
    var_ids = data['var_id']      # (8114,) str
    scale   = float(data['scale'])
    pearson = float(data['pearson'])
    mad     = float(data['mad'])
    print(f"  coords.shape = {coords.shape}  scale={scale:.5f}  Pearson={pearson:.4f}")

    # ── 加载 AnnData ──────────────────────────────────────────────────────────
    print(f"\n加载 AnnData: {h5ad_in.name} ...")
    t0 = time.time()
    adata = ad.read_h5ad(h5ad_in)
    print(f"  shape = {adata.shape}  耗时 {time.time()-t0:.1f}s")

    # ── 按 var_names 对齐坐标 ─────────────────────────────────────────────────
    adata_var_set = set(adata.var_names)
    npz_var_set   = set(var_ids.tolist())

    only_adata = adata_var_set - npz_var_set
    only_npz   = npz_var_set - adata_var_set

    if only_adata:
        print(f"  ⚠ AnnData 独有 var（在 npz 中缺失）: {len(only_adata)}")
        for v in list(only_adata)[:5]:
            print(f"    {v}")
    if only_npz:
        print(f"  ⚠ npz 独有 var（在 AnnData 中缺失）: {len(only_npz)}")

    varid_to_npz_idx = {v: i for i, v in enumerate(var_ids.tolist())}

    aligned = np.zeros((adata.n_vars, coords.shape[1]), dtype=np.float32)
    n_matched = 0
    for i, var_name in enumerate(adata.var_names):
        if var_name in varid_to_npz_idx:
            aligned[i] = coords[varid_to_npz_idx[var_name]]
            n_matched += 1

    print(f"  对齐: {n_matched}/{adata.n_vars} var 匹配成功")
    assert n_matched == adata.n_vars, f"对齐失败: 只有 {n_matched}/{adata.n_vars} 匹配"

    # ── 写入 varm ─────────────────────────────────────────────────────────────
    adata.varm['position_encoding'] = aligned
    print(f"  varm['position_encoding'].shape = {aligned.shape}")

    # ── 更新 provenance ───────────────────────────────────────────────────────
    if 'provenance' not in adata.uns:
        adata.uns['provenance'] = {}
    method = str(data['method']) if 'method' in data.files else 'Adam + Euclidean L2 distance'
    adata.uns['provenance']['position_encoding'] = {
        'source':  'scripts/Poincare/02_euclidean_embed.py',
        'method':  method,
        'space':   'euclidean',
        'dim':     args.dim,
        'scale':   scale,
        'pearson': pearson,
        'mad':     mad,
        'npz':     str(npz_path),
        'note':    'd_eucl(coords[i], coords[j]) ≈ phylo_dist[i,j] × scale',
    }

    # ── 保存 ──────────────────────────────────────────────────────────────────
    print(f"\n写出: {h5ad_out.name} ...")
    t0 = time.time()
    adata.write_h5ad(h5ad_out, compression="gzip")
    print(f"  完成  耗时 {time.time()-t0:.1f}s")
    print(f"\n✓ MCFCorpusV2 已保存: {h5ad_out}")


if __name__ == '__main__':
    main()
