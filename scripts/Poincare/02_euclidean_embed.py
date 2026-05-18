#!/usr/bin/env python3
"""
欧氏空间嵌入：GG2 genus → R^D 普通向量。

动机：双曲嵌入虽然 Pearson 高，但下游"MLP + 加到 token embedding"管线是天然欧氏的，
mismatch 大。欧氏嵌入在 D=16~64 时 Pearson 同样能 > 0.95，且管线全程欧氏无 mismatch。

实现：普通 Adam + 欧氏距离 + 绝对 MSE loss
  - d_eucl = ||x_i - x_j||_2
  - Loss = mean((d_eucl - d_tree_norm)²)
  - d_tree 归一化到 [0, TARGET_MAX]（只影响 lr 数值，不影响 Pearson）

用法：
    conda run -n poincare python 02_euclidean_embed.py [--dim 32] [--steps 20000]
"""

import argparse
import time
from pathlib import Path

import numpy as np
import torch

ROOT     = Path(__file__).resolve().parents[2]
DIST_IN  = ROOT / "results/phylogeny/genus_phylo_dist.npz"
VOCAB_IN = ROOT / "results/phylogeny/genus_vocab.tsv"


def pearsonr_np(x, y):
    x = x - x.mean(); y = y - y.mean()
    return float((x*y).sum() / (np.sqrt((x**2).sum() * (y**2).sum()) + 1e-30))


def validate(coords_np, d_mat_np, d_scale, n_sample=100_000, seed=0):
    N = len(coords_np)
    rng = np.random.default_rng(seed)
    pi = rng.integers(0, N, n_sample*3); pj = rng.integers(0, N, n_sample*3)
    m  = pi != pj;  pi, pj = pi[m][:n_sample], pj[m][:n_sample]

    dt_orig = d_mat_np[pi, pj].astype(np.float64)
    dt_norm = dt_orig * d_scale

    ct = torch.from_numpy(coords_np)
    with torch.no_grad():
        de = (ct[pi] - ct[pj]).norm(dim=-1).numpy()

    nz = dt_orig > 0
    pearson = pearsonr_np(de, dt_orig)
    mad     = float(np.mean(np.abs(de[nz] - dt_norm[nz]) / (dt_norm[nz] + 1e-9)))

    norms = np.linalg.norm(coords_np, axis=1)
    norm_stats = {k: float(v) for k, v in zip(
        ['min','median','p95','max'],
        [norms.min(), np.median(norms), np.percentile(norms,95), norms.max()]
    )}
    return pearson, mad, norm_stats


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dim',        type=int,   default=32)
    parser.add_argument('--steps',      type=int,   default=20_000)
    parser.add_argument('--lr',         type=float, default=0.05)
    parser.add_argument('--batch',      type=int,   default=10_000)
    parser.add_argument('--target_max', type=float, default=10.0,
                        help='归一化后 d_tree 的最大值（只影响 lr 数值，对 Pearson 无影响）')
    args = parser.parse_args()

    D       = args.dim
    STEPS   = args.steps
    LR      = args.lr
    BATCH   = args.batch
    TARGET  = args.target_max

    out_path = ROOT / f"results/poincare/euclidean_d{D}.npz"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    torch.set_default_dtype(torch.float32)

    print(f"=== 欧氏嵌入 (Adam + L2 距离) ===")
    print(f"    dim={D}  steps={STEPS}  lr={LR}  target_max={TARGET}")
    print(f"    torch={torch.__version__}\n")

    # ── 1. 加载距离矩阵 ──────────────────────────────────────────────────────
    print(f"加载距离矩阵: {DIST_IN.name} ...")
    t0 = time.time()
    npz      = np.load(DIST_IN)
    d_mat_np = npz['dist'].astype(np.float32)
    var_ids  = list(npz['var_id'])
    N        = len(var_ids)
    d_max    = float(d_mat_np.max())
    d_scale  = TARGET / d_max
    print(f"  N={N}  d_max={d_max:.1f}  归一化比例={d_scale:.5f}  "
          f"d_tree 归一化到 [0, {TARGET}]  耗时 {time.time()-t0:.1f}s")

    d_mat = torch.from_numpy(d_mat_np * d_scale)

    n_zero = int(((d_mat == 0) & ~torch.eye(N, dtype=torch.bool)).sum().item()) // 2
    print(f"  非对角零距离 pair: {n_zero} 对（不参与 loss）\n")

    # ── 2. 初始化 ────────────────────────────────────────────────────────────
    torch.manual_seed(42)
    coords = torch.nn.Parameter(
        torch.randn(N, D, dtype=torch.float32) * 1.0
    )

    optimizer = torch.optim.Adam([coords], lr=LR)
    scheduler = torch.optim.lr_scheduler.MultiStepLR(
        optimizer, milestones=[10_000, 16_000], gamma=0.3
    )

    # ── 3. 训练 ──────────────────────────────────────────────────────────────
    print(f"开始优化 ({STEPS:,} 步，每步 {BATCH:,} 对) ...")
    rng        = np.random.default_rng(0)
    best_loss  = float('inf')
    best_state = coords.data.clone()
    plateau    = 0
    t_start    = time.time()

    for step in range(1, STEPS + 1):
        pi = torch.from_numpy(rng.integers(0, N, BATCH*2).astype(np.int64))
        pj = torch.from_numpy(rng.integers(0, N, BATCH*2).astype(np.int64))
        mask = pi != pj
        pi, pj = pi[mask][:BATCH], pj[mask][:BATCH]

        dt = d_mat[pi, pj]
        nz = dt > 0
        if nz.sum() < 10:
            continue

        de = (coords[pi[nz]] - coords[pj[nz]]).norm(dim=-1)
        loss = ((de - dt[nz]) ** 2).mean()

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        scheduler.step()

        loss_val = float(loss.item())
        if loss_val < best_loss:
            best_loss  = loss_val
            best_state = coords.data.clone()
            plateau    = 0
        else:
            plateau += 1

        if plateau >= 2000 and step > 5000:
            print(f"\n  step {step}: plateau {plateau} 步，提前停止")
            break

        if step % 1000 == 0:
            elapsed = time.time() - t_start
            nmax = coords.data.norm(dim=-1).max().item()
            print(f"  step {step:6d}  loss={loss_val:.4f}  "
                  f"best={best_loss:.4f}  norm_max={nmax:.4f}  {elapsed:.0f}s")

    print(f"\n优化完成  best_loss={best_loss:.4f}  "
          f"总耗时: {time.time()-t_start:.0f}s\n")

    # ── 4. 验证 ──────────────────────────────────────────────────────────────
    print("完整验证 (100k 对) ...")
    final = best_state.numpy()
    pearson, mad, norm_stats = validate(final, d_mat_np, d_scale)

    print(f"\n{'='*64}")
    print(f"验证报告  dim={D}  steps={STEPS}")
    print(f"{'='*64}")
    print(f"  Pearson(d_eucl, d_tree_orig): {pearson:.4f}  "
          f"{'✓ >0.95' if pearson>0.95 else '△ >0.90' if pearson>0.90 else '✗ <0.90'}")
    print(f"  MAD (归一化空间)             : {mad:.4f}      "
          f"{'✓ <0.10' if mad<0.10 else '△ 偏高'}")
    print(f"  坐标范数:")
    for k, v in norm_stats.items():
        print(f"    {k:8s} = {v:.4f}")
    print(f"{'='*64}")

    # 生物合理性
    import csv
    vocab = {}
    with open(VOCAB_IN, newline='') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            vocab[row['var_id']] = row
    idx = {v: i for i, v in enumerate(var_ids)}
    fam_groups = {}
    for vid, row in vocab.items():
        fam_groups.setdefault(row['Family'], []).append(vid)

    ct = torch.from_numpy(final)
    same_d, cross_d = [], []
    print("\n  生物合理性检查:")
    with torch.no_grad():
        for fam, vlist in list(fam_groups.items()):
            if len(vlist) >= 2 and len(same_d) < 5:
                i1, i2 = idx[vlist[0]], idx[vlist[1]]
                d = float((ct[i1] - ct[i2]).norm().item())
                same_d.append(d)
                print(f"    同 Family ({fam.split(';')[-1]}): "
                      f"{vlist[0].split(';')[-1]} ↔ {vlist[1].split(';')[-1]}  d_eucl={d:.4f}")
        fams = list(fam_groups.keys())
        for i in range(0, len(fams)-1, max(1, len(fams)//10)):
            if len(cross_d) >= 5: break
            g1, g2 = fam_groups[fams[i]][0], fam_groups[fams[i+1]][0]
            i1, i2 = idx[g1], idx[g2]
            d = float((ct[i1] - ct[i2]).norm().item())
            cross_d.append(d)
            print(f"    跨 Family: {g1.split(';')[-1]} ↔ {g2.split(';')[-1]}  d_eucl={d:.4f}")
    if same_d and cross_d:
        ok = max(same_d) < min(cross_d)
        print(f"  同 Family max={max(same_d):.4f}，跨 Family min={min(cross_d):.4f}  "
              f"→ {'✓ 通过' if ok else '⚠ 有重叠（正常：同 Family 内部距离离散度本身就大）'}")

    # ── 5. 保存 ──────────────────────────────────────────────────────────────
    np.savez_compressed(
        out_path,
        coords  = final,
        var_id  = np.array(var_ids),
        dim     = np.int32(D),
        scale   = np.float64(d_scale),
        d_max   = np.float64(d_max),
        pearson = np.float32(pearson),
        mad     = np.float32(mad),
        method  = np.array('Adam + Euclidean L2 distance'),
    )
    print(f"\n✓ 已保存: {out_path}")
    print(f"  coords.shape={final.shape}  scale={d_scale:.5f}")
    print(f"  (d_eucl ≈ 原始 patristic distance × {d_scale:.5f})")


if __name__ == '__main__':
    main()
