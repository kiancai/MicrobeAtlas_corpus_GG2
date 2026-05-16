# SampleDistance — 50k 子集 sample×sample 距离矩阵

挑 50k 代表性样本，算两套距离（**Bray-Curtis** + **Weighted Normalized UniFrac**），
再做 PCoA。1.8M 全样本上算 pairwise distance 既存不下也算不动，是规范做法。

## 流水线

| 脚本 | 输入 | 输出 |
|---|---|---|
| `01_stratified_sample.py` | `merged.gg2.with_phylo.h5ad` (1.83M × 8,114) | `subset_50k_index.tsv` |
| `02_build_subset_anndata.py` | + index TSV | `subset_50k.h5ad` (50k × 8,114，X + obs + varp) |
| `03_collapse_tree.py` | GG2 `phylogeny.id.nwk` (23M tip) + vocab | `genus_tree.nwk` (8,114-tip folded) |
| `04_compute_bc.py` | subset_50k.h5ad | → 写回 `obsp['distance_bc']` (float16) |
| `05_compute_wunifrac.py` | subset_50k.h5ad + genus_tree.nwk | → 写回 `obsp['distance_wunifrac']` (float16) |
| `06_pcoa.py` | subset_50k.h5ad + 两套 obsp | → 写回 `obsm['X_pcoa_*']` + `pcoa_eigenvalues.tsv` |

## 抽样配额（50,000 总量，RANDOM_SEED=42）

MA 30,000，一级桶按 `animal > soil > aquatic > plant` 优先级互斥归属：

| 桶 | 含义 | 配额 | 二级分层 |
|---|---|---:|---|
| MA::Human | `MA_IsHuman == 'Human'` | 8,000 | MA_SampleSite (gut/oral/skin/…) |
| MA::Animal_other | Animal=True 且非 pure Human | 8,099 | MA_Env_Animal_Sub top-10 + Other |
| MA::Soil | Soil=True 且非 Animal | 5,184 | MA_Env_Soil_Sub top-10 + Other |
| MA::Aquatic | Aquatic=True 且非 Animal/Soil | 4,969 | MA_Env_Aquatic_Sub top-10 + Other |
| MA::Plant | Plant=True 且非 Animal/Soil/Aquatic | 2,248 | MA_Env_Plant_Sub top-10 + Other |
| MA::Unknown | 无任何 env flag | 1,500 | 桶内随机 |

RM 20,000，按 `RM_Sample_Site` × `Project_ID` 两层 sqrt(N) 分配。

跨库 Run 重复（32,698 个 Run × 2 副本）不去重，两库独立抽样自然产生 ~150 个跨库重叠。

## 距离方案

- **Bray-Curtis (BC)**：经典组成距离；输入相对丰度，**不依赖**树或 genus 距离矩阵。
- **Weighted Normalized UniFrac**：phylogeny-aware；依赖 `genus_tree.nwk`（GG2 phylogeny 折叠到
  8,114 g__ tip，patristic 距离一字不差）。用 McDonald 2018 的 Striped Fast UniFrac
  (`unifrac.weighted_normalized`)。
- **不**做 rarefaction：用 relative abundance 即可（McMurdie & Holmes 2014）。

## 产物

```
results/sample_distance/
├── subset_50k_index.tsv            ~2 MB
├── subset_50k.h5ad                 ~13 GB
│   ├── X (CSR int32)               ~340 MB
│   ├── obs (56 列)                 ~30 MB
│   ├── varp['taxo_dist'/'phylo_dist']     ~310 MB (从主表搬来)
│   ├── obsp['distance_bc']         50k² float16  ~5 GB
│   ├── obsp['distance_wunifrac']   50k² float16  ~5 GB
│   ├── obsm['X_pcoa_bc']           50k × 10      ~2 MB
│   └── obsm['X_pcoa_wunifrac']     50k × 10      ~2 MB
├── genus_tree.nwk                  0.4 MB
└── pcoa_eigenvalues.tsv            < 1 KB
```

主 `merged.gg2.with_phylo.h5ad` 不变。

## 跑流水线

01/02/03 轻量，直接本地跑；04/05/06 计算量大，走 sbatch。

```bash
cd /hpcdisk1/limk_group/caiqy/project/260428_greengene2

# 轻量步：本地跑
PYBIO=/hpcdisk1/limk_group/caiqy/miniforge3/envs/baseBio/bin/python
$PYBIO scripts/SampleDistance/01_stratified_sample.py     # < 1 min
$PYBIO scripts/SampleDistance/02_build_subset_anndata.py  # ~2 min
$PYBIO scripts/SampleDistance/03_collapse_tree.py         # ~1 min

# 重量步：sbatch（必须显式带资源 flag，即使 #SBATCH 已写）
# 注意 --mem ≤ 256G 才能过 QOSMaxMemoryPerJob；BC/UniFrac/PCoA 实际只用 30-50 GB
sbatch -c 64 --mem=256G -t 02:00:00 scripts/SampleDistance/04_compute_bc.sbatch
sbatch -c 64 --mem=256G -t 02:00:00 scripts/SampleDistance/05_compute_wunifrac.sbatch
sbatch -c 64 --mem=256G -t 02:00:00 scripts/SampleDistance/06_pcoa.sbatch
```

各步均带 sanity / cross-check assert，失败会立即 raise。

## 设计要点（一句话版）

- **不动主 anndata**：所有距离产物在 `results/sample_distance/` 下，子集 anndata standalone。
- **float16 距离矩阵**：5 GB / 份，BC ∈ [0,1] / wUniFrac ∈ [0,~1]，精度足够。
- **树折叠 cross-check**：从新树重算 20 对 patristic 距离 vs `genus_phylo_dist.npz`，必须 diff < 1e-3。
- **嵌套 g__**：GG2 有 1,357 个 g__ 是另一个 g__ 的后代（GTDB 子属重命名），折叠时插
  零枝长 self-tip 处理，所有 8,114 g__ 都成为 tip 且距离不变。
- **抽样固定 seed=42**：可重现。
