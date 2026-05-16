# 50k 分层子集 sample × sample 距离矩阵 详档

`scripts/SampleDistance/`，产物在 `results/sample_distance/`。

## 1. 为什么是子集

1.83M × 1.83M 距离矩阵 float32 = 13 TB，float16 = 6.6 TB；都存不下。
经典 PCoA 还要做 1.83M × 1.83M 矩阵特征分解，根本不可行。
业界 100k+ 微生物组数据默认是**分层抽样代表性子集 + PCoA**（MGnify、AGP meta-analysis 等都是这做法）。

## 2. 抽样配额（01_stratified_sample.py）

总量 50,000 = MA 30,000 + RM 20,000。`RANDOM_SEED = 42`。

### MA 一级桶（优先级互斥）

按 `MA_Env_Animal > MA_Env_Soil > MA_Env_Aquatic > MA_Env_Plant` 顺序判主类；
`MA_IsHuman == 'Human'` 单独切出来作 Human 桶；`HumanMix` 归 Animal_other。

| 桶 | 判定 | 总样本 | 配额 |
|---|---|---:|---:|
| Human          | `MA_IsHuman == 'Human'` | 334,115 | 8,000 (user-fixed) |
| Animal_other   | Animal=True 且非 pure Human (含 HumanMix) | 567,424 | 8,099 (sqrt) |
| Soil           | Soil=True 且非 Animal | 232,434 | 5,184 (sqrt) |
| Aquatic        | Aquatic=True 且非 Animal/Soil | 213,335 | 4,969 (sqrt) |
| Plant          | Plant=True 且非 Animal/Soil/Aquatic | 43,588 | 2,248 (sqrt) |
| Unknown        | 无任何 env flag | 341,805 | 1,500 (user-fixed) |

sqrt 配额：4 个非 fixed 桶共 20,500 名额按 `sqrt(N_i) / Σ sqrt(N)` 加权。

### MA 二级分层

各桶内继续做 sqrt 分层：
- `Human` → `MA_SampleSite` 8 部位（gut/oral/skin/urogenital/lung/gastric/nose/bone）+ NA
- `Animal_other` → `MA_Env_Animal_Sub` top-10 子类 + Other
- `Soil/Aquatic/Plant` → 对应 `_Sub` top-10 + Other
- `Unknown` → 桶内随机

### RM 分层

按 `RM_Sample_Site × Project_ID` 两层。一级 sqrt，二级 sqrt 内项目。RM 主要是呼吸道（Nasopharynx 24k / Sputum 15k / Nasal 14.6k / Oropharynx 13k / BALF 7k / …）。

### 跨库 Run 重复

不去重。MA 和 RM 独立抽样，期望 ~150 个 Run 被两边同时抽到，下游可以用 `Database` 标签区分。

## 3. 子集 anndata 结构（02_build_subset_anndata.py）

```
subset_50k.h5ad  (6.3 GB after 04+05+06)
├── X (CSR int32, 50k × 8,114, 5.9M nnz)            340 MB
├── obs (56 列 = 主表 54 + stratum_id + sub_stratum) 30 MB
├── var (7 列 = 6 rank + observed)
├── varp['taxo_dist']  8114² int8                   63 MB（搬主表）
├── varp['phylo_dist'] 8114² float32               251 MB（搬主表）
├── obsp['distance_bc']        50k² float16        ~5 GB（04 产物）
├── obsp['distance_wunifrac']  50k² float16        ~5 GB（05 产物）
├── obsm['X_pcoa_bc']          50k × 10 float32     2 MB
└── obsm['X_pcoa_wunifrac']    50k × 10 float32     2 MB
```

`obs_names` 沿用主表（`MA_xxxxxxx` / `RM_xxxxxxx`），可反查主表。

## 4. Genus 树折叠（03_collapse_tree.py）

输入：`tools/qiime2_database/greengenes2/tree/2024.09.phylogeny.id.nwk`（23M tip / ~1M 内部节点）。
输出：`genus_tree.nwk` 8,114-tip newick（0.4 MB）。

### 算法

1. 手写 newick 解析器（同 `Phylogeny/03_compute_phylogenetic_distance.py`）：只构建内部节点，跳过 23M tip
2. 标 `is_genus`（label 含 g__）+ 自底向上算 `has_g_desc`
3. 迭代式后序 emit newick：
   - 真叶 g__（无 g__ 后代）→ `g__Xxx:length`
   - 非 g__ 内部 → `(child1,child2,...)label:length`（只递归有 g__ desc 的孩子）
   - **嵌套 g__**（自身是 g__ 且还有 g__ 后代）→ 拆成合成内部节点：`(g__Xxx:0,child1,child2,...):length`
4. 内部节点 label 含 `; `、空格等特殊字符 → 单引号包裹

### 嵌套 g__ 处理（关键 trick）

GG2 树里有 1,357 个 g__ 节点是另一个 g__ 节点的后代（GTDB 子属重命名，如 `g__Bacteroides` 下有 `g__Bacteroides_A`）。
如果遇到 g__ 就把子树全砍掉，会丢这 1,357 个嵌套 g__（实测先版本 tip 数 = 6,757 ≠ 8,114）。
解决：在嵌套 g__ 位置插入合成内部节点，子节点 = `[g__Xxx:0]` + 原 g__ 子树递归结果。
这样：
- 8,114 个 g__ 全部成为新树的 tip
- 合成节点继承原 g__ 节点的枝长 L 到 parent
- self-tip 到合成节点枝长 0
- 任意两 g__ 间 patristic 距离不变（直接走原路径）

### Cross-check

读回新树，用 skbio `TreeNode.read(..., convert_underscores=False)` 防 `g__` 被 skbio 当下划线转空格。
随机 20 对 g__ 算 `node.distance()`，与 `genus_phylo_dist.npz` 对比，diff < 1e-3 才允许通过。

## 5. Bray-Curtis（04_compute_bc.py + sbatch）

### 公式

raw int32 → 行除以 row sum → relative abundance dense float32 (50k × 8,114 = 1.51 GB)。
`BC(u,v) = sum|u-v| / sum(u+v)`，归一后 `sum(u+v) = 2`，所以 `BC = sum|u-v|/2 = 1 − sum_min(u,v)`。

### 并行

`multiprocessing.Pool` fork 模式（父进程 P 共享给子进程，COW 不复制），每个 worker 负责一片 BATCH=500 行调 `scipy.spatial.distance.cdist(Pi, P, 'braycurtis')`。100 块 × 1.5 GB P 共享。

`OMP_NUM_THREADS=1 / OPENBLAS_NUM_THREADS=1` 防止 BLAS 内部多线程跟 Pool 抢 CPU。

### 资源

- 本地 8 核：fork Pool 调度异常（4 worker 始终 idle），~80 min 才能跑完
- sbatch 64 核：fork 初始化 7 min，计算只 5 min；**total ~12 min**

写盘 float16 + gzip：5 GB → 压到 ~2.2 GB；对称 + 对角清零，sanity `BC ≤ 1`。

## 6. Weighted Normalized UniFrac（05_compute_wunifrac.py + sbatch）

用 **Striped Fast UniFrac**（McDonald 2018 + Sfiligoi 2022 mSystems），`unifrac.weighted_normalized()`：
- 输入 1：BIOM-Format v2.1 临时文件（features = 8,114 g__ token，samples = obs_names）
- 输入 2：`genus_tree.nwk`
- `OMP_NUM_THREADS = 64`

### 实测

50k × 50k weighted UniFrac **42 秒**完成（C++ SIMD striped + 64 线程）。同等规模 scikit-bio 单线程要 12+ 小时，差 1000 倍。

### var_id ↔ g__ 映射

anndata `var_names` 是完整 6 级路径（`d__...;p__...;...;g__Foo`），树 tip 是裸 g__ token。
建 BIOM 时用 `vocab.loc[var_ids, 'Genus']` 作 `observation_ids`，让 BIOM feature id 和树 tip 对得上。

### 距离范围

weighted normalized UniFrac ∈ [0, 1]。50k 子集实测：median 0.45，非对角 min 可能 = 0（跨库 Run 重复 + 偶然相同丰度向量）。

## 7. PCoA（06_pcoa.py + sbatch）

### 算法

1. `D` 取 float16 → float32，强制对称、对角清零
2. `A = -0.5 * D^2`
3. Double centering：`B = A - row_mean - col_mean + all_mean`（in-place 节省内存）
4. `trace(B) = sum(所有特征值)` 作分母
5. `scipy.sparse.linalg.eigsh(B, k=15, which='LA')` 取最大 15 个特征值；倒序取前 10 个正的
6. 坐标 = `eigvecs * sqrt(eigvals)`

50k × 50k float32 B = 9.3 GB，单步 in-place 减；eigsh 实测 **2.5 秒**（dense float32 + truncated）。

### 解释方差比

| 轴 | BC | wUniFrac |
|---:|---:|---:|
| PC1 | 10.9% | 20.1% |
| PC2 | 5.5% | 18.1% |
| PC3 | 4.4% | 9.1% |
| PC4 | 3.7% | 7.0% |
| PC5 | 3.0% | 5.4% |
| PC6 | 2.7% | 4.3% |
| PC7 | 2.4% | 3.8% |
| PC8 | 2.2% | 3.5% |
| PC9 | 2.0% | 2.1% |
| PC10 | 1.7% | 1.8% |
| top-3 cum | **20.8%** | **47.3%** |

wUniFrac 把信号集中到少数轴（PC1+PC2 共 38%），符合 phylo-aware 距离的预期：把"虽然 ASV/genus 不同但系统发育近"的方向折叠成同一主方向。

## 8. 可视化（07_plot_pcoa.py）

4 张 PNG 输出到 `results/sample_distance/figures/`：

| 文件 | 内容 |
|---|---|
| `pcoa_bucket.png` | 2×2 panel：BC PC12/PC13 + wUniFrac PC12/PC13；按 7 类大桶（Human/Animal/Soil/Aquatic/Plant/RM_Respiratory/Unknown）着色 |
| `pcoa_scree.png`  | BC vs wUniFrac 前 10 轴 explained variance 折线对比 |
| `pcoa_human_sites.png` | MA Human 8,000 内部按 `MA_SampleSite` 9 部位着色，BC + wUniFrac 各一 |
| `pcoa_rm_sites.png`    | RM 20,000 内部按 `RM_Sample_Site` top-12 + Other 着色 |

绘图细节：
- 大类按数量从大到小排序，小类后画（小类点不被盖住）
- `s=1.5, alpha=0.4, rasterized=True` → 50k 点画在 PNG 不会太大
- 共用图例放在 figure 顶部，markerscale=4 让 legend 点更清楚

### 主图速读

- **BC 散成扇形**：Human / Animal / Soil / Aquatic / Plant / RM 各占象限，PC1 主要是 host-associated（PC1<0）vs free-living（PC1>0）的分离
- **wUniFrac 抱成实心椭圆**：结构压缩在 PC1×PC2，phylo 信号集中
- **gut 子样本**（MA Human 内）在 BC 上聚成 PC1 负 / PC2 偏下；wUniFrac 上 gut vs oral/lung 分得更干净

## 9. 资源 & sbatch

QOSMaxMemoryPerJob：单 job `--mem ≤ 256G`，超会 PD 不调度。
BC/UniFrac/PCoA 实际只用 30-50 GB，256G 富余。

提交模板（必须显式 flag）：
```
sbatch -c 64 --mem=256G -t 02:00:00 scripts/SampleDistance/{04,05,06}_*.sbatch
```

`PYTHONUNBUFFERED=1` 让 print 实时落盘，方便 `tail -f logs/sd_*.log` 看进度。

## 10. 关键决策（一句话版）

- 子集 anndata standalone，**主 anndata 不动**
- float16 距离矩阵足够（BC/wUniFrac 都 ∈ [0,1]），换 5 GB/份 vs float32 10 GB/份
- Bray-Curtis（baseline 无争议）+ Weighted Normalized UniFrac（phylo 金标准）双跑作对照
- 不做 rarefaction，relative abundance 即可（McMurdie & Holmes 2014）
- 嵌套 g__ 用 0 枝长 self-tip 保留全 8,114 g__ 为 tip，patristic 距离不变
- PCoA 用 truncated eigsh top-10，避免完整特征分解 50k × 50k 不可行
- 1.8M 全样本上画 PCoA 本身不可行；想要全样本图后续应走 PCA(64 维) + UMAP 路线
