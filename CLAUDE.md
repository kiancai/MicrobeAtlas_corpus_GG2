# GreenGene2 / MicrobeAtlas 项目

## 项目概述

从 MicrobeAtlas 数据库下载的微生物组数据，用于后续分析。两部分：
- **MicrobeAtlas**: 样本-OTU 丰度矩阵、OTU 注释信息、参考序列库（~269 万样本）
- **ResMicroDb**: 抗性微生物数据库的原始测序 reads（398 项目 / ~10 万样本 / ~164 万 ASV）

## 环境配置

- **Conda 环境**: `baseBio`（路径: `/hpcdisk1/limk_group/caiqy/miniforge3/envs/baseBio`）
- 已安装: `biopython`, `biom-format`, `h5py`, `pandas`, `cutadapt`, `fastp`, `numpy`, `pigz`, `pbzip2`, `scikit-bio`（仅 01 用）, `scipy`
- 代理: `proxy_on` 激活代理（`127.0.0.1:60124`），conda 安装前需先激活

## 数据目录结构

```
rawdata/
├── MicrobeAtlas/
│   ├── sample_info/
│   │   └── samples.env.info.tsv          # 样本元信息 (301 MB)
│   ├── OTU_info/
│   │   ├── otus.info.tsv                  # OTU 注释信息 (954 MB)
│   │   └── mapref-3.0.tar.gz             # 参考序列库 (315 MB)
│   └── OTU_count/
│       ├── otus.97.allinfo                # 97% OTU 详细信息 (197 MB)
│       ├── samples-otus.97.mapped.biom.gz          # 未过滤丰度矩阵 (3.6 GB)
│       └── samples-otus.97.mapped.metag.minfilter.refilt.biom.gz  # 过滤后丰度矩阵 (3.1 GB)
└── ResMicroDb/
    ├── 16S/<PROJECT>/                     # 398 个项目，~10 万样本，~164 万 ASV
    │   ├── results/                       # jxt 流水线产物 + 本项目 GG2 注释
    │   │   ├── asv.fa                     # 项目内 ASV (ID 从 ASV_1 起，跨项目重名)
    │   │   ├── otutab.txt                 # ASV × sample dense TSV
    │   │   ├── rep-seqs.qza, table.qza    # 上述两者的 QZA 形式
    │   │   ├── taxonomy.qza               # SILVA NB 分类 QZA（jxt）
    │   │   ├── taxonomy_silva.txt         # SILVA 注释 8 列 (OTUID + K..S)
    │   │   ├── taxonomy_gg2.qza           # GG2 NB 分类 QZA（本项目新增）
    │   │   └── taxonomy_gg2.txt           # GG2 注释 3 列 (Feature ID / Taxon / Confidence)
    │   ├── log/                           # jxt 原日志
    │   └── log_gg2/                       # GG2 分类日志（本项目新增）
    └── raw_reads/PRJEB13657/              # 原始双端 fastq（仅 1 个项目齐全，其它项目 fastq 不在本机）
```

> 各文件列含义、BIOM 读取代码、mapref 解压结构 → `.claude/data/file_schemas.md`
> Environments / body site 详细计数 → `.claude/data/sample_env_stats.md`
> OTU 嵌套层级（`90_x;96_x;97_x;...`）说明 → `.claude/data/otu_hierarchy.md`

## 数据规模摘要

| 指标 | 数量 |
|------|------|
| 总样本数（MicrobeAtlas） | ~269 万 |
| 过滤后样本数（minfilter） | ~188 万 |
| 97% OTU 数（MicrobeAtlas） | ~10.3 万 |
| 参考全长序列数（mapref-3.0） | ~136 万 |
| ResMicroDb 样本 / ASV | ~10 万 / ~164 万 |
| 主要测序技术 | Amplicon (88%)、WGS (10%) |

---

## 流水线脚本一览

```
scripts/MicrobeAtlas/        # OTU → genus 聚合（输入是策展的全长代表序列）
  01_extract_rep_seqs.sh           从 otus.97.allinfo 提取 97% 代表序列
  02_qiime2_classify.sbatch        GG2 NB 分类（用默认 auto，不要用 both）
  03_build_genus_anndata.py        过滤 mito/chloro/非 BA + 6 级路径聚合
  04_drop_empty_samples.py         剔除零 taxon 样本
  05_qc_filter.py                  删 shallow var + 迭代 QC（MIN_READS=1000, MIN_FEATURES=5）
  05_v2_qc_filter_keep_shallow.py  保留 shallow var 的版本
  06_attach_metadata.py            samples.env.info.tsv → obs 26 列 metadata（py:percent ↔ ipynb）
  07_filter_samples.py             删 Sequencing_Type ∈ {RNAseq, NaN}（保留 AMPLICON + WGS）

scripts/ResMicroDb/          # 398 项目 ASV → 跨项目 sample × genus
  01_qiime2_classify.sbatch        单项目 GG2 NB 分类（已加 --p-read-orientation both）
  01_run_loop.sh                   轮询提交器（CAP=90, INTERVAL=30s, 断点续投）
  02_merge_to_asv_anndata.py       398 项目 → sample × ASV CSR（GG2 + SILVA 双注释）
  03_build_genus_anndata.py        过滤 mito/chloro/非 BA + 6 级路径聚合
  04_drop_empty_samples.py         剔除零 taxon 样本（含项目维度报告）
  05_qc_filter.py                  同 MicrobeAtlas（MIN_FEATURES=5）
  06_standardize_metadata.py       metadata_all.txt → metadata_all.standardized.tsv（135,746 × 36；UTF-16LE → UTF-8；jxt 派生 Age_Group/Case_Or_Control/Is_Healthy + projectTable region）
  06b_export_patches.R             从 ps.16s_0105_new7.rds 导出 Run 级 patch tsv（仅供 Patch 1/2 消费，baseR + phyloseq）
  06b_fix_metadata_errors.py       v2 NCBI 真值方案：5 项目 1,697 Run 修正（P1 Sample_Site 1,257 + P2 PRJNA801796 Phenotype 细分 255 + P3 PRJNA822681 NCBI host disease 152 + P4 PRJNA824137 NCBI title 33）→ .standardized.fixed.parquet
  07_attach_metadata.py            metadata.fixed.parquet → anndata.obs 37 列（Database + Run 列 + 35 metadata；Run 既作 obs_names 又作列）
  08_filter_samples.py             passthrough 骨架（anndata 已 100% 16S；保留 Negative Control 等留给下游训练时筛）

scripts/Merge/               # 两库合并 + 挂距离矩阵 + corpus finalize
  09_merge_databases.py            outer-join var；obs_names 重编号 MA_/RM_ 前缀；Sex 标准化；RM 独有列加 RM_ 前缀；从 var_names 重建 6 级 taxonomy 列
  10_expand_and_attach_phylo.py    anndata var 从 6,857 扩到 GG2 24.09 全 8,114（补 1,257 全 0 列 + observed bool 列）；挂 varp['taxo_dist'] + varp['phylo_dist']
  11_finalize_corpus.py(+.sbatch)  X→相对丰度 (float32) + 原 counts 备份到 layers['counts']；obs +3 列 (total_reads/n_taxa/Run_paired_id)；var +4 列 (n_samples_observed/prevalence/mass_fraction/mean_rel_abundance_when_present)；uns['provenance']
  12_attach_poincare.py            euclidean_d32.npz → varm['position_encoding'] (8114×32 float32)；uns['provenance']['position_encoding'] 记录 method/dim/scale/pearson；输出 MCFCorpusV2.gg2.h5ad（gzip 压缩）

scripts/Poincare/            # 8,114 GG2 genus → 欧氏 R^D 位置编码（基于 patristic 距离矩阵）
  02_euclidean_embed.py            普通 Adam + L2 距离 + 绝对 MSE loss；d_tree 归一化到 [0,10]；FP32 训练；CPU < 15s；输出 euclidean_d{D}.npz
                                   （早期 Sarkar 双曲嵌入 Pearson 0.77 / geoopt Poincaré ball Pearson 0.98 都已废弃：
                                    下游"MLP → token embedding 相加"管线天然欧氏，双曲位置编码 mismatch 大）

scripts/Phylogeny/           # 生成 GG2 24.09 全 genus 距离矩阵（独立于 anndata，输出可发布）
  01_extract_genus_vocab.py            从 phylogeny.id.nwk 解析全 8,114 个 g__ 节点 + walk-up 拼完整 6 级 path
  02_compute_taxonomic_distance.py     分类法层级 hop 距离 (8114² int8, 取值 {0..6})；向量化 rank-by-rank 比对
  03_compute_phylogenetic_distance.py  patristic 距离 (8114² float32)；手写 newick 解析（只保留内部节点）+ scipy.csgraph.shortest_path 分批 Dijkstra；不依赖 skbio

scripts/SampleDistance/      # 50k 分层抽样 → BC + weighted UniFrac → PCoA + 可视化
  01_stratified_sample.py        MA 30k + RM 20k 两层 sqrt 分层（seed=42）
  02_build_subset_anndata.py     主表切子集（含 varp 搬运）
  03_collapse_tree.py            phylogeny.id.nwk 折叠到 8,114-tip genus_tree.nwk（嵌套 g__ 插 0 枝长 self-tip）；20 对 patristic cross-check
  04_compute_bc.py(+.sbatch)     50k² Bray-Curtis（64 核 ~12 min；fork Pool + scipy cdist）
  05_compute_wunifrac.py(+.sbatch) 50k² Striped Fast UniFrac（64 核 < 1 min）
  06_pcoa.py(+.sbatch)           truncated eigsh top-10；输出 obsm + eigenvalues.tsv
  07_plot_pcoa.py                4 张 PCoA 散点图（bucket / scree / human_sites / rm_sites）

scripts/SampleDistance100k/  # 独立 100k 分支：corpus-级 PCoA 可视化（不改主表，不替代 50k）
  01_stratified_sample_100k.py     10k paired Runs（双库各 1 份）+ 70k MA + 10k RM = 80k MA / 20k RM（seed=42）
                                   MA 加额配额：Human 22k / Animal 12k / Soil 13k / Aquatic 12k / Plant 6k / Unknown 5k
  02_build_subset_anndata_100k.py  主表切 100k 子集（含 varp 搬运）
  03_compute_bc_100k.py(+.sbatch)  100k² Bray-Curtis float16
  04_compute_wunifrac_100k.py(+.sbatch) 100k² Striped Fast UniFrac float16
  05_pcoa_100k.py(+.sbatch)        truncated eigsh top-10 + 负特征值审计
  06_plot_pcoa_100k.py(+.sbatch)   4×2 主图：MA env / MA human sites / RM sites / paired Runs × {BC, wUniFrac}
                                   产物在 results/sample_distance_100k/，与 50k 完全分离

scripts/jxt_scripts/         # 上游 ASV 与 SILVA 注释来源（外部，不改）
```

## 集群提交注意事项（本 HPC 限制）

- **QOS 限制**：单用户排队作业上限约 100。array job 提 398 会被 `QOSMaxSubmitJobPerUserLimit` 拒绝；
  必须分批，本项目用 `01_run_loop.sh` 维持 90 并发
- **本地 sbatch wrapper 强制要求命令行带资源 flag**（即使 `#SBATCH` 已写）：
  ```bash
  sbatch -c 8 --mem=128G -t 04:00:00 script.sbatch     # ✓
  sbatch script.sbatch                                  # ✗ 报"必须包含申请核数、内存、任务运行时间"
  ```
  从 nohup / 非交互 shell 提交时尤其严格。`01_run_loop.sh` 已在脚本内部用
  `SBATCH_RES=( -c 8 --mem=128G -t 04:00:00 )` 兜底。

---

## 关键决策（一句话版）

### 数据库与方向

- **MicrobeAtlas → GG2 NB `auto`**（默认）— 代表序列方向已统一，加 `both` 反而退化（−1.96pp）
- **ResMicroDb → GG2 NB `--p-read-orientation both`** — 方向混合 ASV，`both` 比 `auto` 多救 +6.51pp reads
- 版本要求：QIIME2 ≥ 2025.7 才有 `both`（本项目 `qiime2-2026.1` ✓）

> 完整对照实测、jxt 流水线 7 步、SRP515474 案例 → `.claude/data/direction_orientation.md`

### 训练数据选定

| 角色 | 文件 | shape |
|---|---|---|
| Stage 1 预训练 | `results/feature_table/gg2.full.qc.h5ad` | 1,762,635 × 6,306 |
| Stage 2 二次预训练 | `results/feature_table/resmicrodb.gg2.genus.qc.h5ad` | 93,425 × 4,952 |
| **MiCoFormer 语料库 (含位置编码)** | `results/feature_table/MCFCorpusV2.gg2.h5ad` | 1,826,126 × 8,114 |

- 用 **full** 不用 minfilter（minfilter 的 ≥20 OTUs 阈值与表征学习目标不符）
- 用 **qc** 不用 qc_v2（接受 ~20% reads 损失换语义干净 var 空间）

> full vs minfilter 诊断、qc vs qc_v2 取舍、文件体积差异 → `.claude/data/training_decisions.md`

### 合并立场：语料库越宽越好，下游训练时再筛子集

- 07/08 步只删"不可能进入任何训练分析"的样本（MA: RNAseq 9,415 + Sequencing_Type NaN 20,519，共 29,934 行 / 1.7%）
- 保留：WGS（MA 经流水线提取 16S read 后与 AMPLICON 在 OTU 表语义对齐）/ 阴性对照 / 异常坐标 / 各种宿主体外样本（动物/土壤/水/植物）
- 同 Run 在两库重复 32,698 个 → 合并时两份都保留（不同流水线视角；带 `Database` 标签便于下游再筛）

### 合并 obs schema（54 列）

- **obs_names 全局唯一**：MA 行 `MA_0000000..MA_1732700`，RM 行 `RM_0000000..RM_0093424`。Run 字段作为独立列保留并允许重复
- **公共 9 列不前缀**：`Database / Run / BioSample / Project_ID / Sequencing_Type / Sex / Smoking / Latitude / Longitude`
- **MA 独有 17 列保持 `MA_*` 前缀** / **RM 独有 28 列新加 `RM_*` 前缀**（07/08 步骤产物不带 RM 前缀，前缀仅在 09 合并步引入；保证 RM 单库使用时列名干净）
- 跨库 token 取值集（不映射）：`Sequencing_Type ∈ {AMPLICON, WGS, 16S}` / `Smoking ⊇ {Smoker, Non-smoker, Ex-smoker}` / `BioSample` MA=SRS 系 RM=SAMN 系并存
- 标准化动作：`Sex` MA `female/male` → `Female/Male`（按 RM 大小写）
- 语义重合但分类法不同的列（样本部位 / 健康状态 / 年龄分箱 / 地理）**不映射**，分别保留 `MA_*` 与 `RM_*` 两套

### 数据库选择：GG2 而非 SILVA

- GG2 = GTDB phylogeny，分类粒度更准确
- SILVA 复合属（如 `Escherichia-Shigella`）+ 大菌系过度合并（Bacteroides 6+ clade 打包）会污染样本表征
- GG2 的 family-stall 弱点集中在 Pseudomonas / Haemophilus / Escherichia / Enterobacter（开关式：仅约 14% 样本里它们高丰度）

> 全局加权、ASV 一致性矩阵、bipolar 影响、SILVA 质量问题、用例判断 → `.claude/data/gg2_vs_silva.md`

### Genus 距离矩阵选树和算法

- **用 `phylogeny.id.nwk`**（真实系统发育树，~23M tip，带枝长 + 内部节点带 GG2 rank 标签），不要用 `taxonomy.id.nwk`（只是把分类法层级用 newick 表达，**无枝长**）
- `taxonomy.id.tsv.gz` 不是 GG2 完整 genus 空间——只有 6,757 个 g__（GG2 自己 DEPP 在 V4 ASV tip 级别的保守输出）；**真值是 phylogeny.id.nwk 内部节点的 8,114 个 g__**
- **不要用 skbio `shear` + `tip_tip_distances`**：在 GG2 这种"23M tip + 31.7% 0 枝长内部节点 + 巨型多分叉"结构下，shear 的 unifurcation collapse 丢枝长，50% 非对角 pair 算出 0。改用手写 newick 解析（只内部节点 ~100 万）+ `scipy.csgraph.shortest_path` 分批 Dijkstra（17 批 × 500 sources，~32 min）
- 4 个 -1e-6 浮点噪音负枝长 → clamp 到 0
- 落盘前必须 cross-check：scipy 结果 vs 独立"祖先链遍历"算法在 20 个随机抽样 pair 上必须吻合（diff < 1e-3）；phylo 距离按 taxo hop 分桶必须单调递增

### Taxonomic hop vs phylogenetic distance 的关系

- 内部节点方案的语义：`phylo_dist` 表示 **GG2 定义的 genus clade node 到另一个 genus clade node 的 patristic distance**，不是"某属所有 tip 到另一属所有 tip 的平均距离"。这与最终 genus-level feature table 的粒度一致；若以后要解释 ASV/OTU 层面的属内多样性，则需要另做 tip/代表序列层面的距离。
- `taxo_dist` 是 0-6 的离散层级 hop，解释性强；`phylo_dist` 是连续树距离，可以细分同一个 hop 内部的远近。两者不是重复信息：全 GG2 genus pair 抽样相关性约 Pearson 0.52 / Spearman 0.48；随机比较不同 taxo hop 的 pair 时，phylo 顺序与 taxo 顺序约 82.4% 一致，仍有约 17.6% 会被连续树距离重新排序。
- 全 8,114 genus 上，phylo 距离按 taxo hop 的中位数和 5%-95% 范围：

| taxo hop | 含义 | phylo median | phylo 5%-95% |
|---:|---|---:|---:|
| 1 | 同 family 不同 genus | 33.4 | 6.1 - 119.7 |
| 2 | 同 order 不同 family | 77.9 | 23.8 - 174.5 |
| 3 | 同 class 不同 order | 127.9 | 52.9 - 260.4 |
| 4 | 同 phylum 不同 class | 180.2 | 104.4 - 269.5 |
| 5 | 同 domain 不同 phylum | 258.8 | 148.5 - 407.5 |
| 6 | 跨 domain | 316.3 | 240.7 - 444.4 |

- 关键观察：phylo 均值随 taxo hop 单调增加，说明大方向符合分类层级；但每个 hop 内部离散度很大，说明把所有同层级 pair 压成 1/2/3/4/5/6 会丢掉不少结构。例如 `taxo_dist=1` 的同 family pair，phylo 可从接近 0 到 221.3；`taxo_dist=5` 的跨 phylum pair，phylo 可从 5.7 到 656.3。
- 具体例子：同为 `Flavobacteriaceae` 的 `g__Hyunsoonleella_826360` vs `g__Muricauda_A_821532`，`taxo_dist=1` 但 `phylo_dist=221.345`，远高于很多 hop=2/3 的 pair；跨 phylum 的 `g__SZUA-191` vs `g__Entotheonella`，`taxo_dist=5` 但 `phylo_dist=5.706`，说明分类 rank hop 有时会掩盖树上的近邻关系。
- 结论：保留两套距离有价值。`taxo_dist` 适合作为离散、可解释的分类层级 prior；`phylo_dist` 适合作为连续、可细化同层级内部差异的系统发育 prior。下游建议分别测试 `taxo only` / `phylo only` / `taxo + phylo`。
- 注意：GG2 树存在零长度内部枝，当前 `phylo_dist` 有少量非对角 0 距离（审计为 219 个 pair）。一般不影响使用；若下游算法要求不同 genus 距离严格大于 0，可用 `phylo_dist + 1e-6 * (taxo_dist > 0)` 作 tie-break。

---

## 流水线产物速查

### MicrobeAtlas (sample × genus)
- 03 → `gg2.full.h5ad` 2,690,735 × 7,424（保留 98.89% biological reads）
- 04 → `gg2.full.nonzero.h5ad` 2,380,504 × 7,424（删 11.53% 零 taxon 样本）
- 05 → `gg2.full.qc.h5ad` 1,762,635 × 6,306（**Stage 1 训练用**）
- 06 → `gg2.full.qc.with_meta.h5ad` 1,762,635 × 6,306（X 不变，obs 补 26 列 metadata）
- 07 → `gg2.full.qc.with_meta.filtered.h5ad` 1,732,701 × 6,306（删 RNAseq + NaN 共 29,934 行 / -1.70%）

> 03 步 reads 漏斗、Unmapped 75% 解释、零 taxon 样本来源 → `.claude/data/pipeline_audit_microbeatlas.md`
> 06 obs 26 列 schema、col2/col3 解析、缺失值约定、跨库 merge 策略、anndata StringArray 写盘陷阱 → `.claude/data/obs_metadata_schema.md`

### ResMicroDb (sample × genus)
- 02 → `resmicrodb.gg2.asv.h5ad` 100,391 × 1,639,277（带 GG2 + SILVA 双注释）
- 03 → `resmicrodb.gg2.genus.h5ad` 100,391 × 5,891（98.89% reads）
- 04 → `resmicrodb.gg2.genus.nonzero.h5ad` 100,342 × 5,891
- 05 → `resmicrodb.gg2.genus.qc.h5ad` 93,425 × 4,952（**Stage 2 训练用**）
- 06 → `metadata_all.standardized.{tsv,parquet}` 135,746 × 36（all 全表标准化版）
- 06b → `metadata_all.standardized.fixed.parquet` 135,746 × 36（v2 NCBI 真值方案：1,697 Run 修正 / 5 列）
- 07 → `resmicrodb.gg2.genus.qc.with_meta.h5ad` 93,425 × 4,952（X 不变，obs 补 37 列 = Database + Run + 35 metadata，消费 fixed parquet）
- 08 → `resmicrodb.gg2.genus.qc.with_meta.filtered.h5ad` 93,425 × 4,952（passthrough；anndata 已 100% 16S，row filter 无动作）

> 02-05 实测漏斗、both vs auto 收益、05 损失分解 → `.claude/data/pipeline_audit_resmicrodb.md`
> 06 输入四件套（metadata_all UTF-16LE / sampleTable_changed / sampleTable_v5 / projectTable）、36 列 schema、派生列逻辑（Age_Group jxt 7 档 + Case_Or_Control + Is_Healthy）、Region_16S 填充规则、与 jxt 对照 → `.claude/data/metadata_resmicrodb_standardize.md`
> 06b v2 NCBI 真值方案（4 个 patch / 共 1,697 行）：Patch 1 Sample_Site Nasal→Nasopharynx 1,257（0105 源） / Patch 2 PRJNA801796 Influenza 细分 + Phenotype_ID 补齐 255（0105 源） / Patch 3 PRJNA822681 NCBI host disease 152（76 COVID-19↔Health 翻转 + 76 反向，比 0105 多覆 39 行） / Patch 4 PRJNA824137 NCBI title 33（13 HC + 3 TBZ + 17 LTBI=Latent Tuberculosis Infection/MONDO_0040753/case，14 TBM 已对）→ `rawdata/ResMicroDb/supplement data/CHANGES_0105_new7.md`

### 跨数据集合并

中间步：
- 09 → `merged.gg2.h5ad` **1,826,126 × 6,857**（obs 54 列；reads 守恒 47.68 B；X int32 CSR，density 2.1%；size 849 MB）
- 10 → `merged.gg2.with_phylo.h5ad` **1,826,126 × 8,114**（X 前 6,857 列 nnz/dtype 不变 + 新增 1,257 全 0 列；var 增 `observed` bool 列；`varp['taxo_dist']` + `varp['phylo_dist']`；obs 54 列不变）；size 1.05 GB

最终产物：
- 11 → **`MCFCorpus.gg2.h5ad`** 1,826,126 × 8,114（**X 改成 relative_abundance float32 CSR**；原 int32 counts 备份到 `layers['counts']`；obs 57 列 = 原 54 + `total_reads`(int64) + `n_taxa`(int32) + `Run_paired_id`(string)；var 11 列 = 原 7 + `n_samples_observed` + `prevalence` + `mass_fraction` + `mean_rel_abundance_when_present`；`uns['provenance']` 元数据；varp 保留）；size 1.80 GB
  - **下游用法**：`adata.X` 直接拿到 row_sum=1 的相对丰度；想要原始 reads 用 `adata.layers['counts']`
  - `obs['Run_paired_id']`：跨库 Run 重复对的 group id（32,698 对 × 2 = 65,396 行非空）；下游分 train/val/test 时 group split 防泄漏
- 12 → **`MCFCorpusV2.gg2.h5ad`** 1,826,126 × 8,114（V1 + `varm['position_encoding']` 8114×32 float32）；size 1.80 GB
  - 欧氏位置编码：从 `varp['phylo_dist']` 出发，普通 Adam + L2 距离 + 绝对 MSE loss，dim=32，d_tree 归一化到 [0,10]
  - 验证：Pearson(d_eucl, d_tree) = **0.9924**，MAD = 0.0506，坐标范数 median=2.63 / max=6.41
  - `uns['provenance']['position_encoding']` 记录 method/dim/scale/pearson/mad/npz 路径
  - 下游用法：`adata.varm['position_encoding']` 直接喂 MLP → 与 token embedding 相加；想还原 patristic 距离用 `d_eucl / scale`

合并语义（不随 11 改变）：
- var 来源：共享 4,401 + 仅 MA 1,905 + 仅 RM 551
- Run 重复 65,396 行 = 32,698 个 Run × 2 副本（两库流水线视角并存）
- obs 公共 9 列 + `MA_*` 17 + `RM_*` 28（11 又叠 3 列）
- 06b 已落地：下游 07/08/09 已基于 fixed metadata 重跑，shape 与 reads 不变（仅 obs 取值修正）

### Genus 距离矩阵（GG2 24.09 全量，独立副产物 + 挂回 anndata）

副产物（`results/phylogeny/`，可独立发布）：
- 01 → `genus_vocab.tsv` 8,114 行 × 6 列（GG2 24.09 全 genus + 完整 6 级 rank path）
- 02 → `genus_taxo_dist.npz` 8114² int8，取值 {0..6}；362 KB gzip
- 03 → `genus_phylo_dist.npz` 8114² float32 patristic 距离；209 MB

挂回 anndata：见上节 10 步产物（`varp['taxo_dist']` + `varp['phylo_dist']`），11 步沿用不动。

三个 8,114 vs 6,857 数字的关系：
- **8,114** = GG2 24.09 完整 genus 空间（`phylogeny.id.nwk` 内部节点 g__ 计数 = NB classifier 训练标签空间）
- **6,757** = `taxonomy.id.tsv.gz` 唯一 g__（GG2 自己 DEPP 给每 tip 分类，V4-only 信息不够，保守）
- **6,857** = anndata var.Genus（你的 OTU 全长 16S → NB → 输出，跨 6,757 但只覆盖 8,114 的 84.5%）；anndata 独有但 TSV 缺的 900 个 = NB 比 DEPP 大胆，**这就是 vocab 必须从树取不能从 TSV 取的原因**

### Sample × sample 距离矩阵（50k 分层子集）

`MCFCorpus.gg2.h5ad` 1.83M 全样本上算 pairwise 不现实（1.83M² 矩阵存不下也算不动），改在 50k 分层子集上算。产物全在 `results/sample_distance/`，**主 anndata 不变**。（注：子集是 11 finalize 之前从 10 步 `merged.gg2.with_phylo.h5ad` 切出来的，X 仍是 counts；BC/UniFrac 计算内部自己做相对丰度归一化，与最终 corpus 的 X 语义无冲突。）

抽样：MA 30k（Human 8k / Animal_other 8.1k / Soil 5.2k / Aquatic 5.0k / Plant 2.2k / Unknown 1.5k，sqrt 加权）+ RM 20k（按 RM_Sample_Site × Project_ID）。跨库 Run 重复不去重。`RANDOM_SEED = 42`。

- 02 → `subset_50k.h5ad` 50,000 × 8,114（X + obs(56 列：原 54 + stratum_id + sub_stratum) + varp 搬主表的 taxo/phylo_dist）
- 03 → `genus_tree.nwk` 8,114-tip folded tree；嵌套 g__ 插 0 枝长 self-tip（GG2 有 1,357 个 g__ 是另一个 g__ 的 GTDB 子属后代）；新树 patristic 与 `genus_phylo_dist.npz` 20 对 cross-check < 1e-3 通过
- 04 → `obsp['distance_bc']` 50k² **float16**；relative abundance 后 scipy cdist
- 05 → `obsp['distance_wunifrac']` 50k² **float16**；Striped Fast UniFrac (McDonald 2018 / Sfiligoi 2022)；BIOM v2.1 临时文件 + tree 双输入；64 核 < 1 min
- 06 → `obsm['X_pcoa_bc'/'X_pcoa_wunifrac']` 50k × 10 + `pcoa_eigenvalues.tsv`；truncated eigsh 上 dense float32 B；trace(B) 作分母算 explained var ratio
- 07 → `figures/pcoa_{bucket,scree,human_sites,rm_sites}.png` 4 张可视化

PCoA top-3 累计 explained variance：BC 20.8%（散得开）/ wUniFrac 47.3%（phylo 把方差压到少数轴）。

不做 rarefaction（McMurdie & Holmes 2014）；relative abundance 即可。

最终文件 `subset_50k.h5ad` 6.3 GB（X 340 MB + 两个 50k² float16 距离 ~5 GB + varp 310 MB + obs/obsm 30 MB）。

> 抽样配额 sqrt 推导、tree 折叠算法细节、嵌套 g__ self-tip 推导、PCoA 各部位偏移解释 → `.claude/data/sample_distance.md`

### Sample × sample 距离矩阵（100k 独立分支）

`scripts/SampleDistance100k/` 是 corpus-级 PCoA 可视化的独立分支，产物落在 `results/sample_distance_100k/`，**与 50k 完全分离，不替代也不修改 50k**。设计目标是给出大样本量下的 PCoA 总览图，覆盖更多 paired Runs 与各 MA bucket。

抽样设计 100,000 行：
- **10,000 cross-database paired Runs**：在两库 Run 交集里按 RM_Sample_Site × Project_ID 双层 sqrt 分层选 10k 个 Run，每个 Run 同时入 MA 和 RM 两份（共 20k 行）—— 用于跨库 paired 可视化
- **70,000 MA 加额**：按 6 个 bucket 配额 Human 22k / Animal_other 12k / Soil 13k / Aquatic 12k / Plant 6k / Unknown 5k，每个 bucket 内部按 body site / env subcategory 做 sqrt 分层
- **10,000 RM 加额**：按 RM_Sample_Site × Project_ID 双层 sqrt 分层
- 最终库分布 80k MA + 20k RM，`RANDOM_SEED=42`

产物：与 50k 同结构（subset h5ad + genus_tree.nwk + BC/wUniFrac float16 距离 + PCoA top-10 + 负特征值审计 + figures），主图是 `figures/pcoa_100k_4x2.png` —— 4 行（MA env / MA human sites / RM sites / paired Runs）× 2 列（BC / wUniFrac），每面板都画 100k 灰色背景，仅 highlight 层切换，坐标在同度量内可直接比较。

子分支详细 README 在 `scripts/SampleDistance100k/README.md`。
