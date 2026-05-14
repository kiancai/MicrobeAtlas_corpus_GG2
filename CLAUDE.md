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

scripts/Merge/               # 两库合并 + 挂距离矩阵
  09_merge_databases.py            outer-join var；obs_names 重编号 MA_/RM_ 前缀；Sex 标准化；RM 独有列加 RM_ 前缀；从 var_names 重建 6 级 taxonomy 列
  10_expand_and_attach_phylo.py    anndata var 从 6,857 扩到 GG2 24.09 全 8,114（补 1,257 全 0 列 + observed bool 列）；挂 varp['taxo_dist'] + varp['phylo_dist']

scripts/Phylogeny/           # 生成 GG2 24.09 全 genus 距离矩阵（独立于 anndata，输出可发布）
  01_extract_genus_vocab.py            从 phylogeny.id.nwk 解析全 8,114 个 g__ 节点 + walk-up 拼完整 6 级 path
  02_compute_taxonomic_distance.py     分类法层级 hop 距离 (8114² int8, 取值 {0..6})；向量化 rank-by-rank 比对
  03_compute_phylogenetic_distance.py  patristic 距离 (8114² float32)；手写 newick 解析（只保留内部节点）+ scipy.csgraph.shortest_path 分批 Dijkstra；不依赖 skbio

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
| **合并语料库** | `results/feature_table/merged.gg2.h5ad` | 1,826,126 × 6,857 |

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

- 09 → `merged.gg2.h5ad` **1,826,126 × 6,857**（obs 54 列；reads 守恒 47.68 B；X int32 CSR，density 2.1%；size 849 MB）
- var 来源：共享 4,401 + 仅 MA 1,905 + 仅 RM 551
- Run 重复 65,396 行 = 32,698 个 Run × 2 副本（两库流水线视角并存）
- obs 54 列 = 公共 9 + `MA_*` 17 + `RM_*` 28
- 06b 已落地：下游 07/08/09 已基于 fixed metadata 重跑，shape 与 reads 不变（仅 obs 取值修正）

### Genus 距离矩阵（GG2 24.09 全量，独立副产物 + 挂回 anndata）

副产物（`results/phylogeny/`，可独立发布）：
- 01 → `genus_vocab.tsv` 8,114 行 × 6 列（GG2 24.09 全 genus + 完整 6 级 rank path）
- 02 → `genus_taxo_dist.npz` 8114² int8，取值 {0..6}；362 KB gzip
- 03 → `genus_phylo_dist.npz` 8114² float32 patristic 距离；209 MB

挂回合并 anndata：
- 10 → `merged.gg2.with_phylo.h5ad` **1,826,126 × 8,114**（X 前 6,857 列 nnz/dtype 不变 + 新增 1,257 全 0 列；var 增 `observed` bool 列：原 6,857 True / 新 1,257 False；`varp['taxo_dist']` + `varp['phylo_dist']`；obs 54 列不变）；size 1.05 GB

三个 8,114 vs 6,857 数字的关系：
- **8,114** = GG2 24.09 完整 genus 空间（`phylogeny.id.nwk` 内部节点 g__ 计数 = NB classifier 训练标签空间）
- **6,757** = `taxonomy.id.tsv.gz` 唯一 g__（GG2 自己 DEPP 给每 tip 分类，V4-only 信息不够，保守）
- **6,857** = anndata var.Genus（你的 OTU 全长 16S → NB → 输出，跨 6,757 但只覆盖 8,114 的 84.5%）；anndata 独有但 TSV 缺的 900 个 = NB 比 DEPP 大胆，**这就是 vocab 必须从树取不能从 TSV 取的原因**
