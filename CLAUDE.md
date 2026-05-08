# GreenGene2 / MicrobeAtlas 项目

## 项目概述

从 MicrobeAtlas 数据库下载的微生物组数据，用于后续分析。两部分：
- **MicrobeAtlas**: 样本-OTU 丰度矩阵、OTU 注释信息、参考序列库（~269 万样本）
- **ResMicroDb**: 抗性微生物数据库的原始测序 reads（398 项目 / ~10 万样本 / ~164 万 ASV）

## 环境配置

- **Conda 环境**: `baseBio`（路径: `/hpcdisk1/limk_group/caiqy/miniforge3/envs/baseBio`）
- 已安装: `biopython`, `biom-format`, `h5py`, `pandas`, `cutadapt`, `fastp`, `numpy`, `pigz`, `pbzip2`
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

scripts/ResMicroDb/          # 398 项目 ASV → 跨项目 sample × genus
  01_qiime2_classify.sbatch        单项目 GG2 NB 分类（已加 --p-read-orientation both）
  01_run_loop.sh                   轮询提交器（CAP=90, INTERVAL=30s, 断点续投）
  02_merge_to_asv_anndata.py       398 项目 → sample × ASV CSR（GG2 + SILVA 双注释）
  03_build_genus_anndata.py        过滤 mito/chloro/非 BA + 6 级路径聚合
  04_drop_empty_samples.py         剔除零 taxon 样本（含项目维度报告）
  05_qc_filter.py                  同 MicrobeAtlas（MIN_FEATURES=5）

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

- 用 **full** 不用 minfilter（minfilter 的 ≥20 OTUs 阈值与表征学习目标不符）
- 用 **qc** 不用 qc_v2（接受 ~20% reads 损失换语义干净 var 空间）

> full vs minfilter 诊断、qc vs qc_v2 取舍、文件体积差异 → `.claude/data/training_decisions.md`

### 数据库选择：GG2 而非 SILVA

- GG2 = GTDB phylogeny，分类粒度更准确
- SILVA 复合属（如 `Escherichia-Shigella`）+ 大菌系过度合并（Bacteroides 6+ clade 打包）会污染样本表征
- GG2 的 family-stall 弱点集中在 Pseudomonas / Haemophilus / Escherichia / Enterobacter（开关式：仅约 14% 样本里它们高丰度）

> 全局加权、ASV 一致性矩阵、bipolar 影响、SILVA 质量问题、用例判断 → `.claude/data/gg2_vs_silva.md`

---

## 流水线产物速查

### MicrobeAtlas (sample × genus)
- 03 → `gg2.full.h5ad` 2,690,735 × 7,424（保留 98.89% biological reads）
- 04 → `gg2.full.nonzero.h5ad` 2,380,504 × 7,424（删 11.53% 零 taxon 样本）
- 05 → `gg2.full.qc.h5ad` 1,762,635 × 6,306（**Stage 1 训练用**）

> 03 步 reads 漏斗、Unmapped 75% 解释、零 taxon 样本来源 → `.claude/data/pipeline_audit_microbeatlas.md`

### ResMicroDb (sample × genus)
- 02 → `resmicrodb.gg2.asv.h5ad` 100,391 × 1,639,277（带 GG2 + SILVA 双注释）
- 03 → `resmicrodb.gg2.genus.h5ad` 100,391 × 5,891（98.89% reads）
- 04 → `resmicrodb.gg2.genus.nonzero.h5ad` 100,342 × 5,891
- 05 → `resmicrodb.gg2.genus.qc.h5ad` 93,425 × 4,952（**Stage 2 训练用**）

> 02-05 实测漏斗、both vs auto 收益、05 损失分解 → `.claude/data/pipeline_audit_resmicrodb.md`

### 跨数据集合并

- MicrobeAtlas 7,424 var ∩ ResMicroDb 4,952 var = **4,446 共享**（89.8% of ResMicroDb）
- outer join 后 7,930 var / ~280 万样本 — 基础模型训练底盘
