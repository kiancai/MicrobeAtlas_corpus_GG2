# GreenGene2 / MicrobeAtlas 项目

> 本文件是**总入口/导航中枢**，只保留概述与编排级信息。脚本/产物/决策/计划/探索的细节都在 `.claude/` 子文档（见下方「文档架构与导航」）。

## 项目概述

构建 **MiCoFormer 微生物组基础模型语料库**：sample × genus 相对丰度矩阵，词表 = GG2 24.09 全 8,114 genus。**五库已合并并 finalize（merge-then-filter，全留不删 + QC 旗标）**：
- **MicrobeAtlas**（~269 万 16S）· **ResMicroDb**（10 万 ASV）· **MicroBioMap**（16.8 万）· **PRIME**（5.3 万，GG2 24.09 现成）· **Qiita**（32.2 万）
- **最终语料 `MCFCorpusV3.gg2.h5ad`（3,334,541 × 8,114）**：X=相对丰度·layers['counts']·obs 185 列(QC旗标+5库 metadata+GMrepo/mBodyMap 富集)·varp(taxo/phylo 距离)·varm(Poincaré 位置编码)。已回读核验全过。

**当前阶段**：语料 V3 已建成、验证、清洗、并出**自包含数据字典中英双版**（`results/feature_table/DATA_DICTIONARY_MCFCorpusV3{.md,.zh.md}`，HF 出版级；含数据集构成/QC/重复说明，不引用内部脚本）。筛选/去重交下游（用户定不做训练视图）。进度见 `.claude/roadmap.md`。
> 合并架构(merge-then-filter / QC=打标不删 / Run_paired_id) → `.claude/merge_architecture.md`；全局审计 → `.claude/data/corpus_sample_overlap.md`。

## 环境配置

- **Conda 环境**：`baseBio`（路径 `/hpcdisk1/limk_group/caiqy/miniforge3/envs/baseBio`）；R/phyloseq 用 `baseR`
- 已装：`biopython`, `biom-format`, `h5py`, `pandas`, `cutadapt`, `fastp`, `numpy`, `pigz`, `pbzip2`, `scikit-bio`（仅 01 用）, `scipy`
- 代理：`proxy_on` 激活（`127.0.0.1:60124`），conda 安装前先激活

## 数据目录结构

```
rawdata/
├── MicrobeAtlas/
│   ├── sample_info/samples.env.info.tsv          # 样本元信息 (301 MB)
│   ├── OTU_info/{otus.info.tsv, mapref-3.0.tar.gz}
│   └── OTU_count/{otus.97.allinfo, samples-otus.97.mapped[.minfilter].biom.gz}
├── ResMicroDb/
│   ├── 16S/<PROJECT>/results/                     # 398 项目；asv.fa / otutab.txt / *.qza
│   │       ├── taxonomy_silva.txt                 # SILVA 注释（jxt）
│   │       └── taxonomy_gg2.{qza,txt}             # GG2 注释（本项目新增）
│   └── raw_reads/PRJEB13657/                      # 仅 1 个项目齐全
└── (扩增候选，待整合)
    ├── PRIME/         # gg_genus_{absolute,relative}.csv + samples.csv（GG2 24.09 直接合并）
    ├── GMrepo/        # gmrepo.sql.gz + *.txt.gz（肠道，仅元数据，无序列）
    ├── mBodyMap/      # mbodymap_{samples,projects}.zip（全身 22 部位，仅元数据）
    └── MicroBioMap/   # 人类微生物组纲要（168k 人肠道，Zenodo 有 ASV 序列）
```

> 各文件列含义 / BIOM 读取 / mapref 结构 → `.claude/data/file_schemas.md`

## 数据规模摘要

| 指标 | 数量 | | 指标 | 数量 |
|------|------|---|------|------|
| MA 总样本 | ~269 万 | | RM 样本 / ASV | ~10 万 / ~164 万 |
| MA minfilter 后 | ~188 万 | | GG2 24.09 全 genus（词表） | 8,114 |
| MA 97% OTU | ~10.3 万 | | **最终语料行 × 列(V3, 五库)** | **3,334,541 × 8,114** |

## 集群提交注意事项（本 HPC 限制）

- **QOS**：单用户排队上限约 100。array job 提 398 会被 `QOSMaxSubmitJobPerUserLimit` 拒；必须分批，用 `01_run_loop.sh` 维持 90 并发。
- **本地 sbatch wrapper 强制命令行带资源 flag**（即使 `#SBATCH` 已写）：
  ```bash
  sbatch -c 8 --mem=128G -t 04:00:00 script.sbatch     # ✓
  sbatch script.sbatch                                  # ✗ 报"必须包含申请核数、内存、任务运行时间"
  ```
  nohup / 非交互 shell 提交尤其严格。`01_run_loop.sh` 内部用 `SBATCH_RES=(...)` 兜底。

## 核心决策速览

> 一句话版；完整理由与实测见 `.claude/decisions.md`。

- **GG2 不用 SILVA**：GG2=GTDB phylogeny，粒度更准；SILVA 复合属/大菌系过度合并会污染表征。
- **分类方向**：MA → NB `auto`（方向已统一，`both` 退化 −1.96pp）；RM → NB `--p-read-orientation both`（混合方向，`both` 多救 +6.51pp）。新库重注释统一走 `both`。
- **训练数据**：用 **full** 不用 minfilter；用 **qc** 不用 qc_v2。关键文件：

  > 注：下表 Stage1/2 + V2 是**旧 filter-then-merge** 产物（已被 V3 取代，旧码归档 `scripts/archive/Merge/`）。
  > **现行语料 = `MCFCorpusV3.gg2.h5ad`（3,334,541 × 8,114，五库 merge-then-filter，QC=旗标不删）**。

  | 角色（旧） | 文件 | shape |
  |---|---|---|
  | Stage 1 | `results/feature_table/gg2.full.qc.h5ad` | 1,762,635 × 6,306 |
  | Stage 2 | `results/feature_table/resmicrodb.gg2.genus.qc.h5ad` | 93,425 × 4,952 |
  | 旧语料 V2 | `results/feature_table/MCFCorpusV2.gg2.h5ad` | 1,826,126 × 8,114 |

- **合并立场（V3 已落地）**：语料越宽越好，**全留不删**（QC 降级为 obs 旗标 `qc_pass`）；同 Run 跨库重复全留（`Database` 标签 + `Run_paired_id` 串联，下游再去重）。详 `merge_architecture.md`。
- **距离矩阵**：从 `phylogeny.id.nwk` **内部节点**取 8,114 genus（非 taxonomy.id.tsv 的 6,757）；保留 `taxo_dist`(0-6 hop) + `phylo_dist`(patristic) 两套。

## 📁 文档架构与导航

```
CLAUDE.md                  # 本文件：概述 + 导航 + 维护约定（精简）
.claude/
├── pipeline.md            # 流水线脚本一览 + 各步产物速查
├── decisions.md           # 关键决策详解（方向/训练数据/合并/GG2-vs-SILVA/距离）
├── roadmap.md             # 准备做 / 进行中 / 已完成
├── merge_architecture.md  # 合并重构 v2(merge-then-filter)：组成表契约 + Layer2 + 决策
├── explorations/          # 每次探索的带日期记录（README.md = 时间线索引）
└── data/                  # 稳定参考（详细信息）
    ├── file_schemas.md                    # 文件列含义 / BIOM 读取 / mapref
    ├── sample_env_stats.md                # Environments / body site 计数
    ├── otu_hierarchy.md                   # OTU 嵌套层级 90_x;96_x;97_x
    ├── direction_orientation.md           # auto vs both 实测 / jxt 7 步
    ├── training_decisions.md              # full vs minfilter / qc vs qc_v2
    ├── gg2_vs_silva.md                    # GG2 vs SILVA 全量论证
    ├── pipeline_audit_microbeatlas.md     # MA reads 漏斗 + 数据快照时间(≈2021-01)
    ├── obs_metadata_schema.md             # MA obs 26 列 schema
    ├── pipeline_audit_resmicrodb.md       # RM 漏斗
    ├── metadata_resmicrodb_standardize.md # RM metadata 36 列标准化
    ├── metadata_merge_design.md            # 多库 metadata 列级"共同含义"原则 + 逐库逐列清单
    ├── sample_distance.md                 # 50k/100k 抽样 + 树折叠 + PCoA
    ├── corpus_sample_overlap.md           # 全局样本审计：五库规模/两两重叠/净≈300万/∩MGnify 88.7%/GMrepo·mBodyMap 覆盖
    ├── metadata_sources_gmrepo_mbodymap.md# GMrepo/mBodyMap 元数据剖析 + 即用即合
    ├── mgnify_evaluation.md               # MGnify 覆盖率基准(不并入语料) + 全库 run 列表 + 重合实测 + API 配方
    └── qiita.md                           # Qiita 结构/规模/访问/ENA桥/组成表计划(扩样本行第二主力)
```

## 🔧 文档维护约定

当用户说"**落实有效信息到文档中**"时，按下表归位（保持 CLAUDE.md 精简，细节进 `.claude/`）：

| 信息类型 | 落到哪 |
|---|---|
| 一次探索/调研的过程与结论 | 新建 `.claude/explorations/YYYY-MM-DD_<topic>.md`（触发/简介/核心结论/详细信息/产物指针），并在 `explorations/README.md` 加一行时间线 |
| 计划/任务状态变动 | `.claude/roadmap.md`（准备做 ↔ 进行中 ↔ 已完成 三态迁移） |
| 稳定可复用的事实 / schema / 审计 | `.claude/data/<topic>.md`（新建或更新对应参考） |
| 跨领域的关键决策 | `.claude/decisions.md` |
| 流水线脚本 / 产物变化 | `.claude/pipeline.md` |
| 概述 / 目录结构 / 规模摘要 变化 | 才改 CLAUDE.md 本体 |

原则：**CLAUDE.md 只保留导航 + 编排级信息；任何细节一律放 `.claude/` 子文档，避免 CLAUDE.md 变长。**
