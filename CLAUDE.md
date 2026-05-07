# GreenGene2 / MicrobeAtlas 项目

## 项目概述

从 MicrobeAtlas 数据库下载的微生物组数据，用于后续分析。主要包含两部分：
- **MicrobeAtlas**: 样本-OTU 丰度矩阵、OTU 注释信息、参考序列库
- **ResMicroDb**: 抗性微生物数据库的原始测序 reads

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

## 各文件详细说明

### 1. `sample_info/samples.env.info.tsv`

**样本元信息表**，约 **2,690,735 行**，9 列（tab 分隔）：

| 列号 | 列名 | 含义 |
|------|------|------|
| 1 | MAP_SID | 样本 ID（如 ERR/SRR 编号） |
| 2 | Environments | 受控词汇环境标签（见下方详述） |
| 3 | _ | 针对 animal;human 样本的身体部位细分（gut/skin/oral/urogenital 等） |
| 4 | Technology | 测序技术（AMPLICON / WGS / RNAseq） |
| 5 | Keywords | 自由文本关键词标签（与 Environments 互补，非子集关系） |
| 6 | _ | （未命名，通常为空） |
| 7 | Project | 原始项目编号（ERP/SRP） |
| 8 | Institution | 机构 + 地理坐标 |

样本分布：Amplicon ~237 万，WGS ~26 万，RNAseq ~1.2 万。

#### Environments 列要点

- **格式**：`主类;子类`（严格两级），多个标签用 `|` 分隔；约 21% 样本无标签
- **Environments vs Keywords**：独立互补，非子集关系。Environments 是受控词汇，Keywords 是自由文本
- **分类体系**：4 主类（animal/soil/aquatic/plant），共 69 个唯一标签
  - `animal`：37 子类，以 human / mouse / pig / cattle 为主
  - `soil`：9 子类，以 forest / field / agricultural 为主
  - `aquatic`：12 子类，以 marine / sediment / sea / river 为主
  - `plant`：7 子类，以 rhizosphere / leaf / wood 为主
- 详细计数 → `.claude/data/sample_env_stats.md`

#### 第3列（body site）要点

- **非 human 专属**：所有 animal 样本共用同一套词汇（`;` 分隔，最多 3-4 级，可带 healthy/disease/infant 等修饰）
- Human 样本（473,574 个）中 93% 有标注，共 172 种 body site，主要为 gut / oral / skin / urogenital / lung
- 详细分布 → `.claude/data/sample_env_stats.md`

### 2. `OTU_info/otus.info.tsv`

**OTU 注释信息表**，跨所有相似度层级，19 列：

| 列号 | 列名 | 含义 |
|------|------|------|
| 1 | OTU | 层级 OTU ID（见下方说明） |
| 2 | Tax | 分类域（Archaea/Bacteria 等） |
| 3 | SpeciesRep | 物种代表序列 |
| 4 | SeqCount | 该 OTU 的序列数量 |
| 5 | GoldCount | Gold 标准序列数 |
| 6 | GenomeCount | 关联基因组数 |
| 7 | TypeStrains | 模式菌株 |
| 8 | Strains | 菌株 |
| 9 | Genomes | 基因组 |
| 10 | GoldSeqs | Gold 标准序列 |
| 11 | Aliases | 别名 |
| 12-14 | GoldHit/GoldID/GoldScore | Gold 比对结果 |
| 15 | RepSpecies | 代表物种 |
| 16 | Taxaname | 分类名（GTDB 格式） |
| 17 | OrigTax | 原始分类（NCBI 格式） |
| 18 | RepSequenceID | 代表序列 GenBank ID |
| 19 | RepSequence | 代表全长 16S 序列 |

### 3. `OTU_info/mapref-3.0.tar.gz`

**MicrobeAtlas 比对参考数据库**，解压后含：

| 文件 | 含义 |
|------|------|
| `mapref-3.0.fna` | 参考全长序列 FASTA（**1,360,792 条**） |
| `mapref-3.0.fna.otutax` | 每条参考序列的层级 OTU 映射 |
| `mapref-3.0.fna.mscluster` | 序列聚类信息 |
| `mapref-3.0.fna.ncbitax` | NCBI 分类注释 |
| `otus.info` | 同 otus.info.tsv 的副本 |

`mapref-3.0.fna.otutax` 格式示例：
```
KC471280:1..1464    90_17776;96_71281;97_92606;98_125911;99_193128
```

### 4. `OTU_count/otus.97.allinfo`

**97% OTU 详细信息**，每行一个 97% OTU，20 列（tab 分隔，无表头）。列含义通过与 `otus.info.tsv` 交叉验证推断（无官方文档）：

| 列号 | 含义 | 置信度 |
|------|------|--------|
| 1 | OTU ID（格式: `MAPv3;90_x;96_x;97_x`） | 确定 |
| 2 | 某种计数（含义不明，与 SeqCount 不完全一致） | 不确定 |
| 3 | 模式菌株（TypeStrains）基因组登录号 | 较确定 |
| 4 | 菌株列表（Strains），格式 `物种名(登录号)` | 较确定 |
| 5 | 基因组列表（Genomes），比 col4 更宽泛 | 较确定 |
| 6 | 别名/同义词（Aliases） | 较确定 |
| 7 | **代表全长 16S 序列**（已与 otus.info.tsv 交叉验证一致） | 确定 |
| 8 | 代表序列 GenBank 登录号（如 `KF075129:1..1352`） | 确定 |
| 9 | GTDB 完整分类（`d__Bacteria;p__...` 格式） | 确定 |
| 10 | 最佳 GTDB 参考序列登录号 | 推断 |
| 11 | 代表序列对应基因组登录号（如 `RS_GCF_xxx~NZ_xxx`） | 较确定 |
| 13 | 与 GTDB 参考的比对相似度（0-1） | 推断 |
| 14 | 与 GTDB 参考的比对长度（bp） | 推断 |
| 15 | NCBI/SILVA 分类（`Bacteria;Firmicutes;...` 格式） | 确定 |
| 16-17 | 最佳 NCBI/Gold 参考序列登录号 | 推断 |
| 19 | 与 NCBI/Gold 参考的比对相似度（0-1） | 推断 |
| 20 | 与 NCBI/Gold 参考的比对长度（bp） | 推断 |

**提取代表序列（已验证）：**
```bash
awk -F'\t' '{printf ">%s\n%s\n", $1, $7}' OTU_count/otus.97.allinfo > otus97_rep.fasta
```

### 5. BIOM 丰度矩阵文件

两个 HDF5 格式的 BIOM 文件：

| 文件 | 样本数 | OTU 数 | 说明 |
|------|--------|--------|------|
| `samples-otus.97.mapped.biom.gz` | **2,690,735** | 103,271 | 全部样本 |
| `metag.minfilter.refilt.biom.gz` | **1,884,129** | 102,824 | 质量过滤后 |

读取方式（需要先解压，直接读 gzip 会 OOM）：
```python
import h5py
with h5py.File('path.biom', 'r') as f:
    sample_ids = f['sample']['ids'][:]   # 样本 ID 列表
    obs_ids = f['observation']['ids'][:]  # OTU ID 列表
    data = f['observation']['matrix']     # 稀疏矩阵（慎用，内存巨大）
```

## OTU 层级结构说明

MicrobeAtlas 的 OTU 采用**嵌套层级聚类**，ID 格式为：

```
90_17776;96_71281;97_92606;98_125911;99_193128
│       │         │         │         │
90%聚类  96%聚类    97%聚类    98%聚类    99%聚类
```

- 从 90%（宽泛）到 99%（精细），逐级嵌套
- 一个 90% OTU 可包含数千个 97% OTU（如 `90_3` 含 3,681 个 97% OTU）
- 每个 97% OTU 对应**恰好一个代表全长 16S 序列**


## 数据规模摘要

| 指标 | 数量 |
|------|------|
| 总样本数 | ~269 万 |
| 过滤后样本数 | ~188 万 |
| 97% OTU 数 | ~10.3 万 |
| 参考全长序列数 | ~136 万 |
| 主要环境类型 | 人体、动物肠道、土壤、水体 |
| 主要测序技术 | Amplicon (88%)、WGS (10%) |


## 03 步 OTU→Genus 聚合产物校验（2026-05-06）

### 输出
- `results/feature_table/gg2.full.h5ad` — sample × var = `2,690,735 × 7,424`，873 MB
- `results/feature_table/gg2.minfilter.h5ad` — `1,884,129 × 7,424`，749 MB
- 存储格式：`X` 是 `scipy.sparse.csr_matrix`（CSR，只存非零），h5ad 再 gzip 压缩

### OTU 数对账
- `taxonomy.tsv`：111,870（来源 = `otus.97.allinfo` 第 7 列代表序列经 GG2 NB 注释）
- 03 过滤掉 1,835 条：83 非 B/A + 901 mitochondria + 851 chloroplast
- 过滤后 110,035 条 → 与 BIOM 的 103,271 OTU 取交集 = **101,507**
- 8,528 条 tax 不在 BIOM：参考库有但 269 万样本中从未检出，正常
- 1,764 条 BIOM 不在 tax：基本就是被 03 过滤掉的 mito/chloro/non-BA
- 聚合 key（QIIME2 风格 6 级路径 `d__;p__;c__;o__;f__;g__`）唯一值 = **7,424**

### 75% reads 丢弃来源（重要：看着吓人但实际正常）

03 报告 `丢弃 reads: 75.55%`，分解后：

| 类别 | reads 占比 | OTU 数 |
|---|---|---|
| **`Unmapped` 伪 OTU** | **75.15%** | 1 |
| chloroplast | 0.28% | 836 |
| mitochondria | 0.12% | 847 |
| 非 B/A | 0.00% | 80 |

`Unmapped` 是 MicrobeAtlas 在 BIOM 里塞的"未比对到任何参考序列的 reads 计数桶"——单独一行，172.5 G reads。它本来就不是生物特征，不在 `taxonomy.tsv` 里，03 自然把它剔除了。**真正的生物信号丢失只有 ~0.4%**。如果将来要对外报告"过滤损失"，建议把 `Unmapped` 单独列出，避免给读者"75% 生物信号被丢了"的误解。

per-sample 丢失中位数 76% 也是同一原因——MicrobeAtlas 数据本身就有大量未比对 reads，与 03 过滤无关。

### 零 taxon 样本（→ 04 步处理）

03 输出后**未做样本级过滤**，需要后续脚本剔除：

| 数据集 | 样本总数 | 零 taxon | 占比 |
|---|---|---|---|
| `gg2.full.h5ad` | 2,690,735 | **310,231** | 11.53% |
| `gg2.minfilter.h5ad` | 1,884,129 | 0 | 0.00% |

零 taxon = 该样本所有 reads 都落在 `Unmapped` 或被 03 剔除的条目上。`gg2.minfilter` 因 MicrobeAtlas 自己已经做过最小过滤（每样本至少 13 reads / 1 var），不需要这一步。

### 检查代码位置
后续若要复现这套对账，思路：
1. 把 `taxonomy.tsv` 按 mito/chloro/非 B/A 分组，用 OTU id 在 BIOM 中查 `observation/matrix` 各行 reads 总和
2. `BIOM_OTU - tax_OTU` 差集第一个就是 `Unmapped`
3. AnnData 行 nnz 用 `np.diff(X.indptr)`，行和用 `cumsum(data)[indptr]` 求差（避免对每行循环）

---

## ResMicroDb 数据集（2026-05-07 起整理）

### 数据规模
- **398 个 study，100,789 样本**，跨项目 ~164 万 ASV（项目内 ID `ASV_N`，
  跨项目 namespace 化为 `<PROJECT>__ASV_N` 后唯一）
- 项目内 ASV 数中位 1,955，最大 173,864（SRP056779），最小 22（SRP397402）
- 总 reads 约 66 亿

### jxt 流水线（ASV 与 SILVA 注释的来源）

完整代码在 `scripts/jxt_scripts/`，摘要 7 步：

1. 双端 merge（`vsearch fastq_mergepairs`）
2. 切引物 + 质控（`vsearch fastx_filter --fastq_stripleft/right 23 --fastq_maxee_rate 0.01`）
3. **方向矫正**（`vsearch --orient -db silva_16s_v123.udb`）— **关键步骤，问题源头**：
   matched 序列翻为正向，**notmatched 序列保留原方向**
4. 去冗余（`vsearch derep --minuniquesize 10`）
5. 去噪（`usearch unoise3`）→ `asv.fa`
6. 生成特征表（`vsearch usearch_global --id 0.97`）→ `otutab.txt`
7. SILVA 分类（QIIME2 2024.2 `classify-sklearn` + `silva-138-99-nb-classifier.qza`）→ `taxonomy_silva.txt`

**重要认知**：`otutab.txt` 是 vsearch_global @0.97 的 reads 归属，**不是精确 ASV 计数**。
对 genus 级聚合分析无影响，但若做 ASV 级精细分析需注意语义。

### GG2 注释的方向问题（本次发现）

**问题**：jxt 第 3 步用 SILVA v123 做 orient。某些项目（典型如 SRP515474）的引物 / 区段
v123 完全不识别 → 全部走 notmatched 分支 → ASV 集合方向混合。
GG2 NB classifier 默认 `--p-read-orientation auto` 只对前 100 条 ASV 选一向锁定全 batch，
混合方向项目里另一半 ASV 报 Unassigned。

**实测对照**（SRP515474, 73,749 ASV, 5.5 亿 reads）：

| 模式 | reads 到 genus | reads Unassigned | 平均 Confidence |
|---|---:|---:|---:|
| `auto` | 7.10% | 87.00% | 0.78 |
| `reverse-complement` | 7.10% (= auto) | 87.00% | 0.78 |
| **`both`** | **84.41%** | **1.09%** | **0.94** |

**解决方案**：分类 sbatch 加 `--p-read-orientation both`：每条 ASV 独立尝试正反两方向，
取置信度高的一次。计算时间约翻倍；对原本方向已齐的项目无影响。

**版本要求**：QIIME2 ≥ 2025.7（commit `be2c2df` 引入 `both`）。本项目用 `qiime2-2026.1` ✓。
jxt 当年用 2024.2 没这选项，所以 SRP515474 类项目 GG2 注释失败但 SILVA 正常——
推测 SILVA NB 训练时做了双向数据增强，对方向不敏感。

**MicrobeAtlas 不受影响**：代表序列是策展库，方向已统一。auto 模式 Unassigned 仅 0.07%。
本项目重跑 02 时仍统一加 `both` 以保证两数据集参数一致。

### 跨项目合并的不变量

- 样本 ID（ERR / SRR / DRR）跨项目唯一（已在 02 校验）
- ASV ID 跨项目重名 → namespace 为 `<PROJECT>__ASV_N`
- 不同项目的 V 区段不同（V1-V3 / V3-V4 / V4 / V1-V9 都有），跨项目 ASV 序列不可直接比较；
  在 GG2 6 级 genus 路径聚合后大致可比

### 流水线脚本

```
scripts/MicrobeAtlas/        # OTU → genus 聚合（输入是策展的全长代表序列）
  01_extract_rep_seqs.sh           从 otus.97.allinfo 提取 97% 代表序列
  02_qiime2_classify.sbatch        GG2 NB 分类（已加 --p-read-orientation both）
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
```

### 集群提交注意事项（本 HPC 限制）

- **QOS 限制**：单用户排队作业上限约 100。array job 提 398 会被 `QOSMaxSubmitJobPerUserLimit` 拒绝；
  必须分批，本项目用 `01_run_loop.sh` 维持 90 并发
- **本地 sbatch wrapper 强制要求命令行带资源 flag**（即使 `#SBATCH` 已写）：
  ```bash
  sbatch -c 8 --mem=128G -t 04:00:00 script.sbatch     # ✓
  sbatch script.sbatch                                  # ✗ 报"必须包含申请核数、内存、任务运行时间"
  ```
  从 nohup / 非交互 shell 提交时尤其严格。`01_run_loop.sh` 已在脚本内部用
  `SBATCH_RES=( -c 8 --mem=128G -t 04:00:00 )` 兜底。

### 备份位置

- `results/_backup_pre_orient_fix/` — 2026-05-07 启用 `both` 前的 MicrobeAtlas 产物
  （taxonomy/ + 8 个 gg2.*.h5ad + qc_plots/）。重跑后核对无误可删。
