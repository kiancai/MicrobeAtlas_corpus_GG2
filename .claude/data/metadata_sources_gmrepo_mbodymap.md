# 候选 metadata 源：GMrepo + mBodyMap（同组互补库）

> 状态：**仅 metadata 增强源，不贡献新样本行。合并动作暂缓**——等 MicroBioMap / Qiita 等样本行扩增完成后再统一 join（届时可贴的样本会更多）。本文件记录分析结论与"即用即合"规格。
> 原始文件：`rawdata/GMrepo/`（gmrepo.sql.gz + 7 个 txt.gz）、`rawdata/mBodyMap/`（mbodymap_samples.zip / mbodymap_projects.zip）。
> 分析日期：2026-06-04。所有数字由 `sample_to_run_info`（GMrepo）与 `mbodymap_samples.csv`（mBodyMap）实测。

## 1. 两库关系：同组、同 schema、按部位互补（不重叠）

- 同一课题组（陈卫华 / evolgenius），**同一套引擎与 30 列 schema**。28 列完全相同（含顺序）；唯一差异：GMrepo 用 `{checking, more}`，mBodyMap 换成 **`{BodySite, BodySite_info}`**。
- **GMrepo = 肠道专库**（仅 gut/stool）：更深（疾病 marker / 共现 / LEfSe）、更新勤（v1/v2/v3，已 118,865 样本）。
- **mBodyMap = 全身广度库**（22 部位）：一次性 2021 释放，63,148 样本，每部位较浅。
- **两库 run_id 重叠 = 0**（连 mBodyMap 的肠道系 11,464 个也不在 GMrepo）——故意不重叠地把人体样本按部位切开。**互补，不冗余。** 这就是"为什么分两个库"的答案。

## 2. 规模与内容

| | mBodyMap (2021) | GMrepo v3 (2026) |
|---|---:|---:|
| 样本（run 级） | 63,148 | 118,865 |
| Amplicon / Metagenomics | 48,747 / 14,401 | 87,018 / 31,847 |
| run 前缀 | SRR 59,697 / ERR 3,076 / DRR 375 | SRA+ENA 系 |

**mBodyMap BodySite 分布**（最大价值，GMrepo 没有此维度）：
Vagina 16,835 · Oral 9,451 · Rectum 9,060 · Skin 7,103 · 上呼吸道 6,830 · Lung 5,467 · Colon 1,711 · Nose 1,553 · Blood 1,154 · Urethral 1,102 · Cervix 858 · Uterus 511 · Esophagus 274 · Ear 269 · Ileum 266 · Trachea 201 · Peritoneal fluid 140 · Fallopian tube 130 · Stomach 116 · Ovary 80 · Duodenum 22 · Cecum 15

**疾病编码**：MeSH D-码（`disease` 列）+ 人读名（`phenotype` 列）。`D006262 = Health`（两库都约 70% 为健康对照）。解码表可用 gmrepo.sql.gz 里的 `mesh_data` 表，或直接用 `phenotype` 列自建 code→name 映射。
mBodyMap phenotype top：Healthy 23,476 · cystic fibrosis 2,046 · respiratory tract infections 1,966 · Pouchitis 1,905 · Bacterial Vaginosis 1,320 · COPD 1,171 · NEC 1,130 · Asthma 870 · endometrial cancer 660 …

## 3. 字段填充率（互补）

| 字段 | mBodyMap | GMrepo | 备注 |
|---|---:|---:|---|
| BodySite | **100%** | —（默认 gut） | mBodyMap 独有 |
| disease(MeSH) | **100%** | 91% | |
| phenotype | 73% | **91%** | |
| sex | **75%** | 32% | mBodyMap 强 |
| country | 38% | **100%** | GMrepo 强 |
| collection_date | 29% | **64%** | GMrepo 强（时序用） |
| 经纬度 | 17% | **56%** | GMrepo 强 |
| host_age | 15% | 28% | 都偏稀 |
| BMI / diet / 抗生素 | 5% / 1% / 0.1% | 1% / 0% / 0% | **都很差，别指望** |

**可稳定提供**：部位 + 疾病/表型 + 性别 + 地理。年龄/BMI/饮食/抗生素太稀。

## 4. 与现有语料重叠（决定能贴多少）

> **⚠ 本节为旧估（MA-only / 旧 corpus）。2026-06-05 已对真正五库并集重算**：
> mBodyMap **可贴 60,479 / 95.8%**、GMrepo **可贴 39,007 / 32.8%**、两库 run 重叠 = 0、可富集语料样本(去重)= 99,486。
> **以新审计为准 → `.claude/data/corpus_sample_overlap.md` §4。** 下表保留作历史参照。

基于 MA 原始 269 万 run 列表（`/tmp/ma_runs.txt`）估算；落地前需对**最终语料 obs['Run']**重算精确值。

| | 总样本 | ∩ MA | ∩ RM | ∩ PRIME | 覆盖率 |
|---|---:|---:|---:|---:|---:|
| mBodyMap | 63,148 | 60,398 | 10,557 | 0 | **95.6%** |
| GMrepo v3 | 118,865 | 32,392 | 41 | 215 | **27%** |

- **mBodyMap 96% 已在 MA** → 纯 metadata 增强，可给约 6 万已有样本补「部位+疾病+性别」（MA 自身部位标注弱，正好补）。公开下载仅 metadata CSV（无逐 run × taxon 全表），**0 可加净新增行。**
- **GMrepo 仅 27% 在 MA** → 只能给那 ~32k 已有肠道样本补疾病/地理。

> **更新（2026-06-04，对最终 corpus `obs['Run']`=1,793,428 重算）**：GMrepo 带可下载丰度的 16S amplicon run = 41,730，其中 12,275 已在 corpus、**29,455（70.6%）是 net-new**。即"0 净新增行"措辞不准——确有约 2.9 万 net-new 样本，但 GMrepo 丰度 keyed 在 NCBI taxon（仅 2,213 属）+ **无任何代表序列**（SQL dump 已确认无 seq/asv/otu 表），无法确定性对齐到 corpus 的 GG2/GTDB 8,114-genus 空间（NCBI↔GTDB many-to-many）。**否决结论不变**，但理由从"无净新增"修正为"有净新增但分类法不兼容 + 无序列 → unmappable"（与 MGnify-V5、WGS compendia 同类）。详见 `.claude/explorations/2026-06-04_dataset_expansion.md`、`.claude/data/mgnify_evaluation.md`。

### 关键洞察：MA 覆盖率 = 时间快照效应
mBodyMap(2021) 96% → PRIME(2025) 49% → GMrepo v3(2026) 27%，单调下降。**MA 抓取快照实测约停在 2021-01 初（最新公开 run SRR13347130 / first_public 2021-01-04，硬天花板；见 `.claude/data/pipeline_audit_microbeatlas.md`「MA 数据快照时间」）；越新的库 MA 漏得越多。** 推论：净新增集中在 2021 后的新资源 → 主力扩展锁 MicroBioMap(2024–25)/最新 Qiita 是对的。

## 5. 即用即合规格（等行扩增完成后执行）

1. **统一成一张 run 级 metadata 表**：两库 run 不重叠，直接 concat；加 `meta_source ∈ {GMrepo, mBodyMap}` 列。
2. **规整列名**：GMrepo 用下划线、mBodyMap 用点（`Recent.Antibiotics.Use` vs `Recent_Antibiotics_Use`），统一。
3. **携带字段**：`run_id`(主键) · `sample_id`(BioSample,备用键) · `BodySite`/`BodySite_info`(仅 mBodyMap) · `disease`(MeSH)+`phenotype` · `sex` · `country` · `collection_date` · `longitude`/`latitude` · `experiment_type` · `instrument_model` · `disease_stage`。
4. **join 到 obs**：按 `Run` 贴；mBodyMap 主要补非肠道部位，GMrepo 补肠道疾病/地理。加前缀（如 `GMrepo_*` / `mBodyMap_*`）避免与现有列冲突。
5. **解码**：MeSH D-码 → 疾病名（`mesh_data` 表或 phenotype 自建映射）。

> 为何暂缓：后续 MicroBioMap/Qiita 加入更多样本行后，这两库可贴中的样本数会上升，届时一次性 join 才划算。

## 6. 其他已排查的 metadata 资源（2026-06 扫描，结论：16S 已覆盖）

逐个核实后，**16S 逐样本 metadata 空间基本已被 GMrepo + mBodyMap + PRIME 覆盖**，无重要新增：

- **PRIME**（已有）：53k 16S，93 部位 / 101 表型，逐样本临床+协议 metadata 最丰富；作为 Tier1 直接合并源，metadata 随表带入。
- **HumanMetagenomeDB**（NAR 2021，69,822 样本）：curated 人体 metagenome metadata（部位用 BRENDA Tissue Ontology、诊断用 Human Disease Ontology）。**但是 WGS 模态** → 与我们 16S 行不同 run，贴不上；属 WGS 路线（已降级）的元数据源。
- **curatedMetagenomicData**：同上，WGS，临床 schema 金标准；WGS 路线再用。
- **MicrobiomeHD**：~5-6k，16S case/control + OTU 序列；量小，可在行扩增阶段顺带（既补行又补疾病标签）。

**另一类：微生物↔疾病关联库**（Disbiome / MicroPhenoDB / BugSigDB / Peryton / GMMAD / Amadis / HMDAD）= 文献策展的「菌-病」知识库，**无逐样本 accession**，不能贴到语料行；属**下游解释/评测**资源（验证模型学到的 genus-disease 信号），留到评测阶段用，不是语料 metadata。

**最完整的未用 metadata 源其实是 NCBI BioSample/SRA 属性本身**（host_age/sex/disease/body_site…），可对**全部**现有样本回填，而非某个策展子集——这是一个 curation 动作，不是"再加一个库"。

**结论：16S metadata 源已覆盖，可进入样本行扩增阶段。**
