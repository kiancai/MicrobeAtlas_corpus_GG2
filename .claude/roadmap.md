# 项目路线图（准备做 / 进行中 / 已完成）

> 任务级状态板。每次"落实有效信息"时更新三态迁移；探索的过程与结论记到 `.claude/explorations/`，稳定事实记到 `.claude/data/`。
> 最近更新：2026-06-05。

## ✅ 已完成

- **MicrobeAtlas 流水线 01–08**：OTU → GG2 genus anndata → QC → metadata → 过滤。产物 `gg2.full.qc.*`（Stage 1：1,762,635 × 6,306）。
- **ResMicroDb 流水线 01–08**：398 项目 ASV → GG2 genus → QC → metadata 标准化(+06b NCBI 真值修正) → 过滤。产物 `resmicrodb.gg2.genus.qc.*`（Stage 2：93,425 × 4,952）。
- **跨库合并 09–12**：outer-join → 扩到 8,114 GG2 全属 + 挂 taxo/phylo 距离 → finalize 相对丰度 → 挂 Poincaré 欧氏位置编码。**最终语料 `MCFCorpusV2.gg2.h5ad` 1,826,126 × 8,114**。
- **Phylogeny 距离矩阵**：GG2 24.09 全 8,114 genus 的 taxo(0-6 hop) + phylo(patristic) 距离，独立副产物 + 挂回 anndata。
- **Sample×sample 距离**：50k 分层子集（BC + wUniFrac + PCoA）+ 100k 独立分支（corpus 级 PCoA 可视化）。
- **数据库扩增调研（2026-06）**：评估 ~77 候选库，定下扩增策略与顺序。详见 `explorations/2026-06-04_dataset_expansion.md`。
- **MicroBioMap 重注释（2026-06-04）**：Human Microbiome Compendium v1.1.1（168,464 人肠道 16S / 366万全局 ASV）→ GG2 24.09。366万 ASV 分类（31片 array，配方逐字复用 RM）→ genus collapse（RM 03 逐字）。产出**完整组成表 `microbiomap.gg2.genus.h5ad` 168,464 × 5,899**（全样本保留、未筛选未合并）。净新增实测 61,990/168,464(36.8%)。脚本/产物见 `pipeline.md`「MicroBioMap」+ `scripts/databases/MicroBioMap/README.md`。
- **Qiita 重注释（2026-06-05，Layer1）**：整包 99G、**全样本不限人体** 16S。流程（B=全局去重，复用 RM/MBM）：选表(每样本1张90nt 16S，混合study prep-Jaccard锁16S，321,502样本)→ 全局去重 **1670万 ASV** → 67 片 GG2 **full-length NB(both)** 分类（5片长ASV OOM→400G救回）→ per-study collapse(RM03/MBM05逐字)→ 拼表。产出 **`qiita.gg2.genus.h5ad` 321,502 × 6,954**（counts/全样本不筛不去重；134亿 reads；真 genus 5,745）。去重资产 ENA 桥(372,512 run)留 Layer2。脚本/产物/坑 见 `data/qiita.md` + `pipeline.md`「Qiita」+ `explorations/2026-06-05_qiita.md`。

## 🔄 进行中

- **语料合并重构（merge-then-filter）已定，P1 进行中**：改为「各库 → 标准 genus 组成表(counts) → 一次性统一合并 + metadata + 筛选」。架构/组成表契约/已锁定&待定决策 → **`.claude/merge_architecture.md`**。
  - 已锁定：merge-then-filter；genus 粒度；unresolved 存 `obs['unresolved_reads']`(不进 token 矩阵)；metadata 统一层一起合；旧码归档不删。
  - ✅ **QC = 打标不删（2026-06-05 定）**：Layer2 一行不删，QC 降级为 obs 旗标 `qc_pass/n_deep_reads/n_deep_features`，去重为 `Run_paired_id` 标记。依据：旧 QC 会删 100 万行(30%)且 48,526 run 跨库判决分歧 → 全留+旗标避免误删、下游再筛。
  - ✅ **五库组成表全部就位（2026-06-05）**：MA `gg2.full`(2,690,735×7,424) · RM `resmicrodb.gg2.genus`(100,391×5,891) · MBM `microbiomap.gg2.genus`(168,464×5,899) · Qiita `qiita.gg2.genus`(321,502×6,954) · **PRIME `prime.gg2.genus`(53,449×6,793) ← 新建**（Tier1 reshape，脚本 `databases/PRIME/01`）。
  - ✅ P1 + 目录重构已落地：`scripts/merge_v2/` 骨架 + 架构文档；旧 Merge → `archive/Merge/`；MA/RM/MBM/PRIME → `databases/`。
  - ✅ **`01_concat` 已跑（2026-06-09）** → `results/feature_table/merged.counts.h5ad` **3,334,541 × 8,114**（int32 counts，全留不删；qc_pass 70%；unresolved 中位 6.9%；跨库副本 693,701 行/331,125 Run；Qiita 无 run 74,148）。⚠ var 暂带旧统计列待 04 重算。
  - ✅ **各库 metadata 标准化（2026-06-09）** → `results/feature_table/metadata/{ma,rm,mbm,prime,qiita}.std.parquet`（公共 27 列统一 + DB 前缀；脚本 `scripts/merge_v2/metadata/*_standardize.py` + 共享 `common_schema.py`/`iso_country.py`）。Age_Group 中央右闭分桶复刻 RM 100%。
  - ✅ **02_merge_metadata（2026-06-09）** → `metadata/merged.obs.parquet` 3,334,541 × 185（旗标+5库metadata+GMrepo/mBodyMap纯前缀富集）。
  - ✅ **04_finalize（2026-06-09）→ 🎉 最终语料 `MCFCorpusV3.gg2.h5ad` 3,334,541 × 8,114**：X=相对丰度·layers['counts']·obs 185列·varp(taxo/phylo 距离)·varm(Poincaré 位置编码,从 V2 拷贝)·uns 口径。
  - ✅ **回读核验全过（2026-06-09，`/tmp/verify_v3.py`）**：reads 守恒(差=0)·相对丰度行和=1·varp/varm 与 V2 逐值一致·配对样本独立(同 Run 跨 4 库各留一行)。
  - ✅ **公共列 sanity 清洗（2026-06-09）**：剖析发现越界值，回溯判定 bug vs 源脏（详 `metadata_merge_design.md`「sanity」段）→ 修 2 处 gap（Region 大写 / inf）+ Lat/Lon/Age/BMI/Height/Weight 越界裁剪 + (0,0)/999909 哨兵→NA（`common_schema.apply_common_sanity`，02 调用）→ 重跑 02+04 出干净 V3。
  - ✅ **HF 数据字典 中英双版（2026-06-09）** → `results/feature_table/DATA_DICTIONARY_MCFCorpusV3{.md,.zh.md}`。**自包含**(含数据集构成 §1/QC §2/重复 §3/结构/全列/约定/来源署名),不引用内部脚本。QC=旗标(=旧标准,验证逐字复现)、重复=Run_paired_id(画像:331,125 组/693,701 行/36,768 组 QC 分歧)均写入。
  - ✅ **最终 ∩MGnify（2026-06-09）**：V3 唯一 run 2,897,817 ∩ MGnify = **473,988 / 88.3%**（qc_pass 子集 64.5%）；MGnify 独有 62,578(11.7%)。
  - ✅ **Qiita 脚本已并入 `scripts/databases/Qiita/`**（databases/ 现五库齐）。
  - 🟢 **主线 + 交付物完成**（语料 V3 + 双版字典）。用户定：**不做训练视图**（筛选/去重交下游）、暂不删盘。剩纯可选：磁盘清理（清单见 `corpus_sample_overlap`/对话）。
  - ✅ **全局样本审计已出（2026-06-05）** → `.claude/data/corpus_sample_overlap.md`（脚本 `scripts/merge_v2/00_overlap_audit.py`）：五库 Run 并集 2,923,879 + Qiita 无 run 74,148 ≈ **300 万唯一样本(pre-QC)**；**∩MGnify 88.7%**；mBodyMap 可贴 95.8% / GMrepo 32.8%，两库 run 重叠 0。
  - ✅ **全局样本审计已出（2026-06-05）** → `.claude/data/corpus_sample_overlap.md`（脚本 `scripts/merge_v2/00_overlap_audit.py`）：五库 Run 并集 2,923,879 + Qiita 无 run 74,148 ≈ **300 万唯一样本(pre-QC)**；净增 MA→+RM 58k→+MBM 49k→+PRIME 22k→+Qiita 104k；**∩MGnify 88.7%**；mBodyMap 可贴 95.8% / GMrepo 32.8%，两库 run 重叠 0。
- **metadata 合并设计：逐库逐列过 进行中** → **`.claude/data/metadata_merge_design.md`**。
  - 原则：列级「共同含义」才去前缀，否则 `<DB>_` 前缀；自底向上逐列定；公共-派生从公共原始量中央算（如 Age_Group←Age）；id 槽列值加 DB 命名空间。
  - ✅ 已过 **MA(26) + RM(36) + MBM(3 表) + PRIME(84+项目) + Qiita(3 来源)** → 公共列候选集已确立（Run/BioSample/Project_ID/Lat/Lon/Sex/Age(+桶)/BMI/Country/Continent/Platform/Instrument/Library_Layout/Region_16S/Subject_ID/PMID/Smoking/Sequencing_Type{Amplicon,Shotgun}/**Collection_Date/Height/Weight**）。
    - MBM：iso→Country/Continent、External Id→BioSample(仅 ENA 子集)、projects.amplicon→Region_16S、tags 选键→Sex/BMI/Age/Subject_ID，余 MBM_；丢 pubdate/total_bases/library_source。
    - PRIME：两伪等价陷阱（`Sequencing_Type`=读长配置→`PRIME_Read_Config`；`Host_Age` 多分桶→仅数值进 Age）；真 BioSample 全覆盖、Doi→PMID；Smoke_Status 映射；余 PRIME_；丢 AvgSpotLen/Bases/Release/Create_Date。
    - Qiita：三来源（ENA 桥+sample 模板+prep 模板）；obs 主键 `study.sample`（原生,永远有），**run+BioSample 成对的 INSDC 叠加层、同覆盖 ~85.6%**（无号样本未缴 ENA→Layer2 `dedup_uncertain`）；country 剥 GAZ:、host_age 单位换算；环境本体/body site/host 物种全 QIITA_。触发两项升公共：**Collection_Date(3 peer)·Height/Weight(2 peer)**。
    - GMrepo/mBodyMap：同组姊妹库纯前缀富集（`GMREPO_*`/`MBODY_*`，不碰公共列）；2026-06-09 用户确认仅-metadata 库**就这两个**。
  - ✅ **列设计全部收尾** → 进入实现：各库 `standardize` → `02_merge_metadata`（公共 23 列 + DB 前缀 + 富集纯前缀）。🔶 跨库重合已审计（`corpus_sample_overlap.md`）；实现完 → HF 数据字典。
- **Qiita Layer1 ✅ 完成**（2026-06-05，见已完成区）：脚本在 `scripts/databases/Qiita/`，待 Layer2 cleanup 时并入 `databases/Qiita/`。
- **MGnify 已收口**（2026-06-05）：不并入语料，全库 536,566 run 列表已抓存（`rawdata/MGnify/`）转作**覆盖率基准**，见「准备做」末步 + `data/mgnify_evaluation.md`。
- 新原始数据已下载：`rawdata/{PRIME, GMrepo, mBodyMap}/`（尚未整合进语料）。

## 📋 准备做（按顺序）

1. **MicroBioMap 后续（待商议）**：重注释完整组成表已出（`microbiomap.gg2.genus.h5ad` 168,464×5,899）。剩下三件**逐项与用户商议后再做**：① sample 级筛选（剔零=RM 04 / QC=RM 05，用户当年对 MA/RM 做过、去过样本）② metadata 合并（sample_metadata + tags + projects；+ GMrepo/mBodyMap）③ 对齐 8,114 + 并进语料（concat + Run_paired_id + 搬 varp/varm + 重 finalize）。
2. ~~**扩样本行 · Qiita**~~ **✅ Layer1 完成（2026-06-05，见已完成区）** → `qiita.gg2.genus.h5ad` 321,502 × 6,954。剩余交 Layer2：metadata 合并 + 去重（ENA 桥，obs 主键 `study.sample` 非 run）+ 并进语料。
3. **合并 PRIME**（Tier1 直接合并 GG2 24.09 genus 表，零重处理）— 用户要求放到行扩增之后再并。
4. **统一合并 metadata**：GMrepo + mBodyMap（+ NCBI BioSample 回填）→ 按 Run join 到 obs。**暂缓至行扩增完成**（届时可贴样本更多）。规格见 `.claude/data/metadata_sources_gmrepo_mbodymap.md`。
5. **特殊任务数据（纵向 / 同人）**：等用户给 special tasks 清单后细化；两条腿 = 定向收设计即纵向的队列 + 回填现有语料的 subject_id/timepoint。
6. **（later）WGS 路线**（HumanMetagenomeDB / curatedMetagenomicData）+ **关联库评测**（Disbiome / BugSigDB / MicroPhenoDB）——下游评测阶段再用。
7. **最终覆盖率检查（合并完所有库后）**：用最终 corpus `obs['Run']` ∩ MGnify 全库 `rawdata/MGnify/mgnify_run_accessions.sorted.txt`(536,566) 算覆盖率，可按 experiment_type 分层报告。MGnify 仅作外部标尺，不并入。配方/实测见 `data/mgnify_evaluation.md`。

## ❌ 已评估不做

~~环境扩展~~（**2026-06-05 作废此条**：改全样本口径，已纳入源 Qiita 等里的环境/动物样本全留、不再剔；主动找"环境-only 新库"仍低优先）· raw-reads-only（中国 NMDC/CNSA/CNGBdb、GGMP、TwinsUK、LogMPIE、SKIOME）· 受控访问（TEDDY/FINRISK、FGFP/LifeLines-DEEP）· **MGnify 不并入语料**（V1–V5 无 rep-seqs 注释不进 GG2、v6 ASV 仍 beta API 0 条；2026-06-05 全库重合实测 ∩MA 85.3% = INSDC 重处理子集）→ **改作覆盖率基准**，见「准备做」末步 + `data/mgnify_evaluation.md`，非"什么都不做"· 几百样本的小集（Isala/FemMicro16S/MetaIBS 暂缓）· 纯参考库（HOMD/WoL2/GSR-DB）· 纯 WGS。理由见 `explorations/2026-06-04_dataset_expansion.md`。
