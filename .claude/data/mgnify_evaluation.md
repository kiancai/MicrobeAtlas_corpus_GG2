# MGnify（EMBL-EBI）：覆盖率基准（held，不并入语料）

> 状态（2026-06-05）：**不并入语料，转作覆盖率基准**。已抓全库 sample 列表持有，等 MA/RM/Qiita 等合并完成后，仅用来算一次"corpus 对 MGnify 的覆盖率"。
> 不打算重分析 MGnify、不接进 merge-then-filter 流水线。
> 重合实测详表 → `scripts/MGnify/RESULTS.md`；上层调研 → `.claude/explorations/2026-06-04_dataset_expansion.md`。

## 1. 定位（为什么是"基准"而非"数据源"）

MGnify 本质是 **INSDC（ENA/SRA）的一个重处理子集**，与本项目高度重叠（见 §3 实测 ∩MA=85.3%），且 V1–V5 分析**无代表序列**（只给 SILVA-OTU 表 + reads）、能干净注释进 GG2 的 v6 ASV 流水线仍 beta、公开 API 0 条（见 §4）→ **没有可干净并入 GG2 8,114 词表的净新增**。

因此 MGnify 不作为扩增行源，而是**作为外部覆盖率标尺**：持有它全部 run 的 sample 列表，最终衡量"我们合并出来的语料覆盖了 MGnify 多少样本"。

## 2. 持有的资产（sample 列表）

- `rawdata/MGnify/mgnify_runs.tsv` — 全表 536,566 行：`accession, secondary_accession, experiment_type, ena_study, is_private`。
- `rawdata/MGnify/mgnify_run_accessions.sorted.txt` — 纯 run accession 排序去重，536,566 条，**覆盖率检查直接用**。
- 抓取/交集脚本 + 完整报告：`scripts/MGnify/{01_fetch_runs.py, 02_intersect.py, RESULTS.md}`。

规模（MGnify v1 API `/runs`，已按 run 去重，2026-06-05）：
- 总 run **536,566**（`/analyses` 635,107 含多版本流水线重复；`/samples` 435,812）。
- experiment_type：amplicon 439,615 · metagenomic 72,658 · unknown 17,570 · metatranscriptomic 5,107 · metabarcoding 1,494 · 空 122。
- 前缀：SRR 269,240 · ERR 252,809 · DRR 14,502 · MET 13 · SRS 2。

## 3. 重合实测（2026-06-05，基线判断 + 覆盖率预演）

参照：MA 全集 2,690,735 run；最终语料 corpus `obs['Run']` 1,793,428。

| | MGnify 全库 | ∩ MA 全集 | ∩ corpus |
|---|---:|---:|---:|
| 全部类型(536,566) | — | **457,456 = 85.3%** | 328,079 = 61.1% |
| amplicon+metabarcoding(441,109) | — | 382,987 = 86.8% | 282,179 = 64.0% |

- 净新增（不在 MA 全集）：全库 79,110 / amplicon 58,122——且无 rep-seqs 注释不进 GG2 → **可用净新增 ≈ 0**。
- ∩corpus(61%) < ∩MA(85%) 的差(129,377) = 我们爬到了但 QC 阶段主动筛掉的，不算 MGnify 净新增。
- WGS/metaT 重叠偏低（metaT 51%）正常：与我们 16S 是同 biosample 的不同 run，run 级天然对不上。
- 分类型/前缀完整分解 → `scripts/MGnify/RESULTS.md`。

> **覆盖率用法**：合并完所有库后，用最终 corpus 的 `obs['Run']` 对 `mgnify_run_accessions.sorted.txt` 取交集 / 占比，即"语料覆盖了 MGnify 多少"。可按 experiment_type（用 `mgnify_runs.tsv` 第 3 列）分层报告。

## 4. 背景：为什么 V1–V5 注释不进 GG2（不并入语料的技术根因）

实测一条 V5.0 amplicon 分析（`MGYA00462714`）可下载文件：Processed reads / SSU rRNA reads（未去噪）/ MAPseq SSU SILVA 注释 / SILVA OTU 表（TSV·HDF5·JSON Biom）。**全程无 representative-sequence FASTA** → 要么 SILVA→GG2 有损 crosswalk（与"GG2 而非 SILVA"决策冲突，否决），要么拿 reads 自己 de novo 去噪=重跑（用户不重跑 raw）。只有 **v6 的 DADA2 流水线产出 `ASV_seqs.fasta`** 能走本项目"序列→GG2 24.09 NB→genus collapse"管线，但 v6 仍 beta、公开 API 实质 0 条（v1 API 只到 V5.0；v2 采样 3000 分析 V6=0）。

人体 host-associated amplicon 分析数（v1 API，2026-06-04）= 222,637（按版本：V4.1 10.1万 / V5.0 8.5万 / V6.0 = 0）。是 analyses 非 unique run，去重后更少。

> 若将来 v6 退 beta 且把人体 amplicon 回灌公开 API（带 `ASV_seqs.fasta`），届时可重评是否值得并入；但按 §3 实测，其净新增大头已被 MA/Qiita/MicroBioMap 吃掉，预期增量有限。当前**不依赖此条**——MGnify 的角色就是覆盖率基准。

## 5. API 访问配方（复用）

- 走代理：`export http_proxy=http://127.0.0.1:60124 https_proxy=同`。
- 抓全库 run（本次用法）：`GET /metagenomics/api/v1/runs?page_size=250&page=N`，N=1..2147，Header `Accept: application/json`（用 `application/vnd.api+json` 会 406）；run accession = `data[].attributes.accession`（+`secondary-accession`）。不带 lineage 过滤，~4–9s/页。
- 计数：任意端点 `?page_size=1` 读 `.meta.pagination.count`。
- 人体过滤（v1，仅在需要人体子集时）：`?experiment_type=amplicon&lineage=root:Host-associated:Human`（带 lineage 过滤的深查 ~31s/页，慢）。
- v2 API 顶层 `/analyses` 不支持 query 过滤（`experiment_type`/`biome`/`pipeline_version` 被静默忽略），但含全部版本；`curl -L`（301）；用于采样估版本分布。
