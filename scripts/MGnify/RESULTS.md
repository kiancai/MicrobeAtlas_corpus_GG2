# MGnify 全库 × 本项目语料：样本级重合实测（2026-06-05）

> 脚本：`scripts/MGnify/01_fetch_runs.py`（抓取）+ `02_intersect.py`（交集）。sample 列表 → `rawdata/MGnify/`。
> 性质：**MGnify 不并入语料，作覆盖率基准**。目的 = 量化"MGnify 是否只是 INSDC 重处理子集、对我们有无净新增"。
> 定位与用法（合并完所有库后算覆盖率）见 `.claude/data/mgnify_evaluation.md`。

## 0. 一句话结论

**MGnify 全库 536,566 个 run 里 85.3% 已在我们 MA 全集中（amplicon 更高，86.8%）**——MGnify 本质就是 INSDC 的一个重处理子集，与本项目高度重叠。净新增（不在 MA）全库仅 **79,110**、amplicon 仅 **58,122**，且这部分(a)多为 2021 快照后、会被 Qiita/MicroBioMap 覆盖，(b)是 V1–V5 的 SILVA-OTU、无代表序列、**注释不进 GG2 8,114 词表** → **实际可用净新增 ≈ 0**。**否决 MGnify 的决策被实测坐实。**

## 1. 数据规模（MGnify v1 API `/runs` 端点，已按 run 去重）

| 端点 | 数量 |
|---|---:|
| `/runs`（全部，已按 run 去重）| **536,566** |
| `/analyses`（含多版本流水线重复）| 635,107 |
| `/samples` | 435,812 |

experiment_type 分布（536,566）：amplicon 439,615 · metagenomic 72,658 · unknown 17,570 · metatranscriptomic 5,107 · metabarcoding 1,494 · (空) 122。
accession 前缀：SRR 269,240 · ERR 252,809 · DRR 14,502 · MET 13 · SRS 2。

## 2. 重合结果

参照集：MA 全集 `/tmp/ma_runs.txt` = 2,690,735 run；最终语料 corpus `obs['Run']` = 1,793,428 unique run。
匹配键 = MGnify primary accession，SRS/ERS 样本号补用 secondary_accession 兜底。

### 总体（全部 experiment_type）

| | 数量 | 占 MGnify |
|---|---:|---:|
| MGnify ∩ MA 全集 | **457,456** | **85.3%** |
| MGnify ∩ corpus | 328,079 | 61.1% |
| net-new vs MA | 79,110 | 14.7% |
| net-new vs corpus | 208,487 | 38.9% |

### 按 experiment_type 分解

| type | total | ∩MA | ∩MA% | ∩corpus | ∩corp% |
|---|---:|---:|---:|---:|---:|
| amplicon | 439,615 | 381,525 | 86.8% | 281,874 | 64.1% |
| metagenomic | 72,658 | 54,965 | 75.6% | 32,134 | 44.2% |
| unknown | 17,570 | 16,768 | 95.4% | 13,052 | 74.3% |
| metatranscriptomic | 5,107 | 2,615 | 51.2% | 687 | 13.5% |
| metabarcoding | 1,494 | 1,462 | 97.9% | 305 | 20.4% |
| (空) | 122 | 121 | 99.2% | 27 | 22.1% |

### amplicon + metabarcoding（最可比我们的 16S 语料）

| | 数量 | 占比 |
|---|---:|---:|
| 总数 | 441,109 | — |
| ∩ MA 全集 | **382,987** | **86.8%** |
| ∩ corpus | 282,179 | 64.0% |
| net-new vs MA | 58,122 | 13.2% |
| net-new vs corpus | 158,930 | 36.0% |

net-new vs corpus 的前缀分布：SRR 102,399 · ERR 97,829 · DRR 8,244 · MET 13 · SRS 2。

## 3. 关键解读

1. **85% 重叠 = INSDC 重处理子集的铁证**。这与 `2026-06-04_dataset_expansion.md` 核心洞察#4（"不要枚举数据库，它们多是同一 INSDC 池的重叠视图"）完全一致；也与时间快照规律吻合（MA 截至 2021-01，MGnify 大量旧分析落在窗口内）。
2. **∩corpus(61%) < ∩MA(85%) 的差 = 我们自己筛掉的**：457,456 − 328,079 = **129,377** 条 MGnify run 在 MA 爬到了、但在 QC 阶段被主动丢弃（corpus 比 MA 全集少约 90 万 run）。这部分**不算 MGnify 的"净新增"**——是"我们有但没要"。
3. **真·净新增 = 79,110（全库）/ 58,122（amplicon）**，但：
   - 受快照效应，其中后快照部分与 Qiita（后快照 INSDC ~10.8万）/ MicroBioMap（人肠道 16.8万）大面积重叠，不是 MGnify 独有；
   - V1–V5 只给 SILVA-OTU 表 + reads、**无代表序列** → 进不了"序列→GG2 24.09 NB→genus collapse"管线；v6 ASV 仍 beta、公开 API 0 条。
   - → **可干净并入 GG2 的净新增 ≈ 0**，与否决理由一致。
4. **WGS/metaT 重叠偏低（metaT 51%）属正常**：它们与我们 16S 是同一 biosample 的不同 run（不同 accession），run 级自然对不上；不代表"样本不重叠"。我们语料是 run-keyed 16S，run 级是正确粒度。

## 4. 方法与口径（可复现）

- 抓取：`/runs?page_size=250&page=N`，N=1..2147，header `Accept: application/json`，12 并发，每页 6 次重试。2147 页 0 失败，536,566 条全唯一。耗时约 20 分钟（深层页 ~9s/页）。代理 `127.0.0.1:60124`。
- 交集：纯集合 membership，primary accession + secondary_accession 兜底。
- 中间产物（均在 `/tmp`，非 repo）：`mgnify_runs.tsv`、`ma_runs.txt`、`corpus_runs.txt`。
- 口径注意：MGnify 的 SRS/ERS（仅 2 条）是 sample 号非 run 号，与 run-keyed 参照集对不上属预期；量极小可忽略。

## 5. 触发"重新考虑 MGnify"的条件（不变）

仍同 `mgnify_evaluation.md` §5：**v6 退 beta 且把人体 amplicon 回灌公开 API（带 `ASV_seqs.fasta`）**。即便如此，按本次实测，其净新增的大头也已被 MA/Qiita/MicroBioMap 吃掉，预期增量有限。
