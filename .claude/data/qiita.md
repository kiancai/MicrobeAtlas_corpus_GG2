# Qiita 数据库：结构 / 规模 / 访问 / 组成表计划 / 去重资产

> 稳定参考。2026-06-05 实测。叙事过程见 `explorations/2026-06-05_qiita.md`；合并架构见 `merge_architecture.md`。
> **一句话**：Qiita 公开包 = 873 study 的微生物数据，784 含 16S 且基本纯 16S；32.4 万 16S 样本有 ASV 能注释。**Layer1 已完成（2026-06-05）**：`qiita.gg2.genus.h5ad` **321,502 × 6,954**（全注释/不筛不去重，去重交 Layer2）。流水线见 §4。

## 1. 下载物（本地）

- `rawdata/Qiita/QIITA-public.tgz` — **99G 压缩**（解压 ~127G BIOM + 921M templates + 1.1G processed_data）。
- Qiita 官网「Download public BIOM and metadata files」整包；快照 **2026-05-31**；MD5 `8785599b1d0843453fddfaa637306d99`。
- 已解压（占 ~128G，盘有 2.4T）：`rawdata/Qiita/{BIOM/, templates/, processed_data/}`。

### 包结构
```
BIOM/<artifact_id>/
   reference-hit.biom / all.biom   # Deblur ASV 表（feature = DNA 序列）★ 我们要的
   otu_table.biom                  # 闭合参考 OTU（无序列，用不了）
   per-gene/ko/ec/pathway.biom     # WGS Woltka（功能，非 16S）
templates/
   <study>_<ts>.txt                # 样本信息（环境/host/empo/geo/timestamp…）873 个
   <study>_prep_<prepid>_<ts>.txt  # prep 信息（target_gene/platform/run_prefix…）2044 个
processed_data/<id>_otu_table.biom # 旧格式闭合参考
```
- 18,471 个 artifact 目录；其中 **7,625 个是 Deblur 序列表**（reference-hit / all）。
- biom 计数：otu_table 4760 / all 3839 / reference-hit 3790 / per-gene 791 / ko·ec·pathway 各 655 / none 1074 / free 870 / feature-table 580。
- **prep 的 `run_prefix` 不是 SRA/ENA 号**（是测序文件名前缀，如 `Knight16sV4_..._L007_R1_001`）→ 模板内无 accession，去重靠 ENA 桥（§5）。

## 2. 规模与组成（实测）

| | 数 |
|---|---:|
| 有 prep 的 study | 873 |
| 含 16S 的 study | **784**（纯16S 742 / 16S+其它混合 42） |
| 无 16S（纯 ITS/18S/WGS） | 89 |
| **16S 样本（有 16S prep）** | **362,011** |
| └ 有 Deblur ASV 表 = **能注释** | **324,220** ★ |
| └ 仅 OTU 表/无序列 = 注释不了 | 37,791 |

- **方法基本不混**：742/784 纯 16S；42 个混（几乎都是加做真菌 ITS / 真核 18S）。混在**样本/prep 层**（如 American Gut 10317：16S 35,414 + ITS 569 + 18S 569）→ **必须按样本/prep 粒度筛 16S，不能整 study 拿**。
- **trim context**（redbiom 全 host）：Deblur V4 **90nt 320,958 / 6.9M ASV**；100nt 313,690；**150nt 255,239 / 11.5M ASV**；125nt 34,208；非 V4 极边缘（V3-V4 仅几百）。同一样本进多个 trim → ingest 时**每样本选一个 context**（建议 90nt：覆盖最广 + 与 GG2 的 V4 ASV 口径一致）。
- redbiom 口径：public 16S 388,189；可见性 public 435,633 / private 119,074 / sandbox 661,175（private+sandbox 匿名不可得）。

### 环境/部位（全部 16S 样本，empo）
- 大类：**Host-associated ~72%（25.9万）** / Free-living（环境）~14%（4.8万）/ Control 空白 ~5%（1.8万）。
- empo_3 top（注意 EMPO 把**人也归在 "Animal"**）：Animal distal gut（含人肠道）104,496 · Animal(non-saline) 75,224 · Animal secretion 23,127 · Animal surface 21,882 · Animal corpus 15,138 · Soil 9,522 · Water 8,685 · Sterile water blank 11,265(控制)。

## 3. metadata（丰富，远好于 GMrepo/mBodyMap）

覆盖率（分母 = 36.2 万 16S 样本）：

| 字段 | 覆盖 | 字段 | 覆盖 |
|---|---:|---|---:|
| sample_type / scientific_name | 99% | **host_subject_id** | **97.5%** |
| latitude/longitude | 98.7% | collection_timestamp | 89.1% |
| env_material | 91.2% | geo_loc_name | 85.4% |
| empo_1 / empo_3 | ~82% | host_taxid | 78.0% |
| host_body_site | 68.8% | sex | 50.2% |
| country | 42.7% | host_age | 39.9% |

- **同人/时序任务金矿**：`host_subject_id` 97.5% + `collection_timestamp` 89%。
- 偏弱：年龄/性别/国家。
- 注：`scientific_name` 列存的是 `human gut metagenome` 这类（样本本体），不是 host 物种名；host 用 `host_taxid`（9606=人）。

## 4. 组成表（Layer1，**已产出 2026-06-05**，按 merge_architecture 契约）

**产物**：`results/feature_table/qiita.gg2.genus.h5ad`（77M）= **321,502 样本 × 6,954 genus** counts。
- 总 reads 134 亿；真 genus 5,745（+ 1,209 浅层 `g__` 桶）；零行 13（保留）；obs 重名 0。
- **全样本**（不筛/不归一/不去重）、taxon 清理（去 mito/chloro/非 B&A）、var = 6 级路径含浅层 `g__` 空桶 —— 与 MA/RM/MBM **完全一致**。
- 覆盖 321,502 / 362,011 16S 样本：其余 = 37,791（只有 OTU 表无序列）+ 2,718（混合 study 无法确信锁定 16S 表，Jaccard 模糊 → 弃，避免拿到 18S/ITS）。

**obs 主键差异**（关键，影响 Layer2 去重）：obs_names = Qiita `study.sample`（如 `10317.000107576`），**不是 run**。Run 经 ENA 桥（§5）映射，仅 ~30.9 万有 run，其余 Layer2 标 `dedup_uncertain`。

### Layer1 流水线（`scripts/databases/Qiita/`，**B 方案 = 全局去重**，已跑通）

选 B（不选 per-study A）原因：per-study 累计 2567 万 ASV + 最大单 study 233 万（怪兽单作业），B 全局去重后 1670 万、均匀切片、无怪兽。**结果等价**（NB 逐序列确定性）。

| 步 | 脚本 | 作用 / 实测 |
|---|---|---|
| 00a | `00a_select_artifacts.py` | 每 16S 样本选 1 张表（90nt reference-hit 优先）；混合 study 用 **artifact↔prep Jaccard 锁 gene=16S**（纯 16S 的 742 study 全算 16S）。→ 321,502 样本 / 1,259 artifact |
| 00b | `00b_build_per_study.py` | 按 study 出 `work/studies/<study>/{asv.fa, asv.h5ad}`（685 study） |
| 01 | `01_global_asv.py` | 全局去重（md5 头）→ `work/asv_global.fa` = **16,737,797** 唯一 ASV（去重率 34.8%；长度 90nt 1513万 + 100/150/250nt 少量） |
| 02 | `02_shard_fasta.sh` | 切 67 片（250k/片） |
| 03 | `03_classify_shard.sbatch` | GG2 24.09 **full-length NB `--p-read-orientation both`** array 分类（≡RM/MBM）。⚠ **含长 ASV(150/250nt) 的片在 128G 会 OOM；用 400G**（本次 5/67 片 OOM，400G 救回） |
| 04 | `04_merge_taxonomy.py` | 67 片 → `work/taxonomy_gg2.global.tsv`（md5→Taxon，校验 ==asv_global.fa） |
| 05 | `05_build_genus.py` | per-study collapse（parse_gg2 + Domain∈{B,A} + 去 mito/chloro + 6 级 var_id 聚合，**= MBM05/RM03 逐字**）；按 md5(seq) 查全局 taxonomy → `work/studies/<study>/genus.h5ad` |
| 06 | `06_concat_genus.py` | outer-join var + vstack → `qiita.gg2.genus.h5ad` |
| — | `run_finalize.sbatch` | SLURM `--dependency=afterok:<分类array>` 串 04→05→06（256G/计算节点） |

> **GG2 用 full-length 分类器**（区域+长度无关）→ 故 90/100/150nt 混合、V4 与否都无碍，同 RM。`assess/{01_dedup_assess,02_inventory_biom,03_profile}.py` 是探索脚本（净新增/盘点/profile），非流水线。

## 5. ENA 桥（Qiita 样本 → run，给 Layer2 去重用）

Qiita 自动把数据缴存到 ENA（center=UCSDMI），ENA 的 `sample_alias` = `qiita_sid_<study>:<study.sample>` 内嵌 Qiita 样本名，并带 run_accession + first_public + BioSample。这是把 Qiita 接到语料 run 空间的唯一干净桥。

- 已拉全表：`rawdata/Qiita/qiita_ena_bridge.tsv`（**372,512 run / 51M**）。
- 拉法（集群代理可达 ENA）：
  ```bash
  export http_proxy=http://127.0.0.1:60124 https_proxy=同
  curl -sG "https://www.ebi.ac.uk/ena/portal/api/search" \
    --data-urlencode 'result=read_run' \
    --data-urlencode 'query=center_name="UCSDMI"' \
    --data-urlencode 'fields=run_accession,sample_alias,sample_accession,secondary_sample_accession,study_accession,first_public,library_strategy,library_source,instrument_platform,scientific_name' \
    --data-urlencode 'format=tsv' --data-urlencode 'limit=0' -o qiita_ena_bridge.tsv
  ```
- 解析：`sample_name = sample_alias.split(":",1)[1]`（去 `qiita_sid_X:` 前缀）。

### 去重基线（Layer2 用）
| run 集 | 数 | 来源 |
|---|---:|---|
| corpus obs['Run'] | 1,793,428 | 当前语料（MA+RM）；`/tmp/corpus_runs.sorted.txt`（易失，可从 MCFCorpusV2 重抽） |
| MicroBioMap run | 168,464 | `microbiomap.gg2.genus.h5ad` obs_names（∩corpus 106,474 → MBM 净新增 61,990） |
| MA 原始(QC前) | 2,690,735 | `/tmp/ma_runs.sorted.txt`（MA 见过但 QC 剔的判定用） |

### 净新增实测（仅 Layer2 参考，**不作 Layer1 注释门槛**）
- Qiita 16S net-new（不在 corpus∪MBM）= **223,011**；其中**有 ASV 表 = 185,642**。
- 分档：🟢 **铁(MA 2.69M 也没见过, ~全 post-2021) = 81,304** / 🟡 MA见过但QC剔 = 39,372 / ⚪ 无桥(未验证) = 102,335。
- 已落盘逐样本：`rawdata/Qiita/qiita_netnew_16s.tsv`（sample_name / study / status / run / first_public / empo / host_taxid）。

## 6. GG2 注释路径（待 pilot 定）

| | NB 分类器 | exact-match |
|---|---|---|
| 本地文件 | ✅ `tools/qiime2_database/greengenes2/2024.09.backbone.full-length.nb.sklearn-1.4.2.qza` | ✗ 需下 `2024.09.taxonomy.asv.tsv.gz`（ftp.microbio.me/greengenes_release/2024.09/） |
| 一致性 | 与 MA/RM/MBM 同一套 ✅ | 新方法（genus 或略有出入） |
| 算力 | 重（唯一 ASV 数百万 → array job） | 轻（序列哈希查表） |
| 90nt 短 ASV | 一般 | 高（这些 ASV 本就是 GG2 训练输入） |

- 默认倾向 **NB（一致性 + 本地现成）**；exact-match 作加速候选。**建议先 1 个 study pilot 两路对比 genus 一致性再定**。
- 其它 GG2 本地：`tree/2024.09.{taxonomy.id.nwk, taxonomy.id.tsv.gz, phylogeny.id.nwk}`。

## 7. 访问方式备查（复用）

- **redbiom**（跨 study 查询/计数）：`pip install redbiom`；webdis 服务器 **`http://qiita.ucsd.edu:7329`**（设 `REDBIOM_HOST`）。`search metadata`/计数可用；**批量 `fetch sample-metadata` 大批量会 `MGET: b''`**（服务器拒超大批，需切 ≤数千的小块；本项目已改用整包 templates 拿全 metadata）。**项目 ncat 代理(60124) 会弄坏 redbiom keepalive**，须在能直连 7329 的干净外网节点跑（用户 Mac 代理 7890 可用）。
- **HTTP 按 study 下载**（无需登录）：`https://qiita.ucsd.edu/public_download/?data={biom|sample_information|prep_information|raw}&study_id=X[&data_type=16S][&prep_id=Y]`。`raw` 是原始 fastq（巨大，勿下）。
- **整包**：官网「Download public BIOM and metadata files」= 上面那 99G tgz（已下）。
- **ENA Portal API**：见 §5，集群代理可达。

## 8. 本会话产物

- 脚本 `scripts/databases/Qiita/`：`01_dedup_assess.py`（净新增三方去重）/ `02_inventory_biom.py`（Deblur 表盘点）/ `03_profile.py`（study×方法 + metadata 覆盖 + 能注释总量）。
- 中间件 `rawdata/Qiita/`：`qiita_ena_bridge.tsv`(51M) / `qiita_netnew_16s.tsv`(34M) / `qiita_biom_inventory.tsv`（逐 artifact：study/n_samples/n_16s/n_netnew/n_asv/trim）/ `templates/`(921M) / `BIOM/`(127G) / `processed_data/`(1.1G)。
