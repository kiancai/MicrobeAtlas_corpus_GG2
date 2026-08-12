# Qiita / redbiom 作为净新增人体 16S 行源：effort + 风险评估

> 状态：**effort = high（偏 very_high 风险）**，但**结构上比 MGnify/GMrepo 干净**（有 deblur ASV 序列 → 可走既有 NB 管线）。
> 评估日期：2026-06-04。上层调研 → `.claude/explorations/2026-06-04_dataset_expansion.md`；本文件专注 effort 维度。
> 一手核实：redbiom mSystems 2019、Greengenes2 Nat Biotech 2023（PMC10818020）、Qiita stats live（2026-06）、ENA Portal API live、Qiita docs(contexts/reanalysis)。

## 0. 一句话结论

Qiita 人体 16S 在**技术上可并入**（deblur ASV 序列可下载 → 喂本项目 24.09 NB → genus collapse → 对齐 8,114），但 effort 是 high：真正的成本不在算力（V4 ASV 池已被 GG2 dedup 到 ~2,100 万、且本项目只需 genus 不需逐 ASV 精度），而在 **(a) 多 context 碎片化合并、(b) 人体 free-text 筛选、(c) 与 MA/MicroBioMap 的去重——而 Qiita 用 `studyID.sampleID` 不是 SRA run，去重需经 EBI accession 桥接**。净新增体量被快照效应 + MicroBioMap 重叠双重压缩，乐观估计**人体 net-new 仅约 3 万–7 万行**（远低于"10 万+"的表面规模）。

## 1. Qiita 当前实测规模（live, 2026-06，qiita.ucsd.edu/stats）

| 维度 | Public | Private | EBI 提交 |
|---|---:|---:|---:|
| Studies | 878 | 177 | 1,136 |
| Samples | 435,633 | 119,074 | 385,491 |

按数据类型（public，跨 study 计）：**16S = 388,189** · 18S 12,221 · ITS 14,747 · Metagenomic 105,908 · Metatranscriptomic 27,161 · 其它少量。
→ 16S 是主体；但 388k 含全部环境/宿主类型，**人体只是子集**（需 EMPO/host 筛）。

## 2. 净新增体量：被快照 + MicroBioMap 双重压缩

MA 快照截止 ≈2021-01-04（硬天花板，见 pipeline_audit_microbeatlas.md）。corpus `obs['Run']` = 1,793,428 unique run。

ENA Portal API 实测（center_name="UCSDMI" = Qiita 的 ENA 提交 center，library_strategy="AMPLICON"）：
- UCSDMI amplicon runs in ENA 总计 = **256,085**
- first_public ≤ 2021-01-04 = **148,307（58%）** → 落在 MA 快照窗口内，大概率已在 corpus
- first_public ≥ 2021-01-05 = **107,778（42%）** → 后快照，**INSDC 侧 net-new 天花板（含全部环境/宿主类型，未人体筛、未 MicroBioMap 去重）**

参照：ENA 全库 human-gut-metagenome(tax 408170) amplicon first_public ≥2021-01-05 = **331,317**（所有来源 net-new 上限）。Qiita 的 UCSDMI 后快照 107,778 是其中一片，且与 MicroBioMap（NCBI BioSample 人肠道 / amplicon / Oct 2021 截止 / 482 BioProjects / 168,484 样本）**大面积重叠**——American Gut/Microsetta 是两库都含的最大人肠道库。

**净新增推断**：107,778（后快照 INSDC Qiita）→ 人体筛（~60-80%）→ 减 MicroBioMap 已覆盖肠道（并行任务，重叠大）→ 减非肠道部位中已在 MA 的 → **人体 net-new 落点约 3 万–7 万行**（粗估，落地需精确差集）。
另有 119,074 private + 未提交 EBI 的样本 = 唯一真·不在 INSDC 的部分，但**不可下载/access-restricted + 无 INSDC accession 无法去重** → 多数拿不到。

## 3. Effort 分解（用户 5 问）

### (1) 多 context 碎片化合并 —— **medium-high**
- redbiom **context** = `<binning>-<refDB>-<seq strategy>-<region>-<trim length>-<id>`。同一生物样本按 V 区(V2-V3/V3/V4/V4-V6/...) × trim(90/100/125/150/200/250nt) × 工具(Deblur/closed-ref/DADA2) 散在多个 context。
- 后果：**不同 context 的 ASV 序列不可在 ASV 层直接合并**（trim 长度不同 → 序列不同；V 区不同 → 完全不可比）。
- 缓解：本项目终点是 **genus**，collapse 后跨 context 可比。所以策略 = 逐 context 各自 `fetch` → 各自跑 NB → 各自 collapse 到 genus → 再按样本 union。**但 V4 与 V3-V4 的 genus 覆盖深度天然不同**（V4 短，分类分辨率低），混入会引入 batch 维度（与现有 corpus 的 AMPLICON/WGS 异质性同类，可接受但要带 context/region 标签）。
- 工程动作：枚举 human 相关 context、逐个拉 BIOM、统一 prep 元数据。属可脚本化但繁琐。

### (2) GG2 NB 重注释的唯一 ASV 数量级 + 算力 —— **low-medium（被高估）**
- GG2 论文已实测：Qiita 全库 V4 ASV = **23,113,447**，dedup+QC 后 **21,074,442**（snapshot 2021-12-14，Deblur 2021.09，仅 V4，DEPP 入树）。这是**全 Qiita（30 万+样本、含 private、含全环境）**的去重 ASV 上界。
- 本项目只取**人体 net-new 子集**（约 3-7 万样本）→ 唯一 ASV 远小于 21M（按样本占比粗折 ~百万级，且 deblur ASV 已是去重单元）。
- 算力参照：ResMicroDb 164 万 ASV 用 NB(both) 在本 HPC 90 并发已跑通。百万级 V4 短序列 NB 注释完全在能力内（V4 短 → 更快）。**不是瓶颈**。
- 风险点：V4-90nt/100nt 短 fragment 的 NB genus 置信度低于全长（本项目 MA 用的是全长代表序列）→ family-stall 比例上升。

### (3) 人体筛选可靠性（free-text metadata 脏）—— **high（主风险）**
- Qiita 有 **EMPO**(Earth Microbiome Project Ontology) + `host_scientific_name`/`host_taxid`(9606) + body site(UBERON) 字段，比纯 free-text 好。
- 但填充率不均、跨 study 不一致；redbiom metadata 搜索是 word-stem NLP（`antibiot` 模糊匹配），**不是结构化保证**。人体/非人体边界、体外/培养物、模式动物 spike-in 需逐 study 审。
- 与本项目既有 ResMicroDb/MicroBioMap 同款 free-text 脏问题同级 → 有现成清洗经验，但仍是最易出错环节。

### (4) 与 MicroBioMap 去重（并行任务，防双计）—— **high**
- MicroBioMap keyed = **NCBI BioSample**（SAMN 系）；Qiita keyed = **`studyID.sampleID`**（如 `10317.xxx`，American Gut = study 10317）。两套命名空间不同。
- 必须经 **BioSample/ENA accession 桥接**去重：Qiita 每 study 的 sample 元数据里带 EBI accession（提交 EBI 的 study），用它映射到 SAMN/ERS → 与 MicroBioMap 的 BioSample list 取差集。
- 风险：**未提交 EBI 的 Qiita 样本无 INSDC accession** → 与 MicroBioMap 无共同键 → 无法确定性去重，只能靠 (study × host_subject_id × 测序协议) 启发式，易漏/误。American Gut 在两库都在 → 不做好去重会双计数万行。

### (5) 私有数据不可得 —— **已确认**
- Private = 119,074 samples / 177 studies，**不可匿名下载**（redbiom 只给 public）。GG2 的 30 万+含 private 是 Knight lab 内部权限，外部拿不到。
- → 这部分直接排除，不计入 net-new。

### (6) Qiita→SRA run 映射缺失时如何去重 —— **核心痛点**
- Qiita ASV 表的样本 ID 是 `studyID.sampleID`，**不携带 SRA/ENA run accession**。
- 去重路径：① 优先用 study 的 EBI accession（PRJEB/ERP）+ 样本 EBI 字段 → 拿到 ERR/ERS → 与 corpus `obs['Run']`(1,793,428) 和 MicroBioMap BioSample 取差集；② EBI 字段缺失的 study → 无 INSDC 键 → 退化为启发式或直接当作"可能重复，保守计为净新增需打标 `dedup_uncertain`"。
- ENA 已证实可桥：UCSDMI center 在 ENA 有 256,085 amplicon run（sample_accession=SAMEA系），说明提交 EBI 的部分能映射；问题只在**未提交 EBI + private** 的缺口。

## 4. 与既有管线的契合度（正面）

- deblur ASV 是**精确去噪代表序列**，可下载（`redbiom fetch samples --context ... --output table.biom`），可重训/重注释 → **完美契合本项目 "序列 → GG2 24.09 NB → genus collapse → 8,114 词表" 管线**（同 ResMicroDb/MicroBioMap）。比 GMrepo(无序列,NCBI-tax)、MGnify V1-V5(无 rep-seq, SILVA-OTU) 干净一个量级。
- 词表 0 扩展：V4 ASV 的 genus 必落在 GTDB 8,114 内（GG2 NB 标签空间），只填已有列。
- 多 V 区批次 → 带 `region`/`context`/`Database=Qiita` 标签，下游可筛（与 corpus "越宽越好、下游再筛" 立场一致）。

## 5. Effort 等级裁定

**整体 = high**（个别子项 very_high 风险）：
- 算力/重注释：low-medium（被高估，不是瓶颈）
- context 合并：medium-high（繁琐但可脚本化）
- 人体筛选：high
- MicroBioMap 去重：high（命名空间不一致 + 私有/未提交 EBI 缺键）
- SRA 映射缺失去重：high（核心痛点）
- 私有不可得：确定性排除

**ROI 判断**：净新增可能仅 3-7 万人体行（被快照 + MicroBioMap 重叠压缩），却要付 context 枚举 + 重注释 + 双重去重 + free-text 人体筛的全套成本。建议**排在 MicroBioMap 之后**做：先落地 MicroBioMap（拿到其 BioSample list 作去重锚点），再用 redbiom 只拉 **MicroBioMap 未覆盖 + 后快照 + 提交 EBI 可去重** 的人体 context 子集，避免重复劳动与双计。

## 6. 风险清单（落地前必读）

1. **双计风险（最高）**：American Gut/Microsetta、EMP 人体部分在 MA、MicroBioMap、Qiita 三处都可能出现；命名空间不一致（studyID.sampleID vs SAMN vs SRR/ERR）→ 不做 BioSample 级去重会重复数万行，污染训练分布与 train/val/test split。
2. **去重不可达**：未提交 EBI + private 样本无 INSDC 键 → 无法与 corpus/MicroBioMap 确定性去重；要么排除，要么打 `dedup_uncertain` 标签留下游。
3. **人体边界脏**：EMPO/host_taxid 填充不全 + word-stem 模糊匹配 → 误纳非人体/体外/动物样本。
4. **多 V 区异质性**：V4-90nt 与 V3-V4 等 genus 深度不同 → 引入 region batch effect；必须带 region 标签，short-read NB family-stall 比例升高。
5. **快照重叠浪费**：58% 的 Qiita ENA amplicon run 在 MA 快照窗口内（≤2021-01-04），大概率已在 corpus；不先去重会重注释一堆已有样本。
6. **context 维护成本**：context 名含 pipeline 版本(Deblur 2021.09)，Qiita 升级 pipeline 会新增 context；需固定快照版本号保证可复现。
7. **private 不可得**：119,074 private samples 直接排除（GG2 的 30 万统计含 private，不要据此高估可得量）。

## 7. 一手数字溯源

| 数字 | 值 | 来源 |
|---|---|---|
| Qiita public 16S samples | 388,189 | qiita.ucsd.edu/stats live 2026-06 |
| Qiita public/private samples | 435,633 / 119,074 | 同上 |
| Qiita EBI 提交 samples | 385,491 | 同上 |
| GG2 整合 Qiita V4 ASV | 23,113,447 → dedup 21,074,442 | Greengenes2 PMC10818020；snapshot 2021-12-14；仅 V4；含 public+private |
| Qiita 整合覆盖样本 | "over 300,000 public and private"（含 EMP + American Gut/Microsetta 全部） | 同上 |
| UCSDMI(Qiita) ENA amplicon runs | 256,085 | ENA Portal API live |
| └ ≤2021-01-04（MA 窗口内） | 148,307（58%） | ENA Portal API live |
| └ ≥2021-01-05（后快照 net-new 上界） | 107,778（42%） | ENA Portal API live |
| ENA human-gut amplicon ≥2021-01-05（全源上界） | 331,317 | ENA Portal API live (tax 408170) |
| MicroBioMap 样本 / BioProjects / 截止 | 168,484 / 482 / NCBI Oct 2021 | Cell 2025 / PMC10592789 |
| redbiom 论文测试规模 | "over 100,000 samples" | mSystems 2019 PMC6593222 |
