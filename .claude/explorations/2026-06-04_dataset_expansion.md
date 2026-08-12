# 2026-06-04 数据库扩增调研

> ⚠ **前提部分已更新（2026-06-05）**：当时"下游只用于人体 / 不扩环境"的约束**已作废**——改为**全样本**（人体+动物+环境全收）、**能注释全注释、去重交统一层(Layer2)**。本文的候选库评估结论（PRIME/MicroBioMap/Qiita/MGnify/GMrepo… 各自定性）仍有效，只是"人体 only"那条范围限定不再适用。见 `2026-06-05_qiita.md` + `merge_architecture.md`。

## 触发
用户在 `rawdata/` 新加了 PRIME + GMrepo 两个数据集，问：(1) 能不能直接利用 / 重新注释；(2) 还有没有别的公共库能扩增语料、大概能扩多少；(3) 哪些库能补 metadata；(4) 特殊任务（时序/同人）数据。

## 简介
用多智能体工作流检索 ~92 次命中 → 77 个唯一候选库 → 深核 45 + 完整性批判 3。把候选分三类评估：**直接合并(GG2 24.09 表) / 可重注释(提供 ASV·OTU 序列) / 仅元数据**。再用用户的硬约束收敛：**下游只用于人体、不重跑 raw reads、规模上千起**。

---

## 核心结论

### 1. 各候选库定性（决策相关）

| 库 | 角色 | 规模 | 关键事实 |
|---|---|---|---|
| **PRIME** | ✅ 直接合并 (Tier1) | 53,449 / 16S | GG2 **24.09** genus 表，列名=完整6级路径，与我们词表逐字匹配 5,308/6,796，**0 新属**；整数 counts；Sample=Run；∩MA=26,227(49%) |
| **MicroBioMap** | 🚩 主力净新增行 | 168k 人肠道 | Cell 2025；Zenodo 全量含 **ASV 序列**+每项目 ASV 表（R 包只给 SILVA genus 表）；拿序列跑我们 24.09 NB；人肠道 only；自带 GG2 是 2022.10（别直接用，用序列重注释） |
| **Qiita/redbiom** | 🚩 第二主力 | 10万+ 人体 | Deblur ASV 表 feature=序列，V4 多 GG2-native；redbiom 拉表不碰 raw；**保 host_subject_id+collection_timestamp**（特殊任务关键） |
| **GMrepo** | 仅元数据 | v3 87k 16S | NCBI/SILVA 无序列；陈卫华组**肠道**专库；∩MA=32,392(27%) |
| **mBodyMap** | 仅元数据 | 63k / 22 部位 | GMrepo 姊妹库，同组同 schema，**全身**（Vagina/Oral/Skin/呼吸道…）；∩MA=60,398(96%)；与 GMrepo run **0 重叠** |
| MicrobiomeHD | 小·可顺带 | ~5-6k | 16S case/control + OTU 序列 |
| MGnify | ✗（同日深核确认） | 人体 amplicon 22.3万 | V1–V5 全是 SILVA-OTU 表+reads、**无 rep-seqs**；**v6(ASV) 仍 beta、公开 API 0 条**（采样 3000 分析 V6=0）；ENA 系重叠 MA。版本拆分+V5 下载+API 配方+触发重评 → `.claude/data/mgnify_evaluation.md` |
| IMNGS/环境/中国库/受控/小集/参考库/纯WGS | ✗ | — | 见下"不做"理由 |

### 2. 两个关键洞察

- **MA 覆盖率 = 时间快照效应（同日核实：实测停在 2021-01 初）**：2026-06-04 用 ENA `first_public` 钉死——MA 最新公开 run = `SRR13347130`（2021-01-04），185 万 SRR 无一 > SRR15M（硬天花板）→ 单次快照，截止 ≈2021-01。覆盖率随发布年单调下降佐证：mBodyMap(2021) **96%** → PRIME(2025) **49%** → GMrepo v3(2026) **27%**。所以净新增集中在 **2021 之后的新库** → 主力锁 MicroBioMap(2024-25)/最新 Qiita 是对的。方法+边界细节（含 SRR13000000=2023 embargo 个例、2026 Cell 版未下载）→ `.claude/data/pipeline_audit_microbeatlas.md`「MA 数据快照时间」。
- **词表(genera) ≈ 0 扩展**：8,114 已是 GTDB 全属空间；新样本只能填充已有列，几乎不新增列（真正新 lineage 连 GG2 NB 标签都没有）。所以"扩数据集"本质 = 更多行 + 更密覆盖 + 更厚 metadata，不是分类广度。

### 3. metadata 源盘点（结论：16S 已覆盖）
- **能贴到 16S 语料行**（逐样本 SRA-keyed）：GMrepo + mBodyMap + PRIME = 基本穷尽，无重要新增。
- **WGS 模态**（贴不到我们 16S 行，留 WGS 路线）：HumanMetagenomeDB（NAR 2021，69,822）、curatedMetagenomicData。
- **微生物↔疾病关联库**（无逐样本号，属下游评测）：Disbiome / MicroPhenoDB / BugSigDB / Peryton / GMMAD / Amadis / HMDAD。
- **最全的未用源 = NCBI BioSample/SRA 属性本身**（host_age/sex/disease/body_site），可对全部现有样本回填——是 curation 动作不是"加库"。

### 4. 方法论结论（怎么"彻底确定能加多少"）
不要枚举"数据库"（它们多是同一 INSDC 池的重叠视图）。改为：
- **总体一 = 在 INSDC 的数据**：MA≈其覆盖；净新增 = INSDC 16S 总量 − MA，用一次普查求差集回答（穷尽）。
- **总体二 = 不在 INSDC 的数据**（唯一真·净新增）：一份**有限可封口清单**——中国本土档案 / 受控访问 / 不镜像 INSDC 的直投平台 / 补充材料-only / 晚于 MA 快照的新释放。
- **三个目标分开**：① 更多样本行（INSDC 已基本吃满，靠 2021 后新库）② 更厚 metadata（finite 元库，已盘点）③ 更均衡（审计自身语料找洞，人体下游则环境无需扩）。

### 5. 不做（已评估）
环境扩展（人体下游 + MA 已够）· raw-reads-only（中国 NMDC/CNSA/CNGBdb、GGMP、TwinsUK、LogMPIE、SKIOME，用户不重跑 raw）· 受控访问（TEDDY/FINRISK）· MGnify · 几百样本小集 · 纯参考库(HOMD/WoL2/GSR-DB) · 纯 WGS。

---

## 实测数字（本地已算）
- PRIME genus 列 6,796：命中词表 5,308 + 未解析 bin 1,488（全空 genus，0 新属）。
- Run 重叠：PRIME∩MA=26,227 · PRIME∩RM=1,073 · PRIME∩GMrepo=215；GMrepo∩MA=32,392 · GMrepo∩RM=41。
- （2026-06-04 重算 vs 最终 corpus `obs['Run']`=1,793,428）GMrepo 带可下载丰度的 16S run 41,730 中 **29,455(70.6%) 不在 corpus = net-new**，但 NCBI-taxonomy keyed + 无序列 → **unmappable 到 GG2 8,114-genus 空间**（"0 净新增行"措辞修正为"有净新增但分类法不兼容"，否决不变；详见 gmrepo doc）。
- mBodyMap：63,148 行；BodySite 22 类（Vagina 16,835 最大）；Amplicon 48,747+Meta 14,401；∩MA=60,398 · ∩RM=10,557 · ∩GMrepo=**0**。
- MA 唯一 run 列表已生成于 `/tmp/ma_runs.txt`（269 万，prefix-before-dot of MAP_SID）。

## 产物与详细信息指针
- GMrepo/mBodyMap metadata 详剖 + 即用即合规格（含 GMrepo net-new 重算）→ `.claude/data/metadata_sources_gmrepo_mbodymap.md`
- MGnify 深核（人体 amplicon 版本拆分 / V5 下载 / v6=0 / v1·v2 API 配方 / 触发重评）→ `.claude/data/mgnify_evaluation.md`
- MA 数据快照时间验证（2021-01-04 + 方法）→ `.claude/data/pipeline_audit_microbeatlas.md`「MA 数据快照时间」
- 工作流原始 45 篇 dossier（每库 16S/WGS 量、序列可下否、GG2/GTDB、重叠、净新增估计）在本 session 的 workflow transcript（`microbiome-dataset-expansion-scan`）。
- 整合顺序与状态 → `.claude/roadmap.md`

## 待办（落到 roadmap）
先测 `MicroBioMap ∩ MA` 与 `Qiita ∩ MA` 精确净新增 → 扩样本行（MicroBioMap→Qiita）→ 合并 PRIME → 统一 metadata（GMrepo+mBodyMap+BioSample 回填）→ 特殊任务数据。
