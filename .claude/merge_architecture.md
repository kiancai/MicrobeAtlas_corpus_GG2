# 语料合并架构 v2（重构：merge-then-filter）

> 状态（2026-06-05）：架构已定，**目录重构已落地**：旧 Merge 已归档到 `scripts/archive/Merge/`，
> 三个库组成模块已归入 `scripts/databases/`，全部路径引用已更新。**Layer2 实跑合并仍待** Qiita/PRIME
> 组成表到位 + 几个待定决策。`scripts/MGnify/` = 覆盖率基准（sample 列表 + 重合分析，**不并入语料**，见 `data/mgnify_evaluation.md`）；
> Qiita 暂停中、脚本暂留 `scripts/databases/Qiita/`（完成后并入 databases/）。

## 为什么重构

- **旧（filter-then-merge）**：各库各自剔零/QC → 筛完才合并。问题：同一 run 在 MA 被 QC 删、
  却从新库（MicroBioMap/Qiita）进来 → **"删了又进"不连贯**；QC 阈值在各库各自 vocab 上算，不统一。
- **新（merge-then-filter）**：各库 → 标准 **genus 组成表（counts）** → **一次性合并 + metadata + 统一筛选**。
  筛选只做一次、在能看见一个 run 所有副本时做。

## ★ 组成表契约（每个库 Layer1 的输出 —— Qiita 也按这个产）

每库产出一张 sample × genus 表，要求：
- **counts（整数）**，不归一、不丢弃、**不做 sample 级筛选/QC**
- **全样本保留**（含零 taxon 样本）
- **taxon 已清理**：去 mitochondria / chloroplast / 非 Bacteria&Archaea
- **var = 6 级路径** `[Domain,Phylum,Class,Order,Family,Genus]`；深层→genus，
  **浅层留 `g__` 空占位桶（不删）**
- 各库各自方法到这步：MA `OTU→genus`；RM/MBM/Qiita `ASV→GG2 24.09 NB(both)→genus collapse`；PRIME 直接是 GG2 genus 表

> **给跑 Qiita 的 claude**：只要落到这个契约即可无缝并入；**不要在组成表阶段做 sample 级筛选或归一化**，
> 也不要删浅层桶。筛选/归一/去重全部留给 Layer2 统一做。

### 三库现状（2026-06-05 只读核对，均符合契约）

| 库 | 组成表 | shape | X | 真 genus + 浅层桶 | obs 主键 |
|---|---|---|---|---|---|
| MA | `gg2.full.h5ad` | 2,690,735 × 7,424 | int32 | 6,410 + 1,014 | `Run.Sample` 复合 ⚠ |
| RM | `resmicrodb.gg2.genus.h5ad` | 100,391 × 5,891 | int32 | 4,954 + 937 | 干净 Run（SRR/ERR/DRR/CRR） |
| MBM | `microbiomap.gg2.genus.h5ad` | 168,464 × 5,899 | int32 | 4,954 + 945 | 干净 Run |

⚠ **MA obs 主键是 `Run.Sample` 复合**（MAP_SID）→ Layer2 concat 时取点前缀得 Run。RM/MBM 已是干净 Run。

## 实际结构（2026-06-05 重构已落地）

```
scripts/
  databases/         # Layer1: 每库 → 组成表 + metadata（已归入）
    MicrobeAtlas/  ResMicroDb/  MicroBioMap/        ← 已移入
    （Qiita 暂留 scripts/databases/Qiita/，完成后并入；PRIME 待建）
  merge_v2/          # Layer2: 统一合并（全新；最终可重命名为 merge/）
    01_concat / 02_merge_metadata / 03_dedup_filter / 04_finalize   （尚未实现）
  archive/Merge/     # 旧 filter-then-merge 09-12（已移入）
  MGnify/            # 覆盖率基准：sample 列表 + 重合分析（不并入语料）
  Phylogeny/ Poincare/ SampleDistance/ SampleDistance100k/   # 库无关资产，暂留顶层
  jxt_scripts/  micoformer/                                   # 外部上游 / 模型码
```

> 各库 sample 级 QC（MA 05/07、RM 05/08）**暂留各库目录** —— 其产物 Stage1/Stage2 仍是开放决策 #3；
> merge-then-filter 下由 Layer2 统一筛选取代，待决策定后再决定去留/归档。

## 规则：如何新增一个语料库（给未来的 claude / Qiita）

1. **做到「组成表契约」**（见上）——用任意适合该库的方法：
   - 有 ASV/OTU 序列 → GG2 24.09 NB(`both`) → genus collapse（复用 `databases/ResMicroDb` 或 `databases/MicroBioMap` 脚本）。
   - 已是 GG2 genus 表（如 PRIME）→ 只需对齐 var 形态。
2. **脚本放** `scripts/databases/<DB>/`（含该库 metadata 标准化脚本，产出 per-sample metadata 表）。
3. **产物落** `results/feature_table/<db>.gg2.genus.h5ad`（counts / 全样本 / `[Domain..Genus]` / 含浅层桶）。
4. **不在组成表阶段做 sample 级筛选 / 归一 / 去重**——全部交 Layer2。obs 主键=Run（非 Run 的库记录映射关系，Layer2 处理）。
5. Layer2 会自动 concat 它；metadata 在 `02_merge_metadata` 加 DB 前缀并入。

> 这就是"各库自成一体、只需达标契约"的全部规则。改不改 Layer2 取决于待定决策，不取决于加了哪个库。

## Layer2 统一合并（全新，`scripts/merge_v2/`）

1. **01_concat** ✅定：各库 genus 表 → 对齐 8,114 → concat 大 counts（全样本，带 `Database`/`Run`；
   MA 取 Run 点前缀）；**浅层桶按样本求和 → `obs['unresolved_reads']`**；主矩阵纯 8,114（varp/varm 干净）。
2. **02_metadata** ✅定（schema 细节⏳）：各库标准化 metadata + GMrepo/mBodyMap 富集 → 按 Run/BioSample join obs（DB 前缀）。
3. **03_flag（原 dedup_filter）** ✅定（2026-06-05）：**不删任何行**；算 N-way 跨库 `Run_paired_id` + 每行 QC 旗标
   `obs[qc_pass / n_deep_reads / n_deep_features]`（旧 v1 阈值：深层 reads≥1000 且 深层属≥5）。去重/筛选下放下游。
4. **04_finalize** ✅定：`layers['counts']`=genus counts；`X`=相对丰度（口径文档化）；挂 `varp` 距离 + `varm` 位置编码；
   `obs[unresolved_reads/total_reads/...]`；`uns` provenance。

## 已锁定的决策

- **merge-then-filter** 架构（统一层筛选，不再 per-DB 筛完才合）。
- **granularity = genus（8,114 GG2）**；组成表存 **counts、全样本、taxon 清理**。
- **unresolved 处理**：深层→genus；**浅层不进 token 矩阵**，按样本求和存 `obs['unresolved_reads']`；
  per-lineage 细节留在各库组成表归档。模型端用 `counts + unresolved_reads` 自行选相对丰度口径：
  - (a) 丢：`genus/genus.sum`（已解析属内部组成）
  - (b) 真分母：`genus/(genus.sum+unresolved_reads)`（占全样本生物量真实比例）
  - (c) 当 UNK：`[genus,unresolved]/total`
- **metadata 在统一层一起合**（含 GMrepo/mBodyMap 富集）。
- 旧 filter-then-merge 代码**归档不删**。
- ✅ **QC = 打标不删（2026-06-05 定）**：Layer2 **一行不删**，QC 降级为 obs 旗标。理由有实测支撑——
  按旧 v1 QC 跨库会删 100 万行(30%)，且 **48,526 个 run 跨库判决分歧**（某库删、另一库留）；
  全留 + 旗标既不误删这 4.8 万，又把"删不删/怎么去重"下放训练端。详 `.claude/data/corpus_sample_overlap.md` +「QC 实测」。
  - obs 存：`qc_pass`(深层 reads≥1000 且 深层属≥5) · `n_deep_reads` · `n_deep_features` · `unresolved_reads` · `Run_paired_id`(跨库副本串联) · `Database`。

## 公共列「冲突」如何处理（2026-06-05 澄清）

全留不删架构下，公共列**不存在 build 时的格子级冲突**：
- **公共列按行填，每行只来自一个母库**（`Database` 标明）→ MA 行的 `Sex` 来自 MA、MBM 行的来自 MBM；**单行内部自洽**，不会"把 A、B 两库的值塞同一格"。
- 同一样本在两库 = **两行**（同 Run 保留两份，`Run_paired_id` 串联）；两行各带各母库的公共值，**可能相同或矛盾**——**合并层不解决**，差异完整保留可审计，下游去重时再面对。
- **库内**多原始列→一公共列（如 `sex`+`host_sex`）：各库已定**优先级**（host_sex 优先），非跨库冲突。
- GMrepo/mBodyMap 富集 = **纯前缀**，更不碰公共列。
- 可选下游工具：标记"配对行公共列是否分歧"（让重复样本两库标注打架处可见）。

## 待定决策（P2/P3 前要敲）

1. ✅ **去重 × QC 已定（2026-06-05）**：**全留不删 + 打标**（见上「QC=打标不删」）。去重粒度=`Run_paired_id`(Run 字面；BioSample 级可下游再精化)。
2. ✅ **QC 阈值已定**：沿用 v1（深层 reads≥1000 且 深层属≥5），但**只作旗标不作删除**。
3. **staged 训练**：per-DB QC 文件（Stage1=MA / Stage2=RM）去留——合并产物自带 `Database`+`qc_pass` 旗标即可复现任意 stage 子集 → 旧 per-DB QC 文件可归档。
4. ✅ **metadata 合并原则已定**：列级「共同含义」才去前缀 → `.claude/data/metadata_merge_design.md`（MA/RM/MBM/PRIME/Qiita ✅，剩 GMrepo/mBodyMap + 其它 metadata-only）。
5. **时机**：五库组成表齐（PRIME Layer1 待建）后一次合。

## 分期

- **P1（✅完成）**：本文档 + `merge_v2/` 骨架 + 只读核对契约。
- **P1-cleanup（✅基本完成）**：旧 Merge → `archive/Merge/`；MA/RM/MBM → `databases/`；全部路径引用已更新。
  剩：Qiita 脚本归位（待其完成）+ Phylogeny/Poincare/SampleDistance 可选归 shared/analysis。
- **P2**：写 Layer2（01/04 可先写；02 等 metadata schema；03 等去重×QC 决策）。
- **P3**：实跑 + 验证（reads 守恒 / N-way 重叠 / QC 漏斗）→ 新语料。
