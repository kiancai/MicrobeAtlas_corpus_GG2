# merge_v2 — 统一合并层（Layer2，重构后）

> **全新目录。** 把"各库 → genus 组成表"统一 concat → metadata → finalize（QC=打标不删，无独立去重步）。
> 取代旧 `scripts/archive/Merge/`（旧 09-12，已归档）。
>
> 完整架构、组成表契约、已锁定决策 → **`.claude/merge_architecture.md`**

## 脚本

| 脚本 | 作用 | 状态 |
|---|---|---|
| `00_overlap_audit.py` | 五库 Run 重叠 + ∩MGnify + GMrepo/mBodyMap 覆盖审计 | ✅ 已跑（→ `data/corpus_sample_overlap.md`） |
| `01_concat_compositions.py`(+`.sbatch`) | 五库 → 对齐 8,114 → vstack 大 counts；浅层桶求和→`unresolved_reads`；QC 旗标 + `Run_paired_id`（全留不删；MA 取 Run 点前缀，Qiita 经 ENA 桥） | ✅ **已跑** → `results/feature_table/merged.counts.h5ad` **3,334,541 × 8,114**（qc_pass 70%，跨库副本 693,701 行） |
| `02_merge_metadata.py` | 各库标准化 metadata（公共列+DB 前缀）+ GMrepo/mBodyMap 纯前缀富集 → join obs | 设计已定（待写） |
| `04_finalize_corpus.py` | X 相对丰度 + varp 距离 + varm 位置编码 + provenance；**删 var 旧统计列并重算** | 设计已定（待写） |

> ⚠ `merged.counts.h5ad` 的 var 暂带 5 个旧 MCFCorpusV2 统计列（prevalence/mass_fraction… 在旧 1.8M 上算的，已过期）→ 04 删除/重算。
> 原计划的 `03_dedup_filter` **已废**：QC=打标不删，去重旗标(`Run_paired_id`)已在 01 内完成，筛选/去重下放下游。

## 输入（组成表契约，见 merge_architecture.md）
- counts、全样本、taxon 清理、`[Domain..Genus]`、浅层留 `g__` 占位桶、obs 主键=Run
- **五库齐**：`results/feature_table/{gg2.full, resmicrodb.gg2.genus, microbiomap.gg2.genus, prime.gg2.genus, qiita.gg2.genus}.h5ad`

## 当前阶段
P2 进行中：01_concat ✅ 已产 `merged.counts.h5ad`。下一步 = 各库 metadata 标准化脚本 + 02_merge_metadata → 04_finalize。
