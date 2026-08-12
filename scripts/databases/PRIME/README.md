# PRIME — Layer1 组成表（Tier1 零重处理）

PRIME（curated 16S meta-analysis 资源）下载物本身**已是 GG2 24.09 genus counts 表**——
实测 5,308 深层真属全部 ⊆ 8,114 GG2 24.09 词表（零版本漂移；`Collapse=Silva and greengenes2`
仅指它同时提供 Silva+GG2 两套，我们用 `gg_genus_*` 的 GG2 套）。故无需重注释，只需 reshape。

## 脚本
| 脚本 | 作用 |
|---|---|
| `01_build_genus_anndata.py` | `rawdata/PRIME/gg_genus_absolute.csv` → `results/feature_table/prime.gg2.genus.h5ad`。分块读 counts + taxon 清理（去 3 列 Mitochondria/Chloroplast/Unassigned）+ 对齐契约 var/obs schema（对标 MBM）。 |

## 产物（符合组成表契约）
`prime.gg2.genus.h5ad` = **53,449 × 6,793**（5,308 真属 + 1,485 浅层 `g__` 桶）；
X=int32 counts、全样本（含 928 零 taxon）、obs 主键=Run（PRIME `Sample` 列，DRR/SRR/ERR）。

## metadata
`samples.csv`（84 列）+ `projects.csv` 的列设计见 `.claude/data/metadata_merge_design.md`「PRIME」段
（两伪等价陷阱：`Sequencing_Type`=读长配置→`PRIME_Read_Config`、`Host_Age` 多分桶→仅数值进 Age）。
metadata 标准化脚本待 Layer2（02_merge_metadata）阶段写。
