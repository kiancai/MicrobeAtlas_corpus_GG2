# ResMicroDb 流水线产物校验（2026-05-07，both 版本）

## 流水线漏斗

| 阶段 | 输出文件 | shape | reads | %02 | 本步丢 |
|---|---|---|---:|---:|---:|
| 02 sample × ASV | `resmicrodb.gg2.asv.h5ad` | 100,391 × 1,639,277 | 6,633,559,930 | 100.00% | – |
| 03 聚合到 genus_var | `resmicrodb.gg2.genus.h5ad` | 100,391 × 5,891 | 6,559,887,108 | **98.89%** | 1.11% |
| 04 drop empty | `resmicrodb.gg2.genus.nonzero.h5ad` | 100,342 × 5,891 | 6,559,887,108 | 98.89% | 0.00%（仅删 49 个零样本） |
| 05 QC | `resmicrodb.gg2.genus.qc.h5ad` | 93,425 × 4,952 | 5,238,029,532 | **78.96%** | 20.15% |

## both vs auto 收益

| 指标 | auto（旧） | both（当前） |
|---|---:|---:|
| 03 reads 保留率 | 92.38% | **98.89% (+6.51pp)** |
| SRP515474 reads-to-genus | 7.10% | **84.41% (+77pp)** |
| 32 个原"质量异常"项目 reads-to-genus 中位数 | < 50% | **61.73%** |
| 全数据集 reads-to-genus < 10% 的项目 | 32 | **3 (0.75%)** |

## 02 防御性检查（脚本里加的）
1. 跳过 QIIME2 export `taxonomy.tsv` 第二行可能的 `#q2:types` 元数据行
2. 三方 ASV ID 一致性 assert（`asv.fa` ∩ `taxonomy_gg2.txt` ∩ `taxonomy_silva.txt`）—— 全 398 项目通过

## 05 内部细分（拆开 04→05 的 20.15% 损失）
| 子步骤 | reads 丢失 | %04 |
|---|---:|---:|
| step 1：删 shallow var (`g__` 占位符) | 1,317,744,319 | **20.09%** |
| step 2：迭代 min_reads/min_features 阈值 | 4,113,257 | 0.06% |

→ **99.7% 的损失来自 step 1**（删 GG2 family-级停滞的 reads），跟样本质量阈值几乎无关。

## 04→05 损失按 GG2 停滞层级
| GG2 停滞层级 | ASV 数 | reads | %02 |
|---|---:|---:|---:|
| Genus（保留）| 896,034 | 5.24 G | 79.02% |
| Family | 202,212 | 517 M | 7.80% |
| Order | 36,272 | 67 M | 1.01% |
| Class | 33,615 | 43 M | 0.65% |
| Phylum | 3,762 | 0.8 M | 0.01% |
| **仅到 Domain** (d__Bacteria) | **423,117** | **696 M** | **10.48%** |
| GG2 Unassigned | 44,265 | 68 M | 1.02% |

注意"仅到 Domain"（10.48%）比 Family 级停滞（7.80%）还大——是 GG2 的最大弱项。
