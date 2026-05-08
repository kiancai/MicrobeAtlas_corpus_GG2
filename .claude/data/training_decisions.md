# 训练数据最终决策（2026-05-08）

## 选定文件

| 角色 | 文件 | shape | 说明 |
|---|---|---|---|
| Stage 1 预训练 | `results/feature_table/gg2.full.qc.h5ad` | 1,762,635 × 6,306 | MicrobeAtlas full→05（删 shallow） |
| Stage 2 二次预训练 | `results/feature_table/resmicrodb.gg2.genus.qc.h5ad` | 93,425 × 4,952 | ResMicroDb 05（删 shallow） |

## 决策一：MicrobeAtlas 用 full 而非 minfilter

`full→05` 是 `minfilter→05` 的**严格超集**（minfilter 子集 0 独有），增量 191,892 样本：
- 93.8% AMPLICON + 1.5% WGS（合法增量），4.2% RNAseq（小头噪声），0.6% NaN
- 21.6% 是 genus < 20 的"低多样性真实样本"（婴儿肠道、感染、极端环境）—— 分布尾部
- 增量样本质量画像反而更好（mean reads 59k vs minfilter 子集 23k）
- minfilter 的 ≥20 OTUs 阈值是为经典 α 多样性分析设计的，与表征学习目标不一致

诊断证据见对话记录；增量构成不是 bug，是 MicrobeAtlas（生态学 QC）vs 我们（表征学习 QC）的哲学差异。

## 决策二：两数据集都用 qc 而非 qc_v2（接受 ~20% reads 损失）

**取舍**：训练阶段丢失 family 级停滞 reads（MicrobeAtlas 20%、ResMicroDb 20%，主要是 d__Bacteria 占位 + family-stall），换取：
1. 每个 var 都是真实 GG2 genus，**语义干净**
2. 与下游经典分析口径一致——别人用模型 embed 的样本表征做 α/β 多样性、丰度比较时不需要先过滤"假 var"
3. 两数据集 var 空间统一，便于 outer join 到 7,930 维合并空间

**自觉接受的代价**：模型不学 family 级占位符信号；Pseudomonas / Haemophilus / Escherichia 这类高 family-stall 菌系在卡停样本中失真（影响 ~14% 高丰度样本，详见 `.claude/data/gg2_vs_silva.md`）。

ResMicroDb 不再补 `05_v2_keep_shallow` 版本。

## 文件体积差异说明（660M vs 18M）

CSR 稀疏存储，体积随 nnz 走：

| | MicrobeAtlas | ResMicroDb | 倍数 |
|---|---:|---:|---:|
| 样本数 | 1.76M | 93k | 18.9× |
| nnz 总数 | 266M | 5.1M | 52× |
| 每样本非零 genus 均值 | 151 | 54 | 2.8× |
| 文件大小 | 660M | 18M | 36.7× |

每样本 genus 数 2.8× 差距来自数据生成方式：MicrobeAtlas 是和 136 万策展全长 16S 参考库比对（高覆盖），ResMicroDb 是项目内 ~2000 ASV 短 V 区段 de novo 聚类（低覆盖天花板）。这是**数据生成方式的天花板，不是质量问题**——stage 2 batch 密度本来就比 stage 1 低，符合预期。
