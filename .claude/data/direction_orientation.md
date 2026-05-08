# 方向矫正与 GG2 NB read-orientation 决策

## jxt 流水线（ASV 与 SILVA 注释的来源）

完整代码在 `scripts/jxt_scripts/`，摘要 7 步：

1. 双端 merge（`vsearch fastq_mergepairs`）
2. 切引物 + 质控（`vsearch fastx_filter --fastq_stripleft/right 23 --fastq_maxee_rate 0.01`）
3. **方向矫正**（`vsearch --orient -db silva_16s_v123.udb`）— **关键步骤，问题源头**：
   matched 序列翻为正向，**notmatched 序列保留原方向**
4. 去冗余（`vsearch derep --minuniquesize 10`）
5. 去噪（`usearch unoise3`）→ `asv.fa`
6. 生成特征表（`vsearch usearch_global --id 0.97`）→ `otutab.txt`
7. SILVA 分类（QIIME2 2024.2 `classify-sklearn` + `silva-138-99-nb-classifier.qza`）→ `taxonomy_silva.txt`

**重要认知**：`otutab.txt` 是 vsearch_global @0.97 的 reads 归属，**不是精确 ASV 计数**。
对 genus 级聚合分析无影响，但若做 ASV 级精细分析需注意语义。

## GG2 注释的方向问题

**问题**：jxt 第 3 步用 SILVA v123 做 orient。某些项目（典型如 SRP515474）的引物 / 区段
v123 完全不识别 → 全部走 notmatched 分支 → ASV 集合方向混合。
GG2 NB classifier 默认 `--p-read-orientation auto` 只对前 100 条 ASV 选一向锁定全 batch，
混合方向项目里另一半 ASV 报 Unassigned。

**实测对照**（SRP515474, 73,749 ASV, 5.5 亿 reads）：

| 模式 | reads 到 genus | reads Unassigned | 平均 Confidence |
|---|---:|---:|---:|
| `auto` | 7.10% | 87.00% | 0.78 |
| `reverse-complement` | 7.10% (= auto) | 87.00% | 0.78 |
| **`both`** | **84.41%** | **1.09%** | **0.94** |

**ResMicroDb 解决方案**：分类 sbatch 加 `--p-read-orientation both`：每条 ASV 独立尝试
正反两方向，取置信度高的一次。计算时间约翻倍；对原本方向已齐的 ASV 无影响。

**版本要求**：QIIME2 ≥ 2025.7（commit `be2c2df` 引入 `both`）。本项目用 `qiime2-2026.1` ✓。
jxt 当年用 2024.2 没这选项，所以 SRP515474 类项目 GG2 注释失败但 SILVA 正常——
推测 SILVA NB 训练时做了双向数据增强，对方向不敏感。

## MicrobeAtlas 不要用 both（重要更正，2026-05-07 实测）

代表序列是策展库，方向已统一（保守锚点检测：91.33% 正向占优，仅 0.01% / 7 条反向）。
在已统一方向数据上用 `both` 反而有害——RC 方向给出"高数值低意义"的 NB Confidence，
覆盖正向浅注释。实测 4,166 个 OTU 退化，reads 到 genus 覆盖率 88.30% → 86.35%（−1.96%）。

## 结论

- **MicrobeAtlas**：用 **`auto`**（默认），不要加 `--p-read-orientation`
- **ResMicroDb**：用 **`both`**，必要

两数据集分别用各自最优参数。`scripts/MicrobeAtlas/02_qiime2_classify.sbatch` 顶部已注释
说明这条决策。

## 跨项目合并的不变量（ResMicroDb）

- 样本 ID（ERR / SRR / DRR）跨项目唯一（已在 02 校验）
- ASV ID 跨项目重名 → namespace 为 `<PROJECT>__ASV_N`
- 不同项目的 V 区段不同（V1-V3 / V3-V4 / V4 / V1-V9 都有），跨项目 ASV 序列不可直接比较；
  在 GG2 6 级 genus 路径聚合后大致可比
