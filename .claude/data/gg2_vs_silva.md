# GG2 vs SILVA 深度对比（2026-05-07 实测）

## 全局 reads 加权
| 类别 | GG2 (both) | SILVA |
|---|---:|---:|
| to_genus | **79.02%** | **88.46%** |
| shallow | 19.95% | 7.17% |
| unassigned | 1.02% | 4.37% |

## ASV 级一致性矩阵（百万级）
```
              SILVA→
GG2↓        shallow  to_genus  unassigned
shallow     217,553  249,216   232,209
to_genus     38,187  850,992     6,855  ← 主对角线 85 万 ASV 双方都到 genus
unassigned    9,117   22,576    12,572
```

## Top 30 (collapsed view) 重叠度
- GG2 4,955 raw genus → 剥 GTDB 后缀后 4,281（**仅 1.16× 膨胀**）
- 拆得最厉害：Clostridium 39 子类（但只占 0.02% reads，几乎没影响）；Pseudomonas 16 子类（1.60%）
- **Top 30 collapsed 与 SILVA top 30 重叠 25/30**，量级一致（Streptococcus 14.82% vs 15.66%）

## Genus 命名差异（两套都到 genus 的 85 万 ASV）
- 名字字面一致：**44.18%**
- 主要差异类型：
  - GTDB 加字母后缀：`Stenotrophomonas_A_615274` vs `Stenotrophomonas`
  - GTDB 加数字后缀：`Bifidobacterium_388775` vs `Bifidobacterium`
  - GTDB 重新分类：`g__Desulfonema_C` vs `Sva0081_sediment_group`

## 核心洞察：GG2 的 family-stall 在 SILVA 下是哪些菌

按 SILVA top 50 反查每个 SILVA genus 的 ASV 在 GG2 下到什么级别：

| SILVA genus | reads（全） | GG2 卡停率 | GG2 退到 |
|---|---:|---:|---|
| Streptococcus / Prevotella / Veillonella / Dolosigranulum | 数十亿 | < 6% | 几乎全部到 genus ✓ |
| Staphylococcus / Neisseria / Rothia / Porphyromonas | 数千万–数亿 | 11-15% | 大部分到 genus ✓ |
| **Haemophilus** | 3.47 亿 | **60.1%** | f__Pasteurellaceae |
| **Pseudomonas** | 2.30 亿 | **53.8%** | f__Pseudomonadaceae |
| **Escherichia-Shigella** | 4,604 万 | **96.6%** | f__Enterobacteriaceae |
| **Enterobacter** | 1,178 万 | **100.0%** | f__Enterobacteriaceae |
| Granulicatella | 3,602 万 | 38.2% | – |
| Atopobium | 1,381 万 | 68.6% | – |

→ Pseudomonas / Haemophilus / Escherichia / Enterobacter 这几个**肠杆菌科 + 巴氏杆菌科**菌系是 GG2 的"软肋"，因为 GTDB 对它们做了 major reclassification（拆成多个 _A/_B/_E 子分支），NB 分类器在 ASV 边缘 case 上判不清子分支就退到 family。

---

## 关键发现：相对丰度的"bipolar"影响（2026-05-07）

**绝对 reads 卡停率 ≠ 相对丰度偏差**——核心在于这些菌只在少数样本里高丰度。

### Haemophilus 样本级别（按 SILVA 视角分层）

| SILVA 相对丰度区间 | 样本数 | SILVA rel% | GG2 rel% | 偏差 |
|---|---:|---:|---:|---:|
| < 1% | 72,286 | 0.01% | 0.00% | +0.00pp |
| 1–5% | 13,513 | 2.17% | 1.55% | +0.06pp |
| 5–10% | 4,436 | 6.89% | 5.81% | +0.22pp |
| **10–30%** | **4,655** | **16.34%** | **2.30%** | **+10.57pp** |
| **30–60%** | **2,809** | **41.09%** | **0.00%** | **+35.89pp** |
| **60–100%** | **2,692** | **83.41%** | **1.03%** | **+64.70pp** |

### 三个重灾区菌的高丰度样本数（≥5% rel）

| 菌 | ≥5% 样本 | ≥10% | ≥30% | 高丰度样本 GG2 中位 vs SILVA 中位 |
|---|---:|---:|---:|---|
| Haemophilus | 14,592 (14.5%) | 10,156 | 5,501 | **2.55% vs 18.70%（−85.5%）** |
| Pseudomonas | 11,574 (11.5%) | 8,896 | 5,346 | **7.36% vs 26.43%（−27.0%）** |
| Escherichia-Shigella | 3,164 (3.2%) | 2,325 | 1,081 | **0.00% vs 16.95%（−100%）** |

→ **86% 样本里这些菌本来就 < 5%，影响 0**。**~14% 样本里它们高丰度，GG2 严重低估**（特别是 Escherichia-Shigella 几乎完全丢失）。

### 这意味着 GG2 影响是"开关式"的
- 健康人鼻咽 / 口腔 / 皮肤 / 阴道 / 健康肠道 → **GG2 完全可用**（主导菌都能到 genus）
- 上呼吸道感染（Haemophilus 主导） → ⚠ GG2 显著低估
- 囊性纤维化 / VAP（Pseudomonas 主导） → ⚠⚠ GG2 显著低估
- 婴儿/感染肠道（Escherichia 主导） → ⚠⚠⚠ GG2 几乎全卡

---

## SILVA 数据库自身的质量问题（2026-05-07 实测）

### SILVA "to_genus" 88.46% 的水分构成
| 类别 | reads | %02 | 性质 |
|---|---:|---:|---|
| **real_genus**（真实有意义） | 5.74 G | **86.55%** | 有效 |
| composite（连字符复合属） | 67.5 M | 1.02% | ⚠ 拼盘 |
| group / clade | 13.6 M | 0.21% | ⚠ 占位 |
| uncultured | 26.3 M | 0.40% | ⚠ catch-all |
| cluster_id（UCG-X / RC9 等） | 19.2 M | 0.29% | ⚠ 编号 |
| **小计 SILVA "to_genus"** | 5.87 G | 88.46% | – |

→ **真正可用 86.55%，含糊命名 1.92%**。**严格口径 GG2 vs SILVA 差距 = 7.52pp，不是 9.44pp**。

### SILVA 复合属（虚假 genus）实例
```
Escherichia-Shigella                                4,604 万 reads
Burkholderia-Caballeronia-Paraburkholderia         1,327 万 (3 合 1)
Methylobacterium-Methylorubrum                       503 万
Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium   261 万 (4 合 1)
```

→ SILVA 把 2-4 个真实 genus 强行拼成一个名字。**Escherichia-Shigella 是最严重的——你研究里若想区分大肠杆菌 vs 志贺氏菌，SILVA 给不了答案。GG2 反而拆成 g__Escherichia / g__Shigella（虽然部分卡 family）**。

### SILVA 群组占位（top 实例）
```
uncultured                                       2,631 万 (catch-all 最大)
TM7x                                                665 万
Clostridia_UCG-014                                  404 万
[Eubacterium]_nodatum_group                         273 万
Lachnospiraceae_NK4A136_group                       175 万
Rikenellaceae_RC9_gut_group                          99.7 万
[Ruminococcus]_torques_group                         60 万
```

→ 这些都不是真 genus：方括号 `[Eubacterium]` 表示"放在 Eubacterium 但不一定是"；`_NK4A136_group` / `_UCG-014` 是研究编号，不是分类阶元。

### SILVA "过度合并"问题
- SILVA 把 1 万多 Bacteroides ASV 都标 1 个 genus
- GTDB 把这些拆成 6+ 个**系统发育上独立**的 genus 类群（Bacteroides_E / _F / _G / _H 等）
- SILVA 的"高覆盖"部分是用过度合并换来的——把不互通的 clade 打包

→ SILVA 用 LPSN（传统命名法）受历史命名约束；GTDB 基于全基因组系统发育重组。**SILVA 在 Bacteroides / Pseudomonas / Lachnospiraceae 等大菌系上把不互通的 clade 当成同一个**。

### "SILVA 独到"（GG2 不到 genus）的命名质量
SILVA 独到 27.2 万 ASV / 8.09 亿 reads 中：
- 90.5% 是真 genus（有效信号）
- 5.8% 是复合属
- 2.3% uncultured
- 1.4% group / cluster_id

→ "SILVA 独到"大部分确实是真有效，GG2 的 family-stall 是真实损失，不能完全甩锅给 SILVA。但 SILVA 自身的 catch-all 也吃掉了它一部分"漂亮覆盖率"。

---

## 用例驱动的最终判断（2026-05-07）

### 核心结论：**GG2 是正确选择**，特别是对基础模型 / sample 表征训练

1. **GG2 是 GTDB phylogeny**——基于全基因组的系统发育，分类粒度更准确
2. **SILVA 在你的语境下不是中性 alternative**：复合属 + 过度合并会**主动污染样本表征**
3. **GG2 family-stall 不是 deal-breaker**：在 03 输出层这部分 reads 仍是合法占位符 var，承载 family-级正确信号
4. **跨数据集一致性**：MicrobeAtlas + ResMicroDb 都用 GG2，可联合分析；用 SILVA 就丢了 MicrobeAtlas 兼容性

### 但 05 的 shallow 删除对 LLM/基础模型训练**过度激进**

| 阶段 | 形态 | 适用场景 |
|---|---|---|
| 03 / 04 | 含 family 占位符 var，**98.89% reads 保留** | **基础模型 / sample 表征训练（推荐）** |
| 05 (qc) | 仅真 genus var，78.96% reads | 经典统计分析（α/β diversity / abundance 比较） |
| 05_v2_keep_shallow | 保留占位符 + 阈值过滤 | LLM + sample QC 折中（建议为 ResMicroDb 也加这个版本） |

→ **对 LLM 训练**：用 04 输出（或仿 MicrobeAtlas 写一份 05_v2_keep_shallow）。05 默认版本会丢 20% reads，这部分在 GTDB phylogeny 里有 family-级真值，不该简单删除。

### 真正影响重大的样本子集（GG2 不可信的场景）
- 上呼吸道感染（Haemophilus 主导）
- 囊性纤维化 / VAP（Pseudomonas 主导）
- 婴儿肠道感染 / IBD / NICU（Escherichia / Klebsiella 主导）

→ 对这些场景，可在 02 输出 (`resmicrodb.gg2.asv.h5ad`) 的 silva_* 列做 fallback 统计，但跨样本仍以 GG2 为主体特征。

### 不建议反向（SILVA → GG2）的根本原因
即使 SILVA 能"救" Haemophilus / Pseudomonas / Escherichia 这几个菌的 reads，但 SILVA 把 Bacteroides / Lachnospiraceae 等更大菌系做"过度合并"会**破坏大量样本的 phylogeny 信号**——影响面比 GG2 family-stall 更广，只是更隐蔽。

---

## 跨数据集合并潜力（MicrobeAtlas + ResMicroDb，2026-05-07）

| 维度 | 数 |
|---|---:|
| MicrobeAtlas (gg2.full) var | 7,424 |
| ResMicroDb (qc) var | 4,952 |
| **交集（var_id 完全一致）** | **4,446 (89.8% 的 ResMicroDb)** |
| ResMicroDb 独有 | 506（多是 Archaea + 项目特化菌） |
| MicrobeAtlas 独有 | 2,978（环境多样性更广） |
| **合并后 var 总数** | 7,930 |

→ 90% var 共享，可直接 outer join。两数据集合计 **~280 万样本**（MicrobeAtlas 269 万 + ResMicroDb 10 万）的 sample × genus 矩阵，是基础模型训练的优质底盘。
