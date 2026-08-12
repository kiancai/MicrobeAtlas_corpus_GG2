# 关键决策详解

> 从 CLAUDE.md 拆出。CLAUDE.md 只保留"核心决策速览"（一句话版）+ 指针；完整理由在此，更细的实测/对照在 `.claude/data/`。

## 数据库与方向

- **MicrobeAtlas → GG2 NB `auto`**（默认）— 代表序列方向已统一，加 `both` 反而退化（−1.96pp）
- **ResMicroDb → GG2 NB `--p-read-orientation both`** — 方向混合 ASV，`both` 比 `auto` 多救 +6.51pp reads
- 版本要求：QIIME2 ≥ 2025.7 才有 `both`（本项目 `qiime2-2026.1` ✓）

> 完整对照实测、jxt 流水线 7 步、SRP515474 案例 → `.claude/data/direction_orientation.md`

## 训练数据选定

| 角色 | 文件 | shape |
|---|---|---|
| Stage 1 预训练 | `results/feature_table/gg2.full.qc.h5ad` | 1,762,635 × 6,306 |
| Stage 2 二次预训练 | `results/feature_table/resmicrodb.gg2.genus.qc.h5ad` | 93,425 × 4,952 |
| **MiCoFormer 语料库 (含位置编码)** | `results/feature_table/MCFCorpusV2.gg2.h5ad` | 1,826,126 × 8,114 |

- 用 **full** 不用 minfilter（minfilter 的 ≥20 OTUs 阈值与表征学习目标不符）
- 用 **qc** 不用 qc_v2（接受 ~20% reads 损失换语义干净 var 空间）

> full vs minfilter 诊断、qc vs qc_v2 取舍、文件体积差异 → `.claude/data/training_decisions.md`

## 合并立场：语料库越宽越好，下游训练时再筛子集

- 07/08 步只删"不可能进入任何训练分析"的样本（MA: RNAseq 9,415 + Sequencing_Type NaN 20,519，共 29,934 行 / 1.7%）
- 保留：WGS（MA 经流水线提取 16S read 后与 AMPLICON 在 OTU 表语义对齐）/ 阴性对照 / 异常坐标 / 各种宿主体外样本（动物/土壤/水/植物）
- 同 Run 在两库重复 32,698 个 → 合并时两份都保留（不同流水线视角；带 `Database` 标签便于下游再筛）

## 合并 obs schema（54 列）

- **obs_names 全局唯一**：MA 行 `MA_0000000..MA_1732700`，RM 行 `RM_0000000..RM_0093424`。Run 字段作为独立列保留并允许重复
- **公共 9 列不前缀**：`Database / Run / BioSample / Project_ID / Sequencing_Type / Sex / Smoking / Latitude / Longitude`
- **MA 独有 17 列保持 `MA_*` 前缀** / **RM 独有 28 列新加 `RM_*` 前缀**（07/08 步骤产物不带 RM 前缀，前缀仅在 09 合并步引入；保证 RM 单库使用时列名干净）
- 跨库 token 取值集（不映射）：`Sequencing_Type ∈ {AMPLICON, WGS, 16S}` / `Smoking ⊇ {Smoker, Non-smoker, Ex-smoker}` / `BioSample` MA=SRS 系 RM=SAMN 系并存
- 标准化动作：`Sex` MA `female/male` → `Female/Male`（按 RM 大小写）
- 语义重合但分类法不同的列（样本部位 / 健康状态 / 年龄分箱 / 地理）**不映射**，分别保留 `MA_*` 与 `RM_*` 两套

## 数据库选择：GG2 而非 SILVA

- GG2 = GTDB phylogeny，分类粒度更准确
- SILVA 复合属（如 `Escherichia-Shigella`）+ 大菌系过度合并（Bacteroides 6+ clade 打包）会污染样本表征
- GG2 的 family-stall 弱点集中在 Pseudomonas / Haemophilus / Escherichia / Enterobacter（开关式：仅约 14% 样本里它们高丰度）

> 全局加权、ASV 一致性矩阵、bipolar 影响、SILVA 质量问题、用例判断 → `.claude/data/gg2_vs_silva.md`

## Genus 距离矩阵选树和算法

- **用 `phylogeny.id.nwk`**（真实系统发育树，~23M tip，带枝长 + 内部节点带 GG2 rank 标签），不要用 `taxonomy.id.nwk`（只是把分类法层级用 newick 表达，**无枝长**）
- `taxonomy.id.tsv.gz` 不是 GG2 完整 genus 空间——只有 6,757 个 g__（GG2 自己 DEPP 在 V4 ASV tip 级别的保守输出）；**真值是 phylogeny.id.nwk 内部节点的 8,114 个 g__**
- **不要用 skbio `shear` + `tip_tip_distances`**：在 GG2 这种"23M tip + 31.7% 0 枝长内部节点 + 巨型多分叉"结构下，shear 的 unifurcation collapse 丢枝长，50% 非对角 pair 算出 0。改用手写 newick 解析（只内部节点 ~100 万）+ `scipy.csgraph.shortest_path` 分批 Dijkstra（17 批 × 500 sources，~32 min）
- 4 个 -1e-6 浮点噪音负枝长 → clamp 到 0
- 落盘前必须 cross-check：scipy 结果 vs 独立"祖先链遍历"算法在 20 个随机抽样 pair 上必须吻合（diff < 1e-3）；phylo 距离按 taxo hop 分桶必须单调递增

## Taxonomic hop vs phylogenetic distance 的关系

- 内部节点方案的语义：`phylo_dist` 表示 **GG2 定义的 genus clade node 到另一个 genus clade node 的 patristic distance**，不是"某属所有 tip 到另一属所有 tip 的平均距离"。这与最终 genus-level feature table 的粒度一致；若以后要解释 ASV/OTU 层面的属内多样性，则需要另做 tip/代表序列层面的距离。
- `taxo_dist` 是 0-6 的离散层级 hop，解释性强；`phylo_dist` 是连续树距离，可以细分同一个 hop 内部的远近。两者不是重复信息：全 GG2 genus pair 抽样相关性约 Pearson 0.52 / Spearman 0.48；随机比较不同 taxo hop 的 pair 时，phylo 顺序与 taxo 顺序约 82.4% 一致，仍有约 17.6% 会被连续树距离重新排序。
- 全 8,114 genus 上，phylo 距离按 taxo hop 的中位数和 5%-95% 范围：

| taxo hop | 含义 | phylo median | phylo 5%-95% |
|---:|---|---:|---:|
| 1 | 同 family 不同 genus | 33.4 | 6.1 - 119.7 |
| 2 | 同 order 不同 family | 77.9 | 23.8 - 174.5 |
| 3 | 同 class 不同 order | 127.9 | 52.9 - 260.4 |
| 4 | 同 phylum 不同 class | 180.2 | 104.4 - 269.5 |
| 5 | 同 domain 不同 phylum | 258.8 | 148.5 - 407.5 |
| 6 | 跨 domain | 316.3 | 240.7 - 444.4 |

- 关键观察：phylo 均值随 taxo hop 单调增加，说明大方向符合分类层级；但每个 hop 内部离散度很大，说明把所有同层级 pair 压成 1/2/3/4/5/6 会丢掉不少结构。例如 `taxo_dist=1` 的同 family pair，phylo 可从接近 0 到 221.3；`taxo_dist=5` 的跨 phylum pair，phylo 可从 5.7 到 656.3。
- 具体例子：同为 `Flavobacteriaceae` 的 `g__Hyunsoonleella_826360` vs `g__Muricauda_A_821532`，`taxo_dist=1` 但 `phylo_dist=221.345`，远高于很多 hop=2/3 的 pair；跨 phylum 的 `g__SZUA-191` vs `g__Entotheonella`，`taxo_dist=5` 但 `phylo_dist=5.706`，说明分类 rank hop 有时会掩盖树上的近邻关系。
- 结论：保留两套距离有价值。`taxo_dist` 适合作为离散、可解释的分类层级 prior；`phylo_dist` 适合作为连续、可细化同层级内部差异的系统发育 prior。下游建议分别测试 `taxo only` / `phylo only` / `taxo + phylo`。
- 注意：GG2 树存在零长度内部枝，当前 `phylo_dist` 有少量非对角 0 距离（审计为 219 个 pair）。一般不影响使用；若下游算法要求不同 genus 距离严格大于 0，可用 `phylo_dist + 1e-6 * (taxo_dist > 0)` 作 tie-break。
