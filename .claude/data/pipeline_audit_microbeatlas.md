# MicrobeAtlas 流水线产物校验（2026-05-06）

## 03 步 OTU→Genus 聚合产物

### 输出
- `results/feature_table/gg2.full.h5ad` — sample × var = `2,690,735 × 7,424`，873 MB
- `results/feature_table/gg2.minfilter.h5ad` — `1,884,129 × 7,424`，749 MB
- 存储格式：`X` 是 `scipy.sparse.csr_matrix`（CSR，只存非零），h5ad 再 gzip 压缩

### OTU 数对账
- `taxonomy.tsv`：111,870（来源 = `otus.97.allinfo` 第 7 列代表序列经 GG2 NB 注释）
- 03 过滤掉 1,835 条：83 非 B/A + 901 mitochondria + 851 chloroplast
- 过滤后 110,035 条 → 与 BIOM 的 103,271 OTU 取交集 = **101,507**
- 8,528 条 tax 不在 BIOM：参考库有但 269 万样本中从未检出，正常
- 1,764 条 BIOM 不在 tax：基本就是被 03 过滤掉的 mito/chloro/non-BA
- 聚合 key（QIIME2 风格 6 级路径 `d__;p__;c__;o__;f__;g__`）唯一值 = **7,424**

### 75% reads 丢弃来源（重要：看着吓人但实际正常）

03 报告 `丢弃 reads: 75.55%`，分解后：

| 类别 | reads 占比 | OTU 数 |
|---|---|---|
| **`Unmapped` 伪 OTU** | **75.15%** | 1 |
| chloroplast | 0.28% | 836 |
| mitochondria | 0.12% | 847 |
| 非 B/A | 0.00% | 80 |

`Unmapped` 是 MicrobeAtlas 在 BIOM 里塞的"未比对到任何参考序列的 reads 计数桶"——单独一行，172.5 G reads。它本来就不是生物特征，不在 `taxonomy.tsv` 里，03 自然把它剔除了。**真正的生物信号丢失只有 ~0.4%**。如果将来要对外报告"过滤损失"，建议把 `Unmapped` 单独列出，避免给读者"75% 生物信号被丢了"的误解。

per-sample 丢失中位数 76% 也是同一原因——MicrobeAtlas 数据本身就有大量未比对 reads，与 03 过滤无关。

### 零 taxon 样本（→ 04 步处理）

03 输出后**未做样本级过滤**，需要后续脚本剔除：

| 数据集 | 样本总数 | 零 taxon | 占比 |
|---|---|---|---|
| `gg2.full.h5ad` | 2,690,735 | **310,231** | 11.53% |
| `gg2.minfilter.h5ad` | 1,884,129 | 0 | 0.00% |

零 taxon = 该样本所有 reads 都落在 `Unmapped` 或被 03 剔除的条目上。`gg2.minfilter` 因 MicrobeAtlas 自己已经做过最小过滤（每样本至少 13 reads / 1 var），不需要这一步。

### 检查代码位置
后续若要复现这套对账，思路：
1. 把 `taxonomy.tsv` 按 mito/chloro/非 B/A 分组，用 OTU id 在 BIOM 中查 `observation/matrix` 各行 reads 总和
2. `BIOM_OTU - tax_OTU` 差集第一个就是 `Unmapped`
3. AnnData 行 nnz 用 `np.diff(X.indptr)`，行和用 `cumsum(data)[indptr]` 求差（避免对每行循环）

---

## MA 数据快照时间 = 2021-01 初（2026-06-04 验证）

本地 MA 是**一次性快照爬取**，截止约 2021 年 1 月初——决定了"净新增样本只能来自 2021-01 之后首次公开的资源"（语料扩增的核心前提，见 `explorations/2026-06-04_dataset_expansion.md`）。

**方法**：对 MA 全量 run（2,690,735）取各前缀最大编号，向 ENA Portal API 查 `first_public`。

- 前缀计数：SRR 1,852,203 / ERR 766,634 / DRR 71,898。
- **最大 SRR = `SRR13347130`，`first_public = 2021-01-04`**；185 万 SRR **无一 > SRR15M**（硬天花板，分布到 13.35M 戛然而止）。
- 最大 ERR ≈ `ERR5013241`；`ERR5000000` first_public = 2021-01-05。
- **坑**：run 编号→日期非严格单调（`SRR13000000` first_public 反而 2023-06-05，embargo 个例）；但结论靠"MA 含的最大编号 + 硬天花板"，不受影响——若爬于 2023 必含大量 SRR15M–23M 的 run，而它一个都没有。
- **澄清**：上游 "MicrobeAtlas 2026 (Cell)" 更新版（2.39M 样本）未下载；本地是 2.69M 早期快照。
- 覆盖率佐证（各库 run ∩ MA，随发布年单调下降）：mBodyMap(2021) 96% → PRIME(2025) 49% → GMrepo v3(2026) 27%。
- **净新增审计基线** = 最终 corpus `obs['Run']` **1,793,428** unique run（= 1,826,126 行 − 32,698 跨库 Run 重复对；可从 `MCFCorpusV2.gg2.h5ad` 的 `obs['Run']` 提取）。
