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
