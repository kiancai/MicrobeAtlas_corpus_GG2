# MicroBioMap (Human Microbiome Compendium) → GG2 24.09 重注释

> 把 MicroBioMap v1.1.1（Zenodo `15122187`，168,464 人肠道 16S 样本 / 3,663,370 全局 ASV）
> **重新注释到 GG2 24.09**，产出与 MA/ResMicroDb 同口径的 sample × genus 物种组成表。
>
> **核心承诺：重注释流程与 ResMicroDb 严格一致**（用户当年严格监督过 RM）。下游
> genus collapse 逻辑是 RM 脚本的逐字复用；分类配方逐字相同。

## 当前交付物（仅到这里）

**`results/feature_table/microbiomap.gg2.genus.h5ad`** —— GG2 24.09 重注释的**完整物种
组成表**：168,464 样本（SRR 为 obs 主键）× genus（6 级路径 var）。**全样本保留，不做
sample 级筛选(QC/剔零)，不 padding 到 8,114 词表。**

> sample 级筛选、metadata 合并、并进 MCFCorpus —— **均待与用户商议后再做**（见
> `.claude/roadmap.md`）。用户原话："合并的事情都不要着急"。

## 步骤与 ResMicroDb 对应

| MicroBioMap 脚本 | 作用 | ≡ ResMicroDb |
|---|---|---|
| `01_shard_fasta.sh` | 全局 ASV fasta(366万) 切 31 片 | （RM 无需：已有 per-project asv.fa） |
| `02_classify_shard.sbatch` | GG2 24.09 NB **classify-sklearn `--p-read-orientation both`** array 分类 | `01_qiime2_classify.sbatch`（同 classifier / 同 env / 同参数） |
| `03_merge_taxonomy.py` | 31 片 taxonomy → 全局 `taxonomy_gg2.global.tsv` | （RM 无需：per-project taxonomy 已就地） |
| `04_build_asv_anndata.py` | 项目表 + 全局 taxonomy → sample × ASV（同 `parse_gg2` / 同 var schema） | `02_merge_to_asv_anndata.py` |
| `05_build_genus_anndata.py` | ASV → genus collapse（Domain∈{B,A}+去mito/chloro+6级聚合） | `03_build_genus_anndata.py`（**逐字逻辑**） |
| `run_build_composition.sbatch` | 驱动 03→04→05 一键出组成表 | — |
| （待商议）剔零样本 | sample 级 | `04_drop_empty_samples.py` |
| （待商议）QC 筛选 | sample 级(MIN_READS=1000/MIN_FEATURES=5/删shallow) | `05_qc_filter.py` |

## 一致性要点

- **分类 sharding-invariant**：NB classify-sklearn 逐序列独立，"全局分片"与"逐项目分类"
  对每条 ASV 给出完全相同标签 → 与 RM 等价。
- **唯一结构差异（不影响数据）**：MBM 的 ASV 是全局共享 id（同序列跨项目同号、只分类
  一次），故不做 `PROJECT__ASV` 命名空间化，`var['project']` 设常量 `'MicroBioMap'`；
  此列仅用于诊断打印，不进入丰度/聚合。05 的"按项目丢失率"诊断因此只显示一组。
- **counts**：项目表是 float 文本（整数值），04 读入断言整性后 `astype(int32)`。
- **样本主键**：项目表列名是 SRS，04 按 `sample_metadata.tsv` 映射成 SRR（= corpus Run 键）。

## 运行

```bash
# 1) 分类（已提交）：array 0-30，每片 ≤4h，~31 并发
sbatch -c 8 --mem=128G -t 04:00:00 --array=0-30 scripts/databases/MicroBioMap/02_classify_shard.sbatch

# 2) 分类全部完成后，一键出组成表
sbatch -c 8 --mem=256G -t 08:00:00 scripts/databases/MicroBioMap/run_build_composition.sbatch
```

## 数据来源文件（zip 内，不全量解压）

`rawdata/MicroBioMap/15122187.zip`（775MB，CRC 校验完整）。04 用 `zipfile`+`tarfile`
**流式逐项目读** `project_asv_tables.tar.gz`，不落盘。关键成员：
- `obs_md.txt.gz` `#OTUID→sequence`（366万 ASV 序列，→ `work/asv_global.fa`）
- `project_asv_tables.tar.gz` 483 个 `asv/PRJ*.txt`（行=ASV全局id, 列=SRS, 值=count）
- `sample_metadata.tsv`（srs/srr/project/region/iso/pubdate…）
- `classifications.tsv`/`taxonomic_table.csv.gz`（**它自带 SILVA-138.2 / GG2-2022.10-OTU 注释，不用**）
