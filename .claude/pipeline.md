# 流水线脚本与产物速查

> 从 CLAUDE.md 拆出，集中记录"已建成的流水线脚本"与"各步产物漏斗"。更细的审计/schema 见 `.claude/data/`。

## 流水线脚本一览

```
scripts/databases/MicrobeAtlas/        # OTU → genus 聚合（输入是策展的全长代表序列）
  01_extract_rep_seqs.sh           从 otus.97.allinfo 提取 97% 代表序列
  02_qiime2_classify.sbatch        GG2 NB 分类（用默认 auto，不要用 both）
  03_build_genus_anndata.py        过滤 mito/chloro/非 BA + 6 级路径聚合
  04_drop_empty_samples.py         剔除零 taxon 样本
  05_qc_filter.py                  删 shallow var + 迭代 QC（MIN_READS=1000, MIN_FEATURES=5）
  05_v2_qc_filter_keep_shallow.py  保留 shallow var 的版本
  06_attach_metadata.py            samples.env.info.tsv → obs 26 列 metadata（py:percent ↔ ipynb）
  07_filter_samples.py             删 Sequencing_Type ∈ {RNAseq, NaN}（保留 AMPLICON + WGS）

scripts/databases/ResMicroDb/          # 398 项目 ASV → 跨项目 sample × genus
  01_qiime2_classify.sbatch        单项目 GG2 NB 分类（已加 --p-read-orientation both）
  01_run_loop.sh                   轮询提交器（CAP=90, INTERVAL=30s, 断点续投）
  02_merge_to_asv_anndata.py       398 项目 → sample × ASV CSR（GG2 + SILVA 双注释）
  03_build_genus_anndata.py        过滤 mito/chloro/非 BA + 6 级路径聚合
  04_drop_empty_samples.py         剔除零 taxon 样本（含项目维度报告）
  05_qc_filter.py                  同 MicrobeAtlas（MIN_FEATURES=5）
  06_standardize_metadata.py       metadata_all.txt → metadata_all.standardized.tsv（135,746 × 36；UTF-16LE → UTF-8；jxt 派生 Age_Group/Case_Or_Control/Is_Healthy + projectTable region）
  06b_export_patches.R             从 ps.16s_0105_new7.rds 导出 Run 级 patch tsv（仅供 Patch 1/2 消费，baseR + phyloseq）
  06b_fix_metadata_errors.py       v2 NCBI 真值方案：5 项目 1,697 Run 修正（P1 Sample_Site 1,257 + P2 PRJNA801796 Phenotype 细分 255 + P3 PRJNA822681 NCBI host disease 152 + P4 PRJNA824137 NCBI title 33）→ .standardized.fixed.parquet
  07_attach_metadata.py            metadata.fixed.parquet → anndata.obs 37 列（Database + Run 列 + 35 metadata；Run 既作 obs_names 又作列）
  08_filter_samples.py             passthrough 骨架（anndata 已 100% 16S；保留 Negative Control 等留给下游训练时筛）

scripts/databases/MicroBioMap/         # Human Microbiome Compendium → GG2 24.09 重注释（流程严格复用 RM；用户监督一致性）
  01_shard_fasta.sh                366万全局 ASV fasta 切 31 片（120k/片）
  02_classify_shard.sbatch         GG2 24.09 NB(both) array 分类（≡RM 01；classifier/env/参数逐字相同；NB逐序列独立→分片与逐项目等价）
  03_merge_taxonomy.py             31 片 taxonomy → 全局 taxonomy_gg2.global.tsv（校验 366万行/唯一/==fasta ASV集）
  04_build_asv_anndata.py          482真项目表(跳过1个._macOS垃圾)+全局taxonomy → sample×ASV（≡RM 02；同parse_gg2；zip→tar BytesIO流式不落盘；SRS→SRR；float count断言整性→int32；ASV全局共享故var['project']恒='MicroBioMap'）
  05_build_genus_anndata.py        ASV→genus collapse（RM 03 逐字逻辑：Domain∈{B,A}+去mito/chloro+6级聚合；浅层留g__占位；全样本保留不剔）
  run_build_composition.sbatch     驱动 03→04→05（输出存在则跳过，断点续跑；256G/8h）
  README.md                        MBM↔RM 步骤映射表
  （待商议，未做）剔零样本=RM 04 / QC=RM 05 / metadata / 并进语料

scripts/databases/Qiita/               # Qiita 整包 99G → GG2 24.09 重注释（B=全局去重；复用 RM/MBM）；详 .claude/data/qiita.md
  assess/{01_dedup_assess,02_inventory_biom,03_profile}.py  探索：净新增三方去重 / Deblur表盘点 / study×方法+metadata+能注释总量
  00a_select_artifacts.py    每 16S 样本选 1 张表（90nt ref-hit 优先；混合 study 用 prep Jaccard 锁 16S）→ 321,502 样本 / 1,259 表
  00b_build_per_study.py     按 study 出 work/studies/<study>/{asv.fa, asv.h5ad}（685 study）
  01_global_asv.py           全局去重（md5 头）→ asv_global.fa = 16,737,797 唯一 ASV（去重率 34.8%）
  02_shard_fasta.sh          切 67 片（250k/片）
  03_classify_shard.sbatch   GG2 full-length NB(both) array 分类（≡RM/MBM）；⚠ 含长 ASV(150/250nt) 片需 ≥256G（本次 5/67 片 128G OOM→400G 救回）
  04_merge_taxonomy.py       67 片 → taxonomy_gg2.global.tsv（md5→Taxon，校验 ==asv_global.fa）
  05_build_genus.py          per-study collapse（RM03/MBM05 逐字逻辑）；按 md5(seq) 查全局 taxonomy → <study>/genus.h5ad
  06_concat_genus.py         outer-join var + vstack → qiita.gg2.genus.h5ad
  run_finalize.sbatch        SLURM afterok:<分类array> 串 04→05→06（256G/计算节点）
  （待商议，未做）sample 级筛选 / metadata / 并进语料 —— 交 Layer2

scripts/archive/Merge/               # 两库合并 + 挂距离矩阵 + corpus finalize
  09_merge_databases.py            outer-join var；obs_names 重编号 MA_/RM_ 前缀；Sex 标准化；RM 独有列加 RM_ 前缀；从 var_names 重建 6 级 taxonomy 列
  10_expand_and_attach_phylo.py    anndata var 从 6,857 扩到 GG2 24.09 全 8,114（补 1,257 全 0 列 + observed bool 列）；挂 varp['taxo_dist'] + varp['phylo_dist']
  11_finalize_corpus.py(+.sbatch)  X→相对丰度 (float32) + 原 counts 备份到 layers['counts']；obs +3 列 (total_reads/n_taxa/Run_paired_id)；var +4 列 (n_samples_observed/prevalence/mass_fraction/mean_rel_abundance_when_present)；uns['provenance']
  12_attach_poincare.py            euclidean_d32.npz → varm['position_encoding'] (8114×32 float32)；uns['provenance']['position_encoding'] 记录 method/dim/scale/pearson；输出 MCFCorpusV2.gg2.h5ad（gzip 压缩）

scripts/Poincare/            # 8,114 GG2 genus → 欧氏 R^D 位置编码（基于 patristic 距离矩阵）
  02_euclidean_embed.py            普通 Adam + L2 距离 + 绝对 MSE loss；d_tree 归一化到 [0,10]；FP32 训练；CPU < 15s；输出 euclidean_d{D}.npz
                                   （早期 Sarkar 双曲嵌入 Pearson 0.77 / geoopt Poincaré ball Pearson 0.98 都已废弃：
                                    下游"MLP → token embedding 相加"管线天然欧氏，双曲位置编码 mismatch 大）

scripts/Phylogeny/           # 生成 GG2 24.09 全 genus 距离矩阵（独立于 anndata，输出可发布）
  01_extract_genus_vocab.py            从 phylogeny.id.nwk 解析全 8,114 个 g__ 节点 + walk-up 拼完整 6 级 path
  02_compute_taxonomic_distance.py     分类法层级 hop 距离 (8114² int8, 取值 {0..6})；向量化 rank-by-rank 比对
  03_compute_phylogenetic_distance.py  patristic 距离 (8114² float32)；手写 newick 解析（只保留内部节点）+ scipy.csgraph.shortest_path 分批 Dijkstra；不依赖 skbio

scripts/SampleDistance/      # 50k 分层抽样 → BC + weighted UniFrac → PCoA + 可视化
  01_stratified_sample.py        MA 30k + RM 20k 两层 sqrt 分层（seed=42）
  02_build_subset_anndata.py     主表切子集（含 varp 搬运）
  03_collapse_tree.py            phylogeny.id.nwk 折叠到 8,114-tip genus_tree.nwk（嵌套 g__ 插 0 枝长 self-tip）；20 对 patristic cross-check
  04_compute_bc.py(+.sbatch)     50k² Bray-Curtis（64 核 ~12 min；fork Pool + scipy cdist）
  05_compute_wunifrac.py(+.sbatch) 50k² Striped Fast UniFrac（64 核 < 1 min）
  06_pcoa.py(+.sbatch)           truncated eigsh top-10；输出 obsm + eigenvalues.tsv
  07_plot_pcoa.py                4 张 PCoA 散点图（bucket / scree / human_sites / rm_sites）

scripts/SampleDistance100k/  # 独立 100k 分支：corpus-级 PCoA 可视化（不改主表，不替代 50k）
  01_stratified_sample_100k.py     10k paired Runs（双库各 1 份）+ 70k MA + 10k RM = 80k MA / 20k RM（seed=42）
                                   MA 加额配额：Human 22k / Animal 12k / Soil 13k / Aquatic 12k / Plant 6k / Unknown 5k
  02_build_subset_anndata_100k.py  主表切 100k 子集（含 varp 搬运）
  03_compute_bc_100k.py(+.sbatch)  100k² Bray-Curtis float16
  04_compute_wunifrac_100k.py(+.sbatch) 100k² Striped Fast UniFrac float16
  05_pcoa_100k.py(+.sbatch)        truncated eigsh top-10 + 负特征值审计
  06_plot_pcoa_100k.py(+.sbatch)   4×2 主图：MA env / MA human sites / RM sites / paired Runs × {BC, wUniFrac}
                                   产物在 results/sample_distance_100k/，与 50k 完全分离

scripts/jxt_scripts/         # 上游 ASV 与 SILVA 注释来源（外部，不改）
```

## 流水线产物速查

### Qiita (sample × genus, Layer1 完整组成表, 2026-06-05)
- 06 → **`results/feature_table/qiita.gg2.genus.h5ad` 321,502 × 6,954**（counts；总 134 亿 reads；真 genus 5,745 + 浅层桶 1,209；零行 13；全样本不筛不去重，待 Layer2）
- 中间件 `rawdata/Qiita/work/`：`asv_global.fa`(1670万)、`shards/`(67)、`tax/`(67片) + `taxonomy_gg2.global.tsv`、`studies/<study>/{asv.fa,asv.h5ad,genus.h5ad}`、`sample_artifact_map.tsv`
- 去重资产（Layer2 用）：`rawdata/Qiita/qiita_ena_bridge.tsv`(372,512 run)、`qiita_netnew_16s.tsv`

### MicrobeAtlas (sample × genus)
- 03 → `gg2.full.h5ad` 2,690,735 × 7,424（保留 98.89% biological reads）
- 04 → `gg2.full.nonzero.h5ad` 2,380,504 × 7,424（删 11.53% 零 taxon 样本）
- 05 → `gg2.full.qc.h5ad` 1,762,635 × 6,306（**Stage 1 训练用**）
- 06 → `gg2.full.qc.with_meta.h5ad` 1,762,635 × 6,306（X 不变，obs 补 26 列 metadata）
- 07 → `gg2.full.qc.with_meta.filtered.h5ad` 1,732,701 × 6,306（删 RNAseq + NaN 共 29,934 行 / -1.70%）

> 03 步 reads 漏斗、Unmapped 75% 解释、零 taxon 样本来源 → `.claude/data/pipeline_audit_microbeatlas.md`
> 06 obs 26 列 schema、col2/col3 解析、缺失值约定、跨库 merge 策略、anndata StringArray 写盘陷阱 → `.claude/data/obs_metadata_schema.md`

### ResMicroDb (sample × genus)
- 02 → `resmicrodb.gg2.asv.h5ad` 100,391 × 1,639,277（带 GG2 + SILVA 双注释）
- 03 → `resmicrodb.gg2.genus.h5ad` 100,391 × 5,891（98.89% reads）
- 04 → `resmicrodb.gg2.genus.nonzero.h5ad` 100,342 × 5,891
- 05 → `resmicrodb.gg2.genus.qc.h5ad` 93,425 × 4,952（**Stage 2 训练用**）
- 06 → `metadata_all.standardized.{tsv,parquet}` 135,746 × 36（all 全表标准化版）
- 06b → `metadata_all.standardized.fixed.parquet` 135,746 × 36（v2 NCBI 真值方案：1,697 Run 修正 / 5 列）
- 07 → `resmicrodb.gg2.genus.qc.with_meta.h5ad` 93,425 × 4,952（X 不变，obs 补 37 列 = Database + Run + 35 metadata，消费 fixed parquet）
- 08 → `resmicrodb.gg2.genus.qc.with_meta.filtered.h5ad` 93,425 × 4,952（passthrough；anndata 已 100% 16S，row filter 无动作）

> 02-05 实测漏斗、both vs auto 收益、05 损失分解 → `.claude/data/pipeline_audit_resmicrodb.md`
> 06 输入四件套（metadata_all UTF-16LE / sampleTable_changed / sampleTable_v5 / projectTable）、36 列 schema、派生列逻辑（Age_Group jxt 7 档 + Case_Or_Control + Is_Healthy）、Region_16S 填充规则、与 jxt 对照 → `.claude/data/metadata_resmicrodb_standardize.md`
> 06b v2 NCBI 真值方案（4 个 patch / 共 1,697 行）：Patch 1 Sample_Site Nasal→Nasopharynx 1,257（0105 源） / Patch 2 PRJNA801796 Influenza 细分 + Phenotype_ID 补齐 255（0105 源） / Patch 3 PRJNA822681 NCBI host disease 152（76 COVID-19↔Health 翻转 + 76 反向，比 0105 多覆 39 行） / Patch 4 PRJNA824137 NCBI title 33（13 HC + 3 TBZ + 17 LTBI=Latent Tuberculosis Infection/MONDO_0040753/case，14 TBM 已对）→ `rawdata/ResMicroDb/supplement data/CHANGES_0105_new7.md`

### MicroBioMap (Human Microbiome Compendium → sample × genus)

源 `rawdata/MicroBioMap/15122187.zip`（v1.1.1，CRC 校验完整）：168,464 人肠道 16S 样本 / **482 真项目** / **3,663,370 全局 ASV**。净新增 vs corpus `obs['Run']`：**61,990 / 168,464 (36.8%)** net-new（其余 63.2% 已在 MA）。

- 01-02 → 366万 ASV 全局 GG2 24.09 NB(both) 分类（31片 array；genus 解析率 ASV计数 57.7% / **reads加权 71.57%**）
- 03 → `work/taxonomy_gg2.global.tsv` 3,663,370 行（Feature ID/Taxon/Confidence；Confidence 中位 0.98）
- 04 → `microbiomap.gg2.asv.h5ad` **168,464 × 3,663,370**（CSR int32，reads 守恒逐项目校验；226 MB）
- 05 → **`microbiomap.gg2.genus.h5ad` 168,464 × 5,899**（**交付物=完整组成表**；var=[Domain..Genus] 同 MA/RM schema；obs主键=SRR；X int32 nnz 1,010万 / 12.49B reads；reads守恒 97.89%，丢 mito/chloro/非BA 共 2.11%；**真genus列 4,954**≈RM 4,952 + 浅层占位 945；零taxon样本 14 **保留**）
  - **到此为止**：sample级筛选(剔零=RM04 / QC=RM05)、metadata、并进语料 = **待商议，未做**
- 详细：`scripts/databases/MicroBioMap/README.md`（MBM↔RM 映射 + 一致性要点）

### 跨数据集合并

中间步：
- 09 → `merged.gg2.h5ad` **1,826,126 × 6,857**（obs 54 列；reads 守恒 47.68 B；X int32 CSR，density 2.1%；size 849 MB）
- 10 → `merged.gg2.with_phylo.h5ad` **1,826,126 × 8,114**（X 前 6,857 列 nnz/dtype 不变 + 新增 1,257 全 0 列；var 增 `observed` bool 列；`varp['taxo_dist']` + `varp['phylo_dist']`；obs 54 列不变）；size 1.05 GB

最终产物：
- 11 → **`MCFCorpus.gg2.h5ad`** 1,826,126 × 8,114（**X 改成 relative_abundance float32 CSR**；原 int32 counts 备份到 `layers['counts']`；obs 57 列 = 原 54 + `total_reads`(int64) + `n_taxa`(int32) + `Run_paired_id`(string)；var 11 列 = 原 7 + `n_samples_observed` + `prevalence` + `mass_fraction` + `mean_rel_abundance_when_present`；`uns['provenance']` 元数据；varp 保留）；size 1.80 GB
  - **下游用法**：`adata.X` 直接拿到 row_sum=1 的相对丰度；想要原始 reads 用 `adata.layers['counts']`
  - `obs['Run_paired_id']`：跨库 Run 重复对的 group id（32,698 对 × 2 = 65,396 行非空）；下游分 train/val/test 时 group split 防泄漏
- 12 → **`MCFCorpusV2.gg2.h5ad`** 1,826,126 × 8,114（V1 + `varm['position_encoding']` 8114×32 float32）；size 1.80 GB
  - 欧氏位置编码：从 `varp['phylo_dist']` 出发，普通 Adam + L2 距离 + 绝对 MSE loss，dim=32，d_tree 归一化到 [0,10]
  - 验证：Pearson(d_eucl, d_tree) = **0.9924**，MAD = 0.0506，坐标范数 median=2.63 / max=6.41
  - `uns['provenance']['position_encoding']` 记录 method/dim/scale/pearson/mad/npz 路径
  - 下游用法：`adata.varm['position_encoding']` 直接喂 MLP → 与 token embedding 相加；想还原 patristic 距离用 `d_eucl / scale`

合并语义（不随 11 改变）：
- var 来源：共享 4,401 + 仅 MA 1,905 + 仅 RM 551
- Run 重复 65,396 行 = 32,698 个 Run × 2 副本（两库流水线视角并存）
- obs 公共 9 列 + `MA_*` 17 + `RM_*` 28（11 又叠 3 列）
- 06b 已落地：下游 07/08/09 已基于 fixed metadata 重跑，shape 与 reads 不变（仅 obs 取值修正）

### Genus 距离矩阵（GG2 24.09 全量，独立副产物 + 挂回 anndata）

副产物（`results/phylogeny/`，可独立发布）：
- 01 → `genus_vocab.tsv` 8,114 行 × 6 列（GG2 24.09 全 genus + 完整 6 级 rank path）
- 02 → `genus_taxo_dist.npz` 8114² int8，取值 {0..6}；362 KB gzip
- 03 → `genus_phylo_dist.npz` 8114² float32 patristic 距离；209 MB

挂回 anndata：见上节 10 步产物（`varp['taxo_dist']` + `varp['phylo_dist']`），11 步沿用不动。

三个 8,114 vs 6,857 数字的关系：
- **8,114** = GG2 24.09 完整 genus 空间（`phylogeny.id.nwk` 内部节点 g__ 计数 = NB classifier 训练标签空间）
- **6,757** = `taxonomy.id.tsv.gz` 唯一 g__（GG2 自己 DEPP 给每 tip 分类，V4-only 信息不够，保守）
- **6,857** = anndata var.Genus（你的 OTU 全长 16S → NB → 输出，跨 6,757 但只覆盖 8,114 的 84.5%）；anndata 独有但 TSV 缺的 900 个 = NB 比 DEPP 大胆，**这就是 vocab 必须从树取不能从 TSV 取的原因**

### Sample × sample 距离矩阵（50k 分层子集）

`MCFCorpus.gg2.h5ad` 1.83M 全样本上算 pairwise 不现实（1.83M² 矩阵存不下也算不动），改在 50k 分层子集上算。产物全在 `results/sample_distance/`，**主 anndata 不变**。（注：子集是 11 finalize 之前从 10 步 `merged.gg2.with_phylo.h5ad` 切出来的，X 仍是 counts；BC/UniFrac 计算内部自己做相对丰度归一化，与最终 corpus 的 X 语义无冲突。）

抽样：MA 30k（Human 8k / Animal_other 8.1k / Soil 5.2k / Aquatic 5.0k / Plant 2.2k / Unknown 1.5k，sqrt 加权）+ RM 20k（按 RM_Sample_Site × Project_ID）。跨库 Run 重复不去重。`RANDOM_SEED = 42`。

- 02 → `subset_50k.h5ad` 50,000 × 8,114（X + obs(56 列：原 54 + stratum_id + sub_stratum) + varp 搬主表的 taxo/phylo_dist）
- 03 → `genus_tree.nwk` 8,114-tip folded tree；嵌套 g__ 插 0 枝长 self-tip（GG2 有 1,357 个 g__ 是另一个 g__ 的 GTDB 子属后代）；新树 patristic 与 `genus_phylo_dist.npz` 20 对 cross-check < 1e-3 通过
- 04 → `obsp['distance_bc']` 50k² **float16**；relative abundance 后 scipy cdist
- 05 → `obsp['distance_wunifrac']` 50k² **float16**；Striped Fast UniFrac (McDonald 2018 / Sfiligoi 2022)；BIOM v2.1 临时文件 + tree 双输入；64 核 < 1 min
- 06 → `obsm['X_pcoa_bc'/'X_pcoa_wunifrac']` 50k × 10 + `pcoa_eigenvalues.tsv`；truncated eigsh 上 dense float32 B；trace(B) 作分母算 explained var ratio
- 07 → `figures/pcoa_{bucket,scree,human_sites,rm_sites}.png` 4 张可视化

PCoA top-3 累计 explained variance：BC 20.8%（散得开）/ wUniFrac 47.3%（phylo 把方差压到少数轴）。

不做 rarefaction（McMurdie & Holmes 2014）；relative abundance 即可。

最终文件 `subset_50k.h5ad` 6.3 GB（X 340 MB + 两个 50k² float16 距离 ~5 GB + varp 310 MB + obs/obsm 30 MB）。

> 抽样配额 sqrt 推导、tree 折叠算法细节、嵌套 g__ self-tip 推导、PCoA 各部位偏移解释 → `.claude/data/sample_distance.md`

### Sample × sample 距离矩阵（100k 独立分支）

`scripts/SampleDistance100k/` 是 corpus-级 PCoA 可视化的独立分支，产物落在 `results/sample_distance_100k/`，**与 50k 完全分离，不替代也不修改 50k**。设计目标是给出大样本量下的 PCoA 总览图，覆盖更多 paired Runs 与各 MA bucket。

抽样设计 100,000 行：
- **10,000 cross-database paired Runs**：在两库 Run 交集里按 RM_Sample_Site × Project_ID 双层 sqrt 分层选 10k 个 Run，每个 Run 同时入 MA 和 RM 两份（共 20k 行）—— 用于跨库 paired 可视化
- **70,000 MA 加额**：按 6 个 bucket 配额 Human 22k / Animal_other 12k / Soil 13k / Aquatic 12k / Plant 6k / Unknown 5k，每个 bucket 内部按 body site / env subcategory 做 sqrt 分层
- **10,000 RM 加额**：按 RM_Sample_Site × Project_ID 双层 sqrt 分层
- 最终库分布 80k MA + 20k RM，`RANDOM_SEED=42`

产物：与 50k 同结构（subset h5ad + genus_tree.nwk + BC/wUniFrac float16 距离 + PCoA top-10 + 负特征值审计 + figures），主图是 `figures/pcoa_100k_4x2.png` —— 4 行（MA env / MA human sites / RM sites / paired Runs）× 2 列（BC / wUniFrac），每面板都画 100k 灰色背景，仅 highlight 层切换，坐标在同度量内可直接比较。

子分支详细 README 在 `scripts/SampleDistance100k/README.md`。
