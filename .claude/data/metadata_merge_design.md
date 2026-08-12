# metadata 合并设计：列级「共同含义」原则 + 逐库逐列清单

> 状态：**原则已定（2026-06-05，用户定）**；逐库逐列分类**进行中**。
> 这是 Layer2 `merge_v2/02_merge_metadata` 的依据。逐列结论填在下方各库小节。

## 核心原则

**默认带数据库前缀 `<DB>_`；只有当一列的「含义在所有库中可证明完全相同」时，才去前缀、合并成公共列。**

- 判断在**列级语义**：`age=6` 在哪个库都是"6 岁"，含义零差异 → 可公共。
- `body_site` / `phenotype` 这类各库**各自独立标准化**、看似相近实则不同口径 → 强行合并 = **伪等价** → **必须带前缀，不能合**。
- **举证责任在"去前缀"那一方**；**自底向上逐列测**，不预设 canonical schema。

## 判定细则

1. **公共的两种合法来源**：
   - (a) **原始就同义**：Run / BioSample / Project_ID / Latitude / Longitude / 日期 / Library_Layout / Platform / Instrument / 数值 Age …
   - (b) **从公共原始量中央派生**：用同一套规则从公共原始量算（如 `Age_Group` 用统一分桶从公共数值 `Age` 算）。**绝不合并各库 pre-derived 版本**（各库分桶/定义可能不同）。
2. **本体 ID 列**是唯一"跨库合并可辩护"的特例（如 `disease→MONDO id`，共享外部词表，id 含义客观相同）——但须核对各库**用同一本体、同口径标注**；不确定则前缀。free-text 一律前缀。
3. **必带前缀**：body_site / sample_type / phenotype(文本) / 各库自己的 age_group / 各库自己的 is_healthy·case_control / 各库自己的 region 推断 / 环境标签 / 机构 / 自由文本 …… 凡"独立标准化"产物。

## 公共列候选集（逐库过完才 finalize）

**已确立（MA + RM + MBM + PRIME + Qiita 逐列核验，2026-06-05）**：
`Database` · `Run` · `BioSample` · `Project_ID` · `Latitude` · `Longitude` · `Sex`{Male/Female} · `Sequencing_Type`{Amplicon/Shotgun} · `Smoking`{Smoker/Non-smoker/Ex-smoker/Unknown} · `Library_Layout`{PAIRED/SINGLE} · `Platform` · `Instrument` · `Age`(+start/end) · `Age_Group`(7-bin 中央派生) · `BMI`(+start/end) · `Country` · `Continent` · `PMID` · `Subject_ID`(DB 命名空间) · `Region_16S` · `Collection_Date`(真实采样日期；3 peer 升公共) · `Height`(cm) · `Weight`(kg)（2 peer 升公共）

> **2026-06-05 两项"有 peer 即升公共"**（Qiita 触发,用户拍板）：`Collection_Date`（PRIME `Collection_Date`/Qiita `collection_timestamp`/MBM tags `collection_date`；仅真实采样日期,ENA/SRA 行政日期丢）；`Height`/`Weight`（PRIME `Host_Height`/`Host_Weight` + Qiita `host_height`/`host_weight`，客观数值,带单位换算）。

> obs_names = 库原生稳定 id，加 DB 前缀保证全局唯一（沿用旧 merge `MA_`/`RM_` 前缀做法）。
> **Subject_ID / 任何"id 槽"列**：含义同义（"本库的 subject 分组"）但**值必须加 DB 命名空间**（`RM:100BAL1`），防跨库 id 碰撞造成伪同人。
> **公共-派生**：能从公共原始量中央派生的（如 `Age_Group ← Age`），用统一规则中央算，**不合并各库 pre-derived 版本**。

---

## 逐库逐列清单

### MA（26 列）— ✅ 已过（2026-06-05）

| 列 | 含义 | 建议 | 理由/待议 |
|---|---|---|---|
| Database | 库标签 | **公共** | provenance 常量 |
| Run | SRA run | **公共(a)** | 同义 |
| BioSample | SRA biosample | **公共(a)** | 同义 |
| Project_ID | BioProject | **公共(a)** | "研究"标识同义 |
| MA_Sample_ID | MAP_SID 原生 id | → obs_names(DB前缀) | 原生 id 即索引，不单列 |
| Latitude / Longitude | 坐标 | **公共(a)** | 数值同义 |
| Sex | female/male | **公共** ✅ | 统一 Male/Female；核验各库无其它类 |
| Sequencing_Type | AMPLICON/WGS/RNAseq | **公共** ✅ `{Amplicon, Shotgun}`(+RNAseq) | 用户定：只合"扩增 vs 鸟枪"层，**RM 对标 MA**(16S→Amplicon, Metagenomics→Shotgun)，**不分前缀**。靶基因 16S-vs-SSU 细分不入公共列。RNAseq 留待 Layer2 筛 |
| Smoking | Smoker(419)+NA | **公共** ✅ 以 RM 为准 | 词表 `{Smoker, Non-smoker, Ex-smoker, Unknown}`；**实测 MA 无一被填 Non-smoker(count=0)**；MA 映射：Smoker→Smoker / **其余→Unknown**；🔴 严禁空填 Non-smoker |
| MA_Env_{Animal,Soil,Aquatic,Plant}(+_Sub) ×8 | MA 环境标签体系 | **MA_ 前缀** | MA 独有 |
| MA_IsHuman | MA 严格 human 定义 | **MA_ 前缀** | 各库定义不同 |
| MA_SampleSite_Raw | col3 原值 | **MA_ 前缀** | MA 特有 |
| MA_SampleSite (8 部位) | MA body-site 标准化 | **MA_ 前缀** | ★各库独立标准化，必前缀 |
| MA_Health | MA 健康 token | **MA_ 前缀** | MA 自己的健康标注 |
| MA_AgeGroup (5 档) | MA 年龄分桶 | **MA_ 前缀** | 分桶 ≠ RM 7 档；MA 无数值 Age |
| MA_Keywords / MA_Geo_Raw / MA_Institution | 自由文本/原值 | **MA_ 前缀** | MA 特有 |

> 注：MA 无数值 `Age` 列（只有 binned `MA_AgeGroup`）→ 公共数值 Age 由 RM 等贡献。
> 注：MicrobeAtlas 是 **SSU rRNA**(16S 原核 + 18S 真核，**非 ITS**)库（[bioRxiv 2025](https://www.biorxiv.org/content/10.1101/2025.07.18.665519v1)）；我们 GG2 表 03 步 Domain 过滤只留 Bacteria+Archaea ⇒ **实际内容 = 16S 等价**，18S 真核内容已被删。

### RM（36 列）— ✅ 已过（2026-06-05，实测取值核验）

**公共（实测同义）**：`Run` · `Project_ID` · `BioSample` · `Sex`{Male/Female} · `Sequencing_Type` · `Smoking` · `Latitude` · `Longitude` · `Library_Layout`{PAIRED/SINGLE} · `Platform`{ILLUMINA/LS454/BGISEQ/OXFORD_NANOPORE/ION_TORRENT/PACBIO_SMRT/DNBSEQ} · `Age`(+Age_start/Age_end, 数值 0–107) · `BMI`(+BMI_start/BMI_end) · `Country`(58 标准国名) · `Continent`(6 大洲) · `PMID`(PubMed id, 100% 覆盖)

**公共-派生**：`Age_Group` —— 采纳 RM 7 档为**中央统一分桶**，从公共 `Age` 派生（RM 自己的=中央；MA 5 档→`MA_AgeGroup` 前缀）。

**公共（6 边界，用户 2026-06-05 已拍板采纳）**：
- `Model` → 公共 **`Instrument`**（仪器型号客观同义；46 值如 Illumina MiSeq/NextSeq 550/GridION）
- `Region_16S` → **公共**（V4/V3-V4… 命名通用客观；各库推断方式写进数据字典）
- `Patient_ID` → 公共 **`Subject_ID`**，**值加 DB 命名空间**（`RM:100BAL1`）防跨库误配（纵向同人锚）
- `Recent_Antibiotic_Use` → **前缀 `RM_Recent_Antibiotic_Use`**（"recent"各研究定义不同，非真同义）
- `Time_Point` → **前缀 `RM_Time_Point`**（实测极杂：'pre vaccination'/'-1'/'718' 文本数字混，RM 内部都没标准化）
- `Location` → **前缀 `RM_Location`**（自由文本；地理已由 Country+Lat/Lon 覆盖）

**RM_ 前缀（各自独立标准化，必前缀）**：`RM_Phenotype`(文本) · `RM_Phenotype_ID`(MONDO/EFO，**≠ MBM 的 MeSH**) · `RM_Disease_Stage` · `RM_Complication` · `RM_Intervention` · `RM_Antibiotics_Used` · `RM_Sample_Site`(body site) · `RM_Sample_Type` · `RM_Case_Or_Control`(派生自 RM Phenotype) · `RM_Is_Healthy`(派生自 RM Phenotype)

> 36 列全覆盖：公共 19 + 公共派生 1(Age_Group) + 边界转公共 3(Instrument/Region_16S/Subject_ID) + 边界转前缀 3 + RM_前缀 10 = 36 ✓
### MicroBioMap（3 表）— ✅ 已过（2026-06-05，实测取值核验）

**结构注记**：MBM **无标准化表型表**（区别于 RM/MA 的一张干净表）。三张元数据：
- `sample_metadata.tsv`（168,464 行，**每样本，srr=Run 1:1 对组成表 obs**）— 半标准化，主表
- `projects.csv`（482 行，**项目级**，按 project 广播到样本）— MBM 作者人工策展
- `tags.tsv.gz`（349 万行，**EAV 长表**，原始 SRA/ENA 属性，数百键，**未标准化、噪声大**）— 归档原始档

**公共（原始同义 / 中央派生）**：
- `srr` → **`Run`**(a)；`project` → **`Project_ID`**(a)
- `library_strategy`（全 AMPLICON）→ **`Sequencing_Type`=Amplicon**（中央派生）
- `instrument` → **`Instrument`**（同 RM Model 槽，需值规整：454 GS / unspecified→Unknown）
- `iso`（ISO 3166-1 alpha-2 国码）→ 中央**确定性查表**（pycountry）派生 **`Country`**（US→United States）+ **`Continent`**；28k `unknown`→NA。比 `geo_loc_name` 自由文本可靠 → iso 是 MBM 的 Country 公共原始量
- tags `External Id`（小写 `samea…`=SAMEA ENA BioSample，38,852 条、与 srr **1:1**）→ 大写化 → **`BioSample`**；仅 ENA 子集有值，**SRR(NCBI)/DRR(DDBJ) 子集 NA**（真 BioSample 需回 SRA efetch 解析，标 TODO）
- tags `host_subject_id`（44.6k）→ **`Subject_ID`**，值加 `MBM:` 命名空间
- tags `host_sex`（~16k，较干净）+ `sex`（清洗，丢 not collected/数字噪声）→ **`Sex`**
- tags `host_body_mass_index`+`body_mass_index`（~18k 数值）→ **`BMI`**
- tags `host_age`+`age`（**纯数值部分**）→ **`Age`**；范围/文本（'17-29 yo'/'6 months'）丢
- tags `collection_date`（145k）→ **`Collection_Date`**（2026-06-05 升公共后回补）；tags `host_height`（~5k）→ **`Height`**（同上;MBM 无 host_weight）
- projects `amplicon`（v1-v2 / v3-v4 / v4…）→ **`Region_16S`**（中央，客观命名；MBM 项目级策展）
- `link`（论文 DOI）→ Layer2 用 **DOI→PMID 解析**（NCBI ID Converter / CrossRef）命中的回填 **`PMID`**，命不中 NA；DOI 原值另存 `MBM_DOI`

**MBM_ 前缀（必前缀 / 独立标准化 / 自由文本）**：
- `srs` → `MBM_SRS`（SRA Sample 原生 id，≠BioSample 命名空间；obs 已用 Run）
- `geo_loc_name` → `MBM_Geo_Raw`（自由文本，同 MA_Geo_Raw）
- `region`（UN M49 大区 'Europe and Northern America'）→ `MBM_Geo_Region`（≠我们 Continent 6 大洲口径；★勿与 `Region_16S` 混名）
- projects `sample_type`（stool/saliva）→ **`MBM_Sample_Type`**（★body site 红线，必前缀）
- projects `condition`（healthy/CRC/IBD…自由文本含逗号）→ **`MBM_Condition`**（自由文本表型，必前缀）
- projects `bead_beating`/`kit`/`avg_length`/`subjects`/`notes` → `MBM_Bead_Beating`/`MBM_Kit`/`MBM_Avg_Length`/`MBM_Subjects`/`MBM_Notes`
- tags `host`（物种，主 homo sapiens + 少量鼠/猴/控制污染）→ `MBM_Host`（暂无公共 Host_Species 槽）
- `MBM_DOI`（link 原值）
- tags 其余**数百键**（疾病问卷 / AGP 饮食频率 / env_* / ethnicity…）→ **全留 EAV 原始档，不进任何列**（原始噪声大、各研究口径不一 = 伪等价红线）

**丢弃（无价值，连前缀都不留）**：`pubdate`、`total_bases`、`library_source`（三者皆 SRA 投稿人随手填：16S 全 AMPLICON，METAGENOMIC/GENOMIC 标注不一致无生物学意义；首发日期/碱基数对建模无用）

> **MBM 未给公共集带来新公共列**——所有公共列均落在 MA+RM 已确立集合内，其余全 `MBM_` 前缀或丢弃 → 公共候选集不变。
> **两个延后 enrichment（Layer2 实现）**：① SRR/DRR 子集真 BioSample 经 SRA efetch 解析 ② DOI→PMID 解析回填 PMID。

### Qiita（ENA 桥 + 873 sample 模板 + 2044 prep 模板）— ✅ 已过（2026-06-05，实测取值核验）

**数据来源（关键,影响主键与 INSDC 覆盖）**：Qiita 是 UCSD 的微生物组**管理/分析平台**——研究者上传原始 reads+metadata，Qiita 用 Deblur 统一重处理，发**自己的 id**（`sample_name=study.sample`，如 `10317.000107576`）。公开时 Qiita 缴存中心 **UCSDMI** 把数据缴到 **EBI/ENA**，ENA 才铸 INSDC 号（study→PRJEB / sample→SAMEA+ERS / run→ERR）。
- ⇒ **obs 主键 = `study.sample`（Qiita 原生,永远有）**；**run + BioSample 是 INSDC 叠加层、成对出现**（ENA 桥实测 372,512 run 每行 run 与 SAMEA 都有,0 失配）→ 映射到 **309,924 个样本同时有 run+BioSample**，其余 ~52k（362,011−309,924）**两者都没有**。
- ★ **run 与 BioSample 覆盖率完全相等（~85.6%）**，不存在"有 biosample 没 run"。无号样本 = 未(以 UCSDMI 名义公开)缴存到 ENA 的（待缴/别 center/控制/老 study）→ Layer2 标 `dedup_uncertain`。

**三来源元数据**：① ENA 桥 `rawdata/Qiita/qiita_ena_bridge.tsv`（最干净 INSDC 号）② sample 模板（873，Qiita 标准化核心列 + 各 study 自有 free 列）③ prep 模板（2044，target_gene/subfragment/platform/instrument_model）。

**公共 — 来自 ENA 桥**：
- `run_accession`→ **`Run`**；`sample_accession`(SAMEA)→ **`BioSample`**（二者同覆盖 309,924,成对）
- `study_accession`(PRJEB,Qiita 缴 ENA 自得)→ **`Project_ID`**
- `instrument_platform`(ILLUMINA/LS454/ION/Nanopore)→ **`Platform`**
- `library_strategy` 杂（AMPLICON/OTHER/WGS/RNA-Seq,因桥含非16S）→**不直接用**；`Sequencing_Type` 靠 prep `target_gene=16S`→ **Amplicon** 常量
- `first_public`→ 丢（ENA 行政日期）；`secondary_sample_accession`(ERS)→ 丢（冗余）

**公共 — 来自 sample 模板**：
- `latitude`/`longitude`(98.7%)→ **`Latitude`/`Longitude`**
- `host_subject_id`(97.5% 金矿)→ **`Subject_ID`**（`QIITA:` 命名空间）
- `sex`→ **`Sex`**；`host_body_mass_index`→ **`BMI`**
- `country`→ **`Country`**（★中央重标准化:剥 `GAZ:` 本体前缀、`USA`→United States、`control sample`/`unknown`→NA）
- `host_age`+`host_age_units`→ **`Age`**（★中央**单位换算→years**:weeks/days/months 占 ~25% 必换算）
- `collection_timestamp`/`collection_date`→ **`Collection_Date`**（公共,见下方升级）
- `host_height`(+units)→ **`Height`**；`host_weight`(+units)→ **`Weight`**（公共,见下方升级）

**公共 — 来自 prep 模板**：
- `target_subfragment`(V4/V3/V9/V3-V4…)→ **`Region_16S`**（★中央归一拼写:V3V4/V3-4/V34→V3-V4）
- `instrument_model`(Illumina MiSeq…)→ **`Instrument`**；`platform`→ **`Platform`**（与桥一致）

**`QIITA_` 前缀（独立标准化/环境本体/free text）**：
- ★body site:`sample_type`/`host_body_site`/`host_body_product`/`host_body_habitat`/`qiita_sample_type` → `QIITA_*`
- EMPO/ENVO 环境本体:`empo_0/1/2/3/4`/`qiita_empo_*`/`env_biome`/`env_feature`/`env_material`/`env_package`/`envo_biome_*` → `QIITA_*`
- host 物种:`host_taxid`/`taxon_id`/`scientific_name`/`host_scientific_name`/`host_common_name` → `QIITA_*`（无公共 Host_Species 槽,同 MBM_Host）
- `geo_loc_name`→ `QIITA_Geo_Raw`；`elevation`/`altitude`/`depth`→ `QIITA_*`
- `qiita_study_id`→ `QIITA_Study_ID`（native;Project_ID 用桥 PRJEB）；`sample_name`→ obs_names 源
- `diet`/`life_stage`/`host_life_stage`/`anonymized_name`/`description`/`title` 及**各 study 自有数百 free 列** → `QIITA_*`

**丢弃**：`first_public`/`secondary_sample_accession`、`physical_specimen_*`/`dna_extracted`/`public` 等 Qiita 管理字段。

> **Qiita 触发两项"有 peer 即升公共"（用户 2026-06-05 拍板采纳）**——见公共候选集已更新：
> ① **`Collection_Date`**（采样日期,3 peer：PRIME/Qiita/MBM）② **`Height`/`Weight`**（客观 cm/kg,2 peer：PRIME/Qiita）。仅升**真实采样日期**;ENA/SRA 行政日期仍丢。

### PRIME（samples.csv 84 列 + projects.csv）— ✅ 已过（2026-06-05，实测取值核验）

**结构注记**：PRIME = curated 16S meta-analysis 资源。`samples.csv`（53,449 样本 × 84 列，每样本 `Run` 1:1）是主表，列分类丰富但**藏两个同名异义陷阱**（正是「列名相同 ≠ 含义相同」原则的活案例）。`projects.csv`（128 项目）= 项目级版本，列与 samples 子集重叠，归 PRIME_。

**★ 两个伪等价陷阱（必前缀，绝不盲合）**：
1. `Sequencing_Type`（PRIME）= **读长×单双端配置**（`150`/`2×250`/`2×300`/`2×150`/`300`），**非**公共集的"扩增子 vs 鸟枪" → **`PRIME_Read_Config`**；PRIME 实际全 16S 扩增子 → 贡献**公共 `Sequencing_Type=Amplicon`**（常量）。
2. `Host_Age`（PRIME）大量是**分桶/文本**（`30-39`/`above 70`/`child`，257 值）→ **仅纯数值**进公共 `Age`；分桶值（PRIME 自己分桶 ≠ 中央 7 桶）→ `PRIME_Host_Age_Raw`，不喂公共数值 Age。

**公共（原始同义 / 中央派生）**：
- `Run`(53,449 唯一,1:1)→ **`Run`**(a)；`BioProject`(PRJ)→ **`Project_ID`**(a)
- `BioSample`→ **`BioSample`**(★**真 SAMD/SAMN/SAMEA**,51,038 唯一,全覆盖；比 MBM 强,MBM 只 ENA 子集有真 BioSample)
- `Instrument`(Illumina MiSeq…10 值干净)→ **`Instrument`**
- `Variable_Region`(V4/V3-V4/Unknown…)→ **`Region_16S`**
- `Host_Sex`(Female/Male)→ **`Sex`**；`Host_BMI`(数值,31k 非空)→ **`BMI`**
- `Country`（★特殊命名 `America (U.S.)`/`Britain (U.K.)`,85 值）→ **`Country`**（中央**重标准化**映射→United States/United Kingdom）
- `Continent`（North America/Europe/Asia/Africa/Oceania/South America）→ **`Continent`**（映射;或从 Country 中央派生）
- `Library_Layout`（Single/Paired/**Paired (merged)/Paired (forward only)**）→ **`Library_Layout`**（归一 Paired*→PAIRED;细分→`PRIME_Library_Layout_Raw`）
- `Participant_Id`(39,002 唯一)→ **`Subject_ID`**(值加 `PRIME:` 命名空间)
- `Host_Age`(纯数值部分)→ **`Age`**（见陷阱 2）
- `Doi` → DOI→PMID 解析回填 **`PMID`**;DOI 原值→`PRIME_DOI`
- `Smoke_Status`（10 值,杂）→ **`Smoking`**（映射,用户 2026-06-05 拍板）：`no→Non-smoker`、`Former→Ex-smoker`、`yes/Daily/Current/Regularly/**Rarely/Occasionally→Smoker**`（有暴露即 Smoker）、空→`Unknown`

**`PRIME_` 前缀（独立标准化/自由文本/study-specific/无 peer）**：
- ★`Body_Site`(Stool/feces/Oral cavity…94 值)→ **`PRIME_Body_Site`**（body site 红线）
- ★`Phenotype`(Healthy/Seasonal allergies/HCC…1149 值,含多值)→ **`PRIME_Phenotype`**（自由文本）
- `Systems`(Digestive System/Gut…)→ `PRIME_Systems`；`Sequencing_Type`→ **`PRIME_Read_Config`**（陷阱 1）
- **PRIME 流水线 provenance**：`Primer_Cut`/`Primer_Par`/`Denoise_Par`/`Collapse`/`Sequencing_Quality`/`Filter_Pass`（注 `Collapse=Silva and greengenes2` ⇒ PRIME genus 表用 Silva+GG2 双折叠,**对齐我们 8,114 vocab 待 Layer1 核**）→ `PRIME_*`
- **study 设计**：`Time_series`/`Comparison`/`Matched`/`Study_Group`/`Time_Point`/`Project_name`/`Note` → `PRIME_*`
- **来源 id**（命名空间各异,无 peer）：`SRA_Study`(DRP)/`Experiment`(DRX)/`Library_Name`/`Sample_Name` → `PRIME_*`
- **宿主属性/生活方式**：`Race_or_Ethnicity`/`Diet_Type`/`Drinking_Status`/`Probiotic_Frequency`/`Teethbrushing_Frequency`/`Sleep_Duration`/`Host_Menopausal_Status`/`Host_Birth_Year`/`Host_Country_Birth`/`Dominant_Hand` → `PRIME_*`
- **临床/疾病专项（~30 列）**：`Tumor_Stage`/`Node_Stage`/`Metastasis`/`Antibiotics_Use`/`Antibiotic_Route`/`Antibiotic_Status`/`Anti_Fungal`/`Vaginal_pH`/`ALL_Type`/`Endo_*`(C_Value/M_Value/Location/Stage)/`On_PPI`/`Neoadjuvant`/`Adjuvant`/`Pepsin_*`/`RSI_*`/`Post_Pre_Surgery`/`Confirmed_*`(Adenomyosis/Cysts/Fibroids)/`*_Menstrual_Bleeding`/`Sperm_*`/`Fertility`/`Monobutyl/Monoethyl_Phthalate_Presence`/`Ihmc_medication_code`/`Abnormal_Uterine_Bleeding` → 全 `PRIME_*`
- ~~`Host_Height`/`Host_Weight`→ `PRIME_*`~~ → **升公共 `Height`/`Weight`**（2026-06-05 Qiita 也有 → 触发升公共,用户拍板）
- ~~`Collection_Date`→ `PRIME_Collection_Date`~~ → **升公共 `Collection_Date`**（2026-06-05 3 peer → 升公共,用户拍板）

**丢弃（无价值）**：`AvgSpotLen`、`Bases`（同 MBM total_bases）、`Release_Date`、`Create_Date`（SRA 行政日期）

> **PRIME 未给公共集带来新公共列**——公共列全落已确立集合内（亮点：全有真 BioSample、Doi→PMID）；其余 `PRIME_` 前缀或丢弃。
> **一个延后 enrichment（Layer2）**：Doi→PMID 解析（与 MBM link 同机制）。
> **Layer1 兼容性已核实（2026-06-05）**：PRIME genus 表 = **GG2 24.09，与我们完全一致，零版本漂移**。6,796 列中 **5,308 深层真属全部 ⊆ 8,114 词表（全路径精确匹配）**，落词表外的 1,488 全是浅层桶（深层真属外溢=0）。`Collapse=Silva and greengenes2` 仅指 PRIME 同时提供 Silva+GG2 两套折叠表（我们用 `gg_genus_*` 的 GG2 套）。⇒ **PRIME Layer1 = trivial reshape**：CSV→h5ad、5,308 属对齐 8,114 列、1,488 浅层桶求和进 `unresolved_reads`（符合组成表契约，Tier1 零重处理）。

### GMrepo / mBodyMap（仅富集，不加行）— ✅ 已过（2026-06-05，用户定纯前缀）

**性质**：同组同 schema 姊妹库（各 30 列，28 列相同；run_id **零重叠**）→ 合成**一张富集表**（concat + `meta_source∈{GMrepo,mBodyMap}`），按 **Run join** 到已有 corpus 行（**不加样本行**，无序列）。
- GMrepo 可贴 39,007(32.8%) · mBodyMap 可贴 60,479(95.8%)（详 `corpus_sample_overlap.md`）。

**🔑 用户拍板（2026-06-05）= 纯前缀(A)**：GMrepo/mBodyMap **所有列一律 `GMREPO_*`/`MBODY_*` 旁挂**，**绝不回填/覆盖母库公共列**。provenance 最干净；是否下游用富集值补全空公共列 = 下游另议（不在合并层做）。
- join 键：`run_id` → Run（不存为列）。
- 实测 id：`second_sample_id` = **真 BioSample**(SAMN/SAMEA/SAMD)；`sample_id` = SRA Sample(SRS/ERS/DRS)。`experiment_type` = AMPLICON 87k/Metagenomics 32k。
- 价值字段（均带前缀）：`BodySite`/`BodySite_info`(mBodyMap，★MA 部位弱正好补) · `disease`(MeSH D-码)/`phenotype`/`disease_stage` · `sex`/`country`/`collection_date`/`lat`/`lon` · `instrument_model`/`experiment_type` · `second_sample_id`(BioSample) · `diet`/抗生素系列。
- 丢：`nr_reads_sequenced`(同 total_bases) · 内部 id `our_project_id`/`checking`/`more`。
- 规格沿用 `metadata_sources_gmrepo_mbodymap.md` §5（列名点↔下划线统一、MeSH 解码）。
### 其它"仅采纳 metadata"的数据集 — ✅ 无（2026-06-09 用户确认）
> 仅-metadata 库**就 GMrepo + mBodyMap 两个，没有其它**。富集集合已完整。
> 🔶 **重合(overlap)问题**：仅-metadata 库 + 跨库样本重合（同 run/biosample 出现在多库）怎么贴、怎么去重，**留待后续专门讨论**（与 Layer2 去重×QC 决策相关）。

## 公共数值/地理 sanity 清洗（2026-06-09，写字典时发现并修）

剖析 V3 obs 发现公共数值/地理列有越界值。**回溯判定 bug vs 源脏**：
- **源脏（忠实解析，非 bug）**：Lat/Lon（MA `MA_Geo_Raw` 原值即 `482085 163721`/`3209326 10337221` 垃圾投影坐标；Qiita 自由文本越界）· Age（MBM/Qiita `host_age` 源含负数 -88/-957 与上千；RM/PRIME 干净）· Height/Weight（PRIME/Qiita 单位错 17462cm；**`Weight=999909` 是源缺失哨兵**）。
- **我的 2 处 gap（已修）**：① `Region_16S` 大小写未统一（小写 129,953 全来自 MBM 源 `projects.amplicon`）→ 中央 upper；② `num_or_na` 放过 `float("inf")`（BMI 24 行，全 Qiita）→ sanity 用 `np.isfinite` 拒绝。

**处理 = `common_schema.apply_common_sanity()`，在 `02_merge_metadata` 调用**（越界→NA，非掩盖 bug）：
- Lat [-90,90] · Lon [-180,180] · **(0,0) 坐标→NA**（MA 46万行占位）· Age/start/end [0,122] · BMI/start/end [8,100]（去 0/inf）· Height [30,250]cm · Weight [1,500]kg（去 999909 哨兵）
- Region_16S 统一大写；Age_Group 从**清洗后** Age 重派生。

## 待办（逐库逐列继续）
- ✅ MA（26）· ✅ RM（36）· ✅ MicroBioMap（3 表）· ✅ PRIME（84+项目）· ✅ Qiita（3 来源）· ✅ GMrepo/mBodyMap（纯前缀富集）已过 → 公共列候选集已确立（见上）
  - Qiita 触发两项升公共：`Collection_Date`（3 peer）· `Height`/`Weight`（2 peer）
- ✅ **全部库已过**（用户 2026-06-09 确认：仅-metadata 库就 GMrepo+mBodyMap 两个，无其它）→ 列设计收尾，进入实现（standardize → 02_merge_metadata）。
- ⏳ 跨库**重合/去重**专题讨论
- ⏳ 全部过完 → polish 成 HF 出版级**数据字典**（每列 名称/dtype/含义/取值标准化/来源库/缺失约定/派生逻辑）
