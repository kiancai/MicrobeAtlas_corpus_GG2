# ResMicroDb metadata 标准化方案（06 步设计）

> 阶段：上游 — 把 `metadata_all.txt` 全表（135,746 行）清成 standardized tsv（36 列）。
> 与 anndata `obs` 的对接（07 步 attach）在另一份文档里展开，本文不涉及。

---

## 1. 输入文件四件套

| 文件 | 行数 | 列数 | 编码 | 性质 |
|---|---|---|---|---|
| `rawdata/ResMicroDb/metadata_all.txt` | 135,746 | 34 | **UTF-16LE + BOM + CRLF** | **超集 / 原始**，本步基础 |
| `rawdata/ResMicroDb/sampleTable_changed_250924.tsv` | 106,464 | 34 | ASCII + CRLF | jxt 标准化后子集（不直接用，仅作对照参考） |
| `rawdata/ResMicroDb/sampleTable_16S_region_v5_250924.tsv` | 106,464 | 34 | ASCII + CRLF | 同上 + region 修正版（除 Patient_ID 外完全等同 changed） |
| `rawdata/ResMicroDb/projectTable_changed_250924.tsv` | 514 | 26 | UTF-8 | **Project 级 master 表**，本步用其 `16S_Region` 列 |

### 集合关系
- `sampleTable.Run ⊂ metadata_all.Run`（106,464 ⊂ 135,746），sampleTable 0 行新增
- jxt 丢掉的 29,282 行：测序类型 ITS/Virome/Full-16S/WGTS/18S/RNA-Seq/Nanopore/MeDIP-Seq/Tn-Seq/miWTS/AMPLICON 等 + Platform=OXFORD_NANOPORE/PACBIO_SMRT + Body_Site 空 + 上游 ASV 没产出
- 我们 anndata qc 的 93,425 样本在 metadata_all 里**100% 覆盖**且全是 16S
- projectTable 涵盖 514 PID（jxt 流程过的所有 PID）；metadata_all 的 75 个 PID 在 pt 里没有（被 jxt 完全丢弃的项目）

### 编码读盘
```python
import pandas as pd
ma = pd.read_csv('rawdata/ResMicroDb/metadata_all.txt',
                 sep='\t', encoding='utf-16-le',
                 dtype=str, keep_default_na=False, low_memory=False)
ma.columns = [c.lstrip('﻿') for c in ma.columns]   # 第一列名带 BOM
```

---

## 2. 标准化输出 schema（36 列）

`Run` 是主键。dtype 与缺失约定与 MicrobeAtlas obs schema 对齐（`Smoking`、`Sex`、`Latitude/Longitude` 等通用列名 / 取值与 06 attach 阶段共享）。

| # | 列 | dtype | 来源 | 缺失值含义 |
|---|---|---|---|---|
| 1 | `Run` | string | metadata_all.Run | — 主键，不会缺失 |
| 2 | `Project_ID` | category | metadata_all.Project_ID | — 不缺失 |
| 3 | `BioSample` | string | metadata_all.BioSample | — 不缺失 |
| 4 | `PMID` | string | metadata_all.PMID | — 不缺失（**字符串**，避免前导 0 / 大整数转换问题） |
| 5 | `Sequencing_Type` | category | metadata_all.Sequencing_Type | 原值空 |
| 6 | `Library_Layout` | category(PAIRED/SINGLE) | metadata_all.Library_Layout | 原值空（实测全表 0 缺失） |
| 7 | `Platform` | category | metadata_all.Platform | 原值空（实测全表 0 缺失） |
| 8 | `Model` | string | metadata_all.Model | 原值空 |
| 9 | `Phenotype` | string | metadata_all.Phenotype | 原值空 |
| 10 | `Phenotype_ID` | string | metadata_all.Phenotype_ID | 原值空 |
| 11 | `Disease_Stage` | string | metadata_all.Disease_Stage | 原值空 |
| 12 | `Complication` | string | metadata_all.Complication | 原值空 |
| 13 | `Intervention` | string | metadata_all.Intervention | 原值空 |
| 14 | `Smoking` | category(Smoker/Non-smoker/Ex-smoker) | metadata_all.**Smoke** 改名 | 原值空 |
| 15 | `Recent_Antibiotic_Use` | category(Yes/No) | metadata_all.**Recent_Antibiotics_Use** 改名 | 原值空 |
| 16 | `Antibiotics_Used` | string | metadata_all.Antibiotics_Used | 原值空 |
| 17 | `Sample_Site` | category | metadata_all.**Body_Site** 改名 + `Lung → Lung Tissue` | 原值空 |
| 18 | `Sample_Type` | string | metadata_all.**Body_Site_Raw** 改名 | 原值空 |
| 19 | `Sex` | category(Male/Female) | metadata_all.Sex | 原值空 |
| 20 | `Age` | float64 | metadata_all.Age（数值，区间样本天然空 → NaN） | 原值空 / 区间 |
| 21 | `Age_start` | float64 | metadata_all.**age_start** 改名（大小写对齐） | 原值空 |
| 22 | `Age_end` | float64 | metadata_all.**age_end** 改名 | 原值空 |
| 23 | `Age_Group` | category | **本步派生**：jxt 7 档 case_when（见 §4） | start/end 缺失或区间跨桶 |
| 24 | `BMI` | float64 | metadata_all.BMI | 原值空 / 区间 |
| 25 | `BMI_start` | float64 | metadata_all.BMI_start | 原值空 |
| 26 | `BMI_end` | float64 | metadata_all.BMI_end | 原值空 |
| 27 | `Country` | category | metadata_all.Country | 原值空 |
| 28 | `Continent` | category(6 大洲) | metadata_all.Continent | 原值空 |
| 29 | `Location` | string | metadata_all.Location | 原值空 |
| 30 | `Latitude` | float64 | metadata_all.Latitude（**保留 7 位精度**） | 原值空 / 解析失败 |
| 31 | `Longitude` | float64 | metadata_all.Longitude | 同上 |
| 32 | `Region_16S` | category | **本步派生**：projectTable.16S_Region 按 PID + Sequencing_Type=='16S' join（见 §5） | 非 16S 样本 / pt `-` / 项目不在 pt |
| 33 | `Patient_ID` | string | metadata_all.Patient_ID（**伪 PID 清洗**，见 §3） | 原值空 / 已清洗 |
| 34 | `Time_Point` | string | metadata_all.Time_Point（**HTML escape 修**） | 原值空 |
| 35 | `Case_Or_Control` | category(case/control) | **本步派生**：基于 Phenotype（见 §4） | Phenotype 缺失 |
| 36 | `Is_Healthy` | bool? (True/False/NA) | **本步派生**：基于 Phenotype（见 §4） | Phenotype 是 Control 或缺失 |

### dtype 与缺失约定
- `string` = pandas nullable `string[python]`，缺失 `pd.NA`
- `category` = pandas Categorical，缺失也是 `pd.NA`（不进 categories）
- `float64` 缺失 `np.nan`
- `bool?` = pandas nullable `boolean`（容许 NA），用 `Is_Healthy`
- 所有原 tsv 中的 `''` 全部 → 缺失（不留空字符串，不用 `-` 占位）

### 丢弃的 metadata_all 列（共 4 个）
- `Age_With_Interval` — 信息已在 `Age_start`/`Age_end`（实测 invariants 100% 成立）
- `BMI_With_Interval` — 同上
- `Body_Site` 已被 `Sample_Site` 取代（改名 + 值映射），原列不再保留
- `Body_Site_Raw` 已被 `Sample_Type` 取代，原列不再保留

### 不进本表的 sampleTable 独有派生量（4 个）
- `Reads / Shannon / Observed / Chao1` — alpha 多样性，**07 步合并 anndata 阶段从我们自己的 X 重算**（覆盖率 100%，jxt sampleTable 仅 84.8%）

---

## 3. 字符串清洗规则（应用于全部 string/category 列）

按以下顺序对所有 string/category 列应用：

1. **strip 前后空白**：`s.str.strip()`
   - jxt 仅对 Disease_Stage 8 行做了，我们对所有列做，无害
2. **空字符串 → NA**：`s.replace('', pd.NA)`
3. **特定列定向清洗**：
   - `Body_Site → Sample_Site`：`Lung → Lung Tissue`（精确匹配，1,706 行）
   - `Patient_ID`：精确匹配 `50.2_13.4(mean,sd)` → NA（112 行）。这是某个 project 把"50.2 ± 13.4"统计描述错误塞进 PID 列的产物
   - `Time_Point`：精确匹配 `&gt;48` → `48+`（2 行 HTML escape 残留）

### 不做的清洗（与 jxt 不一致处）
- **Lat/Lon 不截 4 位**（jxt 截位损失精度，我们保留 metadata_all 的 7 位）
- **Model 拼写不归一**（如 `454 GshFLX Titanium` / `Nextseq 550` / `Illumina Nextseq 550` 等不规范写法保留原值，下游用 Platform 4 类做 group-by 即可）
- **Antibiotics_Used / Time_Point / Disease_Stage / Complication / Intervention 自由文本不做大小写或分隔符归一**（jxt 也没做）
- **Negative Control / Positive Control / Cough swab / Oral 等控制/边缘 Sample_Site 不删行**（jxt 通过 Body_Site 白名单删了，我们保留作为合法值）
- **不做行筛选**（jxt 按 Sequencing_Type / Platform / Body_Site / Run 5 步过滤；我们保留全 135,746 行，筛选交给 07 步合 anndata）

---

## 4. 派生列逻辑（3 个）

### 4.1 `Age_Group`（7 档，照搬 jxt `1.2_clean_phyloseq.Rmd:343-356`）

```python
def age_group(start, end):
    if pd.isna(start) or pd.isna(end):
        return pd.NA
    if start >= 0  and end <= 3:  return '0-3'
    if start >  3  and end <= 18: return '3-18'
    if start > 18  and end <= 35: return '18-35'
    if start > 35  and end <= 45: return '35-45'
    if start > 45  and end <= 60: return '45-60'
    if start > 60  and end <= 75: return '60-75'
    if start > 75:                return '75+'
    return pd.NA   # 区间跨桶（如 (0,18) / (18,100)）→ NA
```

**边界处理**：
- `0-3` 桶用 `>=0` 含 0；其它桶用 `>` 上一个桶上界
- `start=3, end=3` 进 `0-3`；`start=3.0222` 才进 `3-18`
- `75+` 仅判 start > 75 无上界

**实测分布**（全表 135,746 行）：
- 标量 + 进桶：约 44,594 行（标量年龄全部能落进某桶）
- 区间 + 进桶：约 25,231 行（窄区间，如 (65,100) 不进 75+ 而进 NA；(60,75) 进 60-75）
- NA：剩余（区间太宽跨桶 / 全空）

categories 顺序：`['0-3','3-18','18-35','35-45','45-60','60-75','75+']`

### 4.2 `Case_Or_Control`（照搬 jxt `1.2_clean_phyloseq.Rmd:365-367`）

```python
def case_or_control(phenotype):
    if pd.isna(phenotype):              return pd.NA
    if phenotype in ('Control','Health'): return 'control'
    return 'case'
```

categories：`['case','control']`

### 4.3 `Is_Healthy`（照搬 jxt `1.2_clean_phyloseq.Rmd:368-371`）

```python
def is_healthy(phenotype):
    if phenotype == 'Health':          return True
    if pd.isna(phenotype):             return pd.NA
    if phenotype == 'Control':         return pd.NA
    return False
```

dtype = pandas nullable `boolean`，三态：True / False / `pd.NA`

注意 jxt 这里逻辑差异：`Case_Or_Control` 把 Control/Health 都视作 control；`Is_Healthy` 仅 Health → True，Control → NA。两者**不冗余**。

---

## 5. `Region_16S` 填充逻辑

**数据源**：`projectTable_changed_250924.tsv`（study 级，与 sample 级 v5/changed 100% 等价）

**填充规则**：

```python
# 1) 从 pt 取所有含 '16S' 的 Sequencing_Type 行
pt = pd.read_csv('projectTable_changed_250924.tsv', sep='\t',
                 dtype=str, keep_default_na=False)
pt_16s = pt[pt.Sequencing_Type.str.contains('16S')]   # 含 '16S;Metagenomics' 复合形式
pid2region = pt_16s.set_index('Project_ID')['16S_Region'].to_dict()

# 2) 按 ma.Project_ID join，但仅 ma.Sequencing_Type=='16S' 才生效
def assign_region(row):
    if row.Sequencing_Type != '16S':
        return pd.NA   # 非 16S 样本强制 NA（避免误填给同 study 的 Metag 样本）
    region = pid2region.get(row.Project_ID, '-')
    if region == '-':
        return pd.NA
    return region
```

**关键约束**：
- pt 是 study 级（同 PID 内 region 唯一，0 个 PID 多取值）
- pt 的复合 `Sequencing_Type='16S;Metagenomics'` 也含 16S 部分，region 有效，所以用 `contains('16S')` 不是 `=='16S'`
- ma 端必须 `Sequencing_Type=='16S'` 才赋值 — 防止误把 region 填给同 study 内的 Metag/Metat 样本（jxt v5/changed 在这点上误填了 1,079 行）
- 12,988 个被 jxt 丢弃但 Sequencing_Type='16S' 的样本：若所属 PID 在 pt 中，可拿到 region；若不在（75 个完全被丢弃的 PID），→ NA
- pt `-` 标记的 16S 项目（共 42 个）→ NA

**列名**：`Region_16S`（不用 `16S Region` / `16S_Region`，避免空格 / 数字开头列名带来 anndata obs / pandas attribute 访问的麻烦）

**categories 顺序**（按频率）：`['V4','V3-V4','V3','V1-V2','V1-V3','V3-V5','V4-V5','V5-V7','V4-V6','V5-V6','V1-V3/V3-V5','V1-V9','V1-V2/V3-V4','V6-V8','V6-V9','V6']`

---

## 6. metadata_all 内部 invariants（已实测，写脚本时可用作 assert）

### Age 4 列只有 3 种非空模式（135,746 行无例外）
| 模式 | Age | Age_With_Interval | age_start | age_end | 行数 |
|---|---|---|---|---|---|
| 全空 | `''` | `''` | `''` | `''` | 46,182 |
| 标量 | `48` | `48` | `48` | `48` | 44,594 |
| 区间 | `''` | `(65,100)` | `65` | `100` | 44,970 |

`Age == age_start == age_end == Age_With_Interval` 当且仅当标量；区间时 `Age == ''`、`age_start/age_end` 拆开，`Age_With_Interval == "(lo,hi)"`。

### BMI 4 列同样的镜像设计
| 模式 | 行数 |
|---|---|
| 全空 | 128,392 |
| 标量 | 6,854 |
| 区间 | 500 |

### 其它强约束
- `Country ↔ Continent`：0 个 Country 映射多个 Continent；0 行 Country 非空但 Continent 空（同步空）
- `Location ⊆ Country`：0 行 Location 非空但 Country 空
- `Body_Site_Raw == sampleTable.Sample_Type`：106,464 行 0 处不一致
- `Phenotype ↔ Phenotype_ID`：0 个 Phenotype 映射多个 PID；仅 28 行 Phenotype 非空但 PID 空
- `Sequencing_Type` 全表 16 类（含空 672 行）；`Platform` 全表 7 类无空

---

## 7. 与 jxt `1.2_clean_phyloseq.Rmd` 处理的对照

### 我们采纳的（信息一致或保留更多）

| jxt 操作 | 我们做法 |
|---|---|
| **行筛选**（Sequencing_Type 三类 / Platform≠ONT/PacBio / Body_Site 非空 / Run ⊆ ASV） | **不做** — 保留 all 全表，下游 07 步合 anndata 时自动过滤 |
| **Body_Site 10 类白名单** | **不做** — 保留 Negative Control / Positive Control / Cough swab / Oral 等所有合法值 |
| Smoke → Smoking | ✓ |
| Recent_Antibiotics_Use → Recent_Antibiotic_Use | ✓ |
| Body_Site → Sample_Site + `Lung → Lung Tissue` | ✓ |
| Body_Site_Raw → Sample_Type | ✓ |
| Disease_Stage strip 8 行 | ✓ 全字符串列都做 strip |
| Time_Point HTML escape `&gt;48 → 48+` | ✓ |
| Patient_ID `50.2_13.4(mean,sd) → -` | ✓ → NA |
| `空 → '-'` 全列 | ✗ 我们用 `pd.NA`（更干净，不与合法值 `-` 混淆） |
| `Age_Group` 7 档 case_when | ✓ 完全照搬 |
| `Case_Or_Control` | ✓ 完全照搬 |
| `Is_Healthy` | ✓ 完全照搬 |
| 加 `16S Region`（study 级 broadcast） | ✓ 命名 `Region_16S`，仅 16S 样本 |

### 我们不采纳的（jxt 做了但损失信息）

| jxt 操作 | 不采纳理由 |
|---|---|
| Lat/Lon 截 4 位 | 损失精度，metadata_all 是 7 位 |
| 删 `Phenotype_ID` | EFO/MONDO/NCIT 本体 ID 是宝贵语义锚点，必留 |
| 删 `Intervention` | 临床干预信息有用 |
| `Age` 列污染（塞区间字符串如 `(0,18)`） | 我们 `Age` 严格 float，区间在 `Age_start/Age_end` |
| 删 `BMI_start/BMI_end` 仅留 `BMI` | 我们保留 3 列（BMI / BMI_start / BMI_end） |
| 删 `Age_With_Interval`、删 `BMI_With_Interval` | ✓ 我们也删（信息已在 start/end） |

### 移到 07 步处理（不在本步）

| 操作 | 07 步做法 |
|---|---|
| `Reads / Shannon / Observed / Chao1` | 从我们 anndata `X` 直接重算（不用 jxt 的稀释到 5000，用原始 `X.sum(axis=1)` + `scipy` 算 Shannon/Chao1）；100,391 样本 100% 覆盖 |
| 行筛选 | left join `Run ∈ obs_names` 自动只保留 ASV 表里有的样本 |

---

## 8. 输出

**路径**：`results/feature_table/metadata_all.standardized.tsv`（或同目录命名约定）
**编码**：UTF-8 + LF（不再用 UTF-16LE）
**行数**：135,746（与输入 metadata_all 行数一致，未做行筛选）
**列数**：36（见 §2）
**主键**：`Run`（已 dedup 实测唯一）
**缺失值**：`pd.NA` / `np.nan` 对应 dtype；不留 `''` 也不用 `-` 占位

写盘示例：
```python
df.to_csv('metadata_all.standardized.tsv', sep='\t', index=False,
          na_rep='', encoding='utf-8', lineterminator='\n')
```
（`na_rep=''` 让 NA 写成空字符串，与 input 风格一致；下次读盘时 `keep_default_na=False` + 后处理转 NA）

---

## 9. 写脚本时的 sanity 检查（assert）

1. 输入读盘：`metadata_all.txt` 必须按 UTF-16LE 解码，第一列名 strip BOM
2. 行数：`len(df) == 135746`
3. `Run` 唯一无重复
4. Age 标量 invariant：`(df.Age.notna()) → (Age == Age_start == Age_end)`
5. Age 区间 invariant：`(df.Age.isna() & df.Age_start.notna()) → (Age_start <= Age_end)`
6. BMI 同样的两条 invariants
7. `Country ↔ Continent` 同步空：`(Country.isna()) == (Continent.isna())`
8. `Location.notna() → Country.notna()`
9. `Region_16S.notna() → Sequencing_Type == '16S'`（强制 16S 才有 region）
10. `Region_16S` 取值 ⊆ 16 类白名单
11. `Age_Group` 取值 ⊆ 7 类白名单
12. `Case_Or_Control` 取值 ⊆ {case, control, NA}
13. `Is_Healthy.notna()` ⇒ Phenotype 非空且非 Control
14. 派生列与 jxt 完全一致（取 sampleTable 里同 Run 样本，逐行核对 `Age_Group` / `Case_Or_Control`）
