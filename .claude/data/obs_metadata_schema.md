# obs metadata schema（MicrobeAtlas / 06 步）

脚本：`scripts/MicrobeAtlas/06_attach_metadata.py`
输入：`results/feature_table/gg2.full.qc.h5ad` (1,762,635 × 6,306, obs 空)
输出：`results/feature_table/gg2.full.qc.with_meta.h5ad`（871 MB）

把 9 列原始 `samples.env.info.tsv` 解析、拆分并按 `obs_names`（MAP_SID）对齐写入 obs。obs index 重置为字符串 `RangeIndex`（不再继承 MAP_SID）。

---

## 1. 26 列 obs schema

`MA_` 前缀 = MicrobeAtlas 特有；无前缀 = 与 ResMicroDb 通用（为后续两库 outer join 预留）。

| # | 列 | dtype | 来源 / 解析 | NA 含义 |
|---|---|---|---|---|
| 1 | `Database` | category | 常量 `MicrobeAtlas`，categories 含 `ResMicroDb` 占位 | 不会有 NA |
| 2 | `MA_Sample_ID` | string | col1 原 MAP_SID | 不会有 NA |
| 3 | `Run` | string | MAP_SID `.` 前段（SRR/ERR/DRR） | 不会有 NA |
| 4 | `BioSample` | string | MAP_SID `.` 后段（SRS/ERS/DRS） | 不会有 NA |
| 5 | `Project_ID` | category | col7 | 不会有 NA |
| 6 | `Sequencing_Type` | category(AMPLICON/WGS/RNAseq) | col4 | 原 col4 空 |
| 7 | `MA_Env_Animal` | bool | col2 含 `animal;*` 标签 | 不会有 NA |
| 8 | `MA_Env_Animal_Sub` | string | 同主类多子类 `\|` 连；裸 `animal` → `""` | 主类 flag=False |
| 9-10 | `MA_Env_Soil` / `_Sub` | bool / string | 同上 | 同上 |
| 11-12 | `MA_Env_Aquatic` / `_Sub` | bool / string | 同上 | 同上 |
| 13-14 | `MA_Env_Plant` / `_Sub` | bool / string | 同上 | 同上 |
| 15 | `MA_IsHuman` | category(Human/HumanMix) | col2 派生 | col2 不含 `animal;human` |
| 16 | `MA_SampleSite_Raw` | string | col3 原值（242 唯一） | 原 col3 空 |
| 17 | `MA_SampleSite` | category(8 部位) | col3 token 落 SampleSite 槽 | col3 无相关 token |
| 18 | `MA_Health` | string | col3 健康类 token 用 `;` 拼回（不解冲突） | col3 无健康类 token |
| 19 | `Sex` | category(female/male) | col3 token | col3 无 |
| 20 | `MA_AgeGroup` | category(infant/baby/toddler/adult/elderly) | col3 token | col3 无 |
| 21 | `Smoking` | category(Smoker/Non-smoker) | col3 含 `smoker` → Smoker | col3 无 (ResMicroDb 阶段会用 Non-smoker 填充) |
| 22 | `MA_Keywords` | string | col5 原值 | 原 col5 空 |
| 23 | `MA_Geo_Raw` | string | col9 原值 | 原 col9 空 |
| 24 | `Latitude` | float64 | col9 split 后 float[0] | 解析失败 NaN |
| 25 | `Longitude` | float64 | col9 split 后 float[1] | 解析失败 NaN |
| 26 | `MA_Institution` | string | col8 原值 | 原 col8 空 |

### dtype 选择
- **string**（pandas nullable `string[python]`）：开放式自由文本，缺失值用 `pd.NA`，比 object 干净
- **category**：取值有限（≤ 几百），节省内存且利于下游 group-by
- **bool**：`MA_Env_*` 主类 flag，无 NA
- **float64**：经纬度，缺失用 `np.nan`（float 没有 NA）

读盘后 anndata 会把大部分 string 列**自动转成 category**（4 个 Sub 也会，因此 `_Sub` 的 categories 设计上含空字符串 `""`）。这是 anndata 默认行为，语义不变。

---

## 2. col2 Environments 解析

格式：`主类;子类`，多标签 `|` 分隔，约 21% 样本无标签。

```
animal;human|animal;mouse  → 4 主类 flag: animal=True, soil/aquatic/plant=False
                              MA_Env_Animal_Sub = "human|mouse"
animal                     → animal=True, MA_Env_Animal_Sub = ""
（空）                       → 全部 False, 4 个 Sub = NA
```

### round-trip 校验
脚本里 assert：从 8 列重建 → 与原 col2 set 相等，1,762,635 行无一不可逆。

### MA_IsHuman 严格定义
- `Human`：col2 中 `animal;*` 子类**只有** human，无其他动物（339,094）
- `HumanMix`：含 `animal;human` 且含其他 animal 子类（11,572）
- `NA`：不含 `animal;human` 标签（1,411,969）

> "纯 human" 而非"含 human"——避免下游 human-only 子集被混入鼠/猪样本。

---

## 3. col3 28 token 白名单

实测 col3 只有 242 个唯一字符串、共 28 个唯一 token（`;` 分隔，不含 `|`）。脚本 assert 全表无未识别 token。28 token 分到 5 个语义槽：

| 槽位 | tokens |
|---|---|
| `MA_SampleSite` (8) | gut, skin, oral, urogenital, lung, gastric, nose, bone |
| `Sex` (2) | female, male |
| `MA_AgeGroup` (5) | infant, baby, toddler, adult, elderly |
| `Smoking` (1 → Smoker) | smoker |
| `MA_Health` (12) | healthy, disease, infection, inflammatory bowel disease, dermatitis, cystic fibrosis, cholera, malaria, pneumonia, typhoid fever, tuberculosis, patient |

### 取值规则
- SampleSite/Sex/AgeGroup/Smoking：白名单内**第一个**命中 token（assert AgeGroup/Sex 互斥成立）
- Health：所有命中 token 用 `;` 拼回（不解冲突，例如 `healthy;disease` 会原样保留）
- 空 col3 → 全部 NA

### 为什么 Health 不解冲突 / 不细分子状态
- `infection` 在 col3 里就是 `disease` 的细化，没必要再起专列
- 这些标签只有几千~几万样本带，相对总样本占比低，再切分意义不大；保留单 string 列由下游按需 query

### 为什么 Smoking 留 Non-smoker 占位
- MicrobeAtlas col3 里只有 `smoker` token，无显式"non-smoker"
- ResMicroDb 阶段（02 步整合该库 metadata 时）会有 Non-smoker 标签，因此 `Smoking` 的 categories 提前包含 `["Smoker", "Non-smoker"]`，方便后续两库 concat 时 dtype 对齐

---

## 4. col9 经纬度解析

格式预期：`lat lon` 空格分隔的两个 float。脏值类型实测：

| 输入 | 处理 |
|---|---|
| `36.7 -98.4` | (36.7, -98.4) ✓ |
| `""` | (NaN, NaN), `MA_Geo_Raw=NA` |
| `0 0` | (0.0, 0.0) — 形式合法，原样保留 |
| `-122.726`（单值） | (NaN, NaN) |
| `002.9167`（单值带前导零） | (NaN, NaN) |
| `see Institution` 等非数 | (NaN, NaN) |

成功率：1,253,260 / 1,762,635 = **71.1%** 有 lat/lon 数值。原始字符串总保留在 `MA_Geo_Raw`，未来若要放宽规则重解析，可从 raw 出发。

### 范围 sanity（原始数据脏）
实测 `lat ∈ [-90, 518240]`、`lon ∈ [-6093417, 1140019]`，**远超合法范围**。脚本忠实保留，**不**做范围过滤——下游若要"地理上有意义的坐标"需自行加 `-90 ≤ lat ≤ 90 and -180 ≤ lon ≤ 180`。

---

## 5. 缺失值表示约定（统一性）

| dtype | 缺失值 | 备注 |
|---|---|---|
| `string` / `category` | `pd.NA` | 原 tsv 的 `""` 全部转 NA |
| `float64` | `np.nan` | float 无 NA，pandas 标准 |
| `bool` | 无缺失 | `MA_Env_*` flag 严格 True/False |

**唯一允许内容为 `""` 的列**：4 个 `MA_Env_*_Sub`，且仅在主类 flag=True 但无具体子类时（裸 `animal` 标签）。这是**有意区分**：
- `Sub = NA` ↔ 主类 flag=False（没有这个主类）
- `Sub = ""` ↔ 主类 flag=True 但无子类细分
- `Sub = "human"` 等 ↔ 主类 flag=True 且子类为 human

读 h5ad 后 anndata 把 Sub 列转 category，categories 里会含 `""`——这是设计如此，不要 strip。

---

## 6. 跨库 merge 策略（为 ResMicroDb 02 步铺路）

| 命名 | 含义 | merge 行为 |
|---|---|---|
| `MA_*` | MicrobeAtlas 特有，ResMicroDb 不会有 | concat 时 ResMicroDb 行该列填 NA |
| 无前缀 | 两库通用 | 两库都填，concat 后无 NA 缺位 |

通用列（10 个）：`Database, Run, BioSample, Project_ID, Sequencing_Type, MA_Sample_ID*, Sex, Smoking, Latitude, Longitude`

> *`MA_Sample_ID` 名字带 `MA_` 前缀但 ResMicroDb 也会有对应字段（届时改为 `Sample_ID` 通用列或保留 `MA_/RD_` 双前缀，留待 02 步决策）。

ResMicroDb 端会单独建立：`RD_Project, RD_*` 等特有列。

---

## 7. 校验

### 脚本内 assert（写盘前）
1. `obs_names` ↔ tsv MAP_SID 100% 命中
2. col2 拆 8 列再拼回 == 原值（set 相等，全 1.76M 行）
3. col3 全表 token ⊆ 28 词白名单
4. AgeGroup / Sex 同 col3 内互斥（max ≤ 1）
5. MAP_SID 拆 Run/BioSample 后两段都非空

### 独立 audit（写盘后，从 h5ad 重读）
独立脚本 `/tmp/audit_06_meta.py`，从输出 h5ad + 原 tsv 重新解析所有列，11 项检查全过：
1. 行序：`MA_Sample_ID == 原 obs_names`
2. tsv 重对齐 1,762,635 行无丢失
3. `Database` 全 MicrobeAtlas、categories 含 ResMicroDb 占位
4. `Run + "." + BioSample == MA_Sample_ID`
5. 6 个 raw 列（Sequencing_Type / Project_ID / Institution / Keywords / Geo_Raw / SampleSite_Raw）NA ↔ 空字符串一对一
6. `MA_Env_*` 4 主类 flag + Sub 与 col2 重解析逐行相等
7. `MA_IsHuman` 重算严格定义一致
8. col3 → 5 槽位重算逐行相等
9. `Latitude/Longitude` 重解析 + NaN 行号一致
10. 非 `_Sub` 字符串列 0 个空字符串泄漏
11. 类别列 categories 不含 `""`（_Sub 豁免）

### 每列 NA 计数（参考）
```
Database              0           Sequencing_Type  20,519
MA_Sample_ID          0           Project_ID            0
Run                   0           MA_Institution  1,594,422
BioSample             0           MA_Keywords      28,248
MA_Env_Animal         0           MA_Geo_Raw      508,163
MA_Env_Animal_Sub  849,819        Latitude        509,375
MA_Env_Soil           0           Longitude       509,375
MA_Env_Soil_Sub  1,499,688        MA_SampleSite_Raw 961,082
MA_Env_Aquatic        0           MA_SampleSite   978,106
MA_Env_Aquatic_Sub 1,486,837      MA_Health     1,622,670
MA_Env_Plant          0           Sex           1,719,107
MA_Env_Plant_Sub  1,594,881       MA_AgeGroup   1,687,577
MA_IsHuman       1,411,969        Smoking       1,762,216
```
（Latitude/Longitude NaN 数完全相同，符合预期）

---

## 8. 写盘陷阱：anndata 0.10.x 与 nullable string

**症状**：用 `pd.array(..., dtype="string")` 建 obs 列后调用 `adata.write_h5ad(...)`，文件被创建（662 MB）但 `n_vars=0`、obs 仅剩 Database + `_index`，**无显式报错**。

**根因**：anndata 0.10.x 默认 `settings.allow_write_nullable_strings = False`——拒绝写 `pd.arrays.StringArray`，因 0.11 之前其磁盘格式未稳定。运行时实际抛了 `RuntimeError`，但被 `conda run` 的 stderr 缓冲吞了。HDF5 已 truncate 旧文件并写到第一个 string 列才崩，造成"文件存在但内容残缺"的假象。

**修复**：脚本顶部 import 后加一行：
```python
ad.settings.allow_write_nullable_strings = True
```
0.11+ 默认开，本 setting 取消。

**调试技巧**：怀疑被 `conda run` 吞 stderr 时，改用 `python script.py > /tmp/log 2>&1` 显式落盘，traceback 才会出现。
