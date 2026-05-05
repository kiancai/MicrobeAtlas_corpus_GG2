# samples.env.info.tsv — 环境标签详细统计

数据来源：`rawdata/MicrobeAtlas/sample_info/samples.env.info.tsv`，统计日期：2026-05-01

## 多标签分布（Environments 列，`|` 分隔）

| 标签数 | 样本数 |
|--------|--------|
| 无标签 | 564,385 |
| 1 个 | 1,632,496 |
| 2 个 | 419,153 |
| 3 个 | 65,986 |
| 4 个 | 5,286 |
| 5 个 | 3,428 |

## 环境分类体系（4 主类，共 69 个唯一标签）

每个标签格式为 `主类;子类`（严格两级），多标签用 `|` 分隔。

### animal（37 个子类，展开计数 1,233,286）

| 标签 | 计数 |
|------|------|
| animal;human | 473,574 |
| animal | 437,117 |
| animal;mouse | 119,563 |
| animal;pig | 36,173 |
| animal;cattle | 34,124 |
| animal;bird | 22,767 |
| animal;fish | 19,370 |
| animal;dog | 13,955 |
| animal;rat | 12,721 |
| animal;bee | 7,953 |
| animal;mosquito | 7,362 |
| animal;macaque | 5,911 |
| animal;horse | 5,935 |
| animal;bat | 5,662 |
| animal;sheep | 4,955 |
| animal;tick | 4,339 |
| animal;termite | 3,402 |
| animal;zebrafish | 2,570 |
| animal;cat | 2,503 |
| animal;bumblebee | 2,484 |
| animal;goat | 2,214 |
| animal;fly | 1,762 |
| animal;baboon | 1,296 |
| animal;hamster | 1,122 |
| animal;swallow | 861 |
| animal;chimpanzee | 798 |
| animal;whale | 603 |
| animal;pigeon | 565 |
| animal;tadpole | 527 |
| animal;sparrow | 345 |
| animal;shark | 320 |
| animal;dolphin | 174 |
| animal;gorilla | 87 |
| animal;roe | 75 |
| animal;seagull | 64 |
| animal;eagle | 27 |
| animal;stork | 4 |
| animal;fruitfly | 2 |

### soil（9 个子类，展开计数 548,141）

| 标签 | 计数 |
|------|------|
| soil | 312,786 |
| soil;forest | 71,555 |
| soil;field | 70,999 |
| soil;agricultural | 29,873 |
| soil;farm | 22,060 |
| soil;tundra | 9,865 |
| soil;shrub | 8,566 |
| soil;desert | 7,962 |
| soil;peatland | 7,325 |
| soil;paddy | 7,150 |

### aquatic（12 个子类，展开计数 540,749）

| 标签 | 计数 |
|------|------|
| aquatic | 92,267 |
| aquatic;marine | 109,435 |
| aquatic;sediment | 68,027 |
| aquatic;sea | 67,431 |
| aquatic;river | 54,596 |
| aquatic;lake | 49,034 |
| aquatic;waste water | 42,427 |
| aquatic;ocean | 30,579 |
| aquatic;estuary | 11,221 |
| aquatic;groundwater | 5,845 |
| aquatic;reservoir | 5,556 |
| aquatic;ice | 3,603 |
| aquatic;brine | 728 |

### plant（7 个子类，展开计数 384,874）

| 标签 | 计数 |
|------|------|
| plant | 228,780 |
| plant;rhizosphere | 83,097 |
| plant;leaf | 35,093 |
| plant;wood | 23,054 |
| plant;seed | 6,223 |
| plant;stem | 4,688 |
| plant;flower | 3,702 |
| plant;sprout | 237 |

## Body site 列（第3列）统计

- **非 human 专属**：human、mouse、pig、bird 等所有 animal 样本共用同一套 body site 词汇
- 格式：`;` 分隔，最多 3-4 级，可带状态修饰词（healthy/disease/infant/male/female）

### Human（`animal;human`）样本 body site Top 分布

总样本：473,574；有标注：442,921（~93%）；唯一值：172 种

| body site | 样本数 |
|-----------|--------|
| gut | 175,911 |
| oral | 37,618 |
| skin | 24,160 |
| urogenital | 22,861 |
| gut;infant | 20,507 |
| gut;healthy | 12,228 |
| urogenital;female;disease | 10,917 |
| lung | 10,716 |
| gut;disease | 9,997 |
| gut;infection | 5,282 |
| gut;male | 5,177 |
| oral;male | 7,493 |
| oral;female | 5,874 |
| skin;baby | 5,077 |
| gut;disease;inflammatory bowel disease | 4,242 |
| gut;female | 4,188 |
| lung;infant | 3,807 |
| urogenital;female | 3,474 |
| skin;male;disease;dermatitis | 2,904 |
| gastric | 2,835 |
| nose | 2,470 |
