# 数据文件 schema 详表

## 1. `sample_info/samples.env.info.tsv`

样本元信息表，约 **2,690,735 行**，9 列（tab 分隔）：

| 列号 | 列名 | 含义 |
|------|------|------|
| 1 | MAP_SID | 样本 ID（如 ERR/SRR 编号） |
| 2 | Environments | 受控词汇环境标签（见下方详述） |
| 3 | _ | 针对 animal;human 样本的身体部位细分（gut/skin/oral/urogenital 等） |
| 4 | Technology | 测序技术（AMPLICON / WGS / RNAseq） |
| 5 | Keywords | 自由文本关键词标签（与 Environments 互补，非子集关系） |
| 6 | _ | （未命名，通常为空） |
| 7 | Project | 原始项目编号（ERP/SRP） |
| 8 | Institution | 机构 + 地理坐标 |

样本分布：Amplicon ~237 万，WGS ~26 万，RNAseq ~1.2 万。

### Environments 列要点

- **格式**：`主类;子类`（严格两级），多个标签用 `|` 分隔；约 21% 样本无标签
- **Environments vs Keywords**：独立互补，非子集关系。Environments 是受控词汇，Keywords 是自由文本
- **分类体系**：4 主类（animal/soil/aquatic/plant），共 69 个唯一标签
  - `animal`：37 子类，以 human / mouse / pig / cattle 为主
  - `soil`：9 子类，以 forest / field / agricultural 为主
  - `aquatic`：12 子类，以 marine / sediment / sea / river 为主
  - `plant`：7 子类，以 rhizosphere / leaf / wood 为主
- 详细计数 → `.claude/data/sample_env_stats.md`

### 第 3 列（body site）要点

- **非 human 专属**：所有 animal 样本共用同一套词汇（`;` 分隔，最多 3-4 级，可带 healthy/disease/infant 等修饰）
- Human 样本（473,574 个）中 93% 有标注，共 172 种 body site，主要为 gut / oral / skin / urogenital / lung
- 详细分布 → `.claude/data/sample_env_stats.md`

---

## 2. `OTU_info/otus.info.tsv`

OTU 注释信息表，跨所有相似度层级，19 列：

| 列号 | 列名 | 含义 |
|------|------|------|
| 1 | OTU | 层级 OTU ID |
| 2 | Tax | 分类域（Archaea/Bacteria 等） |
| 3 | SpeciesRep | 物种代表序列 |
| 4 | SeqCount | 该 OTU 的序列数量 |
| 5 | GoldCount | Gold 标准序列数 |
| 6 | GenomeCount | 关联基因组数 |
| 7 | TypeStrains | 模式菌株 |
| 8 | Strains | 菌株 |
| 9 | Genomes | 基因组 |
| 10 | GoldSeqs | Gold 标准序列 |
| 11 | Aliases | 别名 |
| 12-14 | GoldHit/GoldID/GoldScore | Gold 比对结果 |
| 15 | RepSpecies | 代表物种 |
| 16 | Taxaname | 分类名（GTDB 格式） |
| 17 | OrigTax | 原始分类（NCBI 格式） |
| 18 | RepSequenceID | 代表序列 GenBank ID |
| 19 | RepSequence | 代表全长 16S 序列 |

---

## 3. `OTU_info/mapref-3.0.tar.gz`

MicrobeAtlas 比对参考数据库，解压后含：

| 文件 | 含义 |
|------|------|
| `mapref-3.0.fna` | 参考全长序列 FASTA（**1,360,792 条**） |
| `mapref-3.0.fna.otutax` | 每条参考序列的层级 OTU 映射 |
| `mapref-3.0.fna.mscluster` | 序列聚类信息 |
| `mapref-3.0.fna.ncbitax` | NCBI 分类注释 |
| `otus.info` | 同 otus.info.tsv 的副本 |

`mapref-3.0.fna.otutax` 格式示例：
```
KC471280:1..1464    90_17776;96_71281;97_92606;98_125911;99_193128
```

---

## 4. `OTU_count/otus.97.allinfo`

97% OTU 详细信息，每行一个 97% OTU，20 列（tab 分隔，无表头）。列含义通过与 `otus.info.tsv` 交叉验证推断（无官方文档）：

| 列号 | 含义 | 置信度 |
|------|------|--------|
| 1 | OTU ID（格式: `MAPv3;90_x;96_x;97_x`） | 确定 |
| 2 | 某种计数（含义不明，与 SeqCount 不完全一致） | 不确定 |
| 3 | 模式菌株（TypeStrains）基因组登录号 | 较确定 |
| 4 | 菌株列表（Strains），格式 `物种名(登录号)` | 较确定 |
| 5 | 基因组列表（Genomes），比 col4 更宽泛 | 较确定 |
| 6 | 别名/同义词（Aliases） | 较确定 |
| 7 | **代表全长 16S 序列**（已与 otus.info.tsv 交叉验证一致） | 确定 |
| 8 | 代表序列 GenBank 登录号（如 `KF075129:1..1352`） | 确定 |
| 9 | GTDB 完整分类（`d__Bacteria;p__...` 格式） | 确定 |
| 10 | 最佳 GTDB 参考序列登录号 | 推断 |
| 11 | 代表序列对应基因组登录号（如 `RS_GCF_xxx~NZ_xxx`） | 较确定 |
| 13 | 与 GTDB 参考的比对相似度（0-1） | 推断 |
| 14 | 与 GTDB 参考的比对长度（bp） | 推断 |
| 15 | NCBI/SILVA 分类（`Bacteria;Firmicutes;...` 格式） | 确定 |
| 16-17 | 最佳 NCBI/Gold 参考序列登录号 | 推断 |
| 19 | 与 NCBI/Gold 参考的比对相似度（0-1） | 推断 |
| 20 | 与 NCBI/Gold 参考的比对长度（bp） | 推断 |

**提取代表序列（已验证）：**
```bash
awk -F'\t' '{printf ">%s\n%s\n", $1, $7}' OTU_count/otus.97.allinfo > otus97_rep.fasta
```

---

## 5. BIOM 丰度矩阵文件

两个 HDF5 格式的 BIOM 文件：

| 文件 | 样本数 | OTU 数 | 说明 |
|------|--------|--------|------|
| `samples-otus.97.mapped.biom.gz` | **2,690,735** | 103,271 | 全部样本 |
| `metag.minfilter.refilt.biom.gz` | **1,884,129** | 102,824 | 质量过滤后 |

读取方式（需要先解压，直接读 gzip 会 OOM）：
```python
import h5py
with h5py.File('path.biom', 'r') as f:
    sample_ids = f['sample']['ids'][:]   # 样本 ID 列表
    obs_ids = f['observation']['ids'][:]  # OTU ID 列表
    data = f['observation']['matrix']     # 稀疏矩阵（慎用，内存巨大）
```
