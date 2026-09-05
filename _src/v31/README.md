# V3.1 metadata-only producer

本目录属于独立母库生产仓库，不属于 V3 的下游 split/task 脚本。固定输入路径与 SHA256 在
`build_release.py`；执行合同在项目 `tmp/20260905_mcfcorpus_v31_build/PLAN.md`。

从 MCFProjet 根目录、MiCoFormerV2 环境执行（大任务先核对当时 Slurm 资源）：

```bash
python data/MCFCorpus_incoming/_src/v31/test_release.py
python data/MCFCorpus_incoming/_src/v31/build_release.py build --output tmp/20260905_mcfcorpus_v31_build/staging/MCFCorpusV3.1.gg2.h5ad
python data/MCFCorpus_incoming/_src/v31/build_release.py verify --output tmp/20260905_mcfcorpus_v31_build/staging/MCFCorpusV3.1.gg2.h5ad --report tmp/20260905_mcfcorpus_v31_build/audit/verification.json
```

build/report 拒绝覆盖既有文件；重放请选新的实验内输出路径。验证不修改 h5ad，不自行发布或改默认版本。
先复制固定 V3，再只更新 qc_pass/uns.qc_def，新增 17 个 obs 字段及三个 uns 字段；原 QC 定义和 flag 都保留。
未授权变化的 HDF5 对象逐项比较属性、类型、shape 和全部值（包含原 NaN），不是只比 shape。
还全量比对候选每行的身份、标签、QC；检查词表前缀、深度、逐库计数和 AnnData 可读性。
脚本以源码 SHA256 绑定，不将未提交源码伪称为旧 repo commit 的内容。

实际成功运行采用 `staging/retry/MCFCorpusV3.1.gg2.h5ad`（首轮 verifier 节点拼写修正后重建，原尝试保留）。
独立 verifier PASS 后，`python data/MCFCorpus_incoming/_src/v31/publish_release.py` 读取本实验 verification.json，
再次核实输出/原 V3/源码 hash，再独占创建 `data/MCFCorpusV3.1/`。发布脚本拒绝覆盖既有版本。
该目录附带 release_manifest、verification、原准备包 manifest 和官方词表快照。逻辑重放需重新验证；
构建时间写入 provenance，因此不声称另一次构建具有相同文件级 SHA256。

此脚本只重建 V3→V3.1 metadata 更新，不声称修复旧 raw→V3 管线的环境/分类器缺口。
保留原 V3；没有 biome 不自动 QC false；不生产 train/val 或任务 sidecar，不自动运行模型。
