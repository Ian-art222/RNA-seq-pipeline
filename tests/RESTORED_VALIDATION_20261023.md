# 全部修复恢复及 MAPQ 0 验证（2026-10-23）

新版源码：`/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline`。

已在 GitHub 原版 `8d6a389d1e8b1980333d570c03667bf16b69f900` 基础上恢复本次对话中的全部修复和改进，并将比对 BAM 的默认 MAPQ 过滤值改为 0。**完整回归 164 个测试方法全部通过：0 失败、0 错误、0 跳过，146.271 秒。**

## 恢复完整性

从 `/home/weiyuekuan/tmp/rnaseq-sample-isolation-20261023/source-before/` 及其 `source-changes.diff` 恢复删除前最终状态；33 个文件逐项 SHA256 均与当时的 `source-after.sha256.json` 一致。随后仅调整 MAPQ 默认值、新路径文档及对应测试，并增加本报告。

GitHub 原版仍完整保留在 `/home/weiyuekuan/10.RNA-seq/4.RNA-seq-pipeline`，其工作区无修改。新版的 Git 基线提交仍为上述原版提交，实际改进在本地工作区及本轮归档 diff 中，没有推送 GitHub。

已逐项复核恢复内容：

| 功能/修复 | 当前实现与测试覆盖 |
|---|---|
| raw FASTQ 清洗 | 可选 fastp，HTML/JSON 报告，完整原始配对与计数校验，成功后发布 clean 数据 |
| polyG 开关 | true 显式 `--trim_poly_g`；false 显式 `--disable_trim_poly_g`，避免仪器自动识别意外开启 |
| 清洗后读长传递 | 对每个样本两端 clean FASTQ 全量统计平均单条 read 长度，保留小数，用于对应样本 prepDE 换算 |
| 比对与资源 | 限制 HISAT2/SAMtools 实际线程，支持低线程顺序模式及引号路径，临时输出按样本隔离 |
| 参数与依赖保护 | 阻止额外参数覆盖输入/输出、二次裁剪/跳读、附加 BAM；追踪索引、剪接位点、ptf、脚本及 ctab |
| 输入输出保护 | 建图前验证所有受保护输入与输出/日志的路径、符号链接和硬链接冲突 |
| 参考与 BAM | 两种定量模式共用规范注释检查及 BAM/BAI 兼容检查；GFF3 多 gene Parent 与父子模型冲突明确报错 |
| 样本独立判链 | 每个样本自己的 JSON/env 与 StringTie 方向；合格 FR/RF 混合样本分别定量；全局 TSV 只汇总 |
| 失败继续与恢复 | 默认 profile 启用 keep-going；坏样本失败不阻断独立好样本，修复后复用未变化的成功结果 |
| prepDE 与模型完整性 | 保留原始 gene ID，与规范参考比对，仅凭严格匹配的 Ballgown 明确零证据恢复缺失模型 |
| 矩阵与 Excel | 严格样本/基因集合、数值、TPM 检查；文本 ID、大整数精确回读、工作表大小限制及原子发布 |
| 比对报告 | 完整结构及计数关系校验，拒绝缺失父标题等损坏 summary |

## MAPQ 默认 0

以下三处已同步：

- `config/config.yaml` 的 `alignment.mapq_filter: 0`。
- `workflow/Snakefile` 缺少该配置键时回退到 0。
- `workflow/scripts/align_reads.py` 独立调用省略 `--mapq` 时默认 0。

0 表示此步骤不向 `samtools view` 传 `-q`，保留低 MAPQ 及未比对记录；后续工具仍按各自算法处理这些记录。显式设置 30 时过滤继续有效。**RSeQC 的独立判链阈值 `strandness.mapq: 30` 保留**，只用于选择判链证据，不改变定量 BAM。

新增行为验证使用受控 SAM 来源及真实 SAMtools 1.20，覆盖顺序/管道两种模式、默认/显式 0/显式 30 六个组合：

| 设定 | 保留 MAPQ |
|---|---|
| 默认或显式 0 | 0、29、30、60，以及未比对记录 |
| 显式 30 | 30、60 |

另验证实际 Snakemake 配置默认、缺键回退和显式覆盖；原有关闭 Ballgown 的真实 HISAT2→Excel 测试改用默认 MAPQ 0，端到端通过。其余显式 30 的测试继续覆盖用户自定义阈值。

## 本轮完整测试

完整套件在**新源码目录**重新执行，包含原有 162 个方法及新增 2 个 MAPQ 方法，单次运行全部通过。覆盖 raw→fastp→Excel、polyG 两种设置、过滤边界、不同读长、普通/FR/RF 模式、混合方向、样本失败与续跑、错误参考、GFF 边界、参数绕过、路径引号、线程预算、Ballgown 零模型恢复/关闭、矩阵和 Excel，以及输入/代码/参数变化后的重跑依赖。

原始结果与复现脚本：

- [完整测试日志](/home/weiyuekuan/tmp/rnaseq-final-restore-20261023-115747/all-tests.log)
- [结构化测试结果](/home/weiyuekuan/tmp/rnaseq-final-restore-20261023-115747/test-summary.json)
- [完整运行脚本](/home/weiyuekuan/tmp/rnaseq-final-restore-20261023-115747/run_full_tests.py)
- [实际工具路径](/home/weiyuekuan/tmp/rnaseq-final-restore-20261023-115747/test-tool-paths.json)
- [MAPQ 六组合保留记录](/home/weiyuekuan/tmp/rnaseq-final-restore-20261023-115747/mapq-retention-results.json)
- [MAPQ 工作流专项](/home/weiyuekuan/tmp/rnaseq-final-restore-20261023-115747/mapq-workflow-tests.log)

本轮实际人工流程目录：

- raw/fastp→Excel：`/home/weiyuekuan/tmp/rnaseq-fastq-regressions-fzzcqotj`。
- FR/RF 与失败隔离：`/home/weiyuekuan/tmp/rnaseq-strand-integration-vqtuk884`。
- clean/raw 故障恢复及默认 0：`/home/weiyuekuan/tmp/rnaseq-fastq-regressions-3v6qulgf`。

## 软件与真实参考复验

重新核实实际版本：Python 3.11.9、HISAT2 2.2.1、SAMtools 1.20、StringTie 2.2.3、pandas 2.2.2、openpyxl 3.1.5、fastp 1.3.6、RSeQC 5.0.4、gffread 0.12.9、pysam 0.24.0，控制器 Snakemake 9.19.0。

主流程使用此前按固定 YAML 离线创建的环境；fastp、RSeQC/gffread 复用已有可执行文件。测试中参考准备/BAM校验脚本由现有 Python 3.13.9+pysam 0.24.0 执行，其余业务脚本使用固定 Python 3.11.9。本轮没有重新安装依赖，也没有把三套独立 YAML 环境均重新创建或全部经 `--use-conda` 执行作为已完成事项。

HISAT2 的 `hisat2-inspect` 包装脚本需要其所属 Python：参考复验脚本显式用 Python 3.11 运行，避免系统 Python 3.13 缺少旧 `imp` 模块。正式规则通过各自 conda 环境运行。

当前默认真实参考只读复验通过：39,686 个转录本、201,714 条 exon 记录、12 条序列；注释全部序列名与坐标界限与 HISAT2 索引兼容。此检查不等同于逐碱基确认组装来源相同。

- [软件版本输出](/home/weiyuekuan/tmp/rnaseq-final-restore-20261023-115747/runtime-versions.json)
- [真实参考验证](/home/weiyuekuan/tmp/rnaseq-final-restore-20261023-115747/reference-recheck/summary.json)
- [恢复记录](/home/weiyuekuan/tmp/rnaseq-final-restore-20261023-115747/restoration.json)

## 范围与保留边界

没有清洗、比对或重算真实样本。默认仍从已有 clean FASTQ 开始；使用 raw 数据时设置 `fastp.enabled: true`，按平台选择 polyG 开关。prepDE 输出仍是覆盖度估计的 read counts，不是精确 fragment counts。

严格输入扫描、共享参考错误仍可阻止建图或影响依赖任务；机器 CPU/内存/磁盘并非硬隔离。完整矩阵要求所有纳入样本成功，失败重跑可能保留历史矩阵，应检查退出码与日志。现有人工测试通过说明本次覆盖范围内未发现失败，不证明所有物种、注释和真实数据场景均无问题。
