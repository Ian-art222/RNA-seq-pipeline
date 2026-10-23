# fastp 接入验证（2026-10-23）

已完成 raw FASTQ → 可选 fastp → clean FASTQ → 原有比对、定量与 Excel 流程的源码修改和验证。默认 `fastp.enabled: false` 保持已有 clean 输入方式；设为 true 启用清洗。`trim_poly_g` 为本次运行全部样本的显式布尔开关。

## 验证结果

共覆盖 **130 个测试方法，全部已通过，无待补跑的跳过项**：

- 原有117项：首轮114项通过，3项因 PATH 未包含现有 gffread 而跳过；补上本地 gffread 0.12.9 后，这3项全部通过。
- 新增13项：全部通过，无跳过。覆盖真实 fastp 的过滤边界、polyG 开关、损坏/错配输入、输入输出别名保护、失败清理、零 clean 报告和 raw 到 Excel。
- 在13项通过后增强 polyG 切换的端到端断言，单独补跑该方法并通过，未把这次重复执行另计为新测试。

完整日志：[原117项首轮](/home/weiyuekuan/tmp/rnaseq-fastp-20261023-094653/existing-tests.log)、[3项补跑](/home/weiyuekuan/tmp/rnaseq-fastp-20261023-094653/gffread-tests.log)、[新增13项](/home/weiyuekuan/tmp/rnaseq-fastp-20261023-094653/fastp-tests.log)、[增强后的polyG集成补跑](/home/weiyuekuan/tmp/rnaseq-fastp-20261023-094653/fastp-rerun-tests.log)。

另有独立 Snakemake dry-run 检查：clean/raw 两种入口、目录重叠、非法文件后缀、符号/硬链接及 N 阈值越界，详见 [workflow-dryrun](/home/weiyuekuan/tmp/rnaseq-fastp-20261023-094653/workflow-dryrun)。

## 实际工具测试

固定并实际使用 fastp **1.3.6**，保留 `--detect_adapter_for_pe -q 20 -u 30 -n 5 -l 40`，默认传 `-w 8`；Snakemake 分配3线程时实测传 `-w 3`。`-w` 指工作线程，fastp 另有 I/O 辅助线程。

- 带 NextSeq 仪器标识、20 bp polyG 尾的80 bp合成 reads：显式关闭时80 bp，开启时该测试序列为60 bp，原始输入 SHA256 不变。
- Q20通过/Q19失败；80 bp中24个低质量碱基通过、25个失败；N5通过/N6失败；40 bp通过/39 bp失败。合格 read 内的低质量碱基仍保留，证实没有额外启用质量窗口剪切。
- R1/R2的ID错配、记录数不等、FASTQ截断均在 fastp 前阻断；低质量过滤不会掩盖配对问题。
- 输出或日志与 raw 的符号链接/硬链接冲突会阻断。fastp失败不发布临时 clean 文件；所有 reads 被过滤时保留HTML/JSON，后续读长检查阻断比对。
- 两个样本先在关闭polyG时从raw运行至Excel。随后给部分合成reads追加20G，形成75/95 bp及150/170 bp混合raw长度，开启polyG并实际重跑，恢复75/150 bp clean序列，实测读长正确进入prepDE。两个基因在各样本的估算count均为201，与首轮矩阵逐字节一致。
- 无变动时不重跑；修改fastp脚本或参数、删除HTML均触发对应任务。含空格/单引号路径通过。raw入口与链特异分支的DAG通过；已有clean入口的普通和链特异分支均有实际端到端回归。

最终 raw→Excel 和 polyG 重跑产物：[合成运行目录](/home/weiyuekuan/tmp/rnaseq-fastq-regressions-dgrfimkr)。fastp 单项测试产物：[核心测试](/home/weiyuekuan/tmp/rnaseq-fastp-tests-uwhzlza0)。

## 源码核对与验证范围

[fastp 官方 polyG 说明](https://github.com/OpenGene/fastp#polyg-tail-trimming)确认省略开关仍可能按仪器ID自动启用，因此关闭分支显式使用 `--disable_trim_poly_g`。参数范围按 [v1.3.6 options.cpp](https://raw.githubusercontent.com/OpenGene/fastp/v1.3.6/src/options.cpp)核对。其[polyG算法](https://raw.githubusercontent.com/OpenGene/fastp/v1.3.6/src/polyx.cpp)允许少量非G碱基，因此不能把所有数据的裁剪量假定为末尾连续G的长度；pipeline始终完整测量实际clean读长。

本次复用现有工具，没有下载安装，没有运行真实样本清洗或改写历史分析结果。其他实际测试工具为 Snakemake 9.19.0、HISAT2 2.2.1、StringTie 2.2.3、samtools 1.23.1、RSeQC 5.0.4、pysam 0.24.0、gffread 0.12.9；测试解释器和下游表格解释器沿用原验证环境。新增 fastp Conda YAML 固定 Python 3.11.9、fastp 1.3.6，但未重新求解/安装整套Conda环境，不宣称全部固定依赖组合已经安装测试。

修改前备份：[source-before](/home/weiyuekuan/tmp/rnaseq-fastp-20261023-094653/source-before)。测试重现方式及所需环境变量见 [README.md](../README.md) 的测试部分。
