# RNA-seq pipeline 全面复审（2026-10-23）

审查目录：`/home/weiyuekuan/10.RNA-seq/4.RNA-seq-pipeline`。

结论：本轮从配置、样本发现、raw/clean 入口、fastp、读长、HISAT2/SAMtools、参考转换、RSeQC、StringTie、prepDE、矩阵到 Excel 重新审查，确认并修复下列9类问题。最终155个测试方法全部通过，无跳过或待补跑项；未发现这些测试覆盖范围内尚未修复的流程错误。该结论不等同于证明所有物种、注释和测序数据均无问题。

## 新发现与修复

| 问题 | 实测或源码证据 | 修复后的行为 |
|---|---|---|
| HISAT2 额外参数可绕过受管输入/输出，或在读长统计后再次裁剪、跳读 | 真实2.2.1运行 `-3 10`：clean为75 bp，BAM为65 bp，而prepDE仍取75；别名/长参数缩写/短参数组合可绕过旧检查 | 阻断二次裁剪、跳读、抽样、输入输出覆盖及手动 `--rna-strandness`；允许正常评分、read group和已追踪剪接输入 |
| StringTie extra可混入另一个BAM | 原生2.2.3将两个BAM合并处理，gene count由40变80且退出成功 | 禁止附加BAM位置参数、模式切换及受管参数覆盖 |
| StringTie `--ptf` 文件未登记依赖 | extra只作为字符串参数，不追踪引用文件内容 | 按workdir解析并验证点特征文件，列入input；支持空格和单引号路径，实际StringTie运行通过 |
| StringTie省略部分零估计冗余转录本，使不同样本模型集合不一致 | 表达样本GTF缺T1，无表达样本保留T1；同次运行t_data.ctab明确T1 coverage/FPKM均0；上游标准prepDE也会失败 | 所有样本与规范参考全模型比对；仅用严格匹配的Ballgown明确零证据恢复缺失模型，逐ID记日志，保持原count公式 |
| GFF3多gene Parent静默丢基因 | gffread0.12.9对 `Parent=g1,g2` 只保留g1，旧准备器仍成功 | 转换前拒绝无法无损表示成单一gene_id的多归属转录本 |
| gene与子模型的染色体/链/坐标冲突被掩盖 | 反链gene声明和正链transcript原先可通过，gene行在标准化后消失 | GTF/GFF3均核对有声明的gene与子模型关系；报告无显式exon的原始模型ID |
| 错误配置让日志路径与输入相同，任务准备阶段删除输入 | 在人工FASTQ中将read_lengths日志配置为R1同路径，旧版Snakemake清理后R1消失 | 建立DAG前统一检查所有输出/日志与FASTQ、参考、附加输入、脚本、配置的同路径/链接冲突；也检查输出相互冲突 |
| HISAT2 summary缺少父标题但保留计数时仍获接受 | discordant或mate段残缺格式的反例 | 要求可选段落标题与计数完整配对，损坏报告报错 |
| Excel大整数精度损失被近似比较放过 | `10000000000000001` 写出并读回为 `1e16`，旧校验仍通过 | 整数逐值精确回读；超过Excel行列限制时明确失败，失败不替换已有工作簿，CSV仍可保留 |

详细独立证据：

- [比对参数、5种实际HISAT2运行](/home/weiyuekuan/tmp/rnaseq-full-audit-20261023-104011/alignment-independent/result.json)
- [StringTie额外BAM与ptf实际对照](/home/weiyuekuan/tmp/rnaseq-full-audit-20261023-104011/stringtie-options-proof/results.json)
- [多转录本与Excel专项报告](/home/weiyuekuan/tmp/rnaseq-full-audit-20261023-104011/quant-audit.md)
- [参考转换修前/修后对照](/home/weiyuekuan/tmp/rnaseq-full-audit-20261023-104011/reference-bug-reproduction/summary.json)
- [旧版输入/日志同路径复现](/home/weiyuekuan/tmp/rnaseq-full-audit-20261023-104011/log-collision-before/run.log)

## 最终测试

原有130项，加本轮新增25项，共155项。先执行当时全部154项（108.871秒），再针对新增加的 `use_ballgown:false` 分支补跑1项（7.076秒）；两次均无skip、无失败。

- [完整154项日志](/home/weiyuekuan/tmp/rnaseq-full-audit-20261023-104011/all-tests.log)
- [关闭Ballgown端到端日志](/home/weiyuekuan/tmp/rnaseq-full-audit-20261023-104011/ballgown-disabled-test.log)
- [重跑脚本与工具路径](/home/weiyuekuan/tmp/rnaseq-full-audit-20261023-104011/run_full_tests.py)

覆盖 raw→fastp→Excel、polyG显式开关及实际重跑、完整clean读长传递、普通/FR/RF判链及失败门控、不同读长、同基因多个非零剪接异构体、冗余零模型恢复、关闭Ballgown、ID文本保持、参数边界、坏配对、错误参考、空数据、索引/脚本/报告/ctab重跑依赖和路径引号。两个非零剪接异构体独立实测各20条reads，transcript counts为20和20，gene count为40。

最终合成端到端目录：

- [raw/fastp与polyG切换](/home/weiyuekuan/tmp/rnaseq-fastq-regressions-uq31elbr)
- [clean普通与链特异分支](/home/weiyuekuan/tmp/rnaseq-fastq-regressions-x0uqxlfg)
- [FR/RF和失败门控](/home/weiyuekuan/tmp/rnaseq-strand-integration-axmv9sdf)

当前真实默认参考仅做只读兼容核实：39,686个transcript、201,714个exon、12条染色体通过。[参考验证记录](/home/weiyuekuan/tmp/rnaseq-full-audit-20261023-104011/reference-recheck/report.json)。没有清洗、重比对或重算真实样本。

## 软件环境核实

主流程按原 `workflow/envs/rnaseq.yaml` 使用本地缓存离线创建63包环境，以下版本全部实际启动/导入并用于最终回归：Python3.11.9、HISAT2 2.2.1、SAMtools1.20、StringTie2.2.3、pandas2.2.2、openpyxl3.1.5。间接依赖HTSlib解为1.24，实际包清单已保存。

fastp1.3.6、RSeQC5.0.4、pysam0.24.0、gffread0.12.9复用现有可执行环境；判链/参考测试Python为现有3.13，控制器Snakemake9.19.0。三个YAML均能在本地缓存索引上求解，但fastp独立环境因镜像缓存身份/包存档不齐未重建，strandness独立Python3.11环境因缺相应pysam等缓存未安装。因此不声称三套YAML环境都已完成原样部署。

安装仅使用本地缓存，临时清除代理并使用 `--offline`；没有网络下载，也未改已有环境或原环境YAML。离线创建的主环境在审查目录中，不替换日常运行环境。

[环境核实结果](/home/weiyuekuan/tmp/rnaseq-full-audit-20261023-104011/conda-offline/summary.json) · [主环境完整包清单](/home/weiyuekuan/tmp/rnaseq-full-audit-20261023-104011/conda-offline/rnaseq-installed-packages.json)

## 保留的科学和运行边界

- prepDE仍使用 `ceil(coverage × exon_length / 每样本clean平均单条read长度)`，输出为覆盖度估计read counts，不能当成精确PE fragment counts。
- RSeQC5.0.4按BAM顺序取前N条合格alignment，PE两端分开计数；按transcript span（包含intron）判断，不是全基因组随机exon抽样。JSON已记录这些边界；极少量可用reads按现有策略警告，没有擅自添加新的最低样本量阈值。
- `use_ballgown:false`在参考模型完整时实际跑通；若StringTie省略模型且没有Ballgown明确零证据，继续报错并提示启用Ballgown或检查注释。不会无证据补零。
- 一个运行对应一套参考和统一polyG策略。参考名称/长度/坐标一致不等同于逐碱基确认组装来源；不同平台/参考需按设计正确分组。
- 真实生物样本的测序质量、批次、重复设计和差异表达统计不属于本轮合成数据软件验证结论。

官方行为核对采用本机已安装源码/帮助，以及 [HISAT2手册](https://daehwankimlab.github.io/hisat2/manual/)、[StringTie手册](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual)、[StringTie2.2.3源码](https://github.com/gpertea/stringtie/blob/v2.2.3/stringtie.cpp)、[fastp1.3.6源码](https://github.com/OpenGene/fastp/tree/v1.3.6)。

修改前备份：[source-before](/home/weiyuekuan/tmp/rnaseq-full-audit-20261023-104011/source-before)。本轮详细diff、源码SHA256及验证汇总保存在同级审查目录。
