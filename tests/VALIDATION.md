# 修复与验证记录（2026-10-22）

目标源码：`/home/weiyuekuan/10.RNA-seq/4.RNA-seq-pipeline`。

## 最终结果

统一执行 `python -m unittest discover -s tests -v`，开启真实 FASTQ 回归和所有已有工具覆盖：**117 项全部通过，无跳过项，耗时70.121秒**。

| 测试文件 | 项数 | 验证内容 |
|---|---:|---|
| test_alignment_report.py | 18 | 必需字段、计数/比例一致性、可选 HISAT2 模式、样本 ID |
| test_quantification.py | 25 | 单样本读长、canonical ID、GFF、重复记录、零表达、数值和 Excel |
| test_read_lengths_and_resources.py | 11 | 全量双端均值、配对/截断错误、实际线程预算、剪接输入路径、失败传播 |
| test_reference_validation.py | 6 | BAM/BAI、排序、参考序列/范围、旧成功状态清除 |
| test_strand_bed.py | 25 | BED/GTF/GFF3 标准化、转义 ID、未建模基因、子特征、输出保护 |
| test_strandness.py | 23 | FR/RF、质量门控、字段解析、真实 RSeQC、失败诊断 |
| test_workflow_integration.py | 6 | 预建人工 BAM 的 FR/RF、普通模式、方向冲突/模糊、非法配置 |
| test_workflow_regressions.py | 3 | 真实 FASTQ 两分支、默认重跑、索引/脚本/Ballgown、错误参考阻断 |

[完整测试日志](/home/weiyuekuan/tmp/rnaseq-fix-20261022-221623/all-tests.log)。

## FASTQ 到 Excel 的真实工具回归

使用10 kb人工参考、两个基因（ID为 `001`、`Gene|literal`，都有 gene_name），两个样本（`001` 为75 bp，`NA` 为150 bp），各200对人工 reads。

- 普通模式：`--cores 1`，工作目录包含空格及单引号；HISAT2、view、sort顺序执行，单阶段不超预算。
- 链特异模式：`--cores 4`，真实RSeQC通过FR门控，流水线线程之和不超预算。
- 两种模式最终矩阵均为下表，且与 canonical raw count 逐单元一致；150 bp样本按实际读长独立重算一致，与默认75 bp计算不同。

| GeneID | 样本001（75 bp） | 样本NA（150 bp） |
|---|---:|---:|
| 001 | 201 | 201 |
| Gene\|literal | 201 | 201 |

每个基因每样本100对reads，即200条reads；覆盖度输出精度与ceil使估算值为201。这验证的是相同覆盖度公式、正确单条read长度以及ID映射，不将估算count冒充精确fragment count。

每样本TPM总和约100万；Excel的ID、样本名、列名均保持文本，所有数值与CSV一致。完整测试目录保留在：

- [普通分支](</home/weiyuekuan/tmp/rnaseq-fastq-regressions-iwydbmwr/unstranded work 'quoted'>)
- [链特异分支](</home/weiyuekuan/tmp/rnaseq-fastq-regressions-iwydbmwr/stranded work>)

使用默认 Snakemake rerun triggers 实测：无变动时无需重跑；删除Ballgown文件会重建；修改矩阵脚本的实际代码后会重算矩阵；修改索引副本内容后会重新安排比对；缺少索引在DAG阶段失败。普通模式也会因注释contig不属于BAM而阻断StringTie。

## 其他独立核实

- 原prepDE默认CLI与修复版默认CLI，在三组真实StringTie合成GTF（gene_name、GFF3 ncRNA、分离bundle）上的gene/transcript count逐单元相同。pipeline通过新增选项切换canonical ID和每样本读长，计数公式保持不变。证据：[旧版等价检查](/home/weiyuekuan/tmp/prepde-legacy-equivalence-v1ch6hhl)。
- 真实HISAT2的标准、`--no-mixed`、`--no-discordant`及并用模式，含未比对mate的summary均通过新解析器检查。证据：[summary模式检查](/home/weiyuekuan/tmp/hisat-summary-mixed-3fwyxn_3)。
- 独立验证含空格/单引号目录、cores1/3、索引/脚本依赖；HISAT2内部再次拼接shell的路径问题由安全相对别名修复。证据：[独立流程回归](/home/weiyuekuan/tmp/rnaseq-independent-fix-p071q193)。
- 相对剪接文件 `--known-splicesite-infile "known sites.ss"` 在带单引号workdir中实际比对成功；原workdir路径先解析，再创建临时别名，且列入DAG。证据：[剪接输入检查](/home/weiyuekuan/tmp/rnaseq-fix-20261022-221623/splice-check)。
- 当前默认参考已迁移至 `/home/weiyuekuan/3.hista2_reference/`。FASTA与索引的12条序列名/长度一致，共1,093,552,602 bp；GTF 396,289条记录的序列与坐标全部兼容。新版BED准备与矩阵注释解析均通过：39,686转录本/基因模型，12条染色体。证据：[参考准备报告](/home/weiyuekuan/tmp/rnaseq-fix-20261022-221623/reference-validation/report.json)。这些检查未替代参考来源或逐碱基组装一致性证明。

## 工具版本与边界

本次复用现有工具，没有下载或安装依赖：

- Snakemake 9.19.0。
- HISAT2 2.2.1、StringTie 2.2.3，均使用现有缓存可执行文件，与YAML的这两个版本一致。
- samtools 1.23.1、RSeQC 5.0.4、pysam 0.24.0、gffread 0.12.9。
- 测试/判链Python 3.13；矩阵/Excel子进程Python 3.8.20，openpyxl 3.1.5。矩阵和Excel脚本现不再依赖pandas的自动类型推断。

未新建整套Conda环境，因此不宣称已验证完整YAML依赖求解或全部固定版本组合（YAML中samtools为1.20、Python为3.11.9）。本次仅运行人工数据和参考只读验证，没有重算真实样本，没有修改真实参考文件或历史分析结果。

源代码修改前备份：[source-before](/home/weiyuekuan/tmp/rnaseq-fix-20261022-221623/source-before)。详细修复清单、差异和指纹保存在同级审查目录。

## 重复执行

使用含pysam的Python运行测试，设置已有工具路径和包含openpyxl的量化Python：

```bash
RNASEQ_TEST_FULL_FASTQ=1 \
RNASEQ_TEST_KEEP=1 \
RNASEQ_TEST_SNAKEMAKE=/path/to/snakemake \
RNASEQ_TEST_HISAT2=/path/to/hisat2 \
RNASEQ_TEST_HISAT2_BUILD=/path/to/hisat2-build \
RNASEQ_TEST_SAMTOOLS=/path/to/samtools \
RNASEQ_TEST_STRINGTIE=/path/to/stringtie \
RNASEQ_TEST_RSEQC=/path/to/infer_experiment.py \
RNASEQ_TEST_GFFREAD=/path/to/gffread \
RNASEQ_TEST_WORKFLOW_PYTHON=/path/to/quantification/python \
/path/to/pysam/python -m unittest discover -s tests -v
```

`test_strand_bed.py`和`test_strandness.py`的真实工具测试也要求gffread/infer_experiment.py在PATH。未满足工具条件的测试会显示skip，不能算作验证通过。
