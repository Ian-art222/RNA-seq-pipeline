# GTF 参考解析修复与验证（2026-10-23）

本次基线为 `78fdab631284c9e869ecd0e75f180f6fa4c36fd8`。用户指出的 P1 缺陷成立，且实际默认 CtZJ 参考也会触发。此前 164 项测试通过，但没有覆盖该参考从前置处理到 prepDE 的完整链；旧验证报告已补充更正。

## 原始缺陷与影响

`prepare_strand_bed.read_gtf_models()` 接受合法 exon-only GTF 及附带 CDS/UTR 的 GTF，前置处理随后直接复制源 GTF。`prepDE.parse_gtf(reference=True)` 却复用了严格的 StringTie 输出解析逻辑，要求 transcript 行、拒绝 CDS/UTR/codon。

在基线源码上实际复现：

| 参考结构 | 前置处理 | prepDE 参考解析 |
|---|---|---|
| transcript + exon | 通过 | 通过 |
| gene + transcript + exon + UTR/CDS/codon | 通过 | `unexpected StringTie feature` |
| transcript + exon + CDS | 通过 | `unexpected StringTie feature 'CDS'` |
| exon-only | 通过 | `no quantified transcript models` |

默认 `/home/weiyuekuan/3.hista2_reference/CtZJ_T2T_gene.gtf` 有 201,714 条 exon、194,575 条 CDS，无 transcript 行。旧 prepDE 在第 3 行 CDS 报错。源文件 SHA256 为 `03da9ef7f0cc75b5c2a6ecca225ff936e468d05f771db57d89edbe061da13c87`，修复测试前后不变。此缺陷造成汇总阻断，不是本次已证明的静默错误计数。

StringTie 官方说明：exon 足以定义转录本，GTF2 不要求显式 transcript 行，允许附带 CDS/UTR/codon。[格式说明](https://ccb.jhu.edu/software/stringtie/gff.shtml)

## 修复

- 参考解析与 StringTie 定量输出解析分开；参考解析共用 `read_gtf_models()`，支持 exon-only 和合法 CDS/UTR/codon 元数据。
- 用共同的模型生成 BED12 与规范 transcript/exon GTF，再将同一 GTF 传给 StringTie 和 prepDE。转录本范围和长度来自 exon，原始参考不修改。
- 保留 gene_id、transcript_id、gene_name；非法 CDS/UTR/codon 不得落到原 exon 之外。无 exon 的 GTF、gene_id 跨染色体/链、StringTie 无法保留的引号/分号等属性提前拒绝。
- StringTie 输出仍严格检查命令记录、coverage、非有限值、重复/缺失模型与参考一致性，未因支持普通 GTF 而放松。
- fastp、clean 读长统计和 HISAT2 声明参考预检前置依赖；采用 `ancient()` 表示运行顺序约束，参考更新时间本身不强迫已成功的清洗、读长统计和比对重做。定量/矩阵仍按正常参考依赖重算。prepDE 明确追踪共享解析模块。
- GTF 规范化副本不再保证与原文件字节相同，保证的是模型与 ID 一致；原文件仍逐字节保留。旧的对应测试已按这个契约更新。

## StringTie 短内含子边界

实际固定版本 StringTie 2.2.3 将 1–4 bp exon 间隙合并，5 bp 起保留。该行为由源码的 `mergeCloseExons=true` 和 `GFF_MIN_INTRON=4` 确认：[加载参数](https://github.com/gpertea/stringtie/blob/v2.2.3/stringtie.cpp#L353-L357)、[合并条件](https://github.com/gpertea/stringtie/blob/v2.2.3/gclib/gff.cpp#L2162-L2163)。本次还实际测试了 0–10 bp 间隙。

为使 BED、StringTie 与计数使用同一模型，规范化显式采用这一既有行为，在报告 `short_intron_merges` 中逐项记录，并记录 `short_intron_added_exonic_bases`、输出警告。所有子特征先对原 exon 校验，防止合并掩盖错误注释。这是工具的模型兼容处理，不是对原始生物学注释作出修订。

CtZJ 原参考中恰有以下 5 处 4 bp 间隙，共 20 bp；其他 39,681 个模型无此调整：

| transcript_id | 染色体 | 间隙坐标（1-based inclusive） | 链 |
|---|---|---|---|
| CtZJ02G027260.1 | CtZJ02 | 88443311–88443314 | + |
| CtZJ04G008540.1 | CtZJ04 | 20159408–20159411 | + |
| CtZJ06G026850.1 | CtZJ06 | 73915440–73915443 | - |
| CtZJ08G032140.1 | CtZJ08 | 62031811–62031814 | + |
| CtZJ12G021530.1 | CtZJ12 | 79054266–79054269 | - |

外部 BED 也须与规范化后的模型一致；仍保留上述短间隙的外部 BED 会因 exon blocks 不一致被拒绝。

## 回归与实际工具验证

最终完整回归 **177 项全部通过，0 失败、0 错误、0 跳过**，用时 196.897 秒。25 个源码/测试文件与测试时 SHA256 一致；逐项记录与日志位于 `/home/weiyuekuan/tmp/rnaseq-gtf-fix-20261023-150341/`。

新增 [test_reference_pipeline.py](test_reference_pipeline.py) 共 13 项方法，其中包含：

- 四种参考结构在前置处理、BED、prepDE 中使用相同模型；双链、多 exon、literal ID、gene_name 优先级。
- 重叠/重复/相邻 exon、1/4/5 bp 边界、非法子特征与 gene_id 位置信息冲突、严格 StringTie 输出解析。
- 真实 HISAT2 → StringTie → prepDE → CSV → Excel：完整 GTF 与 exon-only 各在普通/链特异模式运行，实际回读 Excel。
- 错误参考分别在 raw/clean 起点阻断，未启动 fastp/读长统计/HISAT2，原 FASTQ 不变。
- 真实 raw → fastp → Excel 后仅更新参考：重新检查参考，15 个 raw/clean/读长/BAM/BAI 文件内容及时间均不变，表达结果与 Excel 验证通过。

另使用完整真实 CtZJ 注释和人工未比对 BAM，实际执行 prepare → StringTie 2.2.3 → prepDE → 矩阵 → Excel。三者的 39,686 个转录本模型逐项相同；gene/transcript counts 各 39,686 行，均为该人工输入预期的 0。count/FPKM/TPM CSV 各 39,686 个基因，ID 集合与顺序一致。实际回读 Excel，三张表达表各 39,687 行（含表头）× 2 列，ID 为文本，数值为 0。该检查不代表任何真实生物样本的表达结果；Excel 的 alignment 表明确使用人工 BAM 的实际记录统计，并非虚构 HISAT2 汇总。

全参考证据：`/home/weiyuekuan/tmp/rnaseq-gtf-independent-20261023-150126/repaired-full-ctzj/results.json`、`final-excel-results.json`。最终增加 gene 一致性约束后重新生成的 GTF/BED 与此工具链输入逐字节相同，记录在主审计目录 `final-reference/verification.json`。

实际软件沿用已安装的固定主环境：HISAT2 2.2.1、samtools 1.20、StringTie 2.2.3、pandas 2.2.2、openpyxl 3.1.5；另用现有 fastp 1.3.6、RSeQC 5.0.4、gffread 0.12.9、Snakemake 9.19.0。未重新安装或下载环境，不宣称三份 YAML 环境均已从头创建。

第一次完整回归在人工 BAM 起点夹具中发现缺少新参考前置产物，导致尝试用其占位索引重比对；主动终止并保留日志。依赖改为时间不触发重算的前置约束，BAM 夹具按真实流程顺序在发布 BAM 前先准备参考；另用真正 raw 数据验证参考更新不重比对。最终完整运行结果以 `test-results.json` 为准。

## 审查表其余条目的处理状态

以下已核实，但本次没有实现相应新功能，也不能据此宣称整张审查表全部解决：

| 条目 | 本轮确认 | 状态 |
|---|---|---|
| raw→clean→HISAT2 质量编码 | 真实 fastp：人工 Phred+64 Q0 输入未声明编码时全保留，显式 --phred64 后全过滤 | 支持缺口保留；README 明确默认 +33。clean-only 可手动 HISAT2 --phred64 |
| 判链最少有效证据 | 真实 RSeQC 仅 1 条 usable alignment 仍 PASS/--fr | 门控策略不足保留，未擅定科学阈值 |
| 混合建库队列 | 全局 stranded 不能逐样本混合 unstranded/stranded | 支持范围限制保留；逐样本 FR/RF 仍支持 |
| 样本科学 QC | 结构有效全零矩阵允许成功 | 尚无统一科学 QC 状态 |
| 历史矩阵残留 | 上游失败时旧文件可保留 | 尚无当前运行 success manifest，需检查退出码和日志 |
| BAI 坐标范围 | 实际 600,000,000 位置 alignment：BAI 失败，CSI 成功 | workflow 仍固定 BAI |
| 内存/磁盘/I/O 调度 | 已有 CPU 线程预算，无对应资源配额 | 工程优化未实现 |
| FASTQ/Excel 性能 | 全量配对校验/读长统计有额外 I/O；Excel 写入全量载入，回读已 read_only | 性能优化未实现 |

余项复现记录：`/home/weiyuekuan/tmp/rnaseq-review-other-rjs4k3u4/results.json`。没有证据表明用户实际 FASTQ 是 +64，也没有重跑真实生物样本或改写其已有表达矩阵。
