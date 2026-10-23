# 新版历史问题逐项复验（2026-10-23）

结论：前面列出的修复和新增功能均在当前新版中，逐项核查未发现遗漏或残留的已知失败。当前目录：`/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline`。

本轮重新执行 **164 个测试方法，0 失败、0 错误、0 跳过**，用时 145.189 秒。核查表共 33 个条目：最初18项、后续9项、功能与默认值6项；部分条目交叉，不表示33个独立缺陷。77个不同测试方法直接对应这些条目，均已自动核对本轮状态为PASS。默认参考条目另由本轮真实参考只读验证支持。

复测前后，上一轮发布的36个文件逐项SHA256均未变化；测试入口从当前新版目录发现测试并启用真实FASTQ分支。此次没有修改业务源码，本文件是新增验证记录。

## 逐项对应

| 编号 | 历史问题或新增功能 | 本轮结果 | 行为证据 |
|---|---|---|---|
| I01 | gene_id拼接gene_name导致count归零 | 通过 3 个关联方法；[主用例](/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline/tests/test_quantification.py:96) | canonical模式保持literal\|G和数值10；旧拼接ID被拒绝；真实两分支raw与最终count相同且预期基因非零。 |
| I02 | prepDE固定75bp | 通过 4 个关联方法；[主用例](/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline/tests/test_read_lengths_and_resources.py:30) | 全量两端500bases/4reads得到125；149.5小数保留；150.5/75.25读长得到20/40；真实150bp计数对照覆盖度公式且不同于固定75。 |
| I03 | HISAT2索引无依赖 | 通过 2 个关联方法；[主用例](/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline/tests/test_workflow_regressions.py:367) | 缺第8片索引时建图失败且未开始比对；修改索引字节后dry-run必须安排hisat2_align_sort。 |
| I04 | ncRNA_gene丢失 | 通过 1 个关联方法；[主用例](/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline/tests/test_quantification.py:155) | gene/mRNA与ncRNA_gene/ncRNA两类模型均在最终TPM保留G1/G2。 |
| I05 | 重复abundance覆盖 | 通过 2 个关联方法；[主用例](/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline/tests/test_quantification.py:208) | 不重叠同gene行TPM500000+500000=1000000且FPKM=2；重叠重复行拒绝。 |
| I06 | 修改脚本不重算 | 通过 1 个关联方法；[主用例](/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline/tests/test_workflow_regressions.py:256) | 无变更无需重跑；修改矩阵脚本真实代码后必须安排矩阵规则；源码登记包括prepDE等导入模块。 |
| I07 | --cores被固定参数绕过 | 通过 4 个关联方法；[主用例](/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline/tests/test_read_lengths_and_resources.py:69) | 3至64核分配总额及各工具上限；1/2核顺序模式单阶段不超预算；真实1/4核与关闭Ballgown的3核回归验证资源日志。 |
| I08 | Excel ID自动转换 | 通过 2 个关联方法；[主用例](/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline/tests/test_quantification.py:267) | 001、NA、literal\|G、=1+1精确保留为文本类型s，表达数值为n；真实两分支同样检查样本和基因ID。 |
| I09 | gene子特征误认为transcript | 通过 1 个关联方法；[主用例](/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline/tests/test_strand_bed.py:350) | 真实gffread处理带TF_binding_site的gene，断言模型数仍为1且唯一名称为真实transcript。 |
| I10 | GFF转义ID不一致 | 通过 3 个关联方法；[主用例](/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline/tests/test_quantification.py:166) | G%3B1不被解码且矩阵成功；GTF标准化字节不变；真实gffread保留%3B和竖线。 |
| I11 | GTF无模型基因误填零 | 通过 4 个关联方法；[主用例](/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline/tests/test_strand_bed.py:236) | gene-only及缺exon明确拒绝；raw缺基因失败且不写矩阵；真正明确全零模型允许。 |
| I12 | CDS/UTR与exon矛盾 | 通过 2 个关联方法；[主用例](/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline/tests/test_strand_bed.py:290) | CDS及两种UTR越界或跨intron均拒绝，合法子区间通过；错误GFF在转换器运行前拒绝。 |
| I13 | 空summary/abundance及非法数值被补零 | 通过 6 个关联方法；[主用例](/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline/tests/test_alignment_report.py:67) | 空报告、缺失/重复字段、非法百分数/计数、空abundance、NaN、负数、非整数count均实测拒绝。 |
| I14 | 样本/基因集合与TPM未校验 | 通过 5 个关联方法；[主用例](/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline/tests/test_quantification.py:190) | 错/缺/重复样本及缺基因拒绝；TPM500000或表达非零但TPM0拒绝；Excel调换样本顺序或替换ID拒绝且无输出。 |
| I15 | 普通模式缺参考检查 | 通过 4 个关联方法；[主用例](/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline/tests/test_workflow_regressions.py:375) | 普通真实比对后错误contig阻断StringTie且GTF不存在；坏BAI、截断BAM、坐标越界、非坐标排序均拒绝。 |
| I16 | Ballgown副产物不受管理 | 通过 1 个关联方法；[主用例](/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline/tests/test_workflow_regressions.py:256) | 两模式均断言5种ctab存在；删除t_data.ctab后必须调度StringTie重建；5种文件均在同一output和rule all登记。 |
| I17 | 空格/单引号路径失败 | 通过 3 个关联方法；[主用例](/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline/tests/test_workflow_regressions.py:256) | 普通分支在含空格及单引号目录实际FASTQ至Excel成功，链特异分支在含空格目录成功；清单空格路径成功；剪接输入解析至原目录且安全临时别名无空白。 |
| I18 | 配置/文档路径过时 | 通过；[本轮真实参考验证](/home/weiyuekuan/tmp/rnaseq-issue-verification-20261023-123742/reference-recheck/summary.json) | config及README默认参考为/home/weiyuekuan/3.hista2_reference/，当前源码和运行示例目录为/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline。 |
| A01 | HISAT2 extra参数输入/输出别名、二次裁剪跳读与链方向覆盖 | 通过 3 个关联方法；[主用例](/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline/tests/test_read_lengths_and_resources.py:85) | 拒绝受管参数及别名/组合；允许正常评分和read group。当前自动化真实HISAT2端到端使用正常参数，恶意参数集合直接调用解析器；拒绝发生在运行工具之前。 |
| A02 | StringTie extra禁止额外BAM及模式覆盖 | 通过 2 个关联方法；[主用例](/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline/tests/test_workflow_audit.py:24) | 拒绝额外位置BAM、merge/mix、受管参数及输入输出覆盖；解析器和Snakemake建图均断言。 |
| A03 | StringTie --ptf依赖及引号路径 | 通过 2 个关联方法；[主用例](/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline/tests/test_workflow_audit.py:33) | 现套件行为测试确认相对路径按workdir解析、单引号空格保留、缺失文件报错、Snakemake input中出现ptf；本轮新增独立真实StringTie 2.2.3基线/--ptf PATH/--ptf=PATH三次执行均退出0，T1/T2各20，gene40，原合成BAM/参考/PTF哈希不变。 |
| A04 | StringTie省略零估计冗余模型的严格恢复 | 通过 8 个关联方法；[主用例](/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline/tests/test_quant_audit.py:49) | 真实2.2.3产生不同GTF集合后恢复显式零，counts/TPM逐值对照；负数NaNInf、非零缺失、错模型/陈旧coverage、缺ctab均拒绝。Ballgown关闭完整模型真实端到端通过。 |
| A05 | GFF3多gene Parent转换前拒绝 | 通过 1 个关联方法；[主用例](/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline/tests/test_strand_bed.py:259) | Parent=g1,g2明确失败，mock gffread证明转换器未调用，无BED输出。 |
| A06 | gene与子模型染色体/链/坐标冲突 | 通过 4 个关联方法；[主用例](/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline/tests/test_strand_bed.py:243) | GTF/GFF3分别测试gene先后出现、染色体/方向/越界和多子模型，错误GFF3在转换前拒绝；真实CDS-only转换记录来源和无显式exon。 |
| A07 | 输入与输出/日志路径及链接冲突 | 通过 5 个关联方法；[主用例](/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline/tests/test_workflow_audit.py:45) | 保护同路径/软硬链接、跨样本输出、文件父路径；真实Snakemake任务准备前拒绝日志=FASTQ，且原文件字节保持。 |
| A08 | HISAT2 summary可选段标题配对 | 通过 4 个关联方法；[主用例](/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline/tests/test_alignment_report.py:111) | 从有效报告删除各父标题必须失败；no-mixed/no-discordant合法省略仍通过，保持正确分母。 |
| A09 | Excel大整数精度与尺寸限制/失败原子保留 | 通过 3 个关联方法；[主用例](/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline/tests/test_quant_audit.py:175) | 真实openpyxl回读10000000000000001精度损失报错且历史输出字节不变；16385列越界报错不发布。行上限有源码保护但无独立越界测试。 |
| F01 | fastp参数边界/报告和显式polyG开关 | 通过 3 个关联方法；[主用例](/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline/tests/test_fastp.py:109) | 真实fastp核对Q20、30%、N5、40bp边界及-w8；NextSeq命名reads关闭polyG保留80bp、开启60bp。真实raw至Excel、polyG改变后的实际重跑及raw哈希不变。 |
| F02 | clean PE实际平均读长传递prepDE | 通过 5 个关联方法；[主用例](/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline/tests/test_read_lengths_and_resources.py:30) | 两端全量平均125/149.5数值断言；75bp与150bp不同样本实际端到端用独立公式重算counts，并证明150bp样本不同于旧固定75结果。polyG后clean长度/总bases和计数均复核。 |
| F03 | 每样本独立门控及混合FR/RF | 通过 3 个关联方法；[主用例](/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline/tests/test_workflow_integration.py:262) | 真实RSeQC和StringTie检查每样本env/GTF方向、混合方向各自正确；坏样本env/GTF不存在、好样本完成且总体非零、新完整矩阵不生成。 |
| F04 | 坏样本继续与修复复用成功结果 | 通过 2 个关联方法；[主用例](/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline/tests/test_workflow_regressions.py:308) | raw与clean两入口破坏一个样本，实际好样本fastp/BAM/读长/判链/定量完成；修复后矩阵完成，好样本各输出内容与纳秒mtime不变。无CLI keep-going，必须自动profile生效。 |
| F05 | MAPQ默认0配置/回退/CLI及显式30 | 通过 3 个关联方法；[主用例](/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline/tests/test_mapq_defaults.py:23) | 受控SAM+真实SAMtools：线程1/3乘默认/显式0/显式30六组合，0含0、29、30、60及未比对，30仅30/60。Snakemake干跑验证yaml默认、双处缺键回退及覆盖30；RSeQC仍mapq30。默认0真实HISAT2到Excel覆盖在Ballgown关闭测试。 |
| F06 | 线程预算、低线程顺序模式和资源错误传播 | 通过 4 个关联方法；[主用例](/home/weiyuekuan/12.Snakemake/4.RNA-seq-pipeline/tests/test_read_lengths_and_resources.py:69) | 3至64核预算不超额、1/2核顺序分配、上游管道失败RuntimeError。真实1核/4核流程读取资源日志并验证sum/max不超过调度预算。 |

## 补充真实工具实测

- 真实HISAT2构造重复参考及未比对reads，普通/判链两种模式均省略MAPQ配置以使用默认0。每个BAM实际保留低MAPQ记录600条、高MAPQ记录200条、未比对记录40条；两分支均完成StringTie和Excel回读。判链仍使用独立MAPQ30并通过。详见[专项报告](/home/weiyuekuan/tmp/rnaseq-issue-verification-20261023-123742/independent-mapq/REPORT.md)。
- 真实StringTie2.2.3分别执行无ptf基线、`--ptf PATH`、`--ptf=PATH`，含空格和单引号路径均正常；T1/T2估计各20、gene总数40，Ballgown输出齐全，输入哈希不变。详见[原始结果](/home/weiyuekuan/tmp/rnaseq-issue-verification-20261023-123742/ptf-real-recheck/results.json)。
- 当前默认真实注释重新验证：39,686个转录本、201,714条exon记录、12条序列，序列名和坐标界限与HISAT2索引兼容；注释SHA256与上轮一致。此项不证明组装逐碱基来源一致。

## 执行证据与范围

- [每个测试方法的实际状态](/home/weiyuekuan/tmp/rnaseq-issue-verification-20261023-123742/test-results.json)。
- [完整测试日志](/home/weiyuekuan/tmp/rnaseq-issue-verification-20261023-123742/all-tests.log)。
- [问题与完整测试ID的对应表](/home/weiyuekuan/tmp/rnaseq-issue-verification-20261023-123742/issue-verification.json)。
- [复现运行脚本](/home/weiyuekuan/tmp/rnaseq-issue-verification-20261023-123742/run_full_tests.py)。
- [源码快照匹配](/home/weiyuekuan/tmp/rnaseq-issue-verification-20261023-123742/source-match.json)。
- [本轮参考结果](/home/weiyuekuan/tmp/rnaseq-issue-verification-20261023-123742/reference-recheck/summary.json)、[实际软件版本](/home/weiyuekuan/tmp/rnaseq-issue-verification-20261023-123742/runtime-versions.json)。

使用本机已有软件与人工FASTQ/BAM，没有安装或下载软件，没有处理真实测序样本。主流程使用固定Python3.11.9/HISAT2 2.2.1/SAMtools1.20/StringTie2.2.3；参考检查和部分辅助测试使用现有Python3.13.9+pysam0.24.0。没有宣称三套独立YAML环境均在本轮重建或完整经`--use-conda`执行。

仍保留的测试粒度边界：脚本变化重跑以矩阵脚本为实际变更代表，其他导入依赖逐项读源码核实；5个ctab均验证生成，删除重建抽测t_data.ctab；ptf已验证真实执行与DAG输入登记，未单独在完整运行后修改ptf再测调度；Excel列上限实测拒绝，行上限有源码保护但未构造百万行越界数据。这些边界不作为已完成实测宣称。

上述结论针对本次对话已列问题及当前用例覆盖范围，不代表对所有物种、注释和真实数据的无条件保证。
