# 样本独立执行验证（2026-10-23）

按本次要求，将原有全队列判链门控改为逐样本门控。各样本的计算、方向选择和输出相互独立；共享参考准备和最终矩阵仍保留必要依赖。

## 实现

- `workflow/Snakefile`：`strandness_gate` 只读取当前样本 JSON，生成 `strandness/{sample}.env` 和独立诊断日志；StringTie 只依赖自己的 env。
- `workflow/scripts/strandness.py summarize`：单独汇总 PASS/FAIL 和原始证据，不产生 env，不作为单样本定量的前置依赖。合格 FR/RF 混合时记录各自方向并警告，不再要求同一批样本方向完全相同。
- `workflow/profiles/default/config.yaml`：默认 `keep-going: true`，独立任务在其他任务失败后继续；有失败时整体仍返回非零退出码。
- 独立调用 `strandness.py gate` 的旧多样本严格校验仍保留，本 workflow 每次仅传一个样本。

## 实际验证

使用人工 FASTQ/BAM 和小型参考，未修改真实数据或参考。沿用已有软件：Snakemake 9.19.0、fastp 1.3.6、HISAT2 2.2.1、samtools 1.20、StringTie 2.2.3、RSeQC 5.0.4、gffread 0.12.9。没有安装或下载软件；本轮没有重新创建各 conda 环境。

| 场景 | 结果 |
|---|---|
| 同批 FR/RF 合格样本 | 分别使用自己的 `--fr` / `--rf`，两样本矩阵和 Excel 完成 |
| 合格样本 + 模糊判链样本 | 好样本完成定量；坏样本无放行 env/GTF；整体失败，不生成新的完整矩阵 |
| 修改一个已完成样本至判链失败 | 其他样本 GTF、abundance、判链 JSON/env 的内容和纳秒时间戳均不变 |
| clean FASTQ 中一个样本记录损坏 | 坏样本在读长校验处失败，好样本实际完成 HISAT2 和 StringTie |
| raw FASTQ 中一个样本记录损坏 | 坏样本在 fastp 前校验失败，好样本完成 fastp、比对、判链和定量 |
| 修复上述 raw/clean 坏文件并续跑 | 完整矩阵及 Excel 恢复；好样本 BAM、读长、GTF、abundance 等内容和时间戳不变 |
| 不显式传 `--keep-going` | 外部运行目录下自动加载源码内默认 profile，仍完成上述故障隔离 |
| 伪造 PASS、损坏报告、输入路径别名 | 证据继续严格校验；不会绕过样本门控或覆盖输入 |

执行记录：首次全套运行 162 个测试方法，161 个通过；新增 FASTQ 故障测试的两个子场景仅因预期错误文字与实际 `invalid/truncated FASTQ record` 不一致而失败。修正该测试断言后，单独重跑此方法，raw/clean 故障与恢复场景全部通过（1/1，34.145 秒）。业务代码没有为此更改。合计 162 个不同测试方法均有通过记录，没有跳过；没有将首轮失败日志标作全部通过。

证据目录：`/home/weiyuekuan/tmp/rnaseq-sample-isolation-20261023/`。

- `source-before/`：改动前源码备份。
- `source-changes.diff`：本轮变更。
- `all-tests.log`：首次完整套件记录（177.542 秒）。
- `isolation-repair-test.log`：修正断言后的 raw/clean 隔离和恢复验证。
- `run_full_tests.py`、`run_isolation_repair_test.py`、`test-tool-paths.json`：复现入口与软件路径。
- 最新故障恢复人工数据：`/home/weiyuekuan/tmp/rnaseq-fastq-regressions-y_in0uj3/`。
- 最新判链独立运行人工数据：`/home/weiyuekuan/tmp/rnaseq-strand-integration-m4h0upja/`。

## 边界

默认严格扫描仍会在启动前拒绝缺失/空 FASTQ、错误样本名、错误配置和输入输出冲突；共享参考故障会影响依赖它的样本。运行期的单样本错误才由独立任务和 keep-going 隔离。

表达矩阵和 Excel 要求全部纳入样本成功，不自动剔除失败样本。失败重跑可能保留历史完整矩阵，文件存在不能证明本次成功，须查看退出码和日志。只修复工具环境时，已有失败判链 JSON 需要按 README 的 `--forcerun` 示例重新生成。

进程、文件与任务依赖独立不等于机器资源硬隔离。CPU 预算由 Snakemake 管理，内存和磁盘仍共享。混合链方向会记录警告，需结合实际建库信息核查。
