# RNA-seq Snakemake pipeline

2026-10-23 恢复版：在 GitHub 原版基础上恢复本次对话中的全部修复、fastp 清洗及样本独立执行改进，源码移至 `/home/weiyuekuan/12.Snakemake/RNA-seq-pipeline`。比对 BAM 的默认 `alignment.mapq_filter` 改为 `0`，RSeQC 判链的独立阈值 `strandness.mapq` 仍为 `30`。本次恢复与复验记录见 [tests/RESTORED_VALIDATION_20261023.md](tests/RESTORED_VALIDATION_20261023.md)。

本轮完整回归 164 项全部通过，0 失败、0 错误、0 跳过；使用人工 FASTQ/BAM 和现有软件环境验证，真实参考另做只读兼容检查。

双端 raw FASTQ → 可选 fastp 清洗 → clean FASTQ → HISAT2 → 排序/索引 BAM → 参考兼容检查 → StringTie `-G -e` → read count 估计、FPKM、TPM 和 Excel。默认从已有 clean FASTQ 开始；可在配置中启用 fastp 并明确指定是否修剪 polyG。输入必须属于同一个参考基因组；不进行新转录本发现或差异表达统计。

2026-10-23 样本独立执行更新：已将判链改为逐样本放行，默认继续执行其他无失败依赖的任务；故障隔离、修复续跑及结果复用的验证见 [tests/SAMPLE_ISOLATION_20261023.md](tests/SAMPLE_ISOLATION_20261023.md)。

## 2026-10-23 全面复审

本轮修复比对参数绕过、StringTie额外BAM输入/点特征依赖、冗余零转录本汇总、GFF3多归属与gene声明冲突、全DAG输入输出保护、损坏summary及Excel整数校验。155项测试全部通过；主流程固定版本环境已离线创建并用于验证。详细证据和未覆盖的环境/科学边界见 [tests/FULL_AUDIT_20261023.md](tests/FULL_AUDIT_20261023.md)。

## 2026-10-22 修复

- 从每个样本两端完整 clean FASTQ 流式统计实际平均 read 长度，分别用于 prepDE，取消 pipeline 隐含的固定 75 bp。计算同时检查 FASTQ 结构、两端记录数及 read ID 配对。
- pipeline 的 prepDE 输出使用原始 `gene_id`，保留 ID 中的 `|`，不会拼接 gene_name，也不靠截断字符串猜测映射。
- 两种模式均标准化并验证注释，且验证 BAM/BAI 与注释序列和坐标。缺少 exon 模型的基因不能作为零表达进入矩阵。
- 支持 GFF3 的 ncRNA_gene 等类型；忽略不属于转录本的普通 gene 子特征；一致保留转义 ID。
- 汇总严格检查基因/样本集合、数值、TPM 总和；合法的同基因、同链、不同区间 abundance 记录聚合，冲突或重复重叠记录报错。
- 索引、辅助脚本、额外剪接位点输入和 Ballgown 文件列入依赖图，修改或丢失后可触发重建。
- 实际线程分配服从 Snakemake 预算，路径引用和 HISAT2 包装器内部路径处理支持空格及单引号。
- Excel 中基因 ID、样本 ID、列名固定为文本，保留 `001`、`NA` 和 `|`。导出后重新打开检查表结构和数值。
- 默认参考位置已更新为当前 `/home/weiyuekuan/3.hista2_reference/`，运行示例改用本源码目录。

完整修复测试见 [tests/VALIDATION.md](tests/VALIDATION.md)。旧版本已经生成的矩阵不会在此次代码修复中自动改写；需要在选定 workdir 中检查 dry-run 并执行相应重算。

## 源码和运行目录

源码目录：`/home/weiyuekuan/12.Snakemake/RNA-seq-pipeline`。源码包含 workflow、config、tests；数据和分析结果放在独立的运行目录，用 `--directory` 指定。

```text
workdir/
├── rawdata/                         # fastp.enabled=true 时的输入
│   ├── sampleA/
│   │   ├── sampleA_1.fastq.gz
│   │   └── sampleA_2.fastq.gz
│   └── sampleB/
│       ├── sampleB_1.fastq.gz
│       └── sampleB_2.fastq.gz
├── cleandata/
│   ├── sampleA/
│   │   ├── sampleA_1.clean.fastq.gz
│   │   └── sampleA_2.clean.fastq.gz
│   └── sampleB/
│       ├── sampleB_1.clean.fastq.gz
│       └── sampleB_2.clean.fastq.gz
├── qc/fastp/                        # 启用 fastp 时生成 HTML/JSON
├── alignment/
├── quantification/
├── logs/
└── .snakemake/
```

`fastp.enabled: false` 时仅扫描 cleandata 下一级样本子目录；设为 true 时改为扫描 rawdata 下一级样本子目录，并自动生成 cleandata，无需预先创建 clean 文件。目录名就是样本名，允许字母、数字、下划线、点、连字符，首字符不能是点或连字符。每个样本必须有一对非空 FASTQ；多个 lane 应先按相应样本正确合并。两种扫描分别使用 `sample_scan` 和 `raw_sample_scan` 的后缀及 `strict` 设置。`strict: true` 会拒绝不完整子目录；设为 false 时明确报告并跳过。普通文件会被忽略，平铺 FASTQ 需先放入对应样本子目录。

FASTQ 支持标准四行记录、gzip 或纯文本（按后缀选择解压）；R1/R2 后缀可配置。每个 read 的序列和质量长度必须一致，配对 ID 可使用 `/1`、`/2` 或 Illumina mate 字段。全量读长统计使用常量内存，但会增加一次读取 FASTQ 的 I/O。

## 配置

修改 `config/config.yaml` 或使用额外配置文件。相对路径都相对于 `--directory` 的 workdir。

```yaml
stranded: false
fastp:
  enabled: false
  trim_poly_g: false
  qualified_quality_phred: 20
  unqualified_percent_limit: 30
  n_base_limit: 5
  length_required: 40
paths:
  raw_dir: rawdata
  qc_dir: qc/fastp
  clean_dir: cleandata
  align_dir: alignment
  quant_dir: quantification
  logs_dir: logs
reference:
  hisat2_index_prefix: /home/weiyuekuan/3.hista2_reference/CtZJ_T2T_index
  annotation: /home/weiyuekuan/3.hista2_reference/CtZJ_T2T_gene.gtf
alignment:
  mapq_filter: 0
  hisat2_extra: ""
stringtie:
  use_ballgown: true
  extra: ""
raw_sample_scan:
  r1_suffix: _1.fastq.gz
  r2_suffix: _2.fastq.gz
  strict: true
sample_scan:
  r1_suffix: _1.clean.fastq.gz
  r2_suffix: _2.clean.fastq.gz
  strict: true
threads:
  fastp: 8
```

### fastp 清洗与 polyG 开关

将 `fastp.enabled` 改为 true 后，每个样本运行 `--detect_adapter_for_pe -q 20 -u 30 -n 5 -l 40`，生成成对 clean FASTQ 及 HTML、JSON 报告。以上质量阈值均可在 `fastp` 配置项中修改，`threads.fastp` 默认为8；传给 `-w` 的工作线程数不超过 Snakemake 分配的预算，fastp 自身还会使用 I/O 辅助线程。输出名称继续使用 `sample_scan` 后缀，例如 `cleandata/sampleA/sampleA_1.clean.fastq.gz`，后续规则自动使用这些文件。

- `trim_poly_g: true`：显式传入 `--trim_poly_g`。
- `trim_poly_g: false`：显式传入 `--disable_trim_poly_g`，确保关闭。仅省略 `--trim_poly_g` 并不保证关闭，因为 fastp 可以根据 FASTQ 仪器 ID 自动启用；参见 [fastp 官方 polyG 说明](https://github.com/OpenGene/fastp#polyg-tail-trimming)。

polyG 设置在本次运行的所有样本间统一。如果样本来自不同平台且需要不同策略，应按策略分别建立运行目录与配置。raw、clean 和 QC 输出目录必须互不重叠，输入与输出也不得链接到同一文件。

`-q 20 -u 30` 用来判断整条 read 的低质量碱基比例；`-n 5` 和 `-l 40` 分别过滤 N 过多和长度不足的 reads。这些参数不等同于按质量窗口剪切末端碱基。本流程保留给定命令的行为，不额外启用 `--cut_front`、`--cut_tail` 或 `--cut_right`；参见 [fastp 官方质量修剪说明](https://github.com/OpenGene/fastp#per-read-cutting-by-quality-score)。

清洗前会完整校验 raw FASTQ 的结构、read ID 和两端记录数，避免过滤掩盖原始配对错误。fastp 成功且报告读数与输入一致后才发布 clean FASTQ 和报告；若全部 reads 被过滤，保留质控报告并在后续读长检查处停止，不进入比对。

已有 clean 数据继续保持 `fastp.enabled: false` 即可；该模式不要求 rawdata 或 fastp 环境。切换清洗策略前，应检查 dry-run 中列出的重算范围。

索引必须是该参考的完整、非空 8 个 `.ht2` 或 8 个 `.ht2l` 文件，不能把前缀写成目录，也不能混放两种格式。不存在、不完整或混合索引在建立 DAG 时立即报错。参考注释必须存在且非空；GTF 需要显式 exon 模型，GFF3 会经 gffread 转换后检查。有 gene 声明时，其染色体、链和范围必须与子模型一致；不能无损表示为单一 gene_id 的 GFF3 多 gene Parent 会在转换前被拒绝，避免静默丢失基因。注释与索引必须来自同一组装；序列名/长度/坐标通过只能说明结构兼容，不能单独证明组装来源一致。

`alignment.mapq_filter` 默认为 `0`，不按 MAPQ 过滤比对 BAM；需要时可设为正整数，例如 `30` 保留 MAPQ≥30 的记录。该设置会影响多重比对相关的基因表达；HISAT2 summary 是过滤前统计，不代表最终保留 BAM 的有效定量比例。RSeQC 使用独立的 `strandness.mapq`（默认 `30`）选择判链证据，不会据此裁剪用于 StringTie 定量的 BAM。

`hisat2_extra` 和 `stringtie.extra` 只用于未被流程管理的选项。不能覆盖线程、参考、输入、输出、定量模式或链方向；不支持改变 summary 格式以及在 HISAT2 extra 中另行生成未被规则追踪的输出。`--known-splicesite-infile`、`--novel-splicesite-infile` 支持相对/绝对路径和空格，按 workdir 解析并登记依赖。HISAT2 extra 禁止再次裁剪（`-3/-5/--trim3/--trim5`）、跳读/抽样和手动 `--rna-strandness`，使比对输入与统计 clean 读长保持一致。StringTie extra 禁止添加额外 BAM 或切换模式；支持的 `--ptf` 点特征文件会按 workdir 解析并列入依赖，修改后重跑定量。额外参数先按 shell 风格解析，再安全传递，不作为任意 shell 命令执行。

## 两种定量模式

- `stranded: false`：不执行判链，不给 StringTie 添加 `--fr/--rf`。仍执行 FASTQ、参考注释和 BAM 兼容检查。
- `stranded: true`：每个样本独立执行 RSeQC 和质量门控，达标后立即按自己的方向进行链特异定量；失败只阻断该样本，不会退回普通模式。

两条分支共享未预设链方向的 HISAT2 `--dta` 比对。开启自动判链时，HISAT2 配置和已有 BAM 的命令记录不得带 `--rna-strandness`。StringTie 方向由流程统一管理。

```yaml
strandness:
  sample_size: 200000
  mapq: 30
  min_determinable: 0.50
  min_dominance: 0.80
  bed12: ""
```

`bed12` 留空时从注释生成。提供外部 BED12 时必须和同一注释的转录本 ID、位置、链及 exon blocks 一致。

- mode1：`1++,1--,2+-,2-+`，优势方向对应 fr-secondstrand、StringTie `--fr`。
- mode2：`1+-,1-+,2++,2--`，优势方向对应 fr-firststrand、StringTie `--rf`。
- 每个样本的放行要求：`mode1 + mode2 >= 0.50` 且 `max(mode1,mode2)/(mode1+mode2) >= 0.80`。

这些是本流程阈值，不是 RSeQC 通用标准。记录实际 usable reads、原始比例、工具版本和命令；无采样、非双端、字段缺失/重复、非法数值或证据与记录方向冲突均阻断该样本。少于目标但大于零会警告，不额外加入最低读数阈值。RSeQC 5.0.4 按 BAM 顺序取前 N 条合格 alignment，PE 两端分别计数；它用转录本完整区间（包含 intron）判断重叠，不按 exon blocks 采样，也不是全基因组随机抽样。这些边界写入判链 JSON；比例不能解释成全部 reads 的注释率。

合格样本的 FR/RF 方向可以不同，分别使用 `--fr` / `--rf`。汇总 TSV 保留每个样本的原始比例、方向和状态，混合方向会在汇总日志中警告；应结合实际建库信息核查。该汇总不控制单样本放行。命令行脚本 `strandness.py gate` 保留多样本严格一致性校验供独立调用，但本 workflow 每次只向它传入一个样本。

## 样本独立执行

每个样本有独立的 fastp、读长检查、比对、索引、BAM 检查、判链和 StringTie 任务，以及各自的输出、日志和临时文件。一个样本完成前置步骤后即可继续，不必等其他样本；每一步由单独的作业进程执行，比对步骤还包含 HISAT2 与 samtools 子进程。

`workflow/profiles/default/config.yaml` 默认启用 `keep-going: true`。运行时某个样本失败，其他没有失败依赖的任务仍继续，整体运行最终返回非零退出码。启动示例也显式使用 `--keep-going`，以免自定义 profile 覆盖默认设置。修复失败原因后在同一运行目录重跑，已成功且输入、参数和代码未变化的任务可以复用。

样本独立执行不改变以下依赖：

- 启动前仍检查配置、样本名、输入输出冲突及配对文件是否存在且非空。默认 `strict: true` 下缺失 FASTQ 会阻止建立任务图；文件内容错误则由各样本任务独立发现。
- 参考注释和索引是共同输入，共享参考错误会影响依赖它的样本。
- 最终表达矩阵和 Excel 等汇总要求全部纳入样本成功，不会静默丢弃失败样本。失败运行不会发布新的完整汇总；旧运行留下的文件不能证明本次成功，应检查退出码与日志。
- 各作业共享机器的 CPU、内存和磁盘；`--cores` 控制线程调度预算，并不是操作系统级内存或磁盘隔离。

## 读长与矩阵的含义

每个样本的平均单条 read 长度为 `两端总碱基数 / 两端总 read 数`，保留小数；它不是 paired-end 片段长度。启用 fastp 时使用清洗后 clean FASTQ 的实际读长，因此接头或 polyG 修剪产生的读长变化会进入 prepDE 计数换算。记录位于 `alignment/read_lengths/{sample}.json`，汇总为 `alignment/read_lengths.json`，例如：

```json
{"sampleA": 149.35, "sampleB": 74.8}
```

pipeline 的 prepDE 调用示例（默认普通模式、在运行目录执行）：

```bash
python /home/weiyuekuan/12.Snakemake/RNA-seq-pipeline/workflow/scripts/prepDE.py \
  -i quantification/sample_lst.txt --read-lengths alignment/read_lengths.json --gene-id-only \
  --reference quantification/reference/reference.gtf --ballgown \
  -g quantification/gene_count_matrix.raw.csv -t quantification/transcript_count_matrix.csv
```

每转录本计数公式仍是 `ceil(coverage × exon_length / sample_mean_read_length)`，再合并至基因。这是 StringTie 覆盖度推算的 read counts，不是直接计数 fragments。每个样本均与规范参考完整模型比对；非有限或负覆盖度会报错。StringTie 2.2.3 可能省略零估计的冗余转录本：默认 `use_ballgown: true` 时，仅在同次运行 `t_data.ctab` 的 ID、基因、坐标、链、外显子数和长度全部匹配，且缺失模型的 coverage/FPKM 明确为零时恢复零估计模型，并记录名单。已有 GTF 的 coverage 也必须与表一致。非零缺失、错误表或缺少证据均阻断；`use_ballgown: false` 遇到此类缺失时会提示启用 Ballgown 或检查参考。这里的零是 StringTie 的估计分配，不证明转录本生物学上不存在。

独立调用 bundled prepDE 时保留原 `-l` 默认75和 gene symbol 拼接行为，新增 `--read-lengths`、`--gene-id-only` 才启用本流程行为。支持含空格的 GTF 路径；不应把修复前的 `gene_id|gene_name` raw 表直接交给新的严格汇总器，应重跑 prepDE。

FPKM/TPM 来自 StringTie `-A`。正常非零样本 TPM 总和检查约为100万；全零样本必须 count/FPKM/TPM 一致为零。缺少记录、错误 ID 映射或无法解释的额外模型都不会自动补零。未表达但有合法模型的基因由 StringTie 明确输出零值。

## 线程和磁盘

`threads.hisat2` 是 HISAT2 总线程；`samtools_view`、`samtools_sort` 按 samtools `-@` 表示附加线程，调度预算另计这两个进程的主线程。`samtools_index`、`stringtie` 为总线程。

比对流程根据 Snakemake 实际分配的 `{threads}` 限制各工具。当获得至少3个线程时使用流水线；只有1或2个线程时顺序执行 HISAT2、view、sort，暂存 SAM/BAM，避免三个进程同时超出预算。低线程模式需要更多临时磁盘。临时文件放在样本 alignment 目录下，步骤退出后清理；只有成功时发布最终 BAM/summary。日志中记录实际资源分配和原始输入路径。

## 环境与运行

Snakemake 使用单独环境；通过 `--use-conda` 管理业务环境：

- `workflow/envs/rnaseq.yaml`：HISAT2 2.2.1、samtools 1.20、StringTie 2.2.3、Python 和表格依赖。
- `workflow/envs/strandness.yaml`：Python、pysam、gffread、RSeQC 5.0.4。**普通模式也使用此环境完成共用注释/BAM检查**，但不会运行 RSeQC 判链。
- `workflow/envs/fastp.yaml`：Python 3.11.9、fastp 1.3.6；仅在 `fastp.enabled: true` 时使用。

以下是 Linux/bash 示例。先把 `RNASEQ_WORKDIR` 改为自己的运行目录，并按所选模式准备其中的 cleandata 或 rawdata；不要把源码目录作为运行目录。要使用 rawdata 清洗，先在配置文件中设置 `fastp.enabled: true`，同时设置 `fastp.trim_poly_g`。

```bash
RNASEQ_SOURCE=/home/weiyuekuan/12.Snakemake/RNA-seq-pipeline
RNASEQ_WORKDIR=/path/to/rnaseq-workdir

snakemake \
  --snakefile "$RNASEQ_SOURCE/workflow/Snakefile" \
  --configfile "$RNASEQ_SOURCE/config/config.yaml" \
  --directory "$RNASEQ_WORKDIR" \
  --use-conda --cores 32 --keep-going --dry-run --printshellcmds \
  --config stranded=true
```

检查样本、参考、运行分支和重跑范围后，去掉 `--dry-run` 正式执行。普通模式改为 `--config stranded=false`。中断后可加 `--rerun-incomplete`；运行目录迁移或已有结果升级时先检查 dry-run，不使用 `--touch` 将未核实的旧文件标记为新结果。

首次缺少环境时会安装依赖，后续复用缓存。本工作区要求下载/安装先使用临时取消代理的直接连接，连接失败才使用 `http://127.0.0.1:17897`，不修改永久代理环境。可先单独准备环境：

```bash
env -u HTTP_PROXY -u HTTPS_PROXY -u ALL_PROXY \
  -u http_proxy -u https_proxy -u all_proxy \
  snakemake --snakefile "$RNASEQ_SOURCE/workflow/Snakefile" \
  --configfile "$RNASEQ_SOURCE/config/config.yaml" \
  --directory "$RNASEQ_WORKDIR" --use-conda --cores 2 --conda-create-envs-only
```

如果上述命令因网络连接失败，再用该命令级代理重试：

```bash
env -u ALL_PROXY -u all_proxy -u http_proxy -u https_proxy \
  HTTP_PROXY=http://127.0.0.1:17897 HTTPS_PROXY=http://127.0.0.1:17897 \
  snakemake --snakefile "$RNASEQ_SOURCE/workflow/Snakefile" \
  --configfile "$RNASEQ_SOURCE/config/config.yaml" \
  --directory "$RNASEQ_WORKDIR" --use-conda --cores 2 --conda-create-envs-only
```

索引文件、FASTQ、注释、脚本和相关参数变化按依赖触发重算。若只修复判链工具环境而之前已生成 FAIL JSON，需强制重新判链：对单个样本可在运行命令中加 `--forcerun "$RNASEQ_WORKDIR/quantification/stranded/strandness/sampleA.json"`（替换样本名和自定义 quant_dir）；对全部样本可加 `--forcerun infer_strandness strandness_gate`。修复前的结果如缺少完整 provenance，不能仅凭文件存在认定可复用。

## 输出

普通模式的定量根目录为 `paths.quant_dir`；链特异模式在其下加 `stranded/`。Excel 同样在配置输出目录下加 stranded 子目录。

```text
cleandata/{sample}/{sample}_1.clean.fastq.gz  # 启用 fastp 时生成
cleandata/{sample}/{sample}_2.clean.fastq.gz
qc/fastp/{sample}/{sample}.fastp.html
qc/fastp/{sample}/{sample}.fastp.json
logs/{sample}/fastp.log

alignment/{sample}/{sample}.sorted.bam[.bai]
alignment/{sample}/{sample}.align_summary.txt
alignment/alignment_report.csv
alignment/read_lengths/{sample}.json
alignment/read_lengths.json

quantification/[stranded/]reference/reference.gtf
quantification/[stranded/]reference/reference.bed12
quantification/[stranded/]reference/bam_validation/{sample}.json
quantification/[stranded/]{sample}/{sample}.gtf
quantification/[stranded/]{sample}/{sample}_abundance.tab
quantification/[stranded/]{sample}/*.ctab        # use_ballgown=true
quantification/[stranded/]sample_lst.txt
quantification/[stranded/]gene_count_matrix.raw.csv
quantification/[stranded/]transcript_count_matrix.csv
quantification/[stranded/]gene_count_matrix.csv
quantification/[stranded/]gene_fpkm_matrix.csv
quantification/[stranded/]gene_tpm_matrix.csv
quantification/[stranded/]gene_expression_matrices.xlsx
quantification/stranded/strandness/{sample}.json
quantification/stranded/strandness/{sample}.env
logs/stranded/strandness/{sample}/strandness_diagnostic.tsv
logs/stranded/strandness/strandness_diagnostic.tsv  # 汇总，不是放行依赖
```

Excel 有 `alignment_report`、`gene_count`、`gene_fpkm`、`gene_tpm` 四个 sheet。它们的样本和基因顺序经过核对，ID 为文本、表达值为数值。日志在配置的 logs 目录；判链失败的诊断 TSV 和原始 RSeQC 输出仍保留。建立 DAG 前会检查所有声明的输出/日志与原始 FASTQ、参考、附加输入及脚本是否存在同路径或链接冲突，避免任务清理/日志打开时覆盖输入。Excel 导出要求整数回读精确一致，超过工作表行列上限时明确失败，CSV仍保留。

## 测试

```bash
cd /home/weiyuekuan/12.Snakemake/RNA-seq-pipeline
python -m unittest discover -s tests -v
```

缺少真实工具时部分测试会 skip；测试结果不能把 skip 当成通过。完整 FASTQ 集成显式设置 `RNASEQ_TEST_FULL_FASTQ=1`，可通过 `RNASEQ_TEST_SNAKEMAKE/HISAT2/HISAT2_BUILD/SAMTOOLS/STRINGTIE/RSEQC/GFFREAD` 指定已有工具；`RNASEQ_TEST_WORKFLOW_PYTHON` 指向包含 openpyxl 的 Python，运行测试的 Python 需含 pysam。fastp 实测使用 `RNASEQ_TEST_FASTP` 指定已有可执行文件，否则从 PATH 查找；raw FASTQ 经 fastp 到 Excel 的完整集成也需要上述完整 FASTQ 开关和工具。设置 `RNASEQ_TEST_KEEP=1` 保留人工输入与所有日志。

这些测试不下载依赖，不修改真实样本；包括 polyG 开关与实际保留读长、低质量/N/短 read 过滤、raw 配对错误阻断、raw 到 Excel、FR/RF与失败门控、含gene_name/竖线ID、每样本不同读长、Excel ID、GFF边界、路径引号、线程上限、Ballgown恢复和默认依赖重跑。

2026-10-23 fastp 接入的具体实测结果和日志见 [tests/FASTP_VALIDATION.md](tests/FASTP_VALIDATION.md)。
