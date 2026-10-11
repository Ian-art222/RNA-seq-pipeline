# RNA-seq Snakemake Pipeline

## 链特异性定量开关

输入仍为同一参考下的双端 clean FASTQ。运行前设置配置文件中的 `stranded`，或在原运行命令末尾添加 `--config stranded=true` / `--config stranded=false`。省略该字段时默认 `false`，字符串 `false` 也会正确解析；其他取值报错。

- `stranded: false`：HISAT2 → 排序/索引 → StringTie 普通定量 → 原有矩阵与 Excel。
- `stranded: true`：HISAT2 → 排序/索引 → RSeQC 逐样本判链 → 全样本质量检查 → StringTie `--fr` / `--rf` 定量 → 矩阵与 Excel。

两条路共享不预设链方向的 `HISAT2 --dta` 比对。开启判链时，`alignment.hisat2_extra` 不得包含 `--rna-strandness`；已有 BAM 的 HISAT2 命令记录也会检查。StringTie 的 `--fr` / `--rf` 统一由开关与判链结果管理，不要写入 `stringtie.extra`。每个样本独立判链，所有样本通过且方向一致后才允许任何样本进入链特异性定量。混合建库协议请按协议拆成独立任务。

```yaml
stranded: true
strandness:
  sample_size: 200000
  mapq: 30
  min_determinable: 0.50
  min_dominance: 0.80
  bed12: ""  # 空值：从 reference.annotation 自动生成
```

`mode1` 为 `1++,1--,2+-,2-+` 的比例，`mode2` 为 `1+-,1-+,2++,2--` 的比例。保留原判定规则：`mode1 + mode2 >= 0.50` 且 `max(mode1, mode2)/(mode1 + mode2) >= 0.80`；mode1 占优用 `--fr`，mode2 占优用 `--rf`。这是本流程采用的阈值，不是 RSeQC 的通用质量标准。

`-s 200000` 是目标采样量。报告会记录 RSeQC stderr 中实际 usable reads 数；实际数量为零、数量缺失或超出目标时判链失败，正数但不足目标时记录警告，不另加未约定的读数阈值。接近判定边界或可用 reads 很少时，应结合原始报告复核。RSeQC 5.0.4 按 BAM 顺序读取符合条件的记录，不是全基因组随机抽样；两种模式之和也不是全部 reads 的注释率。

### 注释与输入检查

开启链特异性模式后，流程从同一 `reference.annotation` 生成 BED12：GTF 逐转录本使用 exon 坐标；GFF/GFF3 先经 gffread 转为 GTF，保留转换文件和日志。缺少 exon 的模型、无效坐标或不一致的转录本信息会明确报错，避免静默丢失模型。如果设置 `strandness.bed12`，会校验 BED12 并与当前注释的转录本模型核对。注释和 BAM 必须来自同一组装；序列名称、坐标检查通过不能单独证明组装版本相同。

判链前检查 BAM/BAI 可读性、排序信息和 BED12 对 BAM 序列名称/长度的兼容性。RSeQC 必须正常退出、识别为 PairEnd，并完整输出两种比例、无法判定比例和实际采样数；缺失/重复字段、非有限数值或非法范围均失败。`NA`、比例不达标或样本间方向冲突会使质量检查规则以非零状态退出，停止后续 StringTie，不自动改用普通定量。

### 输出位置与重跑

普通模式继续输出到配置的 `paths.quant_dir` 和 `output.xlsx`。链特异性模式自动在两者所属位置增加 `stranded/` 子目录，避免覆盖普通定量结果。默认路径为：

```text
quantification/stranded/{sample}/{sample}.gtf
quantification/stranded/{sample}/{sample}_abundance.tab
quantification/stranded/gene_count_matrix.csv
quantification/stranded/gene_fpkm_matrix.csv
quantification/stranded/gene_tpm_matrix.csv
quantification/stranded/gene_expression_matrices.xlsx
quantification/stranded/strandness/reference.bed12
quantification/stranded/strandness/{sample}.json
quantification/stranded/strandness/strandness.env
logs/stranded/strandness/prepare_bed.json
logs/stranded/strandness/{sample}/
logs/stranded/strandness/strandness_diagnostic.tsv
logs/stranded/strandness/gate.log
```

每样本 JSON 记录状态、原始比例、实际/目标采样量、方向、警告和失败原因；日志目录保留命令、版本、RSeQC stdout/stderr。质量检查的 TSV 作为日志保存，失败后仍保留；只有通过才生成 `strandness.env`，StringTie 显式依赖该文件。单样本判链步骤输出诊断记录，即使判链失败也先收集其它样本的结果，最终由质量检查规则阻断定量。

链特异性步骤使用独立的 `workflow/envs/strandness.yaml`（RSeQC 5.0.4、pysam、gffread）；普通模式只需要原来的 `rnaseq.yaml`。继续使用 `--use-conda` 自动管理环境。

修复注释、BAM 或配置后，Snakemake 会按依赖重算。若只修复了运行环境、工具故障，之前的失败 JSON 仍是已完成的诊断文件，需要在原命令中添加 `--forcerun infer_strandness strandness_gate` 重新判链。若旧 BAM 曾以错误或未知的链参数比对，应先核实来源，再按需要显式重跑比对；本升级不会自动重跑用户现有样本。

示例（先检查再正式运行）：

```bash
snakemake \
  --snakefile /home/pipeline/1.Snakemake-pipeline/1.RNA-seq-pipeline/workflow/Snakefile \
  --configfile /home/pipeline/1.Snakemake-pipeline/1.RNA-seq-pipeline/config/config.yaml \
  --directory /home/pipeline/2.Snakemake-workdir/1.RNA-seq-workdir \
  --use-conda --cores 32 --dry-run --printshellcmds \
  --config stranded=true
```

检查通过后去掉 `--dry-run` 执行。源码测试使用 `python -m unittest discover -s tests -v`；真实工具集成测试的环境配置见 `tests/test_workflow_integration.py`，本次验证范围与结果见 [tests/VALIDATION.md](tests/VALIDATION.md)。

## 1. 项目简介
这是一个基于 Snakemake 的 RNA-seq 表达定量流程，工具固定为：
- HISAT2
- Samtools
- StringTie
- Python
流程固定为：
1. 自动扫描 `cleandata/`
2. HISAT2 比对
3. `samtools view` 转 BAM，并支持可调 `mapq_filter`
4. `samtools sort`
5. `samtools index`
6. 按 `stranded` 开关选择普通定量或 RSeQC 判链后定量；StringTie 保持 `-e -A`，`-B` 可开关
7. 汇总比对 summary，生成 `alignment_report.csv`
8. `prepDE.py` 生成 raw count matrix
9. 汇总生成：
   - `gene_count_matrix.csv`
   - `gene_fpkm_matrix.csv`
   - `gene_tpm_matrix.csv`
10. 导出 `gene_expression_matrices.xlsx`
---
## 环境管理模式

本流程采用 **两层环境**，避免把 Snakemake 与 RNA-seq 业务软件混在同一 conda 环境中。

1. **Snakemake 专用环境**（本机示例名 `Snakemake`，亦可自建如 `snakemake_env`）：仅用于安装、激活后运行 Snakemake（负责解析 `Snakefile`、调度 rule、拉取 conda 环境）。**不要**在该环境中手动安装 HISAT2、Samtools、StringTie、pandas、openpyxl；这些由 workflow 环境提供即可。
2. **RNA-seq workflow Conda 环境**：由 `envs/rnaseq.yaml` 定义。使用 `snakemake --use-conda` 时，Snakemake 会根据该文件 **自动创建或复用** 独立环境，并在已声明 `conda:` 的 rule 内激活后执行命令。

**第一次运行**可能较慢，因为 Snakemake 需要先创建 `envs/rnaseq.yaml` 对应的 conda 环境；之后同一机器上通常会直接复用缓存环境。

**后续其他 pipeline**：建议各自维护独立的 `envs/*.yaml`，同样通过 `--use-conda` 调用，保持调度环境与业务依赖分离。

---
## 目录布局

本仓库将 **pipeline 源码** 与 **运行数据/结果** 彻底分离。

**Pipeline 源码目录**（仅包含 `Snakefile`、脚本、conda 环境定义与配置文件；不要在此目录下存放 `cleandata`、`alignment`、`quantification`、`.snakemake` 等运行产物）：

`/home/pipeline/1.Snakemake-pipeline/1.RNA-seq-pipeline/`

推荐布局：

```
/home/pipeline/1.Snakemake-pipeline/1.RNA-seq-pipeline/
├── workflow/
│   ├── Snakefile
│   ├── scripts/
│   │   ├── build_alignment_report.py
│   │   ├── build_expression_matrices.py
│   │   ├── export_expression_xlsx.py
│   │   ├── prepare_strand_bed.py
│   │   ├── strandness.py
│   │   └── prepDE.py
│   └── envs/
│       ├── rnaseq.yaml
│       └── strandness.yaml
├── config/
│   └── config.yaml
├── tests/
└── README.md
```

**Workdir 运行目录**（所有输入、输出、日志与 Snakemake 元数据；使用 `--directory` 指向此处）：

`/home/pipeline/2.Snakemake-workdir/1.RNA-seq-workdir/`

```
/home/pipeline/2.Snakemake-workdir/1.RNA-seq-workdir/
├── cleandata/
├── alignment/
├── quantification/
├── logs/
└── .snakemake/          # 首次运行后由 Snakemake 自动生成，无需手工创建
```

运行 Snakemake 时 **必须** 指定 `--snakefile` 与 `--directory`（workdir），以便相对路径 `cleandata`、`alignment`、`quantification` 等均解析到 workdir 下，而非源码目录。默认配置文件为源码目录下的 `config/config.yaml`（由 `Snakefile` 固定加载）；若使用副本或自定义路径，可再加 `--configfile /path/to/config.yaml` 覆盖。

本流程面向 Linux/bash 环境运行；Windows 用户请通过 WSL、Linux 服务器或集群环境执行。

---
## 2. 输入数据组织方式
以下路径均相对于 **workdir**（`--directory` 所设目录）。输入目录结构固定如下：
```bash
cleandata/{sample}/{sample}_1.fastq.gz
cleandata/{sample}/{sample}_2.fastq.gz
例如（本账号 workdir 内已有数据）：
cleandata/R_MBS1/R_MBS1_1.fastq.gz
cleandata/R_MBS1/R_MBS1_2.fastq.gz
cleandata/R_MBS2/R_MBS2_1.fastq.gz
cleandata/R_MBS2/R_MBS2_2.fastq.gz
若为其他命名，请在 config.yaml 的 sample_scan.r1_suffix / r2_suffix 中修改（如 _1.clean.fastq.gz）。
说明：
只扫描 cleandata/ 下一级子目录
子目录名就是 sample 名
样本名使用字母、数字、下划线、点或连字符，首字符为字母、数字或下划线
必须同时存在一对非空 FASTQ 文件才识别为有效样本
普通文件如 nohup.out 会被忽略
strict=true 时若发现不完整样本目录，流程会直接报错
strict=false 时会 warning 并跳过
```

## 3. 与目录布局的关系
- 源码与配置位于 **pipeline 源码目录**（见上文「目录布局」）。
- `cleandata/`、`alignment/`、`quantification/`、`logs/`、`.snakemake/` 仅存在于 **workdir**，不会出现在源码目录中。

## 4. config.yaml 参数说明
paths
clean_dir：输入 clean FASTQ 根目录
align_dir：比对结果输出目录
quant_dir：定量结果输出目录
logs_dir：运行日志输出目录，默认 logs
reference
hisat2_index_prefix：HISAT2 索引前缀
annotation：参考注释文件，可为 GTF 或 GFF
sample_scan
r1_suffix：R1 文件后缀（相对子目录名 sample），本账号默认 _1.fastq.gz
r2_suffix：R2 文件后缀，本账号默认 _2.fastq.gz
strict：扫描是否严格检查不完整样本目录
alignment
mapq_filter：samtools view 的 MAPQ 过滤阈值
hisat2_extra：传递给 HISAT2 的额外参数
stringtie
use_ballgown：是否加 -B
extra：传递给 StringTie 的额外参数
threads
hisat2
samtools_view
samtools_sort
samtools_index
stringtie
output
xlsx：最终 Excel 输出路径

## 5. 准备 HISAT2 index
先准备参考基因组 FASTA，然后建立索引，例如：
hisat2-build reference.fa reference.fa
之后在 config.yaml 中填写：
reference:
  hisat2_index_prefix: "/your/path/reference.fa"
注意：这里填写的是 索引前缀，不是随便一个目录说明文字。

本账号 `pipeline` 已在 `/home/pipeline/3.Ct_ZJ_genome/hisat2_index/` 下建好 HISAT2 索引，前缀为 `Ct_ZJ_T2T_index`（完整路径见 `config/config.yaml`）。

`workflow/envs/rnaseq.yaml` 中的 **HISAT2 版本应与建索引时所用版本一致**（当前约定为 **2.2.2**），否则运行时可能提示索引与程序版本不匹配。

## 6. 准备注释文件
将 GTF 或 GFF/GFF3 注释文件放到服务器某个路径，然后在 config.yaml 中填写：
reference:
  annotation: "/your/path/reference_annotation.gtf"
本账号（`pipeline`）默认参考与注释位于目录 `/home/pipeline/3.Ct_ZJ_genome/`（HISAT2 索引在 `hisat2_index/`，基因为 `Ct_ZJ_T2T_gene.gff3`），已在 `config/config.yaml` 中配置。

## 7. prepDE.py
Python3 兼容版 `prepDE.py` 已内置在源码目录 `workflow/scripts/prepDE.py`，无需再手动放到 workdir 或 `quantification/` 目录。

## 8. 准备 Snakemake 专用环境
请单独创建一个只用于运行 Snakemake 的 conda 环境，例如：
```bash
conda create -n Snakemake -c conda-forge -c bioconda snakemake
```
（名称可自定；本账号使用环境名 **`Snakemake`**。）
激活后 **无需** 再向该环境安装 HISAT2、Samtools、Stringtie、pandas、openpyxl；这些由 `envs/rnaseq.yaml` 在 `--use-conda` 时自动管理。

### `conda info --envs` 里出现一长串路径？
执行 `conda info --envs` 时，可能会多出一行 **无环境名**、路径形如  
`.../workdir/.snakemake/conda/<哈希>_`。那是 Snakemake 为 `rnaseq.yaml` 创建的 **workflow 业务环境**（HISAT2 / samtools / stringtie 等），由 Snakemake 在跑相应 rule 时自动激活；**不需要**也不建议手动 `conda activate` 到该路径。日常只需激活上方的 **Snakemake**（调度）环境即可。

## 9. 如何运行（两层环境 + `--use-conda` + 源码/workdir 分离）
在已激活 Snakemake 专用环境的前提下，**从任意目录**执行下列命令均可；须指定 `--snakefile` 与 `--directory`；默认使用源码目录中的 `config/config.yaml`（若需其它配置文件，可加 `--configfile`）。

**提醒**：

- 务必先 `conda activate Snakemake`（或你本机为 Snakemake 单独创建的其它环境名）。
- 首次运行会在 **workdir** 下生成 `.snakemake/`，并可能耗时创建 `rnaseq` conda 环境缓存。
- `cleandata`、比对与定量结果、Excel 等均在 workdir 下，**不会**写入源码目录。

1. dry-run（检查规则和依赖，不执行分析；环境安装使用下一步命令）：

```bash
snakemake \
  --snakefile /home/pipeline/1.Snakemake-pipeline/1.RNA-seq-pipeline/workflow/Snakefile \
  --directory /home/pipeline/2.Snakemake-workdir/1.RNA-seq-workdir \
  --use-conda \
  -n
```

2. 只创建 conda 环境并下载依赖（可选，适合第一次运行前先准备环境）：

```bash
snakemake \
  --snakefile /home/pipeline/1.Snakemake-pipeline/1.RNA-seq-pipeline/workflow/Snakefile \
  --directory /home/pipeline/2.Snakemake-workdir/1.RNA-seq-workdir \
  --use-conda \
  --conda-create-envs-only \
  --cores 32
```

3. 正式运行：

```bash
snakemake \
  --snakefile /home/pipeline/1.Snakemake-pipeline/1.RNA-seq-pipeline/workflow/Snakefile \
  --directory /home/pipeline/2.Snakemake-workdir/1.RNA-seq-workdir \
  --use-conda \
  --cores 32 -p
```

4. 若中断后重跑：

```bash
snakemake \
  --snakefile /home/pipeline/1.Snakemake-pipeline/1.RNA-seq-pipeline/workflow/Snakefile \
  --directory /home/pipeline/2.Snakemake-workdir/1.RNA-seq-workdir \
  --use-conda \
  --cores 32 --rerun-incomplete -p
```

**说明**：务必加上 `--use-conda`，否则 Snakemake 不会按 `workflow/envs/rnaseq.yaml` 为各 rule 准备业务软件环境。

运行过程中，Snakemake 自身调度日志会保存在 `.snakemake/log/`；各步骤的详细日志会保存在 `logs/`：

```bash
logs/{sample}/hisat2_align_sort.log
logs/{sample}/samtools_index.log
logs/{sample}/stringtie_quant.log
logs/build_sample_list.log
logs/build_alignment_report.log
logs/prepde_raw_counts.log
logs/finalize_expression_matrices.log
logs/export_expression_xlsx.log
```

## 10. mapq_filter=30 和 mapq_filter=0 的区别
mapq_filter = 30
实际执行：
samtools view -b -q 30
含义：
只保留 MAPQ >= 30 的比对结果
更严格
可能会降低部分基因的计数
对多重比对 reads 更保守
mapq_filter = 0
实际执行：
samtools view -b
含义：
不做 MAPQ 过滤
更接近普通 RNA-seq 默认流程
保留更多 reads

## 11. 最终输出文件说明
alignment
alignment/{sample}/{sample}.sorted.bam
alignment/{sample}/{sample}.sorted.bam.bai
alignment/{sample}/{sample}.align_summary.txt
alignment/alignment_report.csv
logs
logs/{sample}/hisat2_align_sort.log
logs/{sample}/samtools_index.log
logs/{sample}/stringtie_quant.log
logs/build_alignment_report.log
logs/prepde_raw_counts.log
logs/finalize_expression_matrices.log
logs/export_expression_xlsx.log
quantification
quantification/{sample}/{sample}.gtf
quantification/{sample}/{sample}_abundance.tab
quantification/sample_lst.txt
quantification/gene_count_matrix.raw.csv
quantification/transcript_count_matrix.csv
quantification/gene_count_matrix.csv
quantification/gene_fpkm_matrix.csv
quantification/gene_tpm_matrix.csv
quantification/gene_expression_matrices.xlsx

## 12. Excel sheet 说明
最终 Excel 文件 quantification/gene_expression_matrices.xlsx 包含以下 sheet：
alignment_report
gene_count
gene_fpkm
gene_tpm
