# RNA-seq Snakemake Pipeline
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
6. `StringTie -e -A` 定量，`-B` 可开关
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
│   │   └── prepDE.py
│   └── envs/
│       └── rnaseq.yaml
├── config/
│   └── config.yaml
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
必须同时存在一对 FASTQ 才识别为有效样本
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

1. dry-run（会解析 workflow；首次可能触发创建 `envs/rnaseq.yaml` 对应环境）：

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
