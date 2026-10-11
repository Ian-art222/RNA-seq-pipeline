# 两分支升级验证记录

本次只使用人工注释、FASTQ 占位文件和合成 BAM；未运行真实样本，未修改已有分析结果，未下载或安装软件。

## 已通过的验证

- 41 项脚本测试：18 项注释/BED12 测试，23 项判链/质量检查测试。包含真实 gffread 转换、真实 RSeQC 采满目标数量时的输出解析，以及 GFF3 子记录在转换前的坐标/链一致性检查。
- 6 项工作流集成测试（10 次隔离的 Snakemake 调用）：普通与链特异性分支，非法/冲突配置，FR 与 RF，方向模糊与样本间方向冲突。
- 成功场景实际运行 StringTie 2.2.3、prepDE、矩阵汇总与 Excel 导出；核对样本列、基因集合、正数表达值、TPM 总和及 Excel sheet。
- 失败场景检查 Snakemake 非零退出、无放行文件、所有样本均未生成 StringTie GTF，且诊断 TSV 和原始 RSeQC 输出仍存在。

| 人工输入 | 实际采样量 / 目标 | 判链结果 | 后续行为 |
|---|---:|---|---|
| FR | 1000 / 1000 | `--fr` | 完成矩阵和 Excel |
| RF | 2400 / 200000 | `--rf`，采样不足警告 | 完成矩阵和 Excel |
| 非链特异性，普通模式 | 不执行 RSeQC | 无链方向参数 | 完成矩阵和 Excel |
| 两种方向约 1:1，链特异性模式 | 2400 / 200000 | FAIL | 阻断全部定量 |
| 同批一个 FR、一个 RF | 每个 2400 / 200000 | 各自可判定，整批方向冲突 | 阻断全部定量 |

脚本测试还包含缺失/重复比例、NaN/Infinity、非法范围、零采样/缺失计数、非双端结果、阈值边界、工具异常退出、BAM 与注释不兼容、旧放行文件失效等情形。

最后另作样本扫描检查：空 FASTQ 被拒绝，`strict=false` 仍能跳过不完整的非样本目录。

## 验证环境和边界

- Snakemake 9.19.0；RSeQC 5.0.4；pysam 0.24.0；gffread 0.12.9。
- 定量与导出使用已有 workflow 缓存环境：Python 3.11.9、StringTie 2.2.3、pandas 2.2.2、openpyxl 3.1.5。
- 判链测试使用本机 RNA-seq 环境的 Python 3.13.14；新增 Conda 定义使用 Python 3.11。本次没有创建/求解新的 Conda 环境，依赖安装不是此次验证范围。
- 集成测试从预构建的人工 BAM 开始，HISAT2 比对没有重新执行；这次改动的两分支及下游路径已实际执行。
- gffread 对仅 CDS 的测试注释未输出 exon，因此流程正确拒绝该注释；没有以 transcript 跨度静默补 exon。

集成测试使用临时目录中的 `rnaseq-strand-integration-*` 工作目录；设置 `RNASEQ_TEST_KEEP=1` 可保留并查看测试打印的实际路径。每个场景含 `snakemake.run.log` 或 `snakemake.dry-run.log`；链特异性场景的详细报告在 `logs/stranded/strandness/`。

## 重复运行

脚本测试（真实工具不在 PATH 时，其对应可选测试会跳过）：

```bash
python -m unittest discover -s tests -p 'test_strand*.py' -v
```

工作流集成测试用含 pysam 的 Python 执行，可通过环境变量选择现有工具。`RNASEQ_TEST_WORKFLOW_PYTHON` 需包含 pandas/openpyxl；测试会用临时命令分派复现两层业务环境，不安装依赖。

```bash
RNASEQ_TEST_SNAKEMAKE=/path/to/snakemake \
RNASEQ_TEST_RSEQC=/path/to/infer_experiment.py \
RNASEQ_TEST_STRINGTIE=/path/to/stringtie \
RNASEQ_TEST_WORKFLOW_PYTHON=/path/to/quantification/python \
RNASEQ_TEST_KEEP=1 \
/path/to/pysam/python -m unittest discover -s tests -p test_workflow_integration.py -v
```
