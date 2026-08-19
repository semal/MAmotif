# MAmotif

Integrating motif information with ChIP-seq binding signal change.

MAmotif 将 [MAnorm](https://github.com/shao-lab/MAnorm) 定量比较得到的结合信号变化（M 值）与 motif 扫描结果整合在一起，用来识别在两种细胞类型或实验条件下，与差异结合相关的转录因子 motif。

核心思路：把一组 peaks 按是否含有某个 motif 分成两类，比较两类 peak 的 M 值分布（t 检验 / Wilcoxon rank-sum），从而判断该 motif 是否与结合强度的细胞特异性变化相关。

## 功能概览

- 将 MAnorm 比较 peaks 划分为 promoter peaks 与 distal peaks
- 按 motif 存在与否对 peaks 分组，并对 M 值做分布差异检验
- 对 p 值做 Benjamini–Hochberg 或 Bonferroni 多重检验校正
- 将 motif 靶基因与差异表达基因做 Fisher 精确检验
- 关联 peak M 值与靶基因的表达 fold change（Pearson 相关）
- 辅助的 peak 集合重叠、基因集合富集、RNA-seq 预处理工具

支持人类 hg19 与小鼠 mm9 参考注释。内置 JASPAR3 中 196 个 motif 名称列表。

## 依赖环境

本仓库代码使用 **Python 2** 语法编写（例如 `print` 语句、`xrange`、`except Exception, e`），建议在 Python 2.7 下运行。

Python 依赖：

- numpy
- scipy
- pandas
- matplotlib（仅部分可视化脚本需要）

RNA-seq 差异表达分析额外需要 R 与 DESeq。

Motif 扫描结果来自 MotifScan（支持 pandas pickle 结果，以及 MATLAB `.mat` 格式）。

将仓库根目录加入 `PYTHONPATH` 后即可调用包内模块：

```bash
export PYTHONPATH=/path/to/MAmotif:$PYTHONPATH
```

## 分析流程

```text
MAnorm peak 文件  +  MotifScan 结果  +  RefGene 注释
                │
                ▼
     按 promoter 区域分类（可选）
     ├── 全部 peaks
     ├── promoter peaks（TSS ± 2 kb）
     └── distal peaks
                │
                ▼
     按每个 motif 的 target number 分组
     ├── 含该 motif 的 peaks
     └── 不含该 motif 的 peaks
                │
                ▼
     t-test / rank-sum test + 多重检验校正
                │
                ▼
     MAmotif 结果表（.xls）
```

可选后续分析：

1. 用 motif 靶基因与差异表达基因做 Fisher 精确检验
2. 将 peak M 值与靶基因表达变化做相关分析

## 快速开始

主入口是 `MAmotif.py`。给定一份 MAnorm 比较结果中的 peak 文件和对应的 MotifScan 结果，脚本会：

1. 若提供 RefGene 文件，将 peaks 分成 promoter / distal 两类
2. 对全部 peaks（以及分类后的两类 peaks）做 MAmotif 检验
3. 默认在与 MAnorm 比较同名的输出目录中写出结果

```bash
python MAmotif.py \
    -p path/to/sample_peak_MAvalues.xls \
    -M path/to/motifscan_peak_result \
    -r data/modified_hg19_refgene.txt \
    -c benjamin
```

常用参数：

| 参数 | 说明 |
| --- | --- |
| `-p` | MAnorm 比较得到的 peak 文件 |
| `-M` | MotifScan 结果（pandas pickle，或 MATLAB `.mat`） |
| `-r` | RefGene 注释文件；缺省则不做 promoter / distal 分类 |
| `-c` | p 值校正方式：`benjamin`（默认）、`bonferroni`，或其他值表示不校正 |
| `-n` | 对 M 值取负后再检验（比较方向对调） |
| `-d` | 不创建与 MAnorm 比较同名的输出目录，直接在当前目录写出 |

输出文件名形如：

```text
<peak文件名去掉 _MAvalues>_MAmotif_jaspar_output.xls
```

结果按校正后的最大 p 值排序，主要列包括：

- Motif 名称
- 含 / 不含该 motif 的 peak 数目、M 值均值与标准差
- t 检验统计量与右尾 p 值、校正后 p 值
- Rank-sum 检验统计量与右尾 p 值、校正后 p 值
- 两种校正 p 值中的较大者（Maximal P-value）

## 输入文件格式

### MAnorm peak 文件

制表符分隔，`#` 开头的行视为注释。示例见 `data/H1hesc_H3K27ac_Broad_Rep2_promoter_peak_MAvalues.xls`：

```text
#chr	start	end	summit	MAnorm_Mvalue	MAnorm_Avalue	pvalue
chr1	12345	13000	400	1.20	5.30	1.2e-05
```

列含义：染色体、起止坐标、summit（相对 start 的偏移或绝对坐标，由文件决定）、M 值、A 值、p 值。summit / A / p 值缺失时仍可读取 chr、start、end。

### MotifScan 结果

- pandas pickle：每行一个 peak（`chr` / `start` / `end` / `summit`），每个 motif 对应一列 `<motif>.tarnum`
- MATLAB `.mat`：Shao lab MotifScan 输出，含 `chr_id`、`bpstart`、`bpend`、`motif_name`、`motif_tarnum`

### RefGene 注释

UCSC RefSeq 基因表格式。仓库已提供：

- `data/modified_hg19_refgene.txt`
- `data/modified_mm9_refgene.txt`

Promoter 定义为 TSS 上游 2 kb 到下游 2 kb。与该区间重叠的 peak 记为 promoter peak，其余为 distal peak。

### 基因表达变化文件

用于将 M 值与表达变化关联，示例见 `data/H1hesc_K562_genes_fold_change.txt`：

```text
	H1hesc_mean	K562_mean	fold_change	log2_fold_change
DEFA5	0.020756	0.009441	2.198577	1.136570
```

程序读取第一列基因名与第五列 `log2_fold_change`。

## 命令行工具

`command_tools/` 提供可单独调用的分析步骤。

| 脚本 | 作用 |
| --- | --- |
| `classify_pks_by_promoter.py` | 按是否与启动子重叠，将 peaks 分成 promoter / distal |
| `classify_pks_by_motif.py` | 按指定 motif 是否存在，将 peaks 分成两类 |
| `classify_pks_by_pvalue.py` | 按 MAnorm p 值将 peaks 分成显著变化 / 其余 |
| `motifs_classified_peaks_ttest_result.py` | 对全部 motif 做 M 值 t 检验 / rank-sum 检验 |
| `other_peaks_classified_peaks_ttest_result.py` | 用另一组 ChIP-seq peaks 的重叠关系分类并检验 |
| `get_target_genes.py` | 输出 peaks 的靶基因列表 |
| `motif_targenes_de_fisher_test.py` | motif 靶基因是否富集于差异表达基因 |
| `motifs_classified_peaks_target_genes_ttest.py` | 按 motif 分组后，比较两组靶基因的表达变化 |
| `assoicate_Mvalue_GE.py` | 计算 peak M 值与靶基因表达变化的 Pearson 相关 |

示例：

```bash
# 按启动子分类
python command_tools/classify_pks_by_promoter.py \
    -p path/to/peaks_MAvalues.xls \
    -r data/modified_hg19_refgene.txt

# 按单个 motif 分类
python command_tools/classify_pks_by_motif.py \
    -p path/to/peaks_MAvalues.xls \
    -M path/to/motifscan_result \
    -n SOX2

# motif 靶基因 vs 差异表达基因的 Fisher 检验
python command_tools/motif_targenes_de_fisher_test.py \
    -p path/to/peaks_MAvalues.xls \
    -M path/to/motifscan_result \
    -r data/modified_hg19_refgene.txt \
    -d data/K562_H1hesc_DEseq_DE_FC_2_p1e-07_genes.txt

# M 值与基因表达变化的相关
python command_tools/assoicate_Mvalue_GE.py \
    -p path/to/peaks_MAvalues.xls \
    -r data/modified_hg19_refgene.txt \
    -f data/H1hesc_K562_genes_fold_change.txt
```

## 目录结构

```text
MAmotif/
├── MAmotif.py                 # 主分析流程
├── MAmotif_test               # 主流程的早期版本
├── MAmotif_fisher_test.py     # motif 靶基因 / peak 重叠的 Fisher 检验
├── MAmotif_pkg/               # 核心数据结构与 I/O
│   ├── MAmotifIO.py           # 读取 MAnorm / MotifScan / RefGene，分类与检验
│   ├── MAmotifPeaks.py        # MAnormPeak、MotifScanPeak
│   ├── MAnormPeaksClassifier.py
│   ├── sequence.py            # Sequence / Gene / Peak 等基础类
│   └── sequence_permutation.py
├── command_tools/             # 独立命令行工具
├── analyses/                  # 下游分析与可视化脚本
├── RNAseq_analyses/           # RNA-seq 计数、RPKM、DESeq 相关脚本
├── peak_set_analysis.py       # peak 集合重叠与富集
├── gene_set_analysis.py       # 基因集合重叠与 Fisher 检验
├── jaspar3.py                 # JASPAR3 motif 名称与扫描结果读取
├── constant.py                # 示例数据路径
└── data/                      # 参考注释与示例数据
    ├── modified_hg19_refgene.txt
    ├── modified_mm9_refgene.txt
    ├── jaspar3_motif_name_list
    ├── H1hesc_H3K27ac_Broad_Rep2_promoter_peak_MAvalues.xls
    ├── H1hesc_K562_genes_fold_change.txt
    └── K562_H1hesc_DEseq_DE_FC_2_p1e-07_genes.txt
```

## 应用场景

- 比较两种细胞类型中组蛋白修饰或转录因子结合的差异，找出与差异相关的 motif
- 区分 promoter 与 distal 区域上的 motif 效应
- 将 ChIP-seq 结合变化与 RNA-seq 差异表达联系起来

## 相关说明

- 上游定量比较依赖 MAnorm；motif 扫描依赖 MotifScan。本仓库完成匹配、分类与统计检验。
- `data/` 中的示例以 H1 hESC vs K562 的 H3K27ac 比较为主。
- 更完整的模块说明见仓库中的 `项目总结报告.md`。
