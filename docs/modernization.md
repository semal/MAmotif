# MAmotif 现代化说明

## 变更摘要

- 全库从 **Python 2.7** 迁移到 **Python 3.8+**
- `print` / `xrange` / `except E, e` / `cPickle` / `file()` / `raw_input` / `time.clock` 等语法与 API 已替换
- 包内改为相对导入（`MAmotif_pkg`）
- 主入口与 `command_tools/` 改为 `argparse`
- 新增 `pyproject.toml`、`requirements.txt`、可编辑安装与控制台命令 `mamotif`
- 新增 `tests/test_smoke.py` 冒烟测试
- `correct_pvalues` 兼容含 `None` 的列表与现代 NumPy
- 示例路径常量改为懒加载（`constant.py`）

## 安装与自检

```bash
pip install -e ".[dev]"
pytest
mamotif -h
```

## MotifScan pickle 说明

MotifScan 的 pandas pickle 若由 **Python 2** 写出，在 Python 3 下可能无法直接 `read_pickle`。建议：

1. 用当前 Python 3 环境重新跑 MotifScan 生成结果；或
2. 在仍可用的 Py2 环境中读入后以协议 2+ / 兼容格式导出，再在 Py3 中加载。

MATLAB `.mat` MotifScan 结果不受此限制。

## 未改动的部分

- 生物学算法与输出列含义保持原设计
- `analyses/`、`RNAseq_analyses/` 中部分脚本仍含本机历史绝对路径（仅用于作者本地实验复现）
- MotifScan / MAnorm / DESeq 等上游工具仍需自行安装
