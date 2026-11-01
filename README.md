# bioinfo-classify

胃癌进展风险分型与多组学预测（TCGA-STAD）。

## 目录

```
.
├── exp.ipynb                        # 实验核心流程（A–F 六阶段，一次运行产出全部结果）
├── 胃癌进展风险分型报告.md            # 完整分析报告
├── extract.sh                       # 递归解压 生物信息学实验数据.rar → dataset/
├── 生物信息学实验数据.rar            # 原始数据归档（约 376 MB）
├── stad_tcga_pub_clinical_data.tsv   # TCGA 2014 官方分子亚型标签
├── dataset/                         # 解压后的数据文件
├── cache/                           # 预过滤中间结果（自动生成）
├── figures/                         # 全部图表 PNG + PDF（自动生成）
├── requirements.txt                 # 依赖版本锁定
└── pyproject.toml / .python-version # 项目元信息与 Python 版本
```

## 启动

### 1. 解压数据

```bash
bash extract.sh
```

脚本会自动探测 `unrar` / `7z` / `7zz` / `unar` / `bsdtar`，递归解开 `生物信息学实验数据.rar` 及其内部嵌套的压缩包到 `dataset/`。

缺少解压工具时：

```bash
brew install unrar            # macOS
sudo apt install unrar-free   # Debian / Ubuntu
```

### 2. 配置环境

需要 Python 3.14（见 `.python-version`）。

```bash
uv venv && uv pip install -r requirements.txt
# 或
python -m venv .venv && .venv/bin/pip install -r requirements.txt
```

### 3. 运行

```bash
jupyter lab exp.ipynb
```

在仓库根目录下按顺序从头运行，一次产出全部结果：中间结果缓存在 `cache/`，图片写入 `figures/`。

脱离 notebook 重跑：

```bash
jupyter nbconvert --to notebook --execute exp.ipynb --output exp.executed.ipynb
```

> 阶段 E 的通路富集调用 Enrichr 在线服务（`gseapy.enrichr`），断网时该节会失败，其余阶段不受影响。
