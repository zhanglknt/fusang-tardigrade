# Fusang: Tardigrade Edition — Reproducibility Package v1.4

This package reproduces all results from:

> **"Spaced K-mer Frequency Vector Cosine Distances Enable Robust Phylogenetic Inference Under Indel-Rich Conditions"**
> Submitted to *Nucleic Acids Research*, 2026.

**v1.4 Update (2026-06-28)**: Updated Mash benchmark with full 30-seed results (was single-seed). Added Jupyter notebook 07 for Mash vs Fusang comparison. Updated manuscript references to v2.2.

---

## 1. Working Environment (工作环境)

### Operating System

| OS | Supported | Notes |
|----|-----------|-------|
| **Linux** (x86_64) | ✅ 完全支持 | Ubuntu 20.04+ / CentOS 7+ 推荐 |
| **macOS** (Intel / Apple Silicon) | ✅ 完全支持 | macOS 12+ 推荐 |
| **Windows** | ✅ 支持（无 GPU 依赖） | Windows 10/11 x64，需 WSL 或 Git Bash 运行 shell 脚本 |

### Python Version

| Python | 状态 | 说明 |
|--------|------|------|
| **3.9 – 3.13** | ✅ 已测试 | 本包在 **Python 3.13.12** 上开发和验证 |
| 3.8 及以下 | ❌ 不支持 | 需要 `numba>=0.58`，不支持 Python 3.8 |
| 3.14+ | ⚠️ 未测试 | 理论可工作，但 `numba`/`scipy` 可能未适配 |

### Hardware Requirements (硬件要求)

| 场景 | RAM | CPU | 磁盘 | 可处理规模 |
|------|-----|-----|------|-----------|
| **最小 (L0-1 only)** | 4 GB | 2 核 | 2 GB | n ≤ 2,000 |
| **推荐 (全管线)** | 16 GB | 4 核 | 5 GB | n ≤ 10,000 |
| **完整 Benchmark** | 32 GB | 8 核 | 10 GB | n ≤ 10,000 (100+ seeds) |
| **IQ-TREE2 对比** | 16 GB | 4 核 | 5 GB | n ≤ 1,000 (可能超时) |

**注意**: 本包不需要 GPU / CUDA。所有计算在 CPU 上完成。

### Time Estimates (耗时预估)

| 任务 | n=200 | n=500 | n=1000 | n=2000 |
|------|-------|-------|--------|--------|
| Fusang L0-1 (1 seed) | ~5s | ~15s | ~40s | ~120s |
| Fusang MHL (1 seed) | ~8s | ~30s | ~90s | ~300s |
| Full benchmark (30 seeds) | ~3 min | ~8 min | ~25 min | ~60 min |
| IQ-TREE2 GTR (1 seed) | ~10 min | ~60 min | >3 h（超时） | >12 h（超时） |

---

## 2. Required Software (所需软件)

### 2.1 Python Packages (必装)

```bash
pip install -r requirements.txt
```

`requirements.txt` 内容：

```
numpy>=1.24          # 矩阵运算
scipy>=1.10          # 层次聚类、统计
scikit-learn>=1.3    # 标量归一化、随机森林分类器
biopython>=1.81      # Newick 解析、NJ 建树
matplotlib>=3.7      # 绘图
seaborn>=0.12        # 统计可视化
pandas>=2.0          # 数据处理
jupyter>=1.0         # Notebooks
numba>=0.58          # k-mer 计数加速
psutil>=5.9          # 内存/CPU 监控
```

**验证安装：**
```bash
python -c "import numpy, scipy, sklearn, Bio, numba, pandas; print('All OK')"
```

### 2.2 External Tools (外部工具 — 按需安装)

以下工具**不出现在 requirements.txt 中**，因为它们不是 Python 包。
请根据你要复现的结果，按需下载安装。

| # | 工具 | 版本 | 下载地址 | 用途 | 必须？ |
|---|------|------|----------|------|--------|
| **T1** | **MAFFT** | v7.505+ | https://mafft.cbrc.jp/alignment/software/ | 多序列比对 (L3 MSA) | 可选（仅 L3 benchmark） |
| **T2** | **FastTree2** | v2.1.11+ | http://www.microbesonline.org/fasttree/ | ML 建树 (L3 ML) | 可选（仅 L3 benchmark） |
| **T3** | **FastME** | v2.1.6+ | http://www.atgc-montpellier.fr/fastme/ | BME/NJ 建树 (L1 加速替代) | 可选 |
| **T4** | **IQ-TREE2** | v2.2.0+ | http://www.iqtree.org/ | 金标准 ML (P4 对比) | 可选（仅 IQ-TREE 对比） |
| **T5** | **INDELible** | v1.03+ | http://abacus.gene.ucl.ac.uk/software/indelible/ | Indel-rich 序列生成 | 可选（仅 repeat indel benchmark） |

### 2.3 Tool-to-Result Mapping (工具 → 复现结果映射)

下表说明复现论文中每个**结果/表格**需要哪些工具：

| 论文结果 | 脚本/Notebook | 需要的工具 |
|----------|---------------|------------|
| **Fusang L0-1 基准** (nRF, 时间) | `scripts/generate_missing_data.py` | **仅 Python** ✅ |
| **多 k-mer ensemble** (DAHP-V3) | `scripts/generate_multik.py` | **仅 Python** ✅ |
| **可扩展性** (n=200–2000, 内存/时间) | `scripts/scalability_nrf.py` | **仅 Python** ✅ |
| **Table 1: L0-1 vs MSA+ML** | `scripts/benchmark_mhl.py --methods nj,mhl,ft2` | Python + **MAFFT** + **FastTree2** |
| **SwissTree 蛋白质验证** | `notebooks/03_swisstree_validation.ipynb` | Python + **MAFFT** |
| **IQ-TREE2 金标准对比** | `scripts/benchmark_mhl.py --methods nj,iqtree` | Python + **IQ-TREE2** |
| **Mash vs Fusang 对比** (MinHash Jaccard vs k-mer cosine, 30 seeds) | `notebooks/07_mash_benchmark.ipynb` | **仅 Python** ✅ (预计算数据) |
| **Boundary 分类器训练** | `notebooks/04_boundary_classifier.ipynb` | **仅 Python** ✅ |
| **Indel-rich 基准** (n=200/500/1000) | 见 bench_tools/README.md | Python + **INDELible** |
| **16S rRNA 验证** | 外部数据 (SILVA/Greengenes) | Python + **MAFFT** + 数据库下载 |

**✅ = 开箱即用，无需外部工具。**

---

## 3. Environment Setup (环境配置流程)

### 3.1 快速开始（仅 Python，无需外部工具）

```bash
# 1. 确认 Python 版本
python --version   # 应输出 3.9 – 3.13

# 2. 创建虚拟环境（推荐）
python -m venv fusang_env
# Linux/macOS:
source fusang_env/bin/activate
# Windows:
# fusang_env\Scripts\activate

# 3. 安装依赖
cd repro_package
pip install -r requirements.txt

# 4. 验证
python -c "from fusang.fusang_v2 import main; print('Fusang OK')"
```

### 3.2 安装外部工具（按需）

#### Linux (Ubuntu/Debian)

```bash
# MAFFT
sudo apt-get install mafft

# FastTree2
sudo apt-get install fasttree

# IQ-TREE2
sudo apt-get install iqtree2

# FastME (手动下载)
wget http://www.atgc-montpellier.fr/download/sources/fastme/fastme-2.1.6.4-linux64.tar.gz
tar xzf fastme-2.1.6.4-linux64.tar.gz
mv fastme-2.1.6.4-linux64/binaries/fastme-2.1.6.4-linux64-64bit fastme
chmod +x fastme
sudo mv fastme /usr/local/bin/

# INDELible (手动下载)
wget http://abacus.gene.ucl.ac.uk/software/INDELibleV1.03.tar.gz
tar xzf INDELibleV1.03.tar.gz
cd INDELibleV1.03 && make
sudo cp indelible /usr/local/bin/
```

#### macOS

```bash
# MAFFT
brew install mafft

# FastTree2
brew install fasttree

# IQ-TREE2
brew install iqtree2

# FastME
brew install fastme
```

#### Windows

| 工具 | 安装方式 |
|------|----------|
| **MAFFT** | 下载 Windows 版: https://mafft.cbrc.jp/alignment/software/windows.html |
| **FastTree2** | 下载 `FastTree.exe`: http://www.microbesonline.org/fasttree/#Install |
| **IQ-TREE2** | 下载 Windows 版: http://www.iqtree.org/#download |
| **INDELible** | 需 WSL 或 Cygwin 编译（仅 Linux/macOS 二进制发布） |

安装到任意目录后，**设置环境变量**或**将工具放在 `PATH` 中**：

```bash
# 方法 1: 放到 PATH 中
# 将 exe 所在目录加入系统 PATH

# 方法 2: 设置环境变量（Python 脚本会自动读取）
export MAFFT="/c/Program Files/MAFFT/mafft.bat"
export FASTTREE="/c/Program Files/FastTree/FastTree.exe"
export FASTTREE_EXE="/c/Program Files/FastTree/FastTree.exe"
export IQTREE="/c/Program Files/IQ-TREE2/bin/iqtree2.exe"
```

### 3.3 验证安装

```bash
# 检查 Python 包
python -c "import numpy, scipy, sklearn, Bio, numba, pandas; print('Python packages: OK')"

# 检查外部工具
mafft --version      2>/dev/null    && echo "MAFFT: OK"      || echo "MAFFT: NOT FOUND"
FastTree             2>/dev/null    && echo "FastTree2: OK"  || echo "FastTree2: NOT FOUND"
fastme --version     2>/dev/null    && echo "FastME: OK"     || echo "FastME: NOT FOUND"
iqtree2 --version    2>/dev/null    && echo "IQ-TREE2: OK"   || echo "IQ-TREE2: NOT FOUND"
```

---

## 4. Package Structure (包结构)

```
repro_package/
|-- README.md                          # ← 本文件
|-- requirements.txt                   # Python 依赖
|
|-- data/                              # 所有基准数据 + 模型
|   |-- benchmark_clean_seeds130_229.csv     # JC69 n=200 clean (100 seeds)
|   |-- benchmark_mhl_n200_30seeds.csv       # MHL L0/L1/FT2 benchmark
|   |-- benchmark_multik_seeds200_229.csv    # Multi-k ensemble (30 seeds)
|   |-- indel_benchmark_seeds100_129.csv     # Indel n=200 (30 seeds)
|   |-- indel_benchmark_n500.csv             # Indel n=500
|   |-- indel_benchmark_n1000.csv            # Indel n=1000
|   |-- indel_benchmark_MASTER_REBUILT.csv   # 全部 seeds 合并
|   |-- scalability_results.json            # 可扩展性 (时间/内存/nRF)
|   |-- scalability_nrf.json                # nRF 详细
|   |-- benchmark_swisstree_results.csv      # SwissTree 蛋白质基准
|   |-- training_data.pkl                   # 分类器训练数据
|   |-- boundary_rf.pkl                     # 预训练 RF 分类器
|   |-- multik_ensemble_n200_30seeds.csv    # Multi-k 30-seed ensemble
|   |-- mash_benchmark_30seeds.csv          # Mash 30-seed benchmark CSV
|   |-- mash_vs_immi_results.json           # Fusang vs Mash benchmark (30 seeds each)
|   |-- e2e_extended_results.json           # Boundary classifier E2E (88 scenarios)
|   |-- l3_validation_n200/                 # L3 validation (L0/L1/FT2, 30 seeds)
|   |-- table8_results_TRUE_reference.csv   # Co-phylog vs Fusang (TRUE ref)
|   |-- DATA_SOURCES.json                   # v2.0 claim-to-data audit map
|
|-- fusang/                            # 核心源码
|   |-- __init__.py
|   |-- fusang_v2.py                       # 主管线 CLI
|   |-- fusang_v4_dahp_v1.py              # DAHP 增强管线
|   |-- fusang_mhl_main.py                # MHL 入口
|   |-- kmer_distance.py                  # k-mer 频率计算
|   |-- kmer_optimized.py                 # 优化 k-mer 例程
|   |-- tree_simulation.py                # 溯祖树 + JC69 模拟器
|   |-- fusang_minhash_wrapper.py         # MinHash 距离
|   |-- calc_nrf_simple.py                # nRF 计算
|   |-- af_competitor_methods.py          # 竞品方法实现
|   |-- ensemble.py                       # Multi-k ensemble
|   |-- fusang.py                         # 核心类型和工具
|   |-- fusang_mhl/                       # MHL 子包 (11 模块)
|       |-- config.py                     # 全局配置
|       |-- mlh_utils.py                  # 共享工具函数
|       |-- level0_kmer.py                # L0: k-mer 特征提取
|       |-- level1_multik.py              # L1: multi-k 距离矩阵
|       |-- level2_dahp.py                # L2: DAHP 分区
|       |-- level3_msa_ml.py              # L3: MSA + ML 建树
|       |-- merger.py                     # 子树合并
|       |-- boundary_classifier.py        # 边界分类器 (ML)
|       |-- training_data_generator.py    # 训练数据生成
|       |-- ml_split.py                   # ML split 决策 (V4b)
|
|-- scripts/                           # 独立脚本
|   |-- generate_missing_data.py           # 生成 clean benchmark 数据
|   |-- generate_multik.py                 # Multi-k ensemble benchmark
|   |-- scalability_nrf.py                 # 可扩展性 + nRF
|   |-- benchmark_mhl.py                   # MHL benchmark 运行器
|   |-- benchmark_competitors.py           # 竞品方法 benchmark
|   |-- benchmark_multik_ensemble.py       # Multi-k ensemble 完整 benchmark
|   |-- benchmark_swisstree.py             # SwissTree 验证
|   |-- train_boundary_classifier.py       # 训练边界分类器
|   |-- merge_benchmarks.py                # 合并 benchmark 结果
|   |-- summarize_results.py               # 结果汇总
|   |-- validate_l3_e2e.py                 # L3 (MSA+ML) 端到端验证
|   |-- run_l3_batches.py                  # L3 分批运行器
|   |-- run_e2e_structured.py              # E2E structured tree 测试
|   |-- compute_mash_vs_immi.py            # Mash vs IMMI benchmark
|   |-- compute_all_nrf.py                 # Mash nRF 计算
|   |-- test_mhl_e2e_v4b.py                # MHL V4b 端到端测试
|
|-- notebooks/                         # Jupyter Notebooks
|   |-- 01_main_benchmark.ipynb            # 主 Benchmark 分析
|   |-- 02_multik_ensemble.ipynb           # Multi-k ensemble 分析
|   |-- 03_swisstree_validation.ipynb      # SwissTree 蛋白质验证
|   |-- 04_boundary_classifier.ipynb       # 边界分类器训练/评估
|   |-- 05_scalability_demo.ipynb          # 可扩展性演示
|   |-- 06_generate_figures.ipynb          # 论文用图生成
|   |-- 07_mash_benchmark.ipynb             # Mash vs Fusang 对比 (30 seeds)
|
|-- bench_tools/                       # 外部工具说明
|   |-- README.md                          # 下载和安装指南
|
`-- figures/                           # 生成图表输出
    `-- .gitkeep
```

---

## 5. Reproducing Results (复现结果)

### 5.1 开箱即用（无需外部工具）

#### Clean Benchmark (n=200, 100 seeds)

```bash
cd repro_package/scripts
python generate_missing_data.py \
    --start-seed 130 --end-seed 229 --n 200 \
    --output ../data/benchmark_clean_seeds130_229.csv
```

**说明**: 使用 JC69 模拟器（`tree_simulation.py`）。这是无 indel 的 clean 数据。
论文中 indel-rich benchmark 使用 INDELible 生成，nRF 数值会有所不同，但相对排名一致。

#### Multi-k Ensemble (n=200, 30 seeds)

```bash
cd repro_package/scripts
python generate_multik.py \
    --start-seed 200 --end-seed 229 --n 200 \
    --output ../data/benchmark_multik_seeds200_229.csv
```

#### Scalability (n=200–2000)

```bash
cd repro_package/scripts
python scalability_nrf.py \
    --scales 200 500 1000 2000 --seeds-per-scale 5 \
    --output ../data/scalability_nrf.json
```

**预计耗时**: 5 seeds × 4 scales ≈ 15–20 分钟。

### 5.2 需要外部工具

#### Table 1: L0-1 vs MSA+ML (需要 MAFFT + FastTree2)

```bash
# 先确认工具可用
mafft --version && FastTree

# 运行 benchmark (30 seeds)
cd repro_package/scripts
python benchmark_mhl.py \
    --seeds 100-129 --no-l3 --methods nj,mhl,ft2 \
    --output ../data/benchmark_mhl_n200_30seeds.csv
```

**预计耗时**: ~15 分钟 (30 seeds, n=200)。

#### IQ-TREE2 Gold-Standard Comparison (需要 IQ-TREE2)

```bash
iqtree2 --version   # 确认安装

cd repro_package/scripts
python benchmark_mhl.py \
    --seeds 100-129 --no-l3 --methods nj,iqtree \
    --output ../data/benchmark_iqtree.csv
```

**警告**: IQ-TREE2 在 n≥500 时可能超时（>1 小时/seed）。建议从 n=200 开始。

#### Indel-rich Benchmark Generation (需要 INDELible)

Indel-rich 数据使用 INDELible 控制文件生成。
详见 `bench_tools/README.md`。

#### SwissTree Validation

1. 从 AFproject 下载 SwissTree 数据: https://afproject.org
2. 运行: `python scripts/benchmark_swisstree.py --datadir <swisstree_dir>`

---

## 6. Running Fusang (命令行使用)

```bash
# 从 FASTA 构建树 (L0-1, 直接 NJ)
python fusang/fusang_v2.py -i input.fasta -o output.nwk -m default \
    --auto-group --kmer-k 5 --kmer-gap gap2

# Multi-k ensemble 模式 (DAHP-V3)
python fusang/fusang_v4_dahp_v1.py -i input.fasta -o output.nwk --v3

# 完整 MHL 管线 (Level 0-1-2-3)
python fusang/fusang_mhl_main.py -i input.fasta -o output.nwk --all-levels
```

---

## 7. Important Notes (重要说明)

### 7.1 Simulator Discrepancy (模拟器差异)

| 数据 | 本包中的模拟器 | 论文中使用的工具 | nRF 是否可比？ |
|------|---------------|-----------------|---------------|
| Clean (no-indel) | `tree_simulation.py` (Kingman + JC69) | `tree_simulation.py` (Kingman + JC69) | ✅ 完全一致 |
| Indel-rich | **不可在本包中生成** | INDELible (birth-death + GTR+Γ + indels) | — |
| 预计算 data/*.csv | **来自 INDELible** | INDELible | ✅ 使用 data/ 中的预计算文件 |

**结论**: Clean benchmark 可直接复现。Indel benchmark 请使用 `data/` 目录下的预计算 CSV。

### 7.2 nRF Normalization (nRF 归一化)

本文使用的归一化因子为 `2 × (n − 3)`，这是最大 RF 距离（两棵树完全不同的情况）。
**请勿使用 `len(total) − 3`** 或 `n − 3` 作为分母，二者都不正确。

### 7.3 Random Seeds (随机种子)

| Seed 范围 | 数据 | 用途 |
|-----------|------|------|
| 100 – 129 | Indel-rich (INDELible) | Table 1, L0-1 vs MSA+ML |
| 130 – 229 | Clean JC69 | Clean benchmark (100 seeds) |
| 200 – 229 | Clean JC69 | Multi-k ensemble (30 seeds) |

---

## 8. Citation

If you use this code or data in your research, please cite:

> [Authors]. Information-Matched Multi-Level Inference: A General Framework
> for Scalable Phylogenetics. *Nucleic Acids Research Methods*, 2026.

DOI: [10.5281/zenodo.20746742](https://doi.org/10.5281/zenodo.20746742)

---

## 9. License

MIT License. See the source code repository for details.
