# IMMI Framework — Fusang: Tardigrade Edition
## Reproducibility Package

This package contains all code, data, and Jupyter notebooks needed to
reproduce the results in:

> **"Information-Matched Multi-Level Inference: A General Framework for Scalable Phylogenetics"**
> Submitted to *Nucleic Acids Research Methods*, 2026.

### Quick Start

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Verify installation
python -c "from fusang.fusang_v2 import main; print('OK')"

# 3. Run the main benchmark
cd notebooks
jupyter notebook 01_main_benchmark.ipynb
```

### Hardware Requirements

- **Minimum**: 4 GB RAM, 2 CPU cores (supports n ≤ 2,000)
- **Recommended**: 16 GB RAM, 4 CPU cores (supports n ≤ 10,000)
- **Full benchmark**: 32 GB RAM recommended
- **OS**: Windows, Linux, or macOS with Python 3.9+

### Software Dependencies

```
numpy>=1.24
scipy>=1.10
scikit-learn>=1.3
biopython>=1.81
matplotlib>=3.7
seaborn>=0.12
pandas>=2.0
jupyter>=1.0
psutil>=5.9
```

Optional (for Level 3 MSA+ML refinement):
- MAFFT v7+ (https://mafft.cbrc.jp/alignment/software/)
- FastTree2 v2.2.0+ (http://www.microbesonline.org/fasttree/)

### Package Structure

```
repro_package/
├── README.md                    # This file
├── requirements.txt             # Python dependencies
├── build_notebooks.py           # Notebook builder script
├── data/                        # All benchmark data + models
│   ├── benchmark_mhl_n200_30seeds.csv
│   ├── indel_benchmark_seeds100_129.csv
│   ├── indel_benchmark_n500.csv
│   ├── indel_benchmark_n1000.csv
│   ├── indel_benchmark_MASTER_REBUILT.csv
│   ├── training_data.pkl         # 844 classifier training samples
│   ├── boundary_rf.pkl           # Pre-trained RF classifier
│   └── scalability_results.json  # Scalability test output
├── notebooks/                   # Jupyter notebooks
│   ├── 01_main_benchmark.ipynb           # Table 1 + Figure 3
│   ├── 02_multik_ensemble.ipynb          # Table 3 + Figure 4
│   ├── 03_swisstree_validation.ipynb     # Table 5 + Figure 5
│   ├── 04_boundary_classifier.ipynb      # Table 6 + Figure 6
│   ├── 05_scalability_demo.ipynb         # Table 7 + Figure 7
│   └── 06_generate_figures.ipynb         # All figures
├── fusang/                      # Core source code
│   ├── fusang_v2.py             # Main IMMI pipeline (CLI)
│   ├── fusang_mhl_main.py       # MHL entry point
│   ├── kmer_distance.py         # k-mer frequency computation
│   ├── kmer_optimized.py        # Optimized k-mer routines
│   ├── tree_simulation.py       # Coalescent tree + sequence simulator
│   ├── af_competitor_methods.py # Competitor method implementations
│   ├── ensemble.py              # Multi-k ensemble methods
│   └── fusang_mhl/              # MHL sub-package (10 modules)
│       ├── config.py
│       ├── mlh_utils.py
│       ├── level0_feature.py
│       ├── level1_distance.py
│       ├── level2_partition.py
│       ├── level3_msa_ml.py
│       ├── merger.py
│       ├── boundary_classifier.py
│       ├── training_data_generator.py
│       └── __init__.py
├── scripts/                     # Standalone scripts
│   ├── benchmark_mhl.py         # MHL benchmark runner
│   ├── train_boundary_classifier.py
│   ├── benchmark_competitors.py
│   ├── benchmark_swisstree.py
│   ├── benchmark_multik_ensemble.py
│   └── merge_benchmarks.py
└── figures/                     # Generated figures (after running notebooks)
```

### Reproducing Key Results

#### Table 1: L0–1 vs MSA+ML accuracy
```bash
cd scripts
python benchmark_mhl.py --seeds 100-129 --no-l3 --methods nj,mhl,ft2 \
    --output ../data/benchmark_n200.csv
```
Then run `notebooks/01_main_benchmark.ipynb`.

#### Table 3: Multi-k ensemble
```bash
cd scripts
python benchmark_multik_ensemble.py --seeds 100-129 --n 200 --output ../data/multik_ensemble.csv
```
Then run `notebooks/02_multik_ensemble.ipynb`.

#### Table 5: SwissTree validation
```bash
cd scripts
python benchmark_swisstree.py
```
Then run `notebooks/03_swisstree_validation.ipynb`.

#### Table 6: Boundary classifier
The pre-trained classifier (`data/boundary_rf.pkl`) and training data (`data/training_data.pkl`) are included. To retrain:
```bash
cd scripts
python train_boundary_classifier.py --data ../data/training_data.pkl --output ../data/ --n-samples 2000
```
Then run `notebooks/04_boundary_classifier.ipynb`.

#### Table 7: Scalability demonstration
```bash
cd fusang
python -c "
from tree_simulation import make_coalescent_tree, simulate_seqs
# Generate data at target sizes and run fusang_v2.py
"
```
Then run `notebooks/05_scalability_demo.ipynb`.

### Running Fusang from Command Line

```bash
# Build a tree from FASTA
python fusang/fusang_v2.py -i input.fasta -o output.nwk -m default \
    --auto-group --kmer-k 5 --kmer-gap gap2

# Multi-k ensemble mode
python fusang/fusang_v2.py -i input.fasta -o output.nwk --v3

# Full MHL pipeline (Level 0–1–2–3)
python fusang/fusang_mhl_main.py -i input.fasta -o output.nwk --all-levels
```

### Data Generation

To generate new simulated datasets for benchmarking:

```python
from fusang.tree_simulation import make_coalescent_tree, simulate_seqs

import numpy as np
code = {0: 'A', 1: 'T', 2: 'C', 3: 'G'}

# Generate n=200 taxa, L=500 bp, substitution rate 0.05
root, leaves = make_coalescent_tree(200, seed=42)
seqs = simulate_seqs(root, 200, 500, 0.05, seed=42)

# Convert to FASTA
for i in range(200):
    seq_str = ''.join(code[b] for b in seqs[i])
    print(f'>t{i+1:04d}')
    print(seq_str)
```

### License

MIT License. See the source code repository for details.

### Citation

If you use this code or data in your research, please cite:

> [Authors]. Information-Matched Multi-Level Inference: A General Framework
> for Scalable Phylogenetics. *Nucleic Acids Research*, 2026.

### Contact

[Author contact information]
[GitHub repository URL]
[Web server URL]
