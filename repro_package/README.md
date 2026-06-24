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
```

### Hardware Requirements

- **Minimum**: 4 GB RAM, 2 CPU cores (supports n <= 2,000)
- **Recommended**: 16 GB RAM, 4 CPU cores (supports n <= 10,000)
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
numba>=0.58
```

Optional external tools (for full reproduction):
- MAFFT v7+ (https://mafft.cbrc.jp/alignment/software/)
- FastTree2 v2.2.0+ (http://www.microbesonline.org/fasttree/)
- FastME v2.1.6+ (http://www.atgc-montpellier.fr/fastme/)
- IQ-TREE2 v2.2+ (http://www.iqtree.org/)
- INDELible v1.03+ (http://abacus.gene.ucl.ac.uk/software/indelible/)

### Package Structure

```
repro_package/
|-- README.md                    # This file
|-- requirements.txt             # Python dependencies
|-- data/                        # All benchmark data + models
|   |-- benchmark_mhl_n200_30seeds.csv       # MHL benchmark
|   |-- indel_benchmark_seeds100_129.csv     # n=200 indel (30 seeds)
|   |-- indel_benchmark_n500.csv             # n=500 indel
|   |-- indel_benchmark_n1000.csv            # n=1000 indel
|   |-- indel_benchmark_MASTER_REBUILT.csv   # All seeds merged
|   |-- benchmark_clean_seeds130_229.csv     # JC69 seeds 130-229
|   |-- benchmark_multik_seeds200_229.csv    # Multi-k ensemble
|   |-- scalability_results.json            # Runtime/memory scaling
|   |-- scalability_nrf.json                # nRF at scale
|   |-- training_data.pkl                   # 844 classifier samples
|   |-- boundary_rf.pkl                     # Pre-trained RF classifier
|-- fusang/                      # Core source code
|   |-- __init__.py
|   |-- fusang_v2.py             # Main IMMI pipeline (CLI)
|   |-- fusang_v4_dahp_v1.py     # DAHP-enhanced pipeline
|   |-- fusang_mhl_main.py       # MHL entry point
|   |-- kmer_distance.py         # k-mer frequency computation
|   |-- kmer_optimized.py        # Optimized k-mer routines
|   |-- tree_simulation.py       # Coalescent tree + JC69 simulator
|   |-- fusang_minhash_wrapper.py # MinHash distance wrapper
|   |-- calc_nrf_simple.py       # nRF computation utility
|   |-- af_competitor_methods.py # Competitor implementations
|   |-- ensemble.py              # Multi-k ensemble methods
|   |-- fusang.py                # Core types and utilities
|   |-- fusang_mhl/              # MHL sub-package (10 modules)
|       |-- config.py
|       |-- mlh_utils.py
|       |-- level0_feature.py
|       |-- level1_distance.py
|       |-- level2_partition.py
|       |-- level3_msa_ml.py
|       |-- merger.py
|       |-- boundary_classifier.py
|       |-- training_data_generator.py
|       |-- __init__.py
|-- scripts/                     # Standalone scripts
|   |-- generate_missing_data.py     # Generate clean benchmark data
|   |-- generate_multik.py           # Multi-k ensemble benchmark
|   |-- scalability_nrf.py           # nRF at scale
|   |-- benchmark_mhl.py             # MHL benchmark runner
|   |-- train_boundary_classifier.py
|   |-- benchmark_competitors.py
|-- bench_tools/                 # External tool setup
|   |-- README.md                # Download instructions
|-- notebooks/                   # Jupyter notebooks (see below)
`-- figures/                     # Generated figures
```

### IMPORTANT: Simulator Discrepancy

**The bundled simulator (`tree_simulation.py`) uses Kingman coalescent + JC69 model.**
The manuscript's indel-rich benchmarks were generated with **INDELible** (birth-death + GTR+Gamma + indels).
This means:

1. **Clean (no-indel) benchmarks** can be fully reproduced with `scripts/generate_missing_data.py`
2. **Indel-rich benchmarks** require INDELible (external tool) — see `bench_tools/README.md`
3. nRF values from JC69 data will differ from the paper's INDELible-based values

### Reproducing Results

#### Generate Clean Benchmark Data (n=200, seeds 130-229)
```bash
cd scripts
python generate_missing_data.py --start-seed 130 --end-seed 229 --n 200 \
    --output ../data/benchmark_clean_seeds130_229.csv
```
Note: uses JC69 simulator. For INDELible-based data, see section below.

#### Multi-k Ensemble Benchmark (n=200, seeds 200-229)
```bash
cd scripts
python generate_multik.py --start-seed 200 --end-seed 229 --n 200 \
    --output ../data/benchmark_multik_seeds200_229.csv
```

#### Scalability with nRF
```bash
cd scripts
python scalability_nrf.py --scales 200 500 1000 2000 --seeds-per-scale 5 \
    --output ../data/scalability_nrf.json
```

#### Table 1: L0-1 vs MSA+ML accuracy
Requires MAFFT + FastTree2 binaries (see `bench_tools/README.md`).
```bash
cd scripts
python benchmark_mhl.py --seeds 100-129 --no-l3 --methods nj,mhl,ft2 \
    --output ../data/benchmark_n200.csv
```

#### SwissTree Validation (Table 5)
Requires SwissTree data from AFproject (https://afproject.org/).
The script at `scripts/benchmark_competitors.py` can process SwissTree datasets
once downloaded.

#### IQ-TREE2 Comparison (Requires external tool)
```bash
# After installing IQ-TREE2 (see bench_tools/README.md):
iqtree2 -s alignment.fasta -m MFP -nt AUTO
```

#### Indel-rich Data Generation (Requires INDELible)
Indel-rich benchmarks (n=200, 500, 1000, indel rate scan) require INDELible.
See `bench_tools/README.md` for installation and the original INDELible control
files used for the paper's benchmarks.

#### 16S rRNA Validation (Requires external data)
16S sequences from SILVA/Greengenes databases with trusted reference trees.
See manuscript Supplementary Materials for data sources.

### Running Fusang from Command Line

```bash
# Build a tree from FASTA (L0-1, direct NJ)
python fusang/fusang_v2.py -i input.fasta -o output.nwk -m default \
    --auto-group --kmer-k 5 --kmer-gap gap2

# Multi-k ensemble mode (DAHP-V3)
python fusang/fusang_v4_dahp_v1.py -i input.fasta -o output.nwk --v3

# Full MHL pipeline (Level 0-1-2-3)
python fusang/fusang_mhl_main.py -i input.fasta -o output.nwk --all-levels
```

### License

MIT License. See the source code repository for details.

### Citation

If you use this code or data in your research, please cite:

> [Authors]. Information-Matched Multi-Level Inference: A General Framework
> for Scalable Phylogenetics. *Nucleic Acids Research*, 2026.
