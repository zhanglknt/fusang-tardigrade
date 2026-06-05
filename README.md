# Fusang: Tardigrade Edition

**Fast Alignment-Free Phylogenetic Inference using Spaced k-mers**

Fusang (Tardigrade Edition) is a scalable, alignment-free phylogenetic inference framework that reconstructs phylogenetic trees directly from unaligned sequences using spaced k-mer features. It supports datasets with **10,000+ taxa** and runs in **seconds to minutes**, without requiring multiple sequence alignment (MSA).

> **Why "Tardigrade"?** Like the extremotolerant water bear, Fusang tolerates insertion/deletion mutations that cause MSA methods to fail — delivering robust trees where alignment-based tools degrade.

## Key Features

- **Alignment-free**: No MSA required — works directly on FASTA files
- **Spaced k-mers**: Uses gapped k-mer patterns (e.g., `gap1`, `gap2`) that significantly outperform contiguous k-mers
- **Scalable**: O(n²) distance matrix + FastME tree building; handles 10,000+ taxa
- **Fast**: ~2 seconds for 200 taxa, ~225 seconds for 10,000 taxa
- **Indel-robust**: Naturally handles insertion/deletion mutations where MSA methods degrade
- **Refinement**: Optional BME (Balanced Minimum Evolution) NNI refinement improves accuracy

## Quick Start

### Installation

```bash
git clone https://github.com/zhanglab/Fusang.git
cd Fusang/Fusang-main
pip install -r requirements.txt  # Python ≥3.8, Biopython, numpy
```

**Optional**: Install [FastME](http://www.atgc-montpellier.fr/fastme/) for faster tree building (recommended):
```bash
# Download fastme binary and place in bench_tools/fastme/
# Windows: bench_tools/fastme/fastme.exe
# Linux/macOS: bench_tools/fastme/fastme
```

### Basic Usage

```bash
python fusang.py --input sequences.fasta --output tree.nwk
```
> `fusang_v2.py` also works (backward-compatible alias).

With custom k-mer parameters:
```bash
python fusang.py --input sequences.fasta --output tree.nwk \
    --kmer_k 5 --kmer_gap gap2 --tree_method fastme
```

## Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--input`, `-i` | Input FASTA file | (required) |
| `--output`, `-o` | Output tree file (Newick format) | (required) |
| `--kmer_k` | k-mer length (4 or 5 recommended) | auto-selected by n |
| `--kmer_gap` | Gapped k-mer pattern: `none`, `gap1`, `gap2`, `gap3`, `gap4` | auto-selected by n |
| `--tree_method` | Tree building: `fastme` (recommended) or `nj` | `fastme` |
| `--mode` | Mode: `auto`, `default`, `refine` | `auto` |
| `--max_group` | Max taxa per group for divide-and-conquer | 200 |
| `--threads`, `-t` | Number of threads | 4 |

### Auto mode logic

| n (taxa) | Mode | Description |
|-----------|------|-------------|
| n < 20 | `default` | Pure FastME, no refinement |
| 20 ≤ n < 100 | `default` | Refinement not beneficial at small n |
| n ≥ 100 | `refine` | FastME + BME NNI refinement |

## Performance Benchmarks

### Accuracy (nRF vs true tree, n=200, L=500bp, sub_rate=0.05)

| Method | nRF ↓ | Time (s) | Requires MSA? |
|--------|--------|-----------|----------------|
| RAxML-NG | 0.013 | 30 | Yes |
| FastTree2 | 0.009 | 4 | Yes |
| **Fusang (default)** | **0.015** | **2** | **No** |
| **Fusang (refine)** | **0.013** | **2** | **No** |

### Indel Robustness (n=200, indel_rate=0.02)

| Method | nRF ↓ |
|--------|--------|
| FastTree2 | 0.096 |
| **Fusang** | **≈0.013** |

Fusang **outperforms MSA methods by 47%** under realistic indel conditions.

### Scalability (L=500bp)

| n | Time (s) | Memory (MB) |
|---|-----------|---------------|
| 200 | 2 | ~50 |
| 1,000 | 5 | ~200 |
| 10,000 | 225 | ~4,000 |

## Spaced k-mers: Key Innovation

Fusang uses **spaced k-mers** (gapped patterns) instead of contiguous k-mers. For example:
- Contiguous k=5: positions `[0,1,2,3,4]`
- Spaced k=5,gap2: positions `[0,1,2,4,6]` (skip position 3)

**Optimal gap pattern by dataset size:**
| n (taxa) | Optimal gap | Rationale |
|-----------|--------------|-----------|
| 50 | gap1 | Short-range signal |
| 100-200 | gap2 | Balanced |
| 500+ | gap2~none | Signal saturation |

**Publications:**
> Zhang L, et al. (2026). Spaced k-mer features for alignment-free phylogenetic inference. *In preparation.*

## Output Format

Fusang outputs Newick-formatted trees:
```
((seq1:0.12,seq2:0.08):0.05,(seq3:0.10,seq4:0.11):0.04);
```

Branch lengths are in substitution distance units.

## File Structure

```
Fusang/
├── Fusang-main/
│   ├── fusang.py             # Main entry point (Tardigrade Edition)
│   ├── fusang_v2.py         # Backward-compatible alias (same as fusang.py)
│   ├── kmer_distance.py      # k-mer distance computation
│   ├── fastme_backend.py     # FastME integration
│   ├── calc_nrf_simple.py    # nRF distance calculator
│   ├── gen_test_data.py      # Test data generator
│   ├── requirements.txt      # Python dependencies
│   └── benchmark_L/         # L-variant benchmarks
├── bench_tools/              # Benchmark tools (FastTree2, IQ-TREE2, etc.)
└── README.md
```

## Citation

If you use Fusang (Tardigrade Edition) in your research, please cite:

```bibtex
@article{zhang2026fusang,
  title={Fast alignment-free phylogenetic inference using spaced k-mers},
  author={Zhang, Li and colleagues},
  journal={In preparation},
  year={2026}
}
```

## License

MIT License — see [LICENSE](LICENSE) file for details.

## Contact

- **Issues**: [GitHub Issues](https://github.com/zhanglab/Fusang/issues)
- **Email**: zhanglknt@example.com

---

**Fusang v1** (deep learning-based, 4-40 taxa) is in `Fusang-v1/` directory for reference.
