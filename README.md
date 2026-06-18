# Fusang: Tardigrade Edition

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/downloads/)
[![GitHub Release](https://img.shields.io/github/v/release/zhanglknt/fusang-tardigrade)](https://github.com/zhanglknt/fusang-tardigrade/releases)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20746742.svg)](https://doi.org/10.5281/zenodo.20746742)

**Fast Alignment-Free Phylogenetic Inference using Spaced k-mers**

Fusang (Tardigrade Edition) is a scalable, alignment-free phylogenetic inference framework that reconstructs phylogenetic trees directly from unaligned sequences using spaced k-mer features. It supports datasets with **10,000+ taxa** and runs in **seconds to minutes**, without requiring multiple sequence alignment (MSA).

> **Why "Tardigrade"?** Like the extremotolerant water bear, Fusang tolerates insertion/deletion mutations that cause MSA methods to fail — delivering robust trees where alignment-based tools degrade.

## Key Features

- 🧬 **Alignment-free**: No MSA required — works directly on FASTA files
- 🎯 **Spaced k-mers**: Uses gapped k-mer patterns (gap1/gap2) outperforming contiguous k-mers under indels
- ⚡ **Scalable**: Handles 10,000+ taxa (~70 seconds, ~609 MB RAM)
- 🛡️ **Indel-robust**: Outperforms IQ-TREE2 GTR by **1.8× (p<0.001)** on indel-rich data
- 🔬 **IMMI framework**: Information-Matched Multi-level Inference — selects the optimal inference level per dataset
- 🌐 **Web server**: Included Flask app for browser-based access

## Quick Start

### Installation

```bash
git clone https://github.com/zhanglknt/fusang-tardigrade.git
cd fusang-tardigrade
pip install -r requirements.txt
```

Or via conda:

```bash
conda env create -f environment.yml
conda activate fusang
```

### Basic Usage

```bash
python fusang_v2.py -i sequences.fasta -o tree.nwk
```

With custom parameters:

```bash
python fusang_v2.py -i sequences.fasta -o tree.nwk \
    -k 5 -g 2 -d cosine -m nj
```

### Web Server

```bash
python fusang_webapp.py
# Open http://localhost:5001 in browser
```

## Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-i` / `--input` | Input FASTA file | (required) |
| `-o` / `--output` | Output tree file (Newick) | (required) |
| `-k` | k-mer length (4–9) | `5` |
| `-g` | Gap pattern: `0`=none, `1`=gap1, `2`=gap2 | `2` |
| `-d` | Distance metric: `cosine`, `euclidean` | `cosine` |
| `-m` | Tree method: `nj`, `fastme` | `nj` |
| `-t` | Number of threads | `4` |

## Performance Benchmarks

### Accuracy (nRF, lower is better — n=200, sub=0.05, indel=0.02)

| Method | nRF ↓ | Time | Requires MSA? |
|--------|--------|------|----------------|
| IQ-TREE2 GTR | 0.147±0.027 | ~2 min | Yes |
| FastTree2 (MAFFT) | 0.084±0.012 | ~5 s | Yes |
| Co-phylog | 0.419±0.025 | ~1 s | No |
| KmerCosine k=5 | 0.099±0.017 | <1 s | No |
| **Fusang (k=5,gap2)** | **0.112±0.020** | **<2 s** | **No** |
| **Fusang multi-k** | **0.105±0.021** | **<3 s** | **No** |

> Fusang outperforms IQ-TREE2 GTR by **1.8× (p<0.001, d=3.1)** on indel-rich data.

### Scalability

| n taxa | Time (s) | Memory (MB) |
|--------|-----------|-------------|
| 200 | ~2 | ~50 |
| 1,000 | ~8 | ~180 |
| 5,000 | ~32 | ~340 |
| 10,000 | **~70** | **~609** |

### SwissTree Protein Gene Trees (AFproject benchmark, 11 families)

| Method | nRF ↓ |
|--------|--------|
| **Fusang (k=4,gap1)** | **0.239±0.118** |
| Co-phylog k=11 | 0.433±0.076 |

Fusang achieves **1.8× better accuracy** (p=0.014, d=1.13) on real protein families.

## IMMI Framework

Fusang implements the **IMMI (Information-Matched Multi-level Inference)** framework, which automatically selects the optimal phylogenetic inference level based on the information content of the input data:

| Level | Method | Best for |
|-------|--------|---------|
| L0 | k-mer cosine distance + NJ | n<200, high indel rate |
| L1 | Multi-k ensemble | n=200–500, moderate indels |
| L2 | DAHP-V1 selective MSA | n=500–2000, mixed signal |
| L3 | MSA + ML (FastTree2) | n<500, low indel rate |

```bash
python fusang_mhl_main.py -i sequences.fasta -o tree.nwk
```

## Multi-k Ensemble (DAHP-V3)

```bash
python fusang_v4_dahp_v1.py --v3 -i sequences.fasta -o tree.nwk
```

Ensemble over k=5,7,9 with cosine distance — improves nRF by **6.5%** (p=0.006) over single-k.

## File Structure

```
fusang-tardigrade/
├── fusang_v2.py              # Main entry point
├── fusang_v4_dahp_v1.py      # DAHP V1+V3 with multi-k ensemble
├── fusang_mhl_main.py        # IMMI framework entry point
├── kmer_distance.py          # k-mer distance computation
├── fastme_backend.py         # FastME integration
├── calc_nrf_simple.py        # nRF accuracy calculator
├── fusang_webapp.py          # Flask web server
├── fusang_mhl/               # IMMI MHL package
│   ├── level0_kmer.py
│   ├── level1_multik.py
│   ├── level2_dahp.py
│   ├── level3_msa_ml.py
│   ├── merger.py
│   ├── boundary_classifier.py
│   └── models/               # Pre-trained boundary RF classifier
├── af_competitor_methods.py  # Competitor implementations (Co-phylog, etc.)
├── benchmark_competitors.py  # Benchmark runner
├── environment.yml           # Conda environment
├── requirements.txt          # pip dependencies
├── run_webapp.sh             # Linux/macOS server startup
├── run_webapp.bat            # Windows server startup
├── DEPLOYMENT.md             # Deployment guide
├── NAR_MANUSCRIPT_*.md       # Manuscript drafts
├── benchmark_n200_indel_30seeds.csv  # Benchmark data
├── benchmark_n500_indel_30seeds.csv
├── benchmark_n1000_indel_30seeds.csv
├── scalability_results.json
└── Figure*.pdf               # All manuscript figures
```

## Reproducibility

All benchmark data and figure-generation scripts are included. To reproduce main results:

```bash
# Generate all figures
python generate_figure1.py
python generate_figure2.py
# ...etc.

# Run benchmark (requires FastTree2)
python benchmark_competitors.py --n 200 --seeds 30 --output benchmark_n200_new.csv
```

## Citation

If you use Fusang in your research, please cite:

```bibtex
@article{kong2026fusang,
  title={Fast alignment-free phylogenetic inference using spaced k-mers and the IMMI framework},
  author={Kong, Lei and Zhang, Li},
  journal={Nucleic Acids Research},
  year={2026},
  note={Under review}
}
```

## License

MIT License — see [LICENSE](LICENSE) for details.

## Contact

- **Issues**: [GitHub Issues](https://github.com/zhanglknt/fusang-tardigrade/issues)
- **Corresponding author**: Li Zhang (ORCID: 0000-0002-0698-0754)
