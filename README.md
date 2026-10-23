# Fusang: Tardigrade Edition

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/downloads/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20746742.svg)](https://doi.org/10.5281/zenodo.20746742)

**Alignment-free phylogenetic inference from unaligned FASTA in seconds — k-mer frequency vector cosine distances that retain accuracy under indel-rich sequence evolution.**

Fusang reconstructs phylogenetic trees directly from unaligned sequences: spaced k-mer frequency vectors → cosine distance → Neighbor-Joining / FastME. No multiple sequence alignment, no model selection, no manual parameter tuning — adaptive defaults pick k-mer size and gap pattern from your dataset size, and a multi-k ensemble removes the need to choose k at all.

> **Why "Tardigrade"?** Like the extremotolerant water bear, Fusang tolerates insertion/deletion mutations that cause MSA-based methods to degrade.

## Features

| Feature | What it means for you |
|---|---|
| Alignment-free | Feed it raw FASTA — no MSA step, no alignment errors propagated into the tree |
| Spaced k-mers | Gapped patterns (e.g. `1001001001001` for k=5, gap2) tolerate indels better than contiguous k-mers |
| Adaptive parameters | k and gap auto-selected from n (n≤100 → k=4,gap1; n>100 → k=5,gap2) |
| Multi-k ensemble (`--v3`) | Averages distances over k=5,7,9 — no manual k selection; auto-rescues k=5 saturation at genomic scale |
| Two-scale pipeline | Simplified direct pipeline for n≤500; DCM decomposition + EPA grafting for n>500 |
| Boundary classifier (research) | Random-forest model (`boundary_rf.pkl`) flags homogeneous vs. structured datasets — part of the `fusang_mhl_main.py` research framework, not of the shipped `fusang_v2.py` inference pipeline |
| FastME backend | BIONJ + balanced NNI; bundled native Windows (`fastme_bin/fastme.exe`) and Linux binaries — no WSL or manual install needed |
| Scalable | 10,000 taxa in ~70 s on a 4-core workstation |

## 30-second quick start

```bash
git clone https://github.com/zhanglknt/fusang-tardigrade.git
cd fusang-tardigrade
pip install -r requirements.txt   # numpy, numba, biopython, scipy

# Build a tree — one command, sensible defaults chosen automatically
python fusang_v2.py -i sequences.fasta -o tree.nwk
```

Input: a FASTA file of unaligned DNA or protein sequences. Output: a Newick tree. That is the whole workflow.

Optional: multi-k ensemble (recommended for genome-scale or mixed-divergence data):

```bash
python fusang_v4_dahp_v1.py sequences.fasta --v3 --output tree.nwk
```

## CLI reference

### `fusang_v2.py` — main entry point

| Argument | Description | Default |
|---|---|---|
| `-i`, `--input` | Input FASTA file (unaligned) | required |
| `-o`, `--output` | Output tree file (Newick) | required |
| `-m`, `--mode` | `auto`, `default` (NJ+DCM), `refine` (NJ+DCM+BME), `full-dl` | `auto` (n≥100 → refine) |
| `-d`, `--distance_method` | `kmer` or `p-distance` | `kmer` |
| `--kmer_k` | k-mer length | auto: n≤100 → 4, else 5 |
| `--kmer_gap` | `none`, `gap1`, `gap2`, `gap3`, `gap4` | auto: n≤100 → gap1, else gap2 |
| `--tree_method` | `nj` (recommended for k-mer distances) or `fastme` (faster; may degrade k-mer distance accuracy) | `nj` |
| `--auto_group_method` | DCM strategy: `simple`, `nj_centroid`, `epa_improved` | `nj_centroid` |
| `--simple` / `--no-simple` | Force simplified pipeline on/off (auto: simplified for n≤500, DCM above) | auto |
| `--max_group` | Max taxa per DCM group (auto-scaled) | 200 |
| `--overlap` | DCM group overlap ratio | 0.15 |
| `-t`, `--threads` | Threads | 4 |
| `--use_minhash` | MinHash LSH coarse clustering (scales toward 50K+ taxa) | off |

### `fusang_v4_dahp_v1.py` — multi-k ensemble / DAHP

| Argument | Description | Default |
|---|---|---|
| `fasta` | Input FASTA (positional) | required |
| `--v3` | Multi-k distance ensemble | off |
| `--v2` | DAHP-V2 backbone refinement | off |
| `--ks` | Ensemble k values (comma-separated) | `5,7,9` |
| `--fusion` | `average` or `weighted` distance fusion | `average` |
| `-k`, `--gap` | Single-k and gap pattern (non-ensemble modes) | k=5, gap2 |
| `--output` | Output Newick file | stdout/default |

### `fusang_mhl_main.py` — IMMI multi-level pipeline (L0–L3)

```bash
python fusang_mhl_main.py sequences.fasta -o tree.nwk [--k 5 --gap gap2] [--no-l2] [--no-l3] [--boundary-model path]
```

The IMMI framework selects an inference level per dataset: L0 k-mer+NJ (small/indel-rich), L1 multi-k ensemble, L2 DAHP-V2, L3 MSA+ML within reliable clusters (requires MAFFT + FastTree2).

## When to use Fusang — and when not to

| Data type | Verdict |
|---|---|
| Indel-rich gene families (DNA or protein) | ✅ Primary use case — matches or beats MSA+ML at a fraction of the cost |
| Protein families (e.g. SwissTree/AFproject benchmarks) | ✅ Strong accuracy vs. alignment-free competitors |
| Organelle genomes (mitochondrial/chloroplast), deep divergence | ✅ Use the multi-k ensemble (`--v3`) |
| Large unaligned sets (1,000–10,000+ taxa) | ✅ DCM pipeline scales near-linearly |
| 16S / other structured rRNA genes | ❌ Secondary-structure covariation violates k-mer independence assumptions — use alignment-based methods |
| Strain-level, recombination-dominated genomes (e.g. *E. coli*/*Shigella*) | ❌ Recombination scrambles k-mer signal; k-mer distances are not additive here |

If you are unsure, run the MHL entry point: the boundary classifier will tell you whether your dataset looks homogeneous or structured.

## Performance

| n taxa | Time | Notes |
|---|---|---|
| 200 | 2.3 s (NJ, optimised) | single gene, 4-core workstation |
| 500 | 3.8 s | simplified pipeline, FastME BIONJ+BNNI |
| 1,000 | 5.2 s | DCM+EPA pipeline (default routing for n>500), FastME |
| 10,000 | 69.8 s | DCM pipeline, 4-core workstation, ~0.6 GB RAM |

Accuracy (TRUE-relative nRF vs. simulated ground truth, indel rate=0.02, n=200, L=500 bp, 130 seeds; lower is better): Fusang 0.199 ± 0.030 vs FastTree2 (MAFFT+ML) 0.205 ± 0.030 — numerically better but borderline (paired Wilcoxon p=0.0507, not significant at α=0.05); IQ-TREE2 GTR 0.210 ± 0.031 (significantly worse than Fusang, p=9.2×10⁻⁵). The multi-k ensemble (k=5,7,9) improves on the default pipeline (0.194 vs 0.199, p=7.7×10⁻³). Note the direction reverses at scale: at n=1000 indel FastTree2 is significantly better (0.212 vs 0.226, p=1.5×10⁻⁴). See the manuscript and `repro_package/` for full benchmark tables.

## Documentation

- [docs/TUTORIAL.md](docs/TUTORIAL.md) — end-to-end tutorial: installation, your first tree, adaptive parameters, multi-k, large datasets, FAQ
- [docs/EXAMPLES.md](docs/EXAMPLES.md) — three worked examples (gene family, mitochondrial genomes, 10,000 taxa)
- [repro_package/](repro_package/) — reproducibility package for all manuscript results
- [DEPLOYMENT.md](DEPLOYMENT.md) — web server deployment (`python fusang_webapp.py` → http://localhost:5001)

## Citation

If you use Fusang, please cite the archived software and the manuscript:

```bibtex
@software{fusang_zenodo,
  title  = {Fusang: Tardigrade Edition},
  author = {Zhang, Li and Wang, Xiaowo and Li, Yixue},
  doi    = {10.5281/zenodo.20746742},
  url    = {https://doi.org/10.5281/zenodo.20746742}
}
```

Manuscript (under review at *Nucleic Acids Research*): "Fusang: Tardigrade Edition — K-mer Frequency Vector Alignment-Free Phylogenetic Inference under Indel-Rich Sequence Evolution".

## License

MIT — see [LICENSE](LICENSE).

## Contact

- Issues: [GitHub Issues](https://github.com/zhanglknt/fusang-tardigrade/issues)
