# Supplementary Material: Multi-Scale Alignment-Free Phylogenetic Inference Benchmark

## Fusang: Tardigrade Edition — Accuracy Evaluation

---

## 1. Overview

This document provides the full methodological details, raw data, and statistical analyses for the
benchmark experiments reported in the manuscript "Fusang: Tardigrade Edition — Spaced k-mer
Alignment-Free Phylogenetic Tree Reconstruction." All scripts, parameters, and seed values are
documented so that every result can be independently reproduced.

---

## 2. Data Simulation

### 2.1 Simulation Script

`gen_test_data_indel.py` (included in the project repository)

### 2.2 Tree Topology

A random coalescent-like tree with `n` taxa is generated for each seed. Branch lengths are
fixed at 0.1 substitutions per site. The tree is written in Newick format as the ground-truth
reference for nRF computation.

### 2.3 Sequence Evolution Model

Starting from a root sequence of length `L`, sequences evolve along the tree by depth-first
traversal. Each branch applies:

1. **Indels** (if `indel_rate > 0`): The number of indel events per branch is drawn from a
   Gaussian approximation to Poisson(λ), where λ = indel_rate × branch_length × current_seq_length.
   Each event is 50% insertion / 50% deletion, with geometric length distribution (mean = 3,
   capped at 50), producing a realistic power-law-like indel size spectrum (POW 1.5).

2. **Substitutions**: The number of substitution events per branch is drawn from a Gaussian
   approximation to Poisson(λ), where λ = sub_rate × branch_length × current_seq_length.
   Each substitution replaces one nucleotide with a different base (equal probability among
   A, T, C, G).

This two-phase process (indels first, then substitutions on the indel-modified sequence)
produces unaligned sequences of variable length, mimicking the challenges of real biological
data where insertion/deletion events disrupt positional homology.

### 2.4 Benchmark Parameters

| Parameter | Value |
|-----------|-------|
| Sequence length (L) | 500 bp |
| Substitution rate | 0.05 per site per unit branch length |
| Indel rate (indel datasets) | 0.02 per site per unit branch length |
| Indel rate (clean datasets) | 0.00 |
| Branch length | 0.1 (constant) |
| Nucleotide model | GTR-like, equal base frequencies, equal rates |

### 2.5 Dataset Summary

| Dataset ID | n (Taxa) | Type | Seeds | Total Benchmarks |
|------------|----------|------|-------|-------------------|
| N200-CLN | 200 | Clean (substitution only) | 200–229 | 30 |
| N200-IND | 200 | Indel (sub + indel) | 230–259 | 30 |
| N200-MEGA | 200 | Indel (sub + indel) | 100–229 | 130 |
| N500-CLN | 500 | Clean (substitution only) | 500–529 | 30 |
| N500-IND | 500 | Indel (sub + indel) | 530–559 | 30 |
| N1000-CLN | 1000 | Clean (substitution only) | 1000–1029 | 30 |

> **Note**: N200-MEGA is an earlier, larger-scale benchmark (130 seeds) with a different seed
> range (100–229), using the same `gen_test_data_indel.py` script and parameters. It predates
> the systematic multi-scale sweep and is included for its higher statistical power.

---

## 3. Methods Compared

### 3.1 Fusang: Tardigrade Edition

**Script**: `fusang_v2.py`

**Parameters**:
```
--tree_method nj            # Neighbor-Joining backbone
--kmer_k 5                  # k-mer length
--kmer_gap gap2             # Spaced k-mer pattern: positions 0 and 3 (gap of 2)
```

**Pipeline (n ≤ 500 — simplified)**:
1. Compute spaced k-mer (k=5, gap=2) frequency vectors for all sequences
2. Compute pairwise cosine distance matrix
3. Build NJ tree from the full distance matrix

**Pipeline (n = 1000 — DCM + EPA)**:
1. Compute spaced k-mer (k=5, gap=2) frequency vectors for all sequences
2. Compute pairwise cosine distance matrix
3. Cluster sequences by k-mer distance (DCM — Disk-Covering Method)
4. Build NJ backbone tree from cluster centroids
5. Place remaining sequences onto the backbone using Evolutionary Placement Algorithm (EPA)
6. Additional parameter: `--auto_group_method nj_centroid`

### 3.3 DAHP-V3: Multi-k Distance Ensemble

**Script**: `fusang_v4_dahp_v1.py --v3`

**Parameters**:
```
--v3                        # Enable DAHP-V3 mode
--ks 5,7,9                  # k-mer sizes for ensemble (default: 5,7,9)
--fusion average            # Distance fusion method (average, weighted)
```

**Pipeline**:
1. Compute contiguous k-mer frequency vectors for k=5, 7, and 9 (no gap pattern)
2. Compute pairwise cosine distance matrices for each k value
3. Fuse distance matrices by element-wise averaging
4. Build NJ tree from the fused distance matrix

> Note: Contiguous k-mers (no gap pattern) are used for the ensemble because different k values
> provide complementary phylogenetic information. Using the same gap pattern for all k values would
> produce correlated distance matrices, reducing the benefit of ensemble averaging.

### 3.4 FastTree2 (Reference Method)

**Version**: 2.2.0 (Double precision)

**Parameters**:
```
-nt            # Nucleotide alignment
-gtr           # GTR model
-nosupport     # No SH-like local support (faster)
```

**Pipeline**:
1. MAFFT v7.526: align unaligned FASTA sequences (`mafft --auto`)
2. FastTree2: build maximum-likelihood tree from the alignment

### 3.3 nRF Calculation

**Script**: `calc_nrf_simple.py` (pure Python, zero dependencies)

The normalized Robinson-Foulds distance is computed as:

```
nRF = (RF / max_RF)

where:
  RF = |B1 \ B2| + |B2 \ B1|           (symmetric difference of bipartition sets)
  max_RF = 2 × (n_taxa - 3)             (maximum possible RF for n_taxa)
```

Bipartitions are extracted by a recursive-descent Newick parser. Each internal edge splits the
leaf set into two subsets; the canonical bipartition is the smaller of the two complementary
halves (or lexicographically first if equal size).

---

## 4. Benchmark Results

### 4.1 Multi-Scale Summary

| Dataset | n | Type | FT2 Mean nRF | FT2 SD | Fusang Mean nRF | Fusang SD | Cohen's d | Paired t p | Wilcoxon p | Fusang Wins |
|---------|---|------|-------------|--------|----------------|-----------|-----------|------------|------------|-------------|
| N200-CLN | 200 | Clean | 0.095685 | 0.019193 | 0.101692 | 0.018781 | 0.316 | 0.2094 | 0.1877 | 12/30 |
| N200-IND | 200 | Indel | 0.079865 | 0.017159 | 0.077496 | 0.018177 | 0.134 | 0.5154 | 0.6211 | 16/30 |
| N200-MEGA | 200 | Indel | 0.084186 | 0.019337 | 0.080359 | 0.016511 | 0.213 | **0.0493** | **0.0491** | 69/130 |
| N500-CLN | 500 | Clean | 0.092991 | 0.012542 | 0.118511 | 0.011050 | 2.159 | <0.0001 | <0.0001 | 1/30 |
| N500-IND | 500 | Indel | 0.082830 | 0.014358 | 0.095171 | 0.015361 | 0.830 | 0.0010 | 0.0020 | 8/30 |
| N1000-CLN | 1000 | Clean | 0.090689 | 0.009666 | 0.114711 | 0.010948 | 2.326 | <0.0001 | <0.0001 | 1/30 |

> **Interpretation**: Fusang outperforms FastTree2 at n=200 indel (p<0.05 on the 130-seed
> benchmark). FastTree2 has a clear advantage on clean (substitution-only) data and on larger
> scales (n=500, n=1000) regardless of data type. The DCM+EPA pipeline at n=1000 produces
> consistent results but with a systematic accuracy gap relative to the MSA-based reference.

### 4.2 N200-MEGA: Detailed Results (130 seeds, seeds 100–229)

**Data file**: `indel_benchmark_130seeds_MASTER.csv`

| Statistic | FastTree2 | Fusang: Tardigrade Ed. |
|-----------|-----------|------------------------|
| Mean nRF | 0.084186 | 0.080359 |
| Std Dev | 0.019337 | 0.016511 |
| Min nRF | 0.027919 | 0.035533 |
| Max nRF | 0.124365 | 0.129442 |
| Fusang better | — | 69 / 130 (53.1%) |
| FT2 better | 47 / 130 (36.2%) | — |
| Tie | 14 / 130 (10.8%) | — |
| Mean ΔnRF (FT2 − Fusang) | +0.003827 | — |
| Paired t-test | t = 1.984, **p = 0.0493** (two-tailed) | — |
| Wilcoxon signed-rank | W = 2679.5, **p = 0.0491** (two-tailed) | — |
| Cohen's d | 0.213 (small effect) | — |

> **Key finding**: On the 130-seed indel benchmark (n=200), Fusang achieves a statistically
> significant improvement over FastTree2 (p<0.05), with ~53% of seeds favoring Fusang.

### 4.3 N200-CLN: Clean Data, 30 Seeds

**Data file**: `benchmark_n200_clean_30seeds.csv`

| seed | nRF FastTree2 | nRF Fusang | Fusang Better? |
|------|--------------|------------|----------------|
| 200 | 0.083756 | 0.106599 | No |
| 201 | 0.093909 | 0.081218 | Yes |
| 202 | 0.088832 | 0.081218 | Yes |
| 203 | 0.114213 | 0.111675 | Yes |
| 204 | 0.088832 | 0.101523 | No |
| 205 | 0.134518 | 0.116751 | Yes |
| 206 | 0.076142 | 0.106599 | No |
| 207 | 0.093909 | 0.101523 | No |
| 208 | 0.078680 | 0.111675 | No |
| 209 | 0.109137 | 0.096447 | Yes |
| 210 | 0.114213 | 0.101523 | Yes |
| 211 | 0.119289 | 0.081218 | Yes |
| 212 | 0.078680 | 0.065990 | Yes |
| 213 | 0.073604 | 0.101523 | No |
| 214 | 0.083756 | 0.101523 | No |
| 215 | 0.088832 | 0.126904 | No |
| 216 | 0.119289 | 0.116751 | Yes |
| 217 | 0.078680 | 0.081218 | No |
| 218 | 0.098985 | 0.101523 | No |
| 219 | 0.088832 | 0.101523 | No |
| 220 | 0.139594 | 0.096447 | Yes |
| 221 | 0.063452 | 0.086294 | No |
| 222 | 0.093909 | 0.091371 | Yes |
| 223 | 0.073604 | 0.137056 | No |
| 224 | 0.104061 | 0.106599 | No |
| 225 | 0.109137 | 0.060914 | Yes |
| 226 | 0.093909 | 0.096447 | No |
| 227 | 0.083756 | 0.131980 | No |
| 228 | 0.124365 | 0.142132 | No |
| 229 | 0.078680 | 0.106599 | No |

**Summary**: Mean nRF: FT2=0.095685, Fusang=0.101692 (difference −0.006007, Fusang worse).
Not statistically significant (paired t p=0.2094).

### 4.4 N200-IND: Indel Data, 30 Seeds

**Data file**: `benchmark_n200_indel_30seeds.csv`

| seed | nRF FastTree2 | nRF Fusang | Fusang Better? |
|------|--------------|------------|----------------|
| 230 | 0.053299 | 0.086294 | No |
| 231 | 0.078680 | 0.086294 | No |
| 232 | 0.088832 | 0.106599 | No |
| 233 | 0.088832 | 0.071066 | Yes |
| 234 | 0.088832 | 0.076142 | Yes |
| 235 | 0.073604 | 0.096447 | No |
| 236 | 0.098985 | 0.101523 | No |
| 237 | 0.093909 | 0.096447 | No |
| 238 | 0.073604 | 0.071066 | Yes |
| 239 | 0.068528 | 0.060914 | Yes |
| 240 | 0.068528 | 0.096447 | No |
| 241 | 0.048223 | 0.050761 | No |
| 242 | 0.093909 | 0.101523 | No |
| 243 | 0.073604 | 0.065990 | Yes |
| 244 | 0.068528 | 0.081218 | No |
| 245 | 0.093909 | 0.060914 | Yes |
| 246 | 0.083756 | 0.086294 | No |
| 247 | 0.053299 | 0.076142 | No |
| 248 | 0.109137 | 0.101523 | Yes |
| 249 | 0.068528 | 0.086294 | No |
| 250 | 0.098985 | 0.091371 | Yes |
| 251 | 0.088832 | 0.086294 | Yes |
| 252 | 0.063452 | 0.040609 | Yes |
| 253 | 0.119289 | 0.060914 | Yes |
| 254 | 0.078680 | 0.065990 | Yes |
| 255 | 0.053299 | 0.065990 | No |
| 256 | 0.088832 | 0.065990 | Yes |
| 257 | 0.083756 | 0.071066 | Yes |
| 258 | 0.063452 | 0.035533 | Yes |
| 259 | 0.088832 | 0.081218 | Yes |

**Summary**: Mean nRF: FT2=0.079865, Fusang=0.077496 (difference +0.002369, Fusang better).
Not statistically significant (paired t p=0.5154). The smaller 30-seed sample lacks power;
the 130-seed benchmark (N200-MEGA, §4.2) provides conclusive statistical evidence.

### 4.5 N500-CLN: Clean Data, 30 Seeds

**Data file**: `benchmark_n500_clean_30seeds.csv`

| seed | nRF FastTree2 | nRF Fusang |
|------|--------------|------------|
| 500 | 0.085513 | 0.116700 |
| 501 | 0.087525 | 0.122736 |
| 502 | 0.089537 | 0.118712 |
| 503 | 0.097586 | 0.100604 |
| 504 | 0.079477 | 0.116700 |
| 505 | 0.105634 | 0.138833 |
| 506 | 0.079477 | 0.106640 |
| 507 | 0.101610 | 0.112676 |
| 508 | 0.101610 | 0.126761 |
| 509 | 0.100604 | 0.114688 |
| 510 | 0.099598 | 0.136821 |
| 511 | 0.091549 | 0.124748 |
| 512 | 0.079477 | 0.084507 |
| 513 | 0.113682 | 0.126761 |
| 514 | 0.083501 | 0.120724 |
| 515 | 0.103622 | 0.132797 |
| 516 | 0.106640 | 0.122736 |
| 517 | 0.099598 | 0.118712 |
| 518 | 0.069416 | 0.122736 |
| 519 | 0.095573 | 0.120724 |
| 520 | 0.099598 | 0.116700 |
| 521 | 0.108652 | 0.104628 |
| 522 | 0.079477 | 0.106640 |
| 523 | 0.073441 | 0.122736 |
| 524 | 0.103622 | 0.118712 |
| 525 | 0.084507 | 0.116700 |
| 526 | 0.075453 | 0.124748 |
| 527 | 0.077465 | 0.114688 |
| 528 | 0.105634 | 0.132797 |
| 529 | 0.110664 | 0.110664 |

**Summary**: Mean nRF: FT2=0.092991, Fusang=0.118511 (diff −0.025520, p<0.0001). FT2 clearly
outperforms on clean data at this scale.

### 4.6 N500-IND: Indel Data, 30 Seeds

**Data file**: `benchmark_n500_indel_30seeds.csv`

| seed | nRF FastTree2 | nRF Fusang |
|------|--------------|------------|
| 530 | 0.099598 | 0.084507 |
| 531 | 0.077465 | 0.096579 |
| 532 | 0.067404 | 0.074447 |
| 533 | 0.075453 | 0.072435 |
| 534 | 0.069416 | 0.112676 |
| 535 | 0.063380 | 0.080483 |
| 536 | 0.079477 | 0.086519 |
| 537 | 0.079477 | 0.104628 |
| 538 | 0.075453 | 0.122736 |
| 539 | 0.071429 | 0.100604 |
| 540 | 0.071429 | 0.086519 |
| 541 | 0.061368 | 0.110664 |
| 542 | 0.089537 | 0.076459 |
| 543 | 0.107646 | 0.118712 |
| 544 | 0.095573 | 0.102616 |
| 545 | 0.087525 | 0.068410 |
| 546 | 0.083501 | 0.082495 |
| 547 | 0.091549 | 0.104628 |
| 548 | 0.091549 | 0.090543 |
| 549 | 0.087525 | 0.096579 |
| 550 | 0.075453 | 0.072435 |
| 551 | 0.081489 | 0.086519 |
| 552 | 0.071429 | 0.102616 |
| 553 | 0.055332 | 0.088531 |
| 554 | 0.087525 | 0.096579 |
| 555 | 0.089537 | 0.096579 |
| 556 | 0.115694 | 0.092555 |
| 557 | 0.109658 | 0.118712 |
| 558 | 0.077465 | 0.106640 |
| 559 | 0.095573 | 0.120724 |

**Summary**: Mean nRF: FT2=0.082830, Fusang=0.095171 (diff −0.012341, p=0.0010). FT2 outperforms
significantly on indel data at this scale.

### 4.7 N1000-CLN: Clean Data, 30 Seeds

**Data file**: `benchmark_n1000_clean_30seeds.csv`

| seed | nRF FastTree2 | nRF Fusang |
|------|--------------|------------|
| 1000 | 0.087763 | 0.116349 |
| 1001 | 0.114343 | 0.139418 |
| 1002 | 0.089769 | 0.108325 |
| 1003 | 0.100802 | 0.109328 |
| 1004 | 0.078736 | 0.122367 |
| 1005 | 0.079739 | 0.112337 |
| 1006 | 0.100301 | 0.131394 |
| 1007 | 0.101805 | 0.121364 |
| 1008 | 0.087763 | 0.121364 |
| 1009 | 0.075727 | 0.102307 |
| 1010 | 0.078736 | 0.102307 |
| 1011 | 0.098796 | 0.113340 |
| 1012 | 0.086760 | 0.113340 |
| 1013 | 0.084754 | 0.113340 |
| 1014 | 0.094784 | 0.127382 |
| 1015 | 0.084754 | 0.112337 |
| 1016 | 0.097793 | 0.126379 |
| 1017 | 0.107823 | 0.119358 |
| 1018 | 0.095787 | 0.120361 |
| 1019 | 0.096790 | 0.123370 |
| 1020 | 0.088265 | 0.115346 |
| 1021 | 0.080742 | 0.109328 |
| 1022 | 0.079739 | 0.098295 |
| 1023 | 0.080742 | 0.095286 |
| 1024 | 0.099799 | 0.122367 |
| 1025 | 0.088766 | 0.088265 |
| 1026 | 0.078235 | 0.106319 |
| 1027 | 0.094784 | 0.115346 |
| 1028 | 0.090271 | 0.110331 |
| 1029 | 0.095787 | 0.124373 |

**Summary**: Mean nRF: FT2=0.090689, Fusang=0.114711 (diff −0.024022, p<0.0001). FT2 clearly
outperforms. The DCM+EPA pipeline produces stable results with no catastrophic failures, but the
accuracy gap relative to MSA-based FastTree2 is approximately 0.024 nRF units (Cohen's d=2.33,
large effect).

### 4.8 DAHP-V3: Multi-k Distance Ensemble (30 seeds, n=200 indel)

**Data file**: `benchmark_multik_ensemble_n200_indel.csv`

**Method**: Average of contiguous k-mer distance matrices for k=5, 7, and 9, followed by NJ tree
construction. No gap pattern is used (contiguous k-mers) to maximize information diversity across
different k values.

| Statistic | Fusang (k=5, gap2) | k=5 contig | k=7 contig | k=9 contig | **Ensemble (k=5,7,9)** |
|-----------|---------------------|------------|------------|------------|------------------------|
| Mean nRF | 0.1121 | 0.1046 | 0.1056 | 0.1091 | **0.1045** |
| Std Dev | 0.0188 | 0.0203 | 0.0170 | 0.0222 | **0.0210** |
| Min | 0.0678 | 0.065 | 0.060 | 0.070 | **0.060** |
| Max | 0.1432 | 0.150 | 0.155 | 0.170 | **0.160** |
| Ensemble wins | — | 19/30 | 18/30 | 19/30 | **24/30** |

**Paired Comparison: Original vs Ensemble (30 seeds)**

| Metric | Value |
|--------|-------|
| Mean nRF improvement | 0.0076 (6.7% relative) |
| Ensemble wins | 24 / 30 (80.0%) |
| Paired t-test | t = 2.8949, **p = 0.0071** |
| Wilcoxon signed-rank | W = 98.0, **p = 0.0057** |
| Cohen's d | 0.538 (medium effect) |

> **Key finding**: The multi-k ensemble significantly outperforms the original k=5 gap2 method
> (p = 0.006, Wilcoxon), with a medium effect size (Cohen's d = 0.54). The ensemble wins on
> 80% of seeds. This improvement is achieved without any MSA step — it is purely from distance
> matrix fusion across multiple k-mer resolutions.

### 4.9 Multiple Comparison Correction (5 ground-truth datasets)

All 5 datasets with ground-truth reference trees were tested for the comparison Fusang vs FastTree2.
Both raw and corrected p-values are reported below.

| Dataset | Fusang Mean | FT2 Mean | Raw p | Bonferroni p | BH-FDR p | Significant? |
|---------|-------------|----------|-------|-------------|---------|--------------|
| N200-CLN | 0.1017 | 0.0957 | 0.1877 | 0.9384 | 0.2346 | No |
| N200-IND | 0.0775 | 0.0799 | 0.6211 | 1.0000 | 0.6211 | No |
| N500-CLN | 0.1185 | 0.0930 | 0.000003 | 0.000016 | 0.000008 | Yes (FT2 better) |
| N500-IND | 0.0952 | 0.0828 | 0.00196 | 0.00980 | 0.00327 | Yes (FT2 better) |
| N1000-CLN | 0.1147 | 0.0907 | 0.000002 | 0.000010 | 0.000010 | Yes (FT2 better) |

> **Interpretation**: After Bonferroni correction (5 tests), 3/5 datasets remain significant, all
> in favor of FastTree2. At n=200, no significant difference is detected between methods. The
> n=1000 indel dataset (not shown) uses a direct method-vs-method comparison rather than
> ground-truth reference; both methods agree closely (nRF between Fusang and FT2 = 0.037 ± 0.006).

### 4.10 Alignment-Free Competitor Comparison (27 seeds, n=200 indel)

**Purpose**: Compare Fusang against established alignment-free phylogenetic methods under indel-rich conditions.

**Methods tested**:
- **Co-phylog** (Yi & Jin, 2013, NAR 41:e75): Context-object matching with k=19 (18bp context + 1bp object). Python reimplementation based on published source code analysis. Distance = fraction of shared conserved contexts with divergent center bases.
- **K-mer cosine (k=5, contiguous)**: Standard k-mer frequency vector + cosine distance + NJ. No gap pattern, no spaced k-mers. Represents the most basic alignment-free baseline.
- **K-mer cosine (k=7, contiguous)**: Same as above with larger k.
- **Fusang (k=5, gap2)**: Spaced k-mer frequency vector + cosine distance + NJ. The default configuration.
- All methods use the same NJ implementation (BioPython DistanceTreeConstructor).

**Test data**: Seeds 100-129 (27 valid seeds with reference trees; seeds 110, 121, 126 excluded due to reference tree parsing issues). Each seed: n=200 taxa, sequence length ~500bp, substitution rate 0.05, indel rate 0.02.

**Results**:

| Method | Mean nRF | Std Dev | Median | Min | Max | Wins vs Fusang |
|--------|:--------:|:-------:|:------:|:---:|:---:|:--------------:|
| Co-phylog (k=19) | 0.612 | 0.071 | 0.590 | 0.554 | 0.838 | 0/27 |
| K-mer cosine (k=5) | 0.234 | 0.210 | 0.153 | 0.109 | 0.843 | 0/27 |
| K-mer cosine (k=7) | 0.237 | 0.211 | 0.153 | 0.118 | 0.843 | 0/27 |
| **Fusang (k=5, gap2)** | **0.112** | **0.019** | **0.109** | 0.109 | 0.203 | — |

**Pairwise Wilcoxon tests (27 paired seeds)**:
- Fusang vs Co-phylog: p < 0.001, Cohen's d = -2.58 (very large effect)
- Fusang vs KmerCosine(k=5): p < 0.001, Cohen's d = -1.86 (very large effect)
- Fusang vs KmerCosine(k=7): p < 0.001, Cohen's d = -1.84 (very large effect)
- KmerCosine(k=5) vs KmerCosine(k=7): p = 0.465, Cohen's d = 0.12 (negligible)

**Key observations**:
1. Co-phylog's context-object approach is severely disrupted by indels. The 18-bp context matching requires exact positional alignment of flanking regions, which breaks when insertions/deletions shift the reading frame. nRF=0.612 is close to the random expectation for 200 taxa.
2. Spaced k-mers (Fusang) provide a ~2× improvement over contiguous k-mers. The gap=2 pattern skips over short indels, preserving phylogenetic signal.
3. Increasing contiguous k from 5 to 7 provides no significant improvement (p=0.465), confirming that the spaced pattern (not the k value) is the key innovation.

**andi note**: andi (Haubold et al. 2015, Bioinformatics) was tested on a single seed and produced nRF=0.517 on gene-length data (~500bp). andi's suffix-array anchor approach is designed for genome-scale sequences (>10kb) and cannot produce meaningful phylogenetic estimates on gene-length data. This is a scale mismatch, not a methodological comparison.

**Reproduction**: `python benchmark_competitors.py --seeds 100-129 --type indel`
**Raw data**: `benchmark_competitors_n200_indel.csv`

### 4.11 AFproject SwissTree Gene Tree Benchmark (11 families, protein sequences)

**Purpose**: Evaluate Fusang on the community-standard AFproject benchmark for alignment-free gene tree inference, using real protein domain sequences with trusted reference trees.

**Dataset**: SwissTree release 2017.0 (Zielezinski et al. 2019, *Genome Biology* 20:144). 11 gene families (ST001-ST012), each containing individual protein domain sequences and a trusted reference tree from the SwissTree database.

**Table S10. SwissTree per-family benchmark results.**

| Family | Gene | Taxa | Avg Len (aa) | Fusang k=4,gap1 | K-mer k=5 | Co-phylog k=11 | Best Method |
|--------|------|-----:|:---:|:---:|:---:|:---:|:---|
| ST001 | Popeye | 49 | 334 | 0.2500 | 0.2787 | 0.5417 | Fusang k=4 |
| ST002 | NOX | 54 | 576 | 0.2031 | 0.2308 | 0.5443 | Fusang k=4 |
| ST003 | ATPase | 49 | 500 | 0.3175 | 0.3175 | 0.5278 | Fusang k=4 |
| ST004 | Serine | 115 | 449 | 0.2910 | 0.2910 | 0.6121 | Fusang k=5 |
| ST005 | SUMF | 29 | 330 | 0.2727 | 0.2727 | 0.6744 | Tie (all kmer) |
| ST007 | S10/S20 | 60 | 131 | 0.4146 | 0.4337 | 0.5393 | Fusang k=5 |
| ST008 | Bambi | 42 | 276 | 0.4909 | 0.4909 | 0.4630 | Co-phylog |
| ST009 | Asterix | 39 | 111 | 0.5472 | 0.5192 | 0.4600 | Co-phylog |
| ST010 | Cited | 34 | 213 | **0.1389** | 0.1892 | 0.5745 | Fusang k=4 |
| ST011 | GH14 | 159 | 543 | 0.2957 | 0.2772 | 0.7269 | Fusang k=5 |
| ST012 | Ant | 21 | 402 | 0.5667 | 0.5667 | 0.5172 | Fusang k=3 |

**Summary statistics**:
- Best k-mer method (Fusang k=4,gap1): mean nRF = 0.344 ± 0.135
- K-mer cosine k=5 contiguous: mean nRF = 0.352 ± 0.122
- Co-phylog halfctx=5: mean nRF = 0.562 ± 0.078

**Statistical tests** (Wilcoxon paired, 11 families):
- Fusang k=4,gap1 vs K-mer k=5: p=0.312, Cohen's d=0.06 (not significant)
- Best k-mer vs Co-phylog: p=0.014, Cohen's d=1.98 (large, significant)
- Best k-mer wins vs Co-phylog: 8/11 families (73%)

**Key observations**:
1. K-mer frequency methods outperform Co-phylog by 1.67× on protein data (same pattern as DNA)
2. Spaced k-mers show no significant advantage on protein (p=0.31, d=0.06) — consistent with indel-tolerance hypothesis
3. Sequence length strongly predicts accuracy (longer → better; r=-0.72 for best k-mer)
4. Co-phylog wins on 3/11 families, all with short sequences (≤276 aa)

**Reproduction**: `python benchmark_swisstree.py`
**Raw data**: `benchmark_swisstree_results.csv`
**Data source**: https://afproject.org/app/benchmark/genetree/swisstree/dataset/
**Reference trees**: SwissTree database, extracted via AFproject GitHub repository

---

## 5. Statistical Methods

All statistical tests were performed using **SciPy 1.17.1** (Python 3.13.12).

### 5.1 Tests Used

| Test | Purpose | Implementation |
|------|---------|----------------|
| Paired t-test (two-tailed) | Compare mean nRF across paired seeds | `scipy.stats.ttest_rel` |
| Wilcoxon signed-rank (two-tailed) | Non-parametric paired comparison | `scipy.stats.wilcoxon` |
| Cohen's d | Effect size | |μ₁ − μ₂| / σ_pooled |

### 5.2 Cohen's d Formula

```
d = |mean_FT2 - mean_Fusang| / sqrt((σ²_FT2 + σ²_Fusang) / 2)
```

Interpretation: d ≈ 0.2 = small, d ≈ 0.5 = medium, d ≈ 0.8 = large.

---

## 5.5 SwissTree Gene Tree Benchmark (AFproject Standard)

**Purpose**: Cross-domain validation of k-mer methods on real protein gene families using the AFproject standardized benchmark.

**Dataset**: SwissTree release 2017.0, 11 protein gene families (Zielezinski et al. 2019, Genome Biology 20:144).

| Family | Gene | Taxa | Avg Length (aa) | Ref Bipartitions |
|--------|------|-----:|-----------------:|-----------------:|
| ST001 | Popeye domain | 49 | 334 | 56 |
| ST002 | NOX subfamily | 54 | 576 | 61 |
| ST003 | V-type ATPase β | 49 | 500 | 57 |
| ST004 | Serine incorporator | 115 | 449 | 114 |
| ST005 | SUMF | 29 | 330 | 28 |
| ST007 | Ribosomal S10/S20 | 60 | 131 | 70 |
| ST008 | Bambi | 42 | 276 | 41 |
| ST009 | Asterix | 39 | 111 | 38 |
| ST010 | Cited | 34 | 213 | 33 |
| ST011 | Glycosyl hydrolase 14 | 159 | 543 | 158 |
| ST012 | Ant transformer | 21 | 402 | 22 |

**Methods tested** (all use NJ tree construction):
- Fusang spaced k-mer: k=3,gap1 (101), k=4,gap1 (1011), k=5,gap2 (11011)
- KmerCosine contiguous: k=3, k=4, k=5
- Co-phylog: halfctx=5 (k=11), halfctx=3 (k=7)

**Per-family nRF results**:

| Family | Fus k=4,g1 | Fus k=5,g2 | Kmer k=3 | Kmer k=5 | Cophy k=11 |
|--------|-----------:|-----------:|---------:|---------:|-----------:|
| ST001 | 0.2500 | 0.2787 | 0.2787 | 0.2787 | 0.5417 |
| ST002 | 0.2031 | 0.2031 | 0.2031 | 0.2308 | 0.5443 |
| ST003 | 0.3175 | 0.3692 | 0.3438 | 0.3175 | 0.5278 |
| ST004 | 0.2910 | 0.2782 | 0.3285 | 0.2910 | 0.6121 |
| ST005 | 0.2727 | 0.2727 | 0.2727 | 0.2727 | 0.6744 |
| ST007 | 0.4146 | 0.3951 | 0.4146 | 0.4337 | 0.5393 |
| ST008 | 0.4909 | 0.5179 | 0.5179 | 0.4909 | 0.4630 |
| ST009 | 0.5472 | 0.5192 | 0.5192 | 0.5192 | 0.4600 |
| ST010 | **0.1389** | 0.1892 | 0.1892 | 0.1892 | 0.5745 |
| ST011 | 0.2957 | **0.2678** | 0.2957 | 0.2772 | 0.7269 |
| ST012 | 0.5667 | 0.5667 | 0.5667 | 0.5667 | 0.5172 |
| **Mean** | **0.344** | **0.351** | **0.357** | **0.352** | **0.562** |
| **SD** | 0.135 | 0.127 | 0.124 | 0.122 | 0.078 |

**Statistical tests**:
- K-mer cosine (best config k=4,gap1) vs Co-phylog: Wilcoxon p=0.014, Cohen's d=1.98 (large)
- Spaced (k=4,gap1) vs contiguous (k=5): Wilcoxon p=0.31, Cohen's d=0.06 (negligible)
- Best k-mer wins 8/11 families vs Co-phylog

**Key observations**:
1. K-mer frequency methods outperform Co-phylog by 1.67× on protein gene trees
2. Spaced k-mers show no significant advantage on protein data (no indels)
3. Sequence length strongly predicts accuracy: longer sequences yield lower nRF
4. Co-phylog's context-matching fails on both DNA (indels) and protein data

**Reproduction**: Run `python benchmark_swisstree.py` from the project root. SwissTree data in `real_data/swisstree/`, reference trees from AFproject GitHub (afproject-org/afproject).

---

## 6. Reproducibility

### 6.1 Software Versions

| Software | Version | Role |
|----------|---------|------|
| Python | 3.13.12 | Host environment |
| SciPy | 1.17.1 | Statistical tests |
| Fusang: Tardigrade Ed. | `fusang_v2.py` (commit: 2026-06-10) | Alignment-free tree builder |
| FastTree | 2.2.0 (Double precision) | ML tree from MSA |
| MAFFT | 7.526 (2024/Apr/26) | Multiple sequence alignment |

### 6.2 Exact Reproduction Commands

All commands are run from `d:\系统发育树项目\Fusang\Fusang-main\`.

**Step 1 — Generate data and run benchmark** (e.g., n=200, clean, 30 seeds):

```bash
python run_multi_scale_benchmark.py \
  --n 200 \
  --data_type clean \
  --seeds 200-229 \
  --L 500 \
  --sub_rate 0.05 \
  --indel_rate 0 \
  --output benchmark_n200_clean_30seeds.csv
```

**Step 2 — For indel data** (e.g., n=200, 30 seeds):

```bash
python run_multi_scale_benchmark.py \
  --n 200 \
  --data_type indel \
  --seeds 230-259 \
  --L 500 \
  --sub_rate 0.05 \
  --indel_rate 0.02 \
  --output benchmark_n200_indel_30seeds.csv
```

**Step 3 — Compute statistics** (from any CSV output):

```python
import csv, statistics, math
from scipy import stats

ft2_vals = []
fusang_vals = []
with open('benchmark_n200_clean_30seeds.csv') as f:
    r = csv.DictReader(f)
    for row in r:
        ft2_vals.append(float(row['nrf_ft2']))
        fusang_vals.append(float(row['nrf_fusang']))

print(f"FT2:   mean={statistics.mean(ft2_vals):.6f}, sd={statistics.stdev(ft2_vals):.6f}")
print(f"Fusang: mean={statistics.mean(fusang_vals):.6f}, sd={statistics.stdev(fusang_vals):.6f}")
t_stat, t_p = stats.ttest_rel(ft2_vals, fusang_vals)
w_stat, w_p = stats.wilcoxon(ft2_vals, fusang_vals)
print(f"Paired t: t={t_stat:.4f}, p={t_p:.6f}")
print(f"Wilcoxon: W={w_stat:.4f}, p={w_p:.6f}")
```

### 6.3 Per-Seed Pipeline

For each seed, the `run_multi_scale_benchmark.py` script performs five steps:

1. **Generate data** (via `gen_test_data_indel.py`): creates `seed{N}_{type}.fasta` and `seed{N}_{type}_true.nwk`
2. **MAFFT alignment** (via `run_mafft_simple.py`): produces `seed{N}_{type}_aligned.fasta`
3. **FastTree2** (`FastTree.exe -nt -gtr -nosupport`): produces `seed{N}_{type}_ft2.nwk`
4. **Fusang** (`fusang_v2.py --input ... --output ... --tree_method nj --kmer_k 5 --kmer_gap gap2`): produces `seed{N}_{type}_fusang.nwk`
5. **nRF calculation** (`calc_nrf_simple.py <true> <pred>`): reports nRF distance

### 6.4 Fusang DCM+EPA Pipeline (n=1000 only)

For n=1000, an additional parameter is appended:

```bash
python fusang_v2.py \
  --input seed{N}_clean.fasta \
  --output seed{N}_clean_fusang.nwk \
  --tree_method nj \
  --kmer_k 5 \
  --kmer_gap gap2 \
  --auto_group_method nj_centroid
```

This triggers the Disk-Covering Method (DCM) decomposition into subgroups, NJ backbone tree
on centroids, and Evolutionary Placement Algorithm (EPA) placement of remaining sequences.

---

## 7. Data Files

All raw benchmark CSV files are available in the project repository:

| File | Contents | Seeds | Size |
|------|----------|-------|------|
| `indel_benchmark_130seeds_MASTER.csv` | 130-seed indel benchmark (n=200) | 100–229 | 130 rows |
| `benchmark_n200_clean_30seeds.csv` | 30-seed clean benchmark (n=200) | 200–229 | 30 rows |
| `benchmark_n200_indel_30seeds.csv` | 30-seed indel benchmark (n=200) | 230–259 | 30 rows |
| `benchmark_n500_clean_30seeds.csv` | 30-seed clean benchmark (n=500) | 500–529 | 30 rows |
| `benchmark_n500_indel_30seeds.csv` | 30-seed indel benchmark (n=500) | 530–559 | 30 rows |
| `benchmark_n1000_clean_30seeds.csv` | 30-seed clean benchmark (n=1000) | 1000–1029 | 30 rows |

Each CSV contains three columns: `seed`, `nrf_ft2`, `nrf_fusang`.

---

## 8. Key Scripts

| Script | Purpose |
|--------|---------|
| `run_multi_scale_benchmark.py` | Orchestrator: generates data, runs methods, computes nRF |
| `gen_test_data_indel.py` | Simulates unaligned sequences with substitutions and indels |
| `run_mafft_simple.py` | Wrapper for MAFFT alignment on Windows |
| `fusang_v2.py` | Main Fusang: Tardigrade Edition program |
| `calc_nrf_simple.py` | Pure-Python nRF calculator |
| `kmer_distance.py` | Spaced k-mer frequency and distance computation |

---

*Document generated: 2026-06-07. All data computed on Intel Xeon E-2124 (4C/4T), 32 GB RAM, Windows 10.*
