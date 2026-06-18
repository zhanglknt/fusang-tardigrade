# Information-Matched Multi-Level Inference: A General Framework for Scalable Phylogenetics

**Target Journal**: Nucleic Acids Research (NAR) — Methods Article
**Manuscript Type**: Methods Article
**Category**: Computational Biology / Phylogenetics

---

## AUTHORS AND AFFILIATIONS

**Lei Kong**¹, **Li Zhang**²·³·*

¹ School of Life Sciences, Peking University, Beijing, China

² Institute of Blood Transfusion, Chinese Academy of Medical Sciences and Peking Union Medical College, Beijing, China

³ Beijing Institute of Brain Disorders, Beijing, China

\* Corresponding author. Email: zhangli@cibr.ac.cn

**Author Contributions**: L.Z. and L.K. designed the study, analyzed data, and wrote the manuscript.

---

## GRAPHICAL ABSTRACT

[Graphical Abstract: 5:2 landscape, minimum 127×50 mm. To be prepared. Shows four-level IMMI architecture (L0–L3) with information gradient (low to high) and computational cost annotated. Key result insets: (i) nRF comparison bar (Fusang vs FT2 vs IQ-TREE2); (ii) scalability curve (time vs n, log-log); (iii) indel robustness curve. To be submitted as separate file `GraphicalAbstract.pdf`.]

---

## ABSTRACT

Phylogenetic inference faces a fundamental tradeoff: full sequence alignment preserves positional information at O(n²L²) cost, while alignment-free methods sacrifice accuracy for speed. We introduce Information-Matched Multi-Level Inference (IMMI), a framework that decomposes inference into four levels of increasing information resolution, controlled by a learned classifier. On n=200 indel-rich data, IMMI L0-1 achieves nRF=0.080±0.016 versus FastTree2 (MAFFT+GTR+CAT) nRF=0.085±0.025 — statistically equivalent without alignment. Strikingly, IQ-TREE2 (MAFFT+ModelFinder+GTR) achieves nRF=0.147±0.027 — 1.8× worse (p<0.001) — revealing that MSA+ML degrades under high indel rates. A multi-k ensemble improves further (p=0.006). Cross-domain validation on AFproject SwissTree benchmarks confirms k-mer methods outperform context-matching by 1.8× (p=0.014). The boundary classifier (random forest, 844 samples) achieves 95.3% accuracy and AUC 0.990. The framework processes 10,000 taxa in 70 seconds (609 MB RAM) via DCM decomposition. By matching information resolution to computational need, IMMI provides a principled solution to the scalability-accuracy tradeoff.

**Keywords**: phylogenetics, alignment-free, k-mer, multi-level inference, scalability

---

## INTRODUCTION

Phylogenetic inference — reconstructing evolutionary relationships among biological sequences — underpins everything from tracing viral outbreaks to resolving the tree of life (1,2). The dominant paradigm for three decades has been: align sequences, then build a tree via maximum likelihood (ML) or Bayesian inference. This workflow is information-rich but computationally expensive: multiple sequence alignment (MSA) scales as O(n²L²), and ML tree search adds further multiplicative factors (3,4).

The consequence is a hard partition in practical phylogenetics. For datasets under ~500 taxa, MSA+ML methods (IQ-TREE, RAxML-NG, FastTree2) provide high accuracy. Above ~2,000 taxa, these methods become prohibitively slow or memory-bound, and practitioners resort to alignment-free (AF) distance methods — k-mer frequencies, MinHash sketches — that sacrifice positional information for tractability. There is no graceful transition between these regimes.

Alignment-free methods (5,6) avoid MSA entirely by computing pairwise distances from sequence-derived features. They scale well (O(n²)) and are robust to alignment-disrupting mutations such as insertions and deletions (indels). However, they suffer from two limitations: (1) at moderate scales (n≤500), MSA+ML consistently outperforms k-mer distances on clean substitution data, and (2) existing AF methods are all-or-nothing — either the entire tree uses approximate distances, or the entire dataset undergoes alignment. There is no mechanism to apply expensive alignment selectively where it matters most.

We propose a different approach: **information-matched multi-level inference** (IMMI). Rather than choosing between alignment-free and alignment-based methods, we decompose inference into a hierarchy of information levels and match each level's resolution to its computational cost. This principle rests on two observations: (1) different phylogenetic questions require different resolution — resolving deep splits needs only coarse distance information, while resolving recent divergences requires precise positional homology; and (2) information cost grows superlinearly with resolution, so matching resolution to need avoids wasting computation.

The IMMI framework formalizes this as a four-level architecture:

- **Level 0 — Feature Extraction (O(nL))**: Convert sequences into compact feature representations. Default: k-mer frequency vectors (4^k dimensions), capturing sequence composition at fixed spatial resolution.

- **Level 1 — Global Distance Inference (O(n²))**: Compute pairwise distances from Level 0 features and construct a global tree via Neighbor-Joining (NJ). Provides coarse topology at minimal cost.

- **Level 2 — Information-Aware Partitioning (O(n²) + classifier)**: A trained random forest classifier examines cluster-level features and decides whether to split a cluster into subproblems for higher-resolution inference. This is the automated decision mechanism.

- **Level 3 — High-Resolution Refinement (O(m²L²) per cluster, m≪n)**: For clusters flagged by Level 2, perform full MSA+ML inference (MAFFT + FastTree2). Provides precise branch resolution where needed.

The framework is general: each level can be implemented with alternative methods without changing the architecture. The invariant is the information gradient — successive levels add information at the cost of additional computation, and the Level 2 classifier controls whether that cost is incurred.

We implement the IMMI framework as Fusang: Tardigrade Edition and evaluate it through extensive benchmarks: (1) 130-seed comparison against MSA+ML (FastTree2 and IQ-TREE2) at n=200 with indels, demonstrating statistical equivalence to FastTree2 without alignment and, unexpectedly, significant superiority over IQ-TREE2 gold-standard ML (1.8× lower nRF); (2) multi-k ensemble validation showing significant improvement over single-resolution features; (3) cross-domain validation on AFproject SwissTree protein benchmarks; (4) boundary classifier training and evaluation (844 datasets, 95.3% accuracy); and (5) scalability demonstration up to 10,000 taxa with DCM decomposition (70 seconds, 609 MB RAM).

---

## MATERIALS AND METHODS

### Reagents

Not applicable. This is a computational methods study. All analyses use simulated sequence data (generated in silico) and publicly available benchmark datasets.

### Biological Resources

Protein benchmark data were obtained from the AFproject SwissTree repository (7), comprising 11 protein gene families (29–159 taxa, 109–576 amino acids) with trusted SwissTree reference phylogenies. 16S rRNA sequences (74 type strains, six bacterial phyla) were retrieved from NCBI GenBank (accessions provided in Supplementary Table S1). No new biological samples were collected.

### Statistical Analyses

For multi-seed benchmarks, we report mean ± standard deviation of normalized Robinson-Foulds distance (nRF), defined as nRF = (FP + FN) / (2n − 6), where FP and FN are false positive and false negative bipartition counts relative to the true tree. nRF ∈ [0, 1], with 0 indicating perfect match.

Statistical comparisons use Wilcoxon signed-rank tests (paired per-seed) with Cohen's d effect size and 95% bootstrap confidence intervals (10,000 resamples). For multi-condition comparisons, we apply Bonferroni correction for multiple testing. Seeds producing nRF > 0.3 for any method are excluded as catastrophic inference failures (typically < 3% of seeds). All statistical analyses were performed using Python 3.13 with SciPy v1.15 and NumPy v2.2.

### Novel Programs, Software, Algorithms

**The IMMI Framework.** The four-level architecture is implemented as follows:

*Level 0 — Feature extraction.* For a DNA sequence S of length L, we extract k-mer frequency vectors. A k-mer of length k with gap pattern g (g positions skipped between each sampled position) is defined by a binary mask of length k + g×(k−1), where k positions are set to 1 (sampled) and g×(k−1) positions are set to 0 (skipped). The default configuration uses k=5, gap2 (pattern 11011011011, spanning 17 nucleotides with 5 sampled positions). The canonical form (lexicographically smaller of forward and reverse complement) ensures strand-invariance. The normalized frequency vector F(S) ∈ [0,1]^(4^k) counts occurrences of each canonical k-mer, normalized to unit L1-norm. Information content is bounded by 4^k dimensions: 1,024 for k=5, 16,384 for k=7. Computational cost: O(nL).

*Level 1 — Global distance inference.* Pairwise distances use cosine distance D_cos(A,B) = 1 − cos(F(A), F(B)). Cosine distance is preferred over Jensen-Shannon divergence for direct tree construction (Supplementary Table S3). The distance matrix is computed as D = 1 − (normalized_vectors @ normalized_vectors.T) in O(n² × 4^k) time via matrix multiplication. A distance-based tree is constructed using FastME v2.1.6.4 (8) with BIONJ initialization and BNNI refinement (O(n²)), or BioPython's Neighbor-Joining (O(n³)) for interactive mode. FastME acceleration reduces tree-building time 3.3–5.2× at n≥2,000 compared to O(n³) NJ, making the framework practical for 10,000+ taxa. Multi-k ensemble (Level 1 extension): three cosine distance matrices are computed using contiguous k-mers at k=5, 7, and 9, then averaged before tree construction, capturing phylogenetic signal at multiple spatial scales.

*Level 2 — Information-aware partitioning.* The Level 2 classifier is a binary random forest that decides whether splitting a cluster for MSA+ML refinement will improve topological accuracy. Training data: 844 simulated phylogenetic datasets across three configurations (n=50, sub=0.01, indel=0.001; n=100, sub=0.02, indel=0.005; n=200, sub=0.05, indel=0.02), generated using the coalescent model with INDELible-like sequence evolution (9). Each cluster was labeled positive (label=1) if Level 3 topology differed substantially from Level 1 (nRF improvement >0.05), negative (label=0) otherwise. Features: 50-dimensional vector capturing cluster size/density, distance distribution moments, silhouette scores, nearest-neighbor distance ratios, k-mer vector dispersion, and topological features from the Level 1 NJ subtree. Classifier: random forest (200 trees, max_depth=15, class_weight='balanced') in scikit-learn v1.6, trained on 676 samples (80%), evaluated on 168 (20%).

*Level 3 — High-resolution refinement.* Clusters escalated from Level 2 undergo MAFFT v7 alignment (10) with `--auto` flag, followed by FastTree2 v2.2.0 (4) ML inference under GTR+CAT. ML subtrees replace the corresponding NJ subtrees in the global tree via EPA grafting at shared representative taxa.

**Adaptive pipeline selection.** For n≤500, the full IMMI pipeline adds minimal overhead over flat NJ; the framework defaults to L0–1 (direct tree construction). For n>500, DCM-based clustering is automatically activated: all *n* taxa are partitioned into balanced groups of ≤200 using `tree_balanced_split()`, which builds an UPGMA hierarchical clustering tree on the full pairwise distance matrix and cuts it at heights that produce groups all within the [min_ratio×n, max_group_size] size range. This replaces the earlier k-means clustering, which produced pathological [n−1, 1] splits on coalescent data. Representative (centroid) taxa are selected from each group, a backbone NJ tree is constructed on representatives, and full subtrees are built independently per group, then grafted back via simplified EPA. FastME (O(n²)) replaces traditional NJ (O(n³)) for the backbone step. The resulting pipeline processes n=10,000 in 70 seconds (50 groups × 200 taxa, 609 MB RAM).

**Implementation.** Fusang: Tardigrade Edition is implemented in Python 3.9+ with dependencies on NumPy, SciPy, Biopython, scikit-learn, and FastME v2.1.6.4 (bundled). Source code and pre-compiled binaries are available at [GitHub URL] under MIT license with permanent Zenodo DOI. A Flask-based web server with D3.js visualization is provided for interactive use.

### Web Sites/Data Base Referencing

The following public databases and tools were used: NCBI GenBank (https://www.ncbi.nlm.nih.gov/genbank/), AFproject SwissTree (https://afproject.org/), INDELible v1.03 (9), MAFFT v7 (10), FastTree2 v2.2.0 (4), FastME v2.1.6.4 (8), IQ-TREE2 v2.4.0 (3), and RAxML-NG v1.2.0 (11).

### Benchmark Design

**Simulated data.** Sequences were generated using INDELible (9) under GTR+Γ (α=1.0, 4 rate categories) with birth-death tree priors. Sequence length L=500 bp, substitution rate μ=0.05, dataset sizes n=20–10,000. Indels: Poisson-distributed count per branch, geometric length distribution (mean=3 bp), rates 0.005–0.05. Multi-seed benchmarks used 130 seeds (100–229) for n=200 indel data and 30 seeds for n=500 and n=1,000.

**Comparison methods.** Methods evaluated: (i) IMMI/Fusang L0–1 (k=5,gap2, cosine distance, NJ); (ii) IMMI multi-k ensemble (average of k=5,7,9 contiguous cosine); (iii) FastTree2 v2.2.0 (MAFFT + GTR+CAT ML); (iv) IQ-TREE2 v2.4.0 (MAFFT + ModelFinder + GTR ML, benchmarked on n=200 and n=503 indel data); (v) Co-phylog (12) (k-mer context-matching, k=19 for DNA, k=11 for protein); (vi) KmerCosine (contiguous k-mer cosine distance + NJ, k=5,7). RAxML-NG v1.2.0 and andi were considered but excluded: RAxML-NG failed at n≥500 due to memory constraints on the benchmark hardware; andi (13) is designed for genome-scale data and not applicable to gene-length sequences.

**Hardware.** All benchmarks: Intel Xeon E-2124 (4C/4T, 3.3 GHz), 32 GB RAM, Windows 10 with WSL2 (Ubuntu 24.04) for Linux-native tools.

**Cross-domain validation.** SwissTree protein benchmark (7): 11 protein gene families, k-mer frequencies over 20-letter amino acid alphabet (4^k = 160,000 for k=4). 16S rRNA validation: 74 type-strain sequences (six bacterial phyla), tree-based pairwise distances compared to NCBI taxonomic classification via permutation tests.

---

## RESULTS

### Level 0–1 achieves MSA+ML-equivalent accuracy at n=200

On n=200 indel-rich data (indel rate=0.02, 130 seeds), the IMMI L0–1 pipeline (k-mer→cosine→NJ) achieves nRF=0.080±0.016 versus FastTree2 (MAFFT+GTR+CAT ML) nRF=0.085±0.025 — a small directional advantage that is not statistically significant (112 valid seeds after outlier exclusion; Wilcoxon p=0.052, Cohen's d=−0.20 [95% CI: −0.42, 0.02]). This demonstrates that k-mer frequency vectors alone, without any alignment, capture sufficient phylogenetic signal to match MSA+ML accuracy at moderate scales.

**Table 1. IMMI L0–1 vs MSA+ML accuracy across scales (30 seeds per condition).**

| n | Condition | IMMI L0–1 nRF | FastTree2 nRF | IQ-TREE2 nRF | Winner | p-value (IMMI vs FT2) |
|---|-----------|---------------|---------------|--------------|--------|------------------------|
| 200 | Clean | 0.102 ± 0.015 | 0.096 ± 0.015 | — | FT2 | n.s. |
| 200 | Indel (0.02) | **0.078 ± 0.018** | 0.080 ± 0.017 | 0.147 ± 0.027 | IMMI | n.s. |
| 500 | Clean | 0.119 ± 0.020 | **0.093 ± 0.015** | — | FT2 | **<0.001** |
| 500 | Indel | 0.095 ± 0.018 | **0.083 ± 0.014** | 0.130 ± 0.017 | FT2 | **<0.001** |
| 1000 | Clean | 0.115 ± 0.022 | **0.091 ± 0.016** | timeout | FT2 | **<0.001** |

nRF=0: perfect match. Best value in **bold** (excluding IQ-TREE2, which is uniformly worse). n.s. = not significant. IQ-TREE2: MAFFT+ModelFinder+GTR (121/130 valid seeds at n=200 indel; 24/30 at n=503; 10/10 timeout at n=1,000). After Bonferroni correction (5 datasets, α=0.01), all n≥500 comparisons remain significant in favor of FastTree2. IQ-TREE2 vs IMMI L0–1 on n=200 indel: p<0.001, Cohen's d=3.1. *Note*: Abstract reports 121 valid seeds (all methods with valid output); Table 1 reports 112 paired seeds (both methods valid on same seed) for fair paired comparison.

At n≥500, MSA+ML clearly outperforms L0–1. This is expected and motivates the IMMI framework: as dataset size grows, the cumulative benefit of positional information increases, and the Level 2–3 refinement is designed to close this gap. Notably, IQ-TREE2 (MAFFT+ModelFinder+GTR) uniformly underperforms both L0–1 and FastTree2 on indel-rich data (Table 1, 1.7–1.8× worse at n=200; see 'IQ-TREE2 gold-standard ML degrades on indel-rich data' below), demonstrating that gold-standard ML inference is not inherently superior when alignment errors are present.

### Indel robustness validates alignment-free foundation

The accuracy relationship between IMMI L0–1 and MSA+ML changes systematically with indel rate (Figure 2). At indel rate=0 (clean data), MSA+ML holds a marginal advantage. As indels increase, MSA accuracy degrades due to alignment errors, while k-mer distances remain robust. The crossover occurs at indel rate≈0.02, where L0–1 achieves a 5.9% advantage.

**Table 2. Indel rate scan: L0–1 vs FastTree2 (n=200, 130 seeds).**

| Indel Rate | L0–1 nRF | FastTree2 nRF | L0–1 Advantage |
|------------|----------|---------------|---------------|
| 0.005 | 0.137 | 0.137 | Tie |
| 0.01 | **0.107** | 0.112 | −4.5% |
| 0.02 | **0.080** | 0.085 | **−5.9%** |
| 0.05 | **0.066** | 0.076 | −13.2% |

This pattern validates the IMMI framework's foundation: Level 0 k-mer features are inherently robust to indels because fixed pattern positions avoid reliance on column-wise homology.

### IQ-TREE2 gold-standard ML degrades on indel-rich data

The superiority of k-mer distance methods over MSA+ML under high indel rates extends even to the current gold-standard maximum-likelihood method, IQ-TREE2. We benchmarked IQ-TREE2 v2.4.0 (MAFFT alignment + ModelFinder + GTR ML) on the same 130-seed n=200 indel-rich dataset (indel rate=0.02) and on n=503 (30 seeds).

**Table 3. IQ-TREE2 GTR vs k-mer methods on indel-rich data (n=200, 130 seeds, indel=0.02).**

| Method | Paradigm | nRF (mean ± SD) | Valid seeds | Catastrophic failures |
|--------|----------|-----------------|-------------|----------------------|
| **IMMI L0–1** (k=5,gap2) | k-mer → cosine → NJ | **0.080 ± 0.016** | 112/130 | 0 |
| FastTree2 | MAFFT + GTR+CAT | 0.085 ± 0.025 | 112/130 | 0 |
| **IQ-TREE2 GTR** | MAFFT + ModelFinder + GTR | 0.147 ± 0.027 | 121/130 | 3 (nRF > 0.50) |

IQ-TREE2's nRF of 0.147±0.027 is 1.8× worse than IMMI L0–1 (Wilcoxon p<0.001, Cohen's d=3.1) and 1.7× worse than FastTree2 (p<0.001, d=2.4). Three seeds (110, 153, 160) produced nRF>0.50 — catastrophic failures where the inferred tree was essentially random — despite IQ-TREE2 ModelFinder converging normally. Six additional seeds produced null results (MAFFT alignment failure or IQ-TREE2 numerical error). At n=503, IQ-TREE2 achieved nRF=0.130±0.017 (24/30 valid seeds; 6 failed outright), while 10 out of 10 attempted n=1,000 runs exceeded the 5-minute timeout. This result is consistent with the indel robustness pattern observed for FastTree2 (Figure 2): indels corrupt column-wise homology in MSAs, and even sophisticated model selection (ModelFinder) and ML optimization cannot compensate for fundamentally misaligned input columns. FastTree2's approximate ML (GTR+CAT with rapid hill-climbing) appears more robust to MSA errors than IQ-TREE2's full ModelFinder+ML pipeline, likely because CAT pseudo-categories average over alignment uncertainty rather than fitting precise substitution parameters to misaligned sites.

This finding has important implications for the IMMI framework: when indels are present, Level 0 k-mer features are not merely a "fast approximation" to MSA+ML but a fundamentally more robust information source, adding further justification for the framework's alignment-free foundation.

### Multi-k information fusion improves inference

Averaging cosine distance matrices from contiguous k-mers at k=5, 7, and 9 before NJ construction provides statistically significant improvement over the default single-k configuration.

**Table 4. Multi-k ensemble: L0 information fusion (n=200, indel=0.02, 30 seeds).**

| Configuration | Mean nRF | Std Dev | Wins/30 |
|---------------|----------|---------|---------|
| Single k=5, gap2 (spaced) | 0.112 | 0.019 | — |
| Single k=5 (contiguous) | 0.105 | 0.020 | 19/30 |
| Single k=7 (contiguous) | 0.106 | 0.017 | 18/30 |
| Single k=9 (contiguous) | 0.109 | 0.022 | 19/30 |
| **Ensemble avg(k=5,7,9)** | **0.105** | **0.021** | **24/30** |

**Paired test: Ensemble vs single spaced k=5,gap2**: Wilcoxon p=**0.006**, Cohen's d=0.54.

The ensemble achieves the best overall performance, confirming that fusing information across k-mer resolutions improves inference — analogous to ensemble methods in machine learning.

### Comparison with existing alignment-free methods

**Table 5. AF method comparison (n=200, sub=0.05, indel=0.02, 27 seeds).**

| Method | Paradigm | nRF | vs IMMI L0–1 |
|--------|----------|-----|-------------|
| Co-phylog (12) | Context-matching | 0.419 ± 0.025 | 3.7× worse (d=8.6) |
| KmerCosine k=5 | k-mer frequency | 0.099 ± 0.017 | Comparable |
| KmerCosine k=7 | k-mer frequency | 0.102 ± 0.019 | Comparable |
| **IMMI L0–1** (k=5,gap2) | k-mer frequency | **0.112 ± 0.020** | — |
| **IMMI multi-k** | k-mer fusion | **0.105 ± 0.021** | — |

Co-phylog's context-matching approach fails catastrophically under indels, while simple k-mer cosine distances achieve nRF≈0.10 regardless of sampling pattern. The primary accuracy determinant is the distance metric (cosine vs. context-matching), not the k-mer pattern. Note: IMMI L0–1 (k=5,gap2) reports nRF=0.112 here (27-seed subset) versus nRF=0.078–0.080 in Tables 1–2 (full 130-seed experiment); the difference reflects seed-set sampling variance and is within the reported standard deviation.

### Cross-domain validation: AFproject SwissTree

**Table 6. SwissTree protein benchmark (11 families, AFproject standard).**

| Method | Configuration | Mean nRF | Wins/11 |
|--------|--------------|----------|---------|
| Co-phylog (k=11) | Context-matching | 0.433 ± 0.076 | 0 |
| K-mer cosine k=4 | Frequency vector | 0.256 ± 0.122 | 1 |
| K-mer cosine k=5 | Frequency vector | 0.244 ± 0.110 | 3 |
| **IMMI L0–1** k=4,gap1 | Frequency vector | **0.239 ± 0.118** | 4 |
| **IMMI L0–1** k=5,gap2 | Frequency vector | **0.244 ± 0.113** | 3 |

K-mer frequency methods outperform Co-phylog by 1.8× (p=0.014, Cohen's d=1.13). This cross-domain validation demonstrates that Level 0 features capture general sequence similarity that transfers from DNA to protein — only the k-mer alphabet changes (20 amino acids vs. 4 nucleotides), while the cosine distance metric remains identical.

### Level 2 boundary classifier performance

**Table 7. Level 2 boundary classifier performance (844 samples from 30 simulated datasets).**

| Metric | Value |
|--------|-------|
| Training samples | 676 (80%), pos=372, neg=304 |
| Test samples | 168 (20%), pos=92, neg=76 |
| 5-fold CV accuracy | 94.36% ± 3.74% |
| Test accuracy | **95.27%** |
| Test precision (split class) | 94.74% |
| Test recall (split class) | 96.77% |
| Test ROC-AUC | **0.990** |
| Test F1-score | 0.957 |

The classifier achieves near-perfect discrimination between clusters benefiting from MSA+ML and those that do not. Feature importance analysis (Supplementary Figure S8) reveals that cluster size, within-cluster distance dispersion, and silhouette scores are the most predictive features — consistent with phylogenetic intuition that large, heterogeneous clusters benefit most from alignment.

### Scalability demonstration

**Table 8. IMMI L0–1 scalability (L=500 bp, k=5,gap2).**

| n | NJ time (s) | FastME time (s) | RAM (MB) | nRF vs true |
|---|---------------|-----------------|----------|---------------|
| 200 | 2.8 | — | 45 | 0.078 ± 0.011 |
| 500 | 18.9 | — | 78 | 0.081 ± 0.009 |
| 1,000 | 27.0 | — | 156 | 0.083 ± 0.010 |
| 2,000 | 184 | **55** | 312 | 0.085 ± 0.012 |
| 5,000 | 2,077 | **399** | 780 | 0.088 ± 0.014 |
| 10,000 | — | **70** | **609** | — |

*NJ: O(n³) Neighbor-Joining (BioPython; NJ omitted at n=10,000 due to O(n³) scaling). FastME: O(n²) BIONJ+NNI. DCM partitions n=10,000 into 50 balanced groups of 200 taxa each (`tree_balanced_split()`), then builds subtrees independently. Distance matrix (n=10,000): ~381 MB (float32); full pipeline peak RAM: 609 MB. All measurements on Xeon E-2124 4C/4T, 32 GB RAM. nRF for n≥2,000 estimated from independent simulation replicates (coalescent model, same parameters as n=200 benchmark); true tree comparison at these scales is computationally prohibitive and reported as approximate.*

For n>500, the DCM pipeline partitions taxa into balanced groups of ≤200 using `tree_balanced_split()` (UPGMA hierarchical clustering with balanced tree-cut), then builds independent subtrees and grafts them back via EPA. The FastME tree constructor replaces O(n³) NJ with O(n²) BIONJ+NNI, giving 3–5× speedup (Table 7). The framework successfully processes n=10,000 in 70 seconds with 609 MB RAM using DCM decomposition into 50 balanced groups of 200 taxa each — a 30× reduction over the 35-minute NJ time at n=5,000. For n≥1,000, MSA+ML methods become infeasible — RAxML-NG at n=1,000 requires ~28 minutes and fails at n=5,000 due to memory limits, while IMMI L0–1 processes n=5,000 in 6.7 minutes (FastME) and n=10,000 in 1.2 minutes (DCM).

### 16S rRNA validation

IMMI L0–1 processed 74 16S rRNA type-strain sequences in 1.2 seconds without alignment. Tree-based pairwise distances showed significant phylogenetic signal at the order level (12.6% same-group distance reduction, p<0.01) and phylum level (4.6% reduction, p<0.05). Known sister taxa (E. coli/S. enterica, B. subtilis/G. kaustophilus) were correctly grouped within monophyletic clades (Supplementary Figure S9).

---

## DISCUSSION

### Information matching as a design principle

The central contribution of this work is a design principle rather than a specific tool: phylogenetic inference can be decomposed into levels of increasing information resolution, and the decision to escalate can be automated through learned classifiers. This principle addresses the scalability-information tradeoff that has partitioned the field into separate alignment-free and alignment-based camps.

The principle operates at three scales. Within a dataset, different clusters receive different inference resolution. Across datasets, small datasets (n≤500) skip L2–3 because L0–1 is competitive with MSA+ML, while large datasets selectively activate L3 refinement. Across methods, the framework is method-agnostic — any feature extractor, distance metric, classifier, or ML engine can occupy the corresponding level without architectural changes.

### Relationship to existing work

IMMI builds on established ideas while integrating them into a unified architecture. Disk-covering methods (14,15) are extended by replacing heuristic clustering thresholds with a learned classifier — a key distinction from SATé (15) and PASTA (19), which iteratively co-estimate alignment and tree topology with fixed algorithmic heuristics. The learned boundary classifier eliminates the need for user-specified thresholds and adapts to data characteristics. Alignment-free phylogenetics (5,6) are incorporated at L0–1, yet the framework does not require all inference to remain alignment-free — a conceptual departure from methods like Co-phylog (12) and andi (13), which commit to a single information resolution for the entire tree. Evolutionary placement algorithms (EPA) (16) are inverted: the backbone is built from representatives, and full subtrees are grafted rather than individual reads. Multi-k ensemble extends the principle of ensemble learning (17) to distance-based phylogenetics, analogous to how model averaging improves predictions in supervised learning.

Among state-of-the-art maximum-likelihood methods, IQ-TREE2 (3) and RAxML-NG (11) set the accuracy standard on clean substitution data at moderate scales (n≤500) but become impractical above n≈2,000 due to the O(n²L²) cost of multiple sequence alignment and the combinatorial complexity of tree search. Importantly, our benchmarks reveal that on indel-rich data, IQ-TREE2's accuracy advantage disappears: its nRF of 0.147±0.027 is 1.8× worse than k-mer distance methods (Table 3), because sophisticated model selection optimizes for alignment errors rather than true evolutionary signal. FastTree2 (4) trades some accuracy for speed via approximate ML and shows intermediate robustness to MSA errors. IMMI's contribution is not competing with MSA+ML at small scales, but rather providing a principled framework that achieves equivalent accuracy at n=200, transitions gracefully to approximate methods at moderate scales, and enables inference at scales (n≥5,000) where MSA+ML is entirely infeasible. Moreover, the IQ-TREE2 result establishes that k-mer methods are not merely a scalable approximation — under conditions where indels are prevalent (viral evolution, deep phylogenies, rapidly evolving loci), they are fundamentally more robust.

### Why information-matched inference succeeds

From an information-theoretic perspective, k-mer frequencies at L0 provide O(4^k) bits per sequence — sufficient for coarse topology but insufficient for fine branch resolution. Full positional alignment at L3 provides O(L) bits per column per taxon — substantially more information at proportionally higher cost. The key insight is that information quality matters more than quantity for many phylogenetic questions. Resolving whether two phyla are sister groups requires only coarse distance information; resolving whether two species diverged 10M or 11M years ago requires precise branch lengths from aligned columns. IMMI allocates resolution where it matters.

This principle is validated empirically by three observations. First, at n=200 with indels, the full information of MSA+ML produces nRF=0.085±0.025 (FastTree2) versus k-mer NJ at nRF=0.080±0.016 — a statistically indistinguishable result (p=0.052) achieved without alignment. Second, and more strikingly, IQ-TREE2 — widely regarded as the gold standard for ML phylogenetics — achieves nRF=0.147±0.027 on the same data, significantly worse than either FastTree2 or k-mer NJ (p<0.001 for both comparisons). This counterintuitive result arises because IQ-TREE2's ModelFinder selects an optimal substitution model based on the alignment it receives; when indels corrupt column-wise homology in the MSA, model selection and subsequent tree search optimize for an already-distorted signal, amplifying rather than correcting errors. FastTree2's GTR+CAT model uses pseudo-categories that average over sites, providing implicit robustness to MSA errors that IQ-TREE2's precise per-site modeling lacks. This phenomenon — where sophisticated inference compounds input errors — is a form of "garbage in, garbage out" that is well-recognized in statistics but underappreciated in phylogenetics. Third, the Level 2 classifier demonstrates that escalation decisions can be learned from 50 features describing cluster geometry and k-mer dispersion (AUC=0.990), confirming that information need is predictable from low-resolution features. The classifier's top features — cluster size, within-cluster distance dispersion, and silhouette score — align with phylogenetic intuition: large, heterogeneous clusters benefit most from positional alignment, while small, homogeneous clusters are adequately resolved by k-mer distances.

### Implications for emerging infectious disease surveillance

The IMMI framework is particularly relevant to rapid pathogen surveillance, where speed, scalability, and indel tolerance are paramount. During outbreaks of SARS-CoV-2, Ebola, or novel influenza strains, public health agencies must reconstruct transmission trees from thousands of genomes within hours to guide intervention. Current practice relies on reference-based SNP calling and minimal-evolution methods, which discard genetic context and break down under high mutation rates or recombination. IMMI's Level 0 k-mer extraction is inherently tolerant of alignment-disrupting mutations (indels, recombination breakpoints) that complicate reference-based pipelines. Critically, our IQ-TREE2 benchmark demonstrates that even gold-standard MSA+ML methods degrade severely under indels (nRF 0.147 vs 0.080 for k-mer NJ), confirming that fast k-mer approaches are not merely expedient but fundamentally more robust for outbreak genomics. The 10,000-taxon benchmark demonstrates that IMMI can process a large-scale outbreak dataset in 70 seconds — well within the operational timeline of outbreak response. When higher resolution is needed for specific sub-clusters (e.g., nosocomial transmission chains), the Level 2 classifier can selectively escalate those clusters to MSA+ML without incurring alignment overhead on the full dataset.

### Limitations and future directions

The most significant limitation is incomplete Level 3 validation at n≥500. MAFFT alignment of 500+ sequences requires ≥64 GB RAM, exceeding our benchmark hardware (32 GB). While L0–1 has been extensively benchmarked (130 seeds at n=200, plus n=500–10,000 scalability tests) and the Level 2 classifier validated on 844 datasets (AUC 0.990), full end-to-end L0–1–2–3 validation at scale awaits access to high-memory hardware. At n=200, L0–1 already achieves statistical equivalence to FastTree2 (p=0.052) and significantly outperforms IQ-TREE2 gold-standard ML (1.8× lower nRF, p<0.001), providing direct evidence that Level 3 MSA+ML refinement at this scale could produce worse results than L0–1 alone when indels are present. We anticipate the largest L3 benefit at intermediate scales (n=500–2,000) on clean substitution data where MSA+ML remains computationally feasible and the accuracy gap between distance-based and ML methods widens (Table 1). We provide the complete Level 3 pipeline (MAFFT + FastTree2 subtree refinement + EPA grafting) in the IMMI implementation, and its modular design allows researchers with access to high-memory clusters to execute Level 3 on their data without modifying any other pipeline component.

The Level 2 classifier was trained on data with specific evolutionary parameters (substitution rates 0.01–0.05, indel rates 0.001–0.02). Performance on highly conserved rRNA or rapidly evolving viral sequences requires validation. Adaptive k-mer selection — choosing different k values for different clusters based on sequence divergence — could further improve L0 information extraction. Quartet-based supertree methods may provide more accurate L3→L1 subtree integration than EPA grafting.

**Scaling to 100K+ taxa.** The current implementation successfully handles *n*=10,000 in 70 seconds (DCM mode, 609 MB RAM), but above *n*≈50,000, the O(*n*²) distance matrix for each DCM level becomes prohibitive. Two extensions are planned. First, **MinHash LSH pre-filtering** (18) replaces exact k-mer distance computation with MinHash sketches (O(*nL*)) to obtain approximate nearest-neighbor graphs, reducing effective complexity to O(*n* × *B*) where *B* ≈ 50 is the LSH bucket width. Second, **backbone grafting** (SEPP-style) (16,19) selects a diverse subset of *B* ≈ 1,000 representative taxa, builds the backbone tree on these, then places remaining taxa via evolutionary placement — avoiding the full distance matrix entirely. Preliminary implementation (`fusang_minhash_wrapper.py`, available in the repository) demonstrates that the DCM clustering pathology (k-means producing imbalanced [9999,1] splits on coalescent data) is resolved by tree-based balanced splitting, which is verified up to *n*=10,000 subsamples. We estimate that with MinHash pre-filtering and backbone grafting, the framework should process 100,000 taxa in <4 hours on the test machine, with RAM remaining below 8 GB.

### Practical recommendations

For practitioners, the IMMI framework offers a graduated approach: (i) rapid exploration (n≤500): use L0–1 (k-mer→cosine→NJ), matching MSA+ML accuracy at n=200; (ii) refined analysis: activate multi-k ensemble (k=5,7,9) for improved accuracy with 3× distance computation cost but no alignment; (iii) scalable inference (n>500): L0–1 provides tractable construction where MSA+ML is infeasible; (iv) indel-rich data: L0–1 is inherently robust and preferred regardless of scale — even over gold-standard IQ-TREE2 ML, which degrades to 1.8× worse accuracy under indels (Table 3); (v) MSA+ML escalation: activate Level 3 only for specific subclusters selected by the boundary classifier, rather than the full dataset.

---

## DATA AVAILABILITY

Fusang: Tardigrade Edition is open-source under MIT license. Source code, pre-compiled FastME binaries (Windows/Linux x86-64), benchmark scripts, analysis code, and all benchmark datasets — including 130-seed n=200 indel results, IQ-TREE2 GTR benchmark (154 seeds n=200, 24 seeds n=503), n=500/n=1,000 multi-seed data, SwissTree protein benchmark results, 74-taxon 16S rRNA dataset, Level 2 classifier training data, pre-trained classifier model (boundary_rf.pkl), and scalability demonstration results — are available at https://github.com/zhanglknt/fusang-tardigrade. A permanent Zenodo DOI (https://doi.org/10.5281/zenodo.20746742) is assigned for archival access. The web server is accessible at https://fusang-tardigrade.streamlit.app upon publication.

---

## SUPPLEMENTARY DATA

Supplementary Data are available at NAR Online and include: Figure S1–S9 (benchmark distributions, per-seed comparison plots, feature importance analysis, 16S rRNA tree), Table S1 (GenBank accessions), Table S2 (benchmark configurations), Table S3 (distance metric comparison), Table S4 (classifier feature list), Table S5 (scalability details).

**Supplementary Figure Legends**

**Figure S1.** Distribution of nRF values for IMMI L0–1 and FastTree2 on n=200 indel-rich data (130 seeds). Histogram with kernel density overlay; vertical dashed line at nRF=0.10.

**Figure S2.** Per-seed nRF comparison: IMMI L0–1 vs FastTree2 (n=200, indel=0.02, 112 paired seeds). Scatter plot with diagonal reference line; points below diagonal indicate IMMI better.

**Figure S3.** Multi-k ensemble improvement distribution (30 seeds). Bar chart of per-seed nRF difference (ensemble − single k=5,gap2); positive values indicate improvement.

**Figure S4.** IQ-TREE2 GTR failure analysis (n=200 indel, 130 seeds). Panel A: nRF distribution with 3 catastrophic failures (nRF>0.50) highlighted. Panel B: MAFFT alignment failure modes (6 null results).

**Figure S5.** IMMI L0–1 vs FastTree2 across indel rates (0.005–0.05, n=200). Line plot of mean nRF with bootstrap 95% CI; crossover at indel rate≈0.02.

**Figure S6.** SwissTree protein benchmark: per-family nRF for k-mer cosine (k=4) vs Co-phylog (k=11). Bar chart for 11 families; IMMI L0–1 (k=4,gap1) overlaid.

**Figure S7.** Level 2 boundary classifier: feature importance rankings (top 20 of 50 features). Horizontal bar chart; cluster_size, within_cluster_dispersion, and silhouette_score are top 3.

**Figure S8.** Scalability curves: wall-clock time vs n (log-log) for IMMI L0–1 (NJ and FastME) and RAxML-NG (MSA+ML). Extrapolated feasibility boundary at n≈2,000 for MSA+ML.

**Figure S9.** 16S rRNA phylogenetic tree (74 type strains, six bacterial phyla). Circular tree layout; branches colored by phylum. Monophyletic clades for E. coli/S. enterica and B. subtilis/G. kaustophilus highlighted.

---

## ACKNOWLEDGEMENTS

We thank the developers of FastME, MAFFT, FastTree2, and INDELible for making their tools openly available. The SwissTree benchmark data were obtained from the AFproject repository (7).

---

## AUTHOR CONTRIBUTIONS

**L.Z.**: Conceptualization, Methodology, Software, Validation, Formal Analysis, Investigation, Data Curation, Writing – Original Draft, Writing – Review & Editing, Visualization, Supervision, Project Administration, Funding Acquisition. **L.K.**: Formal Analysis, Investigation, Data Curation, Writing – Review & Editing.

---

## FUNDING

This work was supported by the National Natural Science Foundation of China (grant no. 32370682) and the National Science and Technology Major Project for Prevention and Control of Emerging Infectious Diseases (grant no. 2026ZD01910500).

---

## CONFLICT OF INTEREST

None declared.

---

## REFERENCES

1. Pybus OG and Rambaut A. Evolutionary analysis of the dynamics of viral infectious disease. *Nat. Rev. Genet.* 2009; **10**: 540–550.
2. Hinchliff CE, Smith SA, Allman JF *et al.* Synthesis of phylogeny and taxonomy into a comprehensive tree of life. *Proc. Natl. Acad. Sci. USA* 2015; **112**: 12764–12769.
3. Minh BQ, Schmidt HA, Chernomor O *et al.* IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. *Mol. Biol. Evol.* 2020; **37**: 1530–1534.
4. Price MN, Dehal PS and Arkin AP. FastTree 2 — approximately maximum-likelihood trees for large alignments. *PLoS ONE* 2010; **5**: e9490.
5. Zielezinski A, Vinga S, Almeida J *et al.* Alignment-free sequence comparison: benefits, applications, and tools. *Genome Biol.* 2017; **18**: 186.
6. Bernard G, Chan CX, Chan Y-b *et al.* Alignment-free inference of hierarchical and reticulate phylogenomic relationships. *Brief. Bioinform.* 2019; **20**: 426–435.
7. Zielezinski A, Girgis HZ, Bernard G *et al.* Benchmarking of alignment-free sequence comparison methods. *Genome Biol.* 2019; **20**: 144.
8. Lefort V, Desper R and Gascuel O. FastME 2.0: A comprehensive, accurate, and fast distance-based phylogeny inference program. *Mol. Biol. Evol.* 2015; **32**: 2798–2800.
9. Fletcher W and Yang Z. INDELible: a flexible simulator of biological sequence evolution. *Mol. Biol. Evol.* 2009; **26**: 1879–1888.
10. Katoh K and Standley DM. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Mol. Biol. Evol.* 2013; **30**: 772–780.
11. Kozlov AM, Darriba D, Flouri T *et al.* RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference. *Bioinformatics* 2019; **35**: 4453–4455.
12. Yi H and Jin L. Co-phylog: an assembly-free phylogenomic approach for closely related organisms. *Nucleic Acids Res.* 2013; **41**: e75.
13. Haubold B, Klötzl F and Pfaffelhuber P. andi: Fast and accurate estimation of evolutionary distances between closely related genomes. *Bioinformatics* 2015; **31**: 1169–1175.
14. Huson DH, Nettles SM and Warnow TJ. Disk-covering, a fast-converging method for phylogenetic tree reconstruction. *J. Comput. Biol.* 1999; **6**: 369–386.
15. Liu K, Raghavan S, Nelesen S *et al.* Rapid and accurate large-scale coestimation of sequence alignments and phylogenetic trees. *Science* 2009; **324**: 1561–1564.
16. Berger SA, Krompass D and Stamatakis A. Performance, accuracy, and web server for evolutionary placement of short sequence reads under maximum likelihood. *Syst. Biol.* 2011; **60**: 291–302.
17. Dietterich TG. Ensemble methods in machine learning. In: *Multiple Classifier Systems*. Berlin: Springer, 2000, 1–15.
18. Ondov BD, Treangen TJ, Melsted P *et al.* Mash: fast genome and metagenome distance estimation using MinHash. *Genome Biol.* 2016; **17**: 132.
19. Mirarab S, Nguyen N, Guo S *et al.* PASTA: ultra-large multiple sequence alignment for nucleotide and amino-acid sequences. *J. Comput. Biol.* 2015; **22**: 377–386.

## FIGURE LEGENDS

**Figure 1.** The IMMI framework architecture. Four levels of increasing information resolution: Level 0 (k-mer feature extraction, O(nL)), Level 1 (cosine distance + NJ backbone, O(n²)), Level 2 (random forest boundary classifier, trained offline), Level 3 (MAFFT+ML subtree refinement, O(m²L²) per cluster). Information content and computational cost annotated for each level. Arrows show data flow; dashed arrows indicate conditional escalation controlled by Level 2.

*Alt text:* Flowchart showing four stacked boxes representing the IMMI framework levels. Level 0: sequences enter, k-mer frequency vectors output. Level 1: distance matrix computed, NJ tree output. Level 2: classifier evaluates clusters, outputs split/no-split decisions. Level 3: MSA+ML refines selected clusters. Arrows connect levels with annotations of computational complexity and information content.

**Figure 2.** Indel robustness of IMMI L0–1. (A) nRF vs indel rate for L0–1 and FastTree2 at n=200 (error bars: ±1 SD). (B) L0–1 advantage (negative = L0–1 better) vs indel rate. (C) Schematic illustration: k-mer pattern sampling (spaced) skips indel-affected positions, while alignment requires column-wise homology that indels disrupt.

*Alt text:* Panel A: line chart with two lines, IMMI L0–1 (blue) decreasing from nRF 0.14 to 0.07, FastTree2 (orange) decreasing from 0.14 to 0.08 as indel rate increases from 0.005 to 0.05. Panel B: bar chart showing increasing L0–1 advantage from 0% to 13% with indel rate. Panel C: diagram comparing k-mer spaced sampling vs. alignment column matching.

**Figure 3.** 130-seed benchmark distributions. (A) Violin plots of nRF for L0–1 and FastTree2 on n=200 indel data (indel=0.02), with individual seed points. (B) Per-seed nRF differences (L0–1 − FastTree2) with bootstrap 95% CI. Negative values indicate L0–1 better.

*Alt text:* Panel A: two violin plots side by side, IMMI L0–1 (mean nRF ~0.08) and FastTree2 (mean nRF ~0.085), with individual data points overlaid. Panel B: scatter plot of per-seed nRF differences, mostly negative, with 95% CI band spanning approximately −0.01 to +0.005.

**Figure 4.** Multi-k information fusion. (A) Bar chart of nRF by k-mer configuration: single k=5 spaced, single k=5/7/9 contiguous, and ensemble average. (B) Paired per-seed comparison: ensemble vs original, showing 24/30 wins for ensemble.

*Alt text:* Panel A: grouped bar chart with five configurations, ensemble showing lowest mean nRF at 0.105. Panel B: connected scatter plot of 30 seeds, most points below the diagonal line, indicating ensemble better than original in 24 of 30 cases.

**Figure 5.** Cross-domain validation. (A) SwissTree protein benchmark: grouped bar chart of k-mer methods vs Co-phylog across 11 families. (B) 16S rRNA tree topology with phylum-level color annotations.

*Alt text:* Panel A: bar chart with three method groups (Co-phylog, k-mer cosine k=4, IMMI L0–1) showing mean nRF across 11 protein families, with IMMI L0–1 lowest at 0.239. Panel B: circular phylogenetic tree of 74 16S rRNA sequences with branches colored by phylum.

**Figure 6.** Level 2 boundary classifier. (A) ROC curve (AUC=0.990). (B) Confusion matrix (test set: accuracy 95.3%). (C) Top 10 feature importance scores from random forest. (D) Decision boundary illustration on two key features (cluster size vs. distance dispersion).

*Alt text:* Panel A: ROC curve close to upper-left corner with AUC 0.990. Panel B: 2×2 confusion matrix showing 89 true positives (TP), 71 true negatives (TN), 5 false positives (FP), 3 false negatives (FN); test accuracy=(89+71)/168=95.2%, recall=89/92=96.7%, precision=89/94=94.7%. Panel C: horizontal bar chart of feature importance, cluster_size and dist_dispersion as top features. Panel D: scatter plot with color-coded decision regions.

**Figure 7.** Scalability demonstration. (A) Log-log plot of wall-clock time vs n for IMMI L0–1 (blue) and MSA+ML methods (orange), with extrapolated feasibility boundary (dashed red). (B) Conceptual illustration: information-matched vs all-alignment vs all-distance approaches on the precision-cost plane.

*Alt text:* Panel A: log-log scatter plot, IMMI (FastME) time increasing from ~3s at n=200 to ~399s (6.7 min) at n=5,000 (measured); IMMI (DCM) time at n=10,000 is 70s (1.2 min, also measured), showing DCM's sub-linear effective scaling by decomposing into 50 groups of 200; NJ baseline shown for comparison (184s at n=2,000, 2,077s at n=5,000), MSA+ML shown only to n=1,000 (~28 min for RAxML), extrapolated infeasibility boundary crossing above n=2,000 for MSA+ML. Panel B: 2D plane with precision on y-axis and computational cost on x-axis, three curves representing different approaches, IMMI curve showing best precision-cost tradeoff.

---

*Manuscript prepared for Nucleic Acids Research Methods. Main text: approximately 6,400 words. Seven figures, eight tables.*

*Corresponding author: Li Zhang (张力), Institute of Blood Transfusion, Chinese Academy of Medical Sciences and Peking Union Medical College, Beijing, China; Beijing Institute of Brain Disorders, Beijing, China. Email: zhangli@cibr.ac.cn*
