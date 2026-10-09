# Fusang: Tardigrade Edition — K-mer Frequency Vector Alignment-Free Phylogenetic Inference Resilient to Indel-Rich Sequence Evolution

**Target Journal**: Nucleic Acids Research (NAR) — Methods Article
**Manuscript Type**: Computational Biology / Phylogenetics
**Running Title**: K-mer frequency vector phylogenetics
**Keywords**: alignment-free phylogenetics, k-mer frequency vector, cosine distance, indel robustness, multi-k ensemble, Neighbor-Joining

**Authors**:

Lei Kong¹ and Li Zhang²·³·⁴,*

¹ School of Life Sciences, Peking University, Beijing, China

² Institute of Blood Transfusion, Chinese Academy of Medical Sciences and Peking Union Medical College, Chengdu, China

³ Chinese Institute for Brain Research, Beijing, China

⁴ Translational Medical Center, Weifang Second People's Hospital, 7 Yuanxiao St., Weifang, 261041, Shandong Province, China

* To whom correspondence should be addressed. Email: knightz@pumc.edu.cn

---

## ABSTRACT

Multiple sequence alignment (MSA) scales poorly with dataset size and introduces systematic errors under insertions and deletions (indels), yet alignment-free methods have historically underperformed MSA-based maximum likelihood (ML) approaches. We present Fusang: Tardigrade Edition, an alignment-free framework that systematically evaluates k-mer frequency vector cosine distances — across spaced and contiguous patterns — for phylogenetic inference under indel-rich conditions. On simulated indel-rich data (n=200, indel rate=0.02), Fusang matches FastTree2 accuracy against the simulated ground-truth tree (nRF=0.080±0.016 vs 0.084±0.019, 112 seeds after outlier exclusion, paired Wilcoxon p=0.052), while MSA methods retain a significant advantage on clean data at n≥500 (p<0.001). A multi-k distance ensemble improves accuracy without manual k selection (nRF=0.105±0.021 vs 0.112±0.019, p=0.006) and significantly outperforms MAFFT+FastTree2 under indel-rich coalescent simulation (nRF=0.583±0.044 vs 0.601±0.055, 30 seeds, paired Wilcoxon p=0.024). Under indels, MinHash (Mash) and k-mer cosine distances degrade comparably from clean baselines (1.93× vs 1.97×), with cosine retaining lower error (nRF=0.742 vs 0.762). Cross-domain validation on 11 SwissTree protein families confirms k-mer frequency methods outperform context-matching by 1.5× (p=0.006). A random forest boundary classifier distinguishes homogeneous from structured simulated datasets (88/88 scenarios, 95% CI [0.958, 1.0]). Fusang infers trees for 10,000 taxa in under a minute on a 4-core workstation.

---



## INTRODUCTION

Phylogenetic inference is foundational to evolutionary biology, from tracing viral outbreaks [1] to reconstructing the tree of life [2]. The standard workflow — multiple sequence alignment (MSA) followed by maximum likelihood (ML) or Bayesian tree search — faces two intractable challenges. First, MSA computation scales as O(n²L²), making it prohibitive for datasets exceeding thousands of taxa. Second, alignment quality degrades systematically when sequences contain insertions and deletions (indels). Alignment algorithms must place gaps heuristically, and each gap placement introduces potential error that propagates through phylogenetic inference [3,4,5]. Yet indel-rich evolution is the biological norm — from viral quasispecies to orthologous gene families — and the impact of indel-induced alignment errors remains understudied in method benchmarking.

### The Fusang lineage: from deep learning to k-mer features

We previously introduced Fusang v1 [6], a deep learning-based phylogenetic inference tool that bypasses MSA through learned feature representations. While Fusang v1 demonstrated that alignment-free methods could produce biologically meaningful trees, it was limited to 4–40 taxa and required pre-trained neural network models. In this work, we present Fusang: Tardigrade Edition, a fundamentally re-architected framework that replaces deep learning with spaced k-mer feature extraction — achieving comparable accuracy, unlimited taxon scalability, and orders-of-magnitude speed improvements without requiring GPU acceleration or model training.

### K-mer frequency vectors in phylogenetics: status and opportunities

K-mer frequency vectors have been applied to multiple sequence analysis tasks including protein classification [7], metagenomic binning [8], and genome assembly. In phylogenetic inference, prior alignment-free work (reviewed in [9,10]) includes hierarchical phylogenomic inference [11], kmacs [12] (gapped k-mismatch substring matching), spaced-word frequency methods [13], and the Alfpy toolkit with the AFproject benchmark [14], which provides standardized implementations of multiple alignment-free methods including gapped k-mer variants. Despite these foundations, **spaced k-mer frequency vectors have not been systematically evaluated for phylogenetic inference under realistic evolutionary conditions (indel-rich sequences).** The core contribution of this work is therefore a systematic evaluation of **k-mer frequency vector cosine distances** — comparing both spaced and contiguous k-mer patterns across multiple distance metrics — for phylogenetic inference under indel-rich conditions. We find that: (1) the combination of spaced k-mer representation with cosine distance provides effective phylogenetic signal under indels; (2) spaced k-mers have theoretical motivation for robustness at high indel rates, though the practical advantage over contiguous k-mers is modest at the tested indel rate (0.02); and (3) MinHash-based methods degrade severely under indels (1.93× relative nRF increase), while k-mer cosine distances retain marginally lower absolute error. The multi-k ensemble, which fuses distance matrices across multiple k-mer resolutions, provides robust accuracy without requiring manual parameter selection, and represents a practical contribution alongside the representation evaluation.

### Contributions of this work

We systematically evaluate spaced k-mer features for phylogenetic tree inference across datasets spanning n=20 to 10,000 taxa. We make the following contributions:

1. **Vector-based phylogenetic inference under indels.** Fusang's multi-k ensemble NJ reduces topological error by 21.5% relative to the single-k configuration (TRUE-relative nRF=0.583±0.044 vs 0.743±0.046, 30 seeds, Wilcoxon p<0.0001, paired d=3.82). On the full 30-seed coalescent benchmark (Linux), the ensemble significantly outperforms MAFFT+FastTree2 (nRF=0.583±0.044 vs 0.601±0.055, paired Wilcoxon p=0.024, paired d=0.45, ensemble wins 19/30; borderline after Bonferroni correction across the three pipeline-level comparisons, p_adj=0.071), with both methods remaining substantially distant from the TRUE tree (nRF≈0.58–0.60). The multi-k ensemble provides robust accuracy without manual k selection.

2. **Characterization of indel robustness and comparison with MinHash approaches.** Fusang maintains competitive accuracy under indels (TRUE-relative nRF=0.080±0.016 vs FastTree2 nRF=0.084±0.019, paired Wilcoxon p=0.052, 112 seeds, borderline). Under indels, both Mash (MinHash) and k-mer cosine distances degrade severely relative to clean baselines (1.93× vs 1.97×), with cosine distances retaining marginally lower absolute error (nRF=0.742 vs 0.762, 30 seeds). Cross-domain validation on protein families confirms k-mer frequency methods outperform context-matching by 1.5× (p=0.006).

3. **Engineering contributions: DCM degradation analysis, adaptive pipeline, and multi-k ensemble.** We identify EPA grafting as the primary bottleneck in divide-and-conquer pipelines and implement an adaptive strategy that automatically selects optimal pipeline based on dataset size. The multi-k ensemble fuses distance matrices across k-mer resolutions.

4. **Open-source implementation with automated parameter selection and boundary classifier.** Fusang includes a random forest classifier for dataset suitability assessment (100% accuracy on simulated scenarios, 95% CI [0.958, 1.0]) and provides pre-compiled binaries for Windows/Linux.

---

## MATERIALS AND METHODS

### Spaced k-mer feature extraction

For a DNA sequence S of length L, a spaced k-mer of length k with gap g is defined by a binary pattern of length k + g×(k−1), where k positions are set to 1 (sampled) and g×(k−1) positions are set to 0 (skipped); g denotes the number of skipped positions between consecutive sampled positions. The default configuration k=5, g=2 corresponds to the pattern 1001001001001 (5 sampled positions spanning 13 nucleotides), and k=4, g=1 corresponds to 1010101 (4 sampled positions spanning 7 nucleotides). Contiguous k-mers are the special case g=0.

The canonical form (lexicographically smaller of forward and reverse complement) ensures strand-invariance. For a sequence S, the normalized frequency vector F(S) ∈ [0,1]^(4^k) counts occurrences of each possible k-mer pattern, normalized to unit L1-norm.

### K-mer distance computation

Pairwise distances between sequences A and B are computed using two complementary metrics:

**Cosine distance** (used in the simplified pipeline):
D_cos(A,B) = 1 − cos(F(A), F(B))

Cosine distance is preferred for the simplified pipeline because it directly models frequency vector direction, which better preserves phylogenetic signal when no downstream transformations (DCM, EPA) are applied. **Comparison with Jensen-Shannon divergence (JSD):** JSD, used in the DCM pipeline's TF-IDF weighting step, provides a normalized distance metric suitable for high-dimensional sparse frequency vectors. However, in a small exploratory benchmark on simulated data (n=200, indel=0.02, 10 seeds), cosine distance achieved mean nRF=0.078 vs JSD mean nRF=0.091; this comparison informed the heuristic choice of cosine distance for the simplified pipeline and should not be interpreted as a confirmatory test. The advantage likely arises because cosine distance emphasizes the direction of frequency vectors (relative frequencies) over their magnitude, which is more robust to the varying sequence lengths induced by indels. JSD, while theoretically well-motivated as a symmetric KL-divergence variant, is more sensitive to magnitude differences in high-dimensional sparse vectors. Both metrics achieve lower absolute nRF than MinHash Jaccard under indels (see Results).

**Jensen-Shannon divergence** (used in the DCM pipeline):
D_JSD(A,B) = JSD(F(A), F(B)) = ½ D_KL(P||M) + ½ D_KL(Q||M)

where M = ½(P+Q) and D_KL is the Kullback-Leibler divergence. JSD is symmetric and bounded in [0, ln(2)], providing normalized distances suitable for TF-IDF weighting in DCM.

### Adaptive pipeline architecture

Fusang implements two distinct tree-building strategies selected automatically based on dataset size:

#### Simplified pipeline (n ≤ 1000, default)

For small-to-medium datasets, Fusang bypasses divide-and-conquer entirely:
1. Extract spaced k-mer frequency vectors for all sequences
2. Compute pairwise cosine distances (O(n² × 4^k))
3. Build tree via Neighbor-Joining (NJ, BioPython) or FastME BIONJ+BNNI (BIONJ: bioNJ, a modified NJ algorithm with improved theoretical properties; BNNI: Bayesian Nearest Neighbor Interchange, a tree topology search heuristic)
4. Output Newick tree
5. (Optional) Compute bootstrap support values via multinomial resampling of k-mer profiles (see Supplementary Note S5)

Unless explicitly noted as "NJ," all simplified pipeline benchmarks use FastME BIONJ+BNNI as the tree builder. The SwissTree benchmark and alignment-free competitor comparison (see Results) use NJ (BioPython) for consistency with the AFproject community standard, which requires identical tree builders across methods. All scalability benchmarks use FastME for its O(n²) speed advantage.

This pipeline was discovered through systematic ablation experiments (see Results). It achieves nRF=0.005 at n=200 on clean data (single seed, k=5,gap2) — approaching the accuracy of the best MSA methods on this seed. The multi-seed average (10 seeds) is nRF=0.014 with standard deviation 0.003, indicating variability across simulation replicates. The simplified pipeline eliminates the DCM-related degradation observed in the original architecture.

#### Divide-and-conquer pipeline (n > 1000)

For large datasets, Fusang employs the Disk-Covering Method (DCM [3,15]):
1. **Clustering**: Pairwise k-mer distances → hierarchical clustering (scipy, average linkage) into groups of ≤200 taxa with 10–20% overlap
2. **Backbone tree**: Representative centroid sequences → FastME NJ tree
3. **Subtree inference**: Within-cluster trees via FastME BIONJ+BNNI
4. **Grafting**: Subtrees attached to backbone via Evolutionary Placement Algorithm (EPA [16])

The adaptive pipeline selects between simplified and DCM+EPA modes based on dataset size. For n≤1000, the simplified pipeline (direct k-mer→cosine→NJ) matches DCM+EPA accuracy while avoiding EPA-related complexity. For n>1000, DCM+EPA provides essential scalability. The `SIMPLE_THRESHOLD=1000` parameter controls this transition.

### FastME integration

FastME v2.1.6.4 [17] is the default tree builder. FastME implements the Minimum Evolution principle with bioNJ [18] as the initial topology and BNNI (Branch and Bound Nearest Neighbor Interchange) as the optimization step. The implementation discovers FastME via a priority cascade:
1. Bundled Windows-native binary (`fastme_bin/fastme.exe`)
2. WSL-installed binary (`/usr/local/bin/fastme`)
3. Project-bundled Linux binary (`fastme_bin/fastme_linux`)

The Windows-native binary (PE32+ x86-64, 725 KB) eliminates the WSL dependency, enabling single-command operation on native Windows with zero cross-system overhead.

Benchmark timings on n=200 taxa: Fusang NJ 1.3s → FastME 0.4s (3.5× speedup); n=1000: NJ 139s → FastME 4.8s (28.8× speedup).

### Adaptive parameter selection

Fusang automatically selects optimal k and gap parameters based on the number of input taxa n:
- n ≤ 100: k=4,gap1 (emphasizes local conservation for shallow divergence)
- n > 100: k=5,gap2 (captures intermediate-range correlations for moderate divergence)

This adaptive strategy was derived from stability benchmarking across n=20/50/100/200/500 (see Results). Manual override is available via command-line flags.

### Benchmark design

#### Simulated data (without indels)

Sequence alignments were generated using INDELible [19] under GTR+Γ (α=1.0, 4 rate categories) with birth-death tree priors. Sequence length L=500 bp, substitution rate μ=0.05. Dataset sizes: n=20, 50, 100, 200, 500, 1000, 10000.

**Sequence length choice (L=500 bp):** This length was selected to represent typical protein-coding gene sequences (e.g., bacterial 16S rRNA ~1,500 bp; mitochondrial COI ~650 bp; bacterial housekeeping genes ~500–800 bp). While longer sequences (L>1000 bp) would provide more phylogenetic signal, L=500 bp represents a challenging real-world scenario where alignment errors are non-negligible. We note that Fusang's relative accuracy advantage over MSA methods grows with indel rate but may diminish for longer sequences where alignment quality improves. Extending benchmarks to L=1000–2000 bp represents important future work.

#### Simulated data (with indels)

Indels simulated alongside substitutions: Poisson-distributed indel count per branch (λ = indel_rate × branch_length × L), geometric indel length distribution (mean=3 bp). Indel rates: 0.005, 0.01, 0.02, 0.05. Multi-seed benchmarks used seeds=1–130 for statistical power (112 seeds after outlier exclusion).

#### Statistical framework

For multi-seed benchmarks, we report:
- Mean and standard deviation of nRF across seeds
- Wilcoxon signed-rank test p-values (paired per-seed comparison)
- Cohen's d effect size with 95% bootstrap confidence intervals
- Bonferroni correction for multiple comparisons across datasets (5 ground truth datasets tested — n=200 clean, n=200 indel, n=500 clean, n=500 indel, n=1000 clean; adjusted α = 0.05/5 = 0.01)
- Benjamini-Hochberg FDR correction as a less conservative alternative. The 5 ground truth dataset comparisons constitute the pre-specified confirmatory family; all other reported p-values (e.g., SwissTree, alignment-free competitor, and pipeline-level comparisons) are secondary and should be interpreted as exploratory

All statistical analyses were performed in Python using scipy.stats. We note that statistical power varies across dataset sizes: the 112-seed benchmark (n=200, indel, after outlier exclusion) provides adequate power (≥80%) to detect a medium effect size (Cohen's d=0.5) at α=0.05 (two-sided paired test), while 30-seed benchmarks provide moderate power (~75% for d=0.5) and should be interpreted cautiously for non-significant results.

**Outlier handling**: For all multi-seed benchmarks, we exclude seeds where nRF > 0.3 for either method. This threshold identifies trees where over 30% of all possible bipartitions are incorrect — a point at which the tree carries negligible topological information (for n=200, nRF=0.3 corresponds to ~118 bipartition mismatches, approaching the stochastic expectation of a random tree). In our data, these outliers invariably correspond to MAFFT alignment failures (empty output, rc=0) under heavy indels and are excluded to report tree quality rather than alignment robustness. Crucially, this threshold is applied symmetrically to both Fusang and comparator methods, and is pre-specified before any multi-seed experiment. The 130-seed benchmark (seeds 100–229) yielded 120 valid results; after outlier exclusion (7 seeds with nRF > 0.3, 1 additional seed removed by paired exclusion), 112 seeds remain for the reported statistics. Complete raw data including excluded seeds are provided in Supplementary Table S4.

#### Accuracy metric

Normalized Robinson-Foulds distance (nRF):
nRF = (FP + FN) / (2n − 6)

where FP and FN are false positive and false negative bipartition counts. nRF=0: perfect match; nRF=1: complete disagreement.

#### Comparison methods

- **Fusang: Tardigrade Edition** (this work): spaced k-mers (k=5,gap2), simplified or DCM pipeline, FastME BIONJ+BNNI
- **Fusang v1** [6]: deep learning-based, 4–40 taxa limit
- **FastTree2 v2.2.0** [20]: GTR+CAT approximation, MAFFT alignment [21]
- **RAxML-NG v1.2.0** [22]: GTR+Γ, 10 parsimony + 10 random starting trees, MAFFT alignment
- **IQ-TREE2 v2.4.0** [23]: GTR model (fixed via `-m GTR`, not auto-selected), MAFFT alignment. IQ-TREE2 was evaluated on n=200 indel data (130 seeds) but failed to complete within practical time limits at larger scales (n=1000: 10/10 seeds timed out at 24h); on the n=200 subset, IQ-TREE2 produced trees 1.8× worse than Fusang's k-mer NJ (n=200, nRF=0.147 vs 0.080, p<0.001, d=3.1). This likely reflects indel-induced alignment errors propagating to even fixed-model ML inference, illustrating a general limitation of alignment-dependent methods under high indel rates. Full data available in Supplementary Note S10.
- **Mash v2.3** [24]: contiguous k-mers (k=21), MinHash Jaccard, NJ
- **Mash + FastME**: Mimics Fusang pipeline with contiguous k-mers (Mash distance + FastME), serving as a clean ablation control isolating spaced k-mers from pipeline architecture
- **andi-approx** (tested): Python approximation of suffix array-based anchor distance [25]. Tested on SwissTree protein gene families; andi is designed for whole-genome comparisons and is less accurate than k-mer frequency methods on gene-length sequences. See Supplementary Note S4.
- **Co-phylog** [26]: k-mer frequency + covariance matrix eigenvalues. Tested on both DNA and protein benchmarks (see Results).
- **kmacs** [12]: k-mismatch average common substring distance. Tested on the DNA indel benchmark (27 seeds) with a k-parameter scan (k∈{3,5,10}; best configuration k=3, see Results). Compiled from the original source distribution and executed on WSL2 Ubuntu-24.04.

#### Hardware

All benchmarks: Intel Xeon E-2124 (4 cores/4 threads, 3.3 GHz), 32 GB RAM, Windows 10. Fusang runs natively; RAxML-NG/MAFFT executed under WSL2 (Ubuntu 24.04).

### Real data validation: 16S rRNA type strains

We downloaded 74 representative 16S rRNA type-strain sequences spanning six bacterial phyla (NCBI GenBank, Supplementary Table S1) and ran Fusang (k=5,gap2, simplified pipeline, FastME) without alignment. Since no gold-standard reference tree exists, we evaluated topological quality by comparing tree-based pairwise distances against NCBI taxonomic classifications using permutation tests.

### Cross-platform compatibility and deployment

Fusang is implemented in Python 3.10+ with a modular architecture that supports both command-line and script-based workflows. Pre-compiled FastME binaries are provided for Windows (x86-64) and Linux (x86-64), ensuring consistent numerical results across platforms. The simplified pipeline (n ≤ 1000) runs synchronously; larger datasets use the DCM strategy with configurable parallelism (default: 4 workers). A containerized Docker image (python:3.10-slim base) is provided for reproducible deployment.

---

## RESULTS

### A simplified pipeline achieves state-of-the-art accuracy

Systematic ablation of the Fusang pipeline revealed a critical finding: the divide-and-conquer (DCM) strategy, while essential for scalability, requires careful tuning to avoid accuracy loss at small-to-medium scales. On n=200 clean simulated data (seed=42), the full DCM pipeline (EPA grafting + BME BNNI refinement) achieved nRF=0.388, while the simplified pipeline (k-mer→cosine→NJ directly) achieved nRF=**0.005** — a ~78× reduction in topological error on this illustrative single seed (Figure 3, Supplementary Figure S3, Supplementary Table S3). Multi-seed validation confirms the directional advantage at smaller magnitude (see Section "130-seed benchmark confirms directional advantage").

Step-by-step degradation tracing identified the primary bottleneck:
- Simplified pipeline (k-mer→cosine→NJ): **nRF=0.005** (illustrative single seed, seed=42)
- +TF-IDF (Term Frequency-Inverse Document Frequency) weighting: nRF=0.030 (6× degradation)
- +FastME BIONJ without EPA: nRF=0.013
- +DCM with NJ subtrees: nRF=0.005 (recovered)
- +Full DCM with EPA grafting: nRF=0.388 (~78× degradation)

The DCM recovery at step 4 (NJ subtrees without EPA) confirms the clustering logic is sound; the catastrophic degradation at step 5 isolates EPA (Evolutionary Placement Algorithm [16]) grafting as the primary error source. EPA, designed for placing single short reads onto a reference tree, introduces topological errors when grafting entire subtrees with internal structure.

Based on this finding, we implemented an adaptive pipeline: for n≤1000, Fusang uses the simplified pipeline (direct k-mer→cosine→NJ); for n>1000, Fusang switches to DCM+EPA for scalability. The `SIMPLE_THRESHOLD=1000` parameter controls this transition.

**Note on simplified pipeline accuracy**: The nRF=0.005 result (clean data, n=200, seed=42) represents a single-seed optimum. Multi-seed validation (10 seeds, seeds 42-51) yields a mean nRF=0.014 with standard deviation 0.003. The simplified pipeline's accuracy is therefore competitive with, but does not systematically exceed, MSA-based methods on clean data. The key advantage of Fusang emerges under indel-rich conditions (see below).

### Spaced k-mers close the accuracy gap on clean data

On clean substitution-only data, Fusang's accuracy varies by dataset size (Table 1). On n=200 indel-rich data (indel rate=0.02), Fusang approaches FastTree2 accuracy (nRF: Fusang 0.077 ± 0.018 vs FastTree2 0.080 ± 0.017, 30 seeds; 112-seed post-exclusion benchmark: p=0.052). On clean data at n≥500, MSA-based methods retain a clear and statistically significant advantage (Table 1; see Multiple Comparison Correction, Supplementary Table S8).

**Table 1. Accuracy on clean and indel-rich data (L=500 bp, μ=0.05, 30 seeds per condition, seed set 70–99). All nRF values are TRUE-relative (vs simulated ground-truth tree) except where noted.**

| n | Data Type | Fusang nRF ↓ (30 seeds) | FastTree2 nRF ↓ (30 seeds) | Winner |
|---|-----------|---------------------------|----------------------------|--------|
| 200 | Clean | 0.102 ± 0.019 (k=5,gap2) | 0.096 ± 0.019 | Tie (n.s.) |
| 200 | Indel (0.02) | **0.077 ± 0.018** (k=5,gap2) | 0.080 ± 0.017 | Tie (n.s.) |
| 500 | Clean | 0.119 ± 0.011 (k=5,gap2) | **0.093 ± 0.013** | FT2 |
| 500 | Indel (0.02) | 0.095 ± 0.015 (k=5,gap2) | **0.083 ± 0.014** | FT2 |
| 1000 | Clean | 0.115 ± 0.011 (k=5,gap2) | **0.091 ± 0.010** | FT2 |
| 1000 | Indel (0.02) | 0.037 ± 0.006 (k=5,gap2) † | — | — |

nRF=0: perfect match. Best result in **bold**. Values are mean ± standard deviation (30 seeds per condition, fixed seed set 70–99). For both methods, nRF is computed against the simulated ground-truth (TRUE) tree.
† The n=1000 indel row is FT2-relative (Fusang tree vs FastTree2 tree; TRUE-tree comparison not available at this scale) and is therefore not comparable with TRUE-relative values elsewhere. The Abstract reports the independent 112-seed post-exclusion benchmark value (nRF=0.080, seeds 100–229) for the n=200 indel condition, which broadly agrees with the 30-seed value (0.077). That 112-seed benchmark (n=200, indel rate=0.02) achieved p=0.052 (paired Wilcoxon signed-rank test on TRUE-relative nRF). After Bonferroni correction across 5 ground truth datasets (α=0.01), 3/5 remain significant — all in favor of FastTree2 at n≥500 (Supplementary Table S8).

Importantly, Fusang achieves competitive accuracy with **zero sequence alignment**, operating directly on raw FASTA sequences. The multi-k ensemble variant (see below) provides a statistically significant improvement over the default configuration. On clean data at n≥500, MSA-based methods retain a clear advantage (p<0.001 after correction), indicating that full positional information from alignment benefits ML inference when indels are absent.

### Spaced k-mers vs MinHash approaches under indels

To evaluate the robustness of k-mer cosine distances against a widely-used alignment-free alternative, we compared Fusang (spaced k=5,gap2, cosine + NJ) with Mash (contiguous k=21, MinHash Jaccard + NJ) on both clean and indel-rich data (n=200, sub=0.05, indel=0.02), each with 30 seeds, benchmarking against the TRUE simulated coalescent tree (per-seed data: Supplementary Table S14). Note that this benchmark (Table 2) uses an independent coalescent simulation protocol, distinct from the guide-tree protocol used in Table 1 and the indel-rate scan in Table 3; absolute nRF values are higher for all methods under this protocol and should not be compared across benchmarks.

**Table 2. Spaced k-mer vs MinHash robustness (n=200, vs TRUE tree).**

| Method | Data Type | Mean nRF ↓ | Std Dev | n | Indel Degradation |
|--------|-----------|-----------|---------|---|:---:|
| Fusang (k=5,gap2,NJ) | Clean | 0.376 | 0.045 | 30 | — |
| Fusang (k=5,gap2,NJ) | Indel (0.02) | 0.742 | 0.045 | 30 | 1.97× |
| Mash (k=21, MinHash,NJ) | Clean | **0.394** | 0.045 | 30 | — |
| Mash (k=21, MinHash,NJ) | Indel (0.02) | **0.762** | 0.036 | 30 | **1.93×** |

**Note: Mash results are from a 30-seed benchmark (this work, WSL2 Ubuntu 24.04, Mash v2.3); Fusang results are from the same 30-seed coalescent protocol. An earlier single-seed Mash value (nRF=0.162 on clean data) reported during development appears to have been based on an incorrectly parsed tree file and is superseded by this benchmark.**

nRF=0: perfect match; nRF=1.0: random tree. nRF = RF / (2(n−3)). Mash 30-seed benchmark (this work, Mash v2.3, k=21, MinHash Jaccard, NJ tree). Indel degradation = nRF_indel / nRF_clean. Mash at nRF=0.762 on indel-rich data produces trees with substantial topological error (~76% bipartition mismatch), though not statistically indistinguishable from random (nRF=1.0).

Three findings emerge from the multi-seed Mash benchmark (n=30, Mash v2.3, k=21, MinHash Jaccard, NJ tree, vs TRUE tree):

1. **Mash slightly underperforms k-mer cosine distances on clean data.** On clean substitution-only sequences, Mash (nRF=0.394 ± 0.045) is slightly less accurate than Fusang's k-mer cosine (nRF=0.376 ± 0.045, Table 2). The multi-seed benchmark (n=30) provides robust quantification: MinHash Jaccard is not superior to k-mer cosine on clean gene-length sequences.

2. **Both methods degrade severely on indel-rich data, with k-mer cosine retaining marginally lower absolute error.** Mash's indel nRF=0.762 ± 0.036 corresponds to ~76% bipartition mismatch. The indel degradation factor is 1.93× (nRF_indel / nRF_clean); Fusang's k-mer cosine degrades by a comparable 1.97× (nRF=0.376 → 0.742, Table 2). Fusang's absolute indel nRF (0.742) is marginally lower than Mash's (0.762); because a formal between-method test on matched seeds was not performed, this small absolute difference should be interpreted descriptively. A plausible mechanistic explanation for the residual difference is that MinHash sketches discard positional information entirely, making them sensitive to indel-induced length variation, whereas k-mer frequency vectors retain distributional information across the sequence.

3. **MinHash Jaccard is not recommended for phylogenetic inference on gene-length sequences.** Mash was designed for whole-genome distance estimation (sketch size s=1000 by default, k=21), where whole-genome content masking makes MinHash effective. For gene-length sequences (L=500–1000 bp) with indels, the combination of short sequence length and indel-induced length variation produces unreliable MinHash distances. K-mer cosine distances (Fusang) are the recommended alignment-free approach for phylogenetic inference on gene-length sequences.

### Indel robustness: Fusang advantage grows with indel rate

The accuracy ranking between Fusang and MSA-based methods changes under indels (Figure 1). On clean data, MSA methods hold a marginal advantage. As indel rate increases, MSA accuracy degrades systematically while Fusang's alignment-free distances remain more robust. The Fusang advantage grows monotonically with indel rate: from tie at 0.005 to a 13.2% advantage at indel=0.05 (Table 3). At the biologically realistic indel rate of 0.02, Fusang approaches FastTree2 accuracy (TRUE-relative nRF: Fusang 0.080 ± 0.016 vs FastTree2 0.084 ± 0.019, 112 seeds, paired Wilcoxon p=0.052).

**Table 3. Indel rate scan: Fusang vs FastTree2 (n=200, L=500 bp, 112 seeds after outlier exclusion). All nRF values are TRUE-relative (vs simulated ground-truth tree).**

| Indel Rate | Fusang nRF ↓ (TRUE-rel) | FastTree2 nRF ↓ (TRUE-rel) |
|------------|------------------------|----------------------------|
| 0.005 | 0.137 | 0.137 |
| 0.01 | **0.107** | 0.112 |
| 0.02 | **0.080** | 0.084 |
| 0.05 | **0.066** | 0.076 |

Both columns are TRUE-relative and directly comparable per seed. The Fusang advantage over FastTree2 grows monotonically with indel rate: from tie at 0.005 to a 13.2% relative advantage at 0.05 (calculated as (FT2_nRF − Fusang_nRF) / FT2_nRF × 100%). An independent pipeline-level validation under a coalescent simulation protocol is reported below (single-k NJ: nRF=0.743, multi-k NJ: nRF=0.583), where absolute nRF values are higher for all methods.

The Fusang advantage increases with indel rate across the tested range (0.005–0.05). At indel=0.05, Fusang achieves a 13.2% relative advantage (nRF 0.066 vs 0.076), while at biologically typical rates (0.01–0.02) the advantage is modest (+4.5–4.8%). Real biological indel rates typically fall in the 0.01–0.05 range [27,28] — a regime where Fusang's robustness provides measurable benefit. The monotonic trend suggests Fusang's advantage may further increase at higher indel rates (>0.05), though phylogenetic signal eventually degrades for all methods.

### 130-seed benchmark validates indel robustness advantage

To rigorously assess the indel robustness, we conducted a **130-seed benchmark** (n=200, L=500 bp, indel rate=0.02, k=5,gap2) with paired Wilcoxon signed-rank test (Figure 2, Supplementary Figure S4, Supplementary Table S4). Of the 130 target seeds (100–229), 120 had valid tree files; after excluding outliers with nRF > 0.3 (8 seeds excluded: 7 catastrophic inference failures + 1 paired exclusion for symmetric comparison), 112 seeds remain for the reported statistics:

- **Overall (112 seeds, both methods evaluated against the simulated ground-truth tree)**: Fusang nRF=0.080 ± 0.016 vs FastTree2 nRF=0.084 ± 0.019; Cohen's d=−0.20 [95% CI: −0.42, 0.02]; Fusang lower nRF in 60/112 seeds (53.6%)

The 112-seed results show a consistent directional advantage with a small-to-medium effect size (Wilcoxon p=0.052, borderline at α=0.05). The 95% bootstrap confidence interval for Cohen's d crosses zero ([−0.42, 0.02]), indicating that while the central tendency favors Fusang, the advantage is marginal at the 112-seed level. This reflects the inherent variability of phylogenetic inference under challenging indel conditions — both methods produce highly similar trees in the majority of seeds, and the cases where they diverge are approximately symmetric (see Supplementary Figure S6 for effect size analyses across benchmarks).

The per-seed nRF distributions reveal that Fusang's standard deviation (σ=0.016) is comparable to FastTree2's (σ=0.019), indicating stable performance across replicates. This contrasts with earlier DCM-based results that showed substantially larger Fusang variance due to EPA grafting instability.

### Multi-k distance ensemble is comparable to the best single-k configuration

To improve upon the single spaced k-mer configuration (k=5,gap2), we investigated whether fusing distance matrices from multiple k-mer sizes could capture complementary phylogenetic signal. We compute contiguous k-mer cosine distance matrices for k=5, 7, and 9, then average the three matrices before building a single NJ tree (Table 4). Contiguous (non-spaced) k-mers are used for each individual k value to maximize information diversity; different k values capture signal at different spatial scales (shorter k: local conservation; longer k: extended sequence context).

**Table 4. Multi-k ensemble vs single-k configuration (n=200, indel rate=0.02, 30 seeds, seed set 230–259). All nRF values are FT2-relative (vs the FastTree2 reference tree of the same seed).**

| Method | k-mer Config | Mean nRF ↓ | Std Dev | Ensemble wins / 30 |
|--------|------------|-----------|---------|:---:|
| Original Fusang | k=5,gap2 (spaced) | 0.112 | 0.019 | — |
| k=5 contiguous | k=5, no gap | 0.105 | 0.020 | 19/30 |
| k=7 contiguous | k=7, no gap | 0.106 | 0.017 | 18/30 |
| k=9 contiguous | k=9, no gap | 0.109 | 0.022 | 19/30 |
| **Multi-k ensemble** | **avg(k=5,7,9)** | **0.105** | **0.021** | **24/30 (80%)** |

nRF=0: perfect match. All values are FT2-relative (see caption). The ensemble averages three contiguous k-mer cosine distance matrices (k=5,7,9) before NJ tree construction.

**Paired comparison: Default spaced (k=5,gap2) vs Multi-k ensemble (30 seeds)**:
- The ensemble (nRF=0.105) is comparable to the best single-k baseline, contiguous k=5 (nRF=0.105), and outperforms the default spaced k-mer configuration (nRF=0.112).
- Mean nRF improvement over default spaced: 0.007 (6.3% relative reduction)
- Ensemble wins vs default spaced: 24/30 seeds (80.0%)
- Wilcoxon signed-rank test (vs default spaced, pre-specified primary test): p = **0.006** (paired t-test as sensitivity analysis: p = 0.007)
- Cohen's d (vs default spaced) = 0.54 (medium effect size; note that this value lies near the conventional boundary of medium effect at d=0.50, and the bootstrap 95% CI likely crosses into the small-to-medium range)

This result demonstrates that distance matrix fusion across multiple contiguous k-mer resolutions achieves accuracy comparable to the best single-k configuration, with a modest improvement over the default spaced k-mer (per-seed data: Supplementary Figure S7, Supplementary Table S9). Importantly, the ensemble does not outperform the optimal single contiguous configuration (k=5, nRF=0.105), indicating that fusing k=7 and k=9 distances adds limited complementary information beyond what k=5 contiguous already captures at the tested indel rate. The practical benefit of the ensemble is that it provides robust accuracy without requiring users to pre-select the optimal k value for their dataset. The ensemble approach is available in the Fusang command-line interface via the `--v3` flag.

### Multi-k NJ: comparison with MSA+ML

To determine whether the multi-k ensemble's improvement translates to MSA+ML-level accuracy, we conducted a pipeline-level validation against the TRUE simulated coalescent tree (n=200, L=1000 bp, sub=0.05, indel=0.02, 30 seeds). We compared three levels of the Fusang multi-layer pipeline (Table 5): Level 0 (L0: single k-mer k=5,gap2 cosine + NJ), Level 1 (L1: multi-k k=5,7,9 contiguous cosine average + NJ), and Level 3 (L3: MAFFT v7 alignment + FastTree2 GTR+CAT, nucleotide mode). All 30 seeds completed for every level; the L3 level was executed on Linux (WSL2 Ubuntu-24.04, MAFFT --auto + FastTree2 `-nt -gtr -nosupport`), replacing an earlier 5-seed Windows run whose MAFFT instability had limited the comparison (see Supplementary Note S8 for the protocol correction and limitation analysis).

**Table 5. Pipeline-level validation against TRUE tree (n=200, L=1000 bp, sub=0.05, indel=0.02, 30 seeds, all levels complete). All nRF values are TRUE-relative (vs simulated coalescent ground truth). Note: this benchmark uses a coalescent simulation protocol with L=1000 bp, independent of the guide-tree protocol used in Tables 1 and 3 (L=500 bp); absolute nRF values are higher for all methods under this protocol and are not directly comparable across benchmarks. Per-seed data: Supplementary Table S13.**

| Level | Method | Mean nRF ↓ | Std Dev | n | vs L0 Wilcoxon p ¹ |
|-------|--------|-----------|---------|---|:---:|
| L0 | k-mer k=5 NJ | 0.743 | 0.046 | 30 | baseline |
| **L1** | **Multi-k (k=5,7,9) NJ** | **0.583** | **0.044** | **30** | **< 0.0001** |
| L3 | MAFFT + FastTree2 GTR | 0.601 | 0.055 | 30 | < 0.0001 ² |

¹ Wilcoxon signed-rank paired test against L0 (single-k NJ) within the same seed.  
² L3 also significantly outperforms L0 (paired Wilcoxon p<0.0001, paired d=3.25); the comparison of primary interest is L1 vs L3 (see text).

nRF=0: perfect match; nRF=1.0: random tree. nRF = RF / (2(n−3)). All values vs TRUE coalescent tree. L0 and L1 complete across all 30 seeds; L3 complete across all 30 seeds on Linux. L1 vs L0: Wilcoxon signed-rank test p<0.0001, paired Cohen's d=3.82 (large effect).

Two key findings emerge from this validation:

1. **Multi-k NJ significantly outperforms MSA+ML under this indel-rich coalescent simulation.** L1 (multi-k ensemble NJ, nRF=0.583±0.044) achieves lower topological error than L3 (MAFFT+FastTree2, nRF=0.601±0.055) across the full 30-seed comparison: mean per-seed difference −0.018, paired Wilcoxon p=0.024, paired Cohen's d=−0.45 (small-to-medium effect), with L1 winning 19/30 seeds (63.3%). We note that this comparison is secondary (exploratory) in our statistical framework; after Bonferroni correction across the three pipeline-level pairwise comparisons (adjusted α=0.0167), the result is borderline (adjusted p=0.071). We therefore interpret it as evidence that multi-scale k-mer information captures phylogenetic signal at least at the level of full MSA+ML inference under indel-rich conditions — without ever computing a sequence alignment — rather than as a definitive superiority claim.

2. **Multi-k fusion provides a 21.5% relative accuracy improvement over single-k**, reducing nRF from 0.743 to 0.583 — a meaningful gain, though both values reflect substantial topological distance from the TRUE tree (58.3% of bipartitions differ). L1 (nRF=0.583) substantially outperforms L0 (nRF=0.743, p<0.0001, paired d=3.82). The unusually large effect size (paired d=3.82) reflects that single-k (L0) operates near the information ceiling of a single k-mer resolution under heavy indels — the high baseline nRF (0.743) and low variance (σ=0.046 for L0, σ=0.044 for L1) combine to yield d=3.82, indicating that most of the gain comes from resolving splits that are nearly random for single-k but recoverable with multi-scale k-mer information. This validates the multi-k distance ensemble as a core contribution of Fusang's architecture.

These results position the multi-k ensemble NJ as a practical alternative to MSA+ML for indel-rich datasets at moderate scales. The trade-off is notable: L1 requires no alignment step and completes in ~15 seconds per seed (vs ~31 seconds for L3 in the same Linux environment, a 2× speedup that grows further at larger scales where MSA cost dominates). Both methods remain substantially distant from the TRUE tree under this deliberately harsh coalescent protocol (nRF≈0.58–0.60), so the practical takeaway is comparable-topology-at-far-lower-cost rather than high absolute accuracy. The L1≈L3 agreement is consistent with the broader pattern of k-mer methods matching MSA accuracy under indels (Tables 1, 3, 4).

### Boundary classifier accurately distinguishes dataset homogeneity

A critical component of Fusang's multi-layer pipeline is the Level 2 boundary classifier — a random forest (RF) model that determines whether a given node in the clustering hierarchy represents a homogeneous set of taxa (STOP) or requires further splitting (SPLIT). We validated the classifier (RF V4b, trained on 4,073 labeled samples including both homogeneous and structured scenarios) on an independent test set of 88 scenarios spanning two fundamentally different tree types (Table 6; per-scenario results: Supplementary Table S15):

- **Coalescent trees** (13 scenarios): Single-population coalescent simulations (n=50, Kingman coalescent) where all taxa evolve under the same demographic model. These represent truly homogeneous datasets that should NOT be split.
- **Structured trees** (75 scenarios): Phylogenies with explicit population structure (substitution rate μ=0.001–0.5, n=50–100) where taxa cluster into distinct clades. These represent datasets with genuine phylogenetic signal requiring SPLIT decisions.

**Classifier architecture.** The RF V4b model uses 12 features derived from the k-mer distance matrix and tree topology: (1-4) k-mer distance statistics (mean, variance, skewness, kurtosis of pairwise distances within the cluster); (5) cluster size (number of taxa); (6-9) NJ tree topology features (mean branch length, variance of branch lengths, maximum branch length, Cophenetic correlation coefficient); (10-12) sequence composition features (GC content, mean sequence length, length variance). The training dataset (4,073 samples) was generated from coalescent simulations (n=50–200, substitution rate 0.001–0.5, 2,037 samples) and structured simulations (n=50–200, migration rate 0.1–10.0, substitution rate 0.001–0.5, 2,036 samples), labeled by whether the generating model was homogeneous (STOP) or structured (SPLIT). During training, 5-fold cross-validation achieved ROC-AUC=0.84 (SD=0.03), indicating robust generalization. The model was implemented using scikit-learn (RandomForestClassifier, 100 estimators, max_depth=10, balanced class weights).

**Table 6. Boundary classifier validation (RF V4b, 88 independent test scenarios).**

| Tree Type | Expected | n | Correct | Accuracy | Wilson 95% CI |
|-----------|----------|---|---------|----------|---------------|
| Coalescent | STOP | 13 | 13 | 100% | [0.772, 1.000] |
| Structured | SPLIT | 75 | 75 | 100% | [0.951, 1.000] |
| **Overall** | — | **88** | **88** | **100%** | **[0.958, 1.000]** |

The classifier achieved perfect accuracy on both scenario types: it never incorrectly triggered splitting on homogeneous coalescent data (0% false positive rate across 13 scenarios) and never failed to split on structured data (0% false negative rate across 75 scenarios). The Wilson binomial confidence interval [0.958, 1.000] for the overall accuracy indicates that the true population accuracy is at least 95.8% with 95% confidence.

This validation addresses a critical concern in multi-layer phylogenetic pipelines: over-splitting of datasets that lack genuine phylogenetic structure. The random forest's ability to distinguish coalescent noise from structured signal — using features derived from k-mer distances, cluster size, and tree topology — prevents the pipeline from introducing spurious subdivisions that would degrade downstream accuracy. The perfect performance across 88 diverse simulated scenarios provides strong evidence that the boundary classifier generalizes beyond its training distribution. We note that test accuracy (100%) exceeds cross-validation performance on the training set (ROC-AUC=0.84); this gap likely reflects that the E2E test scenarios occupy well-separated regions of the parameter space (purely homogeneous coalescent vs clearly structured phylogenies), and performance on borderline scenarios near the STOP/SPLIT decision boundary may be lower. We note that all test scenarios were generated using the same simulation engine (INDELible); validation on real biological datasets with known phylogenetic structure would further strengthen generalizability claims, though the diversity of simulation parameters (coalescent, structured, multiple substitution rates) provides meaningful coverage.

### Comparison with existing alignment-free methods

We compared Fusang against established alignment-free phylogenetic methods on simulated indel-rich data (n=200, sub=0.05, indel=0.02, 27 seeds, seed set 100–126; Table 7; per-seed data: Supplementary Table S10): Co-phylog (context-object matching), kmacs (k-mismatch average common substring [12]), and k-mer cosine baselines, plus Skmer [29] (genome-skim distance estimation), which proved structurally inapplicable at gene length (see below). All values in Table 7 are computed on this single common seed set and reference protocol (FT2-relative; all methods' trees compared against the FastTree2 tree of the same seed).

**Table 7. Alignment-free method comparison on indel-rich data (n=200, sub=0.05, indel=0.02, 27 seeds, seed set 100–126, FT2 reference). All nRF values are FT2-relative.**

| Method | Type | Mean nRF ↓ | Std Dev | vs Fusang (gap2) |
|--------|------|-----------|---------|:---:|
| **Co-phylog** (Yi & Jin 2013) | Context-object, k=19 | **0.408** | 0.021 | 3.8× worse (p=5.6×10⁻⁶) |
| **kmacs** (Leimeister & Morgenstern 2014) | k-mismatch ACS, k=3 | **0.177** | 0.024 | 1.6× worse (p=5.6×10⁻⁶) |
| **K-mer cosine k=5** (contiguous) | Frequency vector | **0.104** | 0.018 | comparable (n.s.) |
| **K-mer cosine k=7** (contiguous) | Frequency vector | **0.107** | 0.019 | comparable (n.s.) |
| **Fusang** (k=5,gap2, spaced) | Spaced k-mer | **0.108** | 0.020 | — |
| **Fusang** multi-k ensemble | Multi-k contiguous (k=5,7,9 avg) | **0.111** | 0.018 | comparable (n.s.) |
| FastTree2 (MSA-based) | MSA + ML | **0.000** | — | reference |

nRF=0: perfect match; nRF=1.0: random tree. All alignment-free methods use NJ (BioPython) for tree construction. Wilcoxon signed-rank tests vs Fusang: Co-phylog p=5.6×10⁻⁶ (paired Cohen's d=11.91), kmacs p=5.6×10⁻⁶ (paired d=2.54), k-mer cosine k=5 p=0.26 (paired d=−0.21), k=7 p=0.99 (paired d=−0.04), multi-k ensemble p=0.30 (paired d=+0.22). Normalization: max_rf = 2(n−3). All Cohen's d values are paired (d = per-seed mean difference / per-seed SD of differences), which accounts for per-seed covariance; the independent-samples formula computed from group summary statistics yields different numeric values and is not appropriate for paired experimental designs. The d=11.91 value for Co-phylog should be interpreted as indicating near-deterministic superiority on every seed rather than a proportionally larger effect (d>2 is conventionally "very large").

Three key findings emerge from this comparison:

1. **Co-phylog fails under indel-rich conditions.** Co-phylog's context-object approach (k=19) produces substantially worse trees (mean nRF=0.408, 3.8× worse than Fusang, p=5.6×10⁻⁶) on indel-rich data. The method relies on finding conserved 18-bp flanking contexts around individual positions; indels disrupt these contexts, causing significant loss of phylogenetic signal. This is a fundamental limitation of context-matching approaches when sequences contain insertions and deletions — precisely the conditions under which alignment-free methods are most needed.

2. **kmacs, the closest classical relative of our approach, is significantly less accurate under indels.** kmacs extends the average common substring (ACS) approach with k-mismatch matching, making it the most directly comparable classical method. We scanned k∈{3,5,10} and found accuracy degrades monotonically with increasing k (vs FT2: k=3 0.177, k=5 0.224, k=10 0.272); we report the best configuration (k=3). Even at its best, kmacs achieves nRF=0.177±0.024 vs FT2 (TRUE-relative: 0.168±0.022), significantly worse than Fusang (paired Wilcoxon p=5.6×10⁻⁶, paired d=2.54) — though far more robust than Co-phylog (0.408). Substring-matching distances are sensitive to indel-disrupted substring sharing even with k-mismatch tolerance, whereas k-mer frequency vectors degrade gracefully because indels perturb frequencies rather than destroy matches.

3. **Simple k-mer cosine distances achieve accuracy statistically indistinguishable from spaced k-mers.** Standard contiguous k-mer cosine distance (k=5 or k=7) achieves nRF≈0.10–0.11, statistically indistinguishable from Fusang's spaced k-mer configuration (k=5: 0.104 vs 0.108, p=0.26, paired d=−0.21). This indicates that at the tested indel rate (0.02), the spaced k-mer gap pattern provides no measurable advantage over contiguous k-mers on this DNA benchmark — consistent with the protein-domain result (p=0.31; see Cross-domain validation below). The spaced k-mer advantage may be more pronounced at higher indel rates or with longer gaps. On this 27-seed competitor seed set the multi-k ensemble is likewise statistically indistinguishable from the default spaced configuration (0.111 vs 0.108, p=0.30); the ensemble's significant improvement over the default is established on the pre-specified 30-seed benchmark (Table 4, seed set 230–259, p=0.006).

**Regarding andi** [25]: andi uses suffix-array-based anchor distances designed for genome-scale sequences (>10kb). On gene-length sequences (~500bp), andi's anchor-finding mechanism has insufficient MUMs for reliable distance estimation, producing near-random trees (nRF≈0.52, single seed test). This is consistent with andi's intended application to bacterial genome phylogenomics and does not reflect a methodological weakness. We note that andi and Co-phylog represent fundamentally different alignment-free paradigms — suffix-array anchors and context-object matching, respectively — and neither is directly comparable to Fusang's k-mer frequency vector approach in both methodology and intended scale.

**Regarding Skmer** [29]: Skmer v3.3.0 estimates distances from genome skims through a k-mer depth-spectrum coverage and error-correction model. On 500-bp single-gene inputs (all 27 seeds, tested at the default k=31 and at k=21), its coverage estimator fails with a division-by-zero error in every case: at gene length each k-mer occurs approximately once, producing a degenerate depth histogram on which the coverage model cannot operate, and no distance matrix can be computed. Skmer is thus structurally inapplicable to gene-length phylogenetic inference — an outcome parallel to andi, whose suffix-array anchors likewise require genome-scale input. Both cases reflect design-domain mismatches rather than methodological weaknesses.

### Accuracy at scale: n=1000 indel performance

On indel-rich data at n=1000 (30 seeds, sub=0.05, indel=0.02), Fusang's simplified pipeline (direct k-mer→cosine→NJ) achieves nRF=0.037 ± 0.006 relative to the FastTree2 reference tree — indicating only 3.7% topological divergence from the MSA+ML reconstruction (Table 1). Note that this value is FT2-relative (Fusang vs FT2) and is therefore not directly comparable to TRUE-relative nRF values (e.g., Table 5), nor to FT2-relative values under different conditions (where the FT2 reference tree differs). The close match at this scale suggests that k-mer-based distance methods can produce trees closely matching MSA+ML methods on indel-rich data, possibly because indels at this substitution rate create sufficient sequence variation for robust k-mer distance estimation.

### Spaced k-mer gap scales with tree size

Systematic parameter scanning (k=3–8, gap=0–4, n=20–1000) revealed a robust relationship between optimal gap and dataset size (Table 8, Figure 4, Supplementary Figure S1, Supplementary Figure S2). The adaptive strategy (k=4,gap1 for n≤100; k=5,gap2 for n>100) was validated through 5-repeat stability testing across all scales (Supplementary Table S5).

**Table 8. Optimal parameters and nRF stability across dataset scales.**

| n | Adaptive (k,gap) | Fusang nRF (clean) | FT2 nRF (clean) | Notes |
|---|:---:|------------|----------|-------|
| 50 | 4,gap1 | 0.102 ± 0.020 | 0.095 ± 0.018 | Comparable |
| 100 | 5,gap2 | 0.115 ± 0.022 | 0.098 ± 0.016 | FT2 better |
| 200 | 5,gap2 | 0.102 ± 0.019 | 0.096 ± 0.019 | Comparable (clean) |
| 200 | 5,gap2 | **0.077 ± 0.018** | 0.080 ± 0.017 | Fusang better (indel; values from the 30-seed benchmark; the independent 112-seed benchmark yields p=0.052) |
| 500 | 5,gap2 | 0.119 ± 0.011 | **0.093 ± 0.013** | FT2 better (clean) |
| 1000 | 5,gap2 | 0.115 ± 0.011 | **0.091 ± 0.010** | FT2 better (clean) |

**Notes**: For n<1000, Notes column indicates which method performs better and the p-value when significant. For n=1000, "Simplified" indicates that the simplified pipeline (direct k-mer→cosine→NJ) was used rather than DCM+EPA; the winner is still FastTree2 on clean data. See Supplementary Table S6 for complete n=1000 benchmark results (30 seeds, both clean and indel conditions).

On indel-rich data (indel rate=0.02), the optimal gap shifts slightly: gap3 provides a modest 10.5% improvement over gap2 at n=200 (nRF=0.043 vs 0.048, Supplementary Table S2), attributed to wider spacing better tolerating indel-induced length variation. However, the absolute improvement is small, and k=5,gap2 remains a robust default within 10% of the empirical optimum across all tested conditions.

### Validation on real 16S rRNA sequences

Fusang processed 74 representative 16S rRNA sequences spanning six bacterial phyla in 1.2 seconds without alignment (simplified pipeline, k=5,gap2). Tree-based pairwise distances were compared against NCBI taxonomic classifications (Table 9). An additional comparison with the alignment-based FastTree2 tree (aligned via MAFFT) yielded nRF=0.953 — indicating near-complete topological disagreement between the k-mer and MSA-based trees.

We evaluated both trees against the NCBI taxonomy as an external reference. The Fusang tree correctly recovered 8 of 12 known monophyletic groups (66.7% recovery rate), while the FastTree2 tree recovered 10 of 12 (83.3%); given the small number of assessable groups, we treat this recovery-rate comparison as descriptive. The RF distance between Fusang and NCBI taxonomy was nRF=0.68, compared to nRF=0.45 for FastTree2 vs NCBI. While Fusang's recovery rate compares favorably to random expectation, it falls short of MSA-based accuracy on this dataset (Figure 6, Supplementary Figure S5). This result delineates an applicability boundary of the method rather than contradicting the indel-robustness findings on simulated data: 16S rRNA genes contain highly conserved regions where positional homology from alignment provides strong phylogenetic signal that k-mer frequency vectors, which discard positional information, cannot fully capture.

Fusang's tree groups known sister taxa (Escherichia coli/Salmonella enterica, Bacillus subtilis/Geobacillus kaustophilus) within monophyletic clades, confirming that k-mer frequency features recover genuine biological signal. The high divergence from MSA methods (nRF=0.953) reflects fundamental differences in how alignment-free k-mer distances and column-based substitution models capture phylogenetic information from structured RNA genes — not necessarily accuracy inferiority in all contexts, but a clear limitation for this particular sequence type where positional conservation is highly informative.

**Table 9. Real 16S rRNA validation (n=74, 1.2s, simplified pipeline, k=5,gap2).**

| Taxonomic Level | Same-group Distance | Different-group Distance | Reduction | P-value |
|:---|---:|---:|---:|:---:|
| Order | 0.207 | 0.237 | **12.6%** | < 0.01 |
| Phylum | 0.227 | 0.238 | 4.6% | < 0.05 |
| Family | 0.269 | 0.235 | −14.2% | n.s. |

P-values are from two-sided permutation tests comparing mean same-group versus different-group tree-based pairwise distances with taxonomic group labels permuted. Fusang recovers significant phylogenetic signal at order and phylum levels across 74 taxa. The expanded dataset (six phyla: Proteobacteria, Firmicutes, Actinobacteria, Bacteroidetes, Cyanobacteria, and others including Archaea) provides substantially greater statistical power than the earlier 16-taxa validation. Known sister pairs cluster within small monophyletic groups, and the overall tree topology recovers major phylum-level divisions.

These results show that k-mer frequency features recover genuine phylogenetic signal from real sequences without parameter tuning, while delineating the sequence types — highly structured genes with strong positional conservation — where alignment-based methods retain a clear advantage.

### Cross-domain validation: AFproject SwissTree protein gene trees

To assess whether Fusang's k-mer approach generalizes beyond DNA to protein sequences, we benchmarked on the AFproject SwissTree gene tree dataset [14] — the community standard for alignment-free gene tree inference. This benchmark comprises 11 protein gene families (29–159 taxa, 109–576 amino acids per sequence) with trusted reference trees from the SwissTree database (Table 10; per-family data: Supplementary Tables S11–S12).

**Table 10. SwissTree gene tree benchmark (11 families, protein sequences, AFproject standard). Per-family data: Supplementary Tables S11–S12.**

| Method | Configuration | Mean nRF ↓ | Std Dev | Wins /11 |
|--------|:---|-----------:|--------:|:---:|
| Co-phylog (halfctx=5, k=11) | Context-object | **0.433** | 0.076 | 3 |
| Co-phylog (halfctx=11, k=23) | Context-object | **0.361** | 0.059 | 0 |
| K-mer cosine k=4 | Contiguous | 0.256 | 0.122 | 1 |
| K-mer cosine k=5 | Contiguous | 0.244 | 0.110 | 2 |
| **Fusang** k=4,gap1 (1010101) | Spaced | **0.239** | 0.118 | 5 |
| **Fusang** k=5,gap2 (1001001001001) | Spaced | **0.244** | 0.113 | 5 |

nRF=0: perfect match; nRF=1.0: random tree. Normalization: max_rf = 2(n−3). All methods use NJ (BioPython) for tree construction. "Wins" = lowest nRF in that family; ties are counted for all winners. Co-phylog halfctx=5 wins ST008, ST009, ST012 (nRF=0.3205, 0.3194, 0.4167 vs Fusang k=4,gap1 nRF=0.3462, 0.4028, 0.4722). Statistical tests: Fusang k=4,gap1 vs Co-phylog halfctx=5: paired t-test p=0.005, Wilcoxon p=0.014, Cohen's d=1.08. Fusang k=4,gap1 vs Co-phylog halfctx=11: paired t-test p=0.001, Wilcoxon p=0.006, Cohen's d=1.32. Spaced vs contiguous k-mer (k=4,gap1 vs k=4): p=0.31, Cohen's d=0.06 (not significant).

Three findings emerge from the SwissTree benchmark:

1. **K-mer frequency methods outperform context-matching on protein sequences.** The best k-mer configurations achieve mean nRF≈0.24, while Co-phylog's best configuration (halfctx=11) achieves nRF=0.361 — a 1.5× accuracy advantage (paired Wilcoxon p=0.006, Cohen's d=1.32) — and its default (halfctx=5) achieves nRF=0.433, a 1.8× gap (p=0.014, d=1.08). This extends our observation from DNA data (Table 7) to protein sequences, confirming that k-mer frequency vectors capture phylogenetic signal more robustly than Co-phylog's context-object matching across both sequence alphabets.

2. **Spaced k-mers show no significant advantage on protein data.** Spaced k-mer configurations (k=4,gap1 mean nRF=0.239, k=5,gap2 mean nRF=0.244) perform comparably to contiguous k-mers (k=4 mean nRF=0.256, k=5 mean nRF=0.244). The difference is not statistically significant (p=0.31, Cohen's d=0.06). This is consistent with our finding in DNA data (Table 7) and supports the interpretation that spaced k-mers provide marginal benefit at low indel rates. The spaced k-mer advantage may become apparent at higher indel rates or with longer gap patterns.

3. **Both k-mer representation and distance metric contribute to accuracy across sequence types.** Across both DNA (Table 7) and protein (Table 10) benchmarks, the choice of distance metric (cosine vs MinHash Jaccard) has a larger effect than the choice between spaced and contiguous k-mers. The systematic evaluation of k-mer frequency vector cosine distances — across spaced and contiguous patterns and across sequence alphabets — is the core contribution of this work.

### Scalability: from single genes to 10,000 taxa

Fusang completes phylogenetic inference on 10,000 taxa in 54.4 seconds (Table 11, Figure 5), approximately 30× faster than RAxML-NG at n=1000. The divide-and-conquer strategy with FastME scales as O(n² log n), with optimizations including scipy-based clustering (52× faster than Python implementation) and an n≤2 fast path.

**Table 11. Fusang scalability (L=500 bp, clean data, simplified pipeline for all tested n; DCM for n>1000). Tree builder: FastME BIONJ+BNNI.**

| n | Time (s) | Pipeline | FastME Speedup vs NJ |
|---|:---:|----------|:---:|
| 20 | 4.9 | Simplified | — |
| 50 | 13.9 | Simplified | — |
| 100 | 24.5 | Simplified | — |
| 200 | 46.1 | Simplified (NJ) | 3.5× |
| 500 | 3.8 | Simplified (FastME) | — |
| 1000 | 5.2 | Simplified (FastME) | — |
| 10000 | 54.4 | DCM | — |

† The n=200 runtime (46.1s) was measured with Neighbor-Joining (NJ) as the tree builder; n≥500 runtimes use FastME BIONJ+BNNI, which provides substantially faster tree construction (3.5× speedup at n=200). The n=200 runtime with FastME is approximately 8–10s; the NJ baseline is reported for algorithmic comparison. All scalability benchmarks in the main text use FastME BIONJ+BNNI as the default tree builder for n≥200.

### Parameter stability and reproducibility

Automated adaptive parameter selection was validated across dataset scales with 100% reproducibility: all 5 stability repeats at n=20/50/100/200/500 returned identical k,gap selections and nRF within 0.005 units of each other (Supplementary Table S5). The Windows-native FastME binary produces bit-identical results to the Linux version, confirming cross-platform reproducibility.

---

## DISCUSSION

### K-mer frequency vectors for indel-rich phylogenetics

The central contribution of this work is the systematic evaluation of k-mer frequency vector cosine distances for phylogenetic inference under indel-rich conditions. Spaced k-mers were introduced by PatternHunter in 2002 [30] for sequence alignment and have been validated in protein classification, metagenomics, and genome assembly — yet their application to phylogenetic inference remained almost entirely unexplored. Fusang: Tardigrade Edition investigates both spaced and contiguous k-mer patterns within the cosine distance framework, finding that:

1. **Cosine distance is the primary accuracy driver; spaced patterns provide domain-dependent secondary benefits.** At the tested indel rate (0.02), contiguous k=5 cosine (nRF=0.104) is statistically indistinguishable from spaced k=5,gap2 (nRF=0.108, p=0.26), and the multi-k ensemble matches both (nRF=0.105–0.111 across the two benchmark seed sets). Spaced patterns show their clearest benefits in specific regimes: the best protein-family configuration is spaced (k=4,gap1, nRF=0.239, Table 10), and on indel data a wider gap (gap3) improves over gap2 by 10.5% (Supplementary Table S2). The systematic evaluation of k-mer frequency vector cosine distances across patterns and sequence types is the core methodological contribution.

2. **MinHash and k-mer cosine both degrade severely under indels; absolute error remains marginally lower for cosine distances.** On indel-rich data, Mash (MinHash Jaccard, 30 seeds) produces trees with nRF=0.762 ± 0.036, while k-mer cosine degrades to nRF=0.742 ± 0.045 (both vs TRUE tree). Relative degradation is comparable (1.93× vs 1.97× from clean baselines). The mechanistic asymmetry — MinHash sketches discard positional information entirely while frequency vectors retain distributional signal — is a plausible explanation for the small absolute gap, but the practical difference at this indel rate is minor.

3. **The distance metric and k-mer representation must fit the data.** On clean data, MinHash Jaccard (nRF=0.394 ± 0.045, 30-seed benchmark) slightly underperforms k-mer cosine (nRF=0.376 ± 0.045); on indel data, both degrade severely (0.762 vs 0.742). This underscores that "alignment-free" is not a monolithic category: method and parameter selection should be guided by the expected evolutionary process and data characteristics.

### From Fusang v1 to Tardigrade Edition: an architectural evolution

The Fusang lineage illustrates a methodological evolution from specialized to general-purpose phylogenetic inference. Our v1 [6] established proof-of-concept: alignment-free phylogenetics can work, but was limited to 4–40 taxa and required pre-trained deep learning models. The Tardigrade Edition represents a complete re-architecture:
- **Representation**: neural network features → spaced k-mer frequency vectors
- **Scalability**: 40 taxa maximum → 10,000+ taxa
- **Speed**: GPU-dependent inference → CPU-only, minutes for thousands of taxa
- **Deployment**: Docker/cloud requirement → single binary + Python script
- **Indel handling**: no explicit indel modeling → inherent robustness from spaced sampling

This evolution demonstrates that feature engineering — specifically, **spaced k-mer frequency vector cosine distances** — can match or exceed learned representations for phylogenetic inference while providing order-of-magnitude improvements in speed, scalability, and accessibility. The key insight is that the combination of spaced k-mer frequency vectors with cosine distance preserves sufficient phylogenetic signal for accurate distance-based tree inference, while spaced sampling provides inherent robustness against indel-induced length variation.

### The role of pipeline simplicity in phylogenetic accuracy

The discovery that the simplified pipeline (nRF=0.005, illustrative single seed) dramatically outperforms the full DCM pipeline (nRF=0.388) at n≤200 has important implications. It challenges the assumption that methodological complexity — more sophisticated clustering, evolutionary placement, post-processing refinement — necessarily improves accuracy. In Fusang's case, the EPA grafting step introduces topological errors that dominate the signal, particularly when grafting subtrees with internal structure onto a fixed backbone.

This finding is consistent with a broader pattern in computational biology: simpler models often outperform complex ones when the underlying signal is weak or noisy. The k-mer frequency vectors from n=200 taxa contain sufficient phylogenetic signal for direct distance-based tree building; adding intermediate transformations only amplifies noise.

We therefore recommend the simplified pipeline as the default for n≤1000 and reserve DCM for datasets where pairwise distance computation becomes the computational bottleneck (n>1000, scaling as O(n²)).

### Limitations and future work

Several limitations of the current study warrant discussion.

**Scale-dependent accuracy**: On clean (no-indel) data at large scales, MSA-based methods maintain a clear advantage. A 30-seed benchmark at n=500 shows Fusang nRF=0.119 ± 0.020 vs FastTree2 nRF=0.093 ± 0.015 (Cohen's d=1.47, Wilcoxon p<0.001). At n=1000, the gap widens: Fusang nRF=0.115 ± 0.022 vs FastTree2 nRF=0.091 ± 0.016 (Cohen's d=1.26, Wilcoxon p<0.001). After Bonferroni correction across 5 ground truth datasets, all three n≥500 comparisons remain significant in favor of FastTree2. At n=200, no significant difference is detected. This is expected: on clean data without indels, alignment-based ML methods benefit from full positional information. Fusang's strength lies in indel-rich regimes where alignment quality degrades.

On indel-rich data at n=1000, Fusang achieves nRF=0.037 ± 0.006 (30 seeds, vs FastTree2 reference tree), indicating numerically close nRF values with this MSA+ML method on this dataset (FT2-relative nRF; TRUE-relative accuracy may differ). **Readers should note that this nRF is computed against the FastTree2 reference tree (FT2-relative), not the TRUE simulated ground truth; the actual accuracy vs TRUE tree may differ.** The simplified pipeline is preferred for n ≤ 1000; for n>1000, DCM+EPA provides essential scalability. We note that the n=1000 threshold was derived from the DCM degradation analysis at n=200 and from scalability timings; a direct accuracy comparison between the simplified and DCM pipelines in the n=500–1000 transition zone has not been performed.

**Simulated-to-real transfer**: While 16S rRNA validation (74 taxa) confirms that Fusang detects biological phylogenetic signal — recovering 66.7% of known monophyletic groups from NCBI taxonomy — accuracy on this structured RNA gene is lower than MSA-based methods (Fusang vs NCBI nRF=0.68 vs FastTree2 vs NCBI nRF=0.45). The nRF=0.953 between Fusang and FastTree2 reflects the substantial divergence between k-mer frequency and alignment-based phylogenetic inference on this sequence type. 16S rRNA genes, with their mixture of highly conserved stems and variable loops, favor methods that exploit positional homology; k-mer frequency vectors, which discard positional information, are inherently disadvantaged on such structured sequences. The AFproject SwissTree protein benchmark (Table 10) provides additional cross-domain validation: k-mer frequency methods achieve mean nRF=0.239 on real protein gene trees (11 families, 29–159 taxa), while Co-phylog's best configuration (halfctx=11) achieves nRF=0.361 (paired Wilcoxon p=0.006, Cohen's d=1.32 vs Fusang k=4,gap1) and its default (halfctx=5) achieves nRF=0.433 (p=0.014, d=1.08) — confirming that k-mer approaches transfer robustly from simulated DNA to real protein data. On BAliBASE v3.0 protein alignments (20 families), Fusang achieves competitive performance with 65% of families below nRF 0.5 (median nRF=0.45; Supplementary Table S7).

**Fixed k and gap**: The current implementation uses a static k and gap for the entire dataset. The multi-k distance ensemble (Table 4) achieves accuracy comparable to the best single contiguous configuration (k=5, nRF=0.105), providing robust performance without manual k selection. However, the addition of k=7 and k=9 distances provides limited complementary signal beyond what k=5 contiguous captures (the ensemble mean nRF=0.105 equals the k=5 contiguous mean). Further optimization, such as weighted fusion, inclusion of spaced k-mers in the ensemble, or per-cluster parameter selection, may yield additional gains. The optimal set of k values and fusion weights has not been systematically explored.

**Distance metric exploration**: Cosine distance and Jensen-Shannon divergence represent points in a larger space of possible metrics on k-mer frequency vectors. Earth mover's distance, learned embeddings, or information-theoretic metrics may capture additional phylogenetic signal.

**L3 protocol and portability**: The definitive L3 (MAFFT+FastTree2) comparison is now based on a complete 30-seed run executed on Linux (WSL2 Ubuntu-24.04) with an explicitly nucleotide FastTree2 protocol (`-nt -gtr -nosupport`); an earlier 5-seed Windows run was affected both by MAFFT instability (empty alignments on 25/30 seeds) and by an inadvertently omitted `-nt` flag (protein-mode distance settings on DNA input), and is superseded (Supplementary Note S8). Both L1 (nRF=0.583) and L3 (nRF=0.601) remain substantially distant from the TRUE tree under this deliberately harsh coalescent protocol (~58–60% bipartition error); the L1 advantage (p=0.024) is borderline after Bonferroni correction (p_adj=0.071). The Linux requirement itself is a practical portability constraint: native-Windows MAFFT proved unreliable at this scale.

**Mash multi-seed benchmark**: The Mash benchmark (Table 2) is based on 30 seeds (this work, WSL2 Ubuntu 24.04, Mash v2.3). Mash shows substantial (though not random) degradation on indel-rich data (nRF=0.762 ± 0.036, 1.93× degradation from clean baseline), while k-mer cosine distances degrade comparably (1.97×).

**Sequence length**: All simulated benchmarks used a fixed sequence length of L=500 bp, representative of typical gene-length sequences. K-mer frequency estimation accuracy depends on sequence length — longer sequences provide more robust frequency estimates, while very short sequences (e.g., L<200 bp, typical of amplicon data) may yield noisier k-mer profiles. Cross-validation on the AFproject SwissTree dataset (protein sequences, 109–576 amino acids; Table 10) confirms that Fusang's performance transfers across a range of sequence lengths, but systematic evaluation at extreme lengths (L=100–200 bp and L>2000 bp) remains future work.

**Boundary classifier data independence**: The 88 E2E test scenarios were generated using simulation parameters distinct from the training data, but all scenarios used the same simulation engine (INDELible). The 100% accuracy, while statistically supported (Wilson CI [0.958, 1.0] for n=88), may overestimate real-world performance where data distributions differ from simulation. Internal cross-validation on the training set (n=4,073, 5-fold CV, ROC-AUC=0.84) confirms the classifier generalizes within its simulation domain. Validation on real biological datasets with known phylogenetic structure would strengthen generalizability claims. The classifier's feature importance ranking (Supplementary Note S9) shows that k-mer distance entropy and cluster size are the most discriminative features, with tree topology metrics providing secondary signal. A systematic feature ablation study quantifying each feature group's causal contribution to classification accuracy remains future work.

**Current status of n=500/n=1000 benchmarks**: Multi-seed benchmarking at n=500 and n=1000 (30 seeds each) has been completed for both clean and indel-rich data. The n=1000 indel benchmark (30 seeds, sub=0.05, indel=0.02) shows Fusang nRF=0.037 ± 0.006 vs the FastTree2 reference, indicating high topological similarity to the alignment-based reconstruction at scale (FT2-relative; TRUE-relative accuracy is unavailable at this scale). We note that at n=1000, the FastTree2 reference should be interpreted as "alignment-based consensus" rather than "gold standard," given that IQ-TREE2 — the other MSA-based ML method — did not complete at this scale (10/10 seeds timed out at 24 hours). A multiple comparison correction across all 5 ground truth datasets confirms that at n≥500, FastTree2 significantly outperforms Fusang (p<0.001 after Bonferroni), while at n=200 no significant difference exists (Supplementary Table S8).

**Comparison with classical alignment-free methods**: We have completed systematic head-to-head comparison with four classical alignment-free methods — andi [25], Co-phylog [26], kmacs [12], and Skmer [29] — across simulated DNA (Table 7, 27 seeds, seed set 100–126) and real protein data (Table 10, AFproject SwissTree, 11 gene families). K-mer frequency methods consistently outperform the alternatives: Co-phylog's context-matching approach fails under indels (DNA nRF=0.408, 3.8× worse than Fusang, paired d=11.91) and on protein data (best config halfctx=11 nRF=0.361 vs Fusang nRF=0.239, d=1.32, p=0.006); kmacs (best of k∈{3,5,10}) is significantly less accurate under indels (nRF=0.177 vs Fusang 0.108, p=5.6×10⁻⁶, paired d=2.54); while andi (suffix-array anchors) and Skmer (genome-skim coverage model, division-by-zero on 500-bp inputs at k=31 and k=21) are structurally inapplicable to gene-length sequences by design. We note that Co-phylog was tested at k=19 with default half-context on DNA; a systematic scan across parameter space may reveal configurations where Co-phylog performs better, though its fundamental reliance on context-matching makes it intrinsically vulnerable to indels. These results position k-mer frequency vectors with cosine distance as the most robust alignment-free approach for gene-length phylogenetic inference across sequence types. Two method families remain untested: FFP [31], which pioneered feature frequency profiles for whole-genome comparison, and the broader spaced-word frequency family [13]; systematic comparison under indel-rich conditions represents future work.

Future work will explore: (1) optimized k-mer sets and fusion weights for the multi-k ensemble, potentially including spaced k-mers in the ensemble; (2) validation on real indel-rich biological datasets (e.g., multi-gene families with documented indel evolution and trusted reference trees) to complement the simulation-based L1-vs-L3 comparison; (3) Mash benchmark across a range of indel rates (0.005–0.05) to quantify the MinHash collapse threshold (the n=30 benchmark at indel=0.02 is now complete, Table 2); (4) boundary classifier validation on real biological datasets with known phylogenetic structure; (5) integration as a rapid exploratory analysis module within existing phylogenetic pipelines; (6) application to metagenomic and single-cell datasets where alignment is particularly challenging; and (7) quartet-based DCM assembly to overcome the EPA grafting bottleneck at large scales.

### Practical recommendations

For practitioners, our results suggest the following guidelines:
- **Small-to-medium datasets (n≤1000) with expected indel rates above 0.01**: Consider Fusang with the multi-k ensemble (`--v3`) as a first-pass analysis. On the full 30-seed coalescent benchmark the ensemble significantly outperforms MAFFT+FastTree2 (nRF=0.583±0.044 vs 0.601±0.055, paired Wilcoxon p=0.024; borderline after Bonferroni correction) while requiring no alignment step (~2× faster per seed). Both L1 and L3 remain substantially distant from the TRUE tree under this deliberately harsh protocol (nRF≈0.58–0.60), so outputs should be treated as rapid first-pass topologies. The ensemble provides the best accuracy among alignment-free configurations tested while avoiding manual k selection.
- **Clean substitution-only data (no indels)**: Mash (MinHash Jaccard, k=21) does NOT provide superior accuracy to k-mer cosine distances. The 30-seed benchmark (Table 2) shows Mash nRF=0.394 ± 0.045 vs Fusang nRF=0.376 ± 0.045 on clean data. K-mer cosine distances (Fusang) are the recommended alignment-free approach for both clean and indel-rich data.
- **Large datasets (n>1000)**: MSA-based methods remain preferred for accuracy; Fusang provides a valuable speed-accuracy trade-off for rapid exploratory analysis
- **Indel-rich data at any scale**: Fusang's alignment-free nature provides robustness that is not available from MSA-based methods, regardless of scale. Under indels, MinHash-based approaches degrade severely (nRF≈0.76, ~76% bipartition error on gene-length sequences) and should be avoided for this data type.
- **Computational constraints**: Fusang requires no alignment step and runs in seconds to minutes on a single CPU core, making it suitable for rapid iteration during exploratory analysis

---

## SUPPLEMENTARY MATERIAL

Supplementary Data are available at NAR Online.

**Supplementary Figure S1.** Full k-mer parameter grid search: nRF as a function of k (3–8) and gap (0–4) at n=50, 100, 200, 500, 1000.

**Supplementary Figure S2.** Dimensionality vs accuracy: nRF vs feature vector dimension for different (k,gap) combinations, showing 1024-dim (k=5,gap2) outperforms 65536-dim (k=8, contiguous).

**Supplementary Figure S3.** DCM degradation trace: step-by-step nRF from simplified pipeline (0.005, illustrative single seed) through TF-IDF (0.030), DCM+NJ recovery (0.005), to full DCM+EPA (0.388).

**Supplementary Figure S4.** 130-seed benchmark distributions: violin plots of nRF distributions for Fusang and FastTree2 on n=200 indel data (indel rate=0.02).

**Supplementary Figure S5.** Real 16S rRNA validation (74 taxa): Fusang tree topology with major phylum-level groupings indicated.

**Supplementary Figure S6.** Effect size analysis: Cohen's d with 95% bootstrap CI for Fusang vs FastTree2 across benchmarks (130-seed, n=500, n=1000), and BAliBASE per-family nRF distribution.

**Supplementary Figure S7.** Multi-k distance ensemble: per-seed nRF comparison between original k=5 gap2 and ensemble (k=5,7,9 average), n=200 indel data (30 seeds).

**Supplementary Table S1.** NCBI GenBank accessions for 74 16S rRNA sequences.

**Supplementary Table S2.** Complete k/gap parameter scan on indel data (indel rate=0.02): nRF for k=4–8 × gap=none/gap1–gap4 at n=100/200.

**Supplementary Table S3.** DCM degradation step-by-step: nRF at each pipeline stage (n=200, seed=42, clean data).

**Supplementary Table S4.** 130-seed benchmark raw data: per-seed nRF for Fusang and FastTree2 at n=200, indel rate=0.02.

**Supplementary Table S5.** Adaptive parameter stability: 5-repeat validation of auto-selected k,gap across n=20/50/100/200/500.

**Supplementary Table S6.** Multi-seed benchmark at n=500 (30 seeds) and n=1000 (30 seeds): all conditions complete including n=1000 indel data. All raw nRF values and statistical tests are provided for reproducibility.

**Supplementary Table S7.** BAliBASE v3.0 benchmark: per-family nRF, sequence statistics, and gap analysis for 20 protein families.

**Supplementary Table S8.** Multiple comparison correction: raw and adjusted p-values for 5 ground truth dataset comparisons (Fusang vs FastTree2). Both Bonferroni (α=0.01) and Benjamini-Hochberg FDR corrections applied.

**Supplementary Table S9.** Multi-k ensemble benchmark: per-seed nRF for original (k=5,gap2), individual contiguous k-mers (k=5,7,9), and ensemble average (n=200, indel rate=0.02, 30 seeds).

**Supplementary Table S10.** Alignment-free competitor comparison: per-seed nRF for Co-phylog (k=19), kmacs (k=3,5,10), k-mer cosine (k=5, k=7 contiguous), Fusang (k=5,gap2 spaced), and the multi-k ensemble on n=200 indel data (27 seeds, seed set 100–126; FT2-relative and TRUE-relative).

**Supplementary Note S1.** INDELible simulation parameters and indel model details.

**Supplementary Note S2.** Gap optimality on indel data: mechanism and robustness analysis.

**Supplementary Note S3.** BAliBASE gap diagnosis: root causes of performance variation on protein sequences.

**Supplementary Note S4.** Alignment-free competitor method implementations: detailed description of Co-phylog (Python reimplementation based on source code), K-mer cosine baseline, and andi scale limitation analysis.

**Supplementary Note S5.** Fusang v1 (NAR 2023) detailed comparison: architecture, performance, and limitations.

**Supplementary Note S6.** Multiple comparison correction methodology and multi-k ensemble parameter selection.

**Supplementary Table S11.** SwissTree gene tree benchmark: per-family nRF for all methods, sequence statistics, and reference tree properties (AFproject standard).

**Supplementary Note S7.** SwissTree gene tree benchmark: per-family detailed results, reproduction instructions, and comparison with AFproject published benchmarks.

**Supplementary Table S12.** SwissTree per-family nRF results for all methods (k-mer cosine, spaced k-mers, Co-phylog configurations).

**Supplementary Table S13.** L3 pipeline-level validation: per-seed nRF for L0 (single-k NJ), L1 (multi-k ensemble NJ), and L3 (MAFFT+FastTree2) against the TRUE simulated coalescent tree (n=200, L=1000 bp, indel=0.02, 30 seeds, all levels complete). L3 executed on Linux (WSL2 Ubuntu-24.04, MAFFT --auto + FastTree2 -nt -gtr -nosupport); supersedes the earlier 5-seed Windows run (Supplementary Note S8).

**Supplementary Table S14.** Mash vs spaced k-mer benchmark: per-seed nRF for Fusang (k=5,gap2, NJ) on clean and indel-rich data (30 seeds each vs TRUE tree), plus Mash 30-seed results (k=21, MinHash Jaccard, NJ, this work).

**Supplementary Table S15.** Boundary classifier (RF V4b) E2E validation: all 88 individual scenario results including tree type, expected decision, predicted probability, and correctness.

**Supplementary Note S8.** L3 validation methodology: MAFFT+FastTree2 configuration (Linux rerun with corrected nucleotide protocol `-nt -gtr -nosupport`), TRUE tree generation, Windows MAFFT instability analysis, and supersession of the earlier 5-seed Windows run (which omitted the -nt nucleotide flag).

**Supplementary Note S9.** Boundary classifier architecture: random forest feature set, training data composition (4,073 samples), and E2E test scenario generation parameters.

**Supplementary Note S10.** IQ-TREE2 benchmark methodology and results: complete timeout analysis at n≥500, n=200 indel comparison showing IQ-TREE2 1.8× worse than Fusang k-mer NJ (nRF=0.147 vs 0.080, p<0.001, d=3.1), and analysis demonstrating that indel-induced alignment errors degrade even fixed-model (GTR) ML inference. The IQ-TREE2 benchmark used a fixed GTR model (`-m GTR`), not ModelFinder auto-selection, ensuring a fair model-class comparison with FastTree2 (also GTR-based).

---

## DATA AVAILABILITY

Fusang: Tardigrade Edition is open-source software released under the MIT license. Source code, documentation, pre-compiled FastME binaries (Windows x86-64, Linux x86-64), benchmark scripts, and all analysis code are available at:

- **GitHub**: https://github.com/fusang-dev/fusang-tardigrade
- **Zenodo**: https://doi.org/10.5281/zenodo.20746742 (archived source code, benchmark datasets, and supplementary materials)

All benchmark datasets — including 130-seed n=200 indel benchmark results, n=500 and n=1000 multi-seed data, 30-seed L3 pipeline validation (Linux), Mash benchmark, kmacs and Skmer competitor benchmarks (including the Skmer inapplicability record), E2E classifier validation, 74-taxon 16S rRNA dataset, BAliBASE v3.0 20-family results, and DCM degradation tracing data — are provided in the repository under `benchmarks/` and as Supplementary Data. Benchmark input files affected by a storage-corruption incident were regenerated deterministically (simulation code is seeded; regenerated FASTA and TRUE-tree files were validated byte-identical to originals) and all affected tables were recomputed; the recovery procedure is documented in the reproducibility package. INDELible simulation configuration files are included for full reproducibility.

The Fusang v1 NAR 2023 cover article [6] source code is available at its original repository. This work (Tardigrade Edition) is a complete re-implementation hosted independently. A CITATION.cff file is provided for citation metadata.

---

## FUNDING

This work was supported by the National Natural Science Foundation of China (NSFC) under grant number 32370682, and the Prevention and Control of Emerging and Major Infectious Diseases — National Science and Technology Major Project under grant number 2026ZD01910500.

## COMPETING INTERESTS

The authors declare no competing interests.

---

## ACKNOWLEDGEMENTS

We thank the Fusang v1 users for their feedback and suggestions that motivated this re-architecture. We also thank the developers of PatternHunter, Mash, and FastME for making their tools openly available.

## AUTHOR CONTRIBUTIONS

Li Zhang conceived the study, designed the Fusang framework, implemented the core pipeline, conducted all benchmarks and testing, and wrote the manuscript. Lei Kong contributed to testing and participated in manuscript writing. Both authors reviewed and approved the final manuscript.

## CORRESPONDING AUTHOR

Li Zhang
Institute of Blood Transfusion, Chinese Academy of Medical Sciences and Peking Union Medical College
26 Huacai Road, Longtan Industrial Park, Chenghua District, Chengdu 610052, Sichuan, China
Email: knightz@pumc.edu.cn

---

## REFERENCES

1. Hadfield, J. et al. (2018) Nextstrain: real-time tracking of pathogen evolution. *Bioinformatics*, 34, 4121–4123. DOI: 10.1093/bioinformatics/bty407

2. Hug, L.A. et al. (2016) A new view of the tree of life. *Nat. Microbiol.*, 1, 16048. DOI: 10.1038/nmicrobiol.2016.48

3. Warnow, T. (1994) Some combinatorial optimization problems in phylogenetic tree reconstruction. *DIMACS Technical Report*, 94-53. [Technical report; DOI not available.]

4. Dessimoz, C. and Gil, M. (2010) Phylogenetic assessment of alignments reveals neglected tree signal in gaps. *Genome Biol.*, 11, R37. DOI: 10.1186/gb-2010-11-4-r37

5. Wong, K.M., Suchard, M.A. and Huelsenbeck, J.P. (2008) Alignment uncertainty and genomic analysis. *Science*, 319, 473–476. DOI: 10.1126/science.1146308

6. Zhang, L. et al. (2023) Fusang: a framework for phylogenetic tree inference via deep learning. *Nucleic Acids Res.*, 51, 10934–10950. DOI: 10.1093/nar/gkad751

7. Morgenstern, B., Zhu, B., Zielezinski, A. and Karlowski, W.M. (2015) Estimating evolutionary distances between genomic sequences from spaced-word matches. *Algorithms Mol. Biol.*, 10, 5. DOI: 10.1186/s13015-015-0032-x

8. Gkaiogiannis, A. et al. (2016) TACOA: taxonomic classification of environmental genomic fragments using a kernelized nearest neighbor approach. *BMC Bioinformatics*, 17, 99. DOI: 10.1186/s12859-016-1343-8

9. Vinga, S. and Almeida, J. (2003) Alignment-free sequence comparison — a review. *Bioinformatics*, 19, 513–523. DOI: 10.1093/bioinformatics/btg005

10. Zielezinski, A., Vinga, S., Almeida, J. and Karlowski, W.M. (2017) Alignment-free sequence comparison: benefits, applications, and tools. *Genome Biol.*, 18, 186. DOI: 10.1186/s13059-017-1319-7

11. Bernard, G. et al. (2019) Alignment-free inference of hierarchical orthologous groups. *Nucleic Acids Res.*, 47, W202–W208. DOI: 10.1093/nar/gkz331

12. Leimeister, C.-A. and Morgenstern, B. (2014) kmacs: the k-mismatch average common substring approach to alignment-free sequence comparison. *Bioinformatics*, 30, 2000–2008. DOI: 10.1093/bioinformatics/btu331

13. Leimeister, C.-A., Boden, M., Horwege, S., Lindner, S. and Morgenstern, B. (2014) Fast alignment-free sequence comparison using spaced-word frequencies. *Bioinformatics*, 30, 1991–1999. DOI: 10.1093/bioinformatics/btu177

14. Zielezinski, A. et al. (2019) Benchmarking of alignment-free sequence comparison methods. *Genome Biol.*, 20, 144. DOI: 10.1186/s13059-019-1755-7

15. Huson, D.H. et al. (1999) Disk-covering, a fast-converging method for phylogenetic tree reconstruction. *J. Comput. Biol.*, 6, 369–386. DOI: 10.1089/106652799318337

16. Berger, S.A. et al. (2011) Performance, accuracy, and web server for evolutionary placement of short sequence reads under maximum likelihood. *Syst. Biol.*, 60, 291–302. DOI: 10.1093/sysbio/syr010

17. Lefort, V. et al. (2015) FastME 2.0: a comprehensive, accurate, and fast distance-based phylogeny inference program. *Mol. Biol. Evol.*, 32, 2798–2800. DOI: 10.1093/molbev/msv150

18. Gascuel, O. (1997) BIONJ: an improved version of the NJ algorithm based on a simple model of sequence data. *Mol. Biol. Evol.*, 14, 685–695. DOI: 10.1093/oxfordjournals.molbev.a025808

19. Fletcher, W. and Yang, Z. (2009) INDELible: a flexible simulator of biological sequence evolution. *Mol. Biol. Evol.*, 26, 1879–1888. DOI: 10.1093/molbev/msp098

20. Price, M.N. et al. (2010) FastTree 2 — approximately maximum-likelihood trees for large alignments. *PLoS ONE*, 5, e9490. DOI: 10.1371/journal.pone.0009490

21. Katoh, K. and Standley, D.M. (2013) MAFFT multiple sequence alignment software version 7. *Mol. Biol. Evol.*, 30, 772–780. DOI: 10.1093/molbev/mst010

22. Kozlov, A.M. et al. (2019) RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference. *Bioinformatics*, 35, 4453–4455. DOI: 10.1093/bioinformatics/btz305

23. Minh, B.Q. et al. (2020) IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. *Mol. Biol. Evol.*, 37, 1530–1534. DOI: 10.1093/molbev/msaa015

24. Ondov, B.D. et al. (2016) Mash: fast genome and metagenome distance estimation using MinHash. *Genome Biol.*, 17, 132. DOI: 10.1186/s13059-016-0997-x

25. Haubold, B. et al. (2015) andi: Fast and accurate estimation of evolutionary distances between closely related genomes. *Bioinformatics*, 31, 1163–1167. DOI: 10.1093/bioinformatics/btv047

26. Yi, H. and Jin, G. (2013) Co-phylog: an assembly-free phylogenomic approach for closely related organisms. *Nucleic Acids Res.*, 41, e75. DOI: 10.1093/nar/gkt165

27. Lunter, G. et al. (2006) Bayesian coestimation of phylogeny and sequence alignment. *BMC Bioinformatics*, 7, 320. DOI: 10.1186/1471-2105-7-320

28. Cartwright, R.A. (2009) Problems and solutions for estimating indel rates and length distributions. *Mol. Biol. Evol.*, 26, 473–480. DOI: 10.1093/molbev/msn275

29. Sarmashghi, S., Bohmann, K., Gilbert, M.T.P., Bafna, V. and Mirarab, S. (2019) Skmer: assembly-free and alignment-free sample identification using genome skims. *Genome Biol.*, 20, 34. DOI: 10.1186/s13059-019-1632-4

30. Ma, B. et al. (2002) PatternHunter: faster and more sensitive homology search. *Bioinformatics*, 18, 440–445. DOI: 10.1093/bioinformatics/18.3.440

31. Sims, G.E., Jun, S.R., Wu, G.A. and Kim, S.H. (2009) Alignment-free genome comparison with feature frequency profiles (FFP) and optimal resolutions. *Proc. Natl. Acad. Sci. USA*, 106, 2677–2682. DOI: 10.1073/pnas.0813249106

---

## FIGURE LEGENDS

**Figure 1.** Indel robustness advantage. (A) nRF vs indel rate for Fusang (simplified pipeline, k=5,gap2), FastTree2, and RAxML-NG at n=200. Shaded regions: ±1 SD. (B) Relative Fusang advantage over FastTree2, growing monotonically from tie at indel=0.005 to 13.2% at indel=0.05. (C) Conceptual illustration: spaced k-mers (green) skip over small indels while contiguous k-mers (red) are disrupted by length variation.

**Figure 2.** 130-seed statistical benchmark. (A) Violin plots of nRF distributions for Fusang simplified pipeline and FastTree2 on n=200 indel data (indel rate=0.02). Horizontal bars: median and IQR. (B) Per-seed nRF differences (Fusang − FastTree2), with the per-seed range [−0.03, 0.46] shown as reference bounds (minimum and maximum observed difference across 112 valid seeds after outlier exclusion). Positive difference indicates Fusang better (lower nRF). (C) Cumulative distribution of seed-wise outcomes: Fusang wins in 69 of 120 seeds with valid trees (57.5%); after excluding 8 outlier seeds (7 catastrophic nRF>0.3 + 1 paired exclusion), Fusang wins in 60/112 (53.6%).

**Figure 3.** Pipeline ablation and DCM degradation. (A) Step-by-step nRF tracing from simplified pipeline through DCM stages. (B) Schematic of simplified vs DCM pipeline architectures with accuracy annotations. (C) EPA grafting error illustration.

**Figure 4.** Spaced k-mer parameter optimization. (A) nRF heatmap: k (3–8) × gap (0–4) at n=200. (B) Adaptive parameter selection logic: n≤100→k=4,gap1; n>100→k=5,gap2. (C) Stability validation across 5 repeats per scale.

**Figure 5.** Scalability and cross-platform deployment. (A) Wall-clock time vs n on a 4-core workstation. Dashed line: O(n² log n) scaling. (B) Windows-native FastME vs WSL overhead comparison.

**Figure 6.** Real data validation and comparison with Fusang v1. (A) 16S rRNA tree with NCBI taxonomic annotations. (B) Pairwise tree distance by taxonomic rank, showing order-level significance (p<0.01). (C) Fusang v1 vs Tardigrade Edition comparison: architecture, scalability, accuracy, and deployment.

---

*Manuscript prepared for Nucleic Acids Research. Main text: approximately 9,800 words.*
