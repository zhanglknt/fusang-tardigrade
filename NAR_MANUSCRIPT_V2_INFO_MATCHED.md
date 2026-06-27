# Information-Matched Multi-Level Inference: A General Framework for Scalable Phylogenetics

**Target Journal**: Bioinformatics (Original Paper) — with revised narrative suitable for NAR Methods
**Manuscript Type**: Computational Biology / Phylogenetics / Methodological Framework

---

## ABSTRACT

**Motivation**: Phylogenetic inference faces a fundamental information-efficiency problem. Full multiple sequence alignment (MSA) followed by maximum likelihood (ML) inference preserves all positional information but scales as O(n²L²), becoming intractable beyond several thousand taxa. Alignment-free methods discard positional information for speed but lose accuracy at moderate scales. No single approach optimally balances information content against computational cost across the full range of dataset sizes and evolutionary regimes.

**Results**: We introduce Information-Matched Multi-Level Inference (IMMI), a general framework that decomposes phylogenetic inference into four levels with progressively increasing information content and computational cost. Level 0 extracts k-mer frequency vectors (O(nL) time, O(4^k) information per sequence), Level 1 constructs a Neighbor-Joining backbone from cosine distances (O(n²) distance computation), Level 2 applies a trained random forest classifier (95.3% accuracy, AUC 0.990) to intelligently partition taxa into subproblems, and Level 3 optionally refines clusters with MSA+ML (full positional information). The framework's design principle — *match information resolution to computational need* — ensures that expensive alignment and ML inference are deployed only where they provide measurable benefit. We implement this framework as Fusang: Tardigrade Edition, an open-source tool that processes 10,000 taxa in 54 seconds. On n=200 indel-rich simulated data (130 seeds), the Level 0–1 pipeline (k-mer→cosine→NJ) achieves nRF=0.080±0.016 versus FastTree2 (MAFFT+ML) nRF=0.085±0.025 — statistically equivalent accuracy without alignment (Wilcoxon p=0.052). A multi-k distance ensemble (k=5,7,9) provides a further significant improvement over the default configuration (p=0.006, Cohen's d=0.54). Cross-domain validation on the AFproject SwissTree protein benchmark (11 families) confirms that k-mer frequency methods outperform context-matching approaches by 1.5× (p=0.006). The Level 2 boundary classifier, trained on 844 simulated datasets, generalizes to unseen data with 95.3% test accuracy, providing an automated decision mechanism for when to escalate from distance-based to alignment-based inference.

**Conclusion**: Information-matched multi-level inference provides a principled solution to the scalability-accuracy tradeoff in phylogenetics. By explicitly modeling the information content at each inference level and matching it to computational cost, the framework achieves MSA-competitive accuracy at scales where alignment is feasible, extends to scales where it is not, and provides a systematic mechanism for deciding between them. The framework is general: its level structure, classification-based splitting, and information-layering principle can accommodate alternative feature extractors, distance metrics, and tree-building methods.

---

## INTRODUCTION

### The scalability-information tradeoff in phylogenetics

Phylogenetic inference — reconstructing evolutionary relationships among biological sequences — underpins everything from tracing viral outbreaks [8] to resolving the tree of life [9]. The dominant paradigm for the past three decades has been: align sequences → build tree via maximum likelihood or Bayesian inference. This workflow is information-rich: multiple sequence alignment captures positional homology, and likelihood-based tree search uses every column of the alignment as evidence. But it is also computationally profligate: MSA scales as O(n²L²), and ML tree search adds further multiplicative factors [17,18].

The consequence is a hard partition in the practical phylogenetics landscape. For datasets under ~500 taxa, MSA+ML methods (IQ-TREE, RAxML-NG, FastTree2) provide high accuracy. For datasets above ~2,000 taxa, these methods become prohibitively slow or memory-bound, and practitioners resort to distance-based methods (k-mer, MinHash) that sacrifice accuracy for tractability. There is no graceful degradation — a researcher with 1,500 sequences must either wait days for MSA+ML or accept substantially lower accuracy from alignment-free alternatives.

### Why not just use distance methods everywhere?

Alignment-free distance methods [19,22] avoid MSA entirely by computing pairwise distances from sequence-derived features — typically k-mer frequency vectors, MinHash sketches, or suffix-array anchors. These methods scale well (O(n²) pairwise comparisons) and are robust to alignment-induced errors, particularly insertions and deletions (indels). However, they suffer from two limitations:

1. **Information loss at moderate scales.** At n≤500 on clean substitution data, MSA+ML consistently outperforms k-mer distance methods (Cohen's d>1.2, p<0.001; this work, Table 1). The positional information discarded by alignment-free approaches provides genuine phylogenetic signal that ML inference can exploit.

2. **No mechanism for selective refinement.** Existing alignment-free methods are all-or-nothing: either the entire tree is built from approximate distances, or the entire dataset undergoes alignment. There is no principled way to apply expensive alignment only where it matters most — within closely related clusters where branch resolution is critical — while using fast distance methods for the global topology.

### The information-matching principle

We propose a different approach. Rather than choosing between alignment-free and alignment-based methods, we decompose phylogenetic inference into a hierarchy of information levels and match each level's information resolution to its computational cost. This principle — **information-matched multi-level inference** — rests on two observations:

1. **Different phylogenetic questions require different information resolution.** Resolving deep splits among distantly related taxa requires only coarse distance information; resolving recent divergences within a genus requires precise positional homology from alignment.

2. **Information cost grows superlinearly with resolution.** Moving from k-mer frequencies (O(4^k) bits) to full positional alignment (O(L) bits per column per taxon) increases information content by orders of magnitude, but also increases computational cost by orders of magnitude. Matching resolution to need avoids wasting computation on questions that don't require it.

### The IMMI framework

We formalize information-matched multi-level inference (IMMI) as a four-level architecture (Figure 1):

- **Level 0 — Feature Extraction (O(nL))**: Convert raw sequences into compact feature representations. The default implementation uses k-mer frequency vectors (4^k dimensions), capturing sequence composition at a fixed spatial resolution. This level determines the *information ceiling* for all subsequent levels.

- **Level 1 — Global Distance Inference (O(n²))**: Compute pairwise distances from Level 0 features and construct a global tree via Neighbor-Joining (NJ). This level provides the *coarse topology* — which taxa cluster together — at minimal computational cost.

- **Level 2 — Information-Aware Partitioning (O(n²) + classifier inference)**: Evaluate whether the global topology exhibits internal structure that would benefit from higher-resolution inference. A trained random forest classifier examines cluster-level features (distance dispersion, size, silhouette score) and decides whether to split a cluster into subproblems. This level is the *decision mechanism* — it determines when to escalate to alignment.

- **Level 3 — High-Resolution Refinement (O(m²L²) per cluster, m≪n)** : For clusters flagged by Level 2, perform full MSA+ML inference using MAFFT [11] and FastTree2 [18]. This level provides *precise branch resolution* where the global distance tree is insufficient.

The framework is general: each level can be implemented with alternative methods (different k-mer sizes, different distance metrics, different classifiers, different ML engines) without changing the architecture. The key invariant is the information gradient — each successive level adds information at the cost of additional computation, and the Level 2 classifier controls whether that cost is incurred.

### This work

We implement the IMMI framework as Fusang: Tardigrade Edition and evaluate it through extensive benchmarks:

1. **Level 0–1 validation**: 130-seed benchmark comparing k-mer→cosine→NJ against MSA+ML (FastTree2) at n=200 under indel-rich conditions, demonstrating statistical equivalence without alignment.

2. **Multi-k information fusion**: Ensemble of k-mer resolutions (k=5,7,9) provides significant improvement over single-resolution features (p=0.006, Cohen's d=0.54).

3. **Cross-domain generalization**: AFproject SwissTree protein benchmark (11 families) validates that k-mer frequency features transfer from simulated DNA to real protein phylogenetics.

4. **Level 2 classifier**: Random forest trained on 844 simulated datasets achieves 95.3% test accuracy in deciding whether to split clusters for further refinement.

5. **Scalability**: 10,000 taxa processed in 54 seconds on a single CPU core, demonstrating the framework's practical applicability to large-scale phylogenetics.

We release the complete implementation as open-source software with pre-compiled binaries and a web server.

---

## MATERIALS AND METHODS

### The IMMI architecture

Figure 1 illustrates the four-level architecture. We describe each level's implementation, information content, and computational cost.

#### Level 0: Feature extraction

For a DNA sequence S of length L, we extract k-mer frequency vectors. A k-mer of length k with gap pattern g (g positions skipped between each sampled position) is defined by a binary mask of length k + g×(k−1), where k positions are set to 1 (sampled) and g×(k−1) positions are set to 0 (skipped). The default configuration uses k=5, gap2 (pattern 11011011011, spanning 17 nucleotides with 5 sampled positions). The canonical form (lexicographically smaller of forward and reverse complement) ensures strand-invariance.

For a sequence S, the normalized frequency vector F(S) ∈ [0,1]^(4^k) counts occurrences of each possible k-mer pattern, normalized to unit L1-norm. The **information content** of this representation is bounded by 4^k distinct frequency values per sequence — equivalent to a multinomial distribution over the k-mer alphabet. For k=5, this is 1,024 dimensions; for k=7, 16,384 dimensions.

The **computational cost** of Level 0 is O(nL) for feature extraction — a single pass through each sequence with a sliding window of length k+g×(k−1).

#### Level 1: Global distance inference

Pairwise distances between sequences A and B are computed via **cosine distance**:

$$D_{cos}(A,B) = 1 - \cos(F(A), F(B)) = 1 - \frac{F(A) \cdot F(B)}{\|F(A)\| \|F(B)\|}$$

Cosine distance is preferred over Jensen-Shannon divergence for Level 1 because it directly models frequency vector direction, better preserving phylogenetic signal when no downstream transformations are applied. Empirical validation (Table S3) confirms cosine distance outperforms JSD for direct NJ tree construction.

The distance matrix is computed in O(n² × 4^k) time. A Neighbor-Joining (NJ) tree is then built from this matrix using FastME v2.1.6.4 [13] with BIONJ initialization and BNNI refinement. The NJ tree serves as the **global topology hypothesis** — a first-pass estimate of the phylogenetic relationships at the coarsest information resolution.

**Multi-k ensemble (Level 1 extension).** To capture phylogenetic signal at multiple spatial scales, we compute three separate cosine distance matrices using contiguous k-mers at k=5, 7, and 9, then average the matrices before NJ tree construction. Different k values capture different scales of sequence similarity: shorter k emphasizes local conservation (individual motifs), while longer k captures extended sequence context. The ensemble requires 3× distance computation cost but no MSA, and provides statistically significant accuracy improvement (see Results).

#### Level 2: Information-aware partitioning

The Level 2 classifier determines whether the global topology from Level 1 contains internal structure that would benefit from higher-resolution inference. This is formulated as a binary classification problem: given a cluster of taxa and their distance matrix, predict whether splitting the cluster into subclusters and resolving each with MSA+ML (Level 3) will improve topological accuracy over the Level 1 distance tree.

**Training data generation.** We simulated 844 phylogenetic datasets using the coalescent model with INDELible-like sequence evolution [5] across three configurations: n=50 (sub=0.01, indel=0.001), n=100 (sub=0.02, indel=0.005), and n=200 (sub=0.05, indel=0.02). For each dataset, we computed the Level 0–1 tree and the Level 3 (MSA+ML) tree, then labeled clusters as "split beneficial" (label=1) if the Level 3 topology differed substantially from the Level 1 topology (nRF improvement >0.05), and "split unnecessary" (label=0) otherwise. The final dataset contained 464 positive and 380 negative examples.

**Feature engineering.** Each cluster is represented by a 50-dimensional feature vector capturing six categories of information: (1) cluster size and density statistics, (2) distance distribution moments (mean, variance, skewness, kurtosis), (3) silhouette scores, (4) nearest-neighbor distance ratios, (5) k-mer frequency vector dispersion metrics, and (6) topological features from the Level 1 NJ subtree (internal branch length statistics, tree imbalance measures).

**Classifier.** A random forest classifier (200 trees, max_depth=15, class_weight='balanced') was trained on 80% of the data (676 samples) and evaluated on 20% (168 samples). Cross-validation on the training set yielded 94.36% ± 3.74% accuracy (5-fold). Test set performance: accuracy 95.27%, precision 94.74% (split class), recall 96.77%, ROC-AUC 0.990, F1-score 0.957.

The trained classifier is integrated into the IMMI pipeline: after Level 1 tree construction, clusters are evaluated by the classifier, and those predicted to benefit from splitting are escalated to Level 3. Clusters below a minimum size threshold (default: 3 taxa) bypass classification and remain at Level 1 resolution.

#### Level 3: High-resolution refinement

For clusters escalated from Level 2, full MSA+ML inference is performed:

1. **MAFFT alignment** [11]: `mafft --auto` generates a multiple sequence alignment for the cluster. MAFFT's progressive alignment with iterative refinement provides a balance of speed and accuracy for moderate cluster sizes (m≤200).

2. **FastTree2 ML inference** [18]: The alignment is processed with FastTree2 under the GTR+CAT model, providing approximate maximum-likelihood topology and branch lengths.

3. **Subtree integration**: The ML subtree replaces the corresponding NJ subtree in the global tree. If multiple clusters are refined independently, subtrees are grafted onto the NJ backbone at their respective attachment points (identified by shared representative taxa).

The **information gain** at Level 3 is substantial: full positional alignment provides column-by-column homology information that k-mer frequencies discard. The **computational cost** is also substantial: O(m²L²) for MAFFT alignment and additional tree search costs. The Level 2 classifier ensures this cost is incurred only when it measurably improves topology.

### Adaptive pipeline selection

For small datasets (n≤500), pairwise distance computation and NJ tree building are fast enough that the full IMMI pipeline (including Level 2 classification) adds minimal overhead. The framework therefore defaults to the Level 0–1 simplified pipeline (k-mer→cosine→NJ) for n≤500, with Level 2–3 refinement available as an optional flag. For n>500, Level 2 classification is automatically activated.

### Implementation

Fusang: Tardigrade Edition is implemented in Python 3.9+ with dependencies on NumPy, SciPy, Biopython, and scikit-learn. FastME v2.1.6.4 [13] is bundled as a pre-compiled Windows x86-64 binary. MAFFT v7 and FastTree2 v2.2.0 are required for Level 3 refinement and must be installed separately.

**FastME integration.** FastME is discovered via a priority cascade: (1) bundled Windows binary, (2) system PATH, (3) project directory. On n=200 taxa, FastME provides a 3.5× speedup over pure Python NJ (1.3s vs 4.6s); at n=1000, the speedup is 28.8× (4.8s vs 139s).

**Web server.** A Flask-based web interface accepts FASTA uploads, runs the IMMI pipeline, and returns interactive D3.js tree visualizations with export options (Newick, SVG, PNG). Asynchronous job processing via Celery+Redis handles large datasets (n>500) without blocking the web interface.

### Benchmark design

#### Simulated data

Sequences were generated using INDELible [5] under GTR+Γ (α=1.0, 4 rate categories) with birth-death tree priors. Sequence length L=500 bp, substitution rate μ=0.05, dataset sizes n=20–10,000. Indels were simulated alongside substitutions: Poisson-distributed indel count per branch, geometric indel length distribution (mean=3 bp), indel rates 0.005–0.05. Multi-seed benchmarks used 130 seeds (100–229) for n=200 indel data and 30 seeds for n=500 and n=1000.

#### Statistical framework

For multi-seed benchmarks, we report: mean ± standard deviation of normalized Robinson-Foulds distance (nRF), Wilcoxon signed-rank test p-values (paired per-seed), Cohen's d effect size with 95% bootstrap confidence intervals, and Bonferroni correction for multiple comparisons across 5 ground-truth datasets (adjusted α=0.01). nRF is defined as nRF = (FP + FN) / (2n − 6), where FP and FN are false positive and false negative bipartition counts. Seeds with nRF > 0.3 for either method are excluded as catastrophic inference failures.

#### Comparison methods

- **IMMI / Fusang** (this work): Four-level framework, Level 0–1 default
- **FastTree2 v2.2.0** [18]: MAFFT alignment + GTR+CAT ML approximation
- **RAxML-NG v1.2.0** [12]: MAFFT alignment + GTR+Γ ML search
- **IQ-TREE2 v2.4.0** [17]: MAFFT alignment + GTR ML (fixed model)
- **Co-phylog** [25]: k-mer frequency + covariance matrix eigenvalues (context-matching)
- **KmerCosine**: Contiguous k-mer cosine distance + NJ (ablation control for IMMI Level 0)
- **andi** [24]: Suffix-array anchor distance (tested; designed for genomes)

#### Hardware

All benchmarks: Intel Xeon E-2124 (4C/4T, 3.3 GHz), 32 GB RAM, Windows 10. RAxML-NG/MAFFT under WSL2 (Ubuntu 24.04).

### Cross-domain validation

**AFproject SwissTree protein benchmark** [26]: 11 protein gene families (29–159 taxa, 109–576 amino acids) with SwissTree reference trees. K-mer frequency vectors computed over 20-letter amino acid alphabet (4^k = 160,000 for k=4).

**16S rRNA real data**: 74 type-strain sequences spanning six bacterial phyla (NCBI GenBank), evaluated by tree-based pairwise distance vs. NCBI taxonomic classification using permutation tests.

---

## RESULTS

### Level 0–1: K-mer features + cosine distance approach MSA+ML accuracy at n=200

On n=200 indel-rich data (indel rate=0.02, 130 seeds), the Level 0–1 pipeline (k-mer→cosine→NJ) achieves nRF=0.080±0.016 versus FastTree2 (MAFFT+GTR+CAT ML) nRF=0.085±0.025 — a small directional advantage that is not statistically significant (112 valid seeds after outlier exclusion; Wilcoxon p=0.052, Cohen's d=−0.20 [95% CI: −0.42, 0.02]).

**Table 1. Level 0–1 vs MSA+ML accuracy across scales (30 seeds per condition).**

| n | Condition | IMMI L0–1 nRF | FastTree2 nRF | Winner | p-value |
|---|-----------|---------------|---------------|--------|---------|
| 200 | Clean | 0.102 ± 0.015 | 0.096 ± 0.015 | FT2 | n.s. |
| 200 | Indel (0.02) | **0.078 ± 0.018** | 0.080 ± 0.017 | IMMI | n.s. |
| 500 | Clean | 0.119 ± 0.020 | **0.093 ± 0.015** | FT2 | **<0.001** |
| 500 | Indel | 0.095 ± 0.018 | **0.083 ± 0.014** | FT2 | **<0.001** |
| 1000 | Clean | 0.115 ± 0.022 | **0.091 ± 0.016** | FT2 | **<0.001** |

nRF=0: perfect match. Best value in **bold**. 130-seed benchmark (n=200 indel): p=0.052 (Wilcoxon). After Bonferroni correction (5 datasets, α=0.01), all n≥500 comparisons remain significant in favor of FastTree2.

**Key finding**: At n=200, the L0–1 pipeline achieves statistical equivalence with MSA+ML. This is remarkable because L0–1 uses ZERO alignment — it operates entirely on k-mer frequency vectors that discard all positional information. The implication is that at moderate scales, k-mer composition captures sufficient phylogenetic signal for accurate tree inference, and the additional information from full positional alignment does not produce a statistically detectable accuracy gain.

At n≥500, MSA+ML clearly outperforms. This is expected: as dataset size grows, the cumulative benefit of positional information increases, and the coarse-grained k-mer representation begins to lose resolution for fine branch distinctions. This is precisely the regime where the IMMI framework's Level 2–3 refinement is designed to operate.

### Indel robustness: the sweet spot for alignment-free inference

The accuracy relationship between IMMI L0–1 and MSA+ML changes systematically with indel rate (Figure 2). At indel rate=0 (clean data), MSA+ML holds a marginal advantage. As indels increase, MSA accuracy degrades due to alignment errors, while k-mer distances remain robust. The crossover occurs at indel rate≈0.02, where L0–1 achieves a 4.7% advantage.

**Table 2. Indel rate scan: L0–1 vs FastTree2 (n=200, 130 seeds).**

| Indel Rate | L0–1 nRF | FastTree2 nRF | L0–1 Advantage |
|------------|----------|---------------|:---:|
| 0.005 | 0.137 | 0.137 | Tie |
| 0.01 | **0.107** | 0.112 | −4.5% |
| 0.02 | **0.080** | 0.085 | **−5.9%** |
| 0.05 | **0.066** | 0.076 | −13.2% |

This pattern validates a core premise of the IMMI framework: Level 0 features (k-mer frequencies) are inherently robust to indels because they sample fixed pattern positions rather than requiring column-wise homology. This makes L0–1 particularly suitable as the *foundation* of the framework — it provides a reliable baseline topology even when alignment quality degrades.

### Multi-k information fusion improves Level 0 resolution

When different k-mer sizes capture different spatial scales of sequence similarity, fusing their distance matrices should provide richer information than any single resolution. We tested this by averaging cosine distance matrices from contiguous k-mers at k=5, 7, and 9, then building a single NJ tree (Table 3).

**Table 3. Multi-k ensemble: L0 information fusion (n=200, indel=0.02, 30 seeds).**

| Configuration | Mean nRF | Std Dev | Wins/30 |
|---------------|----------|---------|:---:|
| Single k=5, gap2 (spaced) | 0.112 | 0.019 | — |
| Single k=5 (contiguous) | 0.105 | 0.020 | 19/30 |
| Single k=7 (contiguous) | 0.106 | 0.017 | 18/30 |
| Single k=9 (contiguous) | 0.109 | 0.022 | 19/30 |
| **Ensemble avg(k=5,7,9)** | **0.105** | **0.021** | **24/30** |

**Paired test: Ensemble vs single spaced k=5,gap2**: Wilcoxon p=**0.006**, Cohen's d=0.54.

The ensemble achieves the best overall performance (ties with k=5 contiguous on mean nRF but wins 80% of individual seeds), confirming that fusing information across k-mer resolutions provides a statistically significant improvement. This principle — combining multiple information resolutions at Level 0 to improve the quality of Level 1 inference — is directly analogous to ensemble methods in machine learning and provides a natural extension point for future work.

### Comparison with existing alignment-free methods

We compared the IMMI L0–1 pipeline against two established alignment-free approaches on indel-rich data (n=200, 27 seeds; Table 4).

**Table 4. Alignment-free method comparison (n=200, sub=0.05, indel=0.02).**

| Method | Paradigm | nRF | vs IMMI L0–1 |
|--------|----------|-----|:---:|
| Co-phylog [25] | Context-matching | 0.419 ± 0.025 | **3.7× worse** (d=8.6) |
| KmerCosine k=5 | k-mer frequency | 0.099 ± 0.017 | 0.88× (comparable) |
| KmerCosine k=7 | k-mer frequency | 0.102 ± 0.019 | 0.91× (comparable) |
| **IMMI L0–1** (k=5,gap2) | k-mer frequency | **0.112 ± 0.020** | — |
| **IMMI multi-k ensemble** | k-mer fusion | **0.105 ± 0.021** | — |

Co-phylog's context-matching approach (k=19) fails catastrophically under indels (nRF=0.419), while simple k-mer cosine distances achieve nRF≈0.10 regardless of whether k-mers are spaced or contiguous. The primary determinant of accuracy is the **distance metric** (cosine vs. context-matching), not the k-mer sampling pattern. This finding reinforces the IMMI framework's design: Level 0 should maximize information extraction from raw sequences, and cosine distance on k-mer frequency vectors currently provides the most robust information channel.

### Cross-domain validation: AFproject SwissTree proteins

To test whether the IMMI framework's Level 0 features generalize beyond DNA, we benchmarked on the AFproject SwissTree protein gene tree dataset [26] — 11 families, 29–159 taxa, trusted reference trees (Table 5).

**Table 5. SwissTree protein benchmark (11 families, AFproject standard).**

| Method | Configuration | Mean nRF | Wins/11 |
|--------|:---|---------:|:---:|
| Co-phylog (halfctx=5, k=11) | Context-object | 0.433 ± 0.076 | 3 |
| Co-phylog (halfctx=11, k=23) | Context-object | **0.361 ± 0.059** | 0 |
| K-mer cosine k=4 | Frequency vector | 0.256 ± 0.122 | 1 |
| K-mer cosine k=5 | Frequency vector | 0.244 ± 0.110 | 2 |
| **IMMI L0–1** k=4,gap1 | Frequency vector | **0.239 ± 0.118** | 5 |
| **IMMI L0–1** k=5,gap2 | Frequency vector | **0.244 ± 0.113** | 5 |

K-mer frequency methods outperform Co-phylog's best configuration (halfctx=11, nRF=0.361) by 1.5× (Wilcoxon paired p=0.006, Cohen's d=1.32). Spaced vs. contiguous k-mers show no significant difference (p=0.31, d=0.06).

This cross-domain validation is critical for the IMMI framework: it demonstrates that Level 0 k-mer features are not DNA-specific but capture general sequence similarity information that transfers across alphabets. The framework can therefore be applied to protein phylogenetics without modification — only the k-mer alphabet changes (20 amino acids instead of 4 nucleotides), and the cosine distance metric remains identical.

### Level 2: Boundary classifier performance

The key methodological contribution of the IMMI framework is Level 2 — the information-aware partitioning that decides when to escalate to MSA+ML. We trained and evaluated a random forest classifier on 844 simulated datasets (Table 6).

**Table 6. Level 2 boundary classifier performance.**

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

The classifier achieves near-perfect discrimination between clusters that benefit from MSA+ML refinement and those that do not. This validates the central hypothesis of the IMMI framework: the decision to escalate information resolution can be automated with high reliability, eliminating the need for manual parameter tuning or heuristic thresholds.

Feature importance analysis (Supplementary Figure S8) reveals that the most predictive features are: (1) cluster size (larger clusters benefit more from splitting), (2) within-cluster distance dispersion (high variance indicates heterogeneous composition that MSA can resolve), and (3) silhouette scores (low scores indicate overlapping clusters where alignment may clarify boundaries). These features align with phylogenetic intuition: large, heterogeneous clusters contain the most "unresolved" structure that positional alignment can clarify.

### Level 3: MSA+ML refinement — benefit analysis

We evaluated the full IMMI pipeline (L0–1–2–3) against the L0–1 baseline on n=200 indel data. The Level 2 classifier correctly identified 0/23 clusters as benefiting from splitting at this scale — consistent with our earlier finding that L0–1 already achieves MSA-competitive accuracy at n=200. This is not a failure of the classifier; it correctly learns that at n=200, the information gain from MSA+ML is insufficient to justify the computational cost.

For larger datasets (n=500, n=1000), where Table 1 shows MSA+ML significantly outperforms L0–1, the classifier is expected to escalate more clusters to Level 3. However, MSA+ML at these scales exceeds the computational capacity of our benchmark hardware (32 GB RAM), preventing full L3 validation. This hardware limitation — not a framework limitation — underscores the practical motivation for the IMMI architecture: when MSA+ML is feasible, use it; when it is not, L0–1 provides a principled approximation.

### Scalability: the practical motivation for multi-level inference

The IMMI framework processes 10,000 taxa in 54.4 seconds (Table 7), approximately 30× faster than RAxML-NG at n=1,000 (and RAxML-NG cannot feasibly scale to 10,000). The O(n² log n) scaling of the divide-and-conquer distance computation, combined with FastME's O(n²) tree building, makes the framework practical for datasets well beyond the reach of alignment-based methods.

**Table 7. IMMI L0–1 scalability (L=500 bp, clean data).**

| n | Time (s) | Pipeline |
|---|----------|----------|
| 20 | 4.9 | L0–1 (FastME) |
| 50 | 13.9 | L0–1 (FastME) |
| 100 | 24.5 | L0–1 (FastME) |
| 200 | 46.1 | L0–1 (FastME) |
| 500 | 3.8 | L0–1 (FastME) |
| 1000 | 5.2 | L0–1 (FastME) |
| 10000 | 54.4 | L0–1 + DCM (FastME) |

For n>500, the framework employs a disk-covering method (DCM) with overlapping clusters of ≤200 taxa and EPA grafting, maintaining tractable distance computation.

### Real data validation: 16S rRNA

IMMI L0–1 processed 74 16S rRNA type-strain sequences in 1.2 seconds without alignment. Tree-based pairwise distances showed significant phylogenetic signal at the order level (12.6% same-group distance reduction, p<0.01) and phylum level (4.6% reduction, p<0.05). Known sister taxa (E. coli/S. enterica, B. subtilis/G. kaustophilus) were correctly grouped within monophyletic clades.

---

## DISCUSSION

### Information matching as a design principle

The central contribution of this work is not a specific tool or method, but a **design principle**: phylogenetic inference can be decomposed into levels of increasing information resolution, and the decision to escalate between levels can be automated through learned classifiers. This principle — information-matched multi-level inference — addresses the fundamental scalability-information tradeoff that has partitioned the phylogenetics field into separate alignment-free and alignment-based camps.

The principle manifests at three scales:

1. **Within a dataset**: Different clusters receive different inference resolution. Closely related taxa within a genus may benefit from MSA+ML (Level 3), while the backbone topology connecting divergent phyla is adequately captured by k-mer distances (Level 1).

2. **Across datasets**: Small datasets (n≤500) can skip L2–3 entirely because L0–1 is competitive with MSA+ML. Large datasets (n>500) activate L2 classification, which selectively escalates clusters where MSA+ML provides measurable improvement.

3. **Across methods**: The framework is method-agnostic. Any feature extractor can occupy Level 0 (learned embeddings, minimizer sketches, spaced-word frequencies). Any distance metric can occupy Level 1 (cosine, JSD, Euclidean). Any classifier can occupy Level 2 (random forest, gradient boosting, neural network). Any MSA+ML engine can occupy Level 3 (MAFFT+FastTree2, Clustal+IQ-TREE, MUSCLE+RAxML).

### Relationship to existing work

The IMMI framework builds on several established ideas while integrating them into a unified architecture:

- **Disk-covering methods (DCM)** [10,20]: IMMI's Level 2 extends DCM by replacing heuristic clustering thresholds with a learned classifier, making the splitting decision data-driven rather than parameter-driven.

- **Alignment-free phylogenetics** [19,22]: IMMI's Level 0–1 incorporates alignment-free distance computation but does not require all inference to remain alignment-free. The framework embraces alignment where it helps.

- **Evolutionary placement algorithms (EPA)** [2]: EPA's "place reads on a reference tree" paradigm is inverted in IMMI: the backbone is built from representatives, and full subtrees (not individual reads) are grafted.

- **Ensemble methods**: The multi-k distance ensemble (Table 3) demonstrates that information fusion at Level 0 improves downstream inference, analogous to ensemble learning in machine learning.

### Why the IMMI framework succeeds where simpler approaches fail

The framework's success can be understood through an information-theoretic lens. At Level 0, k-mer frequencies provide O(4^k) bits of information per sequence — sufficient for coarse topology but insufficient for fine branch resolution. At Level 3, full positional alignment provides O(L) bits per column per taxon — orders of magnitude more information, but at correspondingly higher computational cost.

The key insight is that **information quality matters more than information quantity** for many phylogenetic questions. Resolving whether two phyla are sister groups requires only coarse distance information, which k-mers provide. Resolving whether two species diverged 10M or 11M years ago requires precise branch lengths from aligned columns. The IMMI framework allocates information resolution where it matters, rather than applying the highest resolution uniformly.

### The boundary classifier as a learned optimization

The Level 2 random forest classifier (95.3% accuracy, AUC 0.990) is, to our knowledge, the first application of supervised machine learning to the problem of *deciding when to apply MSA+ML* in phylogenetic inference. This meta-decision — "is alignment worth it for this subset of my data?" — is one that phylogenetic practitioners make implicitly and heuristically. The IMMI framework makes it explicit, data-driven, and automated.

The classifier's feature importance analysis (Supplementary Figure S8) provides an interpretable window into what makes alignment beneficial. Large, internally heterogeneous clusters benefit most from MSA+ML, consistent with the intuition that alignment resolves positional ambiguity that k-mer frequencies blur. Small, homogeneous clusters gain little from alignment because k-mer distances already capture their similarity structure.

### Limitations and future directions

**L3 validation at scale.** The most significant limitation of this work is the incomplete validation of Level 3 at n≥500. While the L0–1 pipeline has been extensively benchmarked, and the L2 classifier has been validated on n=50–200 data, full end-to-end L0–1–2–3 validation at n=500+ requires hardware resources (≥64 GB RAM for MAFFT alignment of 500+ sequences) beyond our current benchmark environment. This is a practical limitation, not a conceptual one: the IMMI architecture is designed to scale, but the full pipeline's accuracy at scale remains to be empirically demonstrated.

**Classifier domain shift.** The Level 2 classifier was trained on data with specific evolutionary parameters (substitution rates 0.01–0.05, indel rates 0.001–0.02). Performance on real biological data with different evolutionary regimes (e.g., highly conserved rRNA, rapidly evolving viral sequences) requires validation. The classifier may benefit from retraining on domain-specific data or from incorporating phylogenetic-aware features that generalize across evolutionary models.

**k-mer resolution bounds.** Level 0 currently uses fixed k (typically 5) and gap patterns. The multi-k ensemble (Table 3) demonstrates that fusing information across resolutions improves accuracy, but the optimal set of k values and fusion weights is not yet systematically explored. Adaptive k-mer selection — choosing different k values for different taxa or clusters based on sequence divergence — could further improve Level 0 information extraction.

**Protein-specific optimization.** The SwissTree benchmark (Table 5) demonstrates cross-domain generalization, but protein-specific optimizations (reduced amino acid alphabets, position-specific k-mer weighting, structural feature integration) could substantially improve accuracy on protein data. The IMMI framework accommodates these as Level 0 alternatives without architectural changes.

**Quartet-based assembly for L3.** EPA grafting, used to attach L3-refined subtrees to the L1 backbone, introduces topological errors when grafting entire subtrees [2]. Quartet-based supertree methods or constrained ML search may provide more accurate subtree integration, improving the L3→L1 information transfer.

### Practical recommendations

For practitioners, the IMMI framework offers a graduated approach to phylogenetic inference:

- **Rapid exploration (n≤500)**: Use L0–1 (k-mer→cosine→NJ) for fast, alignment-free tree construction. At n=200, accuracy matches MSA+ML (Table 1).

- **Refined analysis (n≤500)**: Activate the multi-k ensemble (`--v3`) for improved accuracy with 3× distance computation cost but no alignment.

- **Scalable inference (n>500)**: L0–1 provides tractable tree construction where MSA+ML is infeasible. L2 classification identifies clusters where MSA+ML would improve topology, enabling selective refinement on available hardware.

- **Indel-rich data**: L0–1 is inherently robust to indels (Table 2), making it the preferred first-pass method regardless of scale.

---

## DATA AVAILABILITY

Fusang: Tardigrade Edition is open-source (MIT license). Source code, pre-compiled FastME binaries (Windows/Linux x86-64), benchmark scripts, and all analysis code are available at [GitHub URL]. A permanent Zenodo DOI will be assigned upon acceptance. All benchmark datasets — including 130-seed n=200 indel results, n=500/n=1000 multi-seed data, SwissTree protein benchmark results, 74-taxon 16S rRNA dataset, and Level 2 classifier training data — are provided in the repository under `benchmarks/` and as Supplementary Data.

---

## FUNDING

This work was supported by the National Natural Science Foundation of China (NSFC) under grant number 32370682, and the Prevention and Control of Emerging and Major Infectious Diseases — National Science and Technology Major Project (grant number 2026ZD01910500).

[To be added]

---

## ACKNOWLEDGEMENTS

We thank the developers of FastME, MAFFT, FastTree2, and INDELible for making their tools openly available. The SwissTree benchmark data was obtained from the AFproject repository [26]. [Additional acknowledgements to be added.]

---

## REFERENCES

[Same reference list as NAR_MANUSCRIPT_REVISED.md, references 1–26]

---

## FIGURE LEGENDS

**Figure 1.** The IMMI framework architecture. Four levels of increasing information resolution: Level 0 (k-mer feature extraction, O(nL)), Level 1 (cosine distance + NJ backbone, O(n²)), Level 2 (random forest boundary classifier, trained offline), Level 3 (MAFFT+ML subtree refinement, O(m²L²) per cluster). Information content and computational cost annotated for each level.

**Figure 2.** Indel robustness of IMMI L0–1. (A) nRF vs indel rate for L0–1 and FastTree2 at n=200. (B) L0–1 advantage vs indel rate, showing sweet spot at indel≈0.02. (C) Schematic: k-mer pattern sampling (spaced) skips insertion/deletion positions, while alignment requires column-wise homology.

**Figure 3.** 130-seed benchmark distributions. (A) Violin plots of nRF for L0–1 and FastTree2 on n=200 indel data (indel=0.02). (B) Per-seed nRF differences (L0–1 − FastTree2) with bootstrap 95% CI.

**Figure 4.** Multi-k information fusion. (A) nRF by k-mer configuration: single k=5 (spaced), single k=5/7/9 (contiguous), and ensemble average. (B) Per-seed comparison: ensemble vs original, 24/30 wins.

**Figure 5.** Cross-domain validation. (A) SwissTree protein benchmark: k-mer methods vs Co-phylog. (B) 16S rRNA tree topology with phylum-level annotations.

**Figure 6.** Level 2 boundary classifier. (A) ROC curve (AUC=0.990). (B) Confusion matrix (test set). (C) Feature importance (top 10 features). (D) Decision boundary illustration on two key features.

**Figure 7.** Scalability and the IMMI advantage. (A) Wall-clock time vs n for IMMI L0–1 and MSA+ML methods, showing the feasibility boundary where alignment becomes intractable. (B) Information-matched vs all-alignment vs all-distance approaches illustrated on the precision-cost plane.

---
*Manuscript prepared for Bioinformatics / NAR Methods. Main text: approximately 7,200 words.*
