# Fusang: Tardigrade Edition — K-mer Frequency Vector Alignment-Free Phylogenetic Inference Resilient to Indel-Rich Sequence Evolution

**Target Journal**: Nucleic Acids Research (NAR) — Methods Article
**Manuscript Type**: Computational Biology / Phylogenetics
**Running Title**: K-mer frequency vector phylogenetics
**Keywords**: alignment-free phylogenetics, k-mer frequency vector, cosine distance, indel robustness, multi-k ensemble, Neighbor-Joining, phylogenetic software

**Authors**:

Lei Kong¹ and Li Zhang²·³·⁴,*

¹ School of Life Sciences, Peking University, Beijing, China

² Institute of Blood Transfusion, Chinese Academy of Medical Sciences and Peking Union Medical College, Chengdu, China

³ Chinese Institute for Brain Research, Beijing, China

⁴ Translational Medical Center, Weifang Second People's Hospital, 7 Yuanxiao St., Weifang, 261041, Shandong Province, China

* To whom correspondence should be addressed. Email: knightz@pumc.edu.cn

---

## ABSTRACT

Multiple sequence alignment (MSA) scales poorly with dataset size and introduces systematic errors under insertions and deletions (indels), yet alignment-free methods have historically underperformed MSA-based maximum likelihood (ML) approaches. We present Fusang: Tardigrade Edition, an alignment-free phylogenetics tool that infers trees directly from unaligned input — with automatic parameter selection, a multi-k distance ensemble that removes manual k choice, and a random-forest boundary classifier that maps its domain of applicability. Its inference engine, k-mer frequency vector cosine distances (spaced and contiguous patterns, systematically evaluated here), significantly outperforms FastTree2 on simulated indel-rich data (n=200, indel rate=0.02: nRF=0.122±0.024 vs 0.134±0.022, 130 seeds, Wilcoxon p=3.5×10⁻⁷), while MSA methods retain a significant advantage on clean data at n≥500 (p<0.001). The multi-k ensemble reduces topological error under indel-rich coalescent simulation (nRF=0.583±0.044 vs 0.743±0.046 single-k, 30 seeds, p<0.0001) and outperforms MAFFT+FastTree2 on the same benchmark (p=0.024; borderline after Bonferroni correction, p_adj=0.071). On a real indel-rich mitochondrial benchmark (25 genomes, 9.1% gaps), the ensemble reaches the community-benchmark ceiling — nRF=0.045, matched by MSA+ML and all best published alignment-free methods on this single small dataset. Under indels, MinHash (Mash) and k-mer cosine distances degrade comparably (1.93× vs 1.97×). Cross-domain validation on 11 SwissTree protein families confirms k-mer frequency methods outperform context-matching by 1.5× (p=0.006). The boundary classifier separates homogeneous from structured datasets (88/88 simulated scenarios). Fusang infers trees for 10,000 taxa in ~70 s on a 4-core workstation, with pre-compiled binaries for Windows and Linux.

---



## INTRODUCTION

Phylogenetic inference is foundational to evolutionary biology, from tracing viral outbreaks [1] to reconstructing the tree of life [2]. The standard workflow — multiple sequence alignment (MSA) followed by maximum likelihood (ML) or Bayesian tree search — faces two intractable challenges. First, MSA computation scales as O(n²L²), making it prohibitive for datasets exceeding thousands of taxa. Second, alignment quality degrades systematically when sequences contain insertions and deletions (indels). Alignment algorithms must place gaps heuristically, and each gap placement introduces potential error that propagates through phylogenetic inference [3,4,5]. Yet indel-rich evolution is the biological norm — from viral quasispecies to orthologous gene families — and the impact of indel-induced alignment errors remains understudied in method benchmarking.

### The Fusang lineage: from deep learning to k-mer features

We previously introduced Fusang v1 [6], a deep learning-based phylogenetic inference tool that bypasses MSA through learned feature representations. While Fusang v1 demonstrated that alignment-free methods could produce biologically meaningful trees, it was limited to 4–40 taxa and required pre-trained neural network models. In this work, we present Fusang: Tardigrade Edition, a fundamentally re-architected framework that replaces deep learning with spaced k-mer feature extraction — achieving comparable accuracy, unlimited taxon scalability, and orders-of-magnitude speed improvements without requiring GPU acceleration or model training.

### K-mer frequency vectors in phylogenetics: status and opportunities

K-mer frequency vectors have been applied to multiple sequence analysis tasks including protein classification [7], metagenomic binning [8], and genome assembly. In phylogenetic inference, prior alignment-free work (reviewed in [9,10]) includes hierarchical phylogenomic inference [11], kmacs [12] (gapped k-mismatch substring matching), spaced-word frequency methods [13], and the Alfpy toolkit with the AFproject benchmark [14], which provides standardized implementations of multiple alignment-free methods including gapped k-mer variants. Despite these foundations, **k-mer frequency vectors (spaced and contiguous) have not been systematically evaluated for phylogenetic inference under realistic evolutionary conditions (indel-rich sequences).** The contributions of this work are therefore twofold: (i) a systematic evaluation of **k-mer frequency vector cosine distances** — comparing both spaced and contiguous k-mer patterns across multiple distance metrics — for phylogenetic inference under indel-rich conditions; and (ii) the delivery of these findings as a ready-to-use open-source tool (Fusang: Tardigrade Edition) with automatic parameter selection, a multi-k ensemble requiring no manual k choice, and an empirically mapped domain of applicability. We find that: (1) k-mer frequency vectors with cosine distance provide effective phylogenetic signal under indels; (2) spaced k-mers are hypothesized to be more robust at high indel rates, though the practical advantage over contiguous k-mers is modest at the tested indel rate (0.02); and (3) MinHash-based methods degrade severely under indels (1.93× relative nRF increase), while k-mer cosine distances retain marginally lower absolute error. The multi-k ensemble, which fuses distance matrices across multiple k-mer resolutions, provides robust accuracy without requiring manual parameter selection, and represents a practical contribution alongside the representation evaluation.

### Contributions of this work

We systematically evaluate spaced k-mer features for phylogenetic tree inference across datasets spanning n=20 to 10,000 taxa. We make the following contributions:

1. **Vector-based phylogenetic inference under indels.** Fusang's multi-k ensemble NJ reduces topological error by 21.5% relative to the single-k configuration (TRUE-relative nRF=0.583±0.044 vs 0.743±0.046, 30 seeds, Wilcoxon p<0.0001, paired d=3.82). On the full 30-seed coalescent benchmark (Linux), the ensemble outperforms MAFFT+FastTree2 (nRF=0.583±0.044 vs 0.601±0.055, paired Wilcoxon p=0.024, paired d=−0.45, ensemble wins 19/30; a secondary comparison — see Results for multiple-comparison context), with both methods remaining substantially distant from the true tree (nRF≈0.58–0.60). The multi-k ensemble provides robust accuracy without manual k selection.

2. **Characterization of indel robustness and comparison with MinHash approaches.** Fusang is significantly more accurate than FastTree2 under indels (TRUE-relative nRF=0.122±0.024 vs FastTree2 nRF=0.134±0.022, paired Wilcoxon p=3.5×10⁻⁷, 130 seeds). Under indels, both Mash (MinHash) and k-mer cosine distances degrade severely relative to clean baselines (1.93× vs 1.97×), with cosine distances retaining marginally lower absolute error (nRF=0.742 vs 0.762, 30 seeds). Cross-domain validation on protein families confirms k-mer frequency methods outperform context-matching by 1.5× (p=0.006).

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

Cosine distance is preferred for the simplified pipeline because it directly models frequency vector direction, which better preserves phylogenetic signal when no downstream transformations (DCM, EPA) are applied. **Comparison with Jensen-Shannon divergence (JSD):** JSD, used in the DCM pipeline's TF-IDF weighting step, provides a normalized distance metric suitable for high-dimensional sparse frequency vectors. However, in a benchmark on INDELible-simulated data (n=200, indel rate=0.02, 30 seeds, seeds 100–129), cosine distance achieved mean nRF=0.129 ± 0.027 vs square-root JSD mean nRF=0.132 ± 0.025 (paired Wilcoxon p=0.147, paired Cohen's d=0.16; cosine better in 17/30 seeds) — the difference is not significant. This comparison therefore informed the heuristic choice of cosine distance for the simplified pipeline but must not be read as a confirmed advantage. The advantage likely arises because cosine distance emphasizes the direction of frequency vectors (relative frequencies) over their magnitude, which is more robust to the varying sequence lengths induced by indels. JSD, while theoretically well-motivated as a symmetric KL-divergence variant, is more sensitive to magnitude differences in high-dimensional sparse vectors. Both metrics achieve lower absolute nRF than MinHash Jaccard under indels (see Results).

**Jensen-Shannon divergence** (used in the DCM pipeline):
D_JSD(A,B) = JSD(F(A), F(B)) = ½ D_KL(P||M) + ½ D_KL(Q||M)

where M = ½(P+Q) and D_KL is the Kullback-Leibler divergence. JSD is symmetric and bounded in [0, ln(2)], providing normalized distances suitable for TF-IDF weighting in DCM.

### Adaptive pipeline architecture

Fusang implements two distinct tree-building strategies selected automatically based on dataset size:

#### Simplified pipeline (n ≤ 1000, default)

For small-to-medium datasets, Fusang bypasses divide-and-conquer entirely:
1. Extract spaced k-mer frequency vectors for all sequences
2. Compute pairwise cosine distances (O(n² × 4^k))
3. Build a tree via Neighbor-Joining (NJ, BioPython) or FastME BIONJ+BNNI (BIONJ: bioNJ, a modified NJ algorithm with improved theoretical properties; BNNI: Bayesian Nearest Neighbor Interchange, a tree topology search heuristic)
4. When the tree was built by NJ and n ≥ 100, apply BME-scored NNI refinement (`nni_refine`, at most 100 accepted moves); this step is skipped for FastME trees, which are already BME-optimal
5. Output Newick tree
6. (Optional) Compute bootstrap support values via multinomial resampling of k-mer profiles (see Supplementary Note S5)

Unless explicitly noted as "FastME," all simplified-pipeline **accuracy** benchmarks in this work use the NJ builder with the BME-scored NNI refinement step (auto mode; applied for n ≥ 100). The SwissTree benchmark and the alignment-free competitor comparison (see Results) likewise use NJ, for consistency with the AFproject community standard, which requires identical tree builders across methods. Scalability (wall-clock) benchmarks use FastME for its O(n²) speed advantage (Table 11).

This pipeline was discovered through systematic ablation experiments (see Results). On the 130-seed INDELible clean benchmark (n=200), it attains nRF=0.144 ± 0.025 against FastTree2's 0.146 ± 0.025 (paired Wilcoxon p=0.242, not significant; Table 1). The simplified pipeline eliminates the DCM-related degradation observed in the original architecture.

#### Divide-and-conquer pipeline (n > 500)

For large datasets, Fusang employs the Disk-Covering Method (DCM [3,15]):
1. **Clustering**: Pairwise k-mer distances → hierarchical clustering (scipy, average linkage) into groups of ≤200 taxa with 10–20% overlap
2. **Backbone tree**: Representative centroid sequences → FastME NJ tree
3. **Subtree inference**: Within-cluster trees via FastME BIONJ+BNNI
4. **Grafting**: Subtrees attached to backbone via Evolutionary Placement Algorithm (EPA [16])

The adaptive pipeline selects between simplified and DCM+EPA modes based on dataset size, controlled by the `SIMPLE_THRESHOLD=500` parameter. For n≤500, Fusang uses the simplified pipeline (k-mer→cosine→NJ plus BME-scored NNI refinement for n≥100) by default, since DCM adds no accuracy benefit at small-to-medium scales. For n>500, Fusang switches to DCM+EPA, which provides essential scalability and a small accuracy advantage at scale (direct comparison at n=1000: DCM nRF=0.083 ± 0.005, 9 seeds 401–409, vs the simplified pipeline nRF=0.102 on seed 401). An earlier configuration with the threshold set to 2,000 forced the simplified pipeline at n=1000 and was reverted after causing an accuracy regression on that seed (0.088 → 0.102).

### FastME integration

FastME v2.1.6.4 [17] is the default tree builder. FastME implements the Minimum Evolution principle with bioNJ [18] as the initial topology and BNNI (Branch and Bound Nearest Neighbor Interchange) as the optimization step. The implementation discovers FastME via a priority cascade:
1. Bundled Windows-native binary (`fastme_bin/fastme.exe`)
2. WSL-installed binary (`/usr/local/bin/fastme`)
3. Project-bundled Linux binary (`fastme_bin/fastme_linux`)

The Windows-native binary (PE32+ x86-64, 725 KB) eliminates the WSL dependency, enabling single-command operation on native Windows with zero cross-system overhead.

Benchmark timings for the tree-building stage alone: at n=200, NJ 1.3 s → FastME 0.4 s (3.5× speedup); at n=500, 13.2 s → 0.9 s (15.1×); at n=1000, 104.3 s → 3.6 s (28.8×).

### Adaptive parameter selection

Fusang automatically selects optimal k and gap parameters based on the number of input taxa n:
- n ≤ 100: k=4,gap1 (emphasizes local conservation for shallow divergence)
- n > 100: k=5,gap2 (captures intermediate-range correlations for moderate divergence)

This adaptive strategy was derived from stability benchmarking across n=20/50/100/200/500 (see Results). Manual override is available via command-line flags.

### Benchmark design

#### Simulated data (without indels)

All simulated benchmark datasets were generated with INDELible v1.03 [19]. Substitutions followed a general time-reversible (GTR) model with exchangeabilities a=0.2, b=0.4, c=0.6, d=0.8, e=1.2 (f=1, i.e. the standard INDELible GTR example rates) and equal base frequencies; among-site rate heterogeneity was modelled with a discrete gamma distribution (shape α=1.0, 4 rate categories). Trees were drawn from a birth–death prior (speciation rate λ=1.0, extinction rate μ=0.5) conditioned on n extant tips; extinct lineages were pruned, tips were extended to the sampling time by an exponential waiting time, and branch lengths were rescaled so that the mean branch length corresponds to 0.005 expected substitutions per site — the divergence level of the pilot simulations (substitution rate 0.05 subs/site per unit branch length × branch length 0.1). Root sequence length L=500 bp. Dataset sizes: n=20, 50, 100, 200, 500, 1000, 10000.

**Sequence length choice (L=500 bp):** This length was selected to represent typical protein-coding gene sequences (e.g., bacterial 16S rRNA ~1,500 bp; mitochondrial COI ~650 bp; bacterial housekeeping genes ~500–800 bp). While longer sequences (L>1000 bp) would provide more phylogenetic signal, L=500 bp represents a challenging real-world scenario where alignment errors are non-negligible. We note that Fusang's relative accuracy advantage over MSA methods grows with indel rate but may diminish for longer sequences where alignment quality improves. Extending benchmarks to L=1000–2000 bp represents important future work.

#### Simulated data (with indels)

Insertions and deletions occurred alongside substitutions at equal rates, with lengths drawn from a geometric distribution of mean 3 bp. In INDELible this is specified as a negative binomial `NB q=2/3, r=1`, whose mean is 1 + r·q/(1−q) = 3 bp; we verified the realised indel-length distribution against theory in a two-taxon control experiment (observed SD 3.8 bp vs theoretical 3.9 bp). Indel rates are reported throughout as the absolute per-site indel rate (0.005, 0.01, 0.02, 0.05), i.e. ρ = 0.1, 0.2, 0.4 and 1.0 relative to the substitution rate. Because INDELible's `[indelrate]` parameter sets the insertion rate and the deletion rate *separately*, the control-file value was set to ρ/2 so that the total indel event rate matches the intended ρ. Multi-seed benchmarks used 130 independent seeds (birth–death tree seed and INDELible random seed both derived from the same integer) for statistical power.

**Simulation reproducibility.** The tree generator (`gen_indelible_sim.py`; pure Python standard library, deterministic given a seed), the INDELible control files it emits, and the benchmark driver (`run_bench_indelible.py`) are included in the reproduction package, together with a pinned Conda environment specification (`environment.yml`: Python 3.11, MAFFT 7.526, FastTree 2.2.0, IQ-TREE 2.4.0). Every simulated dataset — FASTA, true tree, alignment, and both inferred trees — is regenerable from its seed with the exact software versions recorded here.

#### Statistical framework

For multi-seed benchmarks, we report:
- Mean and standard deviation of nRF across seeds
- Wilcoxon signed-rank test p-values (paired per-seed comparison)
- Cohen's d effect size with 95% bootstrap confidence intervals
- Bonferroni correction for multiple comparisons across datasets (5 ground truth datasets tested — n=200 clean, n=200 indel, n=500 clean, n=500 indel, n=1000 clean; adjusted α = 0.05/5 = 0.01)
- Benjamini-Hochberg FDR correction as a less conservative alternative. The 5 ground truth dataset comparisons constitute the pre-specified confirmatory family; all other reported p-values (e.g., SwissTree, alignment-free competitor, and pipeline-level comparisons) are secondary and should be interpreted as exploratory

All statistical analyses were performed in Python using scipy.stats. We note that statistical power varies across dataset sizes: the 130-seed benchmark (n=200, indel) provides adequate power (≥80%) to detect a medium effect size (Cohen's d=0.5) at α=0.05 (two-sided paired test), while 30-seed benchmarks provide moderate power (~75% for d=0.5) and should be interpreted cautiously for non-significant results.

**Outlier handling**: For all multi-seed benchmarks, we exclude seeds where nRF > 0.3 for either method. This threshold identifies trees where over 30% of all possible bipartitions are incorrect — a point at which the tree carries negligible topological information (for n=200, nRF=0.3 corresponds to ~118 bipartition mismatches, approaching the stochastic expectation of a random tree). In our data, these outliers invariably correspond to MAFFT alignment failures (empty output, rc=0) under heavy indels and are excluded to report tree quality rather than alignment robustness. Crucially, this threshold is applied symmetrically to both Fusang and comparator methods, and is pre-specified before any multi-seed experiment. In the 130-seed INDELible benchmark (n=200, indel rate=0.02, seeds 100–229) all 130 seeds produced valid trees and none exceeded the nRF>0.3 threshold, so no exclusion was applied to the reported statistics; the rule is retained for the remaining conditions, where it does remove occasional MAFFT failure seeds. Complete raw data are provided in Supplementary Table S4.

#### Accuracy metric

Normalized Robinson-Foulds distance (nRF):
nRF = (FP + FN) / (2n − 6)

where FP and FN are false positive and false negative bipartition counts. nRF=0: perfect match; nRF=1: complete disagreement.

#### Comparison methods

- **Fusang: Tardigrade Edition** (this work): spaced k-mers (k=5,gap2), simplified or DCM pipeline, FastME BIONJ+BNNI
- **Fusang v1** [6]: deep learning-based, 4–40 taxa limit
- **FastTree2 v2.2.0** [20]: GTR+CAT approximation, MAFFT alignment [21]
- **RAxML-NG** [22] (version not captured in the run logs; the in-repo benchmark reports carry an evidently erroneous "v2.0.1" label, so no version is asserted): GTR+Γ, 10 parsimony + 10 random starting trees, MAFFT alignment
- **IQ-TREE2 v2.4.0** [23]: GTR model (fixed via `-m GTR`, not auto-selected), MAFFT alignment. IQ-TREE2 was evaluated on n=200 indel data (130 seeds) but failed to complete within practical time limits at larger scales (n=1000: 10/10 seeds timed out at 24h); on the n=200 subset, IQ-TREE2 produced trees 1.9× worse than Fusang's k-mer NJ (nRF=0.231 ± 0.032 vs 0.122 ± 0.024, paired Wilcoxon p=4.5×10⁻²³, paired Cohen's d=3.56; Fusang better in 130/130 seeds). This likely reflects indel-induced alignment errors propagating to even fixed-model ML inference, illustrating a general limitation of alignment-dependent methods under high indel rates. Full data available in Supplementary Note S10.
- **Mash v2.3** [24]: contiguous k-mers (k=21), MinHash Jaccard, NJ
- **Mash + FastME**: Mimics Fusang pipeline with contiguous k-mers (Mash distance + FastME), serving as a clean ablation control isolating spaced k-mers from pipeline architecture
- **andi-approx** (tested): Python approximation of suffix array-based anchor distance [25]. Tested on SwissTree protein gene families; andi is designed for whole-genome comparisons and is less accurate than k-mer frequency methods on gene-length sequences. See Supplementary Note S4.
- **Co-phylog** [26]: k-mer frequency + covariance matrix eigenvalues. Tested on both DNA and protein benchmarks (see Results).
- **kmacs** [12]: k-mismatch average common substring distance. Tested on the DNA indel benchmark (130 seeds) with a k-parameter scan (k∈{3,5,10}; best configuration k=3, see Results). Compiled from the original source distribution and executed on WSL2 Ubuntu-24.04.

#### Hardware

All benchmarks: Intel Xeon E-2124 (4 cores/4 threads, 3.3 GHz), 32 GB RAM, Windows 10. Fusang runs natively; RAxML-NG/MAFFT executed under WSL2 (Ubuntu 24.04).

### Real data validation: 16S rRNA type strains

We downloaded 74 representative 16S rRNA sequences spanning 13 phyla (NCBI GenBank, Supplementary Table S1) and ran Fusang (k=5,gap2, simplified pipeline) without alignment. The true composition of the analysed set, established by re-annotating every downloaded accession against NCBI Taxonomy, is 73 sequences across 13 phyla (7 represented by ≥2 taxa), 21 classes (9 with ≥2 taxa) and 40 orders (14 with ≥2 taxa), including 5 archaeal sequences (4 Methanobacteriota, 1 Thermoproteota); one intended taxon (*Sulfolobus solfataricus*) failed to download, giving 73 analysed sequences of which 71 carry usable taxonomy. Since no gold-standard reference tree exists, we evaluated topological quality by comparing tree-based pairwise distances against NCBI taxonomic classifications using permutation tests, and by testing monophyly recovery of taxonomic groups. We note that the accession→organism mapping originally shipped with this dataset was erroneous at the genus/species level for all 73 accessions; all analyses reported here therefore use the corrected NCBI Taxonomy annotation, and the incident is disclosed with the remediation steps in Supplementary Note S17.

### Real genome-scale validation: AFproject community datasets

To test the indel-robustness claims on real data beyond simulations, we benchmarked on two real whole-genome datasets from the AFproject community benchmark [14], each with a published trusted reference tree: (1) Fish mtDNA — 25 fish mitochondrial genomes (~17 kb each, spanning multiple orders; reference tree from Fischer et al. 2013); (2) E. coli/Shigella — 27 genomes (~4.5–5.5 Mb; reference tree from Skippington & Ragan 2011). Sequences were downloaded from NCBI (accessions in Supplementary Table S16). For the fish dataset we additionally computed the MAFFT alignment (n=25, 18,284 columns) to quantify indel content: 9.1% of alignment cells are gaps. For the SwissTree protein benchmark we quantified gap content analogously: median 50% (range 11–81%) across the 11 families, confirming that this benchmark is strongly indel-rich. All methods built trees via Neighbor-Joining (BioPython) on unaligned input, scored with ETE3 Tree.compare (unrooted nRF, the AFproject scoring standard); kmacs and Mash were run at both our working configurations and the exact configurations published in the AFproject benchmark, providing an external cross-validation of the benchmark harness (data, NJ, and ETE3 scoring) rather than of Fusang itself. Whole-genome MSA was not computed for the 27-genome E. coli/Shigella dataset (alignment of 27 × 4.6 Mb is impractical — precisely the regime where alignment-free methods are necessary); for the fish dataset, MAFFT+FastTree2 served as the MSA+ML comparator. We note that the ETE3 nRF used here (informative splits only, max_rf = 2(n−3)) is not numerically comparable with the internal nRF implementation used in Tables 1–11 (which counts trivial splits in the shared set).

### Cross-platform compatibility and deployment

Fusang is implemented in Python 3.10+ with a modular architecture that supports both command-line and script-based workflows. Pre-compiled FastME binaries are provided for Windows (x86-64) and Linux (x86-64), ensuring consistent numerical results across platforms. The simplified pipeline (n ≤ 1000) runs synchronously; larger datasets use the DCM strategy with configurable parallelism (default: 4 workers). A containerized Docker image (python:3.10-slim base) is provided for reproducible deployment.

---

## RESULTS

### A simplified pipeline achieves state-of-the-art accuracy

Systematic ablation of the Fusang pipeline revealed a critical finding: the divide-and-conquer (DCM) strategy, while essential for scalability, requires careful tuning to avoid accuracy loss at small-to-medium scales. On n=200 indel simulated data (indel rate=0.02, seed=42), step-wise tracing from the simplified pipeline to the full DCM+EPA pipeline isolates the error source (Figure 3, Supplementary Figure S3, Supplementary Table S3; raw trace `diag_dcm_loss_n200_s42.json`):

- Simplified pipeline (k-mer→cosine→NJ): **nRF=0.008**
- Swap in the FastME tree builder: nRF=0.015
- Swap in the `distances.py` distance-matrix implementation: nRF=0.013
- Add TF-IDF weighting: nRF=0.013
- Full DCM with EPA grafting: nRF=0.251

For calibration, FastTree2 on the MAFFT alignment of the same data achieves nRF=0.131, and across nine pipeline variants on this seed the full-DCM nRF ranged 0.180–0.251 (mean 0.213 ± 0.083). The first four steps confirm that the clustering and tree-building logic is sound; the large error appears only at the final step, isolating EPA (Evolutionary Placement Algorithm [16]) grafting as the primary error source. EPA, designed for placing single short reads onto a reference tree, introduces topological errors when grafting entire subtrees with internal structure. An earlier informal trace on a different dataset instance (seed=100) reported a larger step-5 error (nRF=0.388); we could not reproduce that value on the canonical seed-42 benchmark and therefore report the reproducible value here. Multi-seed validation confirms the directional advantage at smaller magnitude (see Section "130-seed benchmark confirms directional advantage").

Based on this finding, we implemented an adaptive pipeline: for n≤500, Fusang uses the simplified pipeline (k-mer→cosine→NJ plus BME-scored NNI refinement for n≥100); for n>500, Fusang switches to DCM+EPA, which provides essential scalability and a small accuracy advantage at scale (n=1000: DCM nRF=0.083 ± 0.005, 9 seeds 401–409, vs simplified nRF=0.102 on seed 401). The `SIMPLE_THRESHOLD=500` parameter controls this transition.

**Note on simplified pipeline accuracy**: The simplified pipeline is competitive with, but does not systematically exceed, MSA-based methods on clean data. On 130 INDELible seeds at n=200 clean (Table 1), Fusang attains nRF=0.144 ± 0.025 vs FastTree2 0.146 ± 0.025 (paired Wilcoxon p=0.242, not significant; Fusang better in 68/130 seeds). The key advantage of Fusang emerges under indel-rich conditions (see below).

### Spaced k-mers close the accuracy gap on clean data

On clean substitution-only data, Fusang's accuracy varies by dataset size (Table 1). On n=200 indel-rich data (indel rate=0.02), Fusang significantly outperforms FastTree2 (nRF: Fusang 0.122 ± 0.024 vs FastTree2 0.134 ± 0.022, 130 seeds, paired Wilcoxon p=3.5×10⁻⁷). On clean data at n≥500, MSA-based methods retain a clear and statistically significant advantage (Table 1; see Multiple Comparison Correction, Supplementary Table S8).

**Table 1. Accuracy on clean and indel-rich data (L=500 bp, INDELible birth–death trees; 130 seeds for n=200, 30 seeds for n=500/1000). All nRF values are TRUE-relative (vs simulated ground-truth tree) except where noted.**

| n | Data Type | Fusang nRF ↓ (30 seeds) | FastTree2 nRF ↓ (30 seeds) | Winner |
|---|-----------|---------------------------|----------------------------|--------|
| 200 | Clean | 0.144 ± 0.025 (k=5,gap2) | 0.146 ± 0.025 | Tie (n.s., p=0.242) |
| 200 | Indel (0.02) | **0.122 ± 0.024** (k=5,gap2) | 0.134 ± 0.022 | **Fusang** (p=3.5×10⁻⁷) |
| 500 | Clean | 0.147 ± 0.014 (k=5,gap2) | **0.139 ± 0.011** | FT2 (p=0.004) |
| 500 | Indel (0.02) | **0.125 ± 0.015** (k=5,gap2) | 0.131 ± 0.012 | Fusang (n.s., p=0.056) |
| 1000 | Clean | 0.115 ± 0.011 (k=5,gap2) | **0.091 ± 0.010** | FT2 |
| 1000 | Indel (0.02) | 0.037 ± 0.006 (k=5,gap2) † | — | — |

nRF=0: perfect match. Best result in **bold**. Values are mean ± standard deviation. Seed counts differ by condition: 130 seeds for n=200 (clean: seeds 100–229; indel: seeds 100–229), 30 seeds for n=500 and n=1000. All nRF values are TRUE-relative (vs the simulated ground-truth tree) except where noted.
† The n=1000 indel row is FT2-relative (Fusang tree vs FastTree2 tree; true-tree comparison not available at this scale) and is therefore not comparable with TRUE-relative values elsewhere. The Abstract reports the same 130-seed benchmark (nRF=0.122, seeds 100–229) for the n=200 indel condition that underlies this table row. That benchmark achieved p=3.5×10⁻⁷ (paired Wilcoxon signed-rank test on TRUE-relative nRF). After Bonferroni correction across the 5 ground-truth datasets (α=0.01), the significant comparisons are reported in Supplementary Table S8.

Importantly, Fusang achieves competitive accuracy with **zero sequence alignment**, operating directly on raw FASTA sequences. The multi-k ensemble variant (see below) provides accuracy comparable to the best single-k configuration without manual k selection (its improvement over the default spaced configuration is seed-set-dependent; see the Multi-k section). On clean data at n≥500, MSA-based methods retain a clear advantage (p<0.001 after correction), indicating that full positional information from alignment benefits ML inference when indels are absent.

### Spaced k-mers vs MinHash approaches under indels

To evaluate the robustness of k-mer cosine distances against a widely-used alignment-free alternative, we compared Fusang (spaced k=5,gap2, cosine + NJ) with Mash (contiguous k=21, MinHash Jaccard + NJ) on both clean and indel-rich data (n=200, sub=0.05, indel=0.02), each with 30 seeds, benchmarking against the true simulated coalescent tree (per-seed data: Supplementary Table S14). Note that this benchmark (Table 2) uses an independent coalescent simulation protocol, distinct from the guide-tree protocol used in Table 1 and the indel-rate scan in Table 3; absolute nRF values are higher for all methods under this protocol and should not be compared across benchmarks.

**Table 2. Spaced k-mer vs MinHash robustness (n=200, vs true tree).**

| Method | Data Type | Mean nRF ↓ | Std Dev | n | Indel Degradation |
|--------|-----------|-----------|---------|---|:---:|
| Fusang (k=5,gap2,NJ) | Clean | 0.376 | 0.045 | 30 | — |
| Fusang (k=5,gap2,NJ) | Indel (0.02) | 0.742 | 0.045 | 30 | 1.97× |
| Mash (k=21, MinHash,NJ) | Clean | **0.394** | 0.045 | 30 | — |
| Mash (k=21, MinHash,NJ) | Indel (0.02) | **0.762** | 0.036 | 30 | **1.93×** |

**Note: Mash results are from a 30-seed benchmark (this work, WSL2 Ubuntu 24.04, Mash v2.3); Fusang results are from the same 30-seed coalescent protocol. An earlier single-seed Mash value (nRF=0.162 on clean data) reported during development appears to have been based on an incorrectly parsed tree file and is superseded by this benchmark.**

nRF=0: perfect match; nRF=1.0: random tree. nRF = RF / (2(n−3)). Mash 30-seed benchmark (this work, Mash v2.3, k=21, MinHash Jaccard, NJ tree). Indel degradation = nRF_indel / nRF_clean. Mash at nRF=0.762 on indel-rich data produces trees with substantial topological error (~76% bipartition mismatch), though not statistically indistinguishable from random (nRF=1.0).

Three findings emerge from the multi-seed Mash benchmark (n=30, Mash v2.3, k=21, MinHash Jaccard, NJ tree, vs true tree):

1. **Mash slightly underperforms k-mer cosine distances on clean data.** On clean substitution-only sequences, Mash (nRF=0.394 ± 0.045) is slightly less accurate than Fusang's k-mer cosine (nRF=0.376 ± 0.045, Table 2). The multi-seed benchmark (n=30) provides robust quantification: MinHash Jaccard is not superior to k-mer cosine on clean gene-length sequences.

2. **Both methods degrade severely on indel-rich data, with k-mer cosine retaining marginally lower absolute error.** Mash's indel nRF=0.762 ± 0.036 corresponds to ~76% bipartition mismatch. The indel degradation factor is 1.93× (nRF_indel / nRF_clean); Fusang's k-mer cosine degrades by a comparable 1.97× (nRF=0.376 → 0.742, Table 2). Fusang's absolute indel nRF (0.742) is marginally lower than Mash's (0.762); because a formal between-method test on matched seeds was not performed, this small absolute difference should be interpreted descriptively. A plausible mechanistic explanation for the residual difference is that MinHash sketches discard positional information entirely, making them sensitive to indel-induced length variation, whereas k-mer frequency vectors retain distributional information across the sequence.

3. **MinHash Jaccard is not recommended for phylogenetic inference on gene-length sequences.** Mash was designed for whole-genome distance estimation (sketch size s=1000 by default, k=21), where whole-genome content masking makes MinHash effective. For gene-length sequences (L=500–1000 bp) with indels, the combination of short sequence length and indel-induced length variation produces unreliable MinHash distances. K-mer cosine distances (Fusang) are the recommended alignment-free approach for phylogenetic inference on gene-length sequences.

### Indel robustness: Fusang advantage grows with indel rate

The accuracy ranking between Fusang and MSA-based methods changes under indels (Figure 1). On clean data, MSA methods hold a marginal advantage. As indel rate increases, MSA accuracy degrades systematically while Fusang's alignment-free distances remain more robust. The Fusang advantage grows monotonically with indel rate: from tie at 0.005 to a 13.2% advantage at indel=0.05 (Table 3). At the biologically realistic indel rate of 0.02, Fusang significantly outperforms FastTree2 (TRUE-relative nRF: Fusang 0.122 ± 0.024 vs FastTree2 0.134 ± 0.022, 130 seeds, paired Wilcoxon p=3.5×10⁻⁷).

**Table 3. Indel rate scan: Fusang vs FastTree2 (n=200, L=500 bp; indel rates 0.005, 0.01 and 0.02 use 130 seeds 100–229, 0.05 uses 30 seeds pending completion of the 130-seed scan). All nRF values are TRUE-relative (vs simulated ground-truth tree).**

| Indel Rate | Fusang nRF ↓ (TRUE-rel) | FastTree2 nRF ↓ (TRUE-rel) |
|------------|------------------------|----------------------------|
| 0.005 | 0.137 | 0.141 |
| 0.01 | 0.133 | 0.138 |
| 0.02 | **0.122** | 0.134 |
| 0.05 | **0.066** | 0.076 |

Both columns are TRUE-relative and directly comparable per seed. The Fusang advantage over FastTree2 grows monotonically in magnitude with indel rate — +2.3% at 0.005 (129 valid seeds, Wilcoxon p=0.222, n.s.), +3.3% at 0.01 (130 seeds, p=0.052, borderline), +9.0% at 0.02 (130 seeds, p=3.5×10⁻⁷, significant) and +13.2% at 0.05 (30 seeds) — becoming statistically significant only from indel rate 0.02 upward (calculated as (FT2_nRF − Fusang_nRF) / FT2_nRF × 100%). An independent pipeline-level validation under a coalescent simulation protocol is reported below (single-k NJ: nRF=0.743, multi-k NJ: nRF=0.583), where absolute nRF values are higher for all methods.

The Fusang advantage increases with indel rate across the tested range (0.005–0.05). At indel=0.05, Fusang achieves a 13.2% relative advantage (nRF 0.066 vs 0.076), while at biologically typical rates (0.01–0.02) the advantage is modest (+4.5–5.9%). Real biological indel rates typically fall in the 0.01–0.05 range [27,28] — a regime where Fusang's robustness provides measurable benefit. The monotonic trend suggests Fusang's advantage may further increase at higher indel rates (>0.05), though phylogenetic signal eventually degrades for all methods.

### 130-seed benchmark validates indel robustness advantage

To rigorously assess the indel robustness, we conducted a **130-seed benchmark** (n=200, L=500 bp, indel rate=0.02, k=5,gap2) with paired Wilcoxon signed-rank test (Figure 2, Supplementary Figure S4, Supplementary Table S4). All 130 target seeds (100–229) produced valid tree files, and none exceeded the pre-specified nRF>0.3 outlier threshold, so the reported statistics use all 130 seeds without exclusion:

- **Overall (130 seeds, both methods evaluated against the simulated ground-truth tree)**: Fusang nRF=0.122 ± 0.024 vs FastTree2 nRF=0.134 ± 0.022; Cohen's d=−0.51 [95% CI: −0.69, −0.33]; Fusang lower nRF in 91/130 seeds (70.0%)

The 130-seed results show a consistent and statistically decisive directional advantage with a medium effect size (paired Wilcoxon p=3.5×10⁻⁷). The 95% bootstrap confidence interval for Cohen's d lies entirely below zero ([−0.69, −0.33]), so the advantage is robust rather than marginal. The earlier borderline result reported for this comparison (p=0.052 on 112 seeds) was obtained on data from a legacy JC69 simulator with fixed branch lengths; replacing it with the present INDELible GTR+Γ birth–death protocol both increased the absolute nRF of both methods and strengthened the significance of the Fusang advantage (see Supplementary Figure S6 for effect size analyses across benchmarks).

The per-seed nRF distributions reveal that Fusang's standard deviation (σ=0.024) is comparable to FastTree2's (σ=0.022), indicating stable performance across replicates. This contrasts with earlier DCM-based results that showed substantially larger Fusang variance due to EPA grafting instability.

### Multi-k distance ensemble is comparable to the best single-k configuration

To improve upon the single spaced k-mer configuration (k=5,gap2), we investigated whether fusing distance matrices from multiple k-mer sizes could capture complementary phylogenetic signal. We compute contiguous k-mer cosine distance matrices for k=5, 7, and 9, then average the three matrices before building a single NJ tree (Table 4). Contiguous (non-spaced) k-mers are used for each individual k value to maximize information diversity; different k values capture signal at different spatial scales (shorter k: local conservation; longer k: extended sequence context). All configurations in Table 4 are evaluated on the same 130-seed INDELible benchmark (n=200, indel rate=0.02, seeds 100–229) and, except where noted, use neighbour-joining without NNI refinement so that distance-matrix effects are isolated from the tree-search refinement step.

**Table 4. Multi-k ensemble vs single-k configuration (n=200, indel rate=0.02, 130 seeds 100–229, INDELible). All nRF values are TRUE-relative (vs the simulated ground-truth tree).**

| Method | k-mer Config | Mean nRF ↓ | Std Dev | Ensemble wins / 130 |
|--------|------------|-----------|---------|:---:|
| k=5 contiguous | k=5, no gap | 0.136 | 0.024 | **83/130** (p=1.6×10⁻⁷) |
| k=7 contiguous | k=7, no gap | **0.127** | 0.023 | 44/130 (n.s., p=0.23) |
| k=9 contiguous | k=9, no gap | 0.131 | 0.025 | **73/130** (p=7.3×10⁻⁴) |
| **Multi-k ensemble** | **avg(k=5,7,9)** | **0.128** | **0.024** | — |
| Fusang k=5,gap2 (spaced, NJ only) | k=5, gap2 | 0.132 | 0.023 | 65/130 (p=0.013) ‡ |
| Fusang k=5,gap2 (full pipeline) | spaced + BME-NNI | **0.122** | 0.024 | 45/130 (p=2.2×10⁻⁴) ‡ |

nRF=0: perfect match. All values are TRUE-relative (see caption). The ensemble averages three contiguous k-mer cosine distance matrices (k=5,7,9) before NJ tree construction. All rows except the last use NJ without NNI refinement; the full-pipeline row is the Fusang CLI product (NJ plus BME-scored NNI) and is shown for context. ‡For the spaced NJ-only row, "ensemble wins" means the ensemble achieved lower nRF (65 wins, 22 ties, p=0.013 favouring the ensemble); for the full-pipeline row, the ensemble wins only 45/130 (p=2.2×10⁻⁴ favouring the full pipeline).

**Paired comparison: Multi-k ensemble vs single-k components (130 seeds, TRUE-relative)**:
- The ensemble (nRF=0.128) significantly outperforms contiguous k=5 (0.136; 83/130 wins, Wilcoxon p=1.6×10⁻⁷, paired d=0.49) and contiguous k=9 (0.131; 73/130 wins, 20 ties, p=7.3×10⁻⁴, d=0.30).
- The ensemble is statistically indistinguishable from the best single-k configuration, contiguous k=7 (0.127; 44 wins, 30 ties, p=0.23, d=0.07).
- Against the spaced k=5,gap2 configuration under the same NJ-only protocol, the ensemble is significantly better (0.128 vs 0.132; 65/130 wins, 22 ties, p=0.013, d=0.24) — but the comparison mixes two changes (multi-k fusion and contiguous-vs-spaced patterns); the like-for-like single-k spaced-vs-contiguous comparison (Table 7, Finding 3) shows the spaced pattern itself is significantly better at k=5 (p=0.008).
- The full Fusang pipeline (spaced + BME-NNI, nRF=0.122) significantly outperforms the ensemble (45/130 ensemble wins, p=2.2×10⁻⁴, d=0.36).

This result demonstrates that distance matrix fusion across multiple contiguous k-mer resolutions achieves accuracy comparable to the best single-k configuration (k=7), with significant gains over k=5 and k=9 individually, while not requiring users to pre-select the optimal k value for their dataset (per-seed data: Supplementary Figure S7, Supplementary Table S9). The ensemble does not, however, exceed the full default pipeline at L=500 bp, so we position it as an automatic, robust fallback rather than an accuracy upgrade at gene scale. The ensemble is implemented as a standalone script (`benchmark_multik_ensemble.py`) that averages the k=5/7/9 contiguous cosine distance matrices before NJ tree construction. The ensemble's decisive advantages appear under different conditions: the large single-k→multi-k gain under the L=1000 bp coalescent protocol (Table 5) and the automatic rescue of the saturating gene-scale default at genome scale (see Real genome benchmarks below).

### Multi-k NJ: comparison with MSA+ML

To determine whether the multi-k ensemble's improvement translates to MSA+ML-level accuracy, we conducted a pipeline-level validation against the true simulated coalescent tree (n=200, L=1000 bp, sub=0.05, indel=0.02, 30 seeds). We compared three levels of the Fusang multi-layer pipeline (Table 5): Level 0 (L0: single k-mer k=5,gap2 cosine + NJ), Level 1 (L1: multi-k k=5,7,9 contiguous cosine average + NJ), and Level 3 (L3: MAFFT v7 alignment + FastTree2 GTR+CAT, nucleotide mode). All 30 seeds completed for every level; the L3 level was executed on Linux (WSL2 Ubuntu-24.04, MAFFT --auto + FastTree2 `-nt -gtr -nosupport`), replacing an earlier 5-seed Windows run whose MAFFT instability had limited the comparison (see Supplementary Note S8 for the protocol correction and limitation analysis). The pipeline-level benchmark uses L=1000 bp rather than the L=500 bp of Tables 1 and 3: the coalescent protocol was fixed during multi-layer pipeline development — before the L1-vs-L3 comparison was designed — with L=1000 bp chosen so that the k=9 member of the multi-k set has non-sparse frequency vectors at n=200. The protocol was identical for the superseded 5-seed run and the full 30-seed rerun, and all levels (L0, L1, L3) are compared within this single protocol, so level comparisons are internally fair.

**Table 5. Pipeline-level validation against true tree (n=200, L=1000 bp, sub=0.05, indel=0.02, 30 seeds, all levels complete). All nRF values are TRUE-relative (vs simulated coalescent ground truth). Note: this benchmark uses a coalescent simulation protocol with L=1000 bp, independent of the guide-tree protocol used in Tables 1 and 3 (L=500 bp); absolute nRF values are higher for all methods under this protocol and are not directly comparable across benchmarks. Per-seed data: Supplementary Table S13.**

| Level | Method | Mean nRF ↓ | Std Dev | n | vs L0 Wilcoxon p ¹ |
|-------|--------|-----------|---------|---|:---:|
| L0 | k-mer k=5 NJ | 0.743 | 0.046 | 30 | baseline |
| **L1** | **Multi-k (k=5,7,9) NJ** | **0.583** | **0.044** | **30** | **< 0.0001** |
| L3 | MAFFT + FastTree2 GTR | 0.601 | 0.055 | 30 | < 0.0001 ² |

¹ Wilcoxon signed-rank paired test against L0 (single-k NJ) within the same seed.  
² L3 also significantly outperforms L0 (paired Wilcoxon p<0.0001, paired d=3.25); the comparison of primary interest is L1 vs L3 (see text).

nRF=0: perfect match; nRF=1.0: random tree. nRF = RF / (2(n−3)). All values vs true coalescent tree. L0 and L1 complete across all 30 seeds; L3 complete across all 30 seeds on Linux. L1 vs L0: Wilcoxon signed-rank test p<0.0001, paired Cohen's d=3.82 (large effect).

Two key findings emerge from this validation:

1. **Multi-k NJ significantly outperforms MSA+ML under this indel-rich coalescent simulation.** L1 (multi-k ensemble NJ, nRF=0.583±0.044) achieves lower topological error than L3 (MAFFT+FastTree2, nRF=0.601±0.055) across the full 30-seed comparison: mean per-seed difference −0.018, paired Wilcoxon p=0.024, paired Cohen's d=−0.45 (small-to-medium effect), with L1 winning 19/30 seeds (63.3%). We note that this comparison is secondary (exploratory) in our statistical framework; after Bonferroni correction across the three pipeline-level pairwise comparisons (adjusted α=0.0167), the result is borderline (adjusted p=0.071). We therefore interpret it as evidence that multi-scale k-mer information captures phylogenetic signal at least at the level of full MSA+ML inference under indel-rich conditions — without ever computing a sequence alignment — rather than as a definitive superiority claim.

2. **Multi-k fusion provides a 21.5% relative accuracy improvement over single-k**, reducing nRF from 0.743 to 0.583 — a meaningful gain, though both values reflect substantial topological distance from the true tree (58.3% of bipartitions differ). L1 (nRF=0.583) substantially outperforms L0 (nRF=0.743, p<0.0001, paired d=3.82). The unusually large effect size (paired d=3.82) reflects that single-k (L0) operates near the information ceiling of a single k-mer resolution under heavy indels — the high baseline nRF (0.743) and low variance (σ=0.046 for L0, σ=0.044 for L1) combine to yield d=3.82, indicating that most of the gain comes from resolving splits that are nearly random for single-k but recoverable with multi-scale k-mer information. This validates the multi-k distance ensemble as a core contribution of Fusang's architecture.

These results position the multi-k ensemble NJ as a practical alternative to MSA+ML for indel-rich datasets at moderate scales. The trade-off is notable: L1 requires no alignment step and completes in ~15 seconds per seed (vs ~31 seconds for L3 in the same Linux environment, a 2× speedup that grows further at larger scales where MSA cost dominates). Both methods remain substantially distant from the true tree under this deliberately harsh coalescent protocol (nRF≈0.58–0.60), so the practical takeaway is comparable-topology-at-far-lower-cost rather than high absolute accuracy. The L1≈L3 agreement is consistent with the broader pattern of k-mer methods matching MSA accuracy under indels (Tables 1, 3, 4).

### Boundary classifier accurately distinguishes dataset homogeneity

A critical component of Fusang's multi-layer pipeline is the Level 2 boundary classifier — a random forest (RF) model that determines whether a given node in the clustering hierarchy represents a homogeneous set of taxa (STOP) or requires further splitting (SPLIT). We validated the classifier (RF V4b, trained on 4,073 labeled samples including both homogeneous and structured scenarios) on an independent test set of 88 scenarios spanning two fundamentally different tree types (Table 6; per-scenario results: Supplementary Table S15):

- **Coalescent trees** (13 scenarios): Single-population coalescent simulations (n=50, Kingman coalescent) where all taxa evolve under the same demographic model. These represent truly homogeneous datasets that should not be split.
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

We compared Fusang against established alignment-free phylogenetic methods on simulated indel-rich data (n=200, sub=0.05, indel=0.02, 130 independent seeds 100–229, INDELible birth–death protocol; Table 7; per-seed data: Supplementary Table S10): Co-phylog (context-object matching), kmacs (k-mismatch average common substring [12]), and k-mer cosine baselines, plus Skmer [29] (genome-skim distance estimation), which proved structurally inapplicable at gene length (see below). Table 7 reports TRUE-relative accuracy (all methods' trees compared against the simulated ground-truth tree of the same seed) as the primary frame; FT2-relative values (vs the FastTree2 tree of the same seed) are given in the table footnote. We use the TRUE-relative frame as primary because the FT2-relative frame has poor discriminative power at this scale: FastTree2's own TRUE-relative error (0.134 ± 0.022) is of the same magnitude as every method's FT2-relative distance (~0.13–0.14), so proximity to the FT2 tree does not imply proximity to the truth.

**Table 7. Alignment-free method comparison on indel-rich data (n=200, sub=0.05, indel=0.02, 130 seeds 100–229). All nRF values are TRUE-relative (vs the simulated ground-truth tree).**

| Method | Type | Mean nRF ↓ | Std Dev | vs Fusang (full pipeline) |
|--------|------|-----------|---------|:---:|
| **Co-phylog** (Yi & Jin 2013) | Context-object, k=19 | **0.416** | 0.022 | 3.4× worse (p=4.3×10⁻²³) |
| **kmacs** (Leimeister & Morgenstern 2014) | k-mismatch ACS, k=3 (best) | **0.221** | 0.025 | 1.8× worse (p=4.4×10⁻²³) |
| **K-mer cosine k=5** (contiguous) | Frequency vector, NJ | **0.136** | 0.024 | worse (p=3.7×10⁻¹²) |
| **K-mer cosine k=7** (contiguous) | Frequency vector, NJ | **0.127** | 0.023 | worse (p=1.4×10⁻³) |
| **Fusang** multi-k ensemble | Multi-k contiguous (k=5,7,9 avg), NJ | **0.128** | 0.024 | worse (p=2.2×10⁻⁴) |
| **Fusang** (k=5,gap2, spaced), NJ only | Spaced k-mer, no NNI | **0.132** | 0.023 | worse (p=2.1×10⁻¹²) |
| **Fusang** (k=5,gap2, spaced) + BME-NNI | Spaced k-mer, full pipeline | **0.122** | 0.024 | — |
| FastTree2 (MSA-based) | MSA + ML | **0.134** | 0.022 | 1.10× worse (p=3.5×10⁻⁷) |

nRF=0: perfect match; nRF=1.0: random tree. All alignment-free baselines (Co-phylog, kmacs, contiguous single-k, multi-k ensemble, and the NJ-only Fusang row) use neighbour-joining (BioPython) without NNI refinement; the full-pipeline Fusang row is the CLI product (NJ plus BME-scored NNI refinement, applied at n≥100). Wilcoxon signed-rank tests vs the full-pipeline Fusang configuration (paired Cohen's d, positive d = Fusang better): Co-phylog p=4.3×10⁻²³ (d=9.34), kmacs k=3 p=4.4×10⁻²³ (d=3.84), k-mer cosine k=5 p=3.7×10⁻¹² (d=0.75), k=7 p=1.4×10⁻³ (d=0.30), multi-k ensemble p=2.2×10⁻⁴ (d=0.36), Fusang NJ-only p=2.1×10⁻¹² (d=0.75), FastTree2 p=3.5×10⁻⁷ (d=0.51). In the FT2-relative frame the same ordering holds for the classical methods (Co-phylog 0.404, kmacs k=3 0.220) but the k-mer variants are statistically inseparable (FT2-relative means 0.133–0.138; cosine k=5 vs Fusang p=0.42), consistent with the poor discriminability of that frame noted above. Normalization: max_rf = 2(n−3). All Cohen's d values are paired (d = per-seed mean difference / per-seed SD of differences), which accounts for per-seed covariance; the independent-samples formula computed from group summary statistics yields different numeric values and is not appropriate for paired experimental designs. The d=9.34 value for Co-phylog should be interpreted as indicating near-deterministic superiority on every seed rather than a proportionally larger effect (d>2 is conventionally "very large").

Three key findings emerge from this comparison:

1. **Co-phylog fails under indel-rich conditions.** Co-phylog's context-object approach (k=19) produces substantially worse trees (TRUE-relative mean nRF=0.416, 3.4× worse than Fusang, 0/130 wins, p=4.3×10⁻²³, paired d=9.34) on indel-rich data. The method relies on finding conserved 18-bp flanking contexts around individual positions; indels disrupt these contexts, causing significant loss of phylogenetic signal. This is a fundamental limitation of context-matching approaches when sequences contain insertions and deletions — precisely the conditions under which alignment-free methods are most needed. (Co-phylog results are from our Python reimplementation, Supplementary Note S4; on the fish benchmark — see Real genome benchmarks below — the reimplementation slightly exceeds the published configuration's accuracy, indicating the reimplementation is not disadvantaged relative to the original; the DNA conclusion here is therefore conservative in direction.)

2. **kmacs, the closest classical relative of our approach, is significantly less accurate under indels.** kmacs extends the average common substring (ACS) approach with k-mismatch matching, making it the most directly comparable classical method. We scanned k∈{3,5,10} and found accuracy degrades monotonically with increasing k (TRUE-relative: k=3 0.221, k=5 0.243, k=10 0.262; all 0/130 wins vs Fusang); we report the best configuration (k=3). Even at its best, kmacs achieves nRF=0.221±0.025, 1.8× worse than Fusang (paired Wilcoxon p=4.4×10⁻²³, paired d=3.84) — though far more robust than Co-phylog (0.416). Substring-matching distances are sensitive to indel-disrupted substring sharing even with k-mismatch tolerance, whereas k-mer frequency vectors degrade gracefully because indels perturb frequencies rather than destroy matches.

3. **Spaced k-mers and NNI refinement each contribute significant accuracy gains on DNA.** Comparing like-for-like NJ trees (no NNI on either side), Fusang's spaced k=5,gap2 configuration significantly outperforms the contiguous k=5 cosine baseline (TRUE-relative 0.132±0.023 vs 0.136±0.024; spaced better on 68/130 seeds, paired Wilcoxon p=7.97×10⁻³, paired d=0.234). This reverses the non-significant result (p=0.26) we previously reported from a 27-seed FT2-relative benchmark: the earlier test was underpowered and used the poorly discriminating FT2 reference frame — on the same 130 seeds the FT2-relative comparison remains non-significant (58/130, p=0.48), illustrating how that frame masks real differences. The spaced-pattern advantage on DNA (p=0.008) is directionally consistent with, though weaker than, the protein benchmark (p=0.0078; see Cross-domain validation below), and is consistent with the proposed mechanism: spaced patterns recover k-mers that would otherwise be destroyed by an insertion or deletion falling inside their span, so their benefit grows with alignment gap content. Independently, BME-scored NNI refinement contributes a further significant gain over the pure-NJ spaced tree (0.122±0.024 vs 0.132±0.023; better on 100/130 seeds, p=2.1×10⁻¹², d=0.75), so the full pipeline combines both contributions.

**Multi-k ensemble.** On this 130-seed INDELible benchmark the multi-k ensemble (contiguous k=5,7,9 average, NJ) achieves TRUE-relative nRF=0.128±0.024 — significantly better than the contiguous k=5 baseline (0.136; 83/130 wins, p=1.6×10⁻⁷, d=0.49) and significantly better than the NJ-only spaced configuration (0.132; 65/130 wins, 22 ties, p=1.3×10⁻², d=0.24), but statistically indistinguishable from contiguous k=7 (0.127; p=0.23) and significantly worse than the full Fusang pipeline (0.122; 45/130 wins, p=2.2×10⁻⁴, d=0.36). We therefore do not claim a general accuracy advantage of the multi-k ensemble over the default spaced configuration at L=500 bp; Table 4 examines the ensemble against its single-k components in detail. The ensemble's established benefits are (i) the large single-k→multi-k gain under the L=1000 bp coalescent protocol (Table 5, d=3.82, p<0.0001) and (ii) the automatic rescue of the saturating gene-scale default at genome scale (see Real genome benchmarks below).

**Regarding andi** [25]: andi uses suffix-array-based anchor distances designed for genome-scale sequences (>10kb). On gene-length sequences (~500 bp), andi's anchor-finding mechanism has insufficient maximal unique matches for reliable distance estimation. In our harness andi returned no usable tree on the gene-length families we attempted, so we report **no accuracy figure** for it — the run record (`benchmark_andi_results.json`) captures timings only, with null tree fields for every replicate. We therefore characterise this as a design-domain mismatch rather than a measured accuracy result. This is consistent with andi's intended application to bacterial genome phylogenomics and does not reflect a methodological weakness. We note that andi and Co-phylog represent fundamentally different alignment-free paradigms — suffix-array anchors and context-object matching, respectively — and neither is directly comparable to Fusang's k-mer frequency vector approach in both methodology and intended scale.

**Regarding Skmer** [29]: Skmer v3.3.0 estimates distances from genome skims through a k-mer depth-spectrum coverage and error-correction model. On 500-bp single-gene inputs (tested at the default k=31 and at k=21), its coverage estimator fails with a division-by-zero error in every case: at gene length each k-mer occurs approximately once, producing a degenerate depth histogram on which the coverage model cannot operate, and no distance matrix can be computed. Skmer is thus structurally inapplicable to gene-length phylogenetic inference — an outcome parallel to andi, whose suffix-array anchors likewise require genome-scale input. Both cases reflect design-domain mismatches rather than methodological weaknesses.

### Accuracy at scale: n=1000 indel performance

On indel-rich data at n=1000 (30 seeds, sub=0.05, indel=0.02), Fusang's simplified pipeline (k-mer→cosine→NJ plus BME-scored NNI refinement for n≥100) achieves nRF=0.037 ± 0.006 relative to the FastTree2 reference tree — indicating only 3.7% topological divergence from the MSA+ML reconstruction (Table 1). Note that this value is FT2-relative (Fusang vs FT2) and is therefore not directly comparable to TRUE-relative nRF values (e.g., Table 5), nor to FT2-relative values under different conditions (where the FT2 reference tree differs). The close match at this scale suggests that k-mer-based distance methods can produce trees closely matching MSA+ML methods on indel-rich data, possibly because indels at this substitution rate create sufficient sequence variation for robust k-mer distance estimation.

### Spaced k-mer gap scales with tree size

Systematic parameter scanning (k=3–8, gap=0–4, n=20–1000) revealed a robust relationship between optimal gap and dataset size (Table 8, Figure 4, Supplementary Figure S1, Supplementary Figure S2). The adaptive strategy (k=4,gap1 for n≤100; k=5,gap2 for n>100) was validated through 5-repeat stability testing across all scales (Supplementary Table S5).

**Table 8. Optimal parameters and nRF stability across dataset scales.**

| n | Adaptive (k,gap) | Fusang nRF (clean) | FT2 nRF (clean) | Notes |
|---|:---:|------------|----------|-------|
| 50 | 4,gap1 | 0.102 ± 0.020 | 0.095 ± 0.018 | Comparable |
| 100 | 5,gap2 | 0.115 ± 0.022 | 0.098 ± 0.016 | FT2 better |
| 200 | 5,gap2 | 0.144 ± 0.025 | 0.146 ± 0.025 | Comparable (clean, p=0.242) |
| 200 | 5,gap2 | **0.122 ± 0.024** | 0.134 ± 0.022 | Fusang better (indel, p=3.5×10⁻⁷) ‡ |
| 500 | 5,gap2 | 0.147 ± 0.014 | **0.139 ± 0.011** | FT2 better (clean, p=0.004) |
| 1000 | 5,gap2 | 0.115 ± 0.011 | **0.091 ± 0.010** | FT2 better (clean) |

**Notes**: For n<1000, Notes column indicates which method performs better and the p-value when significant. For n=1000, "Simplified" indicates that the simplified pipeline (k-mer→cosine→NJ plus BME-scored NNI refinement for n≥100) was used rather than DCM+EPA; the winner is still FastTree2 on clean data. See Supplementary Table S6 for complete n=1000 benchmark results (30 seeds, both clean and indel conditions). ‡ Row values are from the 130-seed benchmark (seeds 100–229; all seeds valid, no outlier exclusion), paired Wilcoxon p=3.5×10⁻⁷ (Table 3).

On indel-rich data (indel rate=0.02), the optimal gap shifts slightly: gap3 provides a modest 10.5% improvement over gap2 at n=200 (nRF=0.043 vs 0.048, Supplementary Table S2), attributed to wider spacing better tolerating indel-induced length variation. However, the absolute improvement is small, and k=5,gap2 remains a robust default within 10% of the empirical optimum across all tested conditions.

### Validation on real 16S rRNA sequences

Fusang processed the 73 analysed 16S rRNA sequences without alignment (simplified pipeline, k=5,gap2) in 0.108 s for the k-mer counting plus NJ stage. Tree-based pairwise distances were compared against NCBI taxonomic classifications (Table 9). For reference-tree comparisons we rebuilt the alignment-based comparator end-to-end: the archived alignment file was found to be not genuinely aligned (pairwise identity between the two *Escherichia*/*Bacillus* sequences fell below the 0.25 random baseline), so all sequences were realigned with MAFFT 7.526 (--localpair) and the FastTree2 tree rebuilt from it. Fusang vs this rebuilt FastTree2 tree gives nRF=0.371 — substantial but no longer near-total topological disagreement.

We evaluated both trees against the NCBI taxonomy as an external reference. A taxonomic group is scored as recovered when the MRCA of all tips carrying the label has exactly that tip set. Against the true taxonomy (13 phyla, 7 with ≥2 taxa; 40 orders, 14 with ≥2), the Fusang tree recovers 3/7 assessable phyla and 10/14 assessable orders as monophyletic, while the rebuilt FastTree2 tree recovers 6/7 phyla and 11/14 orders. The RF distance against a taxonomy-derived reference tree (phylum→order→taxon containment hierarchy, polytomous where unresolved) is nRF=0.414 for Fusang vs 0.343 for FastTree2 under the conventional RF/(2n−6) normalisation (0.547 vs 0.453 under the resolution-aware RF/(s₁+s₂) normalisation; both normalisations rank the two trees identically). While Fusang detects substantial genuine signal (Table 9), it falls short of MSA-based accuracy on this dataset (Figure 6, Supplementary Figure S5). This result delineates an applicability boundary of the method rather than contradicting the indel-robustness findings on simulated data: 16S rRNA genes contain highly conserved regions where positional homology from alignment provides strong phylogenetic signal that k-mer frequency vectors, which discard positional information, cannot fully capture.

**Table 9. Real 16S rRNA validation (n=73 taxa spanning 13 phyla and 40 orders under the corrected NCBI Taxonomy annotation; simplified pipeline, k=5,gap2).**

| Taxonomic Level | Same-group Distance | Different-group Distance | Reduction | P-value |
|:---|---:|---:|---:|:---:|
| Phylum | 0.1763 | 0.2570 | **31.4%** | <1×10⁻⁴ |
| Order | 0.1066 | 0.2451 | **56.5%** | <1×10⁻⁴ |
| Family | 0.0709 | 0.2440 | **71.0%** | <1×10⁻⁴ |

P-values are from two-sided permutation tests comparing mean same-group versus different-group tree-based pairwise (patristic) distances with taxonomic group labels permuted (test statistic: difference in group means; 10,000 random permutations per level; `p < 1×10⁻⁴` means no permutation reached the observed statistic). Fusang recovers significant phylogenetic signal at every taxonomic level, with the reduction rising monotonically from phylum (31.4%) to genus (78.6%) — the gradient itself evidencing genuine rather than artefactual signal. On the same sequences the alignment-based FastTree2 tree attains larger reductions (phylum 55.5%, order 80.4%, family 90.1%; all p < 1×10⁻⁴), consistent with the applicability boundary discussed below. Absolute same/different-group distances are not comparable across trees (Fusang NJ branch lengths vs FastTree2 substitution units); only the reductions and p-values are. All values are computed under the corrected NCBI Taxonomy annotation (Supplementary Note S17); an earlier evaluation against the erroneous shipped labels is superseded.

These results show that k-mer frequency features recover genuine phylogenetic signal from real sequences without parameter tuning, while delineating the sequence types — highly structured genes with strong positional conservation — where alignment-based methods retain a clear advantage.

### Cross-domain validation: AFproject SwissTree protein gene trees

To assess whether Fusang's k-mer approach generalizes beyond DNA to protein sequences, we benchmarked on the AFproject SwissTree gene tree dataset [14] — the community standard for alignment-free gene tree inference. This benchmark comprises 11 protein gene families (21–159 taxa, 111–576 amino acids per sequence) with trusted reference trees from the SwissTree database (Table 10; per-family data: Supplementary Tables S11–S12).

**Table 10. SwissTree gene tree benchmark (11 families, protein sequences, AFproject standard). Per-family data: Supplementary Tables S11–S12.**

| Method | Configuration | Mean nRF ↓ | Std Dev | Wins /11 |
|--------|:---|-----------:|--------:|:---:|
| Co-phylog (halfctx=5, k=11) | Context-object | **0.433** | 0.076 | 3 |
| Co-phylog (halfctx=11, k=23) | Context-object | **0.361** | 0.059 | 0 |
| K-mer cosine k=4 | Contiguous | 0.256 | 0.122 | 1 |
| K-mer cosine k=5 | Contiguous | 0.244 | 0.110 | 2 |
| **Fusang** k=4,gap1 (1010101) | Spaced | **0.239** | 0.118 | 5 |
| **Fusang** k=5,gap2 (1001001001001) | Spaced | **0.244** | 0.113 | 5 |

nRF=0: perfect match; nRF=1.0: random tree. Normalization: max_rf = 2(n−3). All methods use NJ (BioPython) for tree construction. "Wins" = lowest nRF in that family; ties are counted for all winners. Co-phylog halfctx=5 wins ST008, ST009, ST012 (nRF=0.3205, 0.3194, 0.4167 vs Fusang k=4,gap1 nRF=0.3462, 0.4028, 0.4722). Statistical tests: Fusang k=4,gap1 vs Co-phylog halfctx=5: paired t-test p=0.005, Wilcoxon p=0.014, Cohen's d=1.08. Fusang k=4,gap1 vs Co-phylog halfctx=11: paired t-test p=0.001, Wilcoxon p=0.006, Cohen's d=1.32. Spaced vs contiguous k-mer (k=4,gap1 vs k=4): exact paired Wilcoxon p=0.0078, Cohen's d=−1.00, with the spaced configuration achieving lower nRF in 8/11 families and never worse in the remaining three (ties).

Three findings emerge from the SwissTree benchmark:

1. **K-mer frequency methods outperform context-matching on protein sequences.** The best k-mer configurations achieve mean nRF≈0.24, while Co-phylog's best configuration (halfctx=11) achieves nRF=0.361 — a 1.5× accuracy advantage (paired Wilcoxon p=0.006, Cohen's d=1.32) — and its default (halfctx=5) achieves nRF=0.433, a 1.8× gap (p=0.014, d=1.08). This extends our observation from DNA data (Table 7) to protein sequences, confirming that k-mer frequency vectors capture phylogenetic signal more robustly than Co-phylog's context-object matching across both sequence alphabets.

2. **Spaced k-mers significantly improve accuracy on protein data.** Spaced k=4,gap1 (mean nRF=0.239) outperforms contiguous k=4 (mean nRF=0.256), achieving lower nRF in 8 of 11 families, tying in the remaining three, and never being worse (exact paired Wilcoxon p=0.0078, Cohen's d=−1.00). Notably, no spaced advantage is detectable on DNA data at low indel rates (Table 7, p=0.26). This domain dependence is consistent with the mechanism we propose: SwissTree alignments are gap-rich (median gap fraction 50%, range 11–81%), so patterns that skip over insertion/deletion sites recover signal that contiguous k-mers lose, whereas at the low indel rates of our DNA benchmarks there is little such signal to recover. We therefore position spaced k-mers as beneficial specifically under indel- and gap-rich conditions rather than as a universal improvement.

3. **Both k-mer representation and distance metric contribute to accuracy across sequence types.** Across both DNA (Table 7) and protein (Table 10) benchmarks, the choice of distance metric (k-mer cosine vs context-object matching vs MinHash Jaccard) has a larger effect than the choice between spaced and contiguous k-mers — but the latter is itself significant on gap-rich protein alignments (finding 2) and on indel-rich DNA (Table 3). The systematic evaluation of k-mer frequency vector cosine distances — across spaced and contiguous patterns and across sequence alphabets — is the core contribution of this work.

### Scalability: from single genes to 10,000 taxa

Fusang completes phylogenetic inference on 10,000 taxa in 69.8 seconds (Table 11, Figure 5). For comparison, RAxML-NG under GTR+Γ with 10 parsimony and 10 random starting trees required 244 s at n=200 (measured; `raxml_bench_results.csv`), against 2.3 s for the Fusang simplified pipeline at the same size. The divide-and-conquer strategy with FastME scales as O(n² log n), with optimizations including scipy-based clustering (52× faster than Python implementation) and an n≤2 fast path.

**Table 11. Fusang scalability (L=500 bp, clean data; simplified pipeline for all n≤1000 timing runs, DCM at n=10,000). Tree builder: FastME BIONJ+BNNI.**

| n | Time (s) | Pipeline | FastME Speedup vs NJ |
|---|:---:|----------|:---:|
| 20 | 4.9 | Simplified | — |
| 50 | 13.9 | Simplified | — |
| 100 | 24.5 | Simplified | — |
| 200 | 46.1 | Simplified (NJ) | 3.5× |
| 500 | 3.8 | Simplified (FastME) | — |
| 1000 | 5.2 | Simplified (FastME) | — |
| 10000 | 69.8 | DCM | — |

† This table mixes two measurement generations and should be read with care. The n=20–200 rows were timed with Neighbor-Joining (NJ) **before** the scipy-based clustering optimisation (reported: 4.9/13.9/24.5/46.1 s; the same NJ builder after optimisation: 1.8/2.0/2.2/2.3 s), whereas the n≥500 rows use FastME BIONJ+BNNI on the optimised code path (n=500/1000 also independently recorded as 3.59/5.08 s). The apparent non-monotonicity between n=200 (46.1 s) and n=500 (3.8 s) therefore reflects the change of tree builder and code version, not an algorithmic discontinuity; a uniform re-timing of the whole table on a single code version is in progress. All scalability claims in the main text use FastME BIONJ+BNNI as the default tree builder for n≥200 (3.5× faster than NJ at n=200).

### Real indel-rich genome benchmarks: AFproject community datasets

The indel-robustness claims rest primarily on simulations; we therefore validated on two real whole-genome datasets with published trusted reference trees from the AFproject community benchmark [14] (Table 12; full results: Supplementary Table S16). Before benchmarking, we quantified the indel content of our real datasets: the fish mitochondrial alignment contains 9.1% gap cells, and the SwissTree protein families (Table 10) contain a median of 50% gap cells (range 11–81%) — establishing that the existing protein benchmark is strongly indel-rich real data.

**Table 12. Real genome benchmarks (AFproject community datasets, ETE3 nRF¹ vs published reference trees). All methods use NJ on unaligned input unless noted.**

| Method | Fish mtDNA (n=25) ↓ | E. coli/Shigella (n=27) ↓ |
|--------|:---:|:---:|
| **Multi-k cosine (k=5,7,9 avg)** | **0.045** | 0.500 |
| **Multi-k cosine (k=7,9,11 avg)** | **0.045** | — |
| Multi-k cosine (k=9,11,13 avg) | — | 0.500 |
| K-mer cosine k=9 / k=11 / k=13 | 0.091 / 0.091 / — | 0.500 / 0.500 / 0.500 |
| K-mer cosine k=5 (gene-scale default) | 0.545 | 0.875 |
| kmacs (k=10, AFproject config) | **0.045** [published 0.05] | — |
| kmacs (k=3) | 0.182 | 0.542 ³ |
| Mash (k=11, s=5000, AFproject config) | **0.045** [published 0.05] | — |
| Mash (k=21, s=10000) | 0.091 | 0.208 |
| Co-phylog (halfctx=5, AFproject config) | **0.045** [published 0.09] | — |
| MAFFT+FastTree2 (MSA+ML) | **0.045** | not feasible ² |
| Published AFproject best (andi/Co-phylog/phylonium) | 0.05 | 0.08 |

¹ nRF computed with ETE3 Tree.compare (unrooted, informative splits only, max_rf = 2(n−3)) — the AFproject scoring standard; not numerically comparable with the internal nRF of Tables 1–11. "—" = not tested (see notes); "n/a" = not applicable. Values in brackets are the scores published in the AFproject benchmark for the same method configuration; kmacs and Mash are original binaries and reproduce the published values (difference ≤ one split, consistent with NJ implementation differences), cross-validating the benchmark harness (data, NJ, and ETE3 scoring) rather than Fusang itself. Co-phylog is our Python reimplementation (as in Table 7), which explains its slightly better-than-published fish score.
² Whole-genome alignment of 27 × ~4.6 Mb is impractical — the regime where alignment-free methods are the only option.
³ kmacs (single-threaded) on the 27 × ~4.6 Mb dataset is far slower than either comparator: an earlier attempt produced no output after >40 min and the run that ultimately yielded a tree was left running overnight, so we report no precise wall-clock figure for it. By contrast the logged runtimes are ~4 min for the full k-mer cosine k-scan (k=5…13; `ecoli_suite_log.txt`) and ~10 s for Mash. The resulting kmacs tree (nRF=0.542) does not beat the k-mer cosine plateau (0.500) on this dataset. Mash and the published AFproject leaderboard remain the stronger genome-scale comparators here.

Three findings emerge:

1. **On deep-divergence indel-rich real data, the multi-k cosine ensemble reaches the community-benchmark ceiling.** On the fish mitochondrial dataset (25 genomes across multiple orders, 9.1% alignment gaps), multi-k cosine achieves nRF=0.045 — identical to MAFFT+FastTree2 and to kmacs, Mash, and Co-phylog at their AFproject-published best configurations. This dataset does not discriminate among these methods: with n=25 taxa (44 informative splits), 0.045 corresponds to a 2-split difference from the trusted tree, and all top methods converge at this floor. The evidential value is therefore not that the ensemble distinguishes itself from MSA+ML, but that an alignment-free method with no k selection and no alignment step (~15s total vs ~75s for MAFFT alone) reaches the community optimum on real deep-divergence indel-rich data.

2. **The k parameter must scale with sequence length, and the multi-k ensemble provides this automatically.** The gene-scale default (k=5) saturates at genome scale (only 1,024 possible 5-mers over 17–4,600 kb of sequence), producing near-random trees (fish nRF=0.545; E. coli nRF=0.875). Single larger k values recover signal on fish (k=9: 0.091; k=11: 0.091), yet the multi-k average — which includes the saturated k=5 and k=7 matrices — still achieves the best result (0.045). A plausible but unverified mechanism (this is a single dataset, and the explanation is post hoc) is that the informative scale-specific variation in the k=7–11 matrices dominates the topology while the ensemble averages out scale-specific noise. The multi-k ensemble thereby removes the need for length-aware k selection across a >2,000-fold range of sequence lengths (500 bp genes to 17 kb mitogenomes).

3. **E. coli/Shigella delineates an applicability boundary: shallow divergence dominated by recombination and gene flux.** On the 27-genome dataset (within-species divergence, recombination and gene gain/loss dominate over substitution), k-mer cosine plateaus at nRF=0.50 for all k≥9 and all multi-k variants, and kmacs (k=3; run left overnight, not instrumented for timing) performs no better (nRF=0.542), whereas Mash reaches 0.208 (our run) / 0.12 (published), and anchor-based methods reach 0.08 (published andi, Co-phylog, phylonium). At within-species identity, the frequency-vector difference between genomes is dominated by the conserved background, while anchor- and sketch-based distances (shared k-mer fractions, average common substrings) remain sensitive to the small divergent fraction — note that kmacs, despite being substring-based, also fails to exploit this fraction at its gene-scale k=3 configuration. We conclude that k-mer cosine distances are the method of choice for moderate-to-deep divergence (gene families across genera to phyla, organelle genomes, the simulated 5%-substitution regime), and recommend anchor-based or MinHash methods for strain-level bacterial collections — a boundary that complements the 16S structured-gene boundary reported above.

### Gene-length validation on real DNA: 13 mitochondrial loci

The whole-genome benchmarks above leave the simulated gene-length regime (L=500 bp) untested on real DNA. We therefore extracted the 13 mitochondrial protein-coding genes (ND1–6, ND4L, COX1–3, ATP6/8, CYTB; 168–1,866 bp, 5 of 13 loci within the 500–1,000 bp band of the simulation benchmarks) from the same 25 fish genomes and ran per-gene benchmarks against a MAFFT+FastTree2 anchor, scoring ETE3 nRF both against the Fischer et al. whole-mitogenome reference tree and — more informatively — as direct Fusang-vs-FastTree2 topological agreement (Supplementary Table S17).

Per-gene gene-tree discordance dominates reference-relative scores: even the MAFFT+FastTree2 anchor differs from the whole-mitogenome reference by nRF=0.332±0.175, so no per-gene tree can be expected to match the full-genome phylogeny. Against this background, the multi-k ensemble (mean of k=5,7,9 cosine distances) agrees with the MSA+ML anchor at nRF=0.451±0.127 across the 13 loci. Reference-relative scores show the same k-dependence as in simulation: the multi-k ensemble (0.535±0.113) and single k=7 (0.532±0.136) perform comparably, and both outperform the gene-scale default k=5 (0.601±0.107). Agreement is length-dependent: loci ≥1 kb show the closest Fusang–FastTree2 agreement (mean nRF=0.37, n=5), while the three shortest loci (ATP8 168 bp, ND4L 297 bp, ND3 350 bp; mean nRF=0.59) are poorly resolved by both method families. Because mitochondrial genes evolve under a low-indel regime (per-locus alignment gap fractions 0–1.7%), this benchmark complements the indel-rich whole-genome test of Table 12: it probes the low-indel gene-length regime on real DNA, and its interpretation is agreement with the MSA+ML reference pipeline rather than absolute accuracy. We note that this is a single locus set from one mitogenome clade; broader gene-family validation remains future work.

### Parameter stability and reproducibility

Automated adaptive parameter selection was validated across dataset scales with 100% reproducibility: all 5 stability repeats at n=20/50/100/200/500 returned identical k,gap selections and nRF within 0.005 units of each other (Supplementary Table S5). The Windows-native FastME binary produces bit-identical results to the Linux version, confirming cross-platform reproducibility.

---

## DISCUSSION

### K-mer frequency vectors for indel-rich phylogenetics

The central contribution of this work is the systematic evaluation of k-mer frequency vector cosine distances for phylogenetic inference under indel-rich conditions, delivered as an open-source tool (Fusang: Tardigrade Edition) with automatic parameter selection, a multi-k ensemble requiring no manual k choice, and an empirically delimited domain of applicability. Spaced k-mers were introduced by PatternHunter in 2002 [30] for sequence alignment and have been validated in protein classification, metagenomics, and genome assembly — yet their application to phylogenetic inference remained almost entirely unexplored. Fusang: Tardigrade Edition investigates both spaced and contiguous k-mer patterns within the cosine distance framework, finding that:

1. **Cosine distance is the primary accuracy driver; spaced patterns and NNI refinement provide significant secondary gains.** On the 130-seed INDELible benchmark, the like-for-like NJ comparison shows spaced k=5,gap2 significantly outperforming contiguous k=5 cosine (TRUE-relative 0.132 vs 0.136, 68/130 seeds, p=7.97×10⁻³, d=0.234) — reversing our earlier underpowered 27-seed FT2-relative result (p=0.26) — and BME-scored NNI refinement adds a further significant gain (0.122 vs 0.132, p=2.1×10⁻¹²). Spaced patterns show their clearest benefits where gap content is high: the best protein-family configuration is spaced (k=4,gap1, nRF=0.239, Table 10; spaced vs contiguous p=0.0078), and on indel data a wider gap (gap3) improves over gap2 by 10.5% (Supplementary Table S2). The systematic evaluation of k-mer frequency vector cosine distances across patterns and sequence types is the core methodological contribution.

2. **MinHash and k-mer cosine both degrade severely under indels; absolute error remains marginally lower for cosine distances.** On indel-rich data, Mash (MinHash Jaccard, 30 seeds) produces trees with nRF=0.762 ± 0.036, while k-mer cosine degrades to nRF=0.742 ± 0.045 (both vs true tree). Relative degradation is comparable (1.93× vs 1.97× from clean baselines). The mechanistic asymmetry — MinHash sketches discard positional information entirely while frequency vectors retain distributional signal — is a plausible explanation for the small absolute gap, but the practical difference at this indel rate is minor.

3. **The distance metric and k-mer representation must fit the data.** On clean data, MinHash Jaccard (nRF=0.394 ± 0.045, 30-seed benchmark) slightly underperforms k-mer cosine (nRF=0.376 ± 0.045); on indel data, both degrade severely (0.762 vs 0.742). This underscores that "alignment-free" is not a monolithic category: method and parameter selection should be guided by the expected evolutionary process and data characteristics.

### From Fusang v1 to Tardigrade Edition: an architectural evolution

The Fusang lineage illustrates a methodological evolution from specialized to general-purpose phylogenetic inference. Our v1 [6] established proof-of-concept: alignment-free phylogenetics can work, but was limited to 4–40 taxa and required pre-trained deep learning models. The Tardigrade Edition represents a complete re-architecture:
- **Representation**: neural network features → spaced k-mer frequency vectors
- **Scalability**: 40 taxa maximum → 10,000+ taxa
- **Speed**: GPU-dependent inference → CPU-only, minutes for thousands of taxa
- **Deployment**: Docker/cloud requirement → single binary + Python script
- **Indel handling**: no explicit indel modeling → hypothesized robustness from spaced sampling (empirically modest at indel rate 0.02; see Results)

This evolution demonstrates that feature engineering — specifically, **k-mer frequency vector cosine distances** — can match or exceed learned representations for phylogenetic inference while providing order-of-magnitude improvements in speed, scalability, and accessibility. The key insight is that k-mer frequency vectors with cosine distance preserve sufficient phylogenetic signal for accurate distance-based tree inference, while spaced sampling is hypothesized to provide additional robustness against indel-induced length variation (an advantage that is empirically modest at the tested indel rate).

### The role of pipeline simplicity in phylogenetic accuracy

The discovery that the simplified pipeline (nRF=0.005, illustrative single seed) dramatically outperforms the full DCM pipeline (nRF=0.388) at n≤200 has important implications. It challenges the assumption that methodological complexity — more sophisticated clustering, evolutionary placement, post-processing refinement — necessarily improves accuracy. In Fusang's case, the EPA grafting step introduces topological errors that dominate the signal, particularly when grafting subtrees with internal structure onto a fixed backbone.

This finding is consistent with a broader pattern in computational biology: simpler models often outperform complex ones when the underlying signal is weak or noisy. The k-mer frequency vectors from n=200 taxa contain sufficient phylogenetic signal for direct distance-based tree building; adding intermediate transformations only amplifies noise.

We therefore recommend the simplified pipeline as the default for n≤500 and DCM for larger datasets, where it provides a modest accuracy advantage (n=1000: DCM nRF=0.083 ± 0.005, 9 seeds 401–409, vs simplified nRF=0.102 on seed 401) in addition to the scalability needed when pairwise distance computation (scaling as O(n²)) becomes the computational bottleneck.

### Limitations and future work

Several limitations of the current study warrant discussion.

**Scale-dependent accuracy**: On clean (no-indel) data at large scales, MSA-based methods maintain a clear advantage. A 30-seed benchmark at n=500 shows Fusang nRF=0.119 ± 0.011 vs FastTree2 nRF=0.093 ± 0.013 (Cohen's d=1.77, Wilcoxon p<0.001). At n=1000, the gap widens: Fusang nRF=0.115 ± 0.011 vs FastTree2 nRF=0.091 ± 0.010 (Cohen's d=2.81, Wilcoxon p<0.001). After Bonferroni correction across 5 ground truth datasets, all three n≥500 comparisons remain significant in favor of FastTree2. At n=200, no significant difference is detected. This is expected: on clean data without indels, alignment-based ML methods benefit from full positional information. Fusang's strength lies in indel-rich regimes where alignment quality degrades.

On indel-rich data at n=1000, Fusang achieves nRF=0.037 ± 0.006 (30 seeds, vs FastTree2 reference tree), indicating numerically close nRF values with this MSA+ML method on this dataset (FT2-relative nRF; TRUE-relative accuracy may differ). **Readers should note that this nRF is computed against the FastTree2 reference tree (FT2-relative), not the true simulated ground truth; the actual accuracy vs true tree may differ.** The simplified pipeline is the default for n≤500; for n>500, DCM+EPA provides essential scalability and a small accuracy advantage (direct comparison at n=1000: DCM nRF=0.083 ± 0.005, 9 seeds 401–409, vs simplified nRF=0.102 on seed 401; an earlier threshold of 2,000 forced the simplified pipeline at that scale, degraded accuracy, and was reverted). A systematic multi-seed accuracy comparison across the n=500–2,000 transition zone remains future work.

**Simulated-to-real transfer**: While 16S rRNA validation (73 analysed taxa, 13 phyla under corrected NCBI Taxonomy annotation; Supplementary Note S17) confirms that Fusang detects biological phylogenetic signal — same-rank tree distances are significantly reduced at every level (phylum 31.4% to genus 78.6%, all p<1×10⁻⁴) and 3/7 assessable phyla plus 10/14 assessable orders are recovered as monophyletic — accuracy on this structured RNA gene is lower than MSA-based methods (Fusang vs taxonomy nRF=0.414 vs FastTree2 vs taxonomy nRF=0.343, RF/(2n−6) normalisation). The nRF=0.371 between Fusang and the rebuilt FastTree2 tree reflects the substantial divergence between k-mer frequency and alignment-based phylogenetic inference on this sequence type. 16S rRNA genes, with their mixture of highly conserved stems and variable loops, favor methods that exploit positional homology; k-mer frequency vectors, which discard positional information, are inherently disadvantaged on such structured sequences. The AFproject SwissTree protein benchmark (Table 10) provides additional cross-domain validation: k-mer frequency methods achieve mean nRF=0.239 on real protein gene trees (11 families, 21–159 taxa; median alignment gap content 50%, range 11–81%), while Co-phylog's best configuration (halfctx=11) achieves nRF=0.361 (paired Wilcoxon p=0.006, Cohen's d=1.32 vs Fusang k=4,gap1) and its default (halfctx=5) achieves nRF=0.433 (p=0.014, d=1.08) — confirming that k-mer approaches transfer robustly from simulated DNA to real, strongly indel-rich protein data. On BAliBASE v3.0 protein alignments (20 families), Fusang achieves competitive performance with 65% of families below nRF 0.5 (median nRF=0.45; Supplementary Table S7). Real genome-scale validation on AFproject community datasets (Table 12) shows the multi-k ensemble reaching the community-benchmark ceiling on deep-divergence indel-rich mitochondrial genomes (nRF=0.045, matched by all top methods), but also delineates a boundary: on shallow-divergence E. coli/Shigella genomes where recombination and gene flux dominate, k-mer cosine (nRF=0.50) is outperformed by Mash (0.12–0.21) and anchor-based methods (published 0.08). The gene-length DNA regime of the simulation benchmarks is now covered by the 13-locus mitochondrial benchmark above (168–1,866 bp), which shows moderate Fusang–FastTree2 agreement (multi-k nRF=0.451±0.127) under a low-indel regime with gene-tree discordance as the dominant error source; real gene-length DNA with documented indel-rich evolution and trusted trees remains an open validation target, with candidate data types including viral quasispecies collections and rapidly evolving gene families with curated reference trees.

**Fixed k and gap**: The current implementation uses a static k and gap for the entire dataset. The multi-k distance ensemble (Table 4, 130-seed INDELible benchmark) achieves accuracy comparable to the best single contiguous configuration (ensemble nRF=0.128 vs k=7 contiguous nRF=0.127, n.s.), providing robust performance without manual k selection. However, it does not exceed the best single-k configuration or the full spaced+NNI pipeline (nRF=0.122) at gene scale, suggesting limited complementary signal beyond what a well-chosen single k captures at L=500 bp. Further optimization, such as weighted fusion, inclusion of spaced k-mers in the ensemble, or per-cluster parameter selection, may yield additional gains. The optimal set of k values and fusion weights has not been systematically explored.

**Distance metric exploration**: Cosine distance and Jensen-Shannon divergence represent points in a larger space of possible metrics on k-mer frequency vectors. Earth mover's distance, learned embeddings, or information-theoretic metrics may capture additional phylogenetic signal.

**L3 protocol and portability**: The definitive L3 (MAFFT+FastTree2) comparison is now based on a complete 30-seed run executed on Linux (WSL2 Ubuntu-24.04) with an explicitly nucleotide FastTree2 protocol (`-nt -gtr -nosupport`); an earlier 5-seed Windows run was affected both by MAFFT instability (empty alignments on 25/30 seeds) and by an inadvertently omitted `-nt` flag (protein-mode distance settings on DNA input), and is superseded (Supplementary Note S8). Both L1 (nRF=0.583) and L3 (nRF=0.601) remain substantially distant from the true tree under this deliberately harsh coalescent protocol (~58–60% bipartition error). The Linux requirement itself is a practical portability constraint: native-Windows MAFFT proved unreliable at this scale.

**Mash multi-seed benchmark**: The Mash benchmark (Table 2) is based on 30 seeds (this work, WSL2 Ubuntu 24.04, Mash v2.3). Mash shows substantial (though not random) degradation on indel-rich data (nRF=0.762 ± 0.036, 1.93× degradation from clean baseline), while k-mer cosine distances degrade comparably (1.97×).

**Sequence length**: All simulated benchmarks used a fixed sequence length of L=500 bp, representative of typical gene-length sequences. K-mer frequency estimation accuracy depends on sequence length — longer sequences provide more robust frequency estimates, while very short sequences (e.g., L<200 bp, typical of amplicon data) may yield noisier k-mer profiles. Cross-validation on the AFproject SwissTree dataset (protein sequences, 111–576 amino acids; Table 10) confirms that Fusang's performance transfers across a range of sequence lengths, but systematic evaluation at extreme lengths (L=100–200 bp and L>2000 bp) remains future work.

**Boundary classifier data independence**: The 88 E2E test scenarios were generated using simulation parameters distinct from the training data, but all scenarios used the same simulation engine (INDELible). The 100% accuracy, while statistically supported (Wilson CI [0.958, 1.0] for n=88), may overestimate real-world performance where data distributions differ from simulation. Internal cross-validation on the training set (n=4,073, 5-fold CV, ROC-AUC=0.84) confirms the classifier generalizes within its simulation domain. Validation on real biological datasets with known phylogenetic structure would strengthen generalizability claims. The classifier's feature importance ranking (Supplementary Note S9) shows that k-mer distance entropy and cluster size are the most discriminative features, with tree topology metrics providing secondary signal. A systematic feature ablation study quantifying each feature group's causal contribution to classification accuracy remains future work.

**Current status of n=500/n=1000 benchmarks**: Multi-seed benchmarking at n=500 and n=1000 (30 seeds each) has been completed for both clean and indel-rich data. The n=1000 indel benchmark (30 seeds, sub=0.05, indel=0.02) shows Fusang nRF=0.037 ± 0.006 vs the FastTree2 reference, indicating high topological similarity to the alignment-based reconstruction at scale (FT2-relative; TRUE-relative accuracy is unavailable at this scale). We note that at n=1000, the FastTree2 reference should be interpreted as "alignment-based consensus" rather than "gold standard," given that IQ-TREE2 — the other MSA-based ML method — did not complete at this scale (10/10 seeds timed out at 24 hours). A multiple comparison correction across all 5 ground truth datasets confirms that at n≥500, FastTree2 significantly outperforms Fusang (p<0.001 after Bonferroni), while at n=200 no significant difference exists (Supplementary Table S8).

**Comparison with classical alignment-free methods**: We have completed systematic head-to-head comparison with four classical alignment-free methods — andi [25], Co-phylog [26], kmacs [12], and Skmer [29] — across simulated DNA (Table 7, 130 INDELible seeds, TRUE-relative) and real protein data (Table 10, AFproject SwissTree, 11 gene families). K-mer frequency methods consistently outperform the alternatives: Co-phylog's context-matching approach fails under indels (DNA TRUE-relative nRF=0.416, 3.4× worse than Fusang, 0/130 wins, paired d=9.34) and on protein data (best config halfctx=11 nRF=0.361 vs Fusang nRF=0.239, d=1.32, p=0.006); kmacs (best of k∈{3,5,10}) is significantly less accurate under indels (TRUE-relative nRF=0.221 vs Fusang 0.122, 0/130 wins, p=4.4×10⁻²³, paired d=3.84); while andi (suffix-array anchors) and Skmer (genome-skim coverage model, division-by-zero on 500-bp inputs at k=31 and k=21) are structurally inapplicable to gene-length sequences by design. We note that Co-phylog was tested at k=19 with default half-context on DNA; a systematic scan across parameter space may reveal configurations where Co-phylog performs better, though its fundamental reliance on context-matching makes it intrinsically vulnerable to indels. These results position k-mer frequency vectors with cosine distance as among the most robust alignment-free approaches for gene-length phylogenetic inference across sequence types. Two method families remain untested — FFP [31], which pioneered feature frequency profiles for whole-genome comparison, and the broader spaced-word frequency family [13] (both were outside the scope of the current benchmark harness); systematic comparison under indel-rich conditions represents future work.

Future work will explore: (1) optimized k-mer sets and fusion weights for the multi-k ensemble, potentially including spaced k-mers in the ensemble and length-aware k selection; (2) validation on real gene-length DNA families with documented indel-rich evolution and trusted trees, extending the real-data validation beyond protein families, mitochondrial genomes and loci (low-indel), and bacterial genomes; (3) Mash benchmark across a range of indel rates (0.005–0.05) to quantify the MinHash collapse threshold (the n=30 benchmark at indel=0.02 is now complete, Table 2); (4) boundary classifier validation on real biological datasets with known phylogenetic structure; (5) integration as a rapid exploratory analysis module within existing phylogenetic pipelines; (6) application to metagenomic and single-cell datasets where alignment is particularly challenging; and (7) quartet-based DCM assembly to overcome the EPA grafting bottleneck at large scales.

### Practical recommendations

For practitioners, our results suggest the following guidelines:
- **Small-to-medium datasets (n≤1000) with expected indel rates above 0.01**: Consider Fusang with the multi-k ensemble (via `benchmark_multik_ensemble.py`) as a first-pass analysis. On the full 30-seed coalescent benchmark the ensemble outperforms MAFFT+FastTree2 (nRF=0.583±0.044 vs 0.601±0.055, paired Wilcoxon p=0.024) while requiring no alignment step (~2× faster per seed). Both L1 and L3 remain substantially distant from the true tree under this deliberately harsh protocol (nRF≈0.58–0.60), so outputs should be treated as rapid first-pass topologies. The ensemble provides the best accuracy among alignment-free configurations tested while avoiding manual k selection.
- **Clean substitution-only data (no indels)**: Mash (MinHash Jaccard, k=21) does not provide superior accuracy to k-mer cosine distances. The 30-seed benchmark (Table 2) shows Mash nRF=0.394 ± 0.045 vs Fusang nRF=0.376 ± 0.045 on clean data. K-mer cosine distances (Fusang) are the recommended alignment-free approach for both clean and indel-rich data.
- **Large datasets (n>1000)**: MSA-based methods remain preferred for accuracy; Fusang provides a valuable speed-accuracy trade-off for rapid exploratory analysis
- **Deep-divergence genome-scale data (organelle genomes, divergent gene families)**: The multi-k cosine ensemble reaches the community-benchmark ceiling on real mitochondrial genomes (Table 12; nRF=0.045, matched by MSA+ML and all top published methods) without alignment and without k selection; the gene-scale default k=5 saturates at genome length — use the multi-k ensemble (via `benchmark_multik_ensemble.py`) instead of any single small k
- **Shallow strain-level divergence (recombination-dominated)**: For within-species bacterial collections, anchor-based methods (andi, Co-phylog, kmacs) or Mash are better suited than k-mer cosine distances (Table 12)
- **Indel-rich data at any scale**: Fusang's alignment-free nature provides robustness that is not available from MSA-based methods, regardless of scale. Under indels, MinHash-based approaches degrade severely (nRF≈0.76, ~76% bipartition error on gene-length sequences) and should be avoided for this data type.
- **Computational constraints**: Fusang requires no alignment step and runs in seconds to minutes on a single CPU core, making it suitable for rapid iteration during exploratory analysis

---

## SUPPLEMENTARY MATERIAL

Supplementary Data are available at NAR Online.

**Supplementary Figure S1.** Full k-mer parameter grid search. (A–E) nRF heatmaps (k=4–8 × gap 0–3) at n=50, 100, 200, 500, and 1000; the best configuration per scale is marked with an asterisk. (F) nRF vs dataset size for the default k=5,gap2 and the best per-n configuration.

**Supplementary Figure S2.** Dimensionality vs accuracy. (A–E) nRF vs feature vector dimension (4^k) for each gap pattern at n=50–1000, showing 1024-dim (k=5,gap2) outperforming 65536-dim (k=8, contiguous). (F) nRF vs dataset size at k=5 for each gap pattern.

**Supplementary Figure S3.** DCM degradation trace and analysis. (A) Stage-by-stage nRF trace on a single illustrative seed (n=200, seed=42, clean): simplified pipeline 0.005, TF-IDF 0.030, DCM+NJ recovery 0.005, DCM+FastME 0.013, full DCM+EPA 0.388. (B) Simplified pipeline multi-seed validation (10 seeds, 0.014±0.003) vs full DCM. (C) Schematic of the DCM pipeline with EPA grafting identified as the primary error source. (D) nRF stability across scales including the DCM regime (clean coalescent, 5 seeds per scale).

**Supplementary Figure S4.** 130-seed indel benchmark: distributions and outlier analysis. (A) nRF violins for all 130 paired seeds. (B) Per-seed Fusang vs FastTree2 scatter with the pre-specified nRF>0.3 exclusion rule (no seed exceeded it in this INDELible dataset). (C) Sorted paired differences for all 130 seeds. (D) Histogram of paired differences (130 seeds; Wilcoxon p=3.5×10⁻⁷, paired d=+0.51, Fusang better in 91/130).

**Supplementary Figure S5.** Real 16S rRNA validation (73 taxa, corrected NCBI Taxonomy annotation). (A) Taxonomic-rank signal (same-rank distance reduction) for the Fusang and rebuilt FastTree2 trees, per rank (Table 9 values). (B) Fusang tree (k=5,gap2) with phylum-coloured tips. (C) Rebuilt FastTree2 tree (MAFFT 7.526 --localpair + FastTree2 -nt -gtr) with identical colouring for direct comparison. (D) Recovery against NCBI taxonomy: monophyletic-group recovery (Fusang 3/7 phyla, 10/14 orders; FastTree2 6/7 phyla, 11/14 orders) and nRF vs a taxonomy-derived reference tree (0.414 vs 0.343, RF/(2n−6); 0.547 vs 0.453, RF/(s₁+s₂)).

**Supplementary Figure S6.** Effect size analysis across benchmarks. (A) Cohen's d with 95% bootstrap CI for Fusang vs FastTree2 (indel n=200, 130 valid seeds; clean n=500 and n=1000, 30 seeds; orientation: FT2 − Fusang, positive = Fusang better). (B) Paired Cohen's d against alignment-free competitors and ablations (Co-phylog, kmacs, IQ-TREE2, SwissTree, multi-k L1 vs L0). (C) BAliBASE v3.0 per-family nRF vs FastTree2 (20 protein families; median 0.445; 65% below 0.5). (D) Bootstrap distribution of the 130-seed paired mean difference (10,000 resamples).

**Supplementary Figure S7.** Multi-k distance ensemble: per-seed analysis (130-seed INDELible benchmark, n=200, indel=0.02, TRUE-relative). (A) Per-seed scatter of the NJ-only spaced k=5,gap2 configuration vs the ensemble (k=5,7,9 average). (B) Configuration comparison across single-k k=5/7/9, ensemble, spaced NJ-only, and the full pipeline (mean ± SD). (C) Sorted per-seed differences: ensemble vs contiguous k=5 (83/130 wins, p=1.6×10⁻⁷, d=0.49). (D) Sorted per-seed differences: ensemble vs full pipeline (45/130 wins, p=2.2×10⁻⁴, d=0.36), showing the ensemble does not exceed the default pipeline at L=500 bp.

**Supplementary Table S1.** NCBI GenBank accessions for 74 16S rRNA sequences.

**Supplementary Table S2.** Complete k/gap parameter scan on indel data (indel rate=0.02): nRF for k=4–8 × gap=none/gap1–gap4 at n=100/200.

**Supplementary Table S3.** DCM degradation step-by-step: nRF at each pipeline stage (n=200, seed=42, clean data).

**Supplementary Table S4.** 130-seed benchmark raw data: per-seed nRF for Fusang and FastTree2 at n=200, indel rate=0.02.

**Supplementary Table S5.** Adaptive parameter stability: 5-repeat validation of auto-selected k,gap across n=20/50/100/200/500.

**Supplementary Table S6.** Multi-seed benchmark at n=500 (30 seeds) and n=1000 (30 seeds): all conditions complete including n=1000 indel data. All raw nRF values and statistical tests are provided for reproducibility.

**Supplementary Table S7.** BAliBASE v3.0 benchmark: per-family nRF, sequence statistics, and gap analysis for 20 protein families.

**Supplementary Table S8.** Multiple comparison correction: raw and adjusted p-values for 5 ground truth dataset comparisons (Fusang vs FastTree2). Both Bonferroni (α=0.01) and Benjamini-Hochberg FDR corrections applied.

**Supplementary Table S9.** Multi-k ensemble benchmark: per-seed nRF for Fusang (k=5,gap2; NJ-only and full pipeline), individual contiguous k-mers (k=5,7,9), and ensemble average (n=200, indel rate=0.02, 130 seeds 100–229, INDELible; TRUE-relative and FT2-relative).

**Supplementary Table S10.** Alignment-free competitor comparison: per-seed nRF for Co-phylog (k=19), kmacs (k=3,5,10), k-mer cosine (k=5, k=7 contiguous), Fusang (k=5,gap2 spaced; NJ-only and full pipeline), and the multi-k ensemble on n=200 indel data (130 seeds 100–229, INDELible; FT2-relative and TRUE-relative).

**Supplementary Table S11.** SwissTree gene tree benchmark: per-family nRF for all methods, sequence statistics, and reference tree properties (AFproject standard).

**Supplementary Table S12.** SwissTree per-family nRF results for all methods (k-mer cosine, spaced k-mers, Co-phylog configurations).

**Supplementary Table S13.** L3 pipeline-level validation: per-seed nRF for L0 (single-k NJ), L1 (multi-k ensemble NJ), and L3 (MAFFT+FastTree2) against the true simulated coalescent tree (n=200, L=1000 bp, indel=0.02, 30 seeds, all levels complete). L3 executed on Linux (WSL2 Ubuntu-24.04, MAFFT --auto + FastTree2 -nt -gtr -nosupport); supersedes the earlier 5-seed Windows run (Supplementary Note S8).

**Supplementary Table S14.** Mash vs spaced k-mer benchmark: per-seed nRF for Fusang (k=5,gap2, NJ) on clean and indel-rich data (30 seeds each vs true tree), plus Mash 30-seed results (k=21, MinHash Jaccard, NJ, this work).

**Supplementary Table S15.** Boundary classifier (RF V4b) E2E validation: all 88 individual scenario results including tree type, expected decision, predicted probability, and correctness.

**Supplementary Table S16.** Real genome benchmarks (AFproject community datasets): complete per-method results for fish mtDNA (25 mitochondrial genomes, Fischer et al. 2013 reference tree) and E. coli/Shigella (27 genomes, Skippington & Ragan 2011 reference tree), including NCBI accessions, k-scan results, multi-k variants, kmacs/Mash/Co-phylog runs at working and AFproject-published configurations, MAFFT+FastTree2 (fish), alignment gap-content quantification, and SwissTree per-family gap fractions; ETE3 unrooted nRF scoring.

**Supplementary Table S17.** Gene-length real DNA benchmark: per-locus results for the 13 mitochondrial protein-coding genes (168–1,866 bp) from the 25 fish genomes — locus lengths, per-locus alignment gap fractions, ETE3 nRF vs the Fischer et al. 2013 whole-mitogenome reference tree, and Fusang-vs-FastTree2 agreement nRF, for the multi-k ensemble (mean of k=5,7,9 cosine distances), single k=5, single k=7, and the MAFFT+FastTree2 anchor.

**Supplementary Note S1.** INDELible simulation parameters and indel model details.

**Supplementary Note S2.** Gap optimality on indel data: mechanism and robustness analysis.

**Supplementary Note S3.** BAliBASE gap diagnosis: root causes of performance variation on protein sequences.

**Supplementary Note S4.** Alignment-free competitor method implementations: detailed description of Co-phylog (Python reimplementation based on source code), K-mer cosine baseline, and andi scale limitation analysis.

**Supplementary Note S5.** Fusang v1 (NAR 2023) detailed comparison: architecture, performance, and limitations.

**Supplementary Note S6.** Multiple comparison correction methodology and multi-k ensemble parameter selection.

**Supplementary Note S7.** SwissTree gene tree benchmark: per-family detailed results, reproduction instructions, and comparison with AFproject published benchmarks.

**Supplementary Note S8.** L3 validation methodology: coalescent protocol parameters (n=200, L=1000 bp, sub=0.05, indel=0.02, identical for the superseded 5-seed run and the 30-seed rerun), MAFFT+FastTree2 configuration (Linux rerun with corrected nucleotide protocol `-nt -gtr -nosupport`), true tree generation, Windows MAFFT instability analysis, and supersession of the earlier 5-seed Windows run (which omitted the -nt nucleotide flag).

**Supplementary Note S9.** Boundary classifier architecture: random forest feature set, training data composition (4,073 samples), and E2E test scenario generation parameters.

**Supplementary Note S10.** IQ-TREE2 benchmark methodology and results: complete timeout analysis at n≥500, n=200 indel comparison showing IQ-TREE2 1.9× worse than Fusang k-mer NJ (nRF=0.231 ± 0.032 vs 0.122 ± 0.024, 130 seeds, paired Wilcoxon p=4.5×10⁻²³, d=3.56), and analysis demonstrating that indel-induced alignment errors degrade even fixed-model (GTR) ML inference. The IQ-TREE2 benchmark used a fixed GTR model (`-m GTR`), not ModelFinder auto-selection, ensuring a fair model-class comparison with FastTree2 (also GTR-based).

**Supplementary Note S17.** 16S rRNA dataset annotation incident and remediation. The accession→organism mapping table (`CURATED_16S` in `download_16s_robust.py`) originally shipped with this dataset was erroneous at the genus/species/strain level for all 73 downloaded accessions. The sequences themselves are genuine GenBank 16S rRNA records — each was verified against NCBI as byte-identical (identity=1.0000 on spot checks) — so the error was in the labels, not the sequences. We re-annotated every accession against NCBI Taxonomy via its TaxID and recomputed every taxonomy-dependent analysis (Table 9 permutation tests, monophyly recovery, RF against the taxonomy-derived reference tree) under the corrected annotation; all values reported in the main text use this corrected annotation. The incident also surfaced a second, independent defect: the archived alignment (`real_16S_full_aligned.fasta`) was not genuinely aligned (cross-species pairwise identity fell below the 0.25 random baseline), so the alignment-based comparator was rebuilt end-to-end (MAFFT 7.526 `--localpair --maxiterate 1000` + FastTree2 `-nt -gtr -nosupport`). The full incident record, true dataset composition, and reproduction steps are provided in the reproducibility package (`real_data/real_16S_TRUE_label_incident.md`, `real_16S_monophyly.json`). The Fusang method itself is unaffected: it is alignment-free and never consumed the erroneous annotation. Earlier evaluations computed against the erroneous labels (including the rank-signal values of the archived analysis) are superseded.

---

## DATA AVAILABILITY

Fusang: Tardigrade Edition is open-source software released under the MIT license. Source code, documentation, pre-compiled FastME binaries (Windows x86-64, Linux x86-64), benchmark scripts, and all analysis code are available at:

- **GitHub**: https://github.com/zhanglknt/fusang-tardigrade
- **Zenodo**: https://doi.org/10.5281/zenodo.20746742 (archived source code, benchmark datasets, and supplementary materials)

All benchmark datasets — including 130-seed n=200 indel benchmark results, n=500 and n=1000 multi-seed data, 30-seed L3 pipeline validation (Linux), Mash benchmark, kmacs and Skmer competitor benchmarks (including the Skmer inapplicability record), AFproject real genome benchmarks (fish mtDNA and E. coli/Shigella, with NCBI accessions and full per-method results), E2E classifier validation, 73-taxon 16S rRNA dataset, BAliBASE v3.0 20-family results, and DCM degradation tracing data — are provided in the repository under `benchmarks/` and as Supplementary Data. A storage-corruption incident affected 21 competitor-benchmark input files (5 FASTA and 16 tree files, seeds 109–126); recovery proceeded via three independent validated paths (reconstruction from intact aligned FASTA, deterministic regeneration from the seeded simulation code — regenerated files validated byte-identical across the independent recovery paths — and FastTree2 reruns on intact alignments). All affected tables were recomputed: 53/54 affected per-seed nRF values match the pre-incident master CSV (one FastTree2 value differs by two bipartitions, 0.0787→0.0838). Corrupted originals are quarantined, and the recovery log with checksums is documented in the reproducibility package. INDELible simulation configuration files are included for full reproducibility.

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

1. Hadfield, J., Megill, C., Bell, S.M., Huddleston, J., Potter, B., Callender, C., Sagulenko, P., Bedford, T. and Neher, R.A. (2018) Nextstrain: real-time tracking of pathogen evolution. *Bioinformatics*, 34, 4121–4123. DOI: 10.1093/bioinformatics/bty407

2. Hug, L.A. et al. (2016) A new view of the tree of life. *Nat. Microbiol.*, 1, 16048. DOI: 10.1038/nmicrobiol.2016.48

3. Warnow, T. (1994) Some combinatorial optimization problems in phylogenetic tree reconstruction. *DIMACS Technical Report*, 94-53. [Technical report; DOI not available.]

4. Dessimoz, C. and Gil, M. (2010) Phylogenetic assessment of alignments reveals neglected tree signal in gaps. *Genome Biol.*, 11, R37. DOI: 10.1186/gb-2010-11-4-r37

5. Wong, K.M., Suchard, M.A. and Huelsenbeck, J.P. (2008) Alignment uncertainty and genomic analysis. *Science*, 319, 473–476. DOI: 10.1126/science.1151532

6. Wang, Z., Sun, J., Gao, Y., Xue, Y., Zhang, Y., Li, K., Zhang, W., Zhang, C., Zu, J. and Zhang, L. (2023) Fusang: a framework for phylogenetic tree inference via deep learning. *Nucleic Acids Res.*, 51, 10909–10923. DOI: 10.1093/nar/gkad805

7. Morgenstern, B., Zhu, B., Horwege, S. and Leimeister, C.-A. (2015) Estimating evolutionary distances between genomic sequences from spaced-word matches. *Algorithms Mol. Biol.*, 10, 5. DOI: 10.1186/s13015-015-0032-x

8. Diaz, N.N., Krause, L., Goesmann, A., Niehaus, K. and Nattkemper, T.W. (2009) TACOA: taxonomic classification of environmental genomic fragments using a kernelized nearest neighbor approach. *BMC Bioinformatics*, 10, 56. DOI: 10.1186/1471-2105-10-56

9. Vinga, S. and Almeida, J. (2003) Alignment-free sequence comparison — a review. *Bioinformatics*, 19, 513–523. DOI: 10.1093/bioinformatics/btg005

10. Zielezinski, A., Vinga, S., Almeida, J. and Karlowski, W.M. (2017) Alignment-free sequence comparison: benefits, applications, and tools. *Genome Biol.*, 18, 186. DOI: 10.1186/s13059-017-1319-7

11. Bernard, G., Chan, C.X., Chan, Y.-b., Chua, X.Y., Cong, Y., Hogan, J.M., Maetschke, S.R. and Ragan, M.A. (2019) Alignment-free inference of hierarchical and reticulate phylogenomic relationships. *Brief. Bioinform.*, 20, 426–435. DOI: 10.1093/bib/bbx067

12. Leimeister, C.-A. and Morgenstern, B. (2014) kmacs: the k-mismatch average common substring approach to alignment-free sequence comparison. *Bioinformatics*, 30, 2000–2008. DOI: 10.1093/bioinformatics/btu331

13. Leimeister, C.-A., Boden, M., Horwege, S., Lindner, S. and Morgenstern, B. (2014) Fast alignment-free sequence comparison using spaced-word frequencies. *Bioinformatics*, 30, 1991–1999. DOI: 10.1093/bioinformatics/btu177

14. Zielezinski, A. et al. (2019) Benchmarking of alignment-free sequence comparison methods. *Genome Biol.*, 20, 144. DOI: 10.1186/s13059-019-1755-7

15. Huson, D.H., Nettles, S.M. and Warnow, T.J. (1999) Disk-covering, a fast-converging method for phylogenetic tree reconstruction. *J. Comput. Biol.*, 6, 369–386. DOI: 10.1089/106652799318337

16. Berger, S.A., Krompass, D. and Stamatakis, A. (2011) Performance, accuracy, and web server for evolutionary placement of short sequence reads under maximum likelihood. *Syst. Biol.*, 60, 291–302. DOI: 10.1093/sysbio/syr010

17. Lefort, V., Desper, R. and Gascuel, O. (2015) FastME 2.0: a comprehensive, accurate, and fast distance-based phylogeny inference program. *Mol. Biol. Evol.*, 32, 2798–2800. DOI: 10.1093/molbev/msv150

18. Gascuel, O. (1997) BIONJ: an improved version of the NJ algorithm based on a simple model of sequence data. *Mol. Biol. Evol.*, 14, 685–695. DOI: 10.1093/oxfordjournals.molbev.a025808

19. Fletcher, W. and Yang, Z. (2009) INDELible: a flexible simulator of biological sequence evolution. *Mol. Biol. Evol.*, 26, 1879–1888. DOI: 10.1093/molbev/msp098

20. Price, M.N., Dehal, P.S. and Arkin, A.P. (2010) FastTree 2 — approximately maximum-likelihood trees for large alignments. *PLoS ONE*, 5, e9490. DOI: 10.1371/journal.pone.0009490

21. Katoh, K. and Standley, D.M. (2013) MAFFT multiple sequence alignment software version 7: improvements in performance and versatility. *Mol. Biol. Evol.*, 30, 772–780. DOI: 10.1093/molbev/mst010

22. Kozlov, A.M., Darriba, D., Flouri, T., Morel, B. and Stamatakis, A. (2019) RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference. *Bioinformatics*, 35, 4453–4455. DOI: 10.1093/bioinformatics/btz305

23. Minh, B.Q., Schmidt, H.A., Chernomor, O., Schrempf, D., Woodhams, M.D., von Haeseler, A. and Lanfear, R. (2020) IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. *Mol. Biol. Evol.*, 37, 1530–1534. DOI: 10.1093/molbev/msaa015

24. Ondov, B.D., Treangen, T.J., Melsted, P., Mallonee, A.B., Bergman, N.H., Koren, S. and Phillippy, A.M. (2016) Mash: fast genome and metagenome distance estimation using MinHash. *Genome Biol.*, 17, 132. DOI: 10.1186/s13059-016-0997-x

25. Haubold, B., Klötzl, F. and Pfaffelhuber, P. (2015) andi: Fast and accurate estimation of evolutionary distances between closely related genomes. *Bioinformatics*, 31, 1169–1175. DOI: 10.1093/bioinformatics/btu815

26. Yi, H. and Jin, L. (2013) Co-phylog: an assembly-free phylogenomic approach for closely related organisms. *Nucleic Acids Res.*, 41, e75. DOI: 10.1093/nar/gkt003

27. Lunter, G., Miklós, I., Drummond, A., Jensen, J.L. and Hein, J. (2005) Bayesian coestimation of phylogeny and sequence alignment. *BMC Bioinformatics*, 6, 83. DOI: 10.1186/1471-2105-6-83

28. Cartwright, R.A. (2009) Problems and solutions for estimating indel rates and length distributions. *Mol. Biol. Evol.*, 26, 473–480. DOI: 10.1093/molbev/msn275

29. Sarmashghi, S., Bohmann, K., Gilbert, M.T.P., Bafna, V. and Mirarab, S. (2019) Skmer: assembly-free and alignment-free sample identification using genome skims. *Genome Biol.*, 20, 34. DOI: 10.1186/s13059-019-1632-4

30. Ma, B., Tromp, J. and Li, M. (2002) PatternHunter: faster and more sensitive homology search. *Bioinformatics*, 18, 440–445. DOI: 10.1093/bioinformatics/18.3.440

31. Sims, G.E., Jun, S.R., Wu, G.A. and Kim, S.H. (2009) Alignment-free genome comparison with feature frequency profiles (FFP) and optimal resolutions. *Proc. Natl. Acad. Sci. USA*, 106, 2677–2682. DOI: 10.1073/pnas.0813249106

---

## FIGURE LEGENDS

**Figure 1.** Indel robustness advantage. (A) nRF vs indel rate for Fusang (simplified pipeline, k=5,gap2) and FastTree2 at n=200 (Table 3 values). (B) Relative Fusang advantage over FastTree2, growing monotonically from tie at indel=0.005 to 13.2% at indel=0.05. (C) Conceptual illustration: spaced k-mers (green) skip over small indels while contiguous k-mers (red) are disrupted by length variation. (D) Method comparison at indel rate=0.02 (nRF vs true simulated tree, mean ± SD): Fusang (0.122±0.024, 130 seeds), FastTree2 (0.134±0.022, 130 seeds), IQ-TREE2 GTR (0.231±0.032, 130 seeds), Mash (0.762±0.036, 30 seeds).

**Figure 2.** 130-seed statistical benchmark. (A) Violin plots of nRF distributions for the Fusang simplified pipeline and FastTree2 on n=200 indel data (indel rate=0.02; all 130 seeds valid, none exceeding the nRF>0.3 threshold), with per-seed points overlaid. (B) Sorted per-seed paired differences (FT2 − Fusang; positive = Fusang better); observed range [−0.038, 0.069], Fusang better in 91/130 seeds. (C) Cumulative win rate across the seed set, converging to 91/130 (70.0%). (D) Histogram of paired differences (FT2 − Fusang); Wilcoxon p=3.5×10⁻⁷, paired Cohen's d=+0.51.

**Figure 3.** Pipeline ablation and DCM degradation. (A) Step-by-step nRF tracing from the simplified pipeline through DCM stages (n=200 indel, seed=42: simplified 0.008 → FastME 0.015 → distances.py 0.013 → +TF-IDF 0.013 → full DCM+EPA 0.251; multi-seed simplified-pipeline accuracy is reported in Table 1). (B) Schematic of simplified vs DCM+EPA pipeline architectures; the adaptive switch is controlled by SIMPLE_THRESHOLD=500. (C) EPA grafting error illustration: subtrees with internal structure grafted onto a fixed backbone accumulate topological errors. (D) Empirical basis of the threshold: at n=1000 (seeds 401–409) DCM outperforms the simplified pipeline (nRF 0.083 ± 0.005 vs 0.102 on seed 401); raising SIMPLE_THRESHOLD to 2,000 forced the simplified pipeline at that scale, degraded accuracy, and was reverted.

**Figure 4.** Spaced k-mer parameter optimization. (A) nRF heatmap: k (4–8) × gap (0–3) at n=200; the red box marks the default k=5,gap2. (B) Adaptive parameter selection logic: n≤100→k=4,gap1; n>100→k=5,gap2; genome scale (>>1 kb)→multi-k ensemble (k=5 saturates its 1,024 possible 5-mers). (C) Adaptive k,gap selection as a function of dataset size n: k=4,gap1 up to n=100 and k=5,gap2 beyond; 5/5 repeats returned identical selections at every scale (n=20–500), with nRF deviation <0.005 across repeats. (D) Multi-k ensemble (k=5,7,9) vs single-k components and the full pipeline on the 130-seed INDELible benchmark (TRUE-relative): the ensemble significantly beats contiguous k=5 (p=1.6×10⁻⁷) and k=9 (p=7.3×10⁻⁴), ties the best single-k (k=7, p=0.23), and does not exceed the full spaced+NNI pipeline (p=2.2×10⁻⁴) — no general advantage over the default pipeline at L=500 bp.

**Figure 5.** Scalability and cross-platform deployment. (A) Wall-clock time vs n on a 4-core workstation (log–log; Table 11), with MAFFT+RAxML-NG runtimes for comparison; dotted line: O(n² log n) reference anchored at n=10,000 (69.8 s). (B) FastME BIONJ+BNNI vs NJ tree-building time: 3.5×, 15.1×, and 28.8× speedup at n=200, 500, and 1000. (C) Peak RAM vs n: below 1 GB even at n=10,000, where DCM caps memory by bounding subtree size (50 subtrees of 200 taxa). (D) FastME binary discovery cascade: bundled native binary → WSL install → system PATH, with an NJ fallback; the Windows-native and Linux binaries produce bit-identical trees.

**Figure 6.** Real-data validation. (A) Fusang 16S rRNA tree (73 taxa, k=5,gap2; 0.11 s without alignment); tip points coloured by phylum (corrected NCBI Taxonomy annotation; Supplementary Note S17). (B) Phylogenetic signal by taxonomic rank under the corrected annotation: same-rank tree-distance reduction for the Fusang tree of 31.4% at phylum (p<1×10⁻⁴), 56.5% at order (p<1×10⁻⁴) and 71.0% at family (p<1×10⁻⁴), rising monotonically towards finer ranks; the rebuilt FastTree2 tree attains 55.5%/80.4%/90.1% at the same ranks. (C) Architectural evolution from Fusang v1 (deep learning, GPU-dependent, ≤40 taxa) to the Tardigrade Edition (k-mer frequency vectors, CPU-only, 10,000+ taxa); inset: maximum dataset size on a log scale. (D) Gene-length DNA benchmark: 13 mitochondrial loci (168–1,866 bp) from 25 fish genomes; Fusang multi-k vs FastTree2 topological agreement (nRF) improves with locus length (mean 0.373 at ≥1 kb loci vs 0.591 at the three shortest loci; overall 0.451±0.127).

---

*Manuscript prepared for Nucleic Acids Research. Main text: approximately 12,000 words.*
