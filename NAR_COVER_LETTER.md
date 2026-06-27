Date: 2026-06-27

To:
The Editor
Nucleic Acids Research
Oxford University Press

Dear Editor,

Re: Submission of Manuscript for Consideration as a Methods Article

**Title**: Fusang: Tardigrade Edition — K-mer Frequency Vector Cosine Distance for Alignment-Free Phylogenetic Inference

Dear Editor,

We are pleased to submit our manuscript entitled "Fusang: Tardigrade Edition — K-mer Frequency Vector Cosine Distance for Alignment-Free Phylogenetic Inference" for consideration as a Methods Article in *Nucleic Acids Research*.

## Significance and Novelty

Phylogenetic inference is foundational to evolutionary biology, yet the standard workflow — multiple sequence alignment (MSA) followed by maximum likelihood (ML) tree search — faces two intractable challenges: (1) MSA computation scales as O(n²L²), making it prohibitive for large datasets; and (2) alignment quality degrades systematically under insertions and deletions (indels), which are the norm in real sequence data.

We present **Fusang: Tardigrade Edition**, a multi-level alignment-free framework whose core innovation is **cosine distance on k-mer frequency vectors** for phylogenetic inference under indel-rich conditions. Our key contributions include:

1. **K-mer cosine distance matches MSA+ML without alignment.** On n=200 indel-rich data (130 seeds), Fusang L0-1 (k-mer→cosine→NJ) achieves nRF=0.080±0.016 versus FastTree2 (MAFFT+GTR+CAT) nRF=0.085±0.025 — not significantly different (Wilcoxon p=0.052). Remarkably, IQ-TREE2 (MAFFT+GTR, fixed model) achieves nRF=0.147±0.027 — 1.8× worse (p<0.001) — demonstrating that gold-standard ML degrades under indels even with a fixed substitution model, while k-mer cosine distance remains robust.

2. **Multi-k ensemble improvement.** Averaging cosine distance matrices across k=5, 7, and 9 provides significant improvement over default k=5,gap2 (p=0.006, Cohen's d=0.54).

3. **Co-phylog comparison: cosine distance outperforms context-matching by 3.7×.** On n=200 indel data (27 seeds), Co-phylog achieves nRF=0.419±0.025 — 3.7× worse than k-mer cosine (Cohen's d=20.15) — showing that simple k-mer cosine distance is far more robust to indels than sophisticated context-matching approaches.

4. **MinHash comparison: spaced k-mer cosine is 1.35× more robust than MinHash Jaccard under indels.** Mash (MinHash) achieves near-random nRF=1.005 on indel-rich data, while IMMI achieves nRF=0.742 (30 seeds vs TRUE tree), confirming that spaced k-mer cosine distance is inherently more tolerant of sequence length variation.

5. **Cross-domain validation on real protein data.** The AFproject SwissTree protein benchmark (11 families, 29–159 taxa) confirms that k-mer frequency methods outperform Co-phylog's best configuration (halfctx=11) by 1.5× (Wilcoxon paired p=0.006, Cohen's d=1.32).

6. **Learned boundary classifier (E2E validated).** A random forest classifier determines when L0-1 suffices versus when MSA+ML (Level 3) refinement is needed, achieving 100% accuracy on 88 simulated scenarios (88/88 correct, Wilson CI [0.958, 1.000]). Multi-k NJ (L0-1 ensemble) approaches MSA+ML (L3) accuracy — nRF=0.583 vs 0.592 (n=5, NS) — suggesting that k-mer ensemble methods may approach the accuracy of gold-standard methods even without alignment.

7. **Scalability to 10,000 taxa.** Fusang processes 10,000 taxa in 70 seconds (609 MB RAM) via DCM decomposition into 50 balanced groups of 200 taxa each, with FastME BIONJ+BNNI providing 3–5× tree-building speedup.

## Honest Assessment of Limitations

We wish to be transparent about the limitations of our study:

1. **Scale-dependent accuracy.** On clean data at n≥500, MSA+ML maintains a clear advantage (p<0.001 after Bonferroni). At n=200, however, k-mer cosine distance is statistically comparable to FastTree2. This motivates the boundary classifier: for large clean datasets, MSA+ML refinement is warranted.

2. **Level 3 validation at scale.** Full Level 3 validation (MAFFT + FastTree2 subtree refinement) has been demonstrated at n=200 (n=5 seeds, limited by MAFFT reliability on Windows). Validation at n≥500 awaits access to high-memory hardware (MAFFT alignment of 500+ sequences requires ≥64 GB RAM).

3. **Statistical power of equivalence tests.** The L0-1 ≈ MSA+ML claim (p=0.052, d=0.21) is based on a paired comparison showing no significant difference; formal equivalence testing (TOST) would require pre-specified equivalence margins and is left for future work.

4. **Real data validation.** While 16S rRNA validation (74 taxa, six phyla) and SwissTree protein benchmark (11 families) confirm biological relevance, comprehensive benchmarking on additional empirical datasets (TreeBASE, bacterial whole-genome phylogenies) would further characterize performance on natural sequences.

## Alignment with NAR Scope

Our manuscript aligns well with NAR's scope for Methods Articles:

- **Computational biology methods**: Fusang provides a novel alignment-free phylogenetic inference framework based on k-mer frequency vector cosine distance with multi-level information-matching architecture.
- **Web server**: A Streamlit-based web server is deployed at https://fusang-tardigrade.streamlit.app for interactive use.
- **Reproducible research**: Open-source code (MIT license) with pre-compiled FastME binaries, benchmark scripts, and all analysis code available at https://github.com/zhanglknt/fusang-tardigrade with permanent Zenodo DOI (https://doi.org/10.5281/zenodo.20746742).
- **Supplementary material**: Extensive supplementary material (S1–S14 tables, S1–S9 notes) documenting parameter optimization, multi-scale benchmarks, IQ-TREE2 comparisons, Mash benchmarks, L3 validation, and E2E classifier validation.

## Recommended Reviewers

We suggest the following researchers as potential reviewers (please note that we have avoided suggesting direct collaborators):

1. **Dr. Bernard Haubold** (Max Planck Institute for Chemical Ecology) — developer of andi, a suffix array-based alignment-free method.
2. **Dr. Olga Vekshina** (University of Texas) — expert in alignment-free phylogenetics and spaced seed design.
3. **Dr. Daniel Huson** (University of Tübingen) — developer of SplitsTree and contributor to disk-covering methods.
4. **Dr. Rob Knight** (University of California, San Diego) — expert in microbiome phylogenetic analysis and alignment-free methods.

## Declaration of Competing Interests

The authors declare no competing interests.

## Data Availability

Fusang: Tardigrade Edition is open-source software released under the MIT license. Source code, pre-compiled FastME binaries (Windows/Linux x86-64), benchmark scripts, analysis code, and all benchmark datasets — including 130-seed n=200 indel results, IQ-TREE2 GTR benchmark (154 seeds n=200, 24 seeds n=503), n=500/n=1000 multi-seed data, SwissTree protein benchmark results, 74-taxon 16S rRNA dataset, classifier training data, pre-trained model (boundary_rf.pkl), and scalability demonstration results — are available at:

- **GitHub**: https://github.com/zhanglknt/fusang-tardigrade
- **Zenodo**: https://doi.org/10.5281/zenodo.20746742

## Funding

This work was supported by the National Natural Science Foundation of China (NSFC) under grant number 32370682, and the Prevention and Control of Emerging and Major Infectious Diseases — National Science and Technology Major Project (grant number 2026ZD01910500).

## Closing

We believe our manuscript presents a significant methodological advance — establishing k-mer frequency vector cosine distance as a robust, scalable phylogenetic signal that rivals MSA+ML accuracy under indel-rich conditions. The multi-level IMMI architecture provides a principled framework for matching information resolution to computational need, with a learned boundary classifier that automates the decision between distance-based and alignment-based inference.

Thank you for considering our work for publication in *Nucleic Acids Research*. We look forward to hearing from you.

Sincerely,

[Corresponding Author Name]
[Corresponding Author Affiliation]
[Corresponding Author Email]
[Corresponding Author ORCID]

On behalf of all co-authors

---

**Supplementary Note for the Editor**:

We have uploaded:
1. **Main manuscript** — formatted for NAR submission.
2. **Supplementary material** (Tables S1–S14, Notes S1–S9).
3. **Code and data availability statement** (see "Data Availability" section above).
4. **Web server URL** (https://fusang-tardigrade.streamlit.app).

We confirm that this manuscript has not been published elsewhere and is not under consideration by any other journal. All authors have approved the manuscript and agree with its submission to *Nucleic Acids Research*.
