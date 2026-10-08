#!/usr/bin/env python3
"""Fusang v2.5 master fix: P0 hard fixes + P1 narrative fixes.
Order: (1) prose replacements (old citation numbers OK, placeholders for new refs/tables)
       (2) table renumber pass  (3) citation renumber pass  (4) placeholder resolution
       (5) reference list replacement  (6) word count + self-checks
"""
import re, sys

PATH = r"D:\系统发育树项目\Fusang\Fusang-main\NAR_MANUSCRIPT_REVISED.md"
text = open(PATH, encoding='utf-8').read()
orig_len = len(text)

R = []  # (label, old, new)
def rep(label, old, new):
    R.append((label, old, new))

# ============ PROSE REPLACEMENTS ============

rep("running_title",
"**Running Title**: Spaced k-mer phylogenetics",
"**Running Title**: K-mer frequency vector phylogenetics")

# --- Abstract: single paragraph, NAR style ---
rep("abstract", """**Background**: Multiple sequence alignment (MSA) scales as O(n²L²) and introduces systematic errors under indels. Alignment-free methods are faster but have historically underperformed MSA-based maximum likelihood (ML) approaches.

**Results**: We present Fusang: Tardigrade Edition, an alignment-free framework that evaluates spaced k-mer frequency vector cosine distances for phylogenetic inference under indel-rich conditions. On simulated data (n=200, indel=0.02), Fusang achieves nRF=0.080±0.016 vs FastTree2 nRF=0.085±0.025 (112 seeds, p=0.052). Multi-k ensemble NJ produces trees with nRF=0.583±0.045 vs MAFFT+FastTree2 nRF=0.592±0.041 on n=5 valid seeds (p=0.24). Mash (MinHash) shows substantial degradation on indel-rich data (nRF=0.762±0.036, 30 seeds) while spaced k-mer cosine distances degrade modestly. A random forest boundary classifier achieves 100% accuracy (88/88 scenarios, 95% CI [0.958, 1.0]). On clean data, Fusang is competitive at n=200; MSA methods retain advantage at n≥500 (p<0.001). Fusang scales to 10,000 taxa in 54 seconds.

**Conclusion**: Spaced k-mer frequency vector cosine distances provide effective phylogenetic signal for indel-rich data without alignment. The multi-k ensemble achieves robust accuracy without manual k selection. MinHash-based approaches collapse under indels, underscoring the importance of k-mer representation and distance metric selection.""",
"""Multiple sequence alignment (MSA) scales poorly with dataset size and introduces systematic errors under insertions and deletions (indels), yet alignment-free methods have historically underperformed MSA-based maximum likelihood (ML) approaches. We present Fusang: Tardigrade Edition, an alignment-free framework that systematically evaluates k-mer frequency vector cosine distances — across spaced and contiguous patterns — for phylogenetic inference under indel-rich conditions. On simulated indel-rich data (n=200, indel rate=0.02), Fusang matches FastTree2 accuracy against the simulated ground-truth tree (nRF=0.080±0.016 vs 0.084±0.019, 112 seeds after pre-specified outlier exclusion, paired Wilcoxon p=0.052, borderline), while MSA methods retain a significant advantage on clean data at n≥500 (p<0.001). A multi-k distance ensemble further improves accuracy without manual k selection (nRF=0.105±0.021 vs 0.112±0.019, p=0.006). Under indels, both MinHash (Mash) and k-mer cosine distances degrade severely relative to their clean baselines (1.93× vs 1.97×), though cosine distances retain marginally lower absolute error (nRF=0.742 vs 0.762, 30 seeds). Cross-domain validation on 11 SwissTree protein families confirms k-mer frequency methods outperform context-matching by 1.5× (p=0.006). A random forest boundary classifier distinguishes homogeneous from structured datasets (88/88 scenarios, 95% CI [0.958, 1.0]). Fusang infers trees for 10,000 taxa in under a minute on a 4-core workstation.""")

# --- L38: add Wong citation ---
rep("intro_wong",
"each gap placement introduces potential error that propagates through phylogenetic inference [21,4].",
"each gap placement introduces potential error that propagates through phylogenetic inference [21,4,22].")

# --- L46: intro k-mer paragraph citations ---
rep("intro_kmer",
"In phylogenetic inference, prior work includes kmacs [15] (Leimeister & Morgenstern 2014, BMC Bioinformatics), which introduced gapped k-mer matching for sequence comparison; SpaMz (2016), which applied spaced word frequencies; and the Alfpy toolkit (Zielezinski et al. 2019, Genome Biology), which provides standardized implementations of multiple alignment-free methods including gapped k-mer variants.",
"In phylogenetic inference, prior alignment-free work (reviewed in [20,23]) includes hierarchical phylogenomic inference [1], kmacs [15] (gapped k-mismatch substring matching), spaced-word frequency methods {{SPACED}}, and the Alfpy toolkit with the AFproject benchmark [27], which provides standardized implementations of multiple alignment-free methods including gapped k-mer variants.")

rep("intro_finding3",
"and (3) MinHash-based methods collapse under indels, underscoring the importance of k-mer representation and distance metric selection.",
"and (3) MinHash-based methods degrade severely under indels (1.93× relative nRF increase), while k-mer cosine distances retain marginally lower absolute error.")

# --- L77: Methods JSD/table mention ---
rep("methods_jsd",
"Both metrics substantially outperform MinHash Jaccard under indels (Table 2).",
"Both metrics achieve lower absolute nRF than MinHash Jaccard under indels (see Results).")

# --- L97: Methods table mentions ---
rep("methods_tables",
"The SwissTree benchmark (Table 9) and alignment-free competitor comparison (Table 8) use NJ (BioPython) for consistency with the AFproject community standard, which requires identical tree builders across methods. All scalability benchmarks (Table 6) use FastME for its O(n²) speed advantage.",
"The SwissTree benchmark and alignment-free competitor comparison (see Results) use NJ (BioPython) for consistency with the AFproject community standard, which requires identical tree builders across methods. All scalability benchmarks use FastME for its O(n²) speed advantage.")

# --- L103: DCM citation fix ---
rep("dcm_cite", "the Disk-Covering Method (DCM [10,20])", "the Disk-Covering Method (DCM [10,21])")

# --- L113: bioNJ citation fix ---
rep("bionj_cite",
"FastME implements the Minimum Evolution principle with bioNJ [11] as the initial topology",
"FastME implements the Minimum Evolution principle with bioNJ {{GASCUEL}} as the initial topology")

# --- L151: power claim fix ---
rep("power",
"while 30-seed benchmarks provide lower power (~55% for d=0.5) and should be interpreted cautiously for non-significant results.",
"while 30-seed benchmarks provide moderate power (~75% for d=0.5) and should be interpreted cautiously for non-significant results.")

# --- L166: MAFFT citation ---
rep("mafft_cite",
"- **FastTree2 v2.2.0** [19]: GTR+CAT approximation, MAFFT alignment",
"- **FastTree2 v2.2.0** [19]: GTR+CAT approximation, MAFFT alignment [11]")

# --- L169: MashTree -> Mash ---
rep("mashtree",
"- **MashTree** [6]: contiguous k-mers (k=21), MinHash Jaccard, NJ",
"- **Mash v2.3** [6]: contiguous k-mers (k=21), MinHash Jaccard, NJ")

# --- L209: n=200 indel 30-seed value ---
rep("l209_nrf",
"(nRF: Fusang 0.078 ± 0.018 vs FastTree2 0.080 ± 0.017, 30 seeds; 112-seed post-exclusion benchmark: p=0.052, borderline)",
"(nRF: Fusang 0.077 ± 0.018 vs FastTree2 0.080 ± 0.017, 30 seeds; 112-seed post-exclusion benchmark: p=0.052, borderline)")

# --- Table 1 block: corrected SDs + single-frame notes ---
rep("table1", """**Table 1. Accuracy on clean data (no indels, L=500 bp, μ=0.05, multi-seed stats).**

**† NOTE: The two nRF columns use different reference trees and are not directly comparable. Fusang column = FT2-relative; FastTree2 column = TRUE-relative. See footnote for details.**

| n | Data Type | Fusang nRF ↓ † (30 seeds) | FastTree2 nRF ↓ ‡ (30 seeds) | Winner |
|---|-----------|---------------------------|----------------------------|--------|
| 200 | Clean | 0.102 ± 0.015 (k=5,gap2) | 0.096 ± 0.015 | FT2 (n.s.) |
| 200 | Indel (0.02) | **0.078 ± 0.018** (k=5,gap2) | 0.080 ± 0.017 | Fusang (n.s.) |
| 500 | Clean | 0.119 ± 0.020 (k=5,gap2) | **0.093 ± 0.015** | FT2 |
| 500 | Indel (0.02) | 0.095 ± 0.018 (k=5,gap2) | **0.083 ± 0.014** | FT2 |
| 1000 | Clean | 0.115 ± 0.022 (k=5,gap2) | **0.091 ± 0.016** | FT2 |
| 1000 | Indel (0.02) | **0.037 ± 0.006** (k=5,gap2) | — †† | — |

nRF=0: perfect match. Best result in **bold**. Values are mean ± standard deviation (30 seeds per condition, fixed seed set 70–99).  
**† Fusang column**: FT2-relative (nRF computed against FastTree2 reference tree).  
**‡ FastTree2 column**: TRUE-relative (nRF computed against simulated ground truth).  
**IMPORTANT**: These two columns use different reference trees and are NOT directly comparable. The n=1000 indel row reports Fusang vs FastTree2 only (TRUE tree unavailable at this scale), and is therefore not directly comparable to TRUE-relative values in Tables 2 and 10. The Abstract reports the 112-seed post-exclusion benchmark value (nRF=0.080, seeds 100–229) for the n=200 indel condition, which broadly agrees with the 30-seed value (0.078). The 112-seed benchmark (n=200, indel rate=0.02) achieved p=0.052 (Wilcoxon signed-rank test, borderline). After Bonferroni correction across 5 ground truth datasets (α=0.01), 3/5 remain significant — all in favor of FastTree2 at n≥500 (Supplementary Table S8).""",
"""**Table 1. Accuracy on clean and indel-rich data (L=500 bp, μ=0.05, 30 seeds per condition, seed set 70–99). All nRF values are TRUE-relative (vs simulated ground-truth tree) except where noted.**

| n | Data Type | Fusang nRF ↓ (30 seeds) | FastTree2 nRF ↓ (30 seeds) | Winner |
|---|-----------|---------------------------|----------------------------|--------|
| 200 | Clean | 0.102 ± 0.019 (k=5,gap2) | 0.096 ± 0.019 | FT2 (n.s.) |
| 200 | Indel (0.02) | **0.077 ± 0.018** (k=5,gap2) | 0.080 ± 0.017 | Fusang (n.s.) |
| 500 | Clean | 0.119 ± 0.011 (k=5,gap2) | **0.093 ± 0.013** | FT2 |
| 500 | Indel (0.02) | 0.095 ± 0.015 (k=5,gap2) | **0.083 ± 0.014** | FT2 |
| 1000 | Clean | 0.115 ± 0.011 (k=5,gap2) | **0.091 ± 0.010** | FT2 |
| 1000 | Indel (0.02) | **0.037 ± 0.006** (k=5,gap2) † | — | — |

nRF=0: perfect match. Best result in **bold**. Values are mean ± standard deviation (30 seeds per condition, fixed seed set 70–99). For both methods, nRF is computed against the simulated ground-truth (TRUE) tree.
† The n=1000 indel row is FT2-relative (Fusang tree vs FastTree2 tree; TRUE-tree comparison not available at this scale) and is therefore not comparable with TRUE-relative values elsewhere. The Abstract reports the independent 112-seed post-exclusion benchmark value (nRF=0.080, seeds 100–229) for the n=200 indel condition, which broadly agrees with the 30-seed value (0.077). That 112-seed benchmark (n=200, indel rate=0.02) achieved p=0.052 (paired Wilcoxon signed-rank test on TRUE-relative nRF, borderline). After Bonferroni correction across 5 ground truth datasets (α=0.01), 3/5 remain significant — all in favor of FastTree2 at n≥500 (Supplementary Table S8).""")

# --- L233: Mash intro rewrite (n=1 fix + protocol distinction) ---
rep("mash_intro", """To evaluate the robustness of spaced k-mer cosine distances against a widely-used alignment-free alternative, we compared Fusang (spaced k=5,gap2, cosine+ NJ) with Mash (contiguous k=21, MinHash Jaccard + NJ) on both clean and indel-rich data (n=200, sub=0.05, indel=0.02), benchmarking against the TRUE simulated coalescent tree. Fusang used 30 seeds; Mash was run on a single representative seed (Windows binary unavailable, Linux-generated trees). **This Mash comparison is preliminary (n=1) and requires multi-seed validation for definitive quantification; the single-seed result is reported to illustrate the qualitative trend.**""",
"""To evaluate the robustness of k-mer cosine distances against a widely-used alignment-free alternative, we compared Fusang (spaced k=5,gap2, cosine + NJ) with Mash (contiguous k=21, MinHash Jaccard + NJ) on both clean and indel-rich data (n=200, sub=0.05, indel=0.02), each with 30 seeds, benchmarking against the TRUE simulated coalescent tree (per-seed data: Supplementary Table S14). Note that this benchmark uses an independent coalescent simulation protocol, distinct from the guide-tree protocol used in Table 1 and Table 3; absolute nRF values are higher for all methods under this protocol and should not be compared across benchmarks.""")

# --- Table 2 note: honest framing ---
rep("table2_note", """**Note: Mash results are from 30-seed benchmark (this work, WSL2 Ubuntu 24.04, Mash v2.3). Multi-seed statistics with standard deviation are reported. Mash on clean data (nRF=0.394) substantially underperforms Fusang (nRF=0.112, Table 1); the previously reported single-seed value (nRF=0.162) [20] appears to have been based on an incorrect tree file (mashtree.pl log misinterpreted as Newick).**""",
"""**Note: Mash results are from a 30-seed benchmark (this work, WSL2 Ubuntu 24.04, Mash v2.3); Fusang results are from the same 30-seed coalescent protocol. An earlier single-seed Mash value (nRF=0.162 on clean data) reported during development appears to have been based on an incorrectly parsed tree file and is superseded by this benchmark.**""")

# --- L250: Mash finding 1 ---
rep("mash_f1", """1. **Mash underperforms k-mer cosine distances on clean data.** On clean substitution-only sequences, Mash (nRF=0.394 ± 0.045) is slightly less accurate than Fusang's k-mer cosine (nRF=0.376 ± 0.045, Table 2). The previously reported single-seed value (nRF=0.162 [20]) appears to have been based on an incorrectly parsed tree file (mashtree.pl log misinterpreted as Newick). The multi-seed benchmark (n=30) provides definitive quantification: Mash's MinHash Jaccard is not superior to k-mer cosine on clean gene-length sequences.""",
"""1. **Mash slightly underperforms k-mer cosine distances on clean data.** On clean substitution-only sequences, Mash (nRF=0.394 ± 0.045) is slightly less accurate than Fusang's k-mer cosine (nRF=0.376 ± 0.045, Table 2). The multi-seed benchmark (n=30) provides robust quantification: MinHash Jaccard is not superior to k-mer cosine on clean gene-length sequences.""")

# --- L252: Mash finding 2 (honest degradation) ---
rep("mash_f2", """2. **Mash shows substantial (though not random) degradation on indel-rich data.** Mash's indel nRF=0.762 ± 0.036 corresponds to ~76% bipartition mismatch — not random (nRF=1.0), but catastrophically inaccurate. The indel degradation factor is 1.93× (nRF_indel / nRF_clean). In comparison, Fusang's k-mer cosine degrades by 1.97× (nRF=0.376 → 0.742, Table 2). Fusang's absolute indel nRF (0.742) is marginally better than Mash's (0.762), and Fusang starts from a comparable clean baseline. The mechanistic explanation is consistent: MinHash sketches discard positional information entirely, making them sensitive to sequence length variation under indels. Spaced k-mer frequency vectors preserve relative positional information through the gap pattern, providing superior indel robustness.""",
"""2. **Both methods degrade severely on indel-rich data, with k-mer cosine retaining marginally lower absolute error.** Mash's indel nRF=0.762 ± 0.036 corresponds to ~76% bipartition mismatch. The indel degradation factor is 1.93× (nRF_indel / nRF_clean); Fusang's k-mer cosine degrades by a comparable 1.97× (nRF=0.376 → 0.742, Table 2). Fusang's absolute indel nRF (0.742) is marginally lower than Mash's (0.762); because a formal between-method test on matched seeds was not performed, this small absolute difference should be interpreted descriptively. A plausible mechanistic explanation for the residual difference is that MinHash sketches discard positional information entirely, making them sensitive to indel-induced length variation, whereas k-mer frequency vectors retain distributional information across the sequence.""")

# --- L258: frame annotation for 112-seed ---
rep("l258_frame",
"At the biologically realistic indel rate of 0.02, Fusang approaches FastTree2 accuracy (nRF: Fusang 0.080 ± 0.016 vs FastTree2 0.084 ± 0.019, 112 seeds, Wilcoxon p=0.052, borderline).",
"At the biologically realistic indel rate of 0.02, Fusang approaches FastTree2 accuracy (TRUE-relative nRF: Fusang 0.080 ± 0.016 vs FastTree2 0.084 ± 0.019, 112 seeds, paired Wilcoxon p=0.052, borderline).")

# --- Table 3 caption + note: single frame ---
rep("table3", """**Table 3. Indel rate scan: Fusang vs FastTree2 (n=200, L=500 bp, 112 seeds after outlier exclusion).**""",
"""**Table 3. Indel rate scan: Fusang vs FastTree2 (n=200, L=500 bp, 112 seeds after outlier exclusion). All nRF values are TRUE-relative (vs simulated ground-truth tree).**""")

rep("table3_note", """**Reference frames**: Fusang column = FT2-relative; FastTree2 column = TRUE-relative. These two columns use different reference trees and are NOT directly comparable. The Fusang advantage over FastTree2 grows monotonically with indel rate: from tie at 0.005 to a 13.3% relative advantage at 0.05 (calculated as (FT2_nRF − Fusang_nRF) / FT2_nRF × 100%, provided in Supplementary Table S2 for reference). See Table 10 for TRUE-relative values (single-k NJ: nRF=0.743, multi-k NJ: nRF=0.583).""",
"""Both columns are TRUE-relative and directly comparable per seed. The Fusang advantage over FastTree2 grows monotonically with indel rate: from tie at 0.005 to a 13.3% relative advantage at 0.05 (calculated as (FT2_nRF − Fusang_nRF) / FT2_nRF × 100%). Table {{T5}} reports an independent pipeline-level validation under a coalescent simulation protocol (single-k NJ: nRF=0.743, multi-k NJ: nRF=0.583), where absolute nRF values are higher for all methods.""")

# --- L275: add Supp Figure S4 ---
rep("l275_supp",
"with paired Wilcoxon signed-rank test (Figure 2, Supplementary Table S4).",
"with paired Wilcoxon signed-rank test (Figure 2, Supplementary Figure S4, Supplementary Table S4).")

# --- L277: FT2 SD fix + frame ---
rep("l277_stats",
"- **Overall (112 seeds)**: Fusang nRF=0.080 ± 0.016 vs FastTree2 nRF=0.085 ± 0.025; Cohen's d=−0.20 [95% CI: −0.42, 0.02]; Fusang lower nRF in 60/112 seeds (53.6%)",
"- **Overall (112 seeds, both methods evaluated against the simulated ground-truth tree)**: Fusang nRF=0.080 ± 0.016 vs FastTree2 nRF=0.084 ± 0.019; Cohen's d=−0.20 [95% CI: −0.42, 0.02]; Fusang lower nRF in 60/112 seeds (53.6%)")

# --- L304: remove duplicate t-test p ---
rep("table7_ttest",
"""- Wilcoxon signed-rank test (vs default spaced): p = **0.006**
- Paired t-test (vs default spaced): p = 0.007
- Cohen's d (vs default spaced) = 0.54""",
"""- Wilcoxon signed-rank test (vs default spaced, pre-specified primary test): p = **0.006** (paired t-test as sensitivity analysis: p = 0.007)
- Cohen's d (vs default spaced) = 0.54""")

# --- L307: add multi-k supp refs ---
rep("multik_supp",
"This result demonstrates that distance matrix fusion across multiple contiguous k-mer resolutions achieves accuracy comparable to the best single-k configuration, with a modest improvement over the default spaced k-mer.",
"This result demonstrates that distance matrix fusion across multiple contiguous k-mer resolutions achieves accuracy comparable to the best single-k configuration, with a modest improvement over the default spaced k-mer (per-seed data: Supplementary Figure S7, Supplementary Table S9).")

# --- Table 10 caption: protocol note + supp ref ---
rep("table10_caption",
"**Table 10. Pipeline-level validation against TRUE tree (n=200, indel=0.02, 30 seeds). All nRF values are TRUE-relative (vs simulated coalescent ground truth).**",
"**Table {{T5}}. Pipeline-level validation against TRUE tree (n=200, indel=0.02, 30 seeds). All nRF values are TRUE-relative (vs simulated coalescent ground truth). Note: this benchmark uses a coalescent simulation protocol, independent of the guide-tree protocol used in Tables 1 and 3; absolute nRF values are higher for all methods under this protocol and are not directly comparable across benchmarks. Per-seed data: Supplementary Table S13.**")

# --- L332: combined table refs ---
rep("l332_tables", "(Tables 1, 3, 7)", "(Tables 1, 3, 4)")

# --- L336: classifier supp ref ---
rep("classifier_supp",
"on an independent test set of 88 scenarios spanning two fundamentally different tree types:",
"on an independent test set of 88 scenarios spanning two fundamentally different tree types (per-scenario results: Supplementary Table S15):")

# --- L357: AF competitor supp ref ---
rep("af_supp",
"We compared Fusang against two established alignment-free phylogenetic methods on simulated indel-rich data (n=200, sub=0.05, indel=0.02, 27 seeds with valid reference trees):",
"We compared Fusang against two established alignment-free phylogenetic methods on simulated indel-rich data (n=200, sub=0.05, indel=0.02, 27 seeds with valid reference trees; per-seed data: Supplementary Table S10):")

# --- Table 8: multi-k type fix ---
rep("table8_multik",
"| **Fusang** multi-k ensemble | Multi-k spaced | **0.105** | 0.021 | — |",
"| **Fusang** multi-k ensemble | Multi-k contiguous (k=5,7,9 avg) | **0.105** | 0.021 | — |")

# --- L386: Figure 4 + Supp Fig S2 citations ---
rep("l386_figs",
"revealed a robust relationship between optimal gap and dataset size (Table 4, Supplementary Figure S1).",
"revealed a robust relationship between optimal gap and dataset size (Table 4, Figure 4, Supplementary Figure S1, Supplementary Figure S2).")

# --- Table 4: verified SDs ---
rep("table4", """| 50 | 4,gap1 | 0.102 ± 0.020 | 0.095 ± 0.018 | Comparable |
| 100 | 5,gap2 | 0.115 ± 0.022 | 0.098 ± 0.016 | FT2 better |
| 200 | 5,gap2 | 0.102 ± 0.015 | 0.096 ± 0.015 | Comparable (clean) |
| 200 | 5,gap2 | **0.078 ± 0.018** | 0.080 ± 0.017 | Fusang better (indel, 112 seeds: p=0.052, borderline) |
| 500 | 5,gap2 | 0.119 ± 0.020 | **0.093 ± 0.015** | FT2 better (clean) |
| 1000 | 5,gap2 | 0.115 ± 0.022 | **0.091 ± 0.016** | FT2 better (clean) |""",
"""| 50 | 4,gap1 | 0.102 ± 0.020 | 0.095 ± 0.018 | Comparable |
| 100 | 5,gap2 | 0.115 ± 0.022 | 0.098 ± 0.016 | FT2 better |
| 200 | 5,gap2 | 0.102 ± 0.019 | 0.096 ± 0.019 | Comparable (clean) |
| 200 | 5,gap2 | **0.077 ± 0.018** | 0.080 ± 0.017 | Fusang better (indel, 112 seeds: p=0.052, borderline) |
| 500 | 5,gap2 | 0.119 ± 0.011 | **0.093 ± 0.013** | FT2 better (clean) |
| 1000 | 5,gap2 | 0.115 ± 0.011 | **0.091 ± 0.010** | FT2 better (clean) |""")

# --- L192: Figure 3 citation ---
rep("fig3_cite",
"achieved nRF=**0.005** — a ~78× reduction in topological error on this illustrative single seed (Supplementary Table S3).",
"achieved nRF=**0.005** — a ~78× reduction in topological error on this illustrative single seed (Figure 3, Supplementary Figure S3, Supplementary Table S3).")

# --- L407: 16S boundary framing + Figure 6 ---
rep("l407_16s",
"While Fusang's recovery rate compares favorably to random expectation, it falls short of MSA-based accuracy on this dataset. This is expected: 16S rRNA genes contain highly conserved regions where positional homology from alignment provides strong phylogenetic signal that k-mer frequency vectors, which discard positional information, cannot fully capture.",
"While Fusang's recovery rate compares favorably to random expectation, it falls short of MSA-based accuracy on this dataset (Figure 6, Supplementary Figure S5). This result delineates an applicability boundary of the method rather than contradicting the indel-robustness findings on simulated data: 16S rRNA genes contain highly conserved regions where positional homology from alignment provides strong phylogenetic signal that k-mer frequency vectors, which discard positional information, cannot fully capture.")

# --- L421: 16S transfer claim soften ---
rep("l421_soften",
"These results demonstrate that k-mer frequency features optimized on simulated data transfer to real sequences without parameter tuning, confirming that alignment-free approaches recover genuine phylogenetic signal rather than simulation artifacts.",
"These results show that k-mer frequency features recover genuine phylogenetic signal from real sequences without parameter tuning, while delineating the sequence types — highly structured genes with strong positional conservation — where alignment-based methods retain a clear advantage.")

# --- Table 9 caption: supp refs ---
rep("table9_caption",
"**Table 9. SwissTree gene tree benchmark (11 families, protein sequences, AFproject standard).**",
"**Table {{T10}}. SwissTree gene tree benchmark (11 families, protein sequences, AFproject standard). Per-family data: Supplementary Tables S11–S12.**")

# --- L446: contribution framing ---
rep("l446_frame",
"The combination of spaced k-mer frequency vectors with cosine distance — systematically evaluated here for phylogenetic inference under indel-rich conditions — is the core contribution of this work.",
"The systematic evaluation of k-mer frequency vector cosine distances — across spaced and contiguous patterns and across sequence alphabets — is the core contribution of this work.")

# --- L450: Figure 5 citation ---
rep("fig5_cite",
"Fusang completes phylogenetic inference on 10,000 taxa in 54.4 seconds (Table 6),",
"Fusang completes phylogenetic inference on 10,000 taxa in 54.4 seconds (Table 6, Figure 5),")

# --- Discussion point 1 (spaced repositioning) ---
rep("disc_p1", """1. **Both spaced k-mer representation and cosine distance contribute to indel robustness.** On clean data at the tested indel rate (0.02), contiguous k=5 cosine (nRF=0.099) slightly outperforms spaced k=5,gap2 (nRF=0.112); the multi-k ensemble matches contiguous k=5 (both nRF=0.105). While cosine distance is the primary accuracy driver on clean data, spaced k-mers provide theoretical robustness against larger or more frequent indels by skipping over length variation. The combination of spaced k-mer frequency vectors with cosine distance — systematically evaluated here for phylogenetic inference under indel-rich conditions — is the core methodological contribution.""",
"""1. **Cosine distance is the primary accuracy driver; spaced patterns provide domain-dependent secondary benefits.** At the tested indel rate (0.02), contiguous k=5 cosine (nRF=0.099) slightly outperforms spaced k=5,gap2 (nRF=0.112), and the multi-k ensemble matches contiguous k=5 (both nRF=0.105). Spaced patterns show their clearest benefits in specific regimes: the best protein-family configuration is spaced (k=4,gap1, nRF=0.239, Table {{T10}}), and on indel data a wider gap (gap3) improves over gap2 by 10.5% (Supplementary Table S2). The systematic evaluation of k-mer frequency vector cosine distances across patterns and sequence types is the core methodological contribution.""")

# --- Discussion point 2 (Mash honest) ---
rep("disc_p2", """2. **Spaced k-mers provide robustness against MinHash collapse under indels.** On indel-rich data, Mash (MinHash Jaccard, n=30) produces trees with nRF=0.762 ± 0.036, while spaced k-mer cosine degrades to nRF=0.742 ± 0.045 (both vs TRUE tree, 30 seeds). The previously reported single-seed value (nRF=1.005 [20]) appears to have been based on an incorrectly parsed tree file. Spaced k-mer frequency vectors tolerate indels by preserving relative positional information through the gap pattern, whereas MinHash sketches discard positions entirely. This result — now validated on 30 seeds — confirms that spaced k-mer frequency vectors provide indel robustness that MinHash cannot match.""",
"""2. **MinHash and k-mer cosine both degrade severely under indels; absolute error remains marginally lower for cosine distances.** On indel-rich data, Mash (MinHash Jaccard, 30 seeds) produces trees with nRF=0.762 ± 0.036, while k-mer cosine degrades to nRF=0.742 ± 0.045 (both vs TRUE tree). Relative degradation is comparable (1.93× vs 1.97× from clean baselines). The mechanistic asymmetry — MinHash sketches discard positional information entirely while frequency vectors retain distributional signal — is a plausible explanation for the small absolute gap, but the practical difference at this indel rate is minor.""")

# --- Discussion point 3 ---
rep("disc_p3", """3. **The distance metric and k-mer representation must fit the data.** On clean data, MinHash Jaccard (nRF=0.394 ± 0.045, 30-seed benchmark) slightly underperforms k-mer cosine (nRF=0.376 ± 0.045). On indel data, MinHash degrades substantially (nRF=0.762 ± 0.036) while k-mer cosine degrades more modestly (nRF=0.742 ± 0.045). This underscores that "alignment-free" is not a monolithic category: method and parameter selection should be guided by the expected evolutionary process and data characteristics.""",
"""3. **The distance metric and k-mer representation must fit the data.** On clean data, MinHash Jaccard (nRF=0.394 ± 0.045, 30-seed benchmark) slightly underperforms k-mer cosine (nRF=0.376 ± 0.045); on indel data, both degrade severely (0.762 vs 0.742). This underscores that "alignment-free" is not a monolithic category: method and parameter selection should be guided by the expected evolutionary process and data characteristics.""")

# --- L486: v1 citation ---
rep("v1_cite",
"Our v1 (NAR 2023) established proof-of-concept",
"Our v1 [24] established proof-of-concept")

# --- L519: remove [20] self-cite ---
rep("l519_mash",
"The previously reported single-seed value (nRF=1.005) [20] appears to have been based on an incorrectly parsed tree file.",
"An earlier single-seed value (nRF=1.005) reported during development appears to have been based on an incorrectly parsed tree file and is superseded by this benchmark.")

# --- L525: n=1000 indel framing ---
rep("l525_n1000",
"The n=1000 indel benchmark (30 seeds, sub=0.05, indel=0.02) shows Fusang nRF=0.037 ± 0.006 vs FastTree2 reference, indicating high accuracy at scale.",
"The n=1000 indel benchmark (30 seeds, sub=0.05, indel=0.02) shows Fusang nRF=0.037 ± 0.006 vs the FastTree2 reference, indicating high topological similarity to the alignment-based reconstruction at scale (FT2-relative; TRUE-relative accuracy is unavailable at this scale).")

# --- L527: AF methods future work with real citations ---
rep("l527_af",
"We note several alignment-free methods not included in our benchmarks: Skmer (Sarmashghi et al. 2019, Genome Biology), which uses k-mer statistics for genome-scale distance estimation and would be a natural comparator for Mash; AAF/FFP (Sims et al. 2009, PNAS) and kr (Leimeister et al. 2014), which pioneered feature frequency profiles and gapped k-mer matching for sequence comparison; and the Alfpy toolkit (Zielezinski et al. 2019), which provides standardized implementations of multiple alignment-free methods including gapped k-mer variants. Systematic comparison with these approaches, particularly Skmer for indel robustness assessment, represents important future work.",
"We note several alignment-free methods not included in our benchmarks: Skmer {{SKMER}}, which uses k-mer statistics for genome-skim distance estimation; FFP {{FFP}}, which pioneered feature frequency profiles for whole-genome comparison; and the broader spaced-word method family {{SPACED}} and kmacs-class gapped-matching approaches [15]. Systematic comparison with these approaches, particularly Skmer and kmacs [15] under indel-rich conditions, represents important future work.")

# --- L535: remove [20] self-cite ---
rep("l535_mash",
"The previously reported single-seed value (nRF=0.162 [20]) appears to have been based on an incorrectly parsed tree file.",
"An earlier single-seed value (nRF=0.162) reported during development appears to have been based on an incorrectly parsed tree file.")

# --- L537: MinHash collapse wording ---
rep("l537_collapse",
"Under indels, MinHash-based approaches collapse to random and should be avoided.",
"Under indels, MinHash-based approaches degrade severely (nRF≈0.76, ~76% bipartition error on gene-length sequences) and should be avoided for this data type.")

# --- Zenodo DOI ---
rep("zenodo",
"- **Zenodo**: DOI to be assigned upon acceptance (archived source code, benchmark datasets, and supplementary materials)",
"- **Zenodo**: https://doi.org/10.5281/zenodo.20746742 (archived source code, benchmark datasets, and supplementary materials)")

# ============ APPLY PROSE REPLACEMENTS ============
missed = []
for label, old, new in R:
    if old in text:
        text = text.replace(old, new, 1)
    else:
        missed.append(label)
print(f"prose replacements: {len(R)-len(missed)}/{len(R)} applied")
if missed:
    print("MISSED:", missed)
    sys.exit(1)

# ============ TABLE RENUMBER (old -> new) ============
# order of first appearance in Results: 1,2,3,7,10,11,8,4,5,9,6
tbl_map = {1:1, 2:2, 3:3, 7:4, 10:5, 11:6, 8:7, 4:8, 5:9, 9:10, 6:11}
# via unique placeholders, longest first
for old_n in sorted(tbl_map, key=lambda x: -x):
    text = re.sub(rf'\bTable {old_n}(?!\d)', f'TBLX{tbl_map[old_n]}X', text)
text = text.replace('TBLX', 'Table ').replace('X\n', '\n')
# fix leftover X marker: 'TBLX5X' -> 'Table 5X'? handle properly:
text = re.sub(r'Table (\d+)X\b', r'Table \1', text)

# ============ CITATION RENUMBER ============
cite_map = {8:1, 9:2, 21:3, 4:4, 22:5, 24:6, 16:7, 7:8, 20:9, 23:10,
            1:11, 15:12, 27:14, 10:15, 2:16, 13:17, 5:19, 19:20, 11:21,
            12:22, 18:23, 6:24, 25:25, 26:26, 14:27, 3:28, 17:29}

# protect reference section (will be replaced wholesale later)
ref_start = text.index('## REFERENCES')
ref_end = text.index('## FIGURE LEGENDS')
body, refsec, tail = text[:ref_start], text[ref_start:ref_end], text[ref_end:]

def renum_cites(m):
    raw = m.group(1)
    if raw.replace(' ', '') == '0,1':
        return m.group(0)  # math interval [0,1], not a citation
    nums = [int(x.strip()) for x in raw.split(',')]
    try:
        new = sorted(cite_map[n] for n in nums)
    except KeyError as e:
        raise SystemExit(f"unmapped citation [{raw}]: missing {e}")
    return '[' + ','.join(str(n) for n in new) + ']'

# [n] or [n,m,...] with optional spaces; decimals excluded by \d+ boundary
body = re.sub(r'\[(\d+(?:\s*,\s*\d+)*)\]', renum_cites, body)

# ============ PLACEHOLDER RESOLUTION ============
placeholders = {
    '{{SPACED}}': '[13]', '{{GASCUEL}}': '[18]', '{{SKMER}}': '[30]', '{{FFP}}': '[31]',
    '{{T5}}': 'Table 5', '{{T10}}': 'Table 10',
}
for ph, val in placeholders.items():
    body = body.replace(ph, val)

# ============ REFERENCE LIST (new citation order) ============
new_refs = """## REFERENCES

1. Hadfield, J., Megill, C., Bell, S.M., Huddleston, J., Potter, B. et al. (2018) Nextstrain: real-time tracking of pathogen evolution. *Bioinformatics*, 34, 4121–4123. DOI: 10.1093/bioinformatics/bty407

2. Hug, L.A., Baker, B.J., Anantharaman, K., Brown, C.T., Probst, A.J. et al. (2016) A new view of the tree of life. *Nat. Microbiol.*, 1, 16048. DOI: 10.1038/nmicrobiol.2016.48

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

14. Zielezinski, A., Girgis, H.Z., Bernard, G., Leimeister, C.-A., Tang, K. et al. (2019) Benchmarking of alignment-free sequence comparison methods. *Genome Biol.*, 20, 144. DOI: 10.1186/s13059-019-1755-7

15. Huson, D.H. et al. (1999) Disk-covering, a fast-converging method for phylogenetic tree reconstruction. *J. Comput. Biol.*, 6, 369–386. DOI: 10.1089/106652799318337

16. Berger, S.A. et al. (2011) Performance, accuracy, and web server for evolutionary placement of short sequence reads under maximum likelihood. *Syst. Biol.*, 60, 291–302. DOI: 10.1093/sysbio/syr010

17. Lefort, V. et al. (2015) FastME 2.0: a comprehensive, accurate, and fast distance-based phylogeny inference program. *Mol. Biol. Evol.*, 32, 2798–2800. DOI: 10.1093/molbev/msv150

18. Gascuel, O. (1997) BIONJ: an improved version of the NJ algorithm based on a simple model of sequence data. *Mol. Biol. Evol.*, 14, 685–695. DOI: 10.1093/oxfordjournals.molbev.a025808

19. Fletcher, W. and Yang, Z. (2009) INDELible: a flexible simulator of biological sequence evolution. *Mol. Biol. Evol.*, 26, 1879–1888. DOI: 10.1093/molbev/msp098

20. Price, M.N. et al. (2010) FastTree 2 — approximately maximum-likelihood trees for large alignments. *PLoS ONE*, 5, e9490. DOI: 10.1371/journal.pone.0009490

21. Katoh, K. and Standley, D.M. (2013) MAFFT multiple sequence alignment software version 7. *Mol. Biol. Evol.*, 30, 772–780. DOI: 10.1093/molbev/mst010

22. Kozlov, A.M. et al. (2019) RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference. *Bioinformatics*, 35, 4453–4455. DOI: 10.1093/bioinformatics/btz305

23. Minh, B.Q. et al. (2020) IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. *Mol. Biol. Evol.*, 37, 1530–1534. DOI: 10.1093/molbev/msaa015

24. Ondov, B.D., Treangen, T.J., Melsted, P., Mallonee, A.B., Bergman, N.H. et al. (2016) Mash: fast genome and metagenome distance estimation using MinHash. *Genome Biol.*, 17, 132. DOI: 10.1186/s13059-016-0997-x

25. Haubold, B. et al. (2015) andi: Fast and accurate estimation of evolutionary distances between closely related genomes. *Bioinformatics*, 31, 1163–1167. DOI: 10.1093/bioinformatics/btv047

26. Yi, H. and Jin, G. (2013) Co-phylog: an assembly-free phylogenomic approach for closely related organisms. *Nucleic Acids Res.*, 41, e75. DOI: 10.1093/nar/gkt165

27. Lunter, G. et al. (2006) Bayesian coestimation of phylogeny and sequence alignment. *BMC Bioinformatics*, 7, 320. DOI: 10.1186/1471-2105-7-320

28. Cartwright, R.A. (2009) Problems and solutions for estimating indel rates and length distributions. *Mol. Biol. Evol.*, 26, 473–480. DOI: 10.1093/molbev/msn275

29. Ma, B. et al. (2002) PatternHunter: faster and more sensitive homology search. *Bioinformatics*, 18, 440–445. DOI: 10.1093/bioinformatics/18.3.440

30. Sarmashghi, S., Bohmann, K., Gilbert, M.T.P., Bafna, V. and Mirarab, S. (2019) Skmer: assembly-free and alignment-free sample identification using genome skims. *Genome Biol.*, 20, 34. DOI: 10.1186/s13059-019-1632-4

31. Sims, G.E., Jun, S.R., Wu, G.A. and Kim, S.H. (2009) Alignment-free genome comparison with feature frequency profiles (FFP) and optimal resolutions. *Proc. Natl. Acad. Sci. USA*, 106, 2677–2682. DOI: 10.1073/pnas.0813249106

---

"""

# ============ WORD COUNT (main text: ABSTRACT .. end of DISCUSSION) ============
main_start = body.index('## ABSTRACT')
main_end = body.index('## SUPPLEMENTARY MATERIAL')
main_text = body[main_start:main_end]
words = len(re.findall(r"[A-Za-z0-9][\w\-\+\.=×%<>/']*", main_text))
body = body.replace(
    '*Manuscript prepared for Nucleic Acids Research. Main text: approximately 6,000 words.*',
    f'*Manuscript prepared for Nucleic Acids Research. Main text: approximately {round(words, -2):,} words.*'
    if '*Manuscript prepared for Nucleic Acids Research. Main text: approximately 6,000 words.*' in body else body)
# the word-count line is actually at the very end (tail); handle there
if '*Manuscript prepared for Nucleic Acids Research. Main text: approximately 6,000 words.*' in tail:
    tail = tail.replace(
        '*Manuscript prepared for Nucleic Acids Research. Main text: approximately 6,000 words.*',
        f'*Manuscript prepared for Nucleic Acids Research. Main text: approximately {round(words, -2):,} words.*')

final = body + new_refs + tail
open(PATH, 'w', encoding='utf-8').write(final)

# ============ SELF-CHECKS ============
print(f"\n=== SELF-CHECKS ===")
print(f"main-text words (Abstract..Discussion): {words}")
# citation sanity
cites = set()
for m in re.finditer(r'\[(\d+(?:,\d+)*)\]', body):
    for x in m.group(1).split(','):
        cites.add(int(x))
print(f"citations used: {sorted(cites)}")
print(f"max citation: {max(cites)} (should be 31)")
uncited = [n for n in range(1, 32) if n not in cites]
print(f"uncited refs: {uncited} (should be [])")
# table sanity
tabs = sorted(set(int(m.group(1)) for m in re.finditer(r'\bTable (\d+)(?!\d)', body)))
print(f"tables referenced: {tabs}")
# abstract word count
abs_start = body.index('## ABSTRACT')
abs_end = body.index('## INTRODUCTION')
abs_words = len(re.findall(r"[A-Za-z0-9][\w\-\+\.=×%<>/']*", body[abs_start:abs_end]))
print(f"abstract words: {abs_words} (limit 250)")
# leftover placeholders
left = re.findall(r'\{\{[^}]+\}\}', body + tail)
print(f"leftover placeholders: {left} (should be [])")
# leftover old-style issues
print(f"'Table 10' count: {body.count('Table 10')}, 'Table 11' count: {body.count('Table 11')}")
print(f"file size: {orig_len} -> {len(final)}")
print("DONE")
