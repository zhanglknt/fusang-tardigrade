#!/usr/bin/env python3
"""fix_v26_master.py — v2.6 batch fixes after 4-reviewer verification round.

Fixes all CLOSED-path residual items flagged in review_v3_verify_{phylo,stats,format,general}.md:
- Table 3 column header (FT2-rel -> TRUE-rel)
- Table first-citation order made monotonic (no table renumbering needed; citation edits only)
- Table 4/6/7/10 missing in-prose citations
- 16S gap1/gap2 contradiction (unified to k=5,gap2 per IMMI L0-1 default provenance)
- sigma=0.017 -> 0.016; 0.008/6.7% -> 0.007/6.3%; 13.3% -> 13.2%; 4.6-4.7% -> 4.5-4.8%
- Methods L62 spaced k-mer notation rewrite; Table 10 pattern labels made consistent
- L42 spaced selling-point softening; L48 contribution 1 rebase to L1 vs L0
- Abstract "simulated" qualifier; Winner -> Tie (n.s.)
- JSD one-tailed softening; confirmatory-family statement; classifier CV/test tension note
- Supplementary orphans S6/S7/S11/S12 cited; reference author truncation unified;
  redundant author-year parens removed; defensive wording reduced
"""

import re
import sys
from pathlib import Path

MD = Path(__file__).parent / "NAR_MANUSCRIPT_REVISED.md"
text = MD.read_text(encoding="utf-8")

# (old, new) exact-string replacements. Each old string must occur exactly once.
REPLACEMENTS = [
    # --- 1. Table 3 column header ---
    (
        "| Indel Rate | Fusang nRF ↓ (FT2-rel) | FastTree2 nRF ↓ (TRUE-rel) |",
        "| Indel Rate | Fusang nRF ↓ (TRUE-rel) | FastTree2 nRF ↓ (TRUE-rel) |",
    ),
    # --- 2. Table 4 reference frame annotation ---
    (
        "**Table 4. Multi-k ensemble vs single-k configuration (n=200, indel rate=0.02, 30 seeds).**",
        "**Table 4. Multi-k ensemble vs single-k configuration (n=200, indel rate=0.02, 30 seeds). All nRF values are FT2-relative (vs the FastTree2 reference tree), consistent with Table 7.**",
    ),
    (
        "nRF=0: perfect match. The ensemble averages three contiguous k-mer cosine distance matrices (k=5,7,9) before NJ tree construction.",
        "nRF=0: perfect match. All values are FT2-relative (see caption). The ensemble averages three contiguous k-mer cosine distance matrices (k=5,7,9) before NJ tree construction.",
    ),
    # --- 3. Unbold FT2-relative value in TRUE-relative Table 1 ---
    (
        "| 1000 | Indel (0.02) | **0.037 ± 0.006** (k=5,gap2) † | — | — |",
        "| 1000 | Indel (0.02) | 0.037 ± 0.006 (k=5,gap2) † | — | — |",
    ),
    # --- 4. Table first-citation order fixes ---
    (
        "- **Co-phylog** [26]: k-mer frequency + covariance matrix eigenvalues. Tested on both DNA (Table 7) and protein (Table 10) benchmarks.",
        "- **Co-phylog** [26]: k-mer frequency + covariance matrix eigenvalues. Tested on both DNA and protein benchmarks (see Results).",
    ),
    (
        "The multi-k ensemble variant (Table 4) provides a statistically significant improvement over the default configuration.",
        "The multi-k ensemble variant (see below) provides a statistically significant improvement over the default configuration.",
    ),
    (
        "Note that this benchmark uses an independent coalescent simulation protocol, distinct from the guide-tree protocol used in Table 1 and Table 3;",
        "Note that this benchmark (Table 2) uses an independent coalescent simulation protocol, distinct from the guide-tree protocol used in Table 1 and the indel-rate scan in Table 3;",
    ),
    (
        "Table 5 reports an independent pipeline-level validation under a coalescent simulation protocol (single-k NJ: nRF=0.743, multi-k NJ: nRF=0.583), where absolute nRF values are higher for all methods.",
        "An independent pipeline-level validation under a coalescent simulation protocol is reported below (single-k NJ: nRF=0.743, multi-k NJ: nRF=0.583), where absolute nRF values are higher for all methods.",
    ),
    (
        "We compute contiguous k-mer cosine distance matrices for k=5, 7, and 9, then average the three matrices before building a single NJ tree.",
        "We compute contiguous k-mer cosine distance matrices for k=5, 7, and 9, then average the three matrices before building a single NJ tree (Table 4).",
    ),
    (
        "We compared three levels of the Fusang multi-layer pipeline: Level 0",
        "We compared three levels of the Fusang multi-layer pipeline (Table 5): Level 0",
    ),
    (
        "on an independent test set of 88 scenarios spanning two fundamentally different tree types (per-scenario results: Supplementary Table S15):",
        "on an independent test set of 88 scenarios spanning two fundamentally different tree types (Table 6; per-scenario results: Supplementary Table S15):",
    ),
    (
        "We compared Fusang against two established alignment-free phylogenetic methods on simulated indel-rich data (n=200, sub=0.05, indel=0.02, 27 seeds with valid reference trees; per-seed data: Supplementary Table S10):",
        "We compared Fusang against two established alignment-free phylogenetic methods on simulated indel-rich data (n=200, sub=0.05, indel=0.02, 27 seeds, seed set 100–126; Table 7; per-seed data: Supplementary Table S10):",
    ),
    (
        "with trusted reference trees from the SwissTree database.",
        "with trusted reference trees from the SwissTree database (Table 10; per-family data: Supplementary Tables S11–S12).",
    ),
    (
        "the community standard for alignment-free gene tree inference (Zielezinski et al. 2019, *Genome Biology*).",
        "the community standard for alignment-free gene tree inference.",
    ),
    # --- 5. Table 8 indel-row note: split 30-seed values from 112-seed p ---
    (
        "| 200 | 5,gap2 | **0.077 ± 0.018** | 0.080 ± 0.017 | Fusang better (indel, 112 seeds: p=0.052, borderline) |",
        "| 200 | 5,gap2 | **0.077 ± 0.018** | 0.080 ± 0.017 | Fusang better (indel; values from the 30-seed benchmark; the independent 112-seed benchmark yields p=0.052) |",
    ),
    # --- 6. sigma fix ---
    (
        "The per-seed nRF distributions reveal that Fusang variance (σ=0.017) is comparable to FastTree2 variance (σ=0.019), indicating stable performance across replicates.",
        "The per-seed nRF distributions reveal that Fusang's standard deviation (σ=0.016) is comparable to FastTree2's (σ=0.019), indicating stable performance across replicates.",
    ),
    # --- 7. ensemble improvement arithmetic ---
    (
        "- Mean nRF improvement over default spaced: 0.008 (6.7% relative reduction)",
        "- Mean nRF improvement over default spaced: 0.007 (6.3% relative reduction)",
    ),
    # --- 9. 4.6-4.7% -> 4.5-4.8% ---
    (
        "while at biologically typical rates (0.01–0.02) the advantage is modest (+4.6–4.7%)",
        "while at biologically typical rates (0.01–0.02) the advantage is modest (+4.5–4.8%)",
    ),
    # --- 10. 16S gap unification (k=5,gap2 per IMMI L0-1 default provenance) ---
    (
        "in 1.2 seconds without alignment (simplified pipeline, k=5,gap1).",
        "in 1.2 seconds without alignment (simplified pipeline, k=5,gap2).",
    ),
    (
        "**Table 9. Real 16S rRNA validation (n=74, 1.2s, simplified pipeline, k=5,gap1).**",
        "**Table 9. Real 16S rRNA validation (n=74, 1.2s, simplified pipeline, k=5,gap2).**",
    ),
    # --- 11. 16S permutation details + recovery-rate caution ---
    (
        "Fusang recovers significant phylogenetic signal at order and phylum levels across 74 taxa.",
        "P-values are from two-sided permutation tests comparing mean same-group versus different-group tree-based pairwise distances with taxonomic group labels permuted. Fusang recovers significant phylogenetic signal at order and phylum levels across 74 taxa.",
    ),
    (
        "while the FastTree2 tree recovered 10 of 12 (83.3%).",
        "while the FastTree2 tree recovered 10 of 12 (83.3%); given the small number of assessable groups, we treat this recovery-rate comparison as descriptive.",
    ),
    # --- 12. Methods spaced k-mer notation rewrite ---
    (
        "For a DNA sequence S of length L, a spaced k-mer of length k with gap g is defined by a binary pattern of length k + g×(k−1), where k positions are set to 1 (sampled) and g×(k−1) positions are set to 0 (skipped). For the default configuration k=5, g=2 (gap1 notation: 10101 with two skipped positions between each sampled position), this yields a pattern spanning 13 nucleotides with 5 sampled positions. For gap2 (11011011011), 3 positions are skipped between each pair of sampled positions, spanning 17 nucleotides.",
        "For a DNA sequence S of length L, a spaced k-mer of length k with gap g is defined by a binary pattern of length k + g×(k−1), where k positions are set to 1 (sampled) and g×(k−1) positions are set to 0 (skipped); g denotes the number of skipped positions between consecutive sampled positions. The default configuration k=5, g=2 corresponds to the pattern 1001001001001 (5 sampled positions spanning 13 nucleotides), and k=4, g=1 corresponds to 1010101 (4 sampled positions spanning 7 nucleotides). Contiguous k-mers are the special case g=0.",
    ),
    # --- 13. Table 10 pattern labels made consistent ---
    (
        "| **Fusang** k=4,gap1 (1011) | Spaced | **0.239** | 0.118 | 5 |",
        "| **Fusang** k=4,gap1 (1010101) | Spaced | **0.239** | 0.118 | 5 |",
    ),
    (
        "| **Fusang** k=5,gap2 (11011) | Spaced | **0.244** | 0.113 | 5 |",
        "| **Fusang** k=5,gap2 (1001001001001) | Spaced | **0.244** | 0.113 | 5 |",
    ),
    # --- 14. L42 selling-point softening ---
    (
        "The core contribution of this work is therefore a systematic evaluation of **spaced k-mer frequency vector cosine distances** — comparing both spaced and contiguous k-mer patterns across multiple distance metrics — for phylogenetic inference under indel-rich conditions.",
        "The core contribution of this work is therefore a systematic evaluation of **k-mer frequency vector cosine distances** — comparing both spaced and contiguous k-mer patterns across multiple distance metrics — for phylogenetic inference under indel-rich conditions.",
    ),
    (
        "(2) spaced k-mers provide theoretical robustness at high indel rates, though the practical advantage over contiguous k-mers is modest at the tested indel rate (0.02);",
        "(2) spaced k-mers have theoretical motivation for robustness at high indel rates, though the practical advantage over contiguous k-mers is modest at the tested indel rate (0.02);",
    ),
    # --- 15. Contribution 1 rebase to 30-seed L1 vs L0 ---
    (
        "1. **Vector-based phylogenetic inference under indels.** Fusang's multi-k ensemble NJ achieves nRF=0.583±0.045 vs MAFFT+FastTree2 nRF=0.592±0.041 on n=5 valid seeds (p=0.24), with both methods remaining substantially distant from the TRUE tree (nRF≈0.58–0.59). The multi-k ensemble provides robust accuracy without manual k selection.",
        "1. **Vector-based phylogenetic inference under indels.** Fusang's multi-k ensemble NJ reduces topological error by 21.5% relative to the single-k configuration (TRUE-relative nRF=0.583±0.044 vs 0.743±0.046, 30 seeds, Wilcoxon p<0.0001, d=3.55). In a preliminary 5-seed comparison, the ensemble also achieves nRF numerically comparable to MAFFT+FastTree2 (0.583 vs 0.592, p=0.24; full 30-seed Linux validation pending), with both methods remaining substantially distant from the TRUE tree (nRF≈0.58–0.59). The multi-k ensemble provides robust accuracy without manual k selection.",
    ),
    # --- 16. Abstract "simulated" qualifier ---
    (
        "A random forest boundary classifier distinguishes homogeneous from structured datasets (88/88 scenarios, 95% CI [0.958, 1.0]).",
        "A random forest boundary classifier distinguishes homogeneous from structured simulated datasets (88/88 scenarios, 95% CI [0.958, 1.0]).",
    ),
    # --- 17. Winner -> Tie (n.s.) ---
    (
        "| 200 | Clean | 0.102 ± 0.019 (k=5,gap2) | 0.096 ± 0.019 | FT2 (n.s.) |",
        "| 200 | Clean | 0.102 ± 0.019 (k=5,gap2) | 0.096 ± 0.019 | Tie (n.s.) |",
    ),
    (
        "| 200 | Indel (0.02) | **0.077 ± 0.018** (k=5,gap2) | 0.080 ± 0.017 | Fusang (n.s.) |",
        "| 200 | Indel (0.02) | **0.077 ± 0.018** (k=5,gap2) | 0.080 ± 0.017 | Tie (n.s.) |",
    ),
    # --- 18. JSD one-tailed softening ---
    (
        "However, in benchmarking on simulated data (n=200, indel=0.02, 10 seeds, preliminary), cosine distance achieved mean nRF=0.078 vs JSD mean nRF=0.091 — a modest but consistent advantage (Wilcoxon p=0.031, one-tailed).",
        "However, in a small exploratory benchmark on simulated data (n=200, indel=0.02, 10 seeds), cosine distance achieved mean nRF=0.078 vs JSD mean nRF=0.091; this comparison informed the heuristic choice of cosine distance for the simplified pipeline and should not be interpreted as a confirmatory test.",
    ),
    # --- 19. Confirmatory-family statement ---
    (
        "- Benjamini-Hochberg FDR correction as a less conservative alternative",
        "- Benjamini-Hochberg FDR correction as a less conservative alternative. The 5 ground truth dataset comparisons constitute the pre-specified confirmatory family; all other reported p-values (e.g., SwissTree, alignment-free competitor, and pipeline-level comparisons) are secondary and should be interpreted as exploratory",
    ),
    # --- 20. Classifier CV vs test tension note ---
    (
        "The perfect performance across 88 diverse simulated scenarios provides strong evidence that the boundary classifier generalizes beyond its training distribution.",
        "The perfect performance across 88 diverse simulated scenarios provides strong evidence that the boundary classifier generalizes beyond its training distribution. We note that test accuracy (100%) exceeds cross-validation performance on the training set (ROC-AUC=0.84); this gap likely reflects that the E2E test scenarios occupy well-separated regions of the parameter space (purely homogeneous coalescent vs clearly structured phylogenies), and performance on borderline scenarios near the STOP/SPLIT decision boundary may be lower.",
    ),
    # --- 21. Feature ablation future work ---
    (
        "with tree topology metrics providing secondary signal.",
        "with tree topology metrics providing secondary signal. A systematic feature ablation study quantifying each feature group's causal contribution to classification accuracy remains future work.",
    ),
    # --- 22. DCM threshold limitation ---
    (
        "The simplified pipeline is preferred for n ≤ 1000; for n>1000, DCM+EPA provides essential scalability.",
        "The simplified pipeline is preferred for n ≤ 1000; for n>1000, DCM+EPA provides essential scalability. We note that the n=1000 threshold was derived from the DCM degradation analysis at n=200 and from scalability timings; a direct accuracy comparison between the simplified and DCM pipelines in the n=500–1000 transition zone has not been performed.",
    ),
    # --- 23. Reference author truncation unification ---
    (
        "Hadfield, J., Megill, C., Bell, S.M., Huddleston, J., Potter, B. et al. (2018)",
        "Hadfield, J. et al. (2018)",
    ),
    (
        "Hug, L.A., Baker, B.J., Anantharaman, K., Brown, C.T., Probst, A.J. et al. (2016)",
        "Hug, L.A. et al. (2016)",
    ),
    (
        "Zielezinski, A., Girgis, H.Z., Bernard, G., Leimeister, C.-A., Tang, K. et al. (2019)",
        "Zielezinski, A. et al. (2019)",
    ),
    (
        "Ondov, B.D., Treangen, T.J., Melsted, P., Mallonee, A.B., Bergman, N.H. et al. (2016)",
        "Ondov, B.D. et al. (2016)",
    ),
    # --- 24. andi redundant author-year ---
    (
        "**Regarding andi** (Haubold et al. 2015, Bioinformatics):",
        "**Regarding andi** [25]:",
    ),
    # --- 25. Supplementary Figure S6 citation ---
    (
        "and the cases where they diverge are approximately symmetric.",
        "and the cases where they diverge are approximately symmetric (see Supplementary Figure S6 for effect size analyses across benchmarks).",
    ),
    # --- 26. Supplementary Table S7 citation (BAliBASE) ---
    (
        "Fusang achieves competitive performance with 65% of families below nRF 0.5 (median nRF=0.45).",
        "Fusang achieves competitive performance with 65% of families below nRF 0.5 (median nRF=0.45; Supplementary Table S7).",
    ),
    # --- 27. Defensive wording reduction (borderline) ---
    (
        "112-seed post-exclusion benchmark: p=0.052, borderline)",
        "112-seed post-exclusion benchmark: p=0.052)",
    ),
    (
        "achieved p=0.052 (paired Wilcoxon signed-rank test on TRUE-relative nRF, borderline)",
        "achieved p=0.052 (paired Wilcoxon signed-rank test on TRUE-relative nRF)",
    ),
    (
        "112 seeds, paired Wilcoxon p=0.052, borderline).",
        "112 seeds, paired Wilcoxon p=0.052).",
    ),
    (
        "a consistent directional advantage with a small-to-medium effect size (Wilcoxon p=0.052, borderline).",
        "a consistent directional advantage with a small-to-medium effect size (Wilcoxon p=0.052, borderline at α=0.05).",
    ),
    # --- 28. Mash error-tree story dedup (keep full note at Table 2 only) ---
    (
        " An earlier single-seed value (nRF=1.005) reported during development appears to have been based on an incorrectly parsed tree file and is superseded by this benchmark.",
        "",
    ),
    (
        "An earlier single-seed value (nRF=0.162) reported during development appears to have been based on an incorrectly parsed tree file. K-mer cosine distances (Fusang) are the recommended alignment-free approach for both clean and indel-rich data.",
        "K-mer cosine distances (Fusang) are the recommended alignment-free approach for both clean and indel-rich data.",
    ),
]

n_applied = 0
for old, new in REPLACEMENTS:
    count = text.count(old)
    if count != 1:
        print(f"ERROR: expected 1 occurrence, found {count}: {old[:90]!r}")
        sys.exit(1)
    text = text.replace(old, new)
    n_applied += 1

# --- 8. 13.3% -> 13.2% (all occurrences: L250, L261, L263, L711) ---
n_133 = text.count("13.3%")
text = text.replace("13.3%", "13.2%")
print(f"13.3% -> 13.2%: {n_133} occurrences")

MD.write_text(text, encoding="utf-8")
print(f"Applied {n_applied} targeted replacements + {n_133} global replacements")

# ================= self-checks =================
errors = []

# 1. Table 3 header consistent
if "(FT2-rel)" in text:
    errors.append("FT2-rel label still present")

# 2. Table first-citation order monotonic
first_cite = {}
for m in re.finditer(r"Table (\d+)", text):
    t = int(m.group(1))
    if 1 <= t <= 11 and t not in first_cite:
        first_cite[t] = m.start()
# exclude caption self-mentions: find first NON-caption occurrence
lines = text.split("\n")
first_prose = {}
for i, line in enumerate(lines, 1):
    if line.startswith("**Table"):
        continue
    for m in re.finditer(r"Table (\d+)", line):
        t = int(m.group(1))
        if 1 <= t <= 11 and t not in first_prose:
            first_prose[t] = i
order = [t for t, _ in sorted(first_prose.items(), key=lambda kv: kv[1])]
print(f"Table first-prose-citation order: {order}")
expected = list(range(1, 12))
if order != expected:
    errors.append(f"Table citation order not monotonic: {order}")

# 3. Every table cited in prose at least once
for t in expected:
    if t not in first_prose:
        errors.append(f"Table {t} never cited in prose")

# 4. Citations 1-31 all present, monotonic first use
cite_first = {}
for i, line in enumerate(lines, 1):
    if i > 640:  # references section
        break
    for m in re.finditer(r"\[(\d+(?:\s*,\s*\d+)*)\]", line):
        for part in m.group(1).split(","):
            n = int(part.strip())
            if n not in cite_first:
                cite_first[n] = i
missing = [n for n in range(1, 32) if n not in cite_first]
if missing:
    errors.append(f"References never cited: {missing}")
cite_order = [n for n, _ in sorted(cite_first.items(), key=lambda kv: kv[1])]
if cite_order != sorted(cite_order):
    errors.append(f"Citation first-use order not monotonic: {cite_order}")

# 5. Contradiction greps
for pat, label in [
    ("σ=0.017", "sigma=0.017 residual"),
    ("0.008 (6.7%", "0.008/6.7% residual"),
    ("13.3%", "13.3% residual"),
    ("4.6–4.7%", "4.6-4.7% residual"),
    ("k=5,gap1).", "16S gap1 residual (Results/Table 9)"),
    ("(1011)", "Table 10 old pattern 1011"),
    ("(11011)", "Table 10 old pattern 11011"),
    ("one-tailed", "JSD one-tailed residual"),
    ("Haubold et al. 2015, Bioinformatics", "andi author-year residual"),
    ("Zielezinski et al. 2019, *Genome Biology*", "Zielezinski author-year residual"),
    ("Fusang (n.s.)", "Winner Fusang (n.s.) residual"),
]:
    if pat in text:
        errors.append(f"{label}: still present")

# 6. 16S gap1 anywhere?
if re.search(r"16S.*gap1|gap1.*16S", text):
    errors.append("16S gap1 cross-reference remains")

# 7. Borderline count
n_borderline = text.count("borderline")
print(f"'borderline' occurrences: {n_borderline}")
if n_borderline > 4:
    errors.append(f"borderline still {n_borderline}x")

# 8. Abstract word count
abstract = text.split("## ABSTRACT")[1].split("---")[0].strip()
n_words = len(abstract.split())
print(f"Abstract words: {n_words}")
if n_words > 200:
    errors.append(f"Abstract {n_words} words > 200")

# 9. Supplementary orphans check: S6, S7, S11, S12 cited outside the SUPPLEMENTARY MATERIAL list
supp_start = text.find("## SUPPLEMENTARY MATERIAL")
main_body = text[:supp_start]
for s in ["Figure S6", "Table S7", "Tables S11–S12"]:
    if s not in main_body:
        errors.append(f"Supplementary {s} still orphaned in main text")

if errors:
    print("\nSELF-CHECK FAILURES:")
    for e in errors:
        print(f"  - {e}")
    sys.exit(1)
print("\nALL SELF-CHECKS PASSED")
