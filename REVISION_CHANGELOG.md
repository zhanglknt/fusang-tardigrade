# Revision Changelog (internal record, not for submission)

## v2.9 (2026-10-09) — v4 expert-panel review fixes (4 reviewers, all Minor Revision)

Resolves all Major and Minor items from the fourth blind-review round
(review_v4_{phylo,stats,format,general}.md; scores 84/85/88/91), plus a
Crossref-audited reference correction discovered during the fix pass.

### Major fixes (cross-reviewer consensus, 7 items)

1. **Abstract p=0.024 overclaim (4/4 reviewers)**: "significantly outperforms
   MAFFT+FastTree2" -> "outperforming MAFFT+FastTree2 on the same benchmark
   (p=0.024; borderline after Bonferroni correction, p_adj=0.071)".
2. **Abstract multi-k claim replaced (stats M1+M2)**: the unstable Table 4
   improvement (p=0.006, not replicated on seed set 100-126) replaced by the
   robust single-k->multi-k gain (0.583 vs 0.743, p<0.0001) from Table 5.
3. **fish "ties MSA+ML" reframed as community-benchmark ceiling (phylo+general+
   stats)**: Abstract, Table 12 finding 1, Limitations, and Practical
   recommendations now state the dataset does not discriminate among top
   methods (n=25, 44 informative splits, all top methods converge at 0.045);
   evidential value = alignment-free + no k selection reaches community optimum.
4. **multi-k seed-set heterogeneity statistically addressed (stats M2)**: new
   "Reproducibility of the ensemble effect across seed sets" paragraph in the
   competitor section: heterogeneity z=-2.71 p=0.007; pooled n=57 analysis
   (mean diff -0.003, paired d=-0.17, 95% CI [-0.45, 0.09], Wilcoxon p=0.19,
   ensemble wins 32/57) computed from per-seed data
   (table7_recomputed_recovered.json + benchmark_multik_ensemble_n200_indel.csv).
   Claim lowered: no general ensemble-vs-default advantage at L=500bp; ensemble
   benefits rest on L=1000bp coalescent gain (d=3.82) and genome-scale rescue.
5. **Discussion n=500/1000 SD+d contradiction resolved (stats M3)**: source data
   (benchmark_n500_clean_30seeds.csv, benchmark_n1000_clean_30seeds.csv)
   confirms Table 1 SDs; Discussion corrected 0.119+-0.020/0.093+-0.015 ->
   0.119+-0.011/0.093+-0.013, 0.115+-0.022/0.091+-0.016 -> 0.115+-0.011/
   0.091+-0.010; d=1.47->1.77 (n=500), d=1.26->2.81 (n=1000); p<0.001 unchanged.
6. **L3 protocol L=1000bp motivation added (general)**: coalescent protocol fixed
   during pipeline development before L1-vs-L3 designed; k=9 non-sparsity
   rationale; identical for superseded 5-seed run and 30-seed rerun; Supp Note
   S8 description now records protocol parameters.
7. **Storage-corruption disclosure precisionized (general+stats m7+format)**:
   Data Availability now lists affected files (5 FASTA + 16 trees, seeds
   109-126), three recovery paths, cross-path byte-identical validation,
   53/54 nRF match vs pre-incident master CSV (one FT2 value 0.0787->0.0838),
   quarantine, checksum log in repro package.

### Reference corrections (Crossref-audited; beyond reviewer findings)

All 31 DOIs verified against Crossref; 7 references had WRONG DOIs/details:
- Ref 5 Wong 2008: DOI 10.1126/science.1146308 (404) -> 10.1126/science.1151532
- Ref 6 Fusang v1: "Zhang, L. et al., 51, 10934-10950, gkad751" (resolves to
  TTD database paper) -> Wang, Z. ... Zhang, L. (10 authors), 51, 10909-10923,
  10.1093/nar/gkad805
- Ref 8 TACOA: "Gkaiogiannis et al. 2016, 17:99" (resolves to Mirnacle) ->
  Diaz, N.N. et al. 2009, 10:56, 10.1186/1471-2105-10-56
- Ref 11: "Bernard et al., NAR 47, W202-W208, gkz331" (resolves to SPADE web
  service) -> Bernard, G. et al., Brief. Bioinform. 20, 426-435,
  10.1093/bib/bbx067 (title corrected: "hierarchical and reticulate
  phylogenomic relationships")
- Ref 25 andi: "31, 1163-1167, btv047" (resolves to KAPPA) -> Haubold, B.,
  Klotzl, F. & Pfaffelhuber, P., 31, 1169-1175, 10.1093/bioinformatics/btu815
- Ref 26 Co-phylog: DOI gkt165 (resolves to Adams hybridization paper) ->
  gkt003
- Ref 27 Lunter: "2006, 7, 320, 10.1186/1471-2105-7-320" (resolves to Ooi et
  al.) -> 2005, 6, 83, 10.1186/1471-2105-6-83
Author truncation unified (format Minor 6): refs 1/15/16/17/20/22/23/24/30
expanded to full author lists (<=10 authors); refs 2/14 remain "et al." (17/19
authors). Audit trail: _ref_audit_crossref.json, _ref_authors_crossref.json.

### Minor fixes

- stats m1: contribution d sign unified to -0.45
- stats m2: Table 9 note now specifies test statistic + 10,000 permutations;
  verification rerun on archived 73-taxon tree reproduces order (7.3%,
  p=0.027) and phylum (3.8%, p=0.003) signal; family labels not archived,
  family row disclosed as from original analysis
- stats m4: cross-validation wording now "benchmark harness ... rather than
  Fusang itself"; NJ attribution softened to "consistent with"
- stats m5: Co-phylog reimplementation direction-conservatism argument added
  (Table 7 finding 1)
- stats m6: "no measurable advantage" -> "no significant difference" + power
  note (n=27, ~73% power for d=0.5)
- stats m8: Table 8 n=200 indel row -> "Fusang better (indel)" + footnote
- general: Abstract classifier sentence + "simulated"; gene-length gap
  candidates added (viral quasispecies, fast-evolving gene families);
  Table 12 kmacs ecoli cell "n/a 3"; fish mechanism marked speculative
- phylo: Introduction L42 spaced-centric bold sentence and finding (1)
  neutralized to "k-mer frequency vectors (spaced and contiguous)";
  "theoretical motivation"/"inherent robustness" -> "hypothesized" (2 sites);
  "the most robust" -> "among the most robust"; FFP/FSWM untested reason added
- format Major 1: word count declaration 9,800 -> 11,200 (actual 11,376)
- format Minor 3: Table 12 footnote 1 anchored to caption "ETE3 nRF^1"
- format Minor 4: supplementary list re-blocked Figure S1-S7 -> Table S1-S16 ->
  Note S1-S10 (was interleaved)
- format Minor 7: TRUE tree -> true tree in running prose (15 sites; definitional
  "(TRUE) tree" and TRUE-relative/TRUE-rel reference-frame labels retained);
  NOT x2 -> not
- format Minor 8: L1-vs-L3 borderline mentions trimmed 5 -> 3 (Abstract,
  contributions pointer, Results); removed from Limitations and Practical recs

### Self-checks (all pass)

Abstract 216 words (<=250); citations 1-31 first-appearance monotonic; Tables
1-12 first-citation monotonic; Supp Figure S1-S7 / Table S1-S16 / Note S1-S10
each grouped and monotonic; no residual superseded values (0.119+-0.020,
d=1.47, btv047, gkt165, gkad751, 1471-2105-7-320, 9,800 words, etc.).
The Table 8 n=100 row value 0.115+-0.022 is a different condition (n=100
clean, seed set 70-99) and is NOT the Discussion n=1000 value - verified
correct.

# Revision Changelog (internal record, not for submission)

### Post-verification fixes (v4-verifier findings, review_v4_verify.md)

1. Ref 7 author list corrected: "Zielezinski, A. and Karlowski, W.M." ->
   "Horwege, S. and Leimeister, C.-A." (Crossref: Morgenstern/Zhu/Horwege/
   Leimeister; Zielezinski/Karlowski are ref 10's authors - copy misplacement).
2. Table 4 d sign unified: "0.54" -> "-0.53" (negative favors ensemble,
   consistent with the paired-d convention and the new L377 paragraph).
3. Heterogeneity stat aligned to precise recompute: z=-2.71/p=0.007 ->
   z~= -2.8 / p~=0.005 (SE-approximation dependent; conclusion unchanged).
4. Ref 21 title completed (": improvements in performance and versatility").
5. Word-count declaration 11,200 -> 11,400 (measured 11,380).
Verifier also independently reproduced pooled n=57 (d=-0.17, p=0.19, 32/57)
and n=500/1000 clean stats (d=1.770/2.807) - all match manuscript values.
Known low-severity residual: Supp Notes S1/S2/S3/S6/S7 have no in-text
citation (pre-existing; supplementary notes are self-standing); classifier
on/off ablation and Co-phylog halfctx=5 DNA scan remain future work.


## v2.8 (2026-10-09) — inherent reservation ③ resolved: real indel-rich genome benchmarks

Reservation ③ (real indel-rich data validation) executed via the AFproject community
benchmark's real whole-genome datasets (trusted published reference trees), replacing the
earlier 2–4 week estimate with a same-day result:

1. **Datasets** (from `real_data/swisstree/afproject_repo/datasets/genome/`, sequences
   downloaded from NCBI eutils by worker; 52/52 accessions validated):
   - Fish mtDNA: 25 mitochondrial genomes (~17 kb), Fischer et al. 2013 reference tree;
     MAFFT alignment gap content 9.1%.
   - E. coli/Shigella: 27 whole genomes (4.4–5.5 Mb), Skippington & Ragan 2011 tree.
2. **Gap-content quantification**: SwissTree 11 protein families measured at median 50%
   gap cells (range 11–81%) — retroactively establishes the existing Table 10 benchmark
   as strongly indel-rich real data.
3. **Fish mtDNA (positive result)**: multi-k cosine ensemble (k=5,7,9 and k=7,9,11
   variants) achieves nRF=0.045 — tied with MAFFT+FastTree2 (MSA+ML) and with kmacs
   (k=10), Mash (k=11,s=5000), Co-phylog (halfctx=5) at their AFproject-published best
   configurations. Gene-scale default k=5 saturates at genome length (nRF=0.545);
   multi-k automatically rescues this (no length-aware k selection needed).
4. **E. coli/Shigella (boundary delineation)**: k-mer cosine plateaus at nRF=0.50 for
   all k≥9 and all multi-k variants (forward-strand k-mers, 64-bit rolling-hash
   implementation cross-validated against the dict implementation); Mash 0.208 (our
   run) / 0.12 (published); published anchors (andi, Co-phylog, phylonium) 0.08.
   Conclusion: shallow recombination-dominated divergence is outside the cosine
   regime — anchor-based or MinHash methods recommended there.
5. **Pipeline cross-validation**: kmacs k=10 and Mash k=11/s=5000 at AFproject-published
   configurations reproduce published scores within one split (0.045 vs 0.05),
   validating our scoring pipeline (ETE3 Tree.compare unrooted nRF, BioPython NJ vs
   AFproject's fneighbor).
6. **Manuscript changes**: new Methods subsection (Real genome-scale validation), new
   Results section + Table 12 (placed after Scalability to keep table order 1–12
   monotonic), Supplementary Table S16 (full results + accessions + gap stats),
   Abstract clause (197 words), Limitations "Simulated-to-real transfer" expanded,
   Practical recommendations +2 bullets, future work item (2) updated, Data
   Availability extended.
7. **Artifacts**: `real_data/afproject_genome/` (fish_mito/, ecoli_shigella/,
   fish_mito_results/, ecoli_results/, ecoli_results_canon/, wsl_out/,
   real_genome_scores.json, Supplementary_Table_S16.md, suite scripts); ete3
   compatibility shim (cgi.py) for Python 3.13.

## v2.7 (2026-10-09) — inherent reservations ① and ② resolved

Reservations ① (L3 Linux 30-seed validation) and ② (kmacs/Skmer head-to-head) from
`review_v3_CLOSURE.md` executed in full; reservation ③ (real indel-rich data) deferred.

1. **L3 Linux rerun (reservation ①)**: MAFFT --auto + FastTree2 `-nt -gtr -nosupport`
   on WSL2 Ubuntu-24.04, all 30 seeds complete (`~/l3run/run_l3.sh`, merged by
   `merge_l3_results.py` into `l3_validation_n200/l3_validation_results.json`;
   pre-Linux backup at `..._pre_linux.json`). **Protocol correction disclosed**: the
   earlier 5-seed Windows run had omitted FastTree2's `-nt` nucleotide flag
   (protein-mode distances on DNA) — superseded and documented in Supplementary
   Note S8. New results: L0=0.743±0.046, L1=0.583±0.044, L3=0.601±0.055 (all n=30);
   **L1 vs L3: paired Wilcoxon p=0.024, paired d=−0.45, L1 wins 19/30** (borderline
   after Bonferroni ×3, p_adj=0.071; L1 vs L0 paired d=3.82, replacing pooled-SD
   d=3.55 per the manuscript's own paired-d convention). Story upgraded from
   "preliminary n=5 numerical similarity" to "full 30-seed significant L1 advantage".
   Timing corrected: L1 ~15s vs L3 ~31s per seed on Linux (old text: 175s Windows).
   Table 5 caption now discloses L=1000 bp (coalescent protocol).
2. **kmacs benchmark (reservation ②)**: original 2014 source recovered from Wayback
   Machine (kmacs.gobics.de dead), compiled on WSL2 (`-Wno-narrowing`); k-scan
   k∈{3,5,10} on 27 seeds (100–126). Best k=3: vs FT2 0.177±0.024, vs TRUE
   0.168±0.022 — significantly worse than Fusang (p=5.6×10⁻⁶, paired d=2.54) but far
   better than Co-phylog. Added to Methods comparison list, Table 7 (new row +
   finding), Limitations, Supp Table S10 caption.
3. **Skmer test (reservation ②)**: Skmer 3.3.0 (conda) structurally inapplicable at
   gene length — ZeroDivisionError in coverage estimation on all 27 seeds at k=31
   and k=21 (`skmer_inapplicability_record.md`). Documented in Results ("Regarding
   Skmer") and Limitations; parallels andi.
4. **Table 7 rebuilt on a single validated seed set**: forensic provenance check
   revealed the v2.6 Table 7 was a chimera — Co-phylog/k5/k7/Fusang values came from
   a `table8_definitive.log` run on seeds 227–253, the multi-k row was copied from
   the Table 4 benchmark (seeds 230–259), while the caption claimed seed set 100–126.
   All rows recomputed on seeds 100–126 (matching the caption and Supp Table S10) via
   `recompute_table7_recovered.py` on recovered benchmark data: Co-phylog 0.408±0.021,
   kmacs 0.177±0.024, k5 0.104±0.018, k7 0.107±0.019, Fusang 0.108±0.020, multi-k
   0.111±0.018. Consequent narrative change: contiguous k=5 vs spaced is now
   n.s. (p=0.26, d=−0.21) instead of "slightly outperforms (p=0.0002)" — more coherent
   with the protein-domain result. Table 4 caption now states its seed set (230–259).
5. **Benchmark data recovery (disk corruption)**: 5 FASTA (seeds 111/114/115/117/121)
   + 16 tree files (FT2 109–126 range, TRUE 110/121/126, Fusang 114/115/120/125) were
   UTF-16/binary garbage. Recovered via three independent validated paths
   (`recover_benchmark_data.py`): ungap-reconstruction from aligned FASTA,
   deterministic regeneration (`gen_test_data_indel.py`, byte-identical), FastTree
   rerun on intact alignments. 53/54 per-seed nRF values match the master CSV
   (only seed109 FT2 differs by 2 bipartitions: 0.0787→0.0838). Corrupted originals
   quarantined in `corrupt_fasta_quarantine/`. Disclosed in Data Availability.
6. **References renumbered**: Skmer [29] (now first cited in Results), PatternHunter
   [30] (first cited in Discussion) — order-preserving swap; citation sequence
   verified monotonic 1–31.
7. **Self-checks**: Abstract 199 words (≤200); citations 1–31 monotonic; tables 1–11
   monotonic; no residual old values (0.592/0.419/0.099/d=3.55/p=0.0002/preliminary).

## v2.6 (2026-10-09) — post-verification residual fixes

Following the 4-reviewer verification round (`review_v3_verify_{phylo,stats,format,general}.md`),
all remaining one-line/one-paragraph items were fixed via `fix_v26_master.py`:

1. **Table 3 column header** (former L254): `(FT2-rel)` → `(TRUE-rel)`. This was a label
   residual from the v2.5 reference-frame correction; the underlying data were always
   TRUE-relative (source-verified in `run_indel_benchmark_v9.py:135-160`, both methods
   scored against `seed{seed}_indel_true.nwk`).
2. **Table first-citation order** made monotonic (1→11) by citation edits only — no table
   renumbering: removed Methods forward references to Tables 7/10; added in-prose citations
   for Tables 2, 4, 5, 6, 7, 10 at their own Results sections.
3. **Table 4 reference-frame annotation** added (FT2-relative, consistent with Table 7).
4. **Table 8 indel row note** split: 30-seed values vs independent 112-seed p=0.052.
5. **σ=0.017 → σ=0.016** (Fusang per-seed SD, now consistent with the five other occurrences).
6. **Arithmetic rounding**: ensemble improvement 0.008/6.7% → 0.007/6.3%;
   indel-scan advantage 13.3% → 13.2% (4 occurrences); 4.6–4.7% → 4.5–4.8%.
7. **16S configuration contradiction** unified to k=5,gap2 (Methods/Results/Table 9).
   Provenance: the 74-taxa IMMI L0-1 run used the pipeline default k=5,gap2; the "gap1"
   in Results/Table 9 was the erroneous instance. Permutation-test description and
   recovery-rate caution added to Table 9.
8. **Methods spaced k-mer notation** rewritten (pattern 1001001001001 for k=5,g=2;
   1010101 for k=4,g=1); Table 10 pattern labels made consistent.
9. **Narrative**: L42 core-contribution bolding changed to "k-mer frequency vector cosine
   distances"; "theoretical robustness" → "theoretical motivation"; contribution 1 rebased
   to the 30-seed L1-vs-L0 result with the n=5 L1-vs-L3 comparison marked preliminary;
   Abstract classifier claim qualified with "simulated"; Winner column "Fusang/FT2 (n.s.)"
   → "Tie (n.s.)".
10. **Statistics reporting**: JSD 10-seed comparison relabeled exploratory/heuristic
    (one-tailed test claim removed); confirmatory family (5 ground-truth dataset
    comparisons) declared, all other p-values labeled exploratory; classifier CV(0.84)
    vs test(100%) gap discussed; feature-ablation and DCM-threshold limitations added;
    27-seed AF benchmark clarified as a designed seed set (100–126), not 30 minus 3.
11. **Format**: references 1, 2, 14, 24 author truncation unified to lead-author et al.;
    redundant author-year parentheticals removed (andi, SwissTree); supplementary orphans
    S6/S7/S11–S12 now cited in main text; "borderline" reduced 8→4 occurrences.
12. **Table 1 n=1000 indel row** unbolded (FT2-relative value in a TRUE-relative table).

## v2.5 (2026-10-08) — major revision after 4-expert panel (avg 60.3/100)

- Reference-frame labels corrected throughout after source-code forensic verification
  proved both Fusang and FastTree2 columns in the 130-seed benchmark are TRUE-relative
  (`calc_nrf_simple.py` recomputation on seed100 matched the master CSV to 6 decimals).
- **SwissTree Wilcoxon p recomputation**: the Fusang k=4,gap1 vs Co-phylog halfctx=5
  Wilcoxon p was corrected from 0.006 (v2.4) to 0.014 after recomputing from per-family
  nRF values (paired t-test p=0.005 unchanged). The 1.5× advantage cited in the Abstract
  corresponds to the halfctx=11 comparison (Wilcoxon p=0.006), which is unchanged.
  Recomputation script output archived with benchmark data.
- Mash narrative rewritten to the honest version (both degrade severely; no matched-seed
  test; earlier single-seed values disclosed as parse-error artifacts and superseded).
- References renumbered to citation order (31 refs); tables renumbered 1–11;
  Figures 3–6 and previously orphaned supplementary items cited.
- 48 prose replacements (P0 mechanical + P1 narrative) via `fix_v25_master.py`.
