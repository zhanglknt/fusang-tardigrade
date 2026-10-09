# Revision Changelog (internal record, not for submission)

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
