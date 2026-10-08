# Revision Changelog (internal record, not for submission)

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
