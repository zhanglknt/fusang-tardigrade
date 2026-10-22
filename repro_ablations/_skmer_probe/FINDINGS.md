# Skmer verification on INDELible data (2026-10-22, ablation-jsd-worker)

> **Metric update (2026-10-23, team lead).** Every nRF number quoted in this
> document (Skmer 0.3198 at k=31 / 0.3147 at k=21, Fusang 0.1624) was computed
> with the **superseded root-dependent metric** `calc_nrf_simple.py`, which
> captures only 64 of the 197 non-trivial bipartitions of an n=200 tree and
> systematically underestimates RF. Those values are **withdrawn**. Under the
> current standard root-independent unrooted RF (`nrf_standard.py`,
> `repro_ablations/skmer_probe_nrfstandard.py`, shipped in this package) the
> same trees score:
>
> | tree | withdrawn (root-dependent) | current (standard unrooted RF) |
> |---|---|---|
> | Skmer k=31 | 0.3198 | **0.7817** |
> | Skmer k=21 | 0.3147 | **0.6142** |
> | Fusang, same seed | 0.1624 | **0.2741** |
>
> The manuscript quotes the current column and explicitly discloses the
> withdrawal of the old one. **Nothing else in this document is affected by
> the metric choice** — the mechanism findings (skmer completes, the coverage
> model is inoperative on every sequence, the fraction of capped sentinel
> distances, the version-label audit, the unreachable division-by-zero
> branch) are metric-independent and stand as written.

Verification of `Fusang-main/skmer_inapplicability_record.md` (dated 2026-10-09),
which claims Skmer v3.3.0 crashes with `ZeroDivisionError` in `estimate_cov` on
500-bp single-gene inputs at k=31 and k=21, making a distance matrix impossible.

## What was re-run

- Input: `indelible_n200_indel/seed100_indel.fasta` (n=200, indel=0.02, INDELible
  v1.03 GTR+G, lengths 476–540 bp), split into 200 per-sequence FASTA files
  (`seed100/seq0.fa` … `seed100/seq199.fa`), as `skmer reference` requires.
- Command: `skmer reference seed100 -t -o seed100_probe` (default k=31),
  WSL2 Ubuntu-24.04, env `/opt/miniconda3/envs/skmer`.
- k=21 rerun: `skmer reference seed100 -t -k 21 -o seed100_probe_k21`.

## Findings

### 1. Version label vs actual code

- `conda list` in both the existing `skmer` env and a freshly created
  `skmer33` env (`conda create -n skmer33 -c conda-forge -c bioconda skmer=3.3.0`)
  reports package `skmer 3.3.0 pyh086e186_0`.
- The installed payload in BOTH envs is `skmer-3.1.0.dist-info`, and the two
  `__main__.py` files are byte-identical (`diff` → IDENTICAL).
- Conclusion: the bioconda package labelled 3.3.0 ships code that reports itself
  as 3.1.0; "3.3.0" is the conda package label, not the code version. The record's
  "Skmer 3.3.0" is consistent with this (same package, same code).

### 2. The recorded crash does NOT reproduce on INDELible data (k=31)

The k=31 run COMPLETED (exit 0, ~46 min). Per-sequence coverage estimation
printed, for every one of the 200 sequences:

    Coverage of seqN is too low, not able to estimate it; no correction applied

and the run went on to sketching and distance estimation, writing
`seed100_probe.txt` (200x200 distance matrix).

### 3. The underlying mechanism IS confirmed (degenerate k-mer depth spectrum)

Every `library/seqN/seqN.hist` file contains exactly one row, e.g.
`seq0.hist` = `1 496` — all 496 k-mers observed exactly once. In
`skmer/__main__.py:96` the guard `if len(count) < 3` catches this case and
returns `coverage=NA` ("too low" branch) instead of crashing.

The record's `ZeroDivisionError` comes from the other branch (`__main__.py:117`):
`cov = ... * count[ind+1] / count[ind]`, which divides by a histogram bin that
can be zero when the depth histogram has >= 2 rows with a gap. On the record's
input (the pre-INDELible "27-seed AF benchmark set, seeds 100–126", which no
longer exists in the repo) at least one sequence evidently had such a histogram.
On the current INDELible data all 200 sequences have single-row histograms, so
the crash branch is never reached.

### 4. The distance matrix is produced but is severely degraded

`seed100_probe.txt` (k=31): 200x200, symmetric, zero diagonal.
- 45.1% of off-diagonal pairs are capped at exactly 5.0 (17962/39800).
- Off-diagonal: min 0.000, median 0.168, mean 2.318 (inflated by the caps).
- NJ tree (fusang_v2.nj_tree) scored against `seed100_indel_true.nwk`:
  **nRF = 0.3198** — ~2x worse than Fusang on the same seed (0.1624).

### 4b. k=21 behaves the same way (empirically confirmed)

`seed100_probe_k21.txt` (k=21): run COMPLETED (exit 0), same all-sequences
"coverage too low / no correction applied" warnings. Matrix:
- 5.6% of off-diagonal pairs capped at 5.0 (fewer caps because k=21 shares more
  k-mers), median 0.1515.
- NJ nRF = **0.3147** — as bad as k=31 despite the mostly non-capped matrix.

A direct Python count of canonical 21-mer depth spectra over all 200 sequences
predicted this: 0/200 sequences have a repeated 21-mer, so every depth histogram
is single-row and the `len(count) < 3` guard catches all of them — the crash
branch is unreachable on this data at either k.

### 5. Implication for the manuscript

The claim "division-by-zero ... no distance matrix can be computed" is NOT
reproducible on the current INDELible benchmark: skmer completes and writes a
matrix, but the coverage/error-correction model is inoperative on every
sequence (all NA), and 45% of pairwise distances hit the 5.0 cap. The
conclusion "structurally inapplicable to gene-length data" survives — and the
evidence is arguably stronger (it runs and produces unusable distances) — but
the mechanism sentence must be rewritten to describe the actual current
behaviour, and the manuscript should state the skmer code version precisely
(bioconda package label 3.3.0 / payload 3.1.0).

## Artifacts

- `seed100_probe.txt`        skmer k=31 distance matrix (200x200, TSV)
- `seed100_probe_k21.txt`    skmer k=21 distance matrix (200x200, TSV)
- `seed100_skmer_k31.npy` / `seed100_skmer_k21.npy`    parsed float64 matrices
- `seed100_skmer_taxa.json`  taxon order for the k=31 matrix
- `seed100/`                 200 per-sequence input FASTAs
- `library/`                 skmer's per-sequence jellyfish/mash intermediates
  (the `*.hist` files are the direct evidence for the degenerate spectra;
  currently holds the k=21 run's intermediates — the k=31 histograms showed
  the same single-row pattern, e.g. `seq0.hist` = `1 496`)
