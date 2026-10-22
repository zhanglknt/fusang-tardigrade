# Skmer inapplicability test record (2026-10-09)

## Setup
- Skmer 3.3.0 (bioconda, conda-forge + bioconda channels), WSL2 Ubuntu-24.04
- Dependencies installed by conda: jellyfish, mash, seqtk
- Input: 27-seed AF benchmark set (seeds 100–126, n=200, L=500 bp, sub=0.05,
  indel=0.02, unaligned FASTA), split into per-sequence files
  (one single-sequence FASTA per taxon, as required by `skmer reference`)

## Commands
```
skmer reference <seed_dir> -t -o seed100            # default k=31
skmer --debug reference <seed_dir> -t -k 21 -o seed100k21   # k=21
```

## Result: FAILS on all inputs (structural inapplicability)
```
[skmer] Estimating coverages using 4 processors...
Traceback (most recent call last):
  File ".../skmer/__main__.py", line 311, in reference
    (name, coverage, genome_length, error_rate, read_length) = result.get(9999999)
  File ".../multiprocessing/pool.py", line 774, in get
    raise self._value
ZeroDivisionError: float division by zero
```

Failure occurs identically at k=31 (default) and k=21, in `estimate_cov` —
Skmer's k-mer depth spectrum coverage model. For 500-bp single-gene assemblies
every k-mer occurs ~once, producing a degenerate depth histogram; the coverage
estimator divides by zero.

## Interpretation
Skmer is designed for **genome skims** (low-coverage whole-genome sequencing
reads), where its coverage/error correction model operates on meaningful k-mer
depth spectra. Gene-length sequences (L=500 bp) are far outside its design
domain, and the tool cannot produce a distance matrix for them. This is a
"tested and structurally inapplicable" outcome, parallel to andi (suffix-array
anchors require genome scale). No distance matrix could be computed, so no nRF
comparison is possible.

## Reproduce
- Split FASTA: see skmer_bench preparation (per-sequence files in seed dirs)
- Env: `conda create -n skmer -c conda-forge -c bioconda skmer`
- Run: `skmer reference <dir> -t -o <prefix>`

---

## Verification on the current INDELible benchmark (2026-10-22, ablation-jsd-worker)

The 2026-10-09 record above used the pre-INDELible "27-seed AF benchmark set"
(seeds 100–126), which no longer exists in the repository. Re-verification on
the CURRENT benchmark data (`indelible_n200_indel/seed100_indel.fasta`,
n=200, indel=0.02, lengths 476–540 bp) found **different behaviour that must
be reflected in the manuscript**:

1. **No crash at k=31 or k=21.** Both runs complete (exit 0). For every one of
   the 200 sequences skmer prints
   `Coverage of seqN is too low, not able to estimate it; no correction applied`
   and proceeds to write a 200x200 distance matrix.
   Mechanism confirmed: every jellyfish depth histogram is a single row
   (e.g. `seq0.hist` = `1 496` at k=31), and `skmer/__main__.py:96`
   (`if len(count) < 3`) catches this gracefully. The record's
   `ZeroDivisionError` comes from `__main__.py:117`
   (`count[ind+1] / count[ind]`), which needs a >= 2-row histogram with a
   zero-count gap; a direct canonical-kmer count shows 0/200 INDELible
   sequences have a repeated 21-mer (or 31-mer), so that branch is unreachable
   on this data.
2. **Version label ambiguity.** `conda list` shows package
   `skmer 3.3.0 pyh086e186_0` in both the old env and a fresh
   `conda create -n skmer33 skmer=3.3.0`, but the installed payload is
   `skmer-3.1.0.dist-info` in both, and the code is byte-identical. The
   bioconda "3.3.0" package ships 3.1.0 code.
3. **The output matrix is severely degraded, so the "inapplicable" conclusion
   stands.** All sequences receive NA coverage (no error correction);
   k=31 caps 45.1% of off-diagonal pairs at 5.0 (k=21: 5.6%).
   NJ trees from these matrices score nRF = 0.3198 (k=31) / 0.3147 (k=21)
   vs the seed-100 true tree — ~2x worse than Fusang (0.1624) and worse than
   every Table 7 competitor except Co-phylog.

Full evidence: `repro_ablations/_skmer_probe/FINDINGS.md` (+ matrices,
histograms, per-sequence inputs in the same directory).

**Manuscript impact**: lines ~384 and ~570 currently say Skmer "fails with a
division-by-zero error" on 500-bp inputs and that "no distance matrix can be
computed". On the INDELible benchmark this is not what happens; the accurate
statement is that the coverage/error-correction model is inoperative on every
sequence (all NA) and the resulting uncorrected distances are dominated by
capped values, giving nRF ≈ 0.32. The conclusion (structural inapplicability
at gene length) is unchanged.
