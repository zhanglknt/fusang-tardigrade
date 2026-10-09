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
