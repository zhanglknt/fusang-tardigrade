# Supplementary Table S16. Real genome benchmarks (AFproject community datasets)

## Datasets

| Dataset | n | Sequence length | Reference tree | Indel content |
|---------|---|-----------------|----------------|---------------|
| Fish mtDNA | 25 | ~16.5-18.3 kb (mitochondrial genomes) | Fischer et al. 2013 (AFproject) | MAFFT alignment: 9.1% gap cells (n=25, 18,284 columns) |
| E. coli/Shigella | 27 | 4.4-5.5 Mb (whole genomes) | Skippington & Ragan 2011 (AFproject) | gene gain/loss dominated (whole-genome; MSA infeasible) |
| SwissTree (Table 10) | 11 families | 259-1,707 aa | SwissTree reference trees | median 50% gap cells (range 11-81%), per family below |

## SwissTree per-family gap fractions (MAFFT --auto)

| Family | n | Aligned length | Gap fraction |
|---------|---|----------------|--------------|
| ST001 | 98 | 524 | 0.362 |
| ST002 | 54 | 1127 | 0.489 |
| ST003 | 49 | 563 | 0.112 |
| ST004 | 115 | 930 | 0.517 |
| ST005 | 29 | 488 | 0.324 |
| ST007 | 60 | 687 | 0.809 |
| ST008 | 42 | 860 | 0.678 |
| ST009 | 39 | 259 | 0.572 |
| ST010 | 34 | 428 | 0.502 |
| ST011 | 159 | 1707 | 0.682 |
| ST012 | 21 | 732 | 0.451 |

## Per-method results (ETE3 Tree.compare, unrooted nRF, max_rf = 2(n-3))

### Fish mtDNA (n=25, 44 informative splits)

| Method | Configuration | nRF | RF splits | Published (AFproject) |
|--------|---------------|-----|-----------|----------------------|
| Multi-k cosine | k=5,7,9 avg | 0.0455 | 2/44 |  |
| Multi-k cosine | k=7,9,11 avg | 0.0455 | 2/44 |  |
| K-mer cosine | k=5 contiguous (gene-scale default) | 0.5455 | 24/44 |  |
| K-mer cosine | k=7 contiguous | 0.2273 | 10/44 |  |
| K-mer cosine | k=9 contiguous | 0.0909 | 4/44 |  |
| K-mer cosine | k=11 contiguous | 0.0909 | 4/44 |  |
| Fusang spaced | k=5,gap2 (default) | 0.5455 | 24/44 |  |
| kmacs | k=3 | 0.1818 | 8/44 |  |
| kmacs | k=10 (AFproject config) | 0.0455 | 2/44 | 0.05 |
| Mash | k=21, s=10000 | 0.0909 | 4/44 |  |
| Mash | k=11, s=5000 (AFproject config) | 0.0455 | 2/44 | 0.05 |
| Co-phylog (Python impl.) | halfctx=5 (AFproject config) | 0.0455 | 2/44 | 0.09 |
| Co-phylog (Python impl.) | halfctx=9 | 0.5455 | 24/44 |  |
| MAFFT+FastTree2 | MSA+ML (MAFFT --auto; FastTree -nt -gtr -nosupport) | 0.0455 | 2/44 |  |

### E. coli/Shigella (n=27, 48 informative splits)

| Method | Configuration | nRF | RF splits | Published (AFproject) |
|--------|---------------|-----|-----------|----------------------|
| Multi-k cosine | k=5,7,9 avg | 0.5000 | 24/48 |  |
| Multi-k cosine | k=5,7,9,11 avg | 0.5000 | 24/48 |  |
| Multi-k cosine | k=7,9,11,13 avg | 0.5000 | 24/48 |  |
| K-mer cosine | k=5 contiguous | 0.8750 | 42/48 |  |
| K-mer cosine | k=7 contiguous | 0.7083 | 34/48 |  |
| K-mer cosine | k=9 contiguous | 0.5000 | 24/48 |  |
| K-mer cosine | k=11 contiguous | 0.5000 | 24/48 |  |
| K-mer cosine | k=13 contiguous | 0.5000 | 24/48 |  |
| Mash | k=21, s=10000 | 0.2083 | 10/48 | 0.12 |
| andi | published | — | — | 0.08 |
| Co-phylog (original binary) | published | — | — | 0.08 |
| phylonium | published | — | — | 0.08 |
| Skmer | published | — | — | 0.17 |
| FSWM | published | — | — | 0.17 |
| FFP | published | — | — | 0.21 |
| spaced words | published | — | — | 0.33 |

## NCBI accessions

Fish mtDNA: NC_009057, NC_009058, NC_009059, NC_009060, NC_009062, NC_009063, NC_009064, NC_009065, NC_009066, NC_009067, NC_009459, NC_010205, NC_011168, NC_011169, NC_011170, NC_011171, NC_011177, NC_011179, NC_012055, NC_013564, NC_013577, NC_013663, NC_013750, NC_018814, NC_018815

E. coli/Shigella: NC_011742, NC_007779, NC_011750, NC_008563, NC_010468, NC_011745, NC_002695, NC_011601, NC_007606, NC_004741, NC_009800, NC_010658, NC_008258, NC_004431, NC_007384, NC_009801, NC_011741, NC_004337, NC_010498, NC_007613, NC_011415, NC_011751, NC_007946, NC_002655, NC_000913, NC_011748, NC_008253

## Protocol notes

- All NJ trees: BioPython DistanceTreeConstructor (AFproject used PHYLIP fneighbor; kmacs k=10 and Mash k=11 reproduce published scores within one split, attributable to the NJ implementation).
- kmacs and Mash are original binaries; Co-phylog is the Python reimplementation used throughout this work (Supplementary Note S4).
- k-mer cosine on E. coli/Shigella used forward-strand k-mers (consistent with all other benchmarks in this work), computed via 64-bit rolling hashes.
- MAFFT+FastTree2 on WSL2 Ubuntu-24.04. Gap fractions measured on MAFFT --auto alignments.
- Data generated 2026-10-09; scripts: real_data/afproject_genome/{fish_mito_suite,ecoli_shigella_suite,ete3_scoring}.py, run logs and distance matrices included in the reproducibility package.- kmacs (k=3) on E. coli/Shigella (27 x ~4.6 Mb, single-threaded) did not complete within our compute budget (>40 min without output) and is therefore not reported for that dataset.
