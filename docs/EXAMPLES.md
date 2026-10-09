# Fusang Worked Examples

Three end-to-end scenarios with commands, expected runtimes, and expected outputs. All timings are from a 4-core workstation; your numbers will scale roughly with core count (`-t`).

---

## Example 1 — Gene family, 200 sequences (~1.5 kb, indel-rich)

The classic use case: homologous gene sequences with enough indels that alignment is unreliable.

```bash
python fusang_v2.py -i gene_family_200.fasta -o gene_family_200.nwk
```

What happens:

```
[MAIN]   read 200 sequences, max length 1531 bp
[MAIN]   auto-selected k=5, gap=gap2 (n=200 taxa)
[MAIN]   auto mode resolved: n=200 → build=default, nni=refine
[MAIN]   auto: n=200 ≤ 500, using simplified pipeline
```

| | |
|---|---|
| **Expected time** | ~1.3 s (NJ) |
| **Output** | `gene_family_200.nwk` — 200-taxon Newick tree |
| **Expected accuracy** | nRF ≈ 0.10–0.11 vs. reference (indel-rate 0.02 simulations) — ties FastTree2-with-MAFFT, without the alignment step |

Optional comparisons:

```bash
# FastME instead of NJ (~0.4 s; slightly lower accuracy on k-mer distances)
python fusang_v2.py -i gene_family_200.fasta -o tree_fastme.nwk --tree_method fastme

# Multi-k ensemble: small consistent accuracy gain, ~3 s
python fusang_v4_dahp_v1.py gene_family_200.fasta --v3 --output tree_multik.nwk
```

---

## Example 2 — 25 mitochondrial genomes (deep divergence, ~16 kb each)

Genome-scale sequences with deep divergence: short k-mers saturate, so use the multi-k ensemble.

```bash
python fusang_v4_dahp_v1.py mito_genomes_25.fasta --v3 --output mito_tree.nwk
```

Why `--v3` here: at ~16 kb with deep divergence, k=5 cosine distances compress toward saturation. The ensemble averages k=5,7,9 distances, and the larger-k components retain resolution — no manual k tuning needed.

| | |
|---|---|
| **Expected time** | ~2–4 s |
| **Output** | `mito_tree.nwk` — 25-taxon Newick tree |
| **Note** | Single-k default (`fusang_v2.py`) also runs, but expect lower resolution on deep internal branches |

To see the saturation effect yourself, compare:

```bash
python fusang_v2.py -i mito_genomes_25.fasta -o mito_single_k.nwk   # k=4,gap1 auto
python fusang_v4_dahp_v1.py mito_genomes_25.fasta --v3 --output mito_multik.nwk
python calc_nrf_simple.py mito_single_k.nwk mito_multik.nwk         # topological difference
```

---

## Example 3 — 10,000 taxa (large-scale benchmark)

The scalability showcase: a full tree of life–sized unaligned dataset in under a minute.

```bash
python fusang_v2.py -i taxa_10000.fasta -o taxa_10000.nwk -t 4
```

What happens:

```
[MAIN]   read 10000 sequences, ...
[MAIN]   auto-selected k=5, gap=gap2 (n=10000 taxa)
[MAIN]   auto-selected max_group=5000 (n=10000 taxa, 2 backbone taxa)
[MAIN]   auto: n=10000 > 500, using DCM pipeline
[DCM] Step 1..6: clustering → group NJ → backbone → EPA grafting → merge
```

| | |
|---|---|
| **Expected time** | ~54 s (4 cores) |
| **Peak memory** | ~0.6 GB |
| **Output** | `taxa_10000.nwk` — 10,000-taxon Newick tree |

Notes:

- The DCM pipeline is automatic; `--max_group` is auto-scaled to keep ≤ 2 backbone taxa (a stability constraint — 3+ backbone taxa make NJ/FastME topology on k-mer distances unstable).
- For finer grafting quality at large n, add `--auto_group_method epa_improved` (NNI refinement after EPA grafting; costs extra time).
- Beyond ~10K taxa, add `--use_minhash` for MinHash-LSH pre-clustering.
- Compare against `--simple` (forced simplified pipeline) only for benchmarking — at n>500 the DCM pipeline is both faster and more accurate.

---

## Reproducing manuscript results

All benchmark data, figure scripts, and frozen code are in [`repro_package/`](../repro_package/). See `repro_package/README.md` for the exact commands behind each manuscript table and figure.
