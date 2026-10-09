# Fusang Tutorial

End-to-end guide: installation → your first tree → adaptive parameters → multi-k ensemble → large datasets → the boundary classifier → FAQ.

---

## 1. Installation

Requirements: Python 3.8+ (3.10+ recommended), pip.

```bash
git clone https://github.com/zhanglknt/fusang-tardigrade.git
cd fusang-tardigrade
pip install -r requirements.txt
```

Dependencies: `numpy`, `numba`, `biopython`, `scipy` (plus `flask`/`streamlit` if you want the web UI).

Optional tools, only needed for specific modes:

| Tool | Needed for | Bundled? |
|---|---|---|
| FastME | `--tree_method fastme` | ✅ Yes — `fastme_bin/fastme.exe` (Windows) and `fastme_bin/fastme_linux` |
| MAFFT + FastTree2 | MHL Level 3 (`fusang_mhl_main.py` without `--no-l3`) | ❌ Install separately |

Verify the install:

```bash
python fusang_v2.py --help
```

---

## 2. Your first tree

Take any FASTA of unaligned sequences — gene-family homologs, organelle genomes, protein sequences:

```
>seq1
ACGTTGCAAGTC...
>seq2
ACGTGCAAGTC...
...
```

Run:

```bash
python fusang_v2.py -i my_sequences.fasta -o my_tree.nwk
```

Typical console output:

```
============================================================
  Fusang: Tardigrade Edition -- Scalable Phylogenetic Inference
============================================================
[MAIN] Step 1: reading FASTA my_sequences.fasta...
[MAIN]   read 50 sequences, max length 1482 bp
[MAIN]   auto-selected k=4, gap=gap1 (n=50 taxa)
[MAIN]   auto: n=50 ≤ 500, using simplified pipeline
...
[MAIN] Done. Tree written to my_tree.nwk
```

The output is a standard Newick tree:

```
((seq1:0.012,seq2:0.014):0.031,(seq3:0.018,...
```

Open it in FigTree, iTOL, MEGA, or `ete3`/`dendropy` in Python.

That's it — no alignment, no substitution model, no k selection.

---

## 3. How adaptive parameters work

Fusang chooses k and gap pattern from the number of taxa (n) unless you override them:

| n taxa | k | gap | Rationale |
|---|---|---|---|
| ≤ 100 | 4 | gap1 | Smaller feature space stays informative when distances are short |
| > 100 | 5 | gap2 (`1001001001001`) | Near-universal optimum across scales in benchmarks |

The pipeline itself also adapts:

| n taxa | Pipeline |
|---|---|
| ≤ 500 | Simplified: k-mer → cosine distance → NJ directly (DCM adds no accuracy at this scale) |
| > 500 | DCM: decompose into overlapping groups, build subtree backbones, graft with EPA-style placement |

And the build mode adapts: n < 100 → plain NJ; n ≥ 100 → NJ + BME refinement (`--mode auto`, the default).

Overriding is one flag each:

```bash
python fusang_v2.py -i in.fasta -o out.nwk --kmer_k 7 --kmer_gap gap3
```

When would you override? Rarely. Indel-heavy data at n≤200 can gain a few percent from `gap3`; very high divergence can favor larger k. Benchmarks first — see `docs/EXAMPLES.md`.

---

## 4. Multi-k ensemble: when and why

Single-k distances can saturate: for long, divergent sequences (genome scale), short k-mers accumulate so many repeated hits that cosine distances compress toward their maximum and lose resolution. Instead of guessing the right k, average several:

```bash
python fusang_v4_dahp_v1.py genomes.fasta --v3 --output tree.nwk
# custom ensemble:
python fusang_v4_dahp_v1.py genomes.fasta --v3 --ks 5,7,9,11 --fusion average --output tree.nwk
```

Rules of thumb:

| Data | Recommendation |
|---|---|
| Gene-length sequences (< ~10 kb), n ≤ 500 | Single-k default (`fusang_v2.py`) is fine |
| Gene families with mixed divergence levels | `--v3` (k=5,7,9) — small but consistent nRF gain over single-k |
| Organelle / genome-scale sequences, deep divergence | `--v3` strongly recommended — the k=7/9 components rescue k=5 saturation automatically |
| Protein sequences | Default k=4,gap1 (auto-selected at small n); `--v3` also applicable |

---

## 5. Large datasets (n > 500, up to 10,000+)

Nothing extra to do — the DCM pipeline engages automatically:

```bash
python fusang_v2.py -i big.fasta -o big_tree.nwk -t 8
```

What happens under the hood:

1. Sequences are clustered into overlapping groups (overlap ratio `--overlap 0.15`), with group size auto-scaled so the backbone has ≤ 2 supertaxa — a stability constraint learned from benchmarking.
2. Each group gets a k-mer NJ subtree.
3. Subtrees are merged into a backbone; remaining taxa are grafted EPA-style (`--auto_group_method epa_improved` adds NNI refinement after grafting).

Reference timing (4-core workstation): n=1,000 ≈ 8 s; n=10,000 ≈ 54 s, ~0.6 GB RAM.

For tens of thousands of taxa, add `--use_minhash` for MinHash-LSH coarse clustering before distance computation.

Practical tips:

- Set `-t` to your physical core count.
- Force the simplified pipeline for comparison runs with `--simple`; force DCM at small n (benchmarking only) with `--no-simple`.
- `--tree_method fastme` is ~3× faster than NJ at n=200 (0.4 s vs 1.3 s) but may degrade k-mer distance accuracy; NJ is the default for a reason. Use FastME when speed matters more than the last few percent of accuracy.

---

## 6. The boundary classifier

`fusang_mhl_main.py` ships a pre-trained random-forest classifier (`fusang_mhl/models/boundary_rf.pkl`) that decides whether a dataset looks **homogeneous** (single cloud in k-mer space — simple methods suffice) or **structured** (distinct clusters — hierarchical decomposition pays off), and routes it through the IMMI levels accordingly.

```bash
# Full IMMI pipeline with boundary detection
python fusang_mhl_main.py sequences.fasta -o tree.nwk

# Explicit model path / skip levels
python fusang_mhl_main.py sequences.fasta -o tree.nwk \
    --boundary-model fusang_mhl/models/boundary_rf.pkl --no-l3 --debug
```

- `--no-l2`: skip DAHP-V2 backbone refinement
- `--no-l3`: skip MSA+ML within reliable clusters (use this if MAFFT/FastTree2 are not installed)
- `--debug`: per-level timing and decisions

If the classifier reports your dataset as structured in a way consistent with recombination or structured genes (rRNA), treat the resulting tree with caution — see the applicability table in the README.

---

## 7. Evaluating your tree

Fusang outputs a topology with branch lengths in cosine-distance units. To sanity-check:

- Compare against a reference method on the same data (`calc_nrf_simple.py` computes normalized Robinson–Foulds distance between two Newick trees).
- For publication-grade support values, bootstrap by resampling columns is not meaningful for alignment-free input; instead, replicate over k and gap settings and check clade stability, or run the multi-k ensemble and compare with single-k results.
- For downstream dating/ML refinement, import the Newick into standard tools.

---

## 8. FAQ

**"FastME not found" / FastME errors on Windows.**
You don't need WSL. Fusang searches in order: `fastme_bin/fastme.exe` (bundled, native Windows) → `fastme` on PATH → `fastme_bin/fastme_linux` via WSL. If you see this error, check that `fastme_bin/fastme.exe` exists in the repo root; it is included in the distribution. If you deliberately want the WSL binary, set `FASTME_WSL_DISTRO=Ubuntu-24.04`.

**Do I need to align my sequences first?**
No — and you shouldn't. Fusang is designed for unaligned input. Aligned input works but wastes the alignment effort.

**DNA or protein?**
Both. The same k-mer machinery applies; adaptive defaults handle typical cases for each.

**How big can I go?**
Benchmarked to 10,000 taxa (~54 s). Beyond that, use `--use_minhash` and more threads.

**When should I NOT use Fusang?**
Structured rRNA genes (16S/23S) and strain-level recombination-dominated genomes (e.g. *E. coli*/*Shigella*). In both cases the k-mer independence/additivity assumptions break down; use alignment-based methods.

**Which entry point should I use?**
`fusang_v2.py` for almost everything. `fusang_v4_dahp_v1.py --v3` when you want the multi-k ensemble (genome-scale or deep divergence). `fusang_mhl_main.py` when you want the full IMMI decision pipeline with boundary classification.

**GUI?**
`python fusang_webapp.py`, then open http://localhost:5001. See DEPLOYMENT.md for production deployment.
