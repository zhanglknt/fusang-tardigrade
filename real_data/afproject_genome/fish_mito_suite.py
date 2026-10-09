#!/usr/bin/env python3
"""
Real indel-rich genome-scale validation, dataset 1: AFproject fish_mito
(25 fish mitochondrial genomes, reference tree Fischer et al. 2013).

Python-side methods: Fusang spaced k=5,gap2 / contiguous k=5,7,9 / multi-k /
Co-phylog (halfctx=9). External tools (kmacs, Mash, MAFFT+FastTree2) run
via WSL separately (run_real_genome_wsl.sh). Scoring vs reference tree via
calc_nrf (nRF = RF / 2(n-3), identical definition to AFproject ete_compare).
"""
import os, sys, json, time
import numpy as np

SRC = r"D:\系统发育树项目\Fusang\Fusang-main"
sys.path.insert(0, SRC)
os.chdir(SRC)

from af_competitor_methods import kmer_cosine_distance_matrix, cophylog_distance_matrix, read_fasta
from fusang_v4_dahp_v1 import build_nj
from calc_nrf_simple import calc_nrf

DATA = os.path.join(SRC, "real_data", "afproject_genome", "fish_mito")
REF_TREE = os.path.join(SRC, "real_data", "swisstree", "afproject_repo",
                        "datasets", "genome", "fish_mito", "tree.newick")
OUT = os.path.join(SRC, "real_data", "afproject_genome", "fish_mito_results")
os.makedirs(OUT, exist_ok=True)

# ---- load sequences (one genome per file, tip name = file stem) ----
seqs, names = {}, []
for fn in sorted(os.listdir(DATA)):
    if fn.endswith(".fasta"):
        rec = read_fasta(os.path.join(DATA, fn))
        assert len(rec) == 1
        name = fn[:-6]
        seqs[name] = list(rec.values())[0]
        names.append(name)
print(f"loaded {len(names)} genomes, lengths "
      f"{min(map(len, seqs.values()))}-{max(map(len, seqs.values()))} bp")

# reference tree tips
ref_tips = set(open(REF_TREE).read().strip().replace("(", " ").replace(")", " ")
               .replace(",", " ").replace(";", " ").split())
ref_tips = {t.split(":")[0] for t in ref_tips if t and not t[0].isdigit()}
missing = ref_tips - set(names)
extra = set(names) - ref_tips
assert not missing, f"missing tips: {missing}"
assert not extra, f"extra tips: {extra}"
print(f"tip sets match reference tree ({len(ref_tips)} tips)")

seq_list = [seqs[n] for n in names]
results = {}

def score(tag, tree_str, t0):
    path = os.path.join(OUT, f"{tag}.nwk")
    open(path, "w").write(tree_str.rstrip() + "\n")
    nrf = calc_nrf(path, REF_TREE)
    results[tag] = {"nrf": nrf, "time_s": round(time.time() - t0, 1)}
    print(f"  {tag}: nRF = {nrf:.4f}  ({results[tag]['time_s']}s)")

# ---- Fusang spaced k=5,gap2 (default) ----
gap_pattern = "1001001001001"  # k=5, g=2
t0 = time.time()
D = kmer_cosine_distance_matrix(seq_list, names, k=5, gap_pattern=gap_pattern)
score("fusang_k5g2", build_nj(D, names), t0)

# ---- contiguous k = 5, 7, 9 ----
mats = {}
for k in (5, 7, 9):
    t0 = time.time()
    D = kmer_cosine_distance_matrix(seq_list, names, k=k, gap_pattern=None)
    mats[k] = D
    score(f"kmer_k{k}", build_nj(D, names), t0)

# ---- multi-k ensemble (average of k=5,7,9) ----
t0 = time.time()
Dm = (mats[5] + mats[7] + mats[9]) / 3.0
score("multik", build_nj(Dm, names), t0)

# ---- Co-phylog (halfctx=9, k=19; DNA benchmark config) ----
t0 = time.time()
D = cophylog_distance_matrix(seq_list, names, halfctx=9)
score("cophylog", build_nj(D, names), t0)

# also save multi-k phylip distance matrix for reference
with open(os.path.join(OUT, "multik_phylip.txt"), "w") as f:
    f.write(f" {len(names)}\n")
    for i, n in enumerate(names):
        row = " ".join(f"{v:.6f}" for v in Dm[i])
        f.write(f"{n:<12s} {row}\n")

json.dump(results, open(os.path.join(OUT, "python_methods.json"), "w"), indent=1)
print("saved python_methods.json")
