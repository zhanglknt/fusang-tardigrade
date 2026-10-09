#!/usr/bin/env python3
"""Score external-tool outputs (kmacs dmat, mash tsv, FastTree2 newick)
against AFproject reference trees."""
import os, sys, json
import numpy as np

SRC = r"D:\系统发育树项目\Fusang\Fusang-main"
sys.path.insert(0, SRC)
os.chdir(SRC)

from fusang_v4_dahp_v1 import build_nj
from calc_nrf_simple import calc_nrf

WSL = os.path.join(SRC, "real_data", "afproject_genome", "wsl_out")
REFS = {
    "fish": os.path.join(SRC, "real_data", "swisstree", "afproject_repo",
                         "datasets", "genome", "fish_mito", "tree.newick"),
    "ecoli": os.path.join(SRC, "real_data", "swisstree", "afproject_repo",
                          "datasets", "genome", "ecoli_shigella", "tree.newick"),
}
results = {}

def parse_phylip(path):
    with open(path, newline="") as f:
        lines = [l.rstrip("\r\n") for l in f if l.strip()]
    n = int(lines[0].split()[0])
    names, M = [], np.zeros((n, n))
    for i, line in enumerate(lines[1:1 + n]):
        parts = line.split()
        names.append(parts[0])
        M[i] = [float(x) for x in parts[1:1 + n]]
    return names, M

def parse_mash_tsv(path):
    D, names = {}, set()
    for line in open(path):
        p = line.split("\t")
        a, b, d = p[0], p[1], float(p[2])
        D[(a, b)] = d
        names.update((a, b))
    names = sorted(names)
    n = len(names)
    M = np.zeros((n, n))
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if i == j:
                continue
            M[i, j] = D.get((a, b), D.get((b, a), 1.0))
    return names, M

def score(tag, tree_str_or_path, ref, is_path=False):
    path = tree_str_or_path if is_path else os.path.join(WSL, tag + ".nwk")
    if not is_path:
        open(path, "w").write(tree_str_or_path.rstrip() + "\n")
    nrf = calc_nrf(path, ref)
    results[tag] = nrf
    print(f"  {tag}: nRF = {nrf:.4f}")

# ---- fish ----
names, M = parse_phylip(os.path.join(WSL, "fish_kmacs_k3.dmat"))
score("fish_kmacs_k3", build_nj(M, names), REFS["fish"])
names, M = parse_mash_tsv(os.path.join(WSL, "fish_mash.tsv"))
score("fish_mash_k21", build_nj(M, names), REFS["fish"])
score("fish_ft2_msa_ml", os.path.join(WSL, "fish_ft2.nwk"), REFS["fish"], is_path=True)

# ---- ecoli ----
names, M = parse_phylip(os.path.join(WSL, "ecoli_kmacs_k3.dmat"))
score("ecoli_kmacs_k3", build_nj(M, names), REFS["ecoli"])
names, M = parse_mash_tsv(os.path.join(WSL, "ecoli_mash.tsv"))
score("ecoli_mash_k21", build_nj(M, names), REFS["ecoli"])

json.dump(results, open(os.path.join(WSL, "external_scores.json"), "w"), indent=1)
print("saved external_scores.json")
