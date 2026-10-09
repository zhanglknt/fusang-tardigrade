#!/usr/bin/env python3
"""Unified scoring for AFproject real-genome benchmarks using ete3 Tree.compare
(unrooted=True) -- identical to AFproject's ete_compare -- plus gap stats."""
import os, sys, json
import numpy as np
from ete3 import Tree

SRC = r"D:\系统发育树项目\Fusang\Fusang-main"
os.chdir(SRC)
sys.path.insert(0, SRC)

REFS = {
    "fish": os.path.join(SRC, "real_data", "swisstree", "afproject_repo",
                         "datasets", "genome", "fish_mito", "tree.newick"),
    "ecoli": os.path.join(SRC, "real_data", "swisstree", "afproject_repo",
                          "datasets", "genome", "ecoli_shigella", "tree.newick"),
}
WSL = os.path.join(SRC, "real_data", "afproject_genome", "wsl_out")
FISH_RES = os.path.join(SRC, "real_data", "afproject_genome", "fish_mito_results")
ECOLI_RES = os.path.join(SRC, "real_data", "afproject_genome", "ecoli_results")

def ete_nrf(tree_path, ref_path):
    qt = Tree(tree_path, format=1)
    rt = Tree(ref_path, format=1)
    res = qt.compare(rt, unrooted=True)
    return float(res["norm_rf"]), float(res["rf"]), float(res["max_rf"])

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
    def norm(n):
        n = n.replace("\\", "/").split("/")[-1]
        return n[:-6] if n.endswith(".fasta") else n
    D, names = {}, set()
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) < 3:
            continue
        a, b, d = norm(p[0]), norm(p[1]), float(p[2])
        if a == b:
            continue
        D[(a, b)] = d
        names.update((a, b))
    names = sorted(names)
    n = len(names)
    M = np.zeros((n, n))
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if i != j:
                M[i, j] = D.get((a, b), D.get((b, a), 1.0))
    return names, M

from fusang_v4_dahp_v1 import build_nj
from af_competitor_methods import cophylog_distance_matrix, read_fasta

out = {"fish": {}, "ecoli": {}}

def add(ds, tag, path):
    nrf, rf, maxrf = ete_nrf(path, REFS[ds])
    out[ds][tag] = {"nrf_ete": round(nrf, 4), "rf": rf, "max_rf": maxrf}
    print(f"  [{ds}] {tag}: nRF(ete3) = {nrf:.4f}")

# ---------- fish: python-method trees ----------
for tag in ["fusang_k5g2", "kmer_k5", "kmer_k7", "kmer_k9", "multik", "cophylog"]:
    p = os.path.join(FISH_RES, tag + ".nwk")
    if os.path.exists(p):
        add("fish", tag, p)

# fish: co-phylog halfctx=5 (AFproject published config, published nRF 0.09)
merged = []
DATA = os.path.join(SRC, "real_data", "afproject_genome", "fish_mito")
seqs, names = {}, []
for fn in sorted(os.listdir(DATA)):
    if fn.endswith(".fasta"):
        rec = read_fasta(os.path.join(DATA, fn))
        name = fn[:-6]
        seqs[name] = list(rec.values())[0]
        names.append(name)
D = cophylog_distance_matrix([seqs[n] for n in names], names, halfctx=5)
p5 = os.path.join(FISH_RES, "cophylog_h5.nwk")
open(p5, "w").write(build_nj(D, names).rstrip() + "\n")
add("fish", "cophylog_halfctx5_published_cfg", p5)

# fish: extended k scan k=11 + genome-scale multi-k (7,9,11)
from af_competitor_methods import kmer_cosine_distance_matrix
seq_list = [seqs[n] for n in names]
D11 = kmer_cosine_distance_matrix(seq_list, names, k=11, gap_pattern=None)
p11 = os.path.join(FISH_RES, "kmer_k11.nwk")
open(p11, "w").write(build_nj(D11, names).rstrip() + "\n")
add("fish", "kmer_k11", p11)
mats = {}
for k in (7, 9, 11):
    mats[k] = (kmer_cosine_distance_matrix(seq_list, names, k=k, gap_pattern=None)
               if k != 11 else D11)
Dm = sum(mats.values()) / 3.0
pm = os.path.join(FISH_RES, "multik_7911.nwk")
open(pm, "w").write(build_nj(Dm, names).rstrip() + "\n")
add("fish", "multik_7911", pm)

# fish: kmacs k=3 / k=10 (published k=10 -> 0.05)
for tag, fn in [("kmacs_k3", "fish_kmacs_k3.dmat"), ("kmacs_k10_published_cfg", "fish_kmacs_k10.dmat")]:
    names_k, M = parse_phylip(os.path.join(WSL, fn))
    p = os.path.join(WSL, "fish_" + tag + ".nwk")
    open(p, "w").write(build_nj(M, names_k).rstrip() + "\n")
    add("fish", tag, p)

# fish: mash k=21 s=10000 / k=11 s=5000 (published k=11 -> 0.05)
for tag, fn in [("mash_k21", "fish_mash_k21.tsv"), ("mash_k11_published_cfg", "fish_mash_k11.tsv")]:
    names_m, M = parse_mash_tsv(os.path.join(WSL, fn))
    p = os.path.join(WSL, "fish_" + tag + ".nwk")
    open(p, "w").write(build_nj(M, names_m).rstrip() + "\n")
    add("fish", tag, p)

# fish: MAFFT+FastTree2 (MSA+ML)
add("fish", "mafft_fasttree2_msa_ml", os.path.join(WSL, "fish_ft2.nwk"))

# fish gap stats (correct per-record parsing)
recs, cur = [], None
for line in open(os.path.join(WSL, "fish_aln.fasta")):
    if line.startswith(">"):
        if cur is not None:
            recs.append(cur)
        cur = ""
    else:
        cur += line.strip()
recs.append(cur)
L = len(recs[0])
gap = sum(s.count("-") for s in recs)
allgap = sum(1 for c in zip(*recs) if all(x == "-" for x in c))
out["fish"]["alignment"] = {"n": len(recs), "length": L,
                            "gap_fraction": round(gap / (L * len(recs)), 4),
                            "all_gap_cols": allgap}
print(f"  [fish] alignment: n={len(recs)} L={L} gap_fraction={gap/(L*len(recs)):.4f}")

# ---------- ecoli: genome cosine trees (if present) ----------
if os.path.isdir(ECOLI_RES):
    for tag in ["kmer_k5", "kmer_k7", "kmer_k9", "kmer_k11", "kmer_k13",
                "multik_579", "multik_57911", "multik_791113"]:
        p = os.path.join(ECOLI_RES, tag + ".nwk")
        if os.path.exists(p):
            add("ecoli", tag, p)

# ecoli: kmacs / mash if present
for tag, fn in [("kmacs_k3", "ecoli_kmacs_k3.dmat")]:
    fp = os.path.join(WSL, fn)
    if os.path.exists(fp):
        names_k, M = parse_phylip(fp)
        p = os.path.join(WSL, "ecoli_" + tag + ".nwk")
        open(p, "w").write(build_nj(M, names_k).rstrip() + "\n")
        add("ecoli", tag, p)
for tag, fn in [("mash_k21", "ecoli_mash_k21.tsv")]:
    fp = os.path.join(WSL, fn)
    if os.path.exists(fp):
        names_m, M = parse_mash_tsv(fp)
        p = os.path.join(WSL, "ecoli_" + tag + ".nwk")
        open(p, "w").write(build_nj(M, names_m).rstrip() + "\n")
        add("ecoli", tag, p)

json.dump(out, open(os.path.join(SRC, "real_data", "afproject_genome",
                                 "real_genome_scores.json"), "w"), indent=1)
print("saved real_genome_scores.json")
