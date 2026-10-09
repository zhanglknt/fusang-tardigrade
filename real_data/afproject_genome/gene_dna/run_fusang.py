#!/usr/bin/env python3
"""Step 3a: per-gene Fusang k-mer cosine distance + BioPython NJ.
k = 5, 7, 9 + multi-k(5,7,9) average. Forward-strand hashing (CANON=False),
identical semantics to af_competitor_methods.kmer_cosine_distance_matrix."""
import os, sys, time
import numpy as np

SRC = r"D:\系统发育树项目\Fusang\Fusang-main"
sys.path.insert(0, SRC)
os.chdir(SRC)
from fusang_v4_dahp_v1 import build_nj

BASE = os.path.join(SRC, "real_data", "afproject_genome", "gene_dna")
FASTA = os.path.join(BASE, "fasta")
TREES = os.path.join(BASE, "trees")
os.makedirs(TREES, exist_ok=True)

GENES = ["ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
         "ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB"]

def canon_kmer_hashes(seq, k):
    s = seq.upper().encode()
    n = len(s) - k + 1
    if n <= 0:
        return np.zeros(0, dtype=np.uint64)
    codes = np.full(len(s), 4, dtype=np.uint8)
    for b, c in ((65, 0), (67, 1), (71, 2), (84, 3)):
        codes[np.frombuffer(s, dtype=np.uint8) == b] = c
    fwd = np.zeros(n, dtype=np.uint64)
    for i in range(k):
        fwd = (fwd << np.uint64(2)) | codes[i:i + n].astype(np.uint64)
    bad = np.zeros(n, dtype=bool)
    for i in range(k):
        bad |= codes[i:i + n] == 4
    return fwd[~bad]

def counts_from_hashes(h):
    u, c = np.unique(h, return_counts=True)
    return u.astype(np.uint64), c.astype(np.float64)

def cosine_dist(counts_list):
    n = len(counts_list)
    D = np.zeros((n, n))
    for i in range(n):
        ui, ci = counts_list[i]
        for j in range(i + 1, n):
            uj, cj = counts_list[j]
            _, ii, jj = np.intersect1d(ui, uj, assume_unique=True,
                                       return_indices=True)
            dot = float(np.dot(ci[ii], cj[jj]))
            na = float(np.dot(ci, ci)) ** 0.5
            nb = float(np.dot(cj, cj)) ** 0.5
            D[i, j] = D[j, i] = 1.0 - dot / (na * nb + 1e-30)
    return D

def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:]
            seqs[name] = ""
        elif name:
            seqs[name] += line
    return seqs

for g in GENES:
    seqs = read_fasta(os.path.join(FASTA, g + ".fasta"))
    names = sorted(seqs)
    seq_list = [seqs[n] for n in names]
    mats = {}
    for k in (5, 7, 9):
        t0 = time.time()
        counts = [counts_from_hashes(canon_kmer_hashes(s, k)) for s in seq_list]
        D = cosine_dist(counts)
        mats[k] = D
        tree = build_nj(D, names)
        open(os.path.join(TREES, f"{g}_fusang_k{k}.nwk"), "w").write(tree.rstrip() + "\n")
    Dm = sum(mats.values()) / 3.0
    tree = build_nj(Dm, names)
    open(os.path.join(TREES, f"{g}_fusang_multik.nwk"), "w").write(tree.rstrip() + "\n")
    print(f"{g:6s} done (L~{len(seq_list[0])})")
print("all fusang trees written to", TREES)
