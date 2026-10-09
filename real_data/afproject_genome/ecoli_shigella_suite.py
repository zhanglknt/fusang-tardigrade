#!/usr/bin/env python3
"""
Real indel-rich genome-scale validation, dataset 2: AFproject ecoli_shigella
(27 E. coli/Shigella whole genomes, reference tree Skippington & Ragan 2011).

Genome-scale hash-based k-mer cosine: k-mers hashed to uint64 (rolling 2-bit
encoding, reverse-complement canonical), counts via np.unique, pairwise cosine
via intersect1d. k scan {5,7,9,11} + multi-k average. Whole-genome MSA is not
feasible (4.6 Mb x 27) -- that is the alignment-free regime by design.
"""
import os, sys, json, time
import numpy as np

SRC = r"D:\系统发育树项目\Fusang\Fusang-main"
sys.path.insert(0, SRC)
os.chdir(SRC)

from fusang_v4_dahp_v1 import build_nj
from calc_nrf_simple import calc_nrf

DATA = os.path.join(SRC, "real_data", "afproject_genome", "ecoli_shigella")
REF_TREE = os.path.join(SRC, "real_data", "swisstree", "afproject_repo",
                        "datasets", "genome", "ecoli_shigella", "tree.newick")
OUT = os.path.join(SRC, "real_data", "afproject_genome", "ecoli_results")
os.makedirs(OUT, exist_ok=True)

COMP = bytes.maketrans(b"ACGTacgtNn", b"TGCAtgcaNn")
CANON = False  # forward-strand k-mers, consistent with af_competitor_methods baseline

def canon_kmer_hashes(seq: str, k: int) -> np.ndarray:
    """uint64 hashes of canonical (min of fwd/revcomp) k-mers."""
    s = seq.upper().encode()
    n = len(s) - k + 1
    if n <= 0:
        return np.zeros(0, dtype=np.uint64)
    # rolling 2-bit codes, N -> sentinel 4
    codes = np.full(len(s), 4, dtype=np.uint8)
    for b, c in ((65, 0), (67, 1), (71, 2), (84, 3)):  # A,C,G,T
        codes[np.frombuffer(s, dtype=np.uint8) == b] = c
    # forward hashes
    fwd = np.zeros(n, dtype=np.uint64)
    for i in range(k):
        fwd = (fwd << np.uint64(2)) | codes[i:i + n].astype(np.uint64)
    bad = np.zeros(n, dtype=bool)
    for i in range(k):
        bad |= codes[i:i + n] == 4
    # reverse-complement hashes: revcomp code = 3 - code, reversed order
    rc_codes = 3 - codes[::-1]
    rc = np.zeros(n, dtype=np.uint64)
    for i in range(k):
        rc = (rc << np.uint64(2)) | rc_codes[i:i + n].astype(np.uint64)
    if CANON:
        h = np.where(fwd <= rc, fwd, rc)
    else:
        h = fwd
    return h[~bad]

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

# ---- load ----
seqs, names = {}, []
for fn in sorted(os.listdir(DATA)):
    if fn.endswith(".fasta"):
        name = fn[:-6]
        parts, cur = [], None
        for line in open(os.path.join(DATA, fn)):
            if line.startswith(">"):
                cur = line[1:].strip()
            else:
                parts.append(line.strip())
        seqs[name] = "".join(parts)
        names.append(name)
print(f"loaded {len(names)} genomes, "
      f"{min(map(len, seqs.values()))/1e6:.1f}-{max(map(len, seqs.values()))/1e6:.1f} Mb")
seq_list = [seqs[n] for n in names]

# reference tree tips check
ref_tips = set(open(REF_TREE).read().strip().replace("(", " ").replace(")", " ")
               .replace(",", " ").replace(";", " ").split())
ref_tips = {t.split(":")[0] for t in ref_tips if t and not t[0].isdigit()}
assert not (ref_tips - set(names)), f"missing: {ref_tips - set(names)}"
assert not (set(names) - ref_tips), f"extra: {set(names) - ref_tips}"
print(f"tip sets match ({len(ref_tips)} tips)")

results = {}
def score(tag, tree_str, t0):
    path = os.path.join(OUT, f"{tag}.nwk")
    open(path, "w").write(tree_str.rstrip() + "\n")
    nrf = calc_nrf(path, REF_TREE)
    results[tag] = {"nrf": nrf, "time_s": round(time.time() - t0, 1)}
    print(f"  {tag}: nRF = {nrf:.4f}  ({results[tag]['time_s']}s)")

mats = {}
for k in (5, 7, 9, 11, 13):
    t0 = time.time()
    t1 = time.time()
    counts = [counts_from_hashes(canon_kmer_hashes(s, k)) for s in seq_list]
    t_hash = time.time() - t1
    D = cosine_dist(counts)
    print(f"  [k={k}] hash {t_hash:.0f}s + cosine {time.time()-t1-t_hash:.0f}s, "
          f"unique k-mers/genome: {counts[0][0].size:,}")
    mats[k] = D
    score(f"kmer_k{k}", build_nj(D, names), t0)

# multi-k variants
def multik_score(tag, ks):
    t0 = time.time()
    Dm = sum(mats[k] for k in ks) / len(ks)
    score(tag, build_nj(Dm, names), t0)

multik_score("multik_579", (5, 7, 9))            # manuscript definition
multik_score("multik_57911", (5, 7, 9, 11))
multik_score("multik_791113", (7, 9, 11, 13))    # genome-scale variant
Dm = sum(mats[k] for k in (7, 9, 11, 13)) / 4

# save phylip for AFproject-style fneighbor if needed
with open(os.path.join(OUT, "multik_phylip.txt"), "w") as f:
    f.write(f" {len(names)}\n")
    for i, n in enumerate(names):
        row = " ".join(f"{v:.6f}" for v in Dm[i])
        f.write(f"{n:<12s} {row}\n")

json.dump(results, open(os.path.join(OUT, "genome_cosine.json"), "w"), indent=1)
print("saved genome_cosine.json")
