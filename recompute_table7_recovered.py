#!/usr/bin/env python3
"""Full Table 7 recomputation on recovered benchmark data.
All methods, 27 seeds (100-126), NJ trees, nRF vs FastTree2 reference trees.
"""
import json
import os
import sys
import tempfile
import numpy as np
from pathlib import Path
from io import StringIO

WORKDIR = Path(r"D:/系统发育树项目/Fusang/Fusang-main")
sys.path.insert(0, str(WORKDIR))
from af_competitor_methods import kmer_cosine_distance_matrix, cophylog_distance_matrix, read_fasta  # noqa: E402
from fusang_v4_dahp_v1 import build_nj  # noqa: E402
from calc_nrf_simple import calc_nrf  # noqa: E402
from scipy.stats import wilcoxon  # noqa: E402
from Bio import Phylo  # noqa: E402
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor  # noqa: E402


def nrf_pair(nwk1, nwk2):
    with tempfile.NamedTemporaryFile("w", suffix=".nwk", delete=False) as f1, \
         tempfile.NamedTemporaryFile("w", suffix=".nwk", delete=False) as f2:
        f1.write(nwk1)
        f2.write(nwk2)
        p1, p2 = f1.name, f2.name
    try:
        return calc_nrf(p1, p2)
    finally:
        os.unlink(p1)
        os.unlink(p2)


def parse_phylip_dmat(path):
    with open(path, newline="") as f:
        lines = [ln for ln in f.read().split("\n") if ln.strip()]
    n = int(lines[0].split()[0])
    names, rows = [], []
    for ln in lines[1:]:
        parts = ln.split()
        names.append(parts[0])
        rows.append([float(x) for x in parts[1:]])
    M = np.array(rows)
    M = (M + M.T) / 2.0
    np.fill_diagonal(M, 0.0)
    return names, M


def nj_newick(names, M):
    lower = [[M[i][j] for j in range(i + 1)] for i in range(len(names))]
    dm = DistanceMatrix(list(names), lower)
    tree = DistanceTreeConstructor().nj(dm)
    buf = StringIO()
    Phylo.write(tree, buf, "newick")
    return buf.getvalue().strip()


def main():
    results = []
    for seed in range(100, 127):
        seqs = read_fasta(str(WORKDIR / f"seed{seed}_indel.fasta"))
        names = list(seqs.keys())
        sequences = [seqs[n] for n in names]
        ft2_nwk = (WORKDIR / f"seed{seed}_indel_ft2.nwk").read_text(encoding="latin-1").strip()
        row = {"seed": seed}

        D = cophylog_distance_matrix(sequences, names, halfctx=9)
        row["cophylog"] = nrf_pair(ft2_nwk, build_nj(D, names))

        kn, M = parse_phylip_dmat(WORKDIR / "kmacs_out" / f"seed{seed}_kmacs_k3.dmat")
        row["kmacs_k3"] = nrf_pair(ft2_nwk, nj_newick(kn, M))

        for k in (5, 7):
            D = kmer_cosine_distance_matrix(sequences, names, k=k, gap_pattern=None)
            row[f"kmer{k}"] = nrf_pair(ft2_nwk, build_nj(D, names))

        fus_nwk = (WORKDIR / f"seed{seed}_indel_fusang.nwk").read_text(encoding="latin-1").strip()
        row["fusang"] = nrf_pair(ft2_nwk, fus_nwk)

        mats = [np.array(kmer_cosine_distance_matrix(sequences, names, k=k, gap_pattern=None))
                for k in (5, 7, 9)]
        avg = np.mean(mats, axis=0)
        row["multik"] = nrf_pair(ft2_nwk, build_nj(avg, names))

        results.append(row)
        print(f"seed {seed}: cophy={row['cophylog']:.3f} kmacs={row['kmacs_k3']:.3f} "
              f"k5={row['kmer5']:.3f} k7={row['kmer7']:.3f} "
              f"fusang={row['fusang']:.3f} multik={row['multik']:.3f}", flush=True)

    print("\n" + "=" * 70)
    print("TABLE 7 (recomputed on recovered data, 27 seeds, FT2-relative)")
    print("=" * 70)
    methods = ["cophylog", "kmacs_k3", "kmer5", "kmer7", "fusang", "multik"]
    summary = {}
    for m in methods:
        v = [r[m] for r in results if r.get(m) is not None]
        summary[m] = {"mean": float(np.mean(v)), "std": float(np.std(v, ddof=1)), "n": len(v)}
        print(f"  {m:10s}: {np.mean(v):.4f} +/- {np.std(v, ddof=1):.4f} (n={len(v)})")

    print("\nPaired vs Fusang:")
    for m in methods:
        if m == "fusang":
            continue
        a = np.array([r[m] for r in results])
        b = np.array([r["fusang"] for r in results])
        stat, p = wilcoxon(a, b)
        d = np.mean(a - b) / (np.std(a - b, ddof=1) + 1e-12)
        print(f"  {m:10s}: Wilcoxon p={p:.4g}, paired d={d:+.2f}")

    out = WORKDIR / "table7_recomputed_recovered.json"
    with open(out, "w") as f:
        json.dump({"summary": summary, "results": results}, f, indent=2)
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    main()
