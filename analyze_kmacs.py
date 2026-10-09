#!/usr/bin/env python3
"""Analyze kmacs distance matrices (phylip DMat) from the 27-seed AF benchmark.

For each seed (100-126) and k in {3,5,10}:
  - parse phylip distance matrix (from WSL kmacs run)
  - build NJ tree (BioPython, same as AF competitor protocol)
  - compute nRF vs FastTree2 tree (FT2-relative, Table 7 frame) and vs TRUE tree
Aggregates mean +/- SD, Wilcoxon paired vs Fusang (precomputed trees).

Usage: python analyze_kmacs.py  (expects kmacs_out/ copied next to this script)
"""
import os
import sys
import json
import tempfile
import numpy as np
from pathlib import Path

WORKDIR = Path(r"D:\系统发育树项目\Fusang\Fusang-main")
sys.path.insert(0, str(WORKDIR))
from calc_nrf_simple import calc_nrf  # noqa: E402
from Bio import Phylo  # noqa: E402
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor  # noqa: E402

SEEDS = list(range(100, 127))
KS = [3, 5, 10]
KMACS_DIR = WORKDIR / "kmacs_out"


def parse_phylip_dmat(path):
    """Parse full-square phylip matrix -> (names, full matrix rows).

    Note: kmacs writes a stray '\\r' after each taxon name; open with
    newline='' so universal newlines do not split rows in half.
    """
    with open(path, newline="") as f:
        lines = [ln for ln in f.read().split("\n") if ln.strip()]
    n = int(lines[0].split()[0])
    names, rows = [], []
    for ln in lines[1:]:
        parts = ln.split()
        names.append(parts[0])
        rows.append([float(x) for x in parts[1:]])
    assert len(names) == n and all(len(r) == n for r in rows), f"bad matrix {path}"
    # symmetrize (kmacs matrix should be symmetric)
    M = np.array(rows)
    M = (M + M.T) / 2.0
    np.fill_diagonal(M, 0.0)
    return names, M


def build_nj_newick(names, M):
    lower = [[M[i][j] for j in range(i + 1)] for i in range(len(names))]
    dm = DistanceMatrix(list(names), lower)
    tree = DistanceTreeConstructor().nj(dm)
    from io import StringIO
    buf = StringIO()
    Phylo.write(tree, buf, "newick", format_branch_length="%0.8f")
    return buf.getvalue().strip()


def nrf_between(nwk1, nwk2):
    with tempfile.NamedTemporaryFile("w", suffix=".nwk", delete=False) as f1:
        f1.write(nwk1)
        p1 = f1.name
    with tempfile.NamedTemporaryFile("w", suffix=".nwk", delete=False) as f2:
        f2.write(nwk2)
        p2 = f2.name
    try:
        return calc_nrf(p1, p2)
    finally:
        os.unlink(p1)
        os.unlink(p2)


def read_tree(path):
    with open(path, encoding="utf-8", errors="replace") as f:
        return f.read().strip()


def main():
    results = []
    for seed in SEEDS:
        ft2_nwk = read_tree(WORKDIR / f"seed{seed}_indel_ft2.nwk")
        true_nwk = read_tree(WORKDIR / f"seed{seed}_indel_true.nwk")
        fusang_nwk = read_tree(WORKDIR / f"seed{seed}_indel_fusang.nwk")
        row = {
            "seed": seed,
            "fusang_vs_ft2": nrf_between(fusang_nwk, ft2_nwk),
            "fusang_vs_true": nrf_between(fusang_nwk, true_nwk),
        }
        for k in KS:
            dm = KMACS_DIR / f"seed{seed}_kmacs_k{k}.dmat"
            if not dm.exists():
                row[f"kmacs_k{k}_vs_ft2"] = None
                row[f"kmacs_k{k}_vs_true"] = None
                continue
            names, M = parse_phylip_dmat(dm)
            nwk = build_nj_newick(names, M)
            row[f"kmacs_k{k}_vs_ft2"] = nrf_between(nwk, ft2_nwk)
            row[f"kmacs_k{k}_vs_true"] = nrf_between(nwk, true_nwk)
        results.append(row)
        print(f"seed {seed}: " + " ".join(
            f"k{k}={row.get(f'kmacs_k{k}_vs_ft2')}" for k in KS), flush=True)

    # Aggregate
    print("\n" + "=" * 70)
    print("kmacs results (27 seeds, n=200 indel, sub=0.05, indel=0.02)")
    print("=" * 70)
    summary = {}
    for k in KS:
        v_ft2 = [r[f"kmacs_k{k}_vs_ft2"] for r in results if r.get(f"kmacs_k{k}_vs_ft2") is not None]
        v_true = [r[f"kmacs_k{k}_vs_true"] for r in results if r.get(f"kmacs_k{k}_vs_true") is not None]
        summary[f"kmacs_k{k}"] = {
            "vs_ft2": {"mean": float(np.mean(v_ft2)), "std": float(np.std(v_ft2, ddof=1)), "n": len(v_ft2)} if v_ft2 else None,
            "vs_true": {"mean": float(np.mean(v_true)), "std": float(np.std(v_true, ddof=1)), "n": len(v_true)} if v_true else None,
        }
    f_ft2 = [r["fusang_vs_ft2"] for r in results]
    f_true = [r["fusang_vs_true"] for r in results]
    summary["fusang_k5gap2"] = {
        "vs_ft2": {"mean": float(np.mean(f_ft2)), "std": float(np.std(f_ft2, ddof=1)), "n": len(f_ft2)},
        "vs_true": {"mean": float(np.mean(f_true)), "std": float(np.std(f_true, ddof=1)), "n": len(f_true)},
    }

    for m, s in summary.items():
        ft2s = s["vs_ft2"]
        trs = s["vs_true"]
        print(f"\n  {m}:")
        print(f"    vs FT2 : {ft2s['mean']:.4f} ± {ft2s['std']:.4f} (n={ft2s['n']})")
        print(f"    vs TRUE: {trs['mean']:.4f} ± {trs['std']:.4f} (n={trs['n']})")

    # Wilcoxon paired vs Fusang (best k by vs_ft2 mean)
    from scipy.stats import wilcoxon
    best_k = min(KS, key=lambda k: np.mean([r[f"kmacs_k{k}_vs_ft2"] for r in results if r.get(f"kmacs_k{k}_vs_ft2") is not None]))
    pairs = [(r[f"kmacs_k{best_k}_vs_ft2"], r["fusang_vs_ft2"]) for r in results if r.get(f"kmacs_k{best_k}_vs_ft2") is not None]
    if len(pairs) >= 5:
        a, b = zip(*pairs)
        stat, p = wilcoxon(a, b)
        d = np.mean(np.array(a) - np.array(b)) / (np.std(np.array(a) - np.array(b), ddof=1) + 1e-12)
        print(f"\n  kmacs(best k={best_k}) vs Fusang (FT2-relative, paired):")
        print(f"    Wilcoxon p={p:.4g}, paired Cohen's d={d:+.2f}")

    out = WORKDIR / "kmacs_benchmark_results.json"
    with open(out, "w") as f:
        json.dump({"summary": summary, "results": results, "best_k": best_k}, f, indent=2)
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    main()
