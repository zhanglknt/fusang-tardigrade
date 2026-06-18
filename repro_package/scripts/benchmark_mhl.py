#!/usr/bin/env python3
"""
Benchmark: Fusang MHL vs baseline methods (NJ, multi-k NJ, FastTree2).

Compares the MHL multi-level hierarchical pipeline against simpler
baseline methods across multiple simulated datasets.

Usage:
    python benchmark_mhl.py --seeds 100-129 --n 200
    python benchmark_mhl.py --seeds 100-129 --n 200 --no-l3

Output: benchmark_mhl_n{N}_{seeds}.csv
"""

import sys
import os
import argparse
import time
import csv
import json

sys.setrecursionlimit(10000)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from fusang_mhl.config import L0_K, L0_GAP, NRF_NORMALIZATION
from fusang_mhl.mlh_utils import read_fasta_simple
from calc_nrf_simple import get_bipartitions_from_newick


def compute_nrf(pred_nwk, ref_nwk):
    """Compute normalized RF distance between two Newick strings.

    Uses 2*(n-3) normalization (standard nRF).
    Returns NaN if either tree is invalid.
    """
    try:
        pred_leaves, pred_bip = get_bipartitions_from_newick(pred_nwk)
        ref_leaves, ref_bip = get_bipartitions_from_newick(ref_nwk)
        if not pred_bip or not ref_bip:
            return float('nan')

        # Restrict to common taxa
        common = pred_leaves & ref_leaves
        if len(common) < 4:
            return float('nan')

        # Filter bipartitions to common taxa
        pred_filtered = set()
        for bip in pred_bip:
            restricted = frozenset(bip & common)
            comp = frozenset(common - restricted)
            if 0 < len(restricted) < len(common):
                pred_filtered.add(frozenset(sorted(restricted)))

        ref_filtered = set()
        for bip in ref_bip:
            restricted = frozenset(bip & common)
            comp = frozenset(common - restricted)
            if 0 < len(restricted) < len(common):
                ref_filtered.add(frozenset(sorted(restricted)))

        if not pred_filtered or not ref_filtered:
            return float('nan')

        shared = pred_filtered & ref_filtered
        rf = len(pred_filtered | ref_filtered) - len(shared)
        n = len(common)
        max_rf = 2 * (n - 3)
        return rf / max_rf if max_rf > 0 else float('nan')
    except Exception as e:
        print(f"    [nRF] Error: {e}", file=sys.stderr)
        return float('nan')


def parse_seed_range(s):
    """Parse '100-129' or '100,101,102' into list of ints."""
    if '-' in s:
        a, b = s.split('-')
        return list(range(int(a), int(b) + 1))
    return [int(x) for x in s.split(',')]


def run_method_nj(sequences, taxon_names, D, verbose=False):
    """Baseline: Single NJ on full k-mer distance matrix."""
    from fusang_v4_dahp_v1 import build_nj
    t0 = time.perf_counter()
    nwk = build_nj(D, taxon_names)
    elapsed = time.perf_counter() - t0
    return nwk, elapsed


def run_method_multik_nj(sequences, taxon_names, verbose=False):
    """Baseline: Multi-k ensemble distance + NJ."""
    from fusang_v4_dahp_v1 import build_nj, compute_multik_distance_matrix
    t0 = time.perf_counter()
    D = compute_multik_distance_matrix(
        sequences, taxon_names,
        ks=(5, 7, 9), gap_pattern=None,
        n_threads=4, fusion_method="average",
    )
    D = np.array(D, dtype=np.float64)
    nwk = build_nj(D, taxon_names)
    elapsed = time.perf_counter() - t0
    return nwk, elapsed


def run_method_mhl(fasta_path, sequences, taxon_names, D,
                   no_l2=False, no_l3=False, verbose=False):
    """Fusang MHL full pipeline."""
    t0 = time.perf_counter()
    from fusang_mhl_main import run_fusang_mhl
    result = run_fusang_mhl(
        fasta_path=fasta_path,
        k=L0_K,
        gap=L0_GAP,
        n_threads=4,
        verbose=False,
        debug=False,
        no_l2=no_l2,
        no_l3=no_l3,
    )
    nwk = result.get("tree")
    elapsed = time.perf_counter() - t0
    return nwk, elapsed


def run_method_ft2(fasta_path, verbose=False):
    """Reference: MAFFT + FastTree2 (requires external tools)."""
    from fusang_mhl.level3_msa_ml import run_l3_on_cluster
    from fusang_mhl.mlh_utils import read_fasta_simple
    t0 = time.perf_counter()
    seqs = read_fasta_simple(fasta_path)
    names = list(seqs.keys())
    seqs_list = [seqs[n] for n in names]
    nwk = run_l3_on_cluster(seqs_list, names, verbose=False)
    elapsed = time.perf_counter() - t0
    return nwk, elapsed


def run_benchmark(
    seed_start, seed_end, data_dir=None,
    no_l2=False, no_l3=False,
    methods=("nj", "multik_nj", "mhl"),
    output_csv=None,
):
    """Run MHL benchmark across seeds."""
    if data_dir is None:
        data_dir = os.path.dirname(os.path.abspath(__file__))

    seeds = list(range(seed_start, seed_end + 1))
    total = len(seeds)
    results = []

    print(f"\n{'='*60}")
    print(f"  MHL Benchmark: {total} seeds")
    print(f"  Methods: {', '.join(methods)}")
    print(f"  L2: {'disabled' if no_l2 else 'enabled'}, "
          f"L3: {'disabled' if no_l3 else 'enabled'}")
    print(f"{'='*60}")

    for idx, seed in enumerate(seeds):
        fasta = os.path.join(data_dir, f'seed{seed}_indel.fasta')
        true_nwk = os.path.join(data_dir, f'seed{seed}_indel_true.nwk')

        if not os.path.exists(fasta):
            print(f'[{idx+1}/{total}] Seed {seed}: SKIP (no fasta)',
                  file=sys.stderr)
            continue

        if not os.path.exists(true_nwk):
            print(f'[{idx+1}/{total}] Seed {seed}: SKIP (no true tree)',
                  file=sys.stderr)
            continue

        # Read reference tree (handle encoding issues gracefully)
        try:
            with open(true_nwk, encoding="utf-8") as f:
                ref_nwk = f.read().strip()
        except UnicodeDecodeError:
            try:
                with open(true_nwk, encoding="latin-1") as f:
                    ref_nwk = f.read().strip()
            except Exception:
                ref_nwk = ""
        # Validate it looks like a Newick tree
        if not ref_nwk.startswith("("):
            ref_nwk = ""

        # Read sequences
        seqs_dict = read_fasta_simple(fasta)
        taxon_names = list(seqs_dict.keys())
        sequences = [seqs_dict[n] for n in taxon_names]
        n = len(taxon_names)

        if n < 4:
            continue

        print(f'\n[{idx+1}/{total}] Seed {seed} (n={n})', file=sys.stderr)

        # Compute k-mer distance (shared by NJ and MHL)
        from fusang_mhl.level0_kmer import compute_l0_distance
        t_kmer = time.perf_counter()
        D = compute_l0_distance(sequences, taxon_names,
                                 k=L0_K, gap=L0_GAP)
        t_kmer = time.perf_counter() - t_kmer
        print(f'  k-mer distance: {t_kmer:.2f}s', file=sys.stderr)

        row = {"seed": seed, "n": n}
        ref_ok = bool(ref_nwk)  # True if reference tree is valid

        # Method: NJ baseline
        if "nj" in methods:
            nwk, t = run_method_nj(sequences, taxon_names, D)
            nrf = compute_nrf(nwk, ref_nwk) if ref_ok else float('nan')
            row["nj_nrf"] = nrf
            row["nj_time"] = round(t, 2)
            print(f'  NJ: nRF={nrf:.4f} ({t:.2f}s)' if ref_ok
                  else f'  NJ: nRF=nan (bad ref) ({t:.2f}s)', file=sys.stderr)

        # Method: Multi-k NJ
        if "multik_nj" in methods:
            nwk, t = run_method_multik_nj(sequences, taxon_names)
            nrf = compute_nrf(nwk, ref_nwk) if ref_ok else float('nan')
            row["multik_nj_nrf"] = nrf
            row["multik_nj_time"] = round(t, 2)
            print(f'  Multi-k NJ: nRF={nrf:.4f} ({t:.2f}s)' if ref_ok
                  else f'  Multi-k NJ: nRF=nan (bad ref) ({t:.2f}s)', file=sys.stderr)

        # Method: MHL
        if "mhl" in methods:
            nwk, t = run_method_mhl(
                fasta, sequences, taxon_names, D,
                no_l2=no_l2, no_l3=no_l3,
            )
            if nwk:
                nrf = compute_nrf(nwk, ref_nwk) if ref_ok else float('nan')
                row["mhl_nrf"] = nrf
                row["mhl_time"] = round(t, 2)
                print(f'  MHL: nRF={nrf:.4f} ({t:.2f}s)' if ref_ok
                      else f'  MHL: nRF=nan (bad ref) ({t:.2f}s)', file=sys.stderr)
            else:
                row["mhl_nrf"] = float('nan')
                row["mhl_time"] = round(t, 2)
                print(f'  MHL: FAILED ({t:.2f}s)', file=sys.stderr)

        # Method: FastTree2 reference
        if "ft2" in methods:
            nwk, t = run_method_ft2(fasta)
            if nwk:
                nrf = compute_nrf(nwk, ref_nwk) if ref_ok else float('nan')
                row["ft2_nrf"] = nrf
                row["ft2_time"] = round(t, 2)
                print(f'  FT2: nRF={nrf:.4f} ({t:.2f}s)' if ref_ok
                      else f'  FT2: nRF=nan (bad ref) ({t:.2f}s)', file=sys.stderr)
            else:
                row["ft2_nrf"] = float('nan')
                row["ft2_time"] = round(t, 2)
                print(f'  FT2: FAILED ({t:.2f}s)', file=sys.stderr)

        results.append(row)

    # Save results
    if not results:
        print("\nNo results collected.", file=sys.stderr)
        return results

    # Determine columns
    base_cols = ["seed", "n"]
    method_cols = []
    for m in methods:
        if m == "nj":
            method_cols.extend(["nj_nrf", "nj_time"])
        elif m == "multik_nj":
            method_cols.extend(["multik_nj_nrf", "multik_nj_time"])
        elif m == "mhl":
            method_cols.extend(["mhl_nrf", "mhl_time"])
        elif m == "ft2":
            method_cols.extend(["ft2_nrf", "ft2_time"])

    all_cols = base_cols + method_cols

    if output_csv is None:
        output_csv = os.path.join(
            data_dir,
            f"benchmark_mhl_n{results[0]['n']}_s{seed_start}-{seed_end}.csv"
        )

    with open(output_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=all_cols, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(results)

    # Print summary
    print(f"\n{'='*60}")
    print(f"  Results saved to: {output_csv}")
    print(f"  Total seeds: {len(results)}")
    print(f"{'='*60}")

    for m in methods:
        col = f"{m}_nrf"
        vals = [r[col] for r in results if col in r and not (r[col] != r[col])]
        if vals:
            mean_val = np.nanmean(vals)
            std_val = np.nanstd(vals)
            print(f"  {m:15s}: nRF = {mean_val:.4f} +/- {std_val:.4f} "
                  f"(n={len(vals)})")

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark Fusang MHL vs baseline methods"
    )
    parser.add_argument("--seeds", type=str, default="100-129",
                        help="Seed range (e.g., 100-129)")
    parser.add_argument("--data-dir", type=str, default=None,
                        help="Directory with seed data")
    parser.add_argument("--no-l2", action="store_true",
                        help="Skip Level 2 (DAHP)")
    parser.add_argument("--no-l3", action="store_true",
                        help="Skip Level 3 (MSA+ML)")
    parser.add_argument("--methods", type=str, default="nj,multik_nj,mhl",
                        help="Methods to benchmark (comma-separated)")
    parser.add_argument("--output", type=str, default=None,
                        help="Output CSV file")
    args = parser.parse_args()

    seed_start, seed_end = parse_seed_range(args.seeds)[0], parse_seed_range(args.seeds)[-1]
    methods = [m.strip() for m in args.methods.split(",")]

    run_benchmark(
        seed_start=seed_start,
        seed_end=seed_end,
        data_dir=args.data_dir,
        no_l2=args.no_l2,
        no_l3=args.no_l3,
        methods=methods,
        output_csv=args.output,
    )


if __name__ == "__main__":
    main()
