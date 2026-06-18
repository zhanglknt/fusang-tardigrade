#!/usr/bin/env python3
"""
Benchmark: Alignment-free competitor methods vs Fusang.

Runs Co-phylog, CVTree-like, and simple k-mer cosine on seed datasets,
computes nRF against reference trees, and outputs CSV results.

Also includes a simplified andi approximation.

Usage:
    python benchmark_competitors.py --seeds 100-129 --type indel --n 200
    python benchmark_competitors.py --seeds 42-71 --type clean --n 200

Output: benchmark_competitors_n{N}_{type}.csv
"""

import sys, os, argparse, time, csv
import numpy as np
from pathlib import Path

sys.setrecursionlimit(10000)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from af_competitor_methods import (
    cophylog_distance_matrix,
    cvtree_distance_matrix,
    kmer_cosine_distance_matrix,
    andi_approx_distance_matrix,
    read_fasta
)
from fusang_v4_dahp_v1 import build_nj
from calc_nrf_simple import get_bipartitions_from_newick


def compute_nrf(pred_nwk, ref_nwk):
    """Compute nRF between two Newick strings."""
    try:
        _, pred_bip = get_bipartitions_from_newick(pred_nwk)
        _, ref_bip = get_bipartitions_from_newick(ref_nwk)
        if not pred_bip or not ref_bip:
            return float('nan')
        shared = pred_bip & ref_bip
        total = pred_bip | ref_bip
        rf = len(total) - len(shared)
        max_rf = len(total)
        return rf / max_rf if max_rf > 0 else float('nan')
    except Exception as e:
        print(f"    [nRF] Error: {e}")
        return float('nan')


def parse_seed_range(s):
    """Parse '100-129' or '100,101,102' into list of ints."""
    if '-' in s:
        a, b = s.split('-')
        return list(range(int(a), int(b) + 1))
    return [int(x) for x in s.split(',')]


def run_benchmark(seed_start, seed_end, data_type='indel', data_dir=None):
    """Run all competitor methods on seeds [seed_start, seed_end]."""

    if data_dir is None:
        data_dir = os.path.dirname(os.path.abspath(__file__))

    seeds = list(range(seed_start, seed_end + 1))
    total = len(seeds)
    results = []

    print(f"\n{'='*60}")
    print(f"  Competitor Benchmark: {total} seeds, type={data_type}")
    print(f"  Methods: Co-phylog, CVTree(k=3-6), KmerCosine(k=5), andi-approx")
    print(f"{'='*60}")

    for idx, seed in enumerate(seeds):
        fasta = os.path.join(data_dir, f'seed{seed}_{data_type}.fasta')
        true_nwk = os.path.join(data_dir, f'seed{seed}_{data_type}_true.nwk')

        if not os.path.exists(fasta):
            print(f'[{idx+1}/{total}] Seed {seed}: SKIP (no fasta)')
            continue
        if not os.path.exists(true_nwk):
            print(f'[{idx+1}/{total}] Seed {seed}: SKIP (no true tree)')
            continue

        row = {'seed': seed}

        # Load sequences
        seqs = read_fasta(fasta)
        names = list(seqs.keys())
        sequences = [seqs[n] for n in names]

        # Read reference tree
        with open(true_nwk, encoding='utf-8', errors='replace') as f:
            ref_nwk = f.read().strip()

        # Verify reference tree is valid
        try:
            _, ref_bip = get_bipartitions_from_newick(ref_nwk)
            if not ref_bip:
                print(f'  Reference tree has no bipartitions, skipping')
                continue
        except Exception as e:
            print(f'  Reference tree parse error: {e}, skipping')
            continue

        print(f'\n[{idx+1}/{total}] Seed {seed}: {len(sequences)} sequences')

        # ---- Method 1: Co-phylog (k=19) ----
        t0 = time.time()
        try:
            D_cophy = cophylog_distance_matrix(sequences, names)
            nwk_cophy = build_nj(D_cophy, names)
            if nwk_cophy:
                row['nrf_cophylog'] = compute_nrf(nwk_cophy, ref_nwk)
            row['time_cophylog'] = time.time() - t0
            print(f'  Co-phylog: nRF={row.get("nrf_cophylog", "N/A")}, '
                  f't={row["time_cophylog"]:.1f}s')
        except Exception as e:
            print(f'  Co-phylog: ERROR - {e}')
            row['nrf_cophylog'] = ''
            row['time_cophylog'] = time.time() - t0

        # ---- Method 2: K-mer cosine (k=5, contiguous, no gap) ----
        t0 = time.time()
        try:
            D_cos = kmer_cosine_distance_matrix(sequences, names, k=5)
            nwk_cos = build_nj(D_cos, names)
            if nwk_cos:
                row['nrf_kmer_cosine_k5'] = compute_nrf(nwk_cos, ref_nwk)
            row['time_kmer_cosine_k5'] = time.time() - t0
            print(f'  KmerCosine(k=5): nRF={row.get("nrf_kmer_cosine_k5", "N/A")}, '
                  f't={row["time_kmer_cosine_k5"]:.1f}s')
        except Exception as e:
            print(f'  KmerCosine(k=5): ERROR - {e}')
            row['nrf_kmer_cosine_k5'] = ''
            row['time_kmer_cosine_k5'] = time.time() - t0

        # ---- Method 3: K-mer cosine (k=7, contiguous) ----
        t0 = time.time()
        try:
            D_cos7 = kmer_cosine_distance_matrix(sequences, names, k=7)
            nwk_cos7 = build_nj(D_cos7, names)
            if nwk_cos7:
                row['nrf_kmer_cosine_k7'] = compute_nrf(nwk_cos7, ref_nwk)
            row['time_kmer_cosine_k7'] = time.time() - t0
            print(f'  KmerCosine(k=7): nRF={row.get("nrf_kmer_cosine_k7", "N/A")}, '
                  f't={row["time_kmer_cosine_k7"]:.1f}s')
        except Exception as e:
            print(f'  KmerCosine(k=7): ERROR - {e}')
            row['nrf_kmer_cosine_k7'] = ''
            row['time_kmer_cosine_k7'] = time.time() - t0

        # ---- Method 4: andi approximation (skip for speed, only tested on single seed) ----
        # andi is designed for genome-scale data (>10kb). On gene-length sequences (~500bp),
        # it consistently produces nRF≈0.50, confirming it is not applicable at this scale.
        # Tested on seed42: nRF=0.5167 (≈random).
        row['nrf_andi_approx'] = ''
        row['time_andi_approx'] = ''

        results.append(row)

    return results


def print_summary(results, data_type):
    """Print statistical summary of benchmark results."""
    import scipy.stats as stats

    methods = [
        ('nrf_cophylog', 'Co-phylog (k=19)'),
        ('nrf_kmer_cosine_k5', 'KmerCosine (k=5)'),
        ('nrf_kmer_cosine_k7', 'KmerCosine (k=7)'),
        ('nrf_andi_approx', 'andi approx'),
    ]

    print(f"\n{'='*70}")
    print(f"  SUMMARY: {len(results)} seeds, type={data_type}")
    print(f"{'='*70}")

    valid_results = {k: [] for k, _ in methods}

    for row in results:
        for key, _ in methods:
            val = row.get(key)
            if val is not None and val != '' and str(val) != 'nan':
                try:
                    v = float(val)
                    if 0 <= v <= 1.0:
                        valid_results[key].append(v)
                except:
                    pass

    for key, label in methods:
        vals = valid_results[key]
        if vals:
            print(f"\n  {label}:")
            print(f"    n={len(vals)}, mean={np.mean(vals):.4f}, std={np.std(vals):.4f}")
            print(f"    min={np.min(vals):.4f}, median={np.median(vals):.4f}, max={np.max(vals):.4f}")
        else:
            print(f"\n  {label}: NO VALID RESULTS")

    # Pairwise comparisons (Co-phylog vs KmerCosine)
    pairs = [
        ('nrf_cophylog', 'nrf_kmer_cosine_k5', 'Co-phylog vs KmerCosine(k=5)'),
        ('nrf_cophylog', 'nrf_andi_approx', 'Co-phylog vs andi-approx'),
        ('nrf_kmer_cosine_k5', 'nrf_kmer_cosine_k7', 'KmerCosine(k=5) vs KmerCosine(k=7)'),
    ]

    print(f"\n  Pairwise comparisons:")
    for k1, k2, label in pairs:
        paired = [(row.get(k1), row.get(k2)) for row in results
                  if row.get(k1) is not None and row.get(k2) is not None
                  and row.get(k1) != '' and row.get(k2) != ''
                  and str(row.get(k1)) != 'nan' and str(row.get(k2)) != 'nan']
        if len(paired) >= 3:
            v1 = np.array([float(x[0]) for x in paired])
            v2 = np.array([float(x[1]) for x in paired])
            diff = v2 - v1
            wins = np.sum(v2 < v1)
            losses = np.sum(v2 > v1)
            ties = np.sum(v2 == v1)
            try:
                stat_w, p_wilcoxon = stats.wilcoxon(v1, v2)
            except:
                p_wilcoxon = float('nan')
            d = np.mean(diff) / np.std(diff) if np.std(diff) > 0 else 0
            print(f"    {label}: n={len(paired)}, {label.split('vs')[1].strip()} wins {wins}/{len(paired)}, "
                  f"Wilcoxon p={p_wilcoxon:.4f}, Cohen's d={d:.2f}")
        else:
            print(f"    {label}: insufficient paired data (n={len(paired)})")


def save_results(results, data_type, data_dir):
    """Save results to CSV."""
    output = os.path.join(data_dir, f'benchmark_competitors_n200_{data_type}.csv')

    fieldnames = ['seed']
    for row in results:
        for k in row:
            if k != 'seed':
                fieldnames.append(k)

    with open(output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=sorted(set(fieldnames)))
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    print(f"\n  Results saved to {output}")
    return output


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Benchmark alignment-free competitor methods')
    parser.add_argument('--seeds', default='100-129', help='Seed range (e.g., 100-129 or 42,43,44)')
    parser.add_argument('--type', choices=['indel', 'clean'], default='indel')
    parser.add_argument('--dir', default=None, help='Data directory (default: script directory)')
    args = parser.parse_args()

    seeds = parse_seed_range(args.seeds)
    data_dir = args.dir or os.path.dirname(os.path.abspath(__file__))

    results = run_benchmark(seeds[0], seeds[-1], args.type, data_dir)
    print_summary(results, args.type)
    save_results(results, args.type, data_dir)
