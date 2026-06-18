#!/usr/bin/env python3
"""
Benchmark: Multi-k ensemble vs original Fusang on n=200 indel data.
Runs ensemble(k=5,7,9 contiguous) on multiple seeds, computes nRF vs FT2 reference.

Usage:
    python benchmark_multik_ensemble.py --seeds 230-259 --dir .
"""

import sys, os, argparse, time, csv
from pathlib import Path

sys.setrecursionlimit(10000)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fusang_v5_multik import (
    read_fasta, build_nj, compute_multik_distance_matrix,
    evaluate_against_reference, resolve_gap_pattern,
    compute_nrf
)
from calc_nrf_simple import get_bipartitions_from_newick
import numpy as np


def parse_seed_range(s):
    """Parse '230-259' or '230,231,232' into list of ints."""
    if '-' in s:
        a, b = s.split('-')
        return list(range(int(a), int(b) + 1))
    return [int(x) for x in s.split(',')]


def run_benchmark(seed_start, seed_end, data_dir, output_csv, threads=4):
    """Run multi-k ensemble benchmark on seeds [seed_start, seed_end]."""

    results = []
    total = seed_end - seed_start + 1

    for i, seed in enumerate(range(seed_start, seed_end + 1)):
        fasta = os.path.join(data_dir, f'seed{seed}_indel.fasta')
        ft2_nwk = os.path.join(data_dir, f'seed{seed}_indel_ft2.nwk')
        fusang_nwk = os.path.join(data_dir, f'seed{seed}_indel_fusang.nwk')

        # Check files exist
        if not os.path.exists(fasta):
            print(f'[{i+1}/{total}] Seed {seed}: SKIP (no fasta)')
            continue
        if not os.path.exists(ft2_nwk):
            print(f'[{i+1}/{total}] Seed {seed}: SKIP (no FT2 reference)')
            continue

        row = {'seed': seed}

        # Read original Fusang nRF (if exists)
        if os.path.exists(fusang_nwk):
            with open(fusang_nwk, 'r', encoding='utf-8', errors='replace') as f:
                fusang_tree = f.read().strip()
            row['nrf_fusang_original'] = evaluate_against_reference(fusang_tree, ft2_nwk)

        # Run multi-k ensemble (k=5,7,9 contiguous)
        print(f'[{i+1}/{total}] Seed {seed}: Running multi-k ensemble...')
        t0 = time.time()

        try:
            seqs = read_fasta(fasta)
            names = list(seqs.keys())
            sequences = [seqs[n] for n in names]

            fused_D, Ds = compute_multik_distance_matrix(
                sequences, names, ks=(5, 7, 9),
                gap_pattern=None, n_threads=threads,
                fusion_method='average'
            )

            nwk, tree = build_nj(fused_D, names)

            if nwk:
                nrf_ensemble = evaluate_against_reference(nwk, ft2_nwk)
                row['nrf_multik_ensemble'] = nrf_ensemble

                # Also compute individual k nRFs for comparison
                for k_val in [5, 7, 9]:
                    single_D, _ = compute_multik_distance_matrix(
                        sequences, names, ks=(k_val,),
                        gap_pattern=None, n_threads=threads,
                        fusion_method='average'
                    )
                    single_nwk, _ = build_nj(single_D, names)
                    if single_nwk:
                        nrf_single = evaluate_against_reference(single_nwk, ft2_nwk)
                        row[f'nrf_k{k_val}_contig'] = nrf_single

                elapsed = time.time() - t0
                row['time'] = round(elapsed, 1)
                print(f'[{i+1}/{total}] Seed {seed}: ensemble nRF={nrf_ensemble:.6f} ({elapsed:.1f}s)')
            else:
                row['nrf_multik_ensemble'] = float('inf')
                print(f'[{i+1}/{total}] Seed {seed}: NJ tree construction failed!')

        except Exception as e:
            row['nrf_multik_ensemble'] = float('inf')
            print(f'[{i+1}/{total}] Seed {seed}: ERROR: {e}')

        results.append(row)

        # Save intermediate results after each seed
        fieldnames = ['seed', 'nrf_fusang_original', 'nrf_k5_contig', 'nrf_k7_contig',
                       'nrf_k9_contig', 'nrf_multik_ensemble', 'time']
        with open(output_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(results)
        print(f'  -> Intermediate results saved to {output_csv}')

    # Print summary
    print(f'\n{"="*60}')
    print(f'Benchmark Summary ({len(results)} seeds)')
    print(f'{"="*60}')

    ensembles = [r['nrf_multik_ensemble'] for r in results if isinstance(r.get('nrf_multik_ensemble'), (int, float)) and r['nrf_multik_ensemble'] != float('inf')]
    originals = [r['nrf_fusang_original'] for r in results if isinstance(r.get('nrf_fusang_original'), (int, float)) and r['nrf_fusang_original'] != float('inf')]

    if ensembles:
        print(f'Multi-k Ensemble (k=5,7,9 contig):')
        print(f'  Mean nRF: {np.mean(ensembles):.6f} ± {np.std(ensembles):.6f}')
        print(f'  Min: {np.min(ensembles):.6f}, Max: {np.max(ensembles):.6f}')
        print(f'  Median: {np.median(ensembles):.6f}')

    if originals:
        print(f'Original Fusang (k=5 gap2):')
        print(f'  Mean nRF: {np.mean(originals):.6f} ± {np.std(originals):.6f}')
        print(f'  Min: {np.min(originals):.6f}, Max: {np.max(originals):.6f}')

    if ensembles and originals:
        # Paired comparison
        paired = [(r['nrf_fusang_original'], r['nrf_multik_ensemble']) for r in results
                  if isinstance(r.get('nrf_fusang_original'), (int, float)) and
                     isinstance(r.get('nrf_multik_ensemble'), (int, float)) and
                     r['nrf_fusang_original'] != float('inf') and
                     r['nrf_multik_ensemble'] != float('inf')]

        if paired:
            orig_nrf = [p[0] for p in paired]
            ens_nrf = [p[1] for p in paired]
            nrf_diff = [o - e for o, e in paired]
            ensemble_wins = sum(1 for o, e in paired if e < o)

            print(f'\nPaired Comparison ({len(paired)} seeds):')
            print(f'  Ensemble wins: {ensemble_wins}/{len(paired)}')
            print(f'  Mean nRF diff (orig - ensemble): {np.mean(nrf_diff):.6f}')
            print(f'  Improvement: {np.mean(nrf_diff)/np.mean(orig_nrf)*100:.1f}%')

            # Statistical test
            from scipy import stats
            t_stat, t_p = stats.ttest_rel(orig_nrf, ens_nrf)
            w_stat, w_p = stats.wilcoxon(orig_nrf, ens_nrf)
            cohens_d = np.mean(nrf_diff) / np.std(nrf_diff) if np.std(nrf_diff) > 0 else 0

            print(f'  Paired t-test: t={t_stat:.4f}, p={t_p:.6f}')
            print(f'  Wilcoxon: W={w_stat:.4f}, p={w_p:.6f}')
            print(f'  Cohen\'s d: {cohens_d:.4f}')

    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Multi-k ensemble benchmark')
    parser.add_argument('--seeds', default='230-259', help='Seed range (e.g., 230-259 or 230,231,232)')
    parser.add_argument('--dir', default='.', help='Data directory')
    parser.add_argument('--output', default='benchmark_multik_ensemble_n200_indel.csv', help='Output CSV')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads')
    args = parser.parse_args()

    seeds = parse_seed_range(args.seeds)
    run_benchmark(min(seeds), max(seeds), args.dir, args.output, args.threads)
