#!/usr/bin/env python3
"""
P1: Generate missing benchmark data for reproducibility package.

NOTE: The current tree_simulation.py uses Kingman coalescent + JC69 model
(NO indel support). Indel-rich data requires INDELible or an upgraded simulator.
This script generates CLEAN (no-indel) benchmark data.

For each seed:
1. Generate coalescent tree + JC69 sequences
2. Build Fusang L0-1 NJ tree (k=5, gap2, cosine)
3. Compute nRF vs true tree
"""
import sys, os, csv, time, statistics as st, tempfile, argparse
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fusang.tree_simulation import make_coalescent_tree, simulate_seqs
from fusang.kmer_distance import compute_kmer_distance_matrix
from fusang.fusang_v4_dahp_v1 import build_nj
from fusang.calc_nrf_simple import calc_nrf

NT_MAP = {0: 'A', 1: 'T', 2: 'C', 3: 'G'}


def seq_array_to_string(arr):
    return ''.join(NT_MAP[b] for b in arr)


def generate_and_benchmark(seed, n=200, L=500, mu=0.05):
    """Generate simulated data and run Fusang L0-1 benchmark for one seed."""
    root, nodes = make_coalescent_tree(n, seed=seed)
    
    # Extract taxon names from coalescent tree leaves (t0001, t0002, ...)
    leaf_nodes = [nd for nd in nodes if nd.is_leaf]
    leaf_nodes.sort(key=lambda nd: nd.idx)
    taxon_names = [nd.name for nd in leaf_nodes]
    
    # simulate_seqs returns (n_taxa, L) int8 array
    seq_array = simulate_seqs(root, n, L, mu=mu, seed=seed + 1000)
    
    sequences = [seq_array_to_string(seq_array[i]) for i in range(n)]
    
    t0 = time.time()
    dist_matrix = compute_kmer_distance_matrix(
        sequences, taxon_names, k=5, metric='cosine',
        gap_pattern=(0, 1, 2, 5, 6),
        n_threads=1  # single-thread to avoid Windows ProcessPoolExecutor issues
    )
    nj_newick = build_nj(dist_matrix, taxon_names)
    fusang_time = time.time() - t0
    
    if nj_newick is None:
        return {'seed': seed, 'n': n, 'L': L, 'mu': mu,
                'fusang_nrf': None, 'fusang_time': round(fusang_time, 2), 'status': 'NJ_FAIL'}
    
    true_nwk = root.to_newick() + ';'
    tmp_true = tempfile.mktemp(suffix='.nwk')
    tmp_nj = tempfile.mktemp(suffix='.nwk')
    with open(tmp_true, 'w') as f: f.write(true_nwk)
    with open(tmp_nj, 'w') as f: f.write(nj_newick + ';\n')
    
    try:
        nrf_val = calc_nrf(tmp_true, tmp_nj)
    except Exception as e:
        print(f"  Seed {seed}: nRF error: {e}")
        nrf_val = None
    
    for p in [tmp_true, tmp_nj]:
        try: os.remove(p)
        except: pass
    
    return {'seed': seed, 'n': n, 'L': L, 'mu': mu,
            'fusang_nrf': nrf_val, 'fusang_time': round(fusang_time, 2),
            'status': 'OK' if nrf_val is not None else 'NRF_FAIL'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--start-seed', type=int, default=130)
    ap.add_argument('--end-seed', type=int, default=229)
    ap.add_argument('--n', type=int, default=200)
    ap.add_argument('--output', default='data/benchmark_clean_seeds130_229.csv')
    args = ap.parse_args()
    
    seeds = list(range(args.start_seed, args.end_seed + 1))
    all_results = []
    
    n_total = len(seeds)
    print(f"Benchmark: seeds {args.start_seed}-{args.end_seed}, n={args.n}, L=500, mu=0.05, JC69")
    print(f"NOTE: Clean (no-indel) data. Simulator uses Kingman coalescent + JC69.")
    t0_total = time.time()
    
    for i, seed in enumerate(seeds):
        t_start = time.time()
        result = generate_and_benchmark(seed, n=args.n)
        elapsed = time.time() - t_start
        all_results.append(result)
        
        nrf_str = f"{result['fusang_nrf']:.4f}" if result['fusang_nrf'] is not None else 'FAIL'
        print(f"  [{i+1}/{n_total}] seed={seed}: nRF={nrf_str}, {elapsed:.1f}s")
        
        if (i + 1) % 10 == 0:
            valid = [r for r in all_results if r['status'] == 'OK']
            if len(valid) >= 2:
                nrf_vals = [r['fusang_nrf'] for r in valid]
                print(f"    --- {len(valid)} valid, nRF={st.mean(nrf_vals):.4f}±{st.stdev(nrf_vals):.4f} ---")
            elif valid:
                print(f"    --- {len(valid)} valid, nRF={valid[0]['fusang_nrf']:.4f} ---")
    
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    fieldnames = ['seed', 'n', 'L', 'mu', 'fusang_nrf', 'fusang_time', 'status']
    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_results)
    
    valid = [r for r in all_results if r['status'] == 'OK']
    if valid:
        nrf_vals = [r['fusang_nrf'] for r in valid]
        if len(nrf_vals) >= 2:
            print(f"\nDone: {len(valid)}/{n_total} valid, nRF={st.mean(nrf_vals):.4f}±{st.stdev(nrf_vals):.4f}")
        else:
            print(f"\nDone: {len(valid)}/{n_total} valid, nRF={st.mean(nrf_vals):.4f}")
    print(f"Total: {time.time() - t0_total:.0f}s  Output: {args.output}")


if __name__ == '__main__':
    main()
