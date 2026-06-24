#!/usr/bin/env python3
"""
D4: Generate multi-k ensemble benchmark data.

For each seed:
1. Build NJ trees with k=5 (gap2), k=7 (gap2), k=9 (contiguous)
2. Build ensemble via average(k=5, k=7, k=9) distance matrix + NJ
3. Compute nRF vs true tree for all methods

NOTE: Uses JC69 simulator (tree_simulation.py), NOT INDELible.
"""
import sys, os, csv, time, statistics as st, tempfile, argparse, json
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fusang.tree_simulation import make_coalescent_tree, simulate_seqs
from fusang.kmer_distance import compute_kmer_distance_matrix
from fusang.fusang_v4_dahp_v1 import build_nj
from fusang.calc_nrf_simple import calc_nrf

NT_MAP = {0: 'A', 1: 'T', 2: 'C', 3: 'G'}


def seq_array_to_string(arr):
    return ''.join(NT_MAP[b] for b in arr)


def benchmark_one(seed, n=200, L=500, mu=0.05):
    """Multi-k benchmark for one seed."""
    root, nodes = make_coalescent_tree(n, seed=seed)
    
    leaf_nodes = [nd for nd in nodes if nd.is_leaf]
    leaf_nodes.sort(key=lambda nd: nd.idx)
    taxon_names = [nd.name for nd in leaf_nodes]
    
    seq_array = simulate_seqs(root, n, L, mu=mu, seed=seed + 1000)
    sequences = [seq_array_to_string(seq_array[i]) for i in range(n)]
    
    # Save true tree
    true_nwk = root.to_newick() + ';'
    
    results = {'seed': seed, 'n': n, 'L': L, 'mu': mu}
    nwk_trees = {}
    
    # k=5 gap2
    t0 = time.time()
    D5 = compute_kmer_distance_matrix(sequences, taxon_names, k=5, metric='cosine',
                                       gap_pattern=(0, 1, 2, 5, 6), n_threads=1)
    nwk5 = build_nj(D5, taxon_names)
    t5 = time.time() - t0
    results['k5_time_s'] = round(t5, 2)
    if nwk5:
        nwk_trees['k5'] = nwk5
    else:
        results['k5_nrf'] = None
    
    # k=7 gap2
    gap7 = (0, 1, 2, 7, 8)  # k=7 with 2 gaps
    t0 = time.time()
    D7 = compute_kmer_distance_matrix(sequences, taxon_names, k=7, metric='cosine',
                                       gap_pattern=gap7, n_threads=1)
    nwk7 = build_nj(D7, taxon_names)
    t7 = time.time() - t0
    results['k7_time_s'] = round(t7, 2)
    if nwk7:
        nwk_trees['k7'] = nwk7
    else:
        results['k7_nrf'] = None
    
    # k=9 contiguous
    t0 = time.time()
    D9 = compute_kmer_distance_matrix(sequences, taxon_names, k=9, metric='cosine',
                                       n_threads=1)
    nwk9 = build_nj(D9, taxon_names)
    t9 = time.time() - t0
    results['k9_time_s'] = round(t9, 2)
    if nwk9:
        nwk_trees['k9'] = nwk9
    else:
        results['k9_nrf'] = None
    
    # Multi-k ensemble: average of k=5(gap2) + k=7(gap2) + k=9(contiguous)
    # Note: compute_multik_distance_matrix applies same gap_pattern to all ks
    # Here we manually average D5, D7, D9 with their respective gap patterns
    t0 = time.time()
    if D5 is not None and D7 is not None and D9 is not None:
        D_ens = (D5 + D7 + D9) / 3.0
    else:
        D_ens = None
    nwk_ens = build_nj(D_ens, taxon_names) if D_ens is not None else None
    t_ens = time.time() - t0
    results['ensemble_time_s'] = round(t_ens, 2)
    if nwk_ens:
        nwk_trees['ensemble'] = nwk_ens
    else:
        results['ensemble_nrf'] = None
    
    # Compute nRF for all
    for key, nwk in nwk_trees.items():
        tmp_true = tempfile.mktemp(suffix='.nwk')
        tmp_nj = tempfile.mktemp(suffix='.nwk')
        with open(tmp_true, 'w') as f: f.write(true_nwk)
        with open(tmp_nj, 'w') as f: f.write(nwk + ';\n')
        try:
            nrf_val = calc_nrf(tmp_true, tmp_nj)
            results[f'{key}_nrf'] = round(nrf_val, 6)
        except Exception as e:
            print(f"  {key} nRF error: {e}")
            results[f'{key}_nrf'] = None
        for p in [tmp_true, tmp_nj]:
            try: os.remove(p)
            except: pass
    
    results['status'] = 'OK'
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--start-seed', type=int, default=200)
    ap.add_argument('--end-seed', type=int, default=229)
    ap.add_argument('--n', type=int, default=200)
    ap.add_argument('--output', default='data/benchmark_multik_seeds200_229.csv')
    args = ap.parse_args()
    
    seeds = list(range(args.start_seed, args.end_seed + 1))
    all_results = []
    n_total = len(seeds)
    
    print(f"Multi-k benchmark: seeds {args.start_seed}-{args.end_seed}, n={args.n}, JC69")
    print(f"Methods: k=5(gap2), k=7(gap2), k=9(contiguous), ensemble(avg)")
    t0_total = time.time()
    
    for i, seed in enumerate(seeds):
        t0 = time.time()
        result = benchmark_one(seed, n=args.n)
        elapsed = time.time() - t0
        all_results.append(result)
        
        nrf_strs = []
        for key in ['k5', 'k7', 'k9', 'ensemble']:
            v = result.get(f'{key}_nrf')
            nrf_strs.append(f"{key}={v:.4f}" if v is not None else f"{key}=FAIL")
        print(f"  [{i+1}/{n_total}] seed={seed}: {', '.join(nrf_strs)}, {elapsed:.1f}s")
        
        if (i + 1) % 5 == 0:
            valid = [r for r in all_results if r['status'] == 'OK']
            if len(valid) >= 2:
                for key in ['k5', 'k7', 'k9', 'ensemble']:
                    vals = [r[f'{key}_nrf'] for r in valid if r.get(f'{key}_nrf') is not None]
                    if vals:
                        print(f"    {key}: nRF={st.mean(vals):.4f}±{st.stdev(vals):.4f} ({len(vals)} valid)")
    
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    fieldnames = ['seed', 'n', 'L', 'mu', 'k5_nrf', 'k7_nrf', 'k9_nrf', 'ensemble_nrf',
                  'k5_time_s', 'k7_time_s', 'k9_time_s', 'ensemble_time_s', 'status']
    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(all_results)
    
    valid = [r for r in all_results if r['status'] == 'OK']
    if valid:
        print(f"\nFinal ({len(valid)} seeds):")
        for key in ['k5', 'k7', 'k9', 'ensemble']:
            vals = [r[f'{key}_nrf'] for r in valid if r.get(f'{key}_nrf') is not None]
            if vals:
                if len(vals) >= 2:
                    print(f"  {key}: nRF={st.mean(vals):.4f}±{st.stdev(vals):.4f}")
                else:
                    print(f"  {key}: nRF={vals[0]:.4f}")
    print(f"Total: {time.time() - t0_total:.0f}s  Output: {args.output}")


if __name__ == '__main__':
    main()
