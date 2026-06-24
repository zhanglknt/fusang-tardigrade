#!/usr/bin/env python3
"""
Generate nRF values for scalability benchmark scales.
D11 fix: Add nRF to scalability results.

For each n in [200, 500, 1000, 2000] (smaller scales for nRF comparison):
- Generate coalescent tree + JC69 sequences
- Build Fusang NJ tree
- Compute nRF vs true tree

NOTE: Uses JC69 simulator (tree_simulation.py), NOT INDELible.
      nRF values will differ from paper's INDELible-based results.
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


def benchmark_one(n, seed, L=500, mu=0.05):
    """Run one scalability+nRF benchmark."""
    root, nodes = make_coalescent_tree(n, seed=seed)
    
    leaf_nodes = [nd for nd in nodes if nd.is_leaf]
    leaf_nodes.sort(key=lambda nd: nd.idx)
    taxon_names = [nd.name for nd in leaf_nodes]
    
    seq_array = simulate_seqs(root, n, L, mu=mu, seed=seed + 1000)
    sequences = [seq_array_to_string(seq_array[i]) for i in range(n)]
    
    t0 = time.time()
    dist_matrix = compute_kmer_distance_matrix(
        sequences, taxon_names, k=5, metric='cosine',
        gap_pattern=(0, 1, 2, 5, 6),
        n_threads=1
    )
    nj_newick = build_nj(dist_matrix, taxon_names)
    fusang_time = time.time() - t0
    
    if nj_newick is None:
        return {'n': n, 'seed': seed, 'status': 'NJ_FAIL', 'time_s': round(fusang_time, 2)}
    
    true_nwk = root.to_newick() + ';'
    tmp_true = tempfile.mktemp(suffix='.nwk')
    tmp_nj = tempfile.mktemp(suffix='.nwk')
    with open(tmp_true, 'w') as f: f.write(true_nwk)
    with open(tmp_nj, 'w') as f: f.write(nj_newick + ';\n')
    
    try:
        nrf_val = calc_nrf(tmp_true, tmp_nj)
    except Exception as e:
        print(f"  n={n} seed={seed}: nRF error: {e}")
        nrf_val = None
    
    for p in [tmp_true, tmp_nj]:
        try: os.remove(p)
        except: pass
    
    return {'n': n, 'seed': seed, 'nrf': nrf_val, 'time_s': round(fusang_time, 2),
            'status': 'OK' if nrf_val is not None else 'NRF_FAIL'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--scales', type=int, nargs='+', default=[200, 500, 1000, 2000])
    ap.add_argument('--seeds-per-scale', type=int, default=5)
    ap.add_argument('--output', default='data/scalability_nrf.json')
    args = ap.parse_args()
    
    all_results = []
    n_total = len(args.scales) * args.seeds_per_scale
    print(f"Scalability nRF benchmark: scales={args.scales}, {args.seeds_per_scale} seeds each")
    print("NOTE: JC69 simulator (not INDELible). Results differ from paper.")
    
    t0_total = time.time()
    count = 0
    
    for n in args.scales:
        for i in range(args.seeds_per_scale):
            seed = 1000 + i  # offset seeds
            count += 1
            t0 = time.time()
            result = benchmark_one(n, seed)
            elapsed = time.time() - t0
            all_results.append(result)
            
            nrf_str = f"{result['nrf']:.4f}" if result.get('nrf') is not None else 'FAIL'
            eta = (time.time() - t0_total) / count * (n_total - count)
            print(f"  [{count}/{n_total}] n={n} seed={seed}: nRF={nrf_str}, {elapsed:.1f}s, ETA {eta:.0f}s")
    
    # Summary by scale
    scale_summary = {}
    for n in args.scales:
        scale_results = [r for r in all_results if r['n'] == n and r['status'] == 'OK']
        if scale_results:
            nrf_vals = [r['nrf'] for r in scale_results]
            scale_summary[str(n)] = {
                'n_seeds': len(scale_results),
                'nrf_mean': round(st.mean(nrf_vals), 4),
                'nrf_std': round(st.stdev(nrf_vals), 4) if len(nrf_vals) >= 2 else 0,
                'time_mean_s': round(st.mean([r['time_s'] for r in scale_results]), 2),
                'comment': 'JC69 simulator (not INDELible). Clean (no-indel) data only.'
            }
    
    output = {'summary': scale_summary, 'details': all_results,
              'note': 'Generated with tree_simulation.py (JC69, coalescent). Different from paper (INDELible, GTR+Γ).'}
    
    with open(args.output, 'w') as f:
        json.dump(output, f, indent=2)
    
    print(f"\nScale summaries:")
    for n_str, summary in sorted(scale_summary.items(), key=lambda x: int(x[0])):
        print(f"  n={n_str}: nRF={summary['nrf_mean']:.4f}±{summary['nrf_std']:.4f} "
              f"({summary['n_seeds']} seeds, {summary['time_mean_s']}s avg)")
    print(f"Total: {time.time() - t0_total:.0f}s  Output: {args.output}")


if __name__ == '__main__':
    main()
