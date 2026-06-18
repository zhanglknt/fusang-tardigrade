#!/usr/bin/env python3
"""
benchmark_swisstree.py — Benchmark Fusang and competitors on AFproject SwissTree gene tree datasets.

This script runs alignment-free phylogenetic methods on the 11 SwissTree protein
gene families and compares against reference trees using nRF distance.

Methods tested:
  - Fusang spaced k-mer (k=5, gap2 = "11011")
  - Fusang spaced k-mer (k=6, gap2 = "110111")
  - KmerCosine contiguous k=3 (protein-appropriate)
  - KmerCosine contiguous k=4 (protein-appropriate)
  - KmerCosine contiguous k=5 (protein-appropriate)
  - Co-phylog (k=19 context-object matching)

Reference: Zielezinski et al. (2019) Genome Biology 20:144
"""

import os
import sys
import json
import time
import csv
from pathlib import Path
from collections import defaultdict
import numpy as np

# Add project root to path
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

from af_competitor_methods import (
    kmer_cosine_distance_matrix,
    cophylog_distance_matrix,
    andi_approx_distance_matrix,
)
from fusang_v4_dahp_v1 import build_nj
from calc_nrf_simple import get_bipartitions_from_newick


def compute_nrf(pred_nwk, ref_nwk):
    """Compute normalized Robinson-Foulds distance between two Newick strings."""
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


def read_multi_fasta(fasta_path, encoding='utf-8', errors='replace'):
    """Read multi-sequence FASTA file. Returns dict name -> sequence."""
    seqs = {}
    name = None
    seq_parts = []
    with open(fasta_path, encoding=encoding, errors=errors) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(seq_parts)
                name = line[1:].split()[0]  # First word after >
                seq_parts = []
            else:
                seq_parts.append(line)
    if name is not None:
        seqs[name] = ''.join(seq_parts)
    return seqs


def load_swisstree_family(swisstree_dir, ref_dir, family_id):
    """
    Load a SwissTree gene family:
    - Read individual FASTA files from swisstree_dir/{family_id}_*.fasta
    - Read reference tree from ref_dir/{family_id}.nwk
    - Read sequence IDs from ref_dir/{family_id}.ids.json

    Returns: (sequences dict, names list, ref_tree string, family_info dict)
    """
    # Load sequence IDs from ids.json
    ids_path = os.path.join(ref_dir, f"{family_id}.ids.json")
    with open(ids_path, encoding='utf-8') as f:
        ids_data = json.load(f)
    expected_ids = ids_data['seqids']

    # Load reference tree
    ref_path = os.path.join(ref_dir, f"{family_id}.nwk")
    with open(ref_path, encoding='utf-8') as f:
        ref_tree = f.read().strip()

    # Load sequences from individual FASTA files
    sequences = {}
    for seq_id in expected_ids:
        fasta_path = os.path.join(swisstree_dir, f"{seq_id}.fasta")
        if os.path.exists(fasta_path):
            with open(fasta_path, encoding='utf-8', errors='replace') as f:
                lines = [l.strip() for l in f if l.strip() and not l.startswith('>')]
            sequences[seq_id] = ''.join(lines).upper()
        else:
            print(f"  [WARN] Missing sequence file: {fasta_path}")

    # Build names list matching the reference tree order
    names = [sid for sid in expected_ids if sid in sequences]

    return sequences, names, ref_tree, {'n_taxa': len(names), 'ids': expected_ids}


def run_method(dist_func, sequences, names, method_name, **kwargs):
    """Run a distance method + NJ tree building. Returns (nrf, time_seconds, tree_string)."""
    t0 = time.time()
    try:
        # Convert dict to list if needed
        seq_list = [sequences[n] for n in names]
        dist_matrix = dist_func(seq_list, names, **kwargs)
        t_dist = time.time() - t0

        # Build NJ tree from distance matrix
        tree_str = build_nj(dist_matrix, names)
        t_total = time.time() - t0

        return t_total, dist_matrix, tree_str
    except Exception as e:
        print(f"    [{method_name}] Error: {e}")
        return None, None, None


def main():
    # Paths
    swisstree_dir = os.path.join(SCRIPT_DIR, "real_data", "swisstree", "swisstree")
    ref_dir = os.path.join(SCRIPT_DIR, "real_data", "swisstree", "ref_trees")

    if not os.path.isdir(swisstree_dir):
        print(f"ERROR: SwissTree directory not found: {swisstree_dir}")
        sys.exit(1)
    if not os.path.isdir(ref_dir):
        print(f"ERROR: Reference trees directory not found: {ref_dir}")
        sys.exit(1)

    families = ['ST001', 'ST002', 'ST003', 'ST004', 'ST005', 'ST007',
                'ST008', 'ST009', 'ST010', 'ST011', 'ST012']

    # Method definitions
    methods = {
        'fusang_k5_gap2': {
            'desc': 'Fusang spaced k=5,gap2 (11011)',
            'func': kmer_cosine_distance_matrix,
            'kwargs': {'k': 5, 'gap_pattern': '11011'},
        },
        'fusang_k4_gap1': {
            'desc': 'Fusang spaced k=4,gap1 (1011)',
            'func': kmer_cosine_distance_matrix,
            'kwargs': {'k': 4, 'gap_pattern': '1011'},
        },
        'fusang_k3_gap1': {
            'desc': 'Fusang spaced k=3,gap1 (101)',
            'func': kmer_cosine_distance_matrix,
            'kwargs': {'k': 3, 'gap_pattern': '101'},
        },
        'kmer_k3': {
            'desc': 'KmerCosine k=3 contiguous',
            'func': kmer_cosine_distance_matrix,
            'kwargs': {'k': 3},
        },
        'kmer_k4': {
            'desc': 'KmerCosine k=4 contiguous',
            'func': kmer_cosine_distance_matrix,
            'kwargs': {'k': 4},
        },
        'kmer_k5': {
            'desc': 'KmerCosine k=5 contiguous',
            'func': kmer_cosine_distance_matrix,
            'kwargs': {'k': 5},
        },
        'cophylog_k5': {
            'desc': 'Co-phylog halfctx=5 (k=11)',
            'func': cophylog_distance_matrix,
            'kwargs': {'halfctx': 5},
        },
        'cophylog_k3': {
            'desc': 'Co-phylog halfctx=3 (k=7)',
            'func': cophylog_distance_matrix,
            'kwargs': {'halfctx': 3},
        },
        'andi_approx_10': {
            'desc': 'andi-approx anchor=10',
            'func': andi_approx_distance_matrix,
            'kwargs': {'min_anchor': 10, 'use_revcomp': False},
        },
        'andi_approx_15': {
            'desc': 'andi-approx anchor=15',
            'func': andi_approx_distance_matrix,
            'kwargs': {'min_anchor': 15, 'use_revcomp': False},
        },
    }

    # Results storage
    results = {}  # family -> method -> {nrf, time, n_taxa}

    print("=" * 80)
    print("SwissTree Gene Tree Benchmark (AFproject)")
    print("=" * 80)
    print(f"Data dir: {swisstree_dir}")
    print(f"Ref dir:  {ref_dir}")
    print(f"Families: {len(families)}")
    print(f"Methods:  {len(methods)}")
    print()

    for family in families:
        print(f"\n{'='*60}")
        print(f"Family: {family}")
        print(f"{'='*60}")

        # Load data
        sequences, names, ref_tree, info = load_swisstree_family(
            swisstree_dir, ref_dir, family
        )
        n_taxa = info['n_taxa']

        # Compute average sequence length
        seq_lens = [len(s) for s in sequences.values()]
        avg_len = np.mean(seq_lens) if seq_lens else 0

        print(f"  Taxa: {n_taxa}, Avg length: {avg_len:.0f} aa")
        print(f"  Ref tree taxa: {len(info['ids'])}")

        if n_taxa < 4:
            print(f"  SKIP: too few taxa ({n_taxa})")
            continue

        # Validate reference tree bipartitions
        try:
            _, ref_bip = get_bipartitions_from_newick(ref_tree)
            if not ref_bip:
                print(f"  SKIP: empty reference bipartitions")
                continue
            print(f"  Ref tree bipartitions: {len(ref_bip)}")
        except Exception as e:
            print(f"  SKIP: cannot parse reference tree: {e}")
            continue

        results[family] = {}
        results[family]['_meta'] = {
            'n_taxa': n_taxa,
            'avg_seq_len': avg_len,
            'ref_bipartitions': len(ref_bip),
        }

        for method_id, method_info in methods.items():
            desc = method_info['desc']
            func = method_info['func']
            kwargs = method_info['kwargs']

            # Skip Co-phylog for very short sequences (context too long)
            if 'cophylog' in method_id and avg_len < 20:
                print(f"  [{desc}] SKIP: sequences too short for context length")
                continue

            print(f"\n  [{desc}]")
            t_total, dist_matrix, tree_str = run_method(
                func, sequences, names, method_id, **kwargs
            )

            if tree_str is None:
                results[family][method_id] = {'nrf': float('nan'), 'time': 0}
                print(f"    FAILED")
                continue

            nrf = compute_nrf(tree_str, ref_tree)
            results[family][method_id] = {'nrf': nrf, 'time': t_total}

            if not np.isnan(nrf):
                print(f"    nRF = {nrf:.4f} ({t_total:.1f}s)")
            else:
                print(f"    nRF = NaN ({t_total:.1f}s)")

    # ---- Summary ----
    print("\n\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    # Build summary table
    method_ids = list(methods.keys())
    header = ["Family", "Taxa", "AvgLen"] + [methods[m]['desc'].split('(')[0].strip() for m in method_ids]
    rows = []

    for family in families:
        if family not in results:
            continue
        meta = results[family]['_meta']
        row = [family, str(meta['n_taxa']), f"{meta['avg_seq_len']:.0f}"]
        for mid in method_ids:
            if mid in results[family]:
                nrf = results[family][mid]['nrf']
                row.append(f"{nrf:.4f}" if not np.isnan(nrf) else "NaN")
            else:
                row.append("SKIP")
        rows.append(row)

    # Print table
    col_widths = [max(len(str(h)), max(len(str(r[i])) for r in rows)) for i, h in enumerate(header)]
    fmt = "  ".join(f"{{:<{w}}}" for w in col_widths)
    print(fmt.format(*header))
    print("-" * sum(col_widths + [2] * (len(col_widths) - 1)))
    for row in rows:
        print(fmt.format(*[str(x) for x in row]))

    # Compute mean nRF across families for each method
    print("\n--- Mean nRF across families ---")
    for mid in method_ids:
        nrf_vals = [results[f][mid]['nrf'] for f in families
                     if f in results and mid in results[f] and
                     not np.isnan(results[f][mid]['nrf'])]
        if nrf_vals:
            mean_nrf = np.mean(nrf_vals)
            std_nrf = np.std(nrf_vals)
            n_families = len(nrf_vals)
            print(f"  {methods[mid]['desc']:40s}  mean={mean_nrf:.4f} +/- {std_nrf:.4f}  (n={n_families} families)")

    # Paired statistical tests
    print("\n--- Paired Wilcoxon tests ---")
    try:
        from scipy.stats import wilcoxon

        # Key comparisons
        comparisons = [
            ('fusang_k4_gap1', 'kmer_k3', 'Spaced k=4,gap1 vs Contiguous k=3'),
            ('fusang_k4_gap1', 'kmer_k4', 'Spaced k=4,gap1 vs Contiguous k=4'),
            ('fusang_k4_gap1', 'cophylog_k5', 'K-mer cosine vs Co-phylog'),
            ('kmer_k3', 'cophylog_k5', 'Contiguous k=3 vs Co-phylog'),
        ]
        for mid1, mid2, label in comparisons:
            pairs = [(results[f][mid1]['nrf'], results[f][mid2]['nrf'])
                     for f in families
                     if f in results and mid1 in results[f] and mid2 in results[f]
                     and not np.isnan(results[f][mid1]['nrf'])
                     and not np.isnan(results[f][mid2]['nrf'])]
            if pairs and len(pairs) >= 5:
                v1 = [p[0] for p in pairs]
                v2 = [p[1] for p in pairs]
                diff = np.array(v1) - np.array(v2)
                cohens_d = np.mean(diff) / np.std(diff) if np.std(diff) > 0 else 0
                try:
                    stat, pval = wilcoxon(v1, v2, alternative='two-sided')
                    sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "n.s."
                    print(f"  {label:50s}  d={cohens_d:+.2f}  p={pval:.4f} {sig}")
                except Exception as e:
                    print(f"  {label:50s}  d={cohens_d:+.2f}  test_error: {e}")
    except ImportError:
        print("  scipy not available, skipping Wilcoxon tests")

    # Save CSV
    csv_path = os.path.join(SCRIPT_DIR, "benchmark_swisstree_results.csv")
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['family', 'n_taxa', 'avg_seq_len'] + method_ids)
        for family in families:
            if family not in results:
                continue
            meta = results[family]['_meta']
            row = [family, meta['n_taxa'], f"{meta['avg_seq_len']:.0f}"]
            for mid in method_ids:
                if mid in results[family]:
                    nrf = results[family][mid]['nrf']
                    row.append(f"{nrf:.6f}" if not np.isnan(nrf) else "NaN")
                else:
                    row.append("SKIP")
            writer.writerow(row)
    print(f"\nResults saved to: {csv_path}")

    # Per-family detail output
    print("\n--- Per-family details ---")
    for family in families:
        if family not in results:
            continue
        meta = results[family]['_meta']
        best_method = None
        best_nrf = 1.0
        for mid in method_ids:
            if mid in results[family] and not np.isnan(results[family][mid]['nrf']):
                if results[family][mid]['nrf'] < best_nrf:
                    best_nrf = results[family][mid]['nrf']
                    best_method = methods[mid]['desc']
        if best_method:
            print(f"  {family} ({meta['n_taxa']} taxa, {meta['avg_seq_len']:.0f} aa): "
                  f"best = {best_method} (nRF={best_nrf:.4f})")


if __name__ == '__main__':
    main()
