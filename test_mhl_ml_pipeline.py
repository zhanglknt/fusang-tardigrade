#!/usr/bin/env python3
"""
End-to-end test: MHL pipeline with ML-driven L0 split decision.

Tests whether the boundary classifier's split/stop decision leads to
better trees when the full MHL pipeline (L0 split → cluster NJ → merge)
is compared against direct NJ on the full dataset.

Usage:
    python test_mhl_ml_pipeline.py [--seeds 20] [--verbose]
"""

import sys
import os
import time
import argparse
import tempfile
import numpy as np
from typing import Dict, List, Tuple, Optional

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from fusang_mhl.level0_kmer import (
    run_level0, compute_l0_distance, cluster_taxa,
    select_centroids, build_backbone_nj,
)
from fusang_mhl.mlh_utils import Timer


# ==========================================================================
# Test Data Generation
# ==========================================================================

DEFAULT_TEST_CONFIGS = [
    # (n, L, sub, indel, n_seeds) — small tests first
    (50,   500,  0.05, 0.02, 10),
    (100,  500,  0.05, 0.02, 10),
    (200,  500,  0.05, 0.02, 10),
    (300,  1000, 0.05, 0.02, 5),
    (500,  1000, 0.05, 0.02, 5),
]


def generate_test_data(n, L, sub, indel, seed):
    """Generate one test dataset with known true tree."""
    from tree_simulation import make_coalescent_tree, simulate_seqs
    root_node, leaves = make_coalescent_tree(n, seed=seed)

    # True tree as Newick
    def _node_to_nwk(node):
        if node.left is None and node.right is None:
            return node.name or f't{node.idx}'
        left_str = _node_to_nwk(node.left)
        right_str = _node_to_nwk(node.right)
        bl = f':{node.branch_length:.6f}' if hasattr(node, 'branch_length') else ''
        return f'({left_str},{right_str}){node.name or ""}{bl}'
    true_nwk = _node_to_nwk(root_node) + ';'

    # Simulate sequences
    leaf_seqs = simulate_seqs(root_node, n, L, sub, seed, indel_rate=indel)

    # Convert to lists
    _CODE_MAP = {0: 'A', 1: 'T', 2: 'C', 3: 'G'}
    seqs_dict = {}
    for i in range(n):
        seqs_dict[f't{i+1:04d}'] = ''.join(_CODE_MAP[b] for b in leaf_seqs[i])

    taxon_names = list(seqs_dict.keys())
    sequences = [seqs_dict[n] for n in taxon_names]
    return true_nwk, sequences, taxon_names


# ==========================================================================
# Tree Building Methods
# ==========================================================================

def build_direct_nj(sequences, taxon_names):
    """Build NJ tree directly on all taxa (baseline)."""
    from fusang_v4_dahp_v1 import build_nj
    D = compute_l0_distance(sequences, taxon_names)
    return build_nj(D, taxon_names)


def build_mhl_split_nj(sequences, taxon_names, verbose=False):
    """MHL pipeline: L0 ML split → per-cluster NJ → merge."""
    from fusang_v4_dahp_v1 import build_nj

    n = len(taxon_names)

    # L0: compute distance + ML split decision
    D = compute_l0_distance(sequences, taxon_names)

    # Get ML split decision
    from fusang_mhl.ml_split import ml_split_decision
    decision = ml_split_decision(D, sequences, taxon_names, verbose=verbose)

    if verbose:
        reason = decision.get('reason', 'unknown')
        p = decision.get('p_split', None)
        p_str = f', p={p:.3f}' if p is not None else ''
        print(f"  [test] ML decision: split={decision['should_split']} "
              f"({reason}{p_str}), groups={decision['target_groups']}",
              file=sys.stderr)

    if not decision['should_split']:
        # No split: build NJ on full dataset
        return build_nj(D, taxon_names), {"split": False, "n_clusters": 1}

    # Split: cluster into groups
    clusters = cluster_taxa(
        D, taxon_names,
        max_group_size=decision['max_group_size'],
        target_groups=decision['target_groups'],
    )

    if len(clusters) <= 1:
        return build_nj(D, taxon_names), {"split": False, "n_clusters": 1}

    if verbose:
        sizes = [len(c) for c in clusters]
        print(f"  [test] Clustered into {len(clusters)} groups: {sizes}",
              file=sys.stderr)

    # Build NJ per cluster
    child_nwks = []
    for ci, cluster_indices in enumerate(clusters):
        sub_seqs = [sequences[i] for i in cluster_indices]
        sub_names = [taxon_names[i] for i in cluster_indices]
        sub_D = D[np.ix_(cluster_indices, cluster_indices)]
        child_nwk = build_nj(sub_D, sub_names)
        if child_nwk is not None:
            child_nwks.append((cluster_indices, child_nwk))
        elif verbose:
            print(f"  [test] Warning: NJ failed for cluster {ci}", file=sys.stderr)

    if len(child_nwks) <= 1:
        return build_nj(D, taxon_names), {"split": False, "n_clusters": 1}

    # Merge children using merge_level
    from fusang_mhl.merger import merge_level
    all_indices = list(range(n))
    merged_nwk = merge_level(
        all_indices, child_nwks,
        D, taxon_names,
        verbose=verbose,
    )

    if merged_nwk is None:
        return build_nj(D, taxon_names), {"split": False, "n_clusters": 1}

    return merged_nwk, {
        "split": True,
        "n_clusters": len(clusters),
        "p_split": decision.get('p_split'),
        "model_used": decision.get('model_used', False),
    }


# ==========================================================================
# nRF Computation
# ==========================================================================

def compute_nrf(tree1_nwk, tree2_nwk, norm="2*(n-3)"):
    """Compute normalized Robinson-Foulds distance."""
    # Parse both trees
    def parse_nwk(nwk_str):
        from fusang_v2 import TreeNode
        return TreeNode.from_newick(nwk_str)

    def get_bipartitions(root):
        """Get all non-trivial bipartitions from a rooted tree."""
        from collections import Counter
        # Assign IDs to leaves
        leaves = {}
        leaf_id = 0
        def assign_ids(node):
            nonlocal leaf_id
            if not node.children:
                if node.name not in leaves:
                    leaves[node.name] = leaf_id
                    leaf_id += 1
                return {leaves[node.name]}
            s = set()
            for ch in node.children:
                s.update(assign_ids(ch))
            # Skip trivial bipartitions (single leaf or all leaves)
            if 1 < len(s) < leaf_id:
                # Canonicalize: use frozenset of the smaller side
                all_leaves = frozenset(range(leaf_id))
                if len(s) <= leaf_id // 2:
                    bipartitions.append(frozenset(s))
                else:
                    bipartitions.append(all_leaves - s)
            return s

        bipartitions = []
        assign_ids(root)
        return bipartitions

    try:
        root1 = parse_nwk(tree1_nwk)
        root2 = parse_nwk(tree2_nwk)
        bp1 = set(get_bipartitions(root1))
        bp2 = set(get_bipartitions(root2))

        # Find leaves
        def get_leaves(root):
            if not root.children:
                return {root.name}
            s = set()
            for ch in root.children:
                s.update(get_leaves(ch))
            return s
        leaves1 = get_leaves(root1)
        leaves2 = get_leaves(root2)
        all_leaves = leaves1 | leaves2
        n_leaves = len(all_leaves)

        # RF distance = |bp1 Δ bp2|
        rf = len(bp1.symmetric_difference(bp2))

        # Normalize
        if norm == "2*(n-3)":
            max_rf = 2 * (n_leaves - 3) if n_leaves > 3 else 1
        else:
            max_rf = len(bp1) + len(bp2) if (len(bp1) + len(bp2)) > 0 else 1

        return rf / max_rf
    except Exception as e:
        print(f"  [nRF] Error: {e}", file=sys.stderr)
        return 1.0


# ==========================================================================
# Main Test
# ==========================================================================

def run_pipeline_test(configs=None, verbose=False):
    """Run full MHL ML pipeline test across multiple configs and seeds."""
    if configs is None:
        configs = DEFAULT_TEST_CONFIGS

    results = []

    for n, L, sub, indel, n_seeds in configs:
        print(f"\n{'='*60}")
        print(f"Testing n={n}, L={L}, sub={sub}, indel={indel}, seeds={n_seeds}")
        print(f"{'='*60}")

        direct_nrfs = []
        mhl_nrfs = []
        ml_decisions = []

        for seed in range(n_seeds):
            try:
                true_nwk, sequences, taxon_names = generate_test_data(
                    n, L, sub, indel, seed=seed * 100,
                )

                # Direct NJ
                direct_nwk = build_direct_nj(sequences, taxon_names)
                if direct_nwk is None:
                    print(f"  Seed {seed}: direct NJ failed, skipping")
                    continue

                # MHL with ML split
                mhl_nwk, meta = build_mhl_split_nj(
                    sequences, taxon_names, verbose=verbose,
                )

                # Compute nRF vs true tree
                nrf_direct = compute_nrf(direct_nwk, true_nwk)
                nrf_mhl = compute_nrf(mhl_nwk, true_nwk)

                direct_nrfs.append(nrf_direct)
                mhl_nrfs.append(nrf_mhl)
                ml_decisions.append(meta)

                if verbose or seed == 0:
                    delta = nrf_mhl - nrf_direct
                    sign = "+" if delta > 0 else ""
                    print(f"  Seed {seed}: direct nRF={nrf_direct:.4f}, "
                          f"mhl nRF={nrf_mhl:.4f} ({sign}{delta:.4f}), "
                          f"split={meta.get('split')}",
                          file=sys.stderr)

            except Exception as e:
                print(f"  Seed {seed}: ERROR: {e}", file=sys.stderr)
                import traceback
                traceback.print_exc(file=sys.stderr)
                continue

        if direct_nrfs:
            d_mean = np.mean(direct_nrfs)
            d_std = np.std(direct_nrfs)
            m_mean = np.mean(mhl_nrfs)
            m_std = np.std(mhl_nrfs)
            n_splits = sum(1 for m in ml_decisions if m.get('split'))
            delta = m_mean - d_mean
            sign = "+" if delta > 0 else ""

            # Statistical test
            from scipy import stats
            t_stat, p_val = stats.ttest_rel(mhl_nrfs, direct_nrfs)
            d_cohen = (m_mean - d_mean) / np.sqrt((np.var(mhl_nrfs) + np.var(direct_nrfs)) / 2)

            print(f"\n  Summary (n={n}, {len(direct_nrfs)} valid seeds):")
            print(f"    Direct NJ nRF: {d_mean:.4f} +/- {d_std:.4f}")
            print(f"    MHL ML   nRF:   {m_mean:.4f} +/- {m_std:.4f}")
            print(f"    Delta:         {sign}{delta:.4f} (p={p_val:.4f}, d={d_cohen:.2f})")
            print(f"    Split decisions: {n_splits}/{len(ml_decisions)} "
                  f"({100*n_splits/len(ml_decisions):.0f}%)")

            n_heuristic = sum(1 for m in ml_decisions if not m.get('model_used', False))
            n_ml = sum(1 for m in ml_decisions if m.get('model_used', False))
            print(f"    Model used: {n_ml}, Heuristic: {n_heuristic}")

            results.append({
                "n": n, "L": L, "sub": sub, "indel": indel,
                "n_valid": len(direct_nrfs),
                "direct_nrf_mean": float(d_mean), "direct_nrf_std": float(d_std),
                "mhl_nrf_mean": float(m_mean), "mhl_nrf_std": float(m_std),
                "p_value": float(p_val), "cohens_d": float(d_cohen),
                "n_splits": n_splits, "n_total": len(ml_decisions),
                "n_model_used": n_ml, "n_heuristic": n_heuristic,
            })

    # Final summary
    print(f"\n{'='*60}")
    print("OVERALL SUMMARY")
    print(f"{'='*60}")
    print(f"{'n':>6} {'Valid':>6} {'Direct':>10} {'MHL':>10} {'Delta':>10} {'p':>8} {'Split%':>8}")
    print("-" * 60)
    for r in results:
        delta = r['mhl_nrf_mean'] - r['direct_nrf_mean']
        sign = "+" if delta > 0 else ""
        split_pct = 100 * r['n_splits'] / r['n_total']
        print(f"{r['n']:>6} {r['n_valid']:>6} "
              f"{r['direct_nrf_mean']:>10.4f} {r['mhl_nrf_mean']:>10.4f} "
              f"{sign}{delta:>9.4f} {r['p_value']:>8.4f} {split_pct:>7.1f}%")

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="MHL ML Pipeline E2E Test")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Verbose output")
    parser.add_argument("--quick", action="store_true",
                        help="Quick test (fewer seeds)")
    args = parser.parse_args()

    if args.quick:
        configs = [
            (50, 500, 0.05, 0.02, 3),
            (200, 500, 0.05, 0.02, 3),
        ]
    else:
        configs = DEFAULT_TEST_CONFIGS

    results = run_pipeline_test(configs, verbose=args.verbose)
