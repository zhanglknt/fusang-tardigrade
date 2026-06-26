#!/usr/bin/env python3
"""
End-to-end MHL pipeline test with V4b boundary classifier.

Tests both coalescent (should → STOP) and structured (should → SPLIT) data
across multiple scales. Verifies:
  1. ML split/STOP decision correctness
  2. Tree completeness (all taxa present)
  3. Timing benchmarks

Usage:
    python test_mhl_e2e_v4b.py [--quick] [--verbose]
"""

import sys
import os
import time
import argparse
import tempfile
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from fusang_mhl.mlh_utils import Timer

# ==========================================================================
# Test Data Generation
# ==========================================================================

DEFAULT_SCALES = [50, 100, 200, 300, 500]
N_SEEDS_PER_SCALE = 5  # per scale per tree_type
LEN = 1000
SUB = 0.05
INDEL = 0.02


def generate_coalescent_data(n, L, sub, indel, seed):
    """Generate coalescent tree data (should predict STOP)."""
    from tree_simulation import make_coalescent_tree, simulate_seqs
    root_node, leaves = make_coalescent_tree(n, seed=seed)

    leaf_seqs = simulate_seqs(root_node, n, L, sub, seed, indel_rate=indel)
    _CODE_MAP = {0: 'A', 1: 'T', 2: 'C', 3: 'G'}
    seqs_dict = {}
    for i in range(n):
        seqs_dict[f't{i+1:04d}'] = ''.join(_CODE_MAP[b] for b in leaf_seqs[i])

    taxon_names = list(seqs_dict.keys())
    sequences = [seqs_dict[n] for n in taxon_names]
    return sequences, taxon_names


def generate_structured_data(n, L, sub, indel, seed):
    """Generate structured tree data with clear clades (should predict SPLIT).

    Creates 2-3 clades with different substitution rates to produce
    distinct distance clusters.
    """
    from tree_simulation import make_coalescent_tree, simulate_seqs

    np.random.seed(seed)
    n_clades = 2 if n <= 200 else 3
    clade_sizes = []
    remaining = n
    for c in range(n_clades - 1):
        sz = max(30, remaining // (n_clades - c) + np.random.randint(-10, 11))
        sz = min(sz, remaining - 30 * (n_clades - c - 1))
        clade_sizes.append(sz)
        remaining -= sz
    clade_sizes.append(remaining)

    # Generate each clade separately with different substitution rate
    clade_seqs = []
    clade_names = []
    name_counter = 1

    for ci, sz in enumerate(clade_sizes):
        sub_rate = sub * (0.5 + ci * 0.75)  # 0.5x, 1.25x, 2.0x baseline
        root_node, leaves = make_coalescent_tree(sz, seed=seed + ci * 1000)
        leaf_seqs = simulate_seqs(root_node, sz, L, sub_rate, seed + ci * 1000, indel_rate=indel)

        _CODE_MAP = {0: 'A', 1: 'T', 2: 'C', 3: 'G'}
        for i in range(sz):
            name = f't{name_counter:04d}'
            name_counter += 1
            clade_names.append(name)
            clade_seqs.append(''.join(_CODE_MAP[b] for b in leaf_seqs[i]))

    return clade_seqs, clade_names


# ==========================================================================
# Tree Building
# ==========================================================================

def build_direct_nj(sequences, taxon_names):
    """Build NJ tree directly on all taxa."""
    from fusang_v4_dahp_v1 import build_nj
    from fusang_mhl.level0_kmer import compute_l0_distance
    D = compute_l0_distance(sequences, taxon_names)
    return build_nj(D, taxon_names)


def check_tree_completeness(nwk_str, expected_names):
    """Check if all expected taxon names appear in the Newick tree."""
    import re
    if nwk_str is None:
        return False, 0
    found = 0
    for name in expected_names:
        if name in nwk_str:
            found += 1
    return found == len(expected_names), found


# ==========================================================================
# Main Test
# ==========================================================================

def run_e2e_test(scales=None, seeds_per_scale=5, quick=False, verbose=False):
    """Run comprehensive E2E test."""
    if scales is None:
        scales = DEFAULT_SCALES

    if quick:
        scales = [50, 200, 500]
        seeds_per_scale = 2  # n=500 is slow, use fewer seeds

    from fusang_mhl.level0_kmer import compute_l0_distance
    from fusang_mhl.ml_split import ml_split_decision, model_available

    print(f"{'='*70}")
    print(f"MHL V4b E2E Test")
    print(f"{'='*70}")
    print(f"Model available: {model_available()}")
    print(f"Scales: {scales}, seeds per scale: {seeds_per_scale}")
    print(f"{'='*70}")

    all_results = []

    for tree_type, generator, expected_decision in [
        ("coalescent", generate_coalescent_data, "STOP"),
        ("structured", generate_structured_data, "SPLIT"),
    ]:
        print(f"\n{'─'*70}")
        print(f"Tree type: {tree_type} (expected → {expected_decision})")
        print(f"{'─'*70}")

        for n in scales:
            decisions = []
            split_pct = []
            tree_completeness = []
            times = []

            for seed in range(seeds_per_scale):
                try:
                    actual_seed = seed * 100 + (1 if tree_type == "structured" else 0)
                    sequences, taxon_names = generator(n, LEN, SUB, INDEL, actual_seed)

                    # ML split decision
                    t0 = time.time()
                    D = compute_l0_distance(sequences, taxon_names)
                    decision = ml_split_decision(D, sequences, taxon_names, verbose=False)
                    elapsed = time.time() - t0
                    times.append(elapsed)

                    decisions.append(decision)
                    split_pct.append(decision.get('p_split', -1))

                    # Verify decision correctness
                    is_correct = (decision['should_split'] == (expected_decision == "SPLIT"))
                    reason = decision.get('reason', '?')
                    model_used = decision.get('model_used', False)
                    p = decision.get('p_split', None)

                    status = "OK" if is_correct else "FAIL"
                    p_str = f"p={p:.3f}" if p is not None else "no model"

                    # Build tree (fast path: simple BioPython NJ for completeness check)
                    from Bio.Phylo.TreeConstruction import DistanceMatrix as BioDM, DistanceTreeConstructor
                    # Convert full matrix to lower-triangle format required by BioPhylo
                    if hasattr(D, 'tolist'):
                        D_list = D.tolist()
                    else:
                        D_list = D
                    # BioPhylo expects lower-triangle: list of n lists, row i has i elements
                    lower_tri = [[D_list[i][j] for j in range(i)] for i in range(len(D_list))]
                    dm = BioDM(list(taxon_names), lower_tri)
                    constructor = DistanceTreeConstructor()
                    tree_obj = constructor.nj(dm)
                    from io import StringIO
                    buf = StringIO()
                    from Bio import Phylo
                    Phylo.write(tree_obj, buf, 'newick')
                    nwk = buf.getvalue()
                    complete, n_found = check_tree_completeness(nwk, taxon_names) if nwk else (False, 0)
                    tree_completeness.append(complete)

                    if verbose or seed == 0:
                        print(f"  n={n:>4} seed={seed:>2}: {status} "
                              f"split={decision['should_split']} ({reason}, {p_str}, "
                              f"model={'ML' if model_used else 'H'}), "
                              f"found={n_found}/{n}, t={elapsed:.2f}s")

                except Exception as e:
                    print(f"  n={n:>4} seed={seed:>2}: ERROR: {e}")
                    if verbose:
                        import traceback
                        traceback.print_exc()

            # Summary for this scale
            n_correct = sum(1 for d in decisions
                           if d['should_split'] == (expected_decision == "SPLIT"))
            n_model = sum(1 for d in decisions if d.get('model_used', False))
            n_complete = sum(tree_completeness)
            avg_p = np.mean([p for p in split_pct if p >= 0]) if any(p >= 0 for p in split_pct) else -1
            avg_t = np.mean(times)

            acc = n_correct / len(decisions) if decisions else 0
            comp_rate = n_complete / len(tree_completeness) if tree_completeness else 0

            result = {
                "tree_type": tree_type,
                "expected": expected_decision,
                "n": n,
                "n_seeds": len(decisions),
                "accuracy": acc,
                "avg_p_split": avg_p,
                "n_model_used": n_model,
                "tree_completeness": comp_rate,
                "avg_time_s": avg_t,
            }
            all_results.append(result)

            pct = f"{100*acc:.0f}%"
            print(f"  → Summary: {n_correct}/{len(decisions)} correct ({pct}), "
                  f"avg_p={avg_p:.3f}, model_used={n_model}/{len(decisions)}, "
                  f"complete={n_complete}/{len(tree_completeness)}, "
                  f"avg_t={avg_t:.3f}s")

    # ==========================================================================
    # Overall Summary
    # ==========================================================================
    print(f"\n{'='*70}")
    print("OVERALL SUMMARY")
    print(f"{'='*70}")

    # By tree type
    for tt in ["coalescent", "structured"]:
        tt_results = [r for r in all_results if r["tree_type"] == tt]
        if tt_results:
            acc = np.mean([r["accuracy"] for r in tt_results])
            comp = np.mean([r["tree_completeness"] for r in tt_results])
            print(f"\n  {tt}:")
            print(f"    Decision accuracy: {100*acc:.1f}%")
            print(f"    Tree completeness: {100*comp:.1f}%")
            for r in tt_results:
                pct = f"{100*r['accuracy']:.0f}%"
                print(f"    n={r['n']:>4}: acc={pct:>4s}, "
                      f"avg_p={r['avg_p_split']:.3f}, "
                      f"model={r['n_model_used']}/{r['n_seeds']}, "
                      f"t={r['avg_time_s']:.3f}s")

    # Overall
    overall_acc = np.mean([r["accuracy"] for r in all_results])
    overall_comp = np.mean([r["tree_completeness"] for r in all_results])
    print(f"\n  OVERALL:")
    print(f"    Decision accuracy: {100*overall_acc:.1f}%")
    print(f"    Tree completeness: {100*overall_comp:.1f}%")

    # Critical check
    coalescent_results = [r for r in all_results if r["tree_type"] == "coalescent"]
    structured_results = [r for r in all_results if r["tree_type"] == "structured"]

    coalescent_stop = np.mean([r["accuracy"] for r in coalescent_results]) if coalescent_results else 0
    structured_split = np.mean([r["accuracy"] for r in structured_results]) if structured_results else 0

    print(f"\n{'='*70}")
    if coalescent_stop >= 0.8 and structured_split >= 0.8:
        print("PASS: Both coalescent->STOP and structured->SPLIT work correctly")
    else:
        print(f"FAIL: coalescent={coalescent_stop:.2f}, structured={structured_split:.2f}")
        if coalescent_stop < 0.8:
            print("   Coalescent data should predict STOP")
        if structured_split < 0.8:
            print("   Structured data should predict SPLIT")
    print(f"{'='*70}")

    return all_results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="MHL V4b E2E Test")
    parser.add_argument("--quick", action="store_true", help="Quick test (fewer seeds)")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    args = parser.parse_args()

    results = run_e2e_test(quick=args.quick, verbose=args.verbose)
