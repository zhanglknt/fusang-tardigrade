#!/usr/bin/env python3
"""
Extended MHL V4b E2E Test — 50+ independent scenarios for P0-2 fix.

Expands the original 12-scenario E2E test to 50+ scenarios with:
  - Multiple evolutionary parameter combinations (indel rate, substitution rate)
  - Varied sequence lengths
  - Both coalescent and structured tree priors
  - Scales: n=50, 100, 200, 300, 500

Reports results with Wilson confidence intervals instead of "100%".

Usage:
    python test_mhl_e2e_extended.py
"""

import sys
import os
import time
import argparse
import json
import pathlib
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from fusang_mhl.mlh_utils import Timer

# ==========================================================================
# Configuration — 50+ scenarios
# ==========================================================================

SCALES = [50, 100, 200, 300, 500]
N_SEEDS_PER = 3  # per (scale, tree_type, param_combo)

# Parameter combinations for varied scenarios
PARAM_COMBOS = [
    # (sub_rate, indel_rate, seq_len, label)
    (0.05, 0.02, 1000, "default"),
    (0.05, 0.01, 1000, "low_indel"),
    (0.05, 0.005, 1000, "vlow_indel"),
    (0.03, 0.02, 500, "short_seq"),
    (0.03, 0.01, 500, "short_lowindel"),
]

TREE_TYPES = ["coalescent", "structured"]  # coalescent → STOP, structured → SPLIT

# Set to a directory of pre-generated INDELible datasets to replace the
# built-in JC69 simulator (see plan_indelible_jobs.py / wsl_run_jobs.py).
INDELIBLE_DIR = None
OUT_JSON = "e2e_extended_results.json"

# Total scenarios = len(SCALES) × N_SEEDS_PER × len(PARAM_COMBOS) × len(TREE_TYPES)
# = 5 × 3 × 5 × 2 = 150 scenarios


def wilson_ci(n_success, n_total, alpha=0.05):
    """Wilson score interval for binomial proportion."""
    if n_total == 0:
        return (0.0, 0.0)
    import math
    from scipy.stats import norm
    z = norm.ppf(1 - alpha / 2)
    p = n_success / n_total
    denom = 1 + z**2 / n_total
    center = (p + z**2 / (2 * n_total)) / denom
    margin = z * math.sqrt((p * (1 - p) + z**2 / (4 * n_total)) / n_total) / denom
    return (max(0, center - margin), min(1, center + margin))


def generate_coalescent_data(n, L, sub, indel, seed,
                             indelible_dir=None, tag=None):
    if indelible_dir is not None:
        from fusang_mhl.indelible_adapter import load_indelible_tree
        nwk, seqs, names = load_indelible_tree(indelible_dir, tag, seed)
        return seqs, names
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


def generate_structured_data(n, L, sub, indel, seed,
                             indelible_dir=None, clade_jobs=None):
    if indelible_dir is not None:
        from fusang_mhl.indelible_adapter import load_indelible_structured
        _nwk, seqs, names = load_indelible_structured(indelible_dir, clade_jobs)
        return seqs, names
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
    clade_seqs = []
    clade_names = []
    name_counter = 1
    for ci, sz in enumerate(clade_sizes):
        sub_rate = sub * (0.5 + ci * 0.75)
        root_node, leaves = make_coalescent_tree(sz, seed=seed + ci * 1000)
        leaf_seqs = simulate_seqs(root_node, sz, L, sub_rate, seed + ci * 1000, indel_rate=indel)
        _CODE_MAP = {0: 'A', 1: 'T', 2: 'C', 3: 'G'}
        for i in range(sz):
            name = f't{name_counter:04d}'
            name_counter += 1
            clade_names.append(name)
            clade_seqs.append(''.join(_CODE_MAP[b] for b in leaf_seqs[i]))
    return clade_seqs, clade_names


def check_tree_completeness(nwk_str, expected_names):
    if nwk_str is None:
        return False, 0
    found = 0
    for name in expected_names:
        if name in nwk_str:
            found += 1
    return found == len(expected_names), found


def run_extended_e2e():
    from fusang_mhl.level0_kmer import compute_l0_distance
    from fusang_mhl.ml_split import ml_split_decision, model_available

    print("=" * 70)
    print("MHL V4b Extended E2E Test — P0-2 Fix")
    print("=" * 70)
    print(f"Model available: {model_available()}")
    print(f"Scales: {SCALES}")
    print(f"Param combos: {len(PARAM_COMBOS)}")
    print(f"Seeds per (scale, type, combo): {N_SEEDS_PER}")
    total = len(SCALES) * N_SEEDS_PER * len(PARAM_COMBOS) * len(TREE_TYPES)
    print(f"Total scenarios: {total}")
    print("=" * 70)

    all_individual = []  # per-scenario results
    summary = []
    failures = []  # explicitly recorded, never silently skipped

    for tree_type in TREE_TYPES:
        expected = "STOP" if tree_type == "coalescent" else "SPLIT"
        generator = generate_coalescent_data if tree_type == "coalescent" else generate_structured_data

        print(f"\n{'─' * 70}")
        print(f"Tree type: {tree_type} (expected → {expected})")
        print(f"{'─' * 70}")

        for sub, indel, L, label in PARAM_COMBOS:
            for n in SCALES:
                decisions_ok = []
                times = []
                complete_ok = []
                p_splits = []
                model_used_count = 0

                for seed_i in range(N_SEEDS_PER):
                    seed = seed_i * 100 + (1 if tree_type == "structured" else 0)
                    kwargs = {}
                    if INDELIBLE_DIR is not None:
                        from fusang_mhl.indelible_adapter import tag_for_e2e
                        if tree_type == "coalescent":
                            kwargs = dict(
                                indelible_dir=INDELIBLE_DIR,
                                tag=tag_for_e2e("coalescent", n, label, seed_i))
                        else:
                            nc = 2 if n <= 200 else 3
                            kwargs = dict(
                                indelible_dir=INDELIBLE_DIR,
                                clade_jobs=[(tag_for_e2e("structured", n, label,
                                                         seed_i, ci),
                                              seed + ci * 1000)
                                            for ci in range(nc)])
                    try:
                        sequences, taxon_names = generator(n, L, sub, indel, seed,
                                                           **kwargs)
                        t0 = time.time()
                        D = compute_l0_distance(sequences, taxon_names)
                        decision = ml_split_decision(D, sequences, taxon_names, verbose=False)
                        elapsed = time.time() - t0
                        times.append(elapsed)
                        p_splits.append(decision.get('p_split', -1))
                        if decision.get('model_used', False):
                            model_used_count += 1
                        correct = (decision['should_split'] == (expected == "SPLIT"))
                        decisions_ok.append(correct)

                        # Build NJ tree and check completeness
                        from Bio.Phylo.TreeConstruction import DistanceMatrix as BioDM, DistanceTreeConstructor
                        if hasattr(D, 'tolist'):
                            D_list = D.tolist()
                        else:
                            D_list = D
                        lower_tri = [[D_list[i][j] for j in range(i+1)] for i in range(len(D_list))]
                        dm = BioDM(list(taxon_names), lower_tri)
                        constructor = DistanceTreeConstructor()
                        tree_obj = constructor.nj(dm)
                        from io import StringIO
                        buf = StringIO()
                        from Bio import Phylo
                        Phylo.write(tree_obj, buf, 'newick')
                        nwk = buf.getvalue()
                        complete, n_found = check_tree_completeness(nwk, taxon_names) if nwk else (False, 0)
                        complete_ok.append(complete)

                        all_individual.append({
                            "tree_type": tree_type, "expected": expected,
                            "n": n, "param": label, "seed": seed_i,
                            "correct": correct, "complete": complete,
                            "p_split": decision.get('p_split', -1),
                            "model_used": decision.get('model_used', False),
                            "time_s": elapsed,
                        })

                    except Exception as e:
                        import traceback
                        failures.append({
                            "tree_type": tree_type, "n": n, "param": label,
                            "seed_i": seed_i,
                            "error": "%s: %s" % (type(e).__name__, e),
                            "traceback": traceback.format_exc()[-1500:],
                        })
                        print(f"  ERROR: n={n} {label} seed={seed_i}: "
                              f"{type(e).__name__}: {e}", flush=True)

                # Per-combination summary
                n_ok = sum(decisions_ok)
                n_total = len(decisions_ok)
                acc = n_ok / n_total if n_total > 0 else 0
                ci_lo, ci_hi = wilson_ci(n_ok, n_total)
                avg_p = np.mean(p_splits) if p_splits else -1
                avg_t = np.mean(times) if times else 0
                comp = sum(complete_ok) / len(complete_ok) if complete_ok else 0

                summary.append({
                    "tree_type": tree_type, "expected": expected,
                    "n": n, "param": label,
                    "n_seeds": n_total, "correct": n_ok,
                    "accuracy": acc, "ci_lower": ci_lo, "ci_upper": ci_hi,
                    "avg_p_split": avg_p, "model_used": model_used_count,
                    "tree_completeness": comp, "avg_time_s": avg_t,
                })

                ci_str = f"[{ci_lo:.1%}, {ci_hi:.1%}]"
                print(f"  n={n:>4} {label:<15s}: {n_ok}/{n_total} correct "
                      f"(acc={acc:.0%}, CI={ci_str}), avg_p={avg_p:.3f}, t={avg_t:.2f}s")

    # ==========================================================================
    # Overall Summary with Wilson CIs
    # ==========================================================================
    print(f"\n{'=' * 70}")
    print("OVERALL SUMMARY")
    print(f"{'=' * 70}")

    # By tree type
    for tt in ["coalescent", "structured"]:
        tt_summary = [s for s in summary if s["tree_type"] == tt]
        tt_indiv = [i for i in all_individual if i["tree_type"] == tt]
        n_correct = sum(s["correct"] for s in tt_summary)
        n_total = sum(s["n_seeds"] for s in tt_summary)
        acc = n_correct / n_total if n_total else 0
        ci_lo, ci_hi = wilson_ci(n_correct, n_total)
        comps = [i["complete"] for i in tt_indiv]
        comp_rate = sum(comps) / len(comps) if comps else 0
        avg_p = np.mean([i["p_split"] for i in tt_indiv if i["p_split"] >= 0])
        print(f"\n  {tt} (expected → {'STOP' if tt == 'coalescent' else 'SPLIT'}):")
        print(f"    Decision accuracy: {acc:.1%} (95% CI: [{ci_lo:.1%}, {ci_hi:.1%}])")
        print(f"    Correct/Total: {n_correct}/{n_total}")
        print(f"    Tree completeness: {comp_rate:.1%}")
        print(f"    Avg p_split: {avg_p:.3f}")

    # Overall
    all_n_correct = sum(s["correct"] for s in summary)
    all_n_total = sum(s["n_seeds"] for s in summary)
    overall_acc = all_n_correct / all_n_total if all_n_total else 0
    overall_ci_lo, overall_ci_hi = wilson_ci(all_n_correct, all_n_total)
    all_comps = [i["complete"] for i in all_individual]
    overall_comp = sum(all_comps) / len(all_comps) if all_comps else 0

    print(f"\n  OVERALL:")
    print(f"    Decision accuracy: {overall_acc:.1%} (95% CI: [{overall_ci_lo:.1%}, {overall_ci_hi:.1%}])")
    print(f"    Correct/Total: {all_n_correct}/{all_n_total}")
    print(f"    Tree completeness: {overall_comp:.1%}")

    # Pass/Fail
    coalescent_acc = sum(s["correct"] for s in summary if s["tree_type"] == "coalescent")
    coalescent_total = sum(s["n_seeds"] for s in summary if s["tree_type"] == "coalescent")
    structured_acc = sum(s["correct"] for s in summary if s["tree_type"] == "structured")
    structured_total = sum(s["n_seeds"] for s in summary if s["tree_type"] == "structured")

    c_acc = coalescent_acc / coalescent_total if coalescent_total else 0
    s_acc = structured_acc / structured_total if structured_total else 0

    print(f"\n{'=' * 70}")
    if c_acc >= 0.8 and s_acc >= 0.8:
        print(f"PASS: Both coalescent→STOP ({c_acc:.1%}) and structured→SPLIT ({s_acc:.1%})")
    else:
        print(f"PARTIAL: coalescent={c_acc:.1%}, structured={s_acc:.1%}")
    print(f"{'=' * 70}")

    # Save results
    output = {
        "total_scenarios": all_n_total,
        "overall_accuracy": overall_acc,
        "overall_ci_lower": overall_ci_lo,
        "overall_ci_upper": overall_ci_hi,
        "tree_completeness": overall_comp,
        "coalescent_accuracy": c_acc,
        "structured_accuracy": s_acc,
        "individual_results": all_individual,
        "per_group_summary": summary,
    }
    with open(OUT_JSON, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved to {OUT_JSON}")
    print(f"Failures recorded: {len(failures)}")

    return output


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--indelible-dir", default=None,
                    help="directory of pre-generated INDELible datasets")
    ap.add_argument("--model", default=None,
                    help="path to a boundary RF .pkl (default: auto-resolve)")
    ap.add_argument("--out", default="e2e_extended_results.json")
    a = ap.parse_args()
    INDELIBLE_DIR = a.indelible_dir
    OUT_JSON = a.out
    if a.model:
        import fusang_mhl.ml_split as _ms
        _ms._V4b_PATH = _ms._V3_PATH = _ms._V2_PATH = _ms._MODEL_PATH = \
            pathlib.Path(a.model)
        _ms._model = None
    run_extended_e2e()
