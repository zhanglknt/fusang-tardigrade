"""
Continue P0-2 E2E test — run structured tree scenarios only.
Loads existing coalescent results and merges.
"""
import sys
import os
import time
import json
import numpy as np
import math
from scipy.stats import norm

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

# Same config as test_mhl_e2e_extended.py
SCALES = [50, 100, 200, 300, 500]
N_SEEDS_PER = 3
PARAM_COMBOS = [
    (0.05, 0.02, 1000, "default"),
    (0.05, 0.01, 1000, "low_indel"),
    (0.05, 0.005, 1000, "vlow_indel"),
    (0.03, 0.02, 500, "short_seq"),
    (0.03, 0.01, 500, "short_lowindel"),
]


def wilson_ci(n_success, n_total, alpha=0.05):
    if n_total == 0:
        return (0.0, 0.0)
    z = norm.ppf(1 - alpha / 2)
    p = n_success / n_total
    denom = 1 + z**2 / n_total
    center = (p + z**2 / (2 * n_total)) / denom
    margin = z * math.sqrt((p * (1 - p) + z**2 / (4 * n_total)) / n_total) / denom
    return (max(0, center - margin), min(1, center + margin))


def generate_structured_data(n, L, sub, indel, seed):
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
    found = sum(1 for name in expected_names if name in nwk_str)
    return found == len(expected_names), found


def main():
    from fusang_mhl.level0_kmer import compute_l0_distance
    from fusang_mhl.ml_split import ml_split_decision, model_available
    from Bio.Phylo.TreeConstruction import DistanceMatrix as BioDM, DistanceTreeConstructor
    from Bio import Phylo
    from io import StringIO

    print("=" * 70)
    print("P0-2 Structured Tree E2E Test (continuation)")
    print("=" * 70)
    print(f"Model available: {model_available()}")
    expected_structured = len(SCALES) * N_SEEDS_PER * len(PARAM_COMBOS)
    print(f"Structured scenarios to run: {expected_structured}")
    print("=" * 70)

    structured_individual = []
    structured_summary = []

    tree_type = "structured"
    expected = "SPLIT"

    for sub, indel, L, label in PARAM_COMBOS:
        for n in SCALES:
            decisions_ok = []
            times = []
            complete_ok = []
            p_splits = []
            model_used_count = 0

            for seed_i in range(N_SEEDS_PER):
                seed = seed_i * 100 + 1  # same seeding as original (+1 for structured)
                try:
                    sequences, taxon_names = generate_structured_data(n, L, sub, indel, seed)
                    t0 = time.time()
                    D = compute_l0_distance(sequences, taxon_names)
                    decision = ml_split_decision(D, sequences, taxon_names, verbose=False)
                    elapsed = time.time() - t0
                    times.append(elapsed)
                    p_splits.append(decision.get('p_split', -1))
                    if decision.get('model_used', False):
                        model_used_count += 1
                    correct = decision['should_split'] == True
                    decisions_ok.append(correct)

                    # NJ tree completeness check
                    if hasattr(D, 'tolist'):
                        D_list = D.tolist()
                    else:
                        D_list = D
                    lower_tri = [[D_list[i][j] for j in range(i+1)] for i in range(len(D_list))]
                    dm = BioDM(list(taxon_names), lower_tri)
                    tree_obj = DistanceTreeConstructor().nj(dm)
                    buf = StringIO()
                    Phylo.write(tree_obj, buf, 'newick')
                    nwk = buf.getvalue()
                    complete, n_found = check_tree_completeness(nwk, taxon_names)
                    complete_ok.append(complete)

                    structured_individual.append({
                        "tree_type": tree_type, "expected": expected,
                        "n": n, "param": label, "seed": seed_i,
                        "correct": correct, "complete": complete,
                        "p_split": decision.get('p_split', -1),
                        "model_used": decision.get('model_used', False),
                        "time_s": elapsed,
                    })

                except Exception as e:
                    print(f"  ERROR: n={n} {label} seed={seed_i}: {e}", flush=True)

            n_ok = sum(decisions_ok)
            n_total = len(decisions_ok)
            acc = n_ok / n_total if n_total > 0 else 0
            ci_lo, ci_hi = wilson_ci(n_ok, n_total)
            avg_p = np.mean(p_splits) if p_splits else -1
            avg_t = np.mean(times) if times else 0
            comp = sum(complete_ok) / len(complete_ok) if complete_ok else 0

            structured_summary.append({
                "tree_type": tree_type, "expected": expected,
                "n": n, "param": label,
                "n_seeds": n_total, "correct": n_ok,
                "accuracy": acc, "ci_lower": ci_lo, "ci_upper": ci_hi,
                "avg_p_split": avg_p, "model_used": model_used_count,
                "tree_completeness": comp, "avg_time_s": avg_t,
            })

            ci_str = f"[{ci_lo:.1%}, {ci_hi:.1%}]"
            print(f"  n={n:>4} {label:<18s}: {n_ok}/{n_total} correct "
                  f"(acc={acc:.0%}, CI={ci_str}), avg_p={avg_p:.3f}, t={avg_t:.2f}s", flush=True)

    # Load existing coalescent results
    print("\nLoading existing coalescent results...")
    existing_file = "e2e_extended_results.json"
    if os.path.exists(existing_file):
        with open(existing_file) as f:
            existing = json.load(f)
        coal_individual = [x for x in existing.get("individual_results", []) if x["tree_type"] == "coalescent"]
        coal_summary = [x for x in existing.get("per_group_summary", []) if x["tree_type"] == "coalescent"]
        print(f"  Loaded {len(coal_individual)} coalescent scenarios")
    else:
        coal_individual = []
        coal_summary = []
        print("  No existing coalescent results found")

    # Merge
    all_individual = coal_individual + structured_individual
    all_summary = coal_summary + structured_summary

    # Compute overall stats
    all_n_correct = sum(s["correct"] for s in all_summary)
    all_n_total = sum(s["n_seeds"] for s in all_summary)
    overall_acc = all_n_correct / all_n_total if all_n_total else 0
    overall_ci_lo, overall_ci_hi = wilson_ci(all_n_correct, all_n_total)
    all_comps = [i["complete"] for i in all_individual]
    overall_comp = sum(all_comps) / len(all_comps) if all_comps else 0

    c_correct = sum(s["correct"] for s in all_summary if s["tree_type"] == "coalescent")
    c_total = sum(s["n_seeds"] for s in all_summary if s["tree_type"] == "coalescent")
    s_correct = sum(s["correct"] for s in all_summary if s["tree_type"] == "structured")
    s_total = sum(s["n_seeds"] for s in all_summary if s["tree_type"] == "structured")
    c_acc = c_correct / c_total if c_total else 0
    s_acc = s_correct / s_total if s_total else 0

    print(f"\n{'=' * 70}")
    print("COMBINED SUMMARY (coalescent + structured)")
    print(f"{'=' * 70}")
    print(f"  Coalescent (→STOP): {c_correct}/{c_total} correct, acc={c_acc:.1%}")
    coal_ci = wilson_ci(c_correct, c_total)
    print(f"    95% CI: [{coal_ci[0]:.1%}, {coal_ci[1]:.1%}]")
    print(f"  Structured (→SPLIT): {s_correct}/{s_total} correct, acc={s_acc:.1%}")
    struct_ci = wilson_ci(s_correct, s_total)
    print(f"    95% CI: [{struct_ci[0]:.1%}, {struct_ci[1]:.1%}]")
    print(f"\n  OVERALL: {all_n_correct}/{all_n_total} correct, acc={overall_acc:.1%}")
    print(f"    95% CI: [{overall_ci_lo:.1%}, {overall_ci_hi:.1%}]")
    print(f"    Tree completeness: {overall_comp:.1%}")

    print(f"\n{'=' * 70}")
    if c_acc >= 0.8 and s_acc >= 0.8:
        print(f"PASS: coalescent={c_acc:.1%}, structured={s_acc:.1%}")
    else:
        print(f"PARTIAL: coalescent={c_acc:.1%}, structured={s_acc:.1%}")
    print(f"{'=' * 70}")

    # Save merged results
    output = {
        "total_scenarios": all_n_total,
        "overall_accuracy": overall_acc,
        "overall_ci_lower": overall_ci_lo,
        "overall_ci_upper": overall_ci_hi,
        "tree_completeness": overall_comp,
        "coalescent_accuracy": c_acc,
        "coalescent_ci_lower": coal_ci[0],
        "coalescent_ci_upper": coal_ci[1],
        "structured_accuracy": s_acc,
        "structured_ci_lower": struct_ci[0],
        "structured_ci_upper": struct_ci[1],
        "individual_results": all_individual,
        "per_group_summary": all_summary,
    }
    with open("e2e_extended_results.json", "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nMerged results saved to e2e_extended_results.json")
    return output


if __name__ == "__main__":
    main()
