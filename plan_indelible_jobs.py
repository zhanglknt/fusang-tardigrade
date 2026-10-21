#!/usr/bin/env python3
"""
Emit the INDELible job list that reproduces the boundary-classifier data.

The tree/clade decomposition is copied verbatim from
``fusang_mhl/training_data_generator.generate_training_data()`` /
``generate_structured_tree()`` and ``test_mhl_e2e_extended.py`` so that the
resulting dataset has exactly the same design (same n, L, clade sizes, same
number of trees) as the old JC69 one -- only the sequence simulator changes.

Usage:
    python plan_indelible_jobs.py --out indelible_bc/jobs.json
"""

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fusang_mhl.training_data_generator import DEFAULT_SIM_CONFIGS
from fusang_mhl.indelible_adapter import tag_for, tag_for_e2e

# --- E2E scenario grid (identical to test_mhl_e2e_extended.py) --------------
SCALES = [50, 100, 200, 300, 500]
N_SEEDS_PER = 3
PARAM_COMBOS = [
    (0.05, 0.02, 1000, "default"),
    (0.05, 0.01, 1000, "low_indel"),
    (0.05, 0.005, 1000, "vlow_indel"),
    (0.03, 0.02, 500, "short_seq"),
    (0.03, 0.01, 500, "short_lowindel"),
]

# Old simulator: expected substitutions/site/branch = sub_rate * bl, bl mean 0.1
SUB_BL = 0.1


def _job(jobs, tag, seed, n, L, sub, indel, meta):
    jobs.append({
        "tag": tag,
        "seed": int(seed),
        "n": int(n),
        "L": int(L),
        "sub": float(sub),
        "indel": float(indel),
        "sub_scale": float(sub) * SUB_BL,
        "meta": meta,
    })


def plan_training(jobs):
    """585 training trees -- same design as generate_training_data()."""
    for ci, (n, L, sub, indel, n_trees) in enumerate(DEFAULT_SIM_CONFIGS):
        for ti in range(n_trees):
            seed = ci * 1000 + ti
            rng_local = np.random.RandomState(seed)
            is_structured = rng_local.uniform() < 0.5

            if is_structured and n >= 8:
                n_clades = rng_local.randint(2, 4)
                rng = np.random.RandomState(seed)  # generate_structured_tree()
                base = max(4, n // n_clades)
                sizes = [base] * n_clades
                rem = n - base * n_clades
                for _ in range(rem):
                    sizes[rng.randint(0, n_clades - 1)] += 1
                for k, size in enumerate(sizes):
                    clade_sub = max(0.001, sub * rng.uniform(0.4, 2.5))
                    clade_indel = max(0.0, indel * rng.uniform(0.4, 2.0))
                    _job(jobs, tag_for(ci, ti, k), seed + k * 10000 + 7,
                         size, L, clade_sub, clade_indel,
                         {"split": "train", "ci": ci, "ti": ti,
                          "clade": k, "n_clades": int(n_clades),
                          "structured": True, "parent_n": n})
            else:
                _job(jobs, tag_for(ci, ti), seed, n, L, sub, indel,
                     {"split": "train", "ci": ci, "ti": ti,
                      "clade": None, "structured": False, "parent_n": n})


def plan_e2e(jobs):
    """150 E2E scenarios (75 coalescent + 75 structured), no silent skips."""
    # --- coalescent (expected STOP) ---
    for sub, indel, L, label in PARAM_COMBOS:
        for n in SCALES:
            for seed_i in range(N_SEEDS_PER):
                seed = seed_i * 100 + 0
                _job(jobs, tag_for_e2e("coalescent", n, label, seed_i),
                     seed, n, L, sub, indel,
                     {"split": "e2e", "tree_type": "coalescent", "n": n,
                      "param": label, "seed_i": seed_i, "clade": None})

    # --- structured (expected SPLIT) ---
    for sub, indel, L, label in PARAM_COMBOS:
        for n in SCALES:
            for seed_i in range(N_SEEDS_PER):
                seed = seed_i * 100 + 1
                np.random.seed(seed)  # generate_structured_data() uses global RNG
                n_clades = 2 if n <= 200 else 3
                clade_sizes = []
                remaining = n
                for c in range(n_clades - 1):
                    sz = max(30, remaining // (n_clades - c)
                             + np.random.randint(-10, 11))
                    sz = min(sz, remaining - 30 * (n_clades - c - 1))
                    clade_sizes.append(sz)
                    remaining -= sz
                clade_sizes.append(remaining)
                for ci, sz in enumerate(clade_sizes):
                    sub_rate = sub * (0.5 + ci * 0.75)
                    _job(jobs, tag_for_e2e("structured", n, label, seed_i, ci),
                         seed + ci * 1000, sz, L, sub_rate, indel,
                         {"split": "e2e", "tree_type": "structured", "n": n,
                          "param": label, "seed_i": seed_i,
                          "clade": ci, "n_clades": n_clades})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="indelible_bc/jobs.json")
    args = ap.parse_args()

    jobs = []
    plan_training(jobs)
    plan_e2e(jobs)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(jobs, fh)

    from collections import Counter
    c = Counter(j["meta"]["split"] for j in jobs)
    print("jobs written: %s -> %s" % (dict(c), args.out))
    print("total: %d" % len(jobs))
    n_trees = sum(1 for j in jobs
                  if j["meta"]["split"] == "train" and j["meta"]["clade"] in (0, None))
    print("training trees: %d" % n_trees)


if __name__ == "__main__":
    main()
