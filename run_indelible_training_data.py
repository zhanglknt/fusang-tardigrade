#!/usr/bin/env python3
"""
Build the boundary-classifier training set from INDELible-simulated data.

This is the INDELible replacement for ``generate_v4b_data.py``.  Everything
downstream (feature extraction, labelling, sub-sampling) is identical -- only
the sequence simulator changes.

Usage:
    python run_indelible_training_data.py \
        --pool indelible_bc/pool \
        --out fusang_mhl/models/training_data_indelible.pkl
"""

import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fusang_mhl.training_data_generator import (  # noqa: E402
    generate_training_data, DEFAULT_SIM_CONFIGS,
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", default="indelible_bc/pool")
    ap.add_argument("--out",
                    default="fusang_mhl/models/training_data_indelible.pkl")
    ap.add_argument("--configs", default="all",
                    help="'all', or comma-separated config indices (for batching)")
    ap.add_argument("--max-sub", type=int, default=15,
                    help="cap on sub-cluster samples per tree (default 15)")
    ap.add_argument("--verbose", action="store_true", default=True)
    args = ap.parse_args()

    idx = None
    configs = DEFAULT_SIM_CONFIGS
    if args.configs != "all":
        idx = [int(x) for x in args.configs.split(",")]
        configs = [DEFAULT_SIM_CONFIGS[i] for i in idx]

    t0 = time.time()
    total_trees = sum(c[4] for c in configs)
    print("=" * 60)
    print("INDELible training data generation")
    print("=" * 60)
    print("configs: %d   trees: %d" % (len(configs), total_trees))
    print("pool:    %s" % os.path.abspath(args.pool))
    print("out:     %s" % os.path.abspath(args.out))
    print()

    failures = []
    data = generate_training_data(
        configs=configs,
        output_pkl=args.out,
        n_samples_target=50000,
        verbose=args.verbose,
        structured_ratio=0.5,
        indelible_dir=args.pool,
        failures=failures,
        config_indices=idx if args.configs != "all" else None,
        max_sub=args.max_sub,
    )

    from collections import Counter
    labels = [s["label"] for s in data]
    n_pos = sum(labels)
    print()
    print("=" * 60)
    print("elapsed: %.0fs (%.1f min)" % (time.time() - t0,
                                         (time.time() - t0) / 60))
    print("samples: %d  SPLIT=%d (%.1f%%)  STOP=%d"
          % (len(data), n_pos, 100.0 * n_pos / max(len(data), 1),
             len(data) - n_pos))
    print("types:   %s" % dict(Counter(
        s.get("meta", {}).get("sample_type", "?") for s in data)))
    print("failures: %d (of %d trees)" % (len(failures), total_trees))

    if failures:
        fpath = os.path.splitext(args.out)[0] + "_failures.json"
        with open(fpath, "w") as fh:
            json.dump(failures, fh, indent=2)
        print("failure log: %s" % fpath)
        for f in failures[:10]:
            print("   ci=%s ti=%s n=%s: %s" % (f["ci"], f["ti"], f["n"],
                                               f["error"]))


if __name__ == "__main__":
    main()
