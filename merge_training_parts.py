#!/usr/bin/env python3
"""Merge the per-batch training-data pickles into one dataset file."""
import argparse
import glob
import os
import pickle
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pattern", default="indelible_bc/part_*.pkl")
    ap.add_argument("--out",
                    default="fusang_mhl/models/training_data_indelible.pkl")
    a = ap.parse_args()

    parts = sorted(glob.glob(a.pattern))
    if not parts:
        raise SystemExit("no parts matched %s" % a.pattern)

    data = []
    for p in parts:
        with open(p, "rb") as fh:
            chunk = pickle.load(fh)
        print("%-32s %5d" % (os.path.basename(p), len(chunk)))
        data.extend(chunk)

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "wb") as fh:
        pickle.dump(data, fh, protocol=pickle.HIGHEST_PROTOCOL)

    n_pos = sum(1 for s in data if s["label"] == 1)
    print()
    print("merged %d parts -> %d samples" % (len(parts), len(data)))
    print("SPLIT=%d (%.1f%%)  STOP=%d"
          % (n_pos, 100.0 * n_pos / len(data), len(data) - n_pos))
    print("types: %s" % dict(Counter(
        s.get("meta", {}).get("sample_type", "?") for s in data)))
    print("n_taxa: %s" % dict(sorted(Counter(
        s.get("meta", {}).get("n_taxa", 0) for s in data).items())))
    print("saved -> %s" % a.out)


if __name__ == "__main__":
    main()
