#!/usr/bin/env python3
"""
Rescore the Skmer probe NJ trees with the STANDARD unrooted RF
(nrf_standard.py), replacing the superseded calc_nrf_simple numbers.

Inputs (already on disk, produced by the 2026-10-22 probe):
  _skmer_probe/seed100_skmer_k31.npy / seed100_skmer_k21.npy  (200x200 Skmer
      distance matrices)
  _skmer_probe/seed100_skmer_taxa.json                        (taxon order)
  ../indelible_n200_indel/seed100_indel_true.nwk              (reference tree)
  ../indelible_n200_indel/seed100_indel_fusang.nwk            (Fusang CLI tree)

NJ trees are rebuilt with fusang_v2.nj_tree (the same builder used by the
original probe) and scored against the reference with nrf_standard.
"""
import os
import sys
import json

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

import nrf_standard as NS   # noqa: E402
import fusang_v2 as F       # noqa: E402

PROBE = os.path.join(HERE, "_skmer_probe")
DATA = os.path.join(ROOT, "indelible_n200_indel")


def main():
    with open(os.path.join(PROBE, "seed100_skmer_taxa.json"), encoding="utf-8") as fh:
        taxa = json.load(fh)
    if isinstance(taxa, dict):
        taxa = taxa.get("taxa") or taxa.get("names") or list(taxa.values())[0]
    with open(os.path.join(DATA, "seed100_indel_true.nwk"), encoding="utf-8") as fh:
        true_nwk = fh.read()

    out = {"scorer": "nrf_standard.py (root-independent unrooted RF, max_rf=2(n-3))"}
    for k in (31, 21):
        mat = np.load(os.path.join(PROBE, "seed100_skmer_k%d.npy" % k))
        assert mat.shape == (len(taxa), len(taxa)), (mat.shape, len(taxa))
        tree = F.nj_tree(mat, list(taxa))
        rf = NS.rf_distance(tree.write(format=1), true_nwk)
        n = len(taxa)
        out["skmer_k%d" % k] = {"rf": rf, "nrf": rf / (2 * (n - 3))}
        print("Skmer k=%d: RF=%d  nRF=%.4f" % (k, rf, rf / (2 * (n - 3))))

    fus = os.path.join(DATA, "seed100_indel_fusang.nwk")
    if os.path.exists(fus):
        with open(fus, encoding="utf-8") as fh:
            fus_nwk = fh.read()
        rf = NS.rf_distance(fus_nwk, true_nwk)
        n = len(NS.leaves(NS.parse_newick(fus_nwk)))
        out["fusang_seed100"] = {"rf": rf, "nrf": rf / (2 * (n - 3))}
        print("Fusang seed100: RF=%d  nRF=%.4f" % (rf, rf / (2 * (n - 3))))

    path = os.path.join(HERE, "skmer_probe_nrfstandard.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2, ensure_ascii=False)
    print("written: %s" % path)


if __name__ == "__main__":
    main()
