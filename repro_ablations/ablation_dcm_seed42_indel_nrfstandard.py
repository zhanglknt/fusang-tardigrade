#!/usr/bin/env python3
"""
Rescore the seed-42 DCM ablation chain (RESULTS, "A simplified pipeline
achieves state-of-the-art accuracy") with the STANDARD unrooted RF
(nrf_standard.py), replacing the superseded calc_nrf_simple numbers.

Dataset: the dataset named in the manuscript text -- n=200, indel rate=0.02,
seed 42 (indelible_n200_indel-independent files seed42_indel.fasta /
seed42_indel_true.nwk at repo root).  NOTE: the old raw trace
(diag_dcm_loss_n200_s42.json) was computed on test_n200.fasta (clean,
substitution-only) even though the manuscript text says indel=0.02; this
rescore uses the dataset the text describes.

Chain (identical construction to ablation_dcm.py, only the scorer changes):
  1. simplified            k-mer -> cosine -> NJ
  2. swap tree builder     k-mer -> cosine -> FastME BIONJ
  3. swap dist-matrix impl kmer_distance.compute_kmer_distance_matrix -> NJ
  4. +TF-IDF weighting     k-mer -> TF-IDF -> cosine -> NJ
  5. full DCM + EPA grafting (divide_and_conquer, NJ subtrees, default mode)
     swept over max_group_size in {10,20,25,30,40,50,60,80,100,120,150};
     the mg=10 configuration is the headline step-5 value.
Calibration: FastTree2 on the MAFFT alignment (seed42_indel_ft2.nwk).

Output: repro_ablations/ablation_dcm_seed42_indel_nrfstandard.json
"""
import os
import sys
import json
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

import nrf_standard as NS                                     # noqa: E402
import fusang_v2 as F                                         # noqa: E402
from kmer_distance import seq_to_kmer_freq, make_gap_pattern, compute_kmer_distance_matrix  # noqa: E402
from ablation_dcm import cosine_distance_matrix, apply_weighting, build_fastme  # noqa: E402

K = 5
GAP_STYLE = "gap2"
FASTA = os.path.join(ROOT, "seed42_indel.fasta")
TRUE_TREE = os.path.join(ROOT, "seed42_indel_true.nwk")
FT2_TREE = os.path.join(ROOT, "seed42_indel_ft2.nwk")


def main():
    t0 = time.time()
    names, seqs = F.read_fasta(FASTA)
    with open(TRUE_TREE, encoding="utf-8") as fh:
        true_nwk = fh.read()
    n = len(names)
    max_rf = 2 * (n - 3)
    print("loaded %d sequences from %s" % (n, FASTA))

    gap_pattern = make_gap_pattern(K, style=GAP_STYLE)
    out = {"scorer": "nrf_standard.py (root-independent unrooted RF, max_rf=2(n-3))",
           "dataset": "seed42_indel.fasta (n=200, indel=0.02, seed 42)",
           "n_taxa": n, "k": K, "gap_style": GAP_STYLE,
           "gap_pattern": list(gap_pattern)}

    def score(tree):
        rf = NS.rf_distance(tree.write(format=1), true_nwk)
        return {"rf": rf, "nrf": rf / max_rf}

    freq = np.array(
        [seq_to_kmer_freq(s, K, normalize=True, gap_pattern=gap_pattern)
         for s in seqs], dtype=np.float32)

    d1 = cosine_distance_matrix(freq)
    out["step1_simplified_nj"] = score(F.nj_tree(d1, names))
    print("step1 simplified NJ:        nRF=%.4f" % out["step1_simplified_nj"]["nrf"])

    try:
        out["step2_fastme"] = score(build_fastme(d1, names))
        print("step2 FastME:               nRF=%.4f" % out["step2_fastme"]["nrf"])
    except Exception as exc:
        out["step2_fastme"] = None
        print("step2 FastME FAILED: %s" % exc)

    d3 = compute_kmer_distance_matrix(seqs, names, k=K, metric="cosine",
                                      n_threads=4, gap_pattern=gap_pattern)
    out["step3_v2_distmat_nj"] = score(F.nj_tree(d3, names))
    out["dist_mat_max_diff"] = float(np.abs(d1 - d3).max())
    print("step3 v2 distmat NJ:        nRF=%.4f" % out["step3_v2_distmat_nj"]["nrf"])

    d4 = cosine_distance_matrix(apply_weighting(freq, "tfidf"))
    out["step4_tfidf_nj"] = score(F.nj_tree(d4, names))
    print("step4 TF-IDF NJ:            nRF=%.4f" % out["step4_tfidf_nj"]["nrf"])

    sweep = {}
    for mg in (10, 20, 25, 30, 40, 50, 60, 80, 100, 120, 150):
        key = "mg%d_nj_default" % mg
        try:
            tree = F.divide_and_conquer(
                seqs, names, max_group_size=mg, overlap=0.15, mode="default",
                distance_method="kmer", k=K, n_threads=4,
                gap_pattern=gap_pattern, _depth=0, tree_method="nj")
            sweep[key] = score(tree)
            print("step5 %-16s        nRF=%.4f" % (key, sweep[key]["nrf"]))
        except Exception as exc:
            sweep[key] = None
            print("step5 %-16s FAILED: %s" % (key, exc))
    out["step5_dcm_sweep_nj_default"] = sweep
    vals = [v["nrf"] for v in sweep.values() if v]
    out["step5_dcm_range"] = {"min": min(vals), "max": max(vals),
                              "mean": float(np.mean(vals)),
                              "sd": float(np.std(vals, ddof=1)),
                              "n_configs": len(vals)}

    if os.path.exists(FT2_TREE):
        with open(FT2_TREE, encoding="utf-8") as fh:
            ft2_nwk = fh.read()
        rf = NS.rf_distance(ft2_nwk, true_nwk)
        out["fasttree2_mafft"] = {"rf": rf, "nrf": rf / max_rf}
        print("FastTree2 (MAFFT):          nRF=%.4f" % out["fasttree2_mafft"]["nrf"])

    out["elapsed_s"] = time.time() - t0
    path = os.path.join(HERE, "ablation_dcm_seed42_indel_nrfstandard.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2, ensure_ascii=False)
    print("written: %s (%.0fs)" % (path, out["elapsed_s"]))


if __name__ == "__main__":
    main()
