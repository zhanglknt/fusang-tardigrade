#!/usr/bin/env python3
"""
Rescore the cosine-vs-JSD comparison (cosine_vs_jsd.py) with the STANDARD
unrooted Robinson-Foulds distance (nrf_standard.py, root-independent split
canonicalisation), replacing the superseded calc_nrf_simple.py numbers.

Trees are rebuilt deterministically from the same inputs
(indelible_n200_indel/seed{100..129}_indel.fasta, k=5, gap2 pattern
(0,1,2,5,6), fusang_v2.nj_tree, no NNI) exactly as in cosine_vs_jsd.py;
only the scorer changes.

Statistics: paired Wilcoxon (scipy, zero_method='zsplit') on INTEGER RF
counts (nRF x 394 = 2(n-3), n=200), paired Cohen's d = mean/sd of per-seed
nRF differences.  diff = sqrtJSD - cosine (positive => cosine better).

Output: repro_ablations/cosine_vs_jsd_n200_indel_30seeds_nrfstandard.json
"""
import os
import sys
import json

import numpy as np
from scipy import stats as st

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

import nrf_standard as NS                                     # noqa: E402
import fusang_v2 as F                                         # noqa: E402
from kmer_distance import seq_to_kmer_freq, make_gap_pattern  # noqa: E402
from cosine_vs_jsd import cosine_distance_matrix, jsd_distance_matrix  # noqa: E402

K = 5
GAP_STYLE = "gap2"
DATA_DIR = os.path.join(ROOT, "indelible_n200_indel")
SEEDS = list(range(100, 130))
MAX_RF = 2 * (200 - 3)  # 394


def main():
    gap_pattern = make_gap_pattern(K, style=GAP_STYLE)
    rows = []
    for seed in SEEDS:
        fa = os.path.join(DATA_DIR, "seed%d_indel.fasta" % seed)
        tr = os.path.join(DATA_DIR, "seed%d_indel_true.nwk" % seed)
        names, seqs = F.read_fasta(fa)
        with open(tr, encoding="utf-8") as fh:
            true_nwk = fh.read()
        freq = np.array(
            [seq_to_kmer_freq(s, K, normalize=True, gap_pattern=gap_pattern)
             for s in seqs], dtype=np.float32)
        d_cos = cosine_distance_matrix(freq)
        _, d_sqrt = jsd_distance_matrix(freq)
        t_cos = F.nj_tree(d_cos, names).write(format=1)
        t_sqrt = F.nj_tree(d_sqrt, names).write(format=1)
        rf_cos = NS.rf_distance(t_cos, true_nwk)
        rf_sqrt = NS.rf_distance(t_sqrt, true_nwk)
        rows.append({"seed": seed,
                     "rf_cosine": rf_cos, "rf_sqrtjsd": rf_sqrt,
                     "nrf_cosine": rf_cos / MAX_RF, "nrf_sqrtjsd": rf_sqrt / MAX_RF})
        print("seed %d: cosine RF=%d (nRF=%.4f)  sqrtJSD RF=%d (nRF=%.4f)"
              % (seed, rf_cos, rf_cos / MAX_RF, rf_sqrt, rf_sqrt / MAX_RF))

    cos = np.array([r["rf_cosine"] for r in rows], dtype=float)
    sqr = np.array([r["rf_sqrtjsd"] for r in rows], dtype=float)
    diff = sqr - cos  # positive => cosine better
    w = st.wilcoxon(cos, sqr, zero_method="zsplit")
    d = float((diff / MAX_RF).mean() / (diff / MAX_RF).std(ddof=1))
    out = {
        "scorer": "nrf_standard.py (root-independent unrooted RF, max_rf=2(n-3))",
        "wilcoxon": "scipy.stats.wilcoxon on integer RF counts, zero_method='zsplit'",
        "n_seeds": len(rows),
        "mean_nrf_cosine": float((cos / MAX_RF).mean()),
        "sd_nrf_cosine": float((cos / MAX_RF).std(ddof=1)),
        "mean_nrf_sqrtjsd": float((sqr / MAX_RF).mean()),
        "sd_nrf_sqrtjsd": float((sqr / MAX_RF).std(ddof=1)),
        "cosine_better_seeds": int((diff > 0).sum()),
        "sqrtjsd_better_seeds": int((diff < 0).sum()),
        "tie_seeds": int((diff == 0).sum()),
        "wilcoxon_statistic": float(w.statistic),
        "wilcoxon_p": float(w.pvalue),
        "cohens_d_paired": d,
        "per_seed": rows,
    }
    path = os.path.join(HERE, "cosine_vs_jsd_n200_indel_30seeds_nrfstandard.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2, ensure_ascii=False)
    print("\ncosine  %.4f +/- %.4f" % (out["mean_nrf_cosine"], out["sd_nrf_cosine"]))
    print("sqrtJSD %.4f +/- %.4f" % (out["mean_nrf_sqrtjsd"], out["sd_nrf_sqrtjsd"]))
    print("cosine better %d/%d, sqrtJSD better %d, ties %d"
          % (out["cosine_better_seeds"], len(rows),
             out["sqrtjsd_better_seeds"], out["tie_seeds"]))
    print("Wilcoxon p = %.4g   paired d = %.3f" % (out["wilcoxon_p"], d))
    print("written: %s" % path)


if __name__ == "__main__":
    main()
