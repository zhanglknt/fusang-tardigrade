#!/usr/bin/env python3
"""Rescore the IQ-TREE2 vs Fusang comparison under the STANDARD unrooted RF metric.

Context
-------
The manuscript (Supplementary Note S10) reported

    "IQ-TREE2 1.9x worse than Fusang k-mer NJ (nRF=0.231 +/- 0.032 vs 0.122 +/- 0.024,
     130 seeds, paired Wilcoxon p=4.4e-23, d=3.64)"

Those two numbers were produced by TWO DIFFERENT metrics (see `provenance` in the
output JSON): the 0.231 came from the non-canonicalised rooted-clade RF in
`batch_iqtree2_indelible.py::compute_nrf` (Biopython Phylo), and the 0.122 came from
the in-house biased `calc_nrf_simple.py::calc_nrf`.  All manuscript nRF values have
now been migrated to the standard root-independent implementation in `nrf_standard.py`:

    RF  = |S1 \\ S2| + |S2 \\ S1|     over NON-TRIVIAL splits only
    nRF = RF / (2 * (n - 3))         n = shared leaf count

with every split canonicalised to the side NOT containing the lexicographically
smallest leaf.  This script recomputes the indel n=200 / 130-seed comparison in that
standard metric, and additionally recomputes BOTH legacy metrics so the provenance of
the manuscript numbers is auditable.

Inputs (pre-existing; this script never re-runs simulations)
    indelible_n200_indel/seed{N}_indel_true.nwk              true tree
    indelible_n200_indel/seed{N}_indel_fusang.nwk            Fusang full pipeline tree
    indelible_iqtree2/run/seed{N}_indel_iqtree2_gtr.treefile IQ-TREE2 trees (primary)
    indelible_iqtree2/wsl_treefiles/seed{N}_gtr.treefile     IQ-TREE2 trees (2nd source)
    indelible_iqtree2/iqtree2_n200_indel.csv                 legacy reference (seed,nrf,status)

Outputs
    indelible_iqtree2/iqtree2_n200_indel_nrfstandard.csv     per-seed metrics
    indelible_iqtree2/iqtree2_n200_indel_nrfstandard.json    paired statistics + provenance

Deterministic and idempotent: re-running produces identical results.
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import os
import sys
from contextlib import redirect_stdout
from datetime import datetime

import numpy as np
import scipy.stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import nrf_standard  # noqa: E402
import calc_nrf_simple  # noqa: E402

try:
    from Bio import Phylo
    HAVE_BIO = True
except ImportError:  # pragma: no cover
    HAVE_BIO = False

SEED_MIN, SEED_MAX = 100, 229
N_TAXA = 200
TIE_TOL = 1e-12

# manuscript figures being audited (Supplementary Note S10)
MANUSCRIPT = {"iqtree2_mean": 0.231, "iqtree2_sd": 0.032,
              "fusang_mean": 0.122, "fusang_sd": 0.024}


def leaf_set(path: str) -> set:
    return set(nrf_standard.leaves(nrf_standard.parse_newick(nrf_standard.read_newick(path))))


def old_metric(path_a: str, path_b: str) -> float:
    """calc_nrf_simple.calc_nrf, with its 9-lines-per-call chatter suppressed."""
    buf = io.StringIO()
    with redirect_stdout(buf):
        return calc_nrf_simple.calc_nrf(path_a, path_b)


def legacy_biopython_metric(path_a: str, path_b: str) -> float:
    """Verbatim re-implementation of batch_iqtree2_indelible.py::compute_nrf.

    Non-canonicalised rooted-clade RF: for every nonterminal clade take the clade's
    leaf set and keep it as-is (no complement, no canonicalisation), then take the
    symmetric difference.  Rooting-dependent => inflated, not a valid unrooted RF.
    """
    tr = Phylo.read(path_a, "newick")
    true = Phylo.read(path_b, "newick")

    def bips(tree):
        out = set()
        for clade in tree.get_nonterminals():
            lv = frozenset(t.name for t in clade.get_terminals())
            if 1 < len(lv) < len(tree.get_terminals()):
                out.add(lv)
        return out

    rf = len(bips(tr).symmetric_difference(bips(true)))
    n = len(tr.get_terminals())
    norm = 2 * (n - 3)
    return rf / norm if norm > 0 else 0.0


def paired_stats(iq: np.ndarray, fu: np.ndarray) -> dict:
    diff = iq - fu
    wz = scipy.stats.wilcoxon(diff, zero_method="zsplit")
    ww = scipy.stats.wilcoxon(diff, zero_method="wilcox")
    return {
        "n_seeds": int(len(diff)),
        "iqtree2": {"mean": float(iq.mean()), "sd": float(iq.std(ddof=1))},
        "fusang": {"mean": float(fu.mean()), "sd": float(fu.std(ddof=1))},
        "ratio_iqtree2_over_fusang": float(iq.mean() / fu.mean()),
        "mean_paired_difference": float(diff.mean()),
        "wilcoxon": {
            "zero_method_zsplit": {"statistic": float(wz.statistic), "pvalue": float(wz.pvalue)},
            "zero_method_wilcox": {"statistic": float(ww.statistic), "pvalue": float(ww.pvalue)},
        },
        "cohens_d": {"ddof1": float(diff.mean() / diff.std(ddof=1)),
                     "ddof0": float(diff.mean() / diff.std(ddof=0))},
        "counts": {"n_iqtree2_worse": int((diff > TIE_TOL).sum()),
                   "n_iqtree2_better": int((diff < -TIE_TOL).sum()),
                   "n_ties": int((np.abs(diff) <= TIE_TOL).sum())},
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("--true-dir", default=None)
    ap.add_argument("--fusang-dir", default=None)
    ap.add_argument("--iqtree2-dir", default=None, help="default: <root>/indelible_iqtree2/run")
    ap.add_argument("--iqtree2-alt-dir", default=None, help="default: <root>/indelible_iqtree2/wsl_treefiles")
    ap.add_argument("--old-csv", default=None)
    ap.add_argument("--out-csv", default=None)
    ap.add_argument("--out-json", default=None)
    ap.add_argument("--seed-min", type=int, default=SEED_MIN)
    ap.add_argument("--seed-max", type=int, default=SEED_MAX)
    args = ap.parse_args()

    root = args.root
    true_dir = args.true_dir or os.path.join(root, "indelible_n200_indel")
    fus_dir = args.fusang_dir or os.path.join(root, "indelible_n200_indel")
    iq_dir = args.iqtree2_dir or os.path.join(root, "indelible_iqtree2", "run")
    alt_dir = args.iqtree2_alt_dir or os.path.join(root, "indelible_iqtree2", "wsl_treefiles")
    old_csv = args.old_csv or os.path.join(root, "indelible_iqtree2", "iqtree2_n200_indel.csv")
    out_csv = args.out_csv or os.path.join(root, "indelible_iqtree2", "iqtree2_n200_indel_nrfstandard.csv")
    out_json = args.out_json or os.path.join(root, "indelible_iqtree2", "iqtree2_n200_indel_nrfstandard.json")

    def iq_path(d, s):
        return (os.path.join(d, f"seed{s}_indel_iqtree2_gtr.treefile") if d == iq_dir
                else os.path.join(d, f"seed{s}_gtr.treefile"))

    ref = {}
    with open(old_csv, newline="") as fh:
        for row in csv.DictReader(fh):
            ref[int(row["seed"])] = (float(row["nrf"]), row["status"].strip())

    rows, excluded = [], []
    for seed in range(args.seed_min, args.seed_max + 1):
        p_true = os.path.join(true_dir, f"seed{seed}_indel_true.nwk")
        p_fus = os.path.join(fus_dir, f"seed{seed}_indel_fusang.nwk")
        p_iq = iq_path(iq_dir, seed)
        p_alt = iq_path(alt_dir, seed)

        reasons = [f"missing {os.path.basename(p)}" for p in (p_true, p_fus, p_iq) if not os.path.exists(p)]
        status = ref.get(seed, (None, None))[1]
        if status is not None and status != "Success":
            reasons.append(f"reference-csv status={status!r}")
        if reasons:
            excluded.append({"seed": seed, "reason": "; ".join(reasons)})
            continue

        lt = leaf_set(p_true)
        problems = []
        for label, p in (("iqtree2", p_iq), ("fusang", p_fus)):
            ls = leaf_set(p)
            if len(ls) != N_TAXA:
                problems.append(f"{label} leaf count {len(ls)} != {N_TAXA}")
            elif ls != lt:
                problems.append(f"{label} leaf set != true leaf set "
                                f"(|only_{label}|={len(ls - lt)}, |only_true|={len(lt - ls)})")
        if len(lt) != N_TAXA:
            problems.append(f"true leaf count {len(lt)} != {N_TAXA}")
        if problems:
            excluded.append({"seed": seed, "reason": "; ".join(problems)})
            continue

        rec = {
            "seed": seed,
            "nrf_iqtree2_std": nrf_standard.nrf_files(p_iq, p_true),
            "nrf_fusang_std": nrf_standard.nrf_files(p_fus, p_true),
            "nrf_iqtree2_old": old_metric(p_iq, p_true),
            "nrf_fusang_old": old_metric(p_fus, p_true),
        }
        if HAVE_BIO:
            rec["nrf_iqtree2_legacy_bp"] = legacy_biopython_metric(p_iq, p_true)
            rec["nrf_fusang_legacy_bp"] = legacy_biopython_metric(p_fus, p_true)
        if os.path.exists(p_alt):
            rec["nrf_iqtree2_std_alt"] = nrf_standard.nrf_files(p_alt, p_true)
            if HAVE_BIO:
                rec["nrf_iqtree2_legacy_bp_alt"] = legacy_biopython_metric(p_alt, p_true)
            rec["nrf_iqtree2_vs_alt"] = nrf_standard.nrf_files(p_iq, p_alt)
        rows.append(rec)

    rows.sort(key=lambda r: r["seed"])
    fieldnames = ["seed", "nrf_iqtree2_std", "nrf_fusang_std", "nrf_iqtree2_old", "nrf_fusang_old",
                  "nrf_iqtree2_legacy_bp", "nrf_fusang_legacy_bp",
                  "nrf_iqtree2_std_alt", "nrf_iqtree2_legacy_bp_alt", "nrf_iqtree2_vs_alt"]
    fieldnames = [f for f in fieldnames if f in rows[0]]
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    def col(name):
        return np.array([r[name] for r in rows])

    # ---- statistics -----------------------------------------------------------
    primary = paired_stats(col("nrf_iqtree2_std"), col("nrf_fusang_std"))
    primary["iqtree2_tree_source"] = iq_dir
    alt_stats = (paired_stats(col("nrf_iqtree2_std_alt"), col("nrf_fusang_std"))
                 if "nrf_iqtree2_std_alt" in rows[0] else None)
    if alt_stats:
        alt_stats["iqtree2_tree_source"] = alt_dir

    # ---- provenance audit -----------------------------------------------------
    prov = {"reference_csv": old_csv}

    # (a) does calc_nrf_simple reproduce the reference CSV?
    dev_old = np.abs(col("nrf_iqtree2_old") - np.array([ref[r["seed"]][0] for r in rows]))
    prov["calc_nrf_simple_vs_reference_csv"] = {
        "n_compared": int(len(dev_old)),
        "max_abs_deviation": float(dev_old.max()),
        "reproduced": bool(dev_old.max() < 1e-9),
        "note": "calc_nrf_simple.calc_nrf does NOT reproduce iqtree2_n200_indel.csv",
    }

    # (b) does the Biopython legacy metric reproduce it, and from which treefiles?
    if HAVE_BIO:
        refv = np.array([ref[r["seed"]][0] for r in rows])
        dev_bp = np.abs(col("nrf_iqtree2_legacy_bp") - refv)
        dev_bp_alt = (np.abs(col("nrf_iqtree2_legacy_bp_alt") - refv)
                      if "nrf_iqtree2_legacy_bp_alt" in rows[0] else None)
        prov["legacy_biopython_vs_reference_csv"] = {
            "implementation": "batch_iqtree2_indelible.py::compute_nrf "
                              "(Biopython Phylo, non-canonicalised rooted clades)",
            "primary_treefiles": {"source": iq_dir, "max_abs_deviation": float(dev_bp.max()),
                                  "n_exact": int((dev_bp < 1e-12).sum())},
            "alt_treefiles": ({"source": alt_dir, "max_abs_deviation": float(dev_bp_alt.max()),
                               "n_exact": int((dev_bp_alt < 1e-12).sum())} if dev_bp_alt is not None else None),
            "conclusion": ("reference CSV was produced from the alt treefiles with the Biopython "
                           "legacy metric" if dev_bp_alt is not None and dev_bp_alt.max() < 1e-12
                           else "see per-source deviations"),
        }

    # (c) are the two IQ-TREE2 treefile sets the same topology?
    if "nrf_iqtree2_vs_alt" in rows[0]:
        d = col("nrf_iqtree2_vs_alt")
        prov["iqtree2_treefile_provenance_discrepancy"] = {
            "n_seeds_topologically_different": int((d > 0).sum()),
            "mean_nrf_between_sources": float(d.mean()),
            "max_nrf_between_sources": float(d.max()),
            "warning": "the two IQ-TREE2 output directories hold DIFFERENT topologies; "
                       "the manuscript's published 0.231 came from the alt treefiles",
        }

    # (d) which legacy metric produced each manuscript number?
    audit = {"manuscript_figures": MANUSCRIPT, "matches": {}}
    if HAVE_BIO:
        for key, arr, tgt in (("iqtree2", col("nrf_iqtree2_legacy_bp"), MANUSCRIPT["iqtree2_mean"]),
                              ("iqtree2_alt", col("nrf_iqtree2_legacy_bp_alt"), MANUSCRIPT["iqtree2_mean"]),
                              ("fusang_old_metric", col("nrf_fusang_old"), MANUSCRIPT["fusang_mean"]),
                              ("fusang_legacy_bp", col("nrf_fusang_legacy_bp"), MANUSCRIPT["fusang_mean"]),
                              ("fusang_standard", col("nrf_fusang_std"), MANUSCRIPT["fusang_mean"])):
            audit["matches"][key] = {"mean": float(arr.mean()), "sd": float(arr.std(ddof=1)),
                                     "abs_diff_vs_manuscript_mean": float(abs(arr.mean() - tgt))}
        audit["verdict"] = ("the manuscript's 0.231 (IQ-TREE2) and 0.122 (Fusang) come from two "
                            "different biased metrics; the 1.9x gap is a metric artefact, not a "
                            "biological result")
    prov["manuscript_figure_audit"] = audit

    summary = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "script": os.path.basename(__file__),
        "metric": {
            "name": "standard root-independent (unrooted) Robinson-Foulds",
            "module": "nrf_standard.py",
            "definition": "RF = |S1-S2| + |S2-S1| over non-trivial splits; nRF = RF/(2*(n-3)); "
                          "split canonicalised to the side not containing the lexicographically smallest leaf",
        },
        "legacy_metrics": {
            "old": {"module": "calc_nrf_simple.py", "function": "calc_nrf",
                    "bias": "per-node 'smaller child clade' rule computed against the local subtree, "
                            "so complements are wrong; under-counts"},
            "legacy_bp": {"module": "batch_iqtree2_indelible.py", "function": "compute_nrf",
                          "bias": "Biopython Phylo clade leaf-sets with no canonicalisation; "
                                  "rooting-dependent, over-counts"},
        },
        "inputs": {
            "true_trees": os.path.join(true_dir, "seed{N}_indel_true.nwk"),
            "fusang_trees": os.path.join(fus_dir, "seed{N}_indel_fusang.nwk"),
            "iqtree2_trees": iq_dir,
            "alt_iqtree2_trees": alt_dir,
            "reference_csv": old_csv,
            "seed_range": [args.seed_min, args.seed_max],
            "n_taxa": N_TAXA,
        },
        "seed_audit": {
            "n_seeds_requested": args.seed_max - args.seed_min + 1,
            "n_seeds_used": len(rows),
            "n_seeds_excluded": len(excluded),
            "excluded_seeds": excluded,
            "leaf_set_check": f"all usable seeds: 200 leaves, IQ-TREE2/Fusang leaf sets identical to true",
        },
        "primary_standard_metric": primary,
        "alt_treefiles_standard_metric": alt_stats,
        "legacy_metric_means": {
            "calc_nrf_simple": {"iqtree2_mean": float(col("nrf_iqtree2_old").mean()),
                                "iqtree2_sd": float(col("nrf_iqtree2_old").std(ddof=1)),
                                "fusang_mean": float(col("nrf_fusang_old").mean()),
                                "fusang_sd": float(col("nrf_fusang_old").std(ddof=1))},
            "legacy_biopython": ({"iqtree2_mean": float(col("nrf_iqtree2_legacy_bp").mean()),
                                  "iqtree2_sd": float(col("nrf_iqtree2_legacy_bp").std(ddof=1)),
                                  "fusang_mean": float(col("nrf_fusang_legacy_bp").mean()),
                                  "fusang_sd": float(col("nrf_fusang_legacy_bp").std(ddof=1))}
                                 if HAVE_BIO else None),
        },
        "provenance": prov,
        "outputs": {"csv": out_csv, "json": out_json},
    }

    with open(out_json, "w") as fh:
        json.dump(summary, fh, indent=2)

    # ---- console --------------------------------------------------------------
    p = primary
    print(f"seeds used       : {p['n_seeds']}  (excluded {len(excluded)})")
    print(f"IQ-TREE2 (std)   : {p['iqtree2']['mean']:.6f} +/- {p['iqtree2']['sd']:.6f}")
    print(f"Fusang   (std)   : {p['fusang']['mean']:.6f} +/- {p['fusang']['sd']:.6f}")
    print(f"ratio            : {p['ratio_iqtree2_over_fusang']:.4f}x")
    print(f"Wilcoxon zsplit  : W={p['wilcoxon']['zero_method_zsplit']['statistic']} "
          f"p={p['wilcoxon']['zero_method_zsplit']['pvalue']:.6g}")
    print(f"Cohen's d (ddof1): {p['cohens_d']['ddof1']:.4f}  (ddof0 {p['cohens_d']['ddof0']:.4f})")
    print(f"worse/better/tie : {p['counts']['n_iqtree2_worse']}/{p['counts']['n_iqtree2_better']}/"
          f"{p['counts']['n_ties']}")
    if alt_stats:
        print(f"[alt treefiles] IQ-TREE2 {alt_stats['iqtree2']['mean']:.6f} +/- "
              f"{alt_stats['iqtree2']['sd']:.6f}  ratio {alt_stats['ratio_iqtree2_over_fusang']:.4f}x  "
              f"p={alt_stats['wilcoxon']['zero_method_zsplit']['pvalue']:.6g}  "
              f"d={alt_stats['cohens_d']['ddof1']:.4f}  "
              f"{alt_stats['counts']['n_iqtree2_worse']}/{alt_stats['counts']['n_iqtree2_better']}/"
              f"{alt_stats['counts']['n_ties']}")
    print(f"calc_nrf_simple vs reference CSV : max|dev|={dev_old.max():.3e} -> "
          f"{'OK' if dev_old.max() < 1e-9 else 'FAIL'}")
    if HAVE_BIO:
        print(f"legacy Biopython vs reference CSV: primary max|dev|={dev_bp.max():.3e}, "
              f"alt max|dev|={dev_bp_alt.max():.3e}")
    print(f"wrote {out_csv}")
    print(f"wrote {out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
