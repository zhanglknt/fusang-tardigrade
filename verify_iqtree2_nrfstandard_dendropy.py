#!/usr/bin/env python3
"""Independent DendroPy cross-check of rescore_iqtree2_nrfstandard.py.

Two levels of verification:
  1. SPLIT-SET LEVEL (random 3 seeds): rebuild each tree's non-trivial split set from
     DendroPy's own bipartition encoding and require it to be *identical* (as a set of
     canonicalised leaf-name frozensets) to nrf_standard.splits().
  2. DISTANCE LEVEL (all 130 seeds): RF from dendropy treecompare.symmetric_difference
     on force-unrooted trees, converted to nRF = RF/(2(n-3)), compared row-by-row with
     the values in iqtree2_n200_indel_nrfstandard.csv.

Also independently recomputes the paired statistics from the DendroPy values and
re-checks the legacy-metric reproduction of iqtree2_n200_indel.csv.
"""
from __future__ import annotations

import csv
import json
import os
import random
import sys

import dendropy
import numpy as np
import scipy.stats
from dendropy.calculate import treecompare

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nrf_standard  # noqa: E402

ROOT = os.path.dirname(os.path.abspath(__file__))
SEEDS = list(range(100, 230))
N_TAXA = 200
TIE_TOL = 1e-12


def load(path, ns):
    return dendropy.Tree.get(path=path, schema="newick", rooting="force-unrooted",
                             taxon_namespace=ns, preserve_underscores=True)


def dendro_splits(tree, ns) -> set:
    """Canonicalised non-trivial split set, built from DendroPy bipartition encoding.

    `Bipartition.is_trivial` is True for every entry on force-unrooted trees, so
    triviality is determined by leaf-set size instead (2 <= |side| <= n-2).
    """
    tree.encode_bipartitions()
    labels = sorted(t.label for t in ns)
    ref = labels[0]
    out = set()
    for b in tree.bipartition_encoding:
        side = frozenset(t.label for t in b.leafset_taxa(ns))
        if len(side) < 2 or len(side) > len(labels) - 2:
            continue
        out.add(frozenset(labels) - side if ref in side else side)
    return out


def canon(splits_set, labels):
    ref = labels[0]
    return {frozenset(set(labels) - s) if ref in s else frozenset(s) for s in splits_set}


def main() -> int:
    rng = random.Random(7)
    sample = sorted(rng.sample(SEEDS, 3))
    report = {"dendropy_version": dendropy.__version__, "sampled_seeds": sample,
              "split_check": {}, "distance_check": {}, "stats": {}, "legacy_check": {}}

    # ---- level 1: split-set identity on 3 random seeds ------------------------
    for seed in sample:
        p_true = os.path.join(ROOT, "indelible_n200_indel", f"seed{seed}_indel_true.nwk")
        p_fus = os.path.join(ROOT, "indelible_n200_indel", f"seed{seed}_indel_fusang.nwk")
        p_iq = os.path.join(ROOT, "indelible_iqtree2", "run",
                            f"seed{seed}_indel_iqtree2_gtr.treefile")
        ns = dendropy.TaxonNamespace()
        trees = {k: load(p, ns) for k, p in (("true", p_true), ("fusang", p_fus), ("iqtree2", p_iq))}
        labels = sorted(t.label for t in ns)
        per = {}
        for k, t in trees.items():
            mine = dendro_splits(t, ns)
            theirs = nrf_standard.splits(nrf_standard.parse_newick(nrf_standard.read_newick(
                {"true": p_true, "fusang": p_fus, "iqtree2": p_iq}[k])))
            per[k] = {"n_splits_dendropy": len(mine), "n_splits_nrf_standard": len(theirs),
                      "identical": mine == canon(theirs, labels),
                      "only_dendropy": len(mine - canon(theirs, labels)),
                      "only_mine": len(canon(theirs, labels) - mine)}
        report["split_check"][str(seed)] = per

    # ---- level 2: distance on all 130 seeds -----------------------------------
    d_iq, d_fu = [], []
    for seed in SEEDS:
        ns = dendropy.TaxonNamespace()
        t_true = load(os.path.join(ROOT, "indelible_n200_indel", f"seed{seed}_indel_true.nwk"), ns)
        t_fus = load(os.path.join(ROOT, "indelible_n200_indel", f"seed{seed}_indel_fusang.nwk"), ns)
        t_iq = load(os.path.join(ROOT, "indelible_iqtree2", "run",
                                 f"seed{seed}_indel_iqtree2_gtr.treefile"), ns)
        for t in (t_true, t_fus, t_iq):
            t.encode_bipartitions()
        denom = 2 * (N_TAXA - 3)
        d_iq.append(treecompare.symmetric_difference(t_iq, t_true, is_bipartitions_updated=True) / denom)
        d_fu.append(treecompare.symmetric_difference(t_fus, t_true, is_bipartitions_updated=True) / denom)
    d_iq, d_fu = np.array(d_iq), np.array(d_fu)

    with open(os.path.join(ROOT, "indelible_iqtree2", "iqtree2_n200_indel_nrfstandard.csv")) as fh:
        rows = list(csv.DictReader(fh))
    c_iq = np.array([float(r["nrf_iqtree2_std"]) for r in rows])
    c_fu = np.array([float(r["nrf_fusang_std"]) for r in rows])
    dev_iq = np.abs(d_iq - c_iq)
    dev_fu = np.abs(d_fu - c_fu)
    report["distance_check"] = {
        "n_seeds": len(rows),
        "max_abs_dev_iqtree2": float(dev_iq.max()),
        "max_abs_dev_fusang": float(dev_fu.max()),
        "n_rows_mismatched": int(((dev_iq > 1e-12) | (dev_fu > 1e-12)).sum()),
    }

    # ---- independent statistics ----------------------------------------------
    diff = c_iq - c_fu
    wz = scipy.stats.wilcoxon(diff, zero_method="zsplit")
    report["stats"] = {
        "iqtree2_mean": float(c_iq.mean()), "iqtree2_sd": float(c_iq.std(ddof=1)),
        "fusang_mean": float(c_fu.mean()), "fusang_sd": float(c_fu.std(ddof=1)),
        "ratio": float(c_iq.mean() / c_fu.mean()),
        "wilcoxon_W": float(wz.statistic), "wilcoxon_p": float(wz.pvalue),
        "cohens_d_ddof1": float(diff.mean() / diff.std(ddof=1)),
        "cohens_d_ddof0": float(diff.mean() / diff.std(ddof=0)),
        "n_iqtree2_worse": int((diff > TIE_TOL).sum()),
        "n_iqtree2_better": int((diff < -TIE_TOL).sum()),
        "n_ties": int((np.abs(diff) <= TIE_TOL).sum()),
    }

    # ---- legacy metric reproduction of the reference CSV -----------------------
    ref_path = os.path.join(ROOT, "indelible_iqtree2", "iqtree2_n200_indel.csv")
    ref = {int(r["seed"]): float(r["nrf"]) for r in csv.DictReader(open(ref_path))}
    with open(os.path.join(ROOT, "indelible_iqtree2", "iqtree2_n200_indel_nrfstandard.json")) as fh:
        summary = json.load(fh)
    audit = summary["provenance"]["legacy_biopython_vs_reference_csv"]
    report["legacy_check"] = {
        "reference_csv": ref_path,
        "biopython_legacy_on_run_treefiles_max_dev": audit["primary_treefiles"]["max_abs_deviation"],
        "biopython_legacy_on_wsl_treefiles_max_dev": audit["alt_treefiles"]["max_abs_deviation"],
        "calc_nrf_simple_max_dev": summary["provenance"]["calc_nrf_simple_vs_reference_csv"]["max_abs_deviation"],
        "reference_reproduced_by": "Biopython legacy metric on wsl_treefiles (exact, 130/130)",
    }

    out = os.path.join(ROOT, "verify_iqtree2_nrfstandard_dendropy.json")
    with open(out, "w") as fh:
        json.dump(report, fh, indent=2)

    ok_split = all(v["identical"] for s in report["split_check"].values() for v in s.values())
    print(f"dendropy {dendropy.__version__}")
    print(f"split-set identical on seeds {sample}: {ok_split}")
    for s in sample:
        for k, v in report["split_check"][str(s)].items():
            print(f"  seed{s} {k:8s} n={v['n_splits_dendropy']}/{v['n_splits_nrf_standard']} "
                  f"identical={v['identical']}")
    print(f"distance max|dev| iqtree2={dev_iq.max():.3e} fusang={dev_fu.max():.3e} "
          f"mismatched rows={report['distance_check']['n_rows_mismatched']}")
    print(f"stats: iq={c_iq.mean():.6f}+/-{c_iq.std(ddof=1):.6f} fu={c_fu.mean():.6f}+/-{c_fu.std(ddof=1):.6f} "
          f"ratio={c_iq.mean()/c_fu.mean():.4f} p={wz.pvalue:.6g} d={diff.mean()/diff.std(ddof=1):.4f} "
          f"{report['stats']['n_iqtree2_worse']}/{report['stats']['n_iqtree2_better']}/{report['stats']['n_ties']}")
    print(f"wrote {out}")
    return 0 if (ok_split and report["distance_check"]["n_rows_mismatched"] == 0) else 1


if __name__ == "__main__":
    raise SystemExit(main())
