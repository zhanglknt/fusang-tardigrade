#!/usr/bin/env python3
"""Step 3c/4: ETE3 scoring (unrooted nRF, AFproject standard) + gap stats.
(a) each inferred tree vs Fischer 2013 whole-mitogenome reference tree
(b) Fusang multi-k tree vs MAFFT+FastTree2 tree (method agreement)"""
import os, sys, json
import numpy as np
from ete3 import Tree

SRC = r"D:\系统发育树项目\Fusang\Fusang-main"
BASE = os.path.join(SRC, "real_data", "afproject_genome", "gene_dna")
TREES = os.path.join(BASE, "trees")
ALN = os.path.join(BASE, "aln")
REF = os.path.join(SRC, "real_data", "swisstree", "afproject_repo",
                   "datasets", "genome", "fish_mito", "tree.newick")

GENES = ["ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
         "ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB"]

def ete_nrf(qpath, rpath):
    qt = Tree(qpath, format=1)
    rt = Tree(rpath, format=1)
    res = qt.compare(rt, unrooted=True)
    return float(res["norm_rf"])

def gap_stats(path):
    recs, cur = [], None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if cur is not None:
                recs.append(cur)
            cur = ""
        else:
            cur += line
    if cur is not None:
        recs.append(cur)
    L = len(recs[0])
    gap = sum(s.count("-") for s in recs)
    return L, gap / (L * len(recs))

extract = json.load(open(os.path.join(BASE, "extract_stats.json")))
results = {}
for g in GENES:
    L_aln, gapfrac = gap_stats(os.path.join(ALN, g + ".fasta"))
    r = {
        "n_taxa": extract[g]["n"],
        "len_mean_bp": extract[g]["len_mean"],
        "len_range_bp": [extract[g]["len_min"], extract[g]["len_max"]],
        "aln_len": L_aln,
        "indel_gap_fraction": round(gapfrac, 4),
    }
    for tag in ("multik", "k5", "k7"):
        p = os.path.join(TREES, f"{g}_fusang_{tag}.nwk")
        r[f"fusang_{tag}_nrf_vs_ref"] = round(ete_nrf(p, REF), 4)
    r["ft2_nrf_vs_ref"] = round(ete_nrf(os.path.join(TREES, f"{g}_ft2.nwk"), REF), 4)
    r["fusang_multik_vs_ft2_nrf"] = round(
        ete_nrf(os.path.join(TREES, f"{g}_fusang_multik.nwk"),
                os.path.join(TREES, f"{g}_ft2.nwk")), 4)
    results[g] = r
    print(f"{g:6s} L={r['len_mean_bp']:7.1f} gap={gapfrac:.3f} "
          f"mk={r['fusang_multik_nrf_vs_ref']:.3f} k5={r['fusang_k5_nrf_vs_ref']:.3f} "
          f"k7={r['fusang_k7_nrf_vs_ref']:.3f} ft2={r['ft2_nrf_vs_ref']:.3f} "
          f"mk_vs_ft2={r['fusang_multik_vs_ft2_nrf']:.3f}")

# summary
def summ(key):
    v = [results[g][key] for g in GENES]
    return {"mean": round(float(np.mean(v)), 4), "sd": round(float(np.std(v, ddof=1)), 4),
            "median": round(float(np.median(v)), 4),
            "min": round(min(v), 4), "max": round(max(v), 4)}

summary = {k: summ(k) for k in
           ("fusang_multik_nrf_vs_ref", "fusang_k5_nrf_vs_ref",
            "fusang_k7_nrf_vs_ref", "ft2_nrf_vs_ref", "fusang_multik_vs_ft2_nrf")}
print("\nsummary:")
for k, s in summary.items():
    print(f"  {k:28s} mean {s['mean']:.4f} +/- {s['sd']:.4f}  "
          f"median {s['median']:.4f}  range [{s['min']}, {s['max']}]")

json.dump({"genes": results, "summary": summary,
           "n_genes": len(GENES), "n_taxa_per_gene": 25,
           "reference": "Fischer et al. 2013 whole-mitogenome tree (AFproject fish_mito)",
           "scoring": "ete3 Tree.compare unrooted=True, norm_rf (max_rf=2(n-3))"},
          open(os.path.join(BASE, "gene_dna_results.json"), "w"), indent=1)
print("\nsaved gene_dna_results.json")
