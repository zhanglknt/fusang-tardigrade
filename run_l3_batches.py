#!/usr/bin/env python3
"""
Run L3 validation in 3 batches of 10 seeds each.
Each batch saves intermediate results; final merge at end.
"""
import subprocess
import sys
import json
import numpy as np

BATCHES = [(1, 10), (11, 10), (21, 10)]

for start, nseeds in BATCHES:
    print(f"\n{'#' * 70}")
    print(f"# BATCH: seeds {start}-{start+nseeds-1}")
    print(f"#{'#' * 70}\n")
    
    result = subprocess.run([
        sys.executable, "validate_l3_e2e.py",
        "--start-index", str(start),
        "--n-seeds", str(nseeds),
        "--n-taxa", "200",
        "--seq-len", "1000",
        "--sub", "0.05",
        "--indel", "0.02",
    ])
    
    if result.returncode != 0:
        print(f"WARNING: Batch {start}-{start+nseeds-1} returned rc={result.returncode}")

# Merge all results
print(f"\n{'=' * 70}")
print("MERGING ALL BATCH RESULTS")
print(f"{'=' * 70}")

all_results = []
for start, _ in BATCHES:
    path = "l3_validation_n200/l3_validation_results.json"
    try:
        with open(path, 'r') as f:
            data = json.load(f)
            batch_results = data.get("results", [])
            all_results.extend(batch_results)
            print(f"  Batch start={start}: {len(batch_results)} results loaded")
    except Exception as e:
        print(f"  WARNING: Could not load batch start={start}: {e}")

# Deduplicate by seed
seen = set()
unique = []
for r in all_results:
    s = r.get("seed")
    if s is not None and s not in seen:
        seen.add(s)
        unique.append(r)
    elif s is not None:
        print(f"  Skipping duplicate seed {s}")

print(f"  Total unique results: {len(unique)}")

# Compute final summary
l0_nrfs = [r.get("l0_nrf") for r in unique if r.get("l0_nrf") is not None]
l1_nrfs = [r.get("l1_nrf") for r in unique if r.get("l1_nrf") is not None]
ft2_nrfs = [r.get("ft2_nrf") for r in unique if r.get("ft2_nrf") is not None]
ft2_failed = sum(1 for r in unique if r.get("ft2_nrf") is None and "skipped" not in str(r.get("ft2_status", "")))

final_summary = {
    "config": {"n": 200, "L": 1000, "sub": 0.05, "indel": 0.02},
    "summary": {
        "l0": {"mean": float(np.mean(l0_nrfs)) if l0_nrfs else None,
               "std": float(np.std(l0_nrfs, ddof=1)) if l0_nrfs else None,
               "n": len(l0_nrfs)},
        "l1": {"mean": float(np.mean(l1_nrfs)) if l1_nrfs else None,
               "std": float(np.std(l1_nrfs, ddof=1)) if l1_nrfs else None,
               "n": len(l1_nrfs)},
        "ft2": {"mean": float(np.mean(ft2_nrfs)) if ft2_nrfs else None,
                "std": float(np.std(ft2_nrfs, ddof=1)) if ft2_nrfs else None,
                "n": len(ft2_nrfs)},
    },
    "ft2_failed": ft2_failed,
    "results": unique,
}

with open("l3_validation_n200/l3_validation_results.json", 'w') as f:
    json.dump(final_summary, f, indent=2, default=str)

print(f"\n{'=' * 70}")
print("FINAL SUMMARY")
print(f"{'=' * 70}")

if l0_nrfs:
    print(f"\n  L0 (k-mer k=5,gap2 NJ): {np.mean(l0_nrfs):.4f} +/- {np.std(l0_nrfs, ddof=1):.4f} (n={len(l0_nrfs)})")
if l1_nrfs:
    print(f"  L1 (multi-k k=5,7,9 NJ): {np.mean(l1_nrfs):.4f} +/- {np.std(l1_nrfs, ddof=1):.4f} (n={len(l1_nrfs)})")
if ft2_nrfs:
    print(f"  L3/FT2 (MAFFT+FT2 GTR): {np.mean(ft2_nrfs):.4f} +/- {np.std(ft2_nrfs, ddof=1):.4f} (n={len(ft2_nrfs)}, {ft2_failed} failed)")

# Comparison: L1 vs FT2
if l1_nrfs and ft2_nrfs:
    common_l1 = [r.get("l1_nrf") for r in unique if r.get("ft2_nrf") is not None]
    from scipy.stats import wilcoxon
    stat, p = wilcoxon(common_l1, ft2_nrfs)
    print(f"\n  L1 vs FT2: delta = {np.mean(common_l1) - np.mean(ft2_nrfs):+.4f}, Wilcoxon p = {p:.4f}")
    if len(l0_nrfs) >= len(ft2_nrfs):
        common_l0 = [r.get("l0_nrf") for r in unique if r.get("ft2_nrf") is not None]
        stat0, p0 = wilcoxon(common_l0, ft2_nrfs)
        print(f"  L0 vs FT2: delta = {np.mean(common_l0) - np.mean(ft2_nrfs):+.4f}, Wilcoxon p = {p0:.4f}")

print(f"\nFinal merged results saved to l3_validation_n200/l3_validation_results.json")
