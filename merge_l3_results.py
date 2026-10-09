#!/usr/bin/env python3
"""Merge Linux-rerun L3 (MAFFT+FastTree2, -nt -gtr -nosupport) trees with existing
L0/L1 results from l3_validation_n200/l3_validation_results.json.

Uses the exact same compute_nrf() as validate_l3_e2e.py for comparability.
Outputs updated summary + paired statistics (L1 vs L3, L0 vs L3, L1 vs L0).
"""
import json
import sys
import numpy as np
from pathlib import Path
from scipy.stats import wilcoxon

WORKDIR = Path(r"D:\系统发育树项目\Fusang\Fusang-main")
sys.path.insert(0, str(WORKDIR))
from validate_l3_e2e import compute_nrf  # noqa: E402  (same function as original run)

VAL_DIR = WORKDIR / "l3_validation_n200"
LINUX_DIR = VAL_DIR / "linux_out"
RESULTS = VAL_DIR / "l3_validation_results.json"


def main():
    data = json.loads(RESULTS.read_text(encoding="utf-8"))
    rows = data["results"]
    by_seed = {r["seed"]: r for r in rows}

    n_done, n_fail = 0, 0
    for i in range(1, 31):
        lin = LINUX_DIR / f"seed{i:03d}_ft2_lin.nwk"
        true_p = VAL_DIR / f"seed{i:03d}_true.nwk"
        if not lin.exists():
            print(f"seed {i}: MISSING linux tree")
            n_fail += 1
            continue
        true_nwk = true_p.read_text(encoding="utf-8", errors="replace").strip()
        ft2_nwk = lin.read_text(encoding="utf-8", errors="replace").strip()
        if not ft2_nwk or len(ft2_nwk) < 20:
            print(f"seed {i}: EMPTY linux tree")
            n_fail += 1
            continue
        nrf = compute_nrf(true_nwk, ft2_nwk)
        if nrf is None:
            print(f"seed {i}: nRF computation failed")
            n_fail += 1
            continue
        if i in by_seed:
            by_seed[i]["ft2_nrf"] = nrf
            by_seed[i]["ft2_status"] = "OK-linux"
            by_seed[i]["ft2_protocol"] = "mafft --auto + FastTree -nt -gtr -nosupport (WSL2)"
        else:
            rows.append({"seed": i, "ft2_nrf": nrf, "ft2_status": "OK-linux"})
        n_done += 1

    print(f"\nLinux L3 trees merged: {n_done} ok, {n_fail} failed")

    # Aggregate
    l0 = np.array([r["l0_nrf"] for r in rows if r.get("l0_nrf") is not None])
    l1 = np.array([r["l1_nrf"] for r in rows if r.get("l1_nrf") is not None])
    ft2 = np.array([r["ft2_nrf"] for r in rows if r.get("ft2_nrf") is not None])
    print(f"\nL0 (k=5,gap2 NJ)      : {l0.mean():.4f} ± {l0.std(ddof=1):.4f} (n={len(l0)})")
    print(f"L1 (multi-k NJ)       : {l1.mean():.4f} ± {l1.std(ddof=1):.4f} (n={len(l1)})")
    print(f"L3 (MAFFT+FT2, Linux) : {ft2.mean():.4f} ± {ft2.std(ddof=1):.4f} (n={len(ft2)})")

    # Paired stats on common seeds (need matched triple)
    matched = [r for r in rows if r.get("l0_nrf") is not None and r.get("l1_nrf") is not None and r.get("ft2_nrf") is not None]
    if len(matched) >= 5:
        a_l0 = np.array([r["l0_nrf"] for r in matched])
        a_l1 = np.array([r["l1_nrf"] for r in matched])
        a_l3 = np.array([r["ft2_nrf"] for r in matched])
        print(f"\nPaired analyses (n={len(matched)} matched seeds):")
        for name, x, y in [("L1 vs L3", a_l1, a_l3), ("L0 vs L3", a_l0, a_l3), ("L1 vs L0", a_l1, a_l0)]:
            diff = x - y
            stat, p = wilcoxon(x, y)
            d = diff.mean() / (diff.std(ddof=1) + 1e-12)
            print(f"  {name}: mean diff = {diff.mean():+.4f}, Wilcoxon p = {p:.4g}, paired Cohen's d = {d:+.2f}")
        # wins
        print(f"  L1 lower than L3 in {(a_l1 < a_l3).sum()}/{len(matched)} seeds")

    # Update summary and save (keep original as _pre_linux backup)
    backup = VAL_DIR / "l3_validation_results_pre_linux.json"
    if not backup.exists():
        backup.write_text(RESULTS.read_text(encoding="utf-8"), encoding="utf-8")
    data["config"]["l3_rerun"] = "WSL2 Ubuntu-24.04, mafft --auto, FastTree -nt -gtr -nosupport (2026-10-09)"
    data["summary"] = {
        "l0": {"mean": float(l0.mean()), "std": float(l0.std(ddof=1)), "n": int(len(l0))},
        "l1": {"mean": float(l1.mean()), "std": float(l1.std(ddof=1)), "n": int(len(l1))},
        "ft2": {"mean": float(ft2.mean()), "std": float(ft2.std(ddof=1)), "n": int(len(ft2))},
    }
    RESULTS.write_text(json.dumps(data, indent=2, default=str), encoding="utf-8")
    print(f"\nUpdated: {RESULTS} (backup: {backup})")


if __name__ == "__main__":
    main()
