"""n=1000 indel simplified-pipeline rerun (forensics F1 follow-up).

Reruns fusang_v2.py --simple on the 30 INDELible seeds (100-129) whose
DCM+EPA trees populate Table 1's n=1000 indel row, then scores simple vs
DCM vs FT2 against the true trees with the standard unrooted RF
(nrf_standard.py). Windows-native (fastme.exe), matching the original
benchmark platform. Idempotent: skips seeds whose output tree exists.
"""
import csv
import json
import os
import subprocess
import sys
import time

ROOT = r"D:/系统发育树项目/Fusang/Fusang-main"
sys.path.insert(0, ROOT)
from nrf_standard import parse_newick, splits  # noqa: E402

PY = r"C:/Users/admin/.workbuddy/binaries/python/envs/default/Scripts/python.exe"
SRC = os.path.join(ROOT, "indelible_n1000_indel")
OUT = os.path.join(ROOT, "_forensics_v5", "n1000_simplified_rerun")
os.makedirs(OUT, exist_ok=True)
N = 1000
DENOM = 2 * (N - 3)


def nrf_of(est_path, true_splits):
    est_splits = splits(parse_newick(open(est_path).read()))
    rf = len(true_splits - est_splits) + len(est_splits - true_splits)
    return rf / DENOM


def main():
    rows = []
    for seed in range(100, 130):
        fasta = os.path.join(SRC, f"seed{seed}_indel.fasta")
        out = os.path.join(OUT, f"seed{seed}_fusang_simple.nwk")
        dt = None
        if not os.path.exists(out):
            t0 = time.time()
            r = subprocess.run(
                [PY, os.path.join(ROOT, "fusang_v2.py"), "--input", fasta,
                 "--output", out, "--simple", "--tree_method", "fastme",
                 "--threads", "4"],
                capture_output=True, text=True, timeout=900, cwd=ROOT)
            if r.returncode != 0 or not os.path.exists(out):
                r = subprocess.run(
                    [PY, os.path.join(ROOT, "fusang_v2.py"), "--input", fasta,
                     "--output", out, "--simple", "--tree_method", "nj",
                     "--threads", "4"],
                    capture_output=True, text=True, timeout=900, cwd=ROOT)
            dt = time.time() - t0
            ok = r.returncode == 0 and os.path.exists(out)
            print(f"[seed {seed}] run {'OK' if ok else 'FAIL'} {dt:.1f}s",
                  flush=True)
            if not ok:
                rows.append((seed, None, None, None, dt))
                continue
        true_splits = splits(parse_newick(
            open(os.path.join(SRC, f"seed{seed}_indel_true.nwk")).read()))
        nrf_simple = nrf_of(out, true_splits)
        nrf_dcm = nrf_of(os.path.join(SRC, f"seed{seed}_indel_fusang.nwk"),
                         true_splits)
        nrf_ft2 = nrf_of(os.path.join(SRC, f"seed{seed}_indel_ft2.nwk"),
                         true_splits)
        rows.append((seed, nrf_simple, nrf_dcm, nrf_ft2, dt))
        print(f"[seed {seed}] simple={nrf_simple:.4f} dcm={nrf_dcm:.4f} "
              f"ft2={nrf_ft2:.4f}", flush=True)

    csv_path = os.path.join(OUT, "n1000_indel_simple_rerun.csv")
    with open(csv_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["seed", "nrf_simple", "nrf_dcm", "nrf_ft2", "seconds"])
        w.writerows(rows)

    import numpy as np
    from scipy.stats import wilcoxon
    good = [r for r in rows if r[1] is not None]
    s = np.array([r[1] for r in good])
    d = np.array([r[2] for r in good])
    f = np.array([r[3] for r in good])

    def paired(a, b):
        wres = wilcoxon(a, b, zero_method="zsplit")
        dz = float(np.mean(a - b) / np.std(a - b, ddof=1))
        return {"p": wres.pvalue, "dz": dz,
                "a_better": int((a < b).sum()), "b_better": int((a > b).sum()),
                "ties": int((a == b).sum())}

    summary = {
        "n_valid": len(good),
        "seed_range": [100, 129],
        "simple_mean": float(s.mean()), "simple_sd": float(s.std(ddof=1)),
        "dcm_mean": float(d.mean()), "dcm_sd": float(d.std(ddof=1)),
        "ft2_mean": float(f.mean()), "ft2_sd": float(f.std(ddof=1)),
        "simple_vs_dcm": paired(s, d),
        "simple_vs_ft2": paired(s, f),
        "dcm_vs_ft2": paired(d, f),
        "metric": "nrf_standard.py, denom=2(n-3)=1994",
        "runner": "fusang_v2.py --simple --tree_method fastme --threads 4 "
                  "(nj fallback), Windows venv python",
    }
    with open(os.path.join(OUT, "n1000_indel_simple_rerun_summary.json"),
              "w") as fh:
        json.dump(summary, fh, indent=2)
    print("SUMMARY:", json.dumps(summary, indent=1), flush=True)


if __name__ == "__main__":
    main()
