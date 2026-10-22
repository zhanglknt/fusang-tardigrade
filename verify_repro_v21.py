"""
Full verification of repro_package_v2.1.zip.

1. Extract the zip to a scratch directory (never the build staging tree).
2. Re-run the PACKAGED nrf_standard.py on a real pair of packaged trees
   (seed 300: L1 vs TRUE) and compare with the packaged per-seed CSV.
3. Recompute every sha256 recorded anywhere in data/DATA_SOURCES.json and
   compare with the extracted bytes.
4. List the zip and check for unwanted content (.git, __pycache__, caches).

Run:  python verify_repro_v21.py
"""
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import zipfile

ROOT = r"D:\系统发育树项目\Fusang\Fusang-main"
os.chdir(ROOT)

ZIP = "repro_package_v2.1.zip"
VER = "repro_v21_verify2"
PY = sys.executable


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


if os.path.exists(VER):
    raise SystemExit(f"{VER}/ exists; `rm -rf {VER}` first")
os.makedirs(VER)
with zipfile.ZipFile(ZIP) as z:
    z.extractall(VER)
print(f"[1] extracted {ZIP} -> {VER}/")

fail = 0

# ---------------------------------------------------------------- 2. recompute
print("\n[2] recompute nRF with the PACKAGED nrf_standard.py")
t_true = os.path.join(VER, "bench_l3_coalescent_indelible", "seed300_tab5_true.nwk")
t_l1 = os.path.join(VER, "bench_l3_coalescent_indelible", "seed300_tab5_l1.nwk")
for p in (t_true, t_l1):
    if not os.path.exists(p):
        raise SystemExit(f"MISSING packaged tree: {p}")

out = subprocess.run([PY, os.path.join(VER, "nrf_standard.py"), t_true, t_l1],
                     capture_output=True, text=True)
if out.returncode != 0:
    print(out.stdout, out.stderr)
    raise SystemExit("packaged nrf_standard.py failed to run")
recomputed = float(out.stdout.strip())
print(f"    nrf_standard.py  seed300 TRUE vs L1 = {recomputed:.12f}  (stdout: {out.stdout.strip()})")

csv_row = None
with open(os.path.join(VER, "data",
                       "bench_l3_coalescent_indelible_30seeds_nrfstandard.csv")) as fh:
    for row in csv.DictReader(fh):
        if row["seed"] == "300":
            csv_row = row
            break
if csv_row is None:
    raise SystemExit("seed 300 not found in the packaged CSV")
# The CSV carries BOTH metrics: unsuffixed columns = superseded calc_nrf_simple,
# "_std" columns = standard root-independent RF (the manuscript values).
import statistics as st

for col, tree in (("l1_true_std", "seed300_tab5_l1.nwk"),
                  ("l0_true_std", "seed300_tab5_l0.nwk")):
    csv_v = float(csv_row[col])
    out2 = subprocess.run([PY, os.path.join(VER, "nrf_standard.py"), t_true,
                           os.path.join(VER, "bench_l3_coalescent_indelible", tree)],
                          capture_output=True, text=True)
    v = float(out2.stdout.strip())
    d = abs(v - csv_v)
    print(f"    [{col:12s}] recomputed {v:.12f} vs CSV {csv_v:.12f}  |d|={d:.3e} "
          f"-> {'MATCH' if d <= 1e-9 else 'MISMATCH'}")
    if d > 1e-9:
        fail += 1

print(f"    (for reference, the superseded-metric column l1_true = "
      f"{float(csv_row['l1_true']):.12f} -- NOT the manuscript value)")

# sanity: the _std column means must reproduce manuscript Table 5
with open(os.path.join(VER, "data",
                       "bench_l3_coalescent_indelible_30seeds_nrfstandard.csv")) as fh:
    rows = list(csv.DictReader(fh))
print("    manuscript Table 5 cross-check (30 seeds):")
for col, expect in (("l0_true_std", 0.150), ("l1_true_std", 0.137), ("l3_true_std", 0.124)):
    m = st.mean(float(r[col]) for r in rows)
    sd = st.stdev(float(r[col]) for r in rows)
    hit = abs(m - expect) < 0.001
    print(f"        {col:12s} mean={m:.4f} sd={sd:.4f}  (manuscript {expect}) "
          f"-> {'OK' if hit else 'MISMATCH'}")
    if not hit:
        fail += 1

# ------------------------- 2b. n=1000 threshold benchmark, packaged trees ----
print("\n[2b] recompute the n=1000 threshold benchmark with the PACKAGED nrf_standard.py")
thr_csv = os.path.join(VER, "table_threshold_n1000_nrfstandard.csv")
with open(thr_csv) as fh:
    thr = {r["seed"]: r for r in csv.DictReader(fh)}
for seed in ("1000", "1007"):
    ref = os.path.join(VER, "bench_n1000_threshold", "true",
                       f"seed{seed}_clean_true.nwk")
    if not os.path.exists(ref):
        print(f"    seed {seed}: true tree MISSING")
        fail += 1
        continue
    for mode, col in (("simple", "nrf_simple"), ("dcm", "nrf_dcm")):
        t = os.path.join(VER, "bench_n1000_threshold",
                         f"seed{seed}_clean_{mode}.nwk")
        o = subprocess.run([PY, os.path.join(VER, "nrf_standard.py"), ref, t],
                           capture_output=True, text=True)
        if o.returncode != 0:
            print(f"    seed {seed} {mode}: nrf_standard.py FAILED: {o.stderr[:200]}")
            fail += 1
            continue
        v = float(o.stdout.strip())
        c = float(thr[seed][col])
        d = abs(v - c)
        print(f"    seed {seed} {mode:6s}: recomputed {v:.9f} vs CSV {c:.9f} "
              f"|d|={d:.3e} -> {'MATCH' if d <= 1e-8 else 'MISMATCH'}")
        if d > 1e-8:
            fail += 1

# ------------------------------------------------------- 3. sha256 full audit
print("\n[3] sha256 audit of data/DATA_SOURCES.json")
ds_path = os.path.join(VER, "data", "DATA_SOURCES.json")
with open(ds_path, encoding="utf-8") as fh:
    ds = json.load(fh)

entries = []


def walk(node, trail):
    if isinstance(node, dict):
        if "path" in node and "sha256" in node:
            entries.append((".".join(str(t) for t in trail), node["path"], node["sha256"]))
        for k, v in node.items():
            walk(v, trail + [k])
    elif isinstance(node, list):
        for i, v in enumerate(node):
            walk(v, trail + [i])


walk(ds, [])
print(f"    {len(entries)} path+sha256 entries found")

missing, mismatch, okc = [], [], 0
for loc, rel, want in entries:
    p = os.path.join(VER, rel.replace("/", os.sep))
    if not os.path.exists(p):
        missing.append((loc, rel))
        continue
    got = sha256(p)
    if got == want:
        okc += 1
    else:
        mismatch.append((loc, rel, want, got))
print(f"    match   : {okc}")
print(f"    missing : {len(missing)}")
for loc, rel in missing[:20]:
    print(f"        MISSING {rel}   ({loc})")
print(f"    mismatch: {len(mismatch)}")
for loc, rel, want, got in mismatch[:20]:
    print(f"        MISMATCH {rel}\n            want {want}\n            got  {got}   ({loc})")
fail += len(missing) + len(mismatch)

# ------------------------------------------------------------- 4. zip hygiene
print("\n[4] zip content hygiene")
bad_dirs = {".git", "__pycache__", ".cache", ".cache_rescore", ".ipynb_checkpoints"}
with zipfile.ZipFile(ZIP) as z:
    names = z.namelist()
    infos = z.infolist()
    total_uncompressed = sum(i.file_size for i in infos)
    bad = [n for n in names
           if any(f"/{d}/" in f"/{n}/" or n.startswith(d + "/") for d in bad_dirs)
           or n.endswith((".pyc", ".pyo", ".tmp", ".bak"))]
    print(f"    entries            : {len(names)}")
    print(f"    zip size           : {os.path.getsize(ZIP)} bytes "
          f"({os.path.getsize(ZIP)/1e6:.2f} MB)")
    print(f"    uncompressed size  : {total_uncompressed} bytes "
          f"({total_uncompressed/1e6:.2f} MB)")
    print(f"    forbidden entries  : {len(bad)} {bad[:10]}")

# presence checks for the v2.1 deliverables
required = [
    "nrf_standard.py",
    "rf_verification_record.md",
    "real_data_rf_recomputed.md",
    "fusang/calc_nrf_simple_SUPERSEDED.py",
    "scripts/gen_test_data_indel_DEFECTIVE_SUPERSEDED.py",
    "scripts/gen_indelible_sim.py",
    "scripts/run_bench_indelible.py",
    "environment.yml",
    "indelible_n200_indel/bench_n200_indel_130seeds_nrfstandard.csv",
    "indelible_n200_indel/competitors_n200_indel_130seeds_stdrf.csv",
    "indelible_n200_indel/competitors_n200_indel_130seeds_stdrf_summary.json",
    "indelible_iqtree2/iqtree2_n200_indel_nrfstandard.csv",
    "indelible_iqtree2/iqtree2_n200_indel_nrfstandard.json",
    "data/bench_l3_coalescent_indelible_30seeds_nrfstandard.csv",
    "data/bench_l3_coalescent_indelible_30seeds_nrfstandard_summary.json",
    "data/bench_mash_indelible_n200_nrfstandard.csv",
    "data/bench_mash_indelible_n200_nrfstandard.json",
    "repro_ablations/cosine_vs_jsd_nrfstandard.py",
    "repro_ablations/cosine_vs_jsd_n200_indel_30seeds_nrfstandard.json",
    "repro_ablations/ablation_dcm_seed42_indel_nrfstandard.py",
    "repro_ablations/ablation_dcm_seed42_indel_nrfstandard.json",
    "repro_ablations/skmer_probe_nrfstandard.py",
    "repro_ablations/skmer_probe_nrfstandard.json",
    "bench_l3_coalescent_indelible/seed300_tab5_true.nwk",
    "bench_l3_coalescent_indelible/seed329_tab5_l1.nwk",
    "swisstree_rf_recomputed/swisstree_rf_standard.csv",
    "balibase_rf_recomputed/balibase_rf_standard.csv",
    "indel_raw_seeds100_129/README_DEFECTIVE.md",
    "table1_nrfstandard/TABLE1_all_summary.json",
    "table1_nrfstandard/n200_indel_nrfstandard.csv",
    "table1_nrfstandard/indel0.05_nrfstandard.csv",
    "scripts/rescore_table1_nrfstandard.py",
    "table_threshold_n1000_nrfstandard.csv",
    "table_threshold_n1000_summary.json",
    "bench_n1000_threshold/ENVIRONMENT.md",
    "bench_n1000_threshold/run_threshold_benchmark.py",
    "bench_n1000_threshold/seed1009_clean_dcm.nwk",
    "bench_n1000_threshold/true/seed1009_clean_true.nwk",
    "figures/Figure3.png",
    "figures/Supplementary_Figure_S7.pdf",
]
# the shipped fusang_v2.py must carry the --simple CLI fix
with zipfile.ZipFile(ZIP) as z:
    fv2 = z.read("fusang/fusang_v2.py").decode("utf-8")
has_fix = "args.force_simple is not None" in fv2
print(f"    fusang/fusang_v2.py --simple fix (4e6bf02) present: "
      f"{'YES' if has_fix else 'NO'}")
if not has_fix:
    fail += 1
absent = [r for r in required if r not in names]
print(f"    required present    : {len(required)-len(absent)}/{len(required)}")
if absent:
    print(f"        ABSENT: {absent}")
fail += len(absent) + len(bad)

# residue scan: any remaining false provenance claim
print("\n[5] residue scan for withdrawn provenance claims")
residue = []
scan_files = [n for n in names if n.endswith((".md", ".py", ".json", ".ipynb"))]
for n in scan_files:
    p = os.path.join(VER, n.replace("/", os.sep))
    try:
        with open(p, encoding="utf-8", errors="ignore") as fh:
            t = fh.read()
    except OSError:
        continue
    for pat in ("deterministic simulator that generated",
                "simulator that generated the indel benchmark data",
                "gen_test_data_indel.py",
                "the actual generator of record",
                "authoritative inputs for the 130-seed"):
        if pat in t:
            # allow the renaming / withdrawal prose itself
            residue.append((n, pat))
# Files that legitimately still name `gen_test_data_indel.py`: they are the
# withdrawal documentation itself, or a historical "replaces X" note.
allow = {"README.md",
         "indel_raw_seeds100_129/README.md",
         "indel_raw_seeds100_129/README_DEFECTIVE.md",
         "scripts/gen_indelible_sim.py",
         "scripts/gen_test_data_indel_DEFECTIVE_SUPERSEDED.py",
         "data/DATA_SOURCES.json"}
real_residue = [(n, p) for n, p in residue if n not in allow]
print(f"    files with old wording (expected, they document the withdrawal): {sorted(set(n for n,_ in residue))}")
print(f"    unexplained residue outside README/DATA_SOURCES: {len(real_residue)}")
for n, p in real_residue:
    print(f"        {n}: {p}")
fail += len(real_residue)

print(f"\n=== VERIFICATION {'PASSED' if fail == 0 else f'FAILED ({fail} problems)'} ===")
