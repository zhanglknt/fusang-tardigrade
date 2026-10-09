#!/usr/bin/env python3
"""Recover corrupted benchmark data files (seeds 100-126 AF benchmark set).

Strategy:
1. Validate determinism of gen_test_data_indel.py by regenerating seed 100
   (intact original) and byte-comparing FASTA + TRUE tree.
2. Regenerate TRUE trees + FASTAs for corrupted seeds (110/121/126 true;
   111/114/115/117/121 fasta already reconstructed from alignments - now
   regenerated exactly).
3. Regenerate FT2 trees (MAFFT mafft-win + FastTree -nt -gtr -nosupport)
   for 9 corrupted seeds.
4. Regenerate Fusang trees (fusang_v2.py, fastme, k=5, gap2) for 4 corrupted seeds.
5. Validate EVERYTHING per-seed against indel_benchmark_MASTER_REBUILT.csv
   (ft2_nrf, fusang_nrf columns, vs TRUE tree).

Corrupted originals are quarantined to corrupt_fasta_quarantine/ (never deleted).
"""
import hashlib
import os
import shutil
import subprocess
import sys
import csv
from pathlib import Path

WORKDIR = Path(r"D:\系统发育树项目\Fusang\Fusang-main")
QUAR = WORKDIR / "corrupt_fasta_quarantine"
PY = r"C:\Users\admin\.workbuddy\binaries\python\envs\default\Scripts\python.exe"
MAFFT_DIR = r"d:\系统发育树项目\Fusang\bench_tools\mafft-win\mafft-win"
FASTTREE = r"d:\系统发育树项目\Fusang\bench_tools\FastTree.exe"
QUAR.mkdir(exist_ok=True)

CORRUPT_TRUE = [110, 121, 126]
CORRUPT_FASTA = [111, 114, 115, 117, 121]     # reconstructed from aligned earlier
CORRUPT_FT2 = [109, 114, 115, 118, 119, 120, 121, 122, 126]
CORRUPT_FUSANG = [114, 115, 120, 125]
ALL_REGEN_FASTA = sorted(set(CORRUPT_TRUE + CORRUPT_FASTA + CORRUPT_FT2 + CORRUPT_FUSANG) | {100})


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]


def run(cmd, timeout=300):
    r = subprocess.run(cmd, shell=True, capture_output=True, timeout=timeout, cwd=WORKDIR)
    return r.returncode, r.stdout, r.stderr


def quarantine(seed, kind, ext="nwk"):
    src = WORKDIR / f"seed{seed}_indel_{kind}.{ext}"
    if src.exists():
        dst = QUAR / f"seed{seed}_indel_{kind}.{ext}.corrupt"
        if not dst.exists():
            shutil.copy2(src, dst)
        print(f"  quarantined {src.name}")


def gen_data(seed):
    """Run gen_test_data_indel.py for a seed; returns (fasta, true_nwk) paths."""
    cmd = f'"{PY}" "{WORKDIR / "gen_test_data_indel.py"}" 200 500 0.05 0.02 {seed}'
    rc, out, err = run(cmd, timeout=120)
    if rc != 0:
        print(f"  GEN FAILED seed {seed}: {err[:200]}")
        return None, None
    fasta = WORKDIR / "test_indel_n200.fasta"
    true = WORKDIR / "test_indel_n200_true.nwk"
    return fasta, true


def main():
    # ---- 1. determinism validation on seed 100 ----
    print("=== Step 1: determinism check (seed 100) ===")
    f100, t100 = gen_data(100)
    orig_f = WORKDIR / "seed100_indel.fasta"
    orig_t = WORKDIR / "seed100_indel_true.nwk"
    f_match = f100.read_bytes() == orig_f.read_bytes()
    t_match = t100.read_bytes() == orig_t.read_bytes()
    print(f"  FASTA byte-identical: {f_match} ({sha(f100)} vs {sha(orig_f)})")
    print(f"  TRUE  byte-identical: {t_match} ({sha(t100)} vs {sha(orig_t)})")
    if not (f_match and t_match):
        print("  !! Non-deterministic generation - inspect before proceeding")
        sys.exit(1)

    # ---- 2. regenerate FASTA + TRUE for corrupted seeds ----
    print("\n=== Step 2: regenerate FASTA + TRUE trees ===")
    regen_seeds = sorted(set(CORRUPT_TRUE + CORRUPT_FASTA))
    for seed in regen_seeds:
        f, t = gen_data(seed)
        if f is None:
            continue
        # fasta
        dst_f = WORKDIR / f"seed{seed}_indel.fasta"
        quarantine(seed, "fasta", "fasta")
        shutil.copy2(f, dst_f)
        # true
        dst_t = WORKDIR / f"seed{seed}_indel_true.nwk"
        quarantine(seed, "true")
        shutil.copy2(t, dst_t)
        # note: keep aligned file for reference; it will be regenerated in step 3 if needed
        print(f"  seed {seed}: fasta+true regenerated (fasta sha {sha(dst_f)})")
        # cross-check reconstruction-from-alignment (for the 5 previously reconstructed)
        aln = WORKDIR / f"seed{seed}_indel_aligned.fasta"
        if aln.exists():
            seqs, order = {}, []
            name = None
            for ln in open(aln, encoding="ascii", errors="ignore"):
                ln = ln.strip()
                if ln.startswith(">"):
                    name = ln[1:].split()[0]
                    seqs[name] = ""
                    order.append(name)
                elif name:
                    seqs[name] += ln
            recon = {k: seqs[k].replace("-", "").upper() for k in order}
            gen_seqs, gorder = {}, []
            name = None
            for ln in open(dst_f):
                ln = ln.strip()
                if ln.startswith(">"):
                    name = ln[1:]
                    gen_seqs[name] = ""
                    gorder.append(name)
                elif name:
                    gen_seqs[name] += ln
            same = gorder == order and all(gen_seqs[k] == recon.get(k) for k in gorder)
            print(f"    ungap-reconstruction matches exact regeneration: {same}")

    # ---- 3. regenerate FT2 trees (FastTree on INTACT original alignments) ----
    print("\n=== Step 3: regenerate FT2 trees (FastTree on intact alignments) ===")
    for seed in CORRUPT_FT2:
        aln = WORKDIR / f"seed{seed}_indel_aligned.fasta"
        raw = aln.read_bytes()
        if not (raw.startswith(b">") and b"\x00" not in raw and raw.count(b">") == 200):
            print(f"  seed {seed}: aligned file NOT intact - needs MAFFT rerun (skipped)")
            continue
        ft2 = WORKDIR / f"seed{seed}_indel_ft2.nwk"
        cmd = f'"{FASTTREE}" -nt -gtr -nosupport < "{aln}"'
        rc, out, err = run(cmd, timeout=300)
        if rc != 0 or len(out) < 100:
            print(f"  seed {seed}: FastTree FAILED rc={rc}")
            continue
        quarantine(seed, "ft2")
        ft2.write_bytes(out)
        print(f"  seed {seed}: FT2 regenerated ({len(out)} bytes)")

    # ---- 4. regenerate Fusang trees ----
    print("\n=== Step 4: regenerate Fusang trees ===")
    for seed in CORRUPT_FUSANG:
        fasta = WORKDIR / f"seed{seed}_indel.fasta"
        out_tree = WORKDIR / f"seed{seed}_indel_fusang.nwk"
        cmd = (f'"{PY}" "{WORKDIR / "fusang_v2.py"}" --input "{fasta}" '
               f'--output "{out_tree}" --tree_method fastme --kmer_k 5 --kmer_gap gap2')
        rc, out, err = run(cmd, timeout=300)
        ok = out_tree.exists() and out_tree.stat().st_size > 1000
        if ok:
            quarantine(seed, "fusang")
            print(f"  seed {seed}: Fusang regenerated ({out_tree.stat().st_size} bytes)")
        else:
            print(f"  seed {seed}: Fusang FAILED rc={rc} {err[:150]}")

    # ---- 5. validate against MASTER CSV ----
    print("\n=== Step 5: validation vs MASTER CSV (nRF vs TRUE) ===")
    sys.path.insert(0, str(WORKDIR))
    from calc_nrf_simple import calc_nrf
    import tempfile
    master = {int(r["seed"]): r for r in csv.DictReader(open(WORKDIR / "repro_package/data/indel_benchmark_MASTER_REBUILT.csv"))}

    def nrf(p1, p2):
        return calc_nrf(str(p1), str(p2))

    n_checked, n_match, mismatches = 0, 0, []
    for seed in range(100, 127):
        true_p = WORKDIR / f"seed{seed}_indel_true.nwk"
        ft2_p = WORKDIR / f"seed{seed}_indel_ft2.nwk"
        fus_p = WORKDIR / f"seed{seed}_indel_fusang.nwk"
        if not (true_p.exists() and ft2_p.exists() and fus_p.exists()):
            print(f"  seed {seed}: missing trees, skip")
            continue
        row = master.get(seed)
        if not row:
            continue
        for kind, tree_p, col in [("ft2", ft2_p, "ft2_nrf"), ("fusang", fus_p, "fusang_nrf")]:
            expected = row.get(col)
            if not expected:
                continue
            got = nrf(true_p, tree_p)
            n_checked += 1
            if got is not None and abs(got - float(expected)) < 0.005:
                n_match += 1
            else:
                mismatches.append((seed, kind, expected, got))
    print(f"  checked {n_checked} values, matched {n_match}")
    for m in mismatches:
        print(f"  MISMATCH: seed {m[0]} {m[1]}: CSV={m[2]} recomputed={m[3]}")
    if n_checked == n_match:
        print("\nALL REGENERATED DATA VALIDATED AGAINST MASTER CSV")


if __name__ == "__main__":
    main()
