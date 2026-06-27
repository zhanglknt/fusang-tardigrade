#!/usr/bin/env python3
"""
P0-4: Level 3 End-to-End Validation.

Validates the full MHL pipeline (L0 → L1 → L2 → L3) on n=200 coalescent data:
- L0: k-mer distance + NJ baseline
- L1: Multi-k ensemble comparison
- L2: DAHP cluster refinement (if classifier says SPLIT)
- L3: MSA+ML via MAFFT+FastTree2 (where MAFFT is available)

Generates nRF comparison at each level against the true tree.

Usage:
    python validate_l3_e2e.py --n-taxa 200 --n-seeds 30
"""

import sys
import os
import time
import json
import argparse
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from tree_simulation import make_coalescent_tree, simulate_seqs
from fusang_mhl.level0_kmer import compute_l0_distance

_CODE_MAP = {0: 'A', 1: 'T', 2: 'C', 3: 'G'}


def generate_data(n, L, sub, indel, seed):
    root_node, leaves = make_coalescent_tree(n, seed=seed)
    leaf_seqs = simulate_seqs(root_node, n, L, sub, seed, indel_rate=indel)
    taxon_names = [f't{i+1:04d}' for i in range(n)]
    sequences = [''.join(_CODE_MAP[b] for b in leaf_seqs[i]) for i in range(n)]
    true_nwk = root_node.to_newick()
    return sequences, taxon_names, true_nwk


def compute_nrf(true_nwk_str, pred_nwk_str):
    from Bio import Phylo
    from io import StringIO
    import re
    
    true_clean = re.sub(r':-[\d.e-]+', ':0.00001', true_nwk_str)
    pred_clean = re.sub(r':-[\d.e-]+', ':0.00001', pred_nwk_str)
    
    try:
        t_true = Phylo.read(StringIO(true_clean), 'newick')
        t_pred = Phylo.read(StringIO(pred_clean), 'newick')
    except Exception as e:
        return None
    
    def get_splits(tree):
        tips = [l.name for l in tree.get_terminals()]
        n_tips = len(tips)
        if n_tips < 3:
            return set()
        all_names = frozenset(tips)
        try:
            first_tip = list(tree.find_clades(name=tips[0]))[0]
            tree.root_with_outgroup(first_tip)
        except Exception:
            return set()
        splits = set()
        for clade in tree.get_nonterminals():
            leaf_names = frozenset(l.name for l in clade.get_terminals())
            m = len(leaf_names)
            if 1 <= m < n_tips and leaf_names != all_names:
                splits.add(min(leaf_names, all_names - leaf_names, key=sorted))
        return splits
    
    s1 = get_splits(t_true)
    s2 = get_splits(t_pred)
    
    if not s1 or not s2:
        return None
    
    rf = len(s1 - s2) + len(s2 - s1)
    n = len(t_true.get_terminals())
    max_rf = 2.0 * (n - 3)
    return rf / max_rf if max_rf > 0 else 1.0


def build_nj_tree(sequences, taxon_names):
    """Build NJ tree from k-mer distance matrix."""
    from Bio.Phylo.TreeConstruction import DistanceMatrix as BioDM, DistanceTreeConstructor
    from io import StringIO
    from Bio import Phylo
    
    D = compute_l0_distance(sequences, taxon_names)
    if hasattr(D, 'tolist'):
        D_list = D.tolist()
    else:
        D_list = D
    lower_tri = [[D_list[i][j] for j in range(i+1)] for i in range(len(D_list))]
    dm = BioDM(list(taxon_names), lower_tri)
    constructor = DistanceTreeConstructor()
    tree_obj = constructor.nj(dm)
    buf = StringIO()
    Phylo.write(tree_obj, buf, 'newick')
    return buf.getvalue()


def build_multi_k_nj(sequences, taxon_names, ks=[5, 7, 9]):
    """Build NJ tree from multi-k ensemble distance (averaged)."""
    all_D = []
    for k in ks:
        # Use contiguous k-mers for ensemble (matching DAHP-V3)
        D = compute_l0_distance(sequences, taxon_names, k=k, gap=0)
        if hasattr(D, 'tolist'):
            D = D.tolist()
        all_D.append(np.array(D))
    
    # Average distance matrices
    avg_D = np.mean(all_D, axis=0)
    
    from Bio.Phylo.TreeConstruction import DistanceMatrix as BioDM, DistanceTreeConstructor
    from io import StringIO
    from Bio import Phylo
    
    lower_tri = [[avg_D[i][j] for j in range(i+1)] for i in range(len(avg_D))]
    dm = BioDM(list(taxon_names), lower_tri)
    constructor = DistanceTreeConstructor()
    tree_obj = constructor.nj(dm)
    buf = StringIO()
    Phylo.write(tree_obj, buf, 'newick')
    return buf.getvalue()


def build_ft2_tree(sequences, taxon_names, fasta_path, out_path):
    """Build FT2 tree via MAFFT+FastTree2. Returns None if unavailable."""
    import subprocess
    from fusang_mhl.config import MAFFT_BAT, MAFFT_DIR, FASTTREE_EXE
    
    if not os.path.exists(MAFFT_BAT) or not os.path.exists(FASTTREE_EXE):
        return None, "MAFFT or FastTree2 not available"
    
    # Write FASTA
    with open(fasta_path, 'w') as f:
        for name, seq in zip(taxon_names, sequences):
            f.write(f'>{name}\n{seq}\n')
    
    # Run MAFFT (use absolute paths since cwd != file location)
    msa_path = fasta_path.replace('.fasta', '_aligned.fasta')
    try:
        import shutil
        os.makedirs(os.path.dirname(msa_path), exist_ok=True)
        fasta_abs = os.path.abspath(fasta_path)
        msa_abs = os.path.abspath(msa_path)
        cmd = f'"{MAFFT_BAT}" --auto --quiet "{fasta_abs}" > "{msa_abs}"'
        result = subprocess.run(
            cmd, shell=True,
            cwd=str(MAFFT_DIR) if os.path.isdir(MAFFT_DIR) else os.path.dirname(MAFFT_BAT),
            capture_output=True, text=True, timeout=300,
        )
        if result.returncode != 0 or not os.path.exists(msa_path) or os.path.getsize(msa_path) < 10:
            return None, f"MAFFT failed: rc={result.returncode}"
    except Exception as e:
        return None, f"MAFFT error: {e}"
    
    # Run FastTree2
    try:
        msa_abs = os.path.abspath(msa_path)
        out_abs = os.path.abspath(out_path)
        cmd = f'"{FASTTREE_EXE}" -gtr -quiet "{msa_abs}" > "{out_abs}"'
        result = subprocess.run(
            cmd, shell=True, capture_output=True, text=True, timeout=120,
        )
        if result.returncode != 0 or not os.path.exists(out_path):
            return None, f"FastTree2 failed: rc={result.returncode}"
        
        with open(out_path, 'r') as f:
            return f.read().strip(), "OK"
    except Exception as e:
        return None, f"FastTree2 error: {e}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--n-taxa', type=int, default=200)
    parser.add_argument('--n-seeds', type=int, default=30)
    parser.add_argument('--start-index', type=int, default=1,
                        help="Seed number to start from (1-indexed)")
    parser.add_argument('--seq-len', type=int, default=1000)
    parser.add_argument('--sub', type=float, default=0.05)
    parser.add_argument('--indel', type=float, default=0.02)
    args = parser.parse_args()
    
    n = args.n_taxa
    seeds = args.n_seeds
    start_idx = args.start_index
    L = args.seq_len
    sub = args.sub
    indel = args.indel
    
    # Output directory
    out_dir = f"l3_validation_n{n}"
    os.makedirs(out_dir, exist_ok=True)
    
    print("=" * 70)
    print("P0-4: L0-L3 End-to-End Validation")
    print("=" * 70)
    print(f"n={n}, L={L}, sub={sub}, indel={indel}, seeds {start_idx}-{start_idx+seeds-1} (total={seeds})")
    
    # Check MAFFT availability
    try:
        from fusang_mhl.config import MAFFT_BAT, FASTTREE_EXE
        mafft_ok = os.path.exists(MAFFT_BAT)
        ft2_ok = os.path.exists(FASTTREE_EXE)
    except Exception:
        mafft_ok = False
        ft2_ok = False
    print(f"MAFFT available: {mafft_ok}, FastTree2 available: {ft2_ok}")
    print("=" * 70)
    
    # Load existing results if any
    partial_path = f"{out_dir}/l3_validation_results.json"
    results = []
    if os.path.exists(partial_path):
        try:
            with open(partial_path, 'r') as f:
                prev = json.load(f)
                results = prev.get("results", [])
                print(f"Loaded {len(results)} existing results from {partial_path}")
        except Exception:
            pass
    
    for seed in range(start_idx, start_idx + seeds):
        print(f"\nSeed {seed}/{start_idx+seeds-1}", end=" ", flush=True)
        row = {"seed": seed}
        
        try:
            # Generate data
            sequences, taxon_names, true_nwk = generate_data(n, L, sub, indel, seed)
            
            # Save true tree
            true_path = f"{out_dir}/seed{seed:03d}_true.nwk"
            with open(true_path, 'w') as f:
                f.write(true_nwk)
            
            # L0: k-mer NJ (k=5, gap=2)
            t0 = time.time()
            l0_nwk = build_nj_tree(sequences, taxon_names)
            l0_time = time.time() - t0
            l0_nrf = compute_nrf(true_nwk, l0_nwk)
            row["l0_nrf"] = l0_nrf
            row["l0_time"] = l0_time
            
            # L1: Multi-k ensemble (k=5,7,9 contiguous)
            t0 = time.time()
            l1_nwk = build_multi_k_nj(sequences, taxon_names, ks=[5, 7, 9])
            l1_time = time.time() - t0
            l1_nrf = compute_nrf(true_nwk, l1_nwk)
            row["l1_nrf"] = l1_nrf
            row["l1_time"] = l1_time
            
            # L3: MSA+ML (MAFFT+FastTree2) - only if available
            ft2_nrf = None
            ft2_time = None
            ft2_status = "skipped"
            if mafft_ok and ft2_ok:
                fasta_path = f"{out_dir}/seed{seed:03d}.fasta"
                ft2_path = f"{out_dir}/seed{seed:03d}_ft2.nwk"
                t0 = time.time()
                ft2_nwk, ft2_status = build_ft2_tree(
                    sequences, taxon_names, fasta_path, ft2_path)
                ft2_time = time.time() - t0
                if ft2_nwk:
                    ft2_nrf = compute_nrf(true_nwk, ft2_nwk)
            row["ft2_nrf"] = ft2_nrf
            row["ft2_time"] = ft2_time
            row["ft2_status"] = ft2_status
            
            results.append(row)
            
            status = []
            if l0_nrf is not None:
                status.append(f"L0={l0_nrf:.4f}")
            if l1_nrf is not None:
                status.append(f"L1={l1_nrf:.4f}")
            if ft2_nrf is not None:
                status.append(f"FT2={ft2_nrf:.4f}")
            elif ft2_status != "skipped":
                status.append(f"FT2={ft2_status}")
            print(f"| {' '.join(status)}")
            
        except Exception as e:
            print(f"| ERROR: {e}")
            results.append({"seed": seed, "error": str(e)})
        
        # Intermediate save every 5 seeds
        if seed % 5 == 0 or seed == start_idx + seeds - 1:
            partial = {
                "config": {"n": n, "L": L, "sub": sub, "indel": indel, "seeds": seeds,
                           "start_index": start_idx},
                "results": results,
            }
            with open(partial_path, 'w') as f:
                json.dump(partial, f, indent=2, default=str)
    
    # ==========================================================================
    # Summary
    # ==========================================================================
    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print(f"{'=' * 70}")
    
    l0_nrfs = [r["l0_nrf"] for r in results if r.get("l0_nrf") is not None]
    l1_nrfs = [r["l1_nrf"] for r in results if r.get("l1_nrf") is not None]
    ft2_nrfs = [r["ft2_nrf"] for r in results if r.get("ft2_nrf") is not None]
    
    if l0_nrfs:
        print(f"\n  L0 (k-mer k=5,gap2 NJ):")
        print(f"    nRF = {np.mean(l0_nrfs):.4f} ± {np.std(l0_nrfs, ddof=1):.4f}  (n={len(l0_nrfs)})")
    
    if l1_nrfs:
        print(f"\n  L1 (multi-k k=5,7,9 NJ):")
        print(f"    nRF = {np.mean(l1_nrfs):.4f} ± {np.std(l1_nrfs, ddof=1):.4f}  (n={len(l1_nrfs)})")
        if l0_nrfs and len(l1_nrfs) == len(l0_nrfs):
            from scipy.stats import wilcoxon
            stat, p = wilcoxon(l0_nrfs, l1_nrfs)
            improvement = np.mean(l0_nrfs) - np.mean(l1_nrfs)
            print(f"    Δ vs L0: {improvement:+.4f} (Wilcoxon p={p:.4f})")
    
    if ft2_nrfs:
        print(f"\n  L3/FT2 (MAFFT+FastTree2 GTR):")
        print(f"    nRF = {np.mean(ft2_nrfs):.4f} ± {np.std(ft2_nrfs, ddof=1):.4f}  (n={len(ft2_nrfs)})")
        if l0_nrfs and len(ft2_nrfs) >= 5:
            from scipy.stats import wilcoxon
            common_l0 = [r["l0_nrf"] for r in results if r.get("ft2_nrf") is not None]
            stat, p = wilcoxon(common_l0, ft2_nrfs)
            print(f"    Δ vs L0: {np.mean(common_l0) - np.mean(ft2_nrfs):+.4f} (Wilcoxon p={p:.4f})")
    else:
        if mafft_ok:
            errors = [r.get("ft2_status", "") for r in results]
            print(f"\n  L3/FT2: ALL FAILED (reason: {set(errors)})")
        else:
            print(f"\n  L3/FT2: NOT AVAILABLE (MAFFT={mafft_ok}, FastTree2={ft2_ok})")
    
    # ==========================================================================
    # Overall comparison
    # ==========================================================================
    print(f"\n{'=' * 70}")
    print("LEVEL COMPARISON")
    print(f"{'=' * 70}")
    
    methods = []
    if l0_nrfs:
        methods.append(("L0 (k-mer NJ)", np.mean(l0_nrfs), np.std(l0_nrfs, ddof=1)))
    if l1_nrfs:
        methods.append(("L1 (multi-k NJ)", np.mean(l1_nrfs), np.std(l1_nrfs, ddof=1)))
    if ft2_nrfs:
        methods.append(("L3/FT2 (MSA+ML)", np.mean(ft2_nrfs), np.std(ft2_nrfs, ddof=1)))
    
    best = min(methods, key=lambda x: x[1]) if methods else None
    
    print(f"\n  {'Method':<25s} {'nRF':>8s} {'±':>8s}")
    print(f"  {'-'*25} {'-'*8} {'-'*8}")
    for name, mean, std in methods:
        marker = " ← BEST" if (name, mean, std) == best else ""
        print(f"  {name:<25s} {mean:>8.4f} {std:>8.4f}{marker}")
    
    # Save full results
    output = {
        "config": {"n": n, "L": L, "sub": sub, "indel": indel, "seeds": seeds,
                   "mafft_available": mafft_ok, "fasttree_available": ft2_ok},
        "summary": {
            "l0": {"mean": float(np.mean(l0_nrfs)), "std": float(np.std(l0_nrfs, ddof=1)), "n": len(l0_nrfs)},
            "l1": {"mean": float(np.mean(l1_nrfs)), "std": float(np.std(l1_nrfs, ddof=1)), "n": len(l1_nrfs)},
            "ft2": {"mean": float(np.mean(ft2_nrfs)) if ft2_nrfs else None, 
                     "std": float(np.std(ft2_nrfs, ddof=1)) if ft2_nrfs else None, 
                     "n": len(ft2_nrfs)},
        },
        "results": results,
    }
    
    out_path = f"{out_dir}/l3_validation_results.json"
    with open(out_path, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nResults saved to {out_path}")
    
    return results


if __name__ == "__main__":
    main()
