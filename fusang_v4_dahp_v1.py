#!/usr/bin/env python3
"""
Fusang v4 - DAHP (Distance-Adaptive Hybrid Phylogenetics)
V1: k-mer NJ + selective MSA refinement
V2: backbone refinement with centroids (experimental)
V3: multi-k distance ensemble (k=5,7,9 contiguous average)
"""

import sys, os, argparse, time, tempfile, subprocess, shutil
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional

sys.setrecursionlimit(10000)  # BioPython NJ uses recursion; n=1000 needs >1000

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from kmer_distance import compute_kmer_distance_matrix
from io import StringIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

MAFFT_DIR = Path(r"d:/系统发育树项目/Fusang/bench_tools/mafft-win/mafft-win")
MAFFT_BAT = str(MAFFT_DIR / "mafft.bat")
MAFFT_TMP = MAFFT_DIR / "tmp"
FASTTREE_EXE = r"d:/系统发育树项目/Fusang/bench_tools/FastTree.exe"

def read_fasta(path):
    seqs, name = {}, None
    VALID = set('ACGTURYSWKMBDHVN.-*')
    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                # Sanitize: keep only valid IUPAC characters
                seqs[name].append(''.join(c for c in line.upper() if c in VALID))
    return {k: ''.join(v) for k, v in seqs.items()}

def resolve_gap_pattern(gap_str):
    """Convert gap shorthand to tuple pattern for spaced k-mer."""
    presets = {
        'gap2': (0, 1, 2, 5, 6),
        'gap3': (0, 1, 2, 6, 7),
        'contiguous': None,
    }
    if gap_str in presets:
        return presets[gap_str]
    raise ValueError(f"Unknown gap pattern: {gap_str}. Use: {list(presets.keys())}")

def build_nj(D, names):
    try:
        n = len(names)
        # BioPython DistanceMatrix requires lower-triangle format:
        # matrix[i] has len i+1 elements, matrix[i][j] = D[i][j] for j<=i
        lower = []
        for i in range(n):
            lower.append([float(D[i][j]) for j in range(i + 1)])
        dm = DistanceMatrix(names, lower)
        tc = DistanceTreeConstructor()
        tree = tc.nj(dm)
        out = StringIO()
        Phylo.write(tree, out, 'newick')
        return out.getvalue().strip()
    except Exception as e:
        print(f"[NJ] Error: {e}", file=sys.stderr)
        import traceback; traceback.print_exc()
        return None


# =========================================================
# Multi-k Distance Ensemble (DAHP-V3)
# =========================================================

def compute_multik_distance_matrix(sequences, taxon_names, ks=(5, 7, 9),
                                  gap_pattern=None, n_threads=4,
                                  fusion_method='average'):
    """
    Compute multi-k fused distance matrix.

    Parameters:
        sequences: list of sequence strings
        taxon_names: list of taxon names
        ks: tuple of k values to use (default: 5,7,9 contiguous)
        gap_pattern: gapped k-mer pattern (None for contiguous k-mers)
        n_threads: parallel threads for k-mer computation
        fusion_method: 'average' (recommended), 'weighted', or 'zscore'

    Returns:
        fused_dist: numpy array (n, n), fused distance matrix
    """
    n = len(sequences)
    Ds = {}
    print(f"[Multi-k] Computing distance matrices for k={ks}")

    for k in ks:
        D = compute_kmer_distance_matrix(
            sequences, taxon_names, k=k, metric='cosine',
            n_threads=n_threads, gap_pattern=gap_pattern
        )
        Ds[k] = D

    # Fuse distance matrices
    if fusion_method == 'average':
        fused = np.zeros((n, n), dtype=np.float64)
        for k in ks:
            fused += Ds[k]
        fused /= len(ks)
    elif fusion_method == 'weighted':
        fused = np.zeros((n, n), dtype=np.float64)
        total_w = 0
        for k in ks:
            w = k
            fused += w * Ds[k]
            total_w += w
        fused /= total_w
    else:
        raise ValueError(f"Unknown fusion method: {fusion_method}")

    fused = (fused + fused.T) / 2
    np.fill_diagonal(fused, 0.0)
    print(f"[Multi-k] Fused distance matrix ({n}x{n}), method={fusion_method}")

    return fused.astype(np.float32)


def dahp_build_v3(fasta_path, ks=(5, 7, 9), n_threads=4, fusion='average'):
    """
    DAHP-V3: Multi-k distance ensemble + NJ tree construction.
    Uses contiguous k-mers (no gap pattern) for maximum information diversity.

    Usage: python fusang_v4_dahp_v1.py input.fasta --v3
    """
    seqs = read_fasta(fasta_path)
    names = sorted(seqs.keys())
    sequences = [seqs[nm] for nm in names]
    n = len(names)

    print(f"\n{'='*60}")
    print(f"[DAHP-V3] Multi-k Ensemble (ks={ks}, fusion={fusion})")
    print(f"[DAHP-V3] {n} sequences, contiguous k-mers")
    print(f"{'='*60}")

    # Step 1: Multi-k distance ensemble (contiguous k-mers)
    t0 = time.time()
    fused_D = compute_multik_distance_matrix(
        sequences, names, ks=ks,
        gap_pattern=None, n_threads=n_threads,
        fusion_method=fusion
    )
    print(f"[DAHP-V3] Distance fusion: {time.time()-t0:.1f}s")

    # Step 2: Build NJ tree from fused distances
    t0 = time.time()
    nwk = build_nj(fused_D, names)
    print(f"[DAHP-V3] NJ tree: {time.time()-t0:.1f}s")

    return nwk


def run_mafft(fasta_path, output_msa):
    """
    Run MAFFT via mafft.bat (verified Windows approach).
    Copy input to MAFFT_DIR/tmp/, run with shell redirect to output file.
    Uses '>' redirect instead of capture_output because mafft.bat->bash.exe
    pipe chain drops stdout when captured via Python PIPE.
    """
    MAFFT_TMP.mkdir(parents=True, exist_ok=True)
    base_name = Path(fasta_path).name
    tmp_input = MAFFT_TMP / base_name
    shutil.copy2(fasta_path, tmp_input)
    rel_path = "tmp/" + base_name   # Use forward slash (bash-compatible)
    # Shell redirect: MSA output goes directly to file (more reliable on Windows)
    cmd = f'"{MAFFT_BAT}" --auto --quiet {rel_path} > "{output_msa}"'
    try:
        r = subprocess.run(
            cmd, shell=True, cwd=str(MAFFT_DIR),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=300
        )
        if r.returncode != 0:
            err = r.stderr.decode('utf-8', errors='replace')[:200] if r.stderr else ''
            print(f"    [MAFFT] rc={r.returncode}, stderr={err}", file=sys.stderr)
            return False
        if not os.path.exists(output_msa) or os.path.getsize(output_msa) < 10:
            err_full = r.stderr.decode('utf-8', errors='replace') if r.stderr else ''
            print(f"    [MAFFT] no output or too small ({os.path.getsize(output_msa) if os.path.exists(output_msa) else 0} bytes)", file=sys.stderr)
            print(f"    [MAFFT] cmd={cmd}", file=sys.stderr)
            print(f"    [MAFFT] stderr_full={err_full[:500]}", file=sys.stderr)
            return False
        return True
    except Exception as e:
        print(f"    [MAFFT] exception: {e}", file=sys.stderr)
        return False
    finally:
        try: tmp_input.unlink(missing_ok=True)
        except: pass

def run_ft2(msa_path):
    """Run FastTree2. Use shell=True + redirect (verified approach)."""
    try:
        tmp = tempfile.gettempdir()
        pid = os.getpid()
        ascii_msa = os.path.join(tmp, f"ft2_{pid}.msa")
        ascii_out = os.path.join(tmp, f"ft2_{pid}.nwk")
        shutil.copy2(msa_path, ascii_msa)
        cmd = f'"{FASTTREE_EXE}" -nt -gtr -quiet "{ascii_msa}" > "{ascii_out}"'
        try:
            result = subprocess.run(cmd, shell=True, capture_output=True, timeout=600)
            if result.returncode != 0 or not os.path.exists(ascii_out):
                err = result.stderr.decode('utf-8', errors='replace')[:200] if result.stderr else ''
                print(f"    [FT2] rc={result.returncode}, stderr={err}", file=sys.stderr)
                return None
            with open(ascii_out, 'r') as f:
                nwk = f.read().strip()
            return nwk if nwk else None
        finally:
            for fp in [ascii_msa, ascii_out]:
                try: os.remove(fp)
                except: pass
    except Exception as e:
        print(f"    [FT2] exception: {e}", file=sys.stderr)
        return None

def refine(seqs_dict, label=""):
    tmp = tempfile.gettempdir()
    pid = os.getpid()
    fa = os.path.join(tmp, f"dahp_{pid}_{label}.fasta")
    ms = os.path.join(tmp, f"dahp_{pid}_{label}.msa")
    try:
        with open(fa, 'w') as f:
            for n, s in seqs_dict.items(): f.write(f">{n}\n{s}\n")
        print(f"    MAFFT: aligning {len(seqs_dict)} seqs...")
        if not run_mafft(fa, ms):
            print(f"    MAFFT failed")
            return None
        msa_size = os.path.getsize(ms) if os.path.exists(ms) else 0
        print(f"    MAFFT ok ({msa_size} bytes)")
        print(f"    FastTree2: building tree...")
        nwk = run_ft2(ms)
        if nwk:
            print(f"    FastTree2 ok ({len(nwk)} chars)")
        else:
            print(f"    FastTree2 failed")
        return nwk
    finally:
        for fp in [fa, ms]:
            try: os.remove(fp)
            except: pass

def find_clade_by_taxa(tree, target_taxa):
    """Find the clade in tree whose terminal names equal target_taxa (as set)."""
    target = set(target_taxa)
    found = []
    def walk(node):
        if not hasattr(node, 'clades') or not node.clades:
            return
        leaves = {l.name for l in node.get_terminals()}
        if leaves == target:
            found.append(node)
            return
        for c in node.clades:
            walk(c)
    walk(tree.root)
    return found[0] if found else None

def merge_refined_into_tree(orig_nwk, refined_list):
    """
    Merge refined subtrees back into the original NJ tree.
    refined_list: list of (taxa_set, refined_newick) tuples.
    Returns merged Newick string.
    """
    if not refined_list:
        return orig_nwk
    tree = Phylo.read(StringIO(orig_nwk), "newick")
    for taxa_set, rnwk in refined_list:
        clade = find_clade_by_taxa(tree, taxa_set)
        if not clade:
            print(f"    [MERGE] WARNING: clade with {len(taxa_set)} taxa not found in original tree, skipping")
            continue
        # Parse refined tree and replace clade's structure
        rtree = Phylo.read(StringIO(rnwk), "newick")
        rroot = rtree.root
        # Replace clade's clades with refined clades
        # Keep the original clade's branch length
        orig_branch = clade.branch_length
        clade.clades = rroot.clades
        clade.name = rroot.name
        # Restore original branch length
        clade.branch_length = orig_branch
        print(f"    [MERGE] Replaced clade ({len(taxa_set)} taxa)")
    out = StringIO()
    Phylo.write(tree, out, "newick")
    return out.getvalue().strip()

def dahp_build(fasta_path, threshold=0.5, max_msa=100, k=5, gap='gap2'):
    t0 = time.time()
    print(f"=== DAHP-V1 ===")
    seqs = read_fasta(fasta_path)
    names = sorted(seqs.keys())
    n = len(names)
    print(f"Read {n} seqs")

    gap_tuple = resolve_gap_pattern(gap)
    print(f"Computing k-mer matrix (k={k}, gap={gap_tuple})...")
    D = compute_kmer_distance_matrix(
        [seqs[nm] for nm in names],
        taxon_names=names,
        k=k,
        gap_pattern=gap_tuple,
        n_threads=1
    )
    
    print("Building NJ tree...")
    nwk = build_nj(D, names)
    if not nwk:
        print("NJ failed")
        return None
    print(f"NJ: {nwk[:50]}...")
    
    print("Identifying clades for refinement...")
    tree = Phylo.read(StringIO(nwk), "newick")
    name_to_idx = {names[i]: i for i in range(n)}
    
    # Step 1: collect all candidate clades (don't refine yet)
    candidates = []   # (n_leaves, taxa_set, max_d)
    def collect(node):
        if not hasattr(node, 'clades') or not node.clades:
            return
        if len(node.clades) >= 2:
            leaves = [l.name for l in node.get_terminals()]
            if 4 <= len(leaves) <= max_msa:
                idx = [name_to_idx[l] for l in leaves]
                mx = max(D[i,j] for i in idx for j in idx if i < j)
                if mx > threshold:
                    candidates.append((len(leaves), frozenset(leaves), mx))
        for c in node.clades:
            collect(c)
    
    collect(tree.root)
    
    # Step 2: greedily select non-overlapping clades (largest first)
    candidates.sort(key=lambda x: -x[0])   # descending by size
    selected = []
    covered = set()
    for size, taxa, mx in candidates:
        if taxa & covered:   # overlaps with already selected
            continue
        selected.append((taxa, mx))
        covered.update(taxa)
    
    print(f"  {len(candidates)} candidate clades, {len(selected)} selected (non-overlapping)")
    
    # Step 3: refine selected clades
    refined_list = []
    for i, (taxa, mx) in enumerate(selected):
        print(f"  Refining clade {i+1}/{len(selected)}: {len(taxa)} leaves, max_d={mx:.4f}")
        seqs_sub = {n: seqs[n] for n in taxa}
        rnwk = refine(seqs_sub, label=f"C{i}")
        if rnwk:
            refined_list.append((taxa, rnwk))
        else:
            print(f"    FAILED")
    
    print(f"Refined {len(refined_list)} clades")
    
    # Merge refined subtrees back into the NJ tree
    if refined_list:
        print("Merging refined subtrees into NJ tree...")
        nwk = merge_refined_into_tree(nwk, refined_list)
        print(f"  Merged {len(refined_list)} subtrees")
    
    print(f"Done in {time.time()-t0:.1f}s")
    return nwk

# ============================================================
# DAHP-V2: Refine inter-clade (backbone) topology using centroids
# ============================================================
# Key insight:
#   - Within clades: k-mer distances are short → k-mer NJ is accurate → no refinement needed
#   - Between clades: centroid distances can be large → k-mer unreliable → refine with MSA
#   - MSA is only done on centroid sequences (<< n), not all sequences

def cut_tree_into_clades(root, D, name_to_idx, threshold):
    """
    Cut NJ tree into clades where max pairwise k-mer distance <= threshold.

    Returns list of frozensets (each a set of taxa names forming one clade).
    Clades are sorted by size descending.
    """
    clades = []

    def cut(node):
        terminals = node.get_terminals()
        if len(terminals) <= 1:
            clades.append(frozenset(t.name for t in terminals))
            return
        idx = [name_to_idx[t.name] for t in terminals]
        max_d = max(D[i][j] for i in idx for j in idx if i < j)
        if max_d <= threshold:
            clades.append(frozenset(t.name for t in terminals))
        else:
            for child in node.clades:
                cut(child)

    cut(root)
    clades.sort(key=lambda c: -len(c))
    return clades


def pick_centroid(taxa_set, D, names):
    """Pick sequence with minimum average k-mer distance to others in taxa_set."""
    name_to_idx = {names[i]: i for i in range(len(names))}
    idx = [name_to_idx[t] for t in taxa_set]
    best = min(idx, key=lambda i: sum(D[i][j] for j in idx if j != i) / max(len(idx) - 1, 1))
    return names[best]


def graft_subtrees(backbone_nwk, clade_info_list):
    """
    Replace centroid leaves in backbone tree with internal clade subtrees.

    Maintains strict binary tree structure by inserting an intermediate
    bifurcating node when grafting multi-leaf subtrees.

    clade_info_list: list of dicts with 'centroid' and 'internal_nwk' keys.
    Returns Newick string of the grafted tree.
    """
    from Bio.Phylo.Newick import Clade

    backbone = Phylo.read(StringIO(backbone_nwk), 'newick')
    centroid_to_sub = {ci['centroid']: ci['internal_nwk'].rstrip(';') for ci in clade_info_list}

    def replace_leaves(node):
        new_children = []
        for child in node.clades:
            if not child.clades and child.name in centroid_to_sub:
                # This leaf is a centroid -> replace with internal subtree
                sub_nwk_str = centroid_to_sub[child.name] + ';'
                sub_tree = Phylo.read(StringIO(sub_nwk_str), 'newick')
                bl = child.branch_length or 0

                if sub_tree.root.clades:
                    # Multi-leaf subtree: maintain binary tree by inserting
                    # an intermediate node that holds the subtree's root clades
                    root_clades = sub_tree.root.clades
                    if len(root_clades) == 1:
                        # Subtree root has single child (rare): attach directly
                        sc = root_clades[0]
                        sc.branch_length = (sc.branch_length or 0) + bl
                        new_children.append(sc)
                    else:
                        # Create intermediate bifurcating node:
                        # Take first clade as left child, rest as right subtree
                        left = root_clades[0]
                        left.branch_length = (left.branch_length or 0) + bl / 2

                        if len(root_clades) == 2:
                            right = root_clades[1]
                            right.branch_length = (right.branch_length or 0) + bl / 2
                        else:
                            # 3+ clades: right side is another intermediate node
                            right = Clade(branch_length=bl / 2)
                            right.clades = root_clades[1:]

                        intermediate = Clade(branch_length=0)
                        intermediate.clades = [left, right]
                        new_children.append(intermediate)
                else:
                    # Size-1 clade: keep as leaf with branch length
                    leaf = Clade(name=sub_tree.root.name or child.name,
                                 branch_length=bl)
                    new_children.append(leaf)
            else:
                replace_leaves(child)
                new_children.append(child)
        node.clades = new_children

    replace_leaves(backbone.root)
    out = StringIO()
    Phylo.write(backbone, out, 'newick')
    return out.getvalue().strip()


def dahp_build_v2(fasta_path, clade_threshold=0.3, k=5, gap='gap2'):
    """
    DAHP-V2: Refine inter-clade (backbone) topology using representative sequences.

    Pipeline:
    1. k-mer distance matrix -> NJ tree
    2. Cut NJ tree into clades (max internal k-mer distance <= clade_threshold)
    3. For each clade: internal NJ tree from k-mer submatrix (accurate for short distances)
    4. Pick centroid for each clade
    5. Backbone: MSA on centroids -> ML tree (FastTree2)
       - Only MSA on centroid_count sequences, NOT on all n sequences
    6. Graft internal trees onto backbone at centroid positions

    Args:
        fasta_path: input FASTA file
        clade_threshold: max pairwise k-mer distance within a clade (default 0.3)
        k: k-mer size (default 5)
        gap: gap pattern (default 'gap2')
    Returns:
        Newick string of final tree, or None on failure
    """
    t0 = time.time()
    seqs = read_fasta(fasta_path)
    names = sorted(seqs.keys())
    n = len(names)
    print(f"=== DAHP-V2 (n={n}, clade_threshold={clade_threshold}) ===")

    gap_tuple = resolve_gap_pattern(gap)

    # Step 1: k-mer distance matrix
    print("Computing k-mer distance matrix...")
    D = compute_kmer_distance_matrix(
        [seqs[nm] for nm in names], taxon_names=names,
        k=k, gap_pattern=gap_tuple, n_threads=1
    )

    # Step 2: NJ tree
    print("Building NJ tree...")
    nwk = build_nj(D, names)
    if not nwk:
        print("NJ failed")
        return None

    tree = Phylo.read(StringIO(nwk), 'newick')
    name_to_idx = {names[i]: i for i in range(n)}

    # Step 3: Cut into clades
    print(f"Cutting tree into clades (max_d <= {clade_threshold})...")
    clades = cut_tree_into_clades(tree.root, D, name_to_idx, clade_threshold)
    nc = len(clades)
    sizes = [len(c) for c in clades]
    print(f"  {nc} clades, sizes: {sizes[:10]}{'...' if nc > 10 else ''}")
    print(f"  min_size={min(sizes)}, max_size={max(sizes)}, avg={sum(sizes)/nc:.1f}")

    # Step 4: Internal trees + centroids
    print("Building internal clade trees and picking centroids...")
    clade_info = []
    for taxa_set in clades:
        clade_names = sorted(taxa_set)
        m = len(clade_names)

        if m == 1:
            clade_info.append({
                'centroid': clade_names[0],
                'internal_nwk': clade_names[0] + ';',
                'taxa': taxa_set,
                'size': m,
            })
            continue

        # Sub-distance matrix for this clade
        idx = [name_to_idx[t] for t in clade_names]
        sub_D = [[D[i][j] for j in idx] for i in idx]

        # Internal NJ tree (k-mer is accurate for short distances within clade)
        internal_nwk = build_nj(sub_D, clade_names)
        if not internal_nwk:
            internal_nwk = ','.join(clade_names) + ';'

        # Pick centroid (sequence with min avg distance to others)
        centroid = pick_centroid(taxa_set, D, names)

        clade_info.append({
            'centroid': centroid,
            'internal_nwk': internal_nwk,
            'taxa': taxa_set,
            'size': m,
        })

    centroids = [c['centroid'] for c in clade_info]
    print(f"  {len(centroids)} centroids selected")
    for ci in clade_info[:5]:
        print(f"    {ci['centroid']}: {ci['size']} seqs")
    if len(clade_info) > 5:
        print(f"    ... and {len(clade_info)-5} more")

    # Step 5: Backbone tree from MSA on centroids
    print(f"\nBackbone: MSA on {len(centroids)} centroids (NOT {n} sequences!)")
    backbone_nwk = None

    if len(centroids) >= 4:
        # Write centroid sequences to temp FASTA
        tmp = tempfile.gettempdir()
        pid = os.getpid()
        tmp_fasta = os.path.join(tmp, f"dahp_v2_cents_{pid}.fasta")
        tmp_msa = os.path.join(tmp, f"dahp_v2_cents_{pid}.msa")
        try:
            with open(tmp_fasta, 'w', encoding='utf-8') as f:
                for c in centroids:
                    f.write(f">{c}\n{seqs[c]}\n")

            t_mafft = time.time()
            print(f"  MAFFT aligning {len(centroids)} centroids...")
            if run_mafft(tmp_fasta, tmp_msa):
                print(f"  MAFFT ok ({time.time()-t_mafft:.1f}s)")

                t_ft2 = time.time()
                print(f"  FastTree2 building backbone...")
                backbone_nwk = run_ft2(tmp_msa)
                if backbone_nwk:
                    print(f"  FT2 ok ({time.time()-t_ft2:.1f}s)")
                else:
                    print(f"  FT2 failed, falling back to k-mer backbone")
            else:
                print(f"  MAFFT failed, falling back to k-mer backbone")
        finally:
            for fp in [tmp_fasta, tmp_msa]:
                try: os.remove(fp)
                except: pass

    if backbone_nwk is None:
        # Fallback: build backbone from k-mer distances between centroids
        if len(centroids) >= 2:
            print("  Building k-mer backbone (NJ on centroids)...")
            cent_idx = [name_to_idx[c] for c in centroids]
            cent_D = [[D[i][j] for j in cent_idx] for i in cent_idx]
            backbone_nwk = build_nj(cent_D, centroids)
        else:
            # Single clade, no backbone
            backbone_nwk = clade_info[0]['internal_nwk']

    if not backbone_nwk:
        print("Backbone construction failed!")
        return None

    # Step 6: Graft internal trees onto backbone
    print("Grafting internal clade trees onto backbone...")
    final_nwk = graft_subtrees(backbone_nwk, clade_info)

    # Verify leaf count
    try:
        final_tree = Phylo.read(StringIO(final_nwk), 'newick')
        n_leaves = len(final_tree.get_terminals())
        print(f"  Final tree: {n_leaves} leaves (expected {n})")
        if n_leaves != n:
            print(f"  WARNING: leaf count mismatch!")
    except Exception as e:
        print(f"  WARNING: could not parse final tree: {e}")

    print(f"\nDone in {time.time()-t0:.1f}s")
    return final_nwk


def main():
    parser = argparse.ArgumentParser(description="Fusang DAHP (Distance-Adaptive Hybrid Phylogenetics)")
    parser.add_argument("fasta", help="Input FASTA")
    parser.add_argument("--v2", action="store_true", help="Use DAHP-V2 (backbone refinement with centroids)")
    parser.add_argument("--v3", action="store_true", help="Use DAHP-V3 (multi-k distance ensemble, k=5,7,9 contiguous)")
    parser.add_argument("--threshold", type=float, default=0.5, help="V1: clade max_d threshold (default 0.5)")
    parser.add_argument("--clade-threshold", type=float, default=0.3, help="V2: max pairwise k-mer distance within clade (default 0.3)")
    parser.add_argument("--max-msa", type=int, default=100, help="V1: max sequences per MSA (default 100)")
    parser.add_argument("--k", type=int, default=5)
    parser.add_argument("--gap", default="gap2")
    parser.add_argument("--ks", default="5,7,9", help="V3: k-mer sizes for ensemble (comma-separated, default: 5,7,9)")
    parser.add_argument("--fusion", choices=["average", "weighted"], default="average", help="V3: distance fusion method")
    parser.add_argument("--output", default=None)
    args = parser.parse_args()

    if args.v3:
        ks = tuple(int(x) for x in args.ks.split(','))
        suffix = '_dahp_v3.nwk'
        nwk = dahp_build_v3(args.fasta, ks=ks, n_threads=4, fusion=args.fusion)
    elif args.v2:
        suffix = '_dahp_v2.nwk'
        nwk = dahp_build_v2(args.fasta, clade_threshold=args.clade_threshold, k=args.k, gap=args.gap)
    else:
        suffix = '_dahp_v1.nwk'
        nwk = dahp_build(args.fasta, threshold=args.threshold, max_msa=args.max_msa, k=args.k, gap=args.gap)

    output = args.output or args.fasta.replace('.fasta', suffix).replace('.fas', suffix)
    if nwk:
        with open(output, 'w') as f:
            f.write(nwk + "\n")
        print(f"Output: {output}")
    else:
        print("Failed")
        sys.exit(1)

if __name__ == "__main__":
    main()
