"""
Level 3: MSA+ML within MRC (Minimum Reliable Cluster).

For clusters where pairwise distance < threshold (d < 0.01),
run full MSA (MAFFT) + ML (FastTree2) to get precise topology.
"""

import sys
import os
import tempfile
import subprocess
from typing import Dict, List, Tuple, Set, Optional, Any
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from .config import (
    MAFFT_BAT, FASTTREE_EXE, MAFFT_DIR, MAFFT_TMP,
    L3_DEFAULTS,
)
from .mlh_utils import Timer, ensure_dir, write_fasta, write_newick


def _run_mafft(fasta_path: str, output_msa: str, verbose: bool = True) -> bool:
    """Run MAFFT on a FASTA file. Returns True on success."""
    ensure_dir(str(MAFFT_TMP))
    tmp_input = MAFFT_TMP / Path(fasta_path).name
    try:
        import shutil
        shutil.copy2(fasta_path, str(tmp_input))
    except Exception:
        pass

    rel_path = Path(fasta_path).name
    cmd = f'"{MAFFT_BAT}" --auto --quiet {rel_path} > "{output_msa}"'
    try:
        result = subprocess.run(
            cmd, shell=True, cwd=str(MAFFT_DIR),
            capture_output=True, text=True, timeout=300,
        )
        if result.returncode != 0:
            if verbose:
                print(f"[L3:MAFFT] FAILED (rc={result.returncode}): {result.stderr[:200]}",
                      file=sys.stderr)
            return False
        # Verify output exists and is non-empty
        if not os.path.exists(output_msa) or os.path.getsize(output_msa) < 10:
            if verbose:
                print(f"[L3:MAFFT] Output empty: {output_msa}", file=sys.stderr)
            return False
        return True
    except subprocess.TimeoutExpired:
        if verbose:
            print("[L3:MAFFT] Timeout", file=sys.stderr)
        return False
    except Exception as e:
        if verbose:
            print(f"[L3:MAFFT] Error: {e}", file=sys.stderr)
        return False


def _run_fasttree2(msa_path: str, output_nwk: str, verbose: bool = True) -> Optional[str]:
    """Run FastTree2 on an MSA file. Returns Newick string or None."""
    cmd = f'"{FASTTREE_EXE}" -nt -gtr -quiet "{msa_path}" > "{output_nwk}"'
    try:
        result = subprocess.run(
            cmd, shell=True,
            capture_output=True, text=True, timeout=600,
        )
        if result.returncode != 0:
            if verbose:
                print(f"[L3:FT2] FAILED (rc={result.returncode})", file=sys.stderr)
            return None
        if not os.path.exists(output_nwk):
            if verbose:
                print("[L3:FT2] Output not found", file=sys.stderr)
            return None
        with open(output_nwk, "r", encoding="utf-8", errors="replace") as f:
            nwk = f.read().strip()
        if not nwk.endswith(";"):
            nwk += ";"
        return nwk
    except subprocess.TimeoutExpired:
        if verbose:
            print("[L3:FT2] Timeout", file=sys.stderr)
        return None
    except Exception as e:
        if verbose:
            print(f"[L3:FT2] Error: {e}", file=sys.stderr)
        return None


def should_run_l3(
    D: np.ndarray,
    cluster_indices: List[int],
    max_pairwise_dist: float = L3_DEFAULTS["max_pairwise_dist"],
) -> bool:
    """Check if a cluster needs Level 3 (MSA+ML).

    Returns True if mean pairwise distance < threshold
    (cluster is close enough that k-mer is unreliable).
    """
    if len(cluster_indices) <= L3_DEFAULTS["min_taxa"]:
        return False
    if len(cluster_indices) > L3_DEFAULTS["max_taxa"]:
        # Too large for MSA+ML — skip L3, stay at L2
        return False
    sub_D = D[np.ix_(cluster_indices, cluster_indices)]
    n = len(cluster_indices)
    if n <= 1:
        return False
    from scipy.spatial.distance import squareform
    pd = squareform(sub_D, checks=False)
    mean_d = float(np.mean(pd))
    return mean_d < max_pairwise_dist


def run_l3_on_cluster(
    sequences: List[str],
    taxon_names: List[str],
    output_dir: Optional[str] = None,
    verbose: bool = True,
) -> Optional[str]:
    """Run Level 3 (MAFFT + FastTree2) on a single cluster.

    Args:
        sequences: List of DNA sequences for this cluster
        taxon_names: List of taxon names (same order)
        output_dir: Directory for temp files (default: system temp)
        verbose: Print progress

    Returns:
        Newick string of ML tree, or None on failure
    """
    if output_dir is None:
        output_dir = os.path.join(tempfile.gettempdir(), "fusang_mhl_l3")
    ensure_dir(output_dir)

    # Write input FASTA
    input_fa = os.path.join(output_dir, "l3_input.fasta")
    msa_out = os.path.join(output_dir, "l3_alignment.fasta")
    nwk_out = os.path.join(output_dir, "l3_tree.nwk")

    write_fasta({n: s for n, s in zip(taxon_names, sequences)}, input_fa)

    if verbose:
        seq_len = len(sequences[0]) if sequences else '?'
        print(f"[L3] Running MAFFT on {len(sequences)} taxa (L={seq_len})",
              file=sys.stderr)

    with Timer("L3: MAFFT", verbose=verbose):
        ok = _run_mafft(input_fa, msa_out, verbose=verbose)
    if not ok:
        return None

    with Timer("L3: FastTree2", verbose=verbose):
        nwk = _run_fasttree2(msa_out, nwk_out, verbose=verbose)

    return nwk


def run_l3_on_clusters(
    cluster_data: List[Tuple[List[str], List[str], str]],
    max_pairwise_dist: float = L3_DEFAULTS["max_pairwise_dist"],
    verbose: bool = True,
) -> List[Optional[str]]:
    """Run Level 3 on multiple clusters (those that need it).

    Args:
        cluster_data: List of (sequences, names, output_dir) per cluster
        max_pairwise_dist: Threshold for triggering L3
        verbose: Print progress

    Returns:
        List of Newick strings (None for clusters that skipped L3 or failed)
    """
    results = []
    for ci, (seqs, names, out_dir) in enumerate(cluster_data):
        if verbose:
            print(f"[L3] Cluster {ci+1}/{len(cluster_data)} (n={len(names)})",
                  file=sys.stderr)

        # Quick distance check to decide if L3 is needed
        if len(seqs) >= 3:
            from .level0_kmer import compute_l0_distance
            try:
                D_small = compute_l0_distance(
                    seqs[:min(10, len(seqs))],
                    names[:min(10, len(names))],
                    verbose=False,
                )
                upper_idx = np.triu_indices(len(D_small), k=1)
                mean_d = float(np.mean(D_small[upper_idx]))
                if mean_d >= max_pairwise_dist * 2:
                    if verbose:
                        print(f"  [L3] Skipped (mean_d={mean_d:.4f} >= threshold, "
                              f"k-mer still reliable)",
                              file=sys.stderr)
                    results.append(None)
                    continue
            except Exception:
                pass  # Run L3 anyway if check fails

        nwk = run_l3_on_cluster(seqs, names, out_dir, verbose=verbose)
        results.append(nwk)

    return results
