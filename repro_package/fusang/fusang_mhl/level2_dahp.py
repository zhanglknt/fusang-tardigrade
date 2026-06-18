"""
Level 2: DAHP-V2 centroid backbone + optional NNI refinement.

Wraps dahp_build_v2 from fusang_v4_dahp_v1.py.
For clusters that pass the boundary classifier (split decision),
refine using centroid-based MSA backbone within each cluster.
"""

import sys
import os
import numpy as np
from typing import Dict, List, Tuple, Set, Optional, Any

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from .config import L2_DEFAULTS
from .mlh_utils import Timer


def run_level2(
    fasta_path: str,
    clade_threshold: float = L2_DEFAULTS["clade_threshold"],
    k: int = L2_DEFAULTS["k"],
    gap: str = L2_DEFAULTS["gap"],
    enable_nni: bool = L2_DEFAULTS["enable_nni"],
    nni_max_rounds: int = L2_DEFAULTS["nni_max_rounds"],
    verbose: bool = True,
) -> Optional[str]:
    """Run Level 2: DAHP-V2 on a single cluster FASTA file.

    Args:
        fasta_path: Path to input FASTA (single cluster)
        clade_threshold: Distance threshold for clade cutting
        k, gap: k-mer parameters for distance
        enable_nni: Whether to run NNI refinement
        nni_max_rounds: Max NNI rounds
        verbose: Print progress

    Returns:
        Newick string of refined tree, or None on error
    """
    if verbose:
        print(f"[L2] DAHP-V2 on {fasta_path}", file=sys.stderr)

    with Timer("L2: DAHP-V2", verbose=verbose):
        from fusang_v4_dahp_v1 import dahp_build_v2
        nwk = dahp_build_v2(
            fasta_path,
            clade_threshold=clade_threshold,
            k=k,
            gap=gap,
        )

    if nwk is None:
        if verbose:
            print("  [L2] dahp_build_v2 returned None", file=sys.stderr)
        return None

    # Optional NNI refinement
    if enable_nni:
        if verbose:
            print("  [L2] Running NNI refinement...", file=sys.stderr)
        with Timer("L2: NNI", verbose=verbose):
            try:
                from fusang_v2 import nni_refine
                refined = nni_refine(nwk, max_rounds=nni_max_rounds)
                if refined:
                    nwk = refined
                    if verbose:
                        print("  [L2] NNI refinement OK", file=sys.stderr)
            except Exception as e:
                if verbose:
                    print(f"  [L2] NNI error (skipped): {e}", file=sys.stderr)

    return nwk


def run_level2_on_clusters(
    cluster_fastas: List[str],
    clade_thresholds: Optional[List[float]] = None,
    k: int = L2_DEFAULTS["k"],
    gap: str = L2_DEFAULTS["gap"],
    enable_nni: bool = False,
    verbose: bool = True,
) -> List[Optional[str]]:
    """Run Level 2 on multiple cluster FASTA files.

    Args:
        cluster_fastas: List of paths to cluster FASTA files
        clade_thresholds: Per-cluster threshold (or None for default)
        k, gap: k-mer parameters
        enable_nni: Whether to run NNI
        verbose: Print progress

    Returns:
        List of Newick strings (one per cluster, None for failures)
    """
    n_clusters = len(cluster_fastas)
    if clade_thresholds is None:
        clade_thresholds = [L2_DEFAULTS["clade_threshold"]] * n_clusters

    results = []
    for ci, (fa_path, thresh) in enumerate(zip(cluster_fastas, clade_thresholds)):
        if verbose:
            print(f"[L2] Cluster {ci+1}/{n_clusters} (threshold={thresh})",
                  file=sys.stderr)
        nwk = run_level2(
            fa_path,
            clade_threshold=thresh,
            k=k,
            gap=gap,
            enable_nni=enable_nni,
            verbose=verbose,
        )
        results.append(nwk)

    return results


def adaptive_clade_threshold(
    cluster_size: int,
    large_threshold: int = L2_DEFAULTS["large_cluster_threshold"],
    strict_thresh: float = L2_DEFAULTS["clade_threshold_large"],
    default_thresh: float = L2_DEFAULTS["clade_threshold"],
) -> float:
    """Get adaptive clade threshold based on cluster size.

    Larger clusters need stricter (smaller) thresholds to ensure
    adequate subdivision.
    """
    if cluster_size >= large_threshold:
        return strict_thresh
    return default_thresh
