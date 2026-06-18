#!/usr/bin/env python3
"""
Fusang MHL (Multi-Level Hierarchical) Phylogenetic Tree Builder.

Main entry point for the 4-level hierarchical phylogenetic inference:
  Level 0: k-mer cosine distance clustering → backbone NJ
  Level 1: multi-k ensemble + boundary classifier decision
  Level 2: DAHP-V2 centroid backbone + optional NNI
  Level 3: MAFFT + FastTree2 within MRC (Minimum Reliable Cluster)
  Merge:   Bridge Taxa constrained NJ supertree assembly

Usage:
    python fusang_mhl.py input.fasta -o output.nwk [--verbose]

    # With custom parameters
    python fusang_mhl.py input.fasta -o output.nwk --k 5 --gap gap2 --threads 4

    # Debug mode (show per-level breakdown)
    python fusang_mhl.py input.fasta -o output.nwk --verbose --debug
"""

import sys
import os
import argparse
import json
from typing import Dict, List, Optional, Any

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fusang_mhl.config import (
    L0_THRESHOLDS, L0_K, L0_GAP, L0_DISTANCE_METHOD,
    L1_DEFAULTS, L2_DEFAULTS, L3_DEFAULTS, MERGE_DEFAULTS,
    TEMP_DIR, BOUNDARY_MODEL_DIR,
)
from fusang_mhl.mlh_utils import (
    Timer, ProgressReporter, setup_logger, ensure_dir,
    write_fasta, write_newick, read_fasta_simple,
    check_leaf_completeness,
)
from fusang_mhl.level0_kmer import run_level0, get_l0_config, cluster_taxa
from fusang_mhl.level1_multik import run_level1, extract_cluster_features
from fusang_mhl.boundary_classifier import BoundaryClassifier
from fusang_mhl.level2_dahp import run_level2, run_level2_on_clusters, adaptive_clade_threshold
from fusang_mhl.level3_msa_ml import run_l3_on_cluster, should_run_l3
from fusang_mhl.merger import merge_level, bottom_up_merge, select_bridge_taxa
from fusang_mhl.level0_kmer import compute_l0_distance

import getpass
USER = getpass.getuser()


def parse_args():
    parser = argparse.ArgumentParser(
        description="Fusang MHL: Multi-Level Hierarchical Phylogenetic Inference",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("input", type=str, help="Input FASTA file")
    parser.add_argument("-o", "--output", type=str, default=None,
                        help="Output Newick file (default: input.nwk)")
    parser.add_argument("--k", type=int, default=L0_K,
                        help=f"k-mer size (default: {L0_K})")
    parser.add_argument("--gap", type=str, default=L0_GAP,
                        help=f"Gap pattern: gap2, gap3, contiguous (default: {L0_GAP})")
    parser.add_argument("--threads", type=int, default=4,
                        help="Number of threads for k-mer (default: 4)")
    parser.add_argument("--verbose", action="store_true", default=True)
    parser.add_argument("--debug", action="store_true",
                        help="Show detailed per-level timing")
    parser.add_argument("--no-l2", action="store_true",
                        help="Skip Level 2 (DAHP-V2)")
    parser.add_argument("--no-l3", action="store_true",
                        help="Skip Level 3 (MSA+ML)")
    parser.add_argument("--boundary-model", type=str, default=None,
                        help="Path to pre-trained boundary classifier model")
    parser.add_argument("--output-dir", type=str, default=None,
                        help="Directory for intermediate files (default: TEMP)")
    return parser.parse_args()


def load_boundary_classifier(model_path: Optional[str] = None):
    """Load pre-trained boundary classifier, or return None."""
    if model_path and os.path.exists(model_path):
        try:
            from fusang_mhl.boundary_classifier import BoundaryClassifier
            clf = BoundaryClassifier()
            clf.load(model_path)
            return clf
        except Exception as e:
            print(f"[MHL] Warning: could not load boundary model: {e}",
                  file=sys.stderr)
    # Try default path
    if BOUNDARY_MODEL_DIR and os.path.exists(BOUNDARY_MODEL_DIR):
        for fname in ["boundary_rf.pkl", "boundary_rf_best.pkl"]:
            path = os.path.join(BOUNDARY_MODEL_DIR, fname)
            if os.path.exists(path):
                try:
                    from fusang_mhl.boundary_classifier import BoundaryClassifier
                    clf = BoundaryClassifier()
                    clf.load(path)
                    return clf
                except Exception:
                    continue
    return None


def heuristic_split_decision(
    D: np.ndarray,
    cluster_indices: List[int],
    cluster_names: List[str],
    n_total: int,
    min_cluster_size: int = L1_DEFAULTS["min_cluster_size"],
) -> bool:
    """Heuristic fallback when no boundary classifier is available.

    Split if:
      1. Cluster is large enough (>= min_cluster_size)
      2. Mean pairwise distance suggests multiple subgroups
      3. Distance distribution has high variance (multiple modes)

    Tuned 2026-06-16:
      - n_total <= 200: never split (NJ is accurate enough)
      - n_total > 200: use conservative thresholds
    """
    # Never split when total taxa is small — NJ is accurate and fast
    if n_total <= 200:
        return False

    n_c = len(cluster_indices)
    if n_c < min_cluster_size:
        return False

    sub_D = D[np.ix_(cluster_indices, cluster_indices)]
    n = len(cluster_indices)
    upper_idx = np.triu_indices(n, k=1)
    pairwise = sub_D[upper_idx]

    mean_d = float(np.mean(pairwise))
    std_d = float(np.std(pairwise))

    if mean_d < 0.005:
        return False  # Very close — stop here
    if mean_d > 0.2:
        return False  # Very diverse — k-mer is sufficient
    if std_d < 0.01:
        return False  # Uniform distances — likely single group

    # Use stricter CV threshold for large n
    cv = std_d / max(mean_d, 1e-6)
    min_cv = 0.5 if n_total > 500 else 0.4
    return cv > min_cv and n_c >= 3 * min_cluster_size


def run_fusang_mhl(
    fasta_path: str,
    output_path: Optional[str] = None,
    k: int = L0_K,
    gap: str = L0_GAP,
    metric: str = L0_DISTANCE_METHOD,
    n_threads: int = 4,
    verbose: bool = True,
    debug: bool = False,
    no_l2: bool = False,
    no_l3: bool = False,
    boundary_model: Optional[str] = None,
    output_dir: Optional[str] = None,
) -> Dict[str, Any]:
    """Run the full Fusang MHL pipeline.

    Args:
        fasta_path: Input FASTA file
        output_path: Output Newick file (None = auto)
        k: k-mer size
        gap: Gap pattern
        metric: Distance metric
        n_threads: Threads for k-mer computation
        verbose: Print progress
        debug: Show detailed timing
        no_l2: Skip Level 2
        no_l3: Skip Level 3
        boundary_model: Path to boundary classifier model
        output_dir: Intermediate files directory

    Returns:
        dict with 'tree', 'stats', 'hierarchy' keys
    """
    import tempfile

    # Setup
    if output_dir is None:
        output_dir = os.path.join(TEMP_DIR, f"run_{os.getpid()}")
    ensure_dir(output_dir)

    if output_path is None:
        output_path = os.path.splitext(fasta_path)[0] + ".nwk"

    if verbose:
        print(f"[MHL] Fusang MHL v1.0", file=sys.stderr)
        print(f"[MHL] Input: {fasta_path}", file=sys.stderr)
        print(f"[MHL] Output: {output_path}", file=sys.stderr)

    # Read input
    with Timer("Read FASTA", verbose=debug):
        seqs_dict = read_fasta_simple(fasta_path)
    taxon_names = list(seqs_dict.keys())
    sequences = [seqs_dict[n] for n in taxon_names]
    n = len(taxon_names)
    if verbose:
        print(f"[MHL] Loaded {n} taxa, L={len(sequences[0]) if sequences else '?'}",
              file=sys.stderr)

    if n < 4:
        if verbose:
            print("[MHL] Too few taxa for tree building", file=sys.stderr)
        return {"tree": None, "stats": {"error": "too few taxa"}, "hierarchy": []}

    # Load boundary classifier
    classifier = load_boundary_classifier(boundary_model)
    if classifier and verbose:
        print("[MHL] Boundary classifier loaded", file=sys.stderr)
    elif verbose:
        print("[MHL] No boundary classifier — using heuristic split decisions",
              file=sys.stderr)

    # ================================================================
    # Level 0: k-mer distance clustering
    # ================================================================
    D, clusters, backbone_nwk, centroids, cluster_stats = run_level0(
        sequences, taxon_names,
        k=k, gap=gap, metric=metric,
        n_threads=n_threads, verbose=verbose,
    )

    l0_config = get_l0_config(n)
    n_clusters = len(clusters)

    if verbose:
        print(f"[MHL] L0: {n_clusters} clusters", file=sys.stderr)

    # ================================================================
    # Level 1+: Process each L0 cluster
    # ================================================================
    hierarchy = []  # List of (level, cluster_indices, newick_or_None)
    level0_entry = []
    for ci, cluster in enumerate(clusters):
        level0_entry.append((cluster, None))
    hierarchy.append(level0_entry)

    all_subtree_results = []  # (cluster_indices, newick)

    for ci, cluster in enumerate(clusters):
        n_c = len(cluster)
        cluster_names = [taxon_names[i] for i in cluster]
        cluster_seqs = [sequences[i] for i in cluster]

        if verbose:
            print(f"\n[MHL] --- Cluster {ci+1}/{n_clusters} (n={n_c}) ---",
                  file=sys.stderr)

        # Decide: split further or stop?
        if n_c <= L1_DEFAULTS["min_cluster_size"]:
            if verbose:
                print(f"[MHL] Cluster too small for further splitting", file=sys.stderr)
            # Build simple NJ tree for this cluster
            from fusang_v4_dahp_v1 import build_nj
            sub_D = D[np.ix_(cluster, cluster)]
            nwk = build_nj(sub_D, cluster_names)
            all_subtree_results.append((cluster, nwk))
            continue

        # Use boundary classifier or heuristic
        should_split = False
        if classifier is not None:
            # For classifier: need multi-k distances for features
            try:
                from fusang_mhl.level1_multik import compute_l1_distance
                sub_D_multi = compute_l1_distance(
                    cluster_seqs, cluster_names, verbose=False,
                )
                features = extract_cluster_features(
                    sub_D_multi, None, None, None,
                    cluster_seqs, list(range(n_c)), 0,
                    None, n_total=n,
                    parent_size=n_c,
                    sibling_sizes=[],
                    current_level=1,
                    ancestor_sizes=[n_c],
                    verbose=False,
                )
                if features is not None:
                    prob = classifier.predict_proba(np.array([features]))[0][1]
                    should_split = prob > L1_DEFAULTS["split_threshold"]
                    if debug:
                        print(f"[MHL]   P(split) = {prob:.3f}",
                              file=sys.stderr)
            except Exception:
                should_split = heuristic_split_decision(
                    D, cluster, cluster_names, n,
                )
        else:
            should_split = heuristic_split_decision(D, cluster, cluster_names, n)

        if not should_split:
            if verbose:
                print(f"[MHL] Decision: STOP (build NJ for this cluster)",
                      file=sys.stderr)
            from fusang_v4_dahp_v1 import build_nj
            sub_D = D[np.ix_(cluster, cluster)]
            nwk = build_nj(sub_D, cluster_names)
            all_subtree_results.append((cluster, nwk))
            continue

        # Split: sub-cluster this L0 group
        if verbose:
            print(f"[MHL] Decision: SPLIT → sub-clustering", file=sys.stderr)

        # Sub-cluster using hierarchical clustering on sub-distance matrix
        sub_D = D[np.ix_(cluster, cluster)]
        n_sub_target = max(2, min(n_c // 30, 8))  # adaptive sub-clusters
        sub_clusters = cluster_taxa(
            sub_D, cluster_names,
            max_group_size=30,
            target_groups=n_sub_target,
        )

        # Convert sub-clusters from local to global indices
        sub_clusters_global = [[cluster[j] for j in sc] for sc in sub_clusters]

        if verbose:
            print(f"[MHL]   Sub-clusters: {len(sub_clusters_global)} "
                  + ", ".join(f"n={len(sc)}" for sc in sub_clusters_global),
                  file=sys.stderr)

        # ================================================================
        # Process each sub-cluster (L2/L3 decisions)
        # ================================================================
        sub_tree_entries = []

        for si, subcluster in enumerate(sub_clusters_global):
            sc_names = [taxon_names[i] for i in subcluster]
            sc_seqs = [sequences[i] for i in subcluster]
            sub_fa = os.path.join(output_dir, f"cluster_{ci}_sub{si}.fasta")
            write_fasta({n: s for n, s in zip(sc_names, sc_seqs)}, sub_fa)

            # Check if L3 is needed (close distance)
            if not no_l3 and should_run_l3(D, subcluster):
                if verbose:
                    print(f"[MHL]   Sub {si+1} (n={len(subcluster)}): close → L3",
                          file=sys.stderr)
                nwk = run_l3_on_cluster(
                    sc_seqs, sc_names,
                    output_dir=os.path.join(output_dir, f"l3_{ci}_{si}"),
                    verbose=verbose,
                )
                if nwk is None:
                    from fusang_v4_dahp_v1 import build_nj
                    sc_D = D[np.ix_(subcluster, subcluster)]
                    nwk = build_nj(sc_D, sc_names)
                sub_tree_entries.append((subcluster, nwk))
            elif not no_l2 and len(subcluster) >= L2_DEFAULTS["min_clade_size"]:
                if verbose:
                    print(f"[MHL]   Sub {si+1} (n={len(subcluster)}): → L2 DAHP-V2",
                          file=sys.stderr)
                thresh = adaptive_clade_threshold(len(subcluster))
                nwk = run_level2(
                    sub_fa,
                    clade_threshold=thresh,
                    verbose=verbose,
                )
                if nwk is None:
                    from fusang_v4_dahp_v1 import build_nj
                    sc_D = D[np.ix_(subcluster, subcluster)]
                    nwk = build_nj(sc_D, sc_names)
                sub_tree_entries.append((subcluster, nwk))
            else:
                # Build simple NJ tree
                from fusang_v4_dahp_v1 import build_nj
                sc_D = D[np.ix_(subcluster, subcluster)]
                nwk = build_nj(sc_D, sc_names)
                sub_tree_entries.append((subcluster, nwk))

        # ================================================================
        # Merge subtrees within this L0 cluster
        # ================================================================
        if len(sub_tree_entries) > 1:
            if verbose:
                print(f"[MHL]   Merging {len(sub_tree_entries)} subtrees",
                      file=sys.stderr)
            merged_nwk = merge_level(
                cluster, sub_tree_entries,
                D, taxon_names,
                verbose=verbose,
            )
            all_subtree_results.append((cluster, merged_nwk))
        else:
            all_subtree_results.append(sub_tree_entries[0])

    # ================================================================
    # Final merge: all L0 clusters into one tree
    # ================================================================
    if n_clusters > 1:
        if verbose:
            print(f"\n[MHL] Final merge: {len(all_subtree_results)} L0 clusters",
                  file=sys.stderr)

        final_nwk = merge_level(
            list(range(n)), all_subtree_results,
            D, taxon_names,
            n_bridge_taxa=MERGE_DEFAULTS["n_bridge_taxa"],
            bridge_strategy=MERGE_DEFAULTS["bridge_selection"],
            verbose=verbose,
        )
    else:
        final_nwk = all_subtree_results[0][1] if all_subtree_results else None

    # ================================================================
    # Validation + Output
    # ================================================================
    if final_nwk is not None:
        # Check leaf completeness
        import re
        leaves_re = set(re.findall(r"[A-Za-z0-9_\.]+(?=:)", final_nwk))
        leaves_in_tree = {l for l in leaves_re if not l.replace(".", "").isdigit()}
        expected = set(taxon_names)
        missing = expected - leaves_in_tree
        n_leaves = len(leaves_in_tree)
        if missing and verbose:
            print(f"[MHL] WARNING: {len(missing)} taxa missing from output tree",
                  file=sys.stderr)

        # Save
        write_newick(final_nwk, output_path)

        if verbose:
            print(f"\n[MHL] Output: {output_path}", file=sys.stderr)
            print(f"[MHL] Leaves in tree: {leaves_in_tree}/{n}", file=sys.stderr)
    else:
        if verbose:
            print("[MHL] ERROR: no tree produced", file=sys.stderr)

    # Collect stats
    stats = {
        "n_taxa": n,
        "n_clusters": n_clusters,
        "k": k,
        "gap": gap,
        "leaves_in_tree": n_leaves if final_nwk else 0,
        "output_path": output_path,
    }

    return {
        "tree": final_nwk,
        "stats": stats,
        "hierarchy": hierarchy,
    }


def main():
    args = parse_args()

    result = run_fusang_mhl(
        fasta_path=args.input,
        output_path=args.output,
        k=args.k,
        gap=args.gap,
        n_threads=args.threads,
        verbose=args.verbose,
        debug=args.debug,
        no_l2=args.no_l2,
        no_l3=args.no_l3,
        boundary_model=args.boundary_model,
        output_dir=args.output_dir,
    )

    if result["tree"] is None:
        sys.exit(1)

    if args.debug:
        import json
        stats = result["stats"]
        print(json.dumps(stats, indent=2, default=str), file=sys.stderr)


if __name__ == "__main__":
    main()
