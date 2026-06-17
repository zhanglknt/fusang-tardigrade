"""
Level 0: k-mer cosine distance clustering and backbone NJ tree construction.

This is the coarsest level, operating on all n taxa with O(n^2) k-mer distances.
Groups sequences into L0 clusters using hierarchical clustering, then builds
a backbone NJ tree from cluster centroids.
"""

import sys
import os
import numpy as np
from typing import Dict, List, Tuple, Set, Optional

# Ensure parent directory is on path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from .config import (
    L0_THRESHOLDS, L0_K, L0_GAP, L0_DISTANCE_METHOD, FUSANG_ROOT
)
from .mlh_utils import Timer


def get_l0_config(n_taxa: int) -> dict:
    """Get adaptive L0 configuration based on total taxa count.

    Returns:
        dict with 'max_group_size', 'target_groups', 'should_split'
    """
    for min_n, max_group, target_groups in reversed(L0_THRESHOLDS):
        if n_taxa >= min_n:
            return {
                "max_group_size": max_group,
                "target_groups": target_groups,
                "should_split": target_groups > 1,
            }
    return {"max_group_size": 200, "target_groups": 1, "should_split": False}


def resolve_gap_pattern(gap_str: str) -> Optional[tuple]:
    """Convert gap shorthand to tuple pattern."""
    if gap_str == "contiguous" or gap_str is None:
        return None
    presets = {
        "gap2": (0, 1, 2, 5, 6),
        "gap3": (0, 1, 2, 6, 7),
    }
    return presets.get(gap_str)


def compute_l0_distance(
    sequences: List[str],
    taxon_names: List[str],
    k: int = L0_K,
    gap: str = L0_GAP,
    metric: str = L0_DISTANCE_METHOD,
    n_threads: int = 4,
) -> np.ndarray:
    """Compute Level 0 k-mer distance matrix.

    Args:
        sequences: List of DNA sequence strings
        taxon_names: List of taxon names (same order)
        k: k-mer length
        gap: Gap pattern ("gap2", "gap3", "contiguous")
        metric: Distance metric ("cosine", "euclidean")
        n_threads: Number of threads for k-mer computation

    Returns:
        n x n distance matrix (numpy array)
    """
    gap_pattern = resolve_gap_pattern(gap)
    from kmer_distance import compute_kmer_distance_matrix
    D = compute_kmer_distance_matrix(
        sequences, taxon_names,
        k=k,
        metric=metric,
        n_threads=n_threads,
        gap_pattern=gap_pattern,
    )
    return np.array(D)


def cluster_taxa(
    D: np.ndarray,
    taxon_names: List[str],
    max_group_size: int = 200,
    target_groups: int = 3,
    method: str = "average",
) -> List[List[int]]:
    """Cluster taxa using hierarchical clustering.

    Args:
        D: n x n distance matrix
        taxon_names: List of taxon names
        max_group_size: Maximum taxa per cluster
        target_groups: Target number of clusters
        method: Linkage method ("average", "complete", "single")

    Returns:
        List of clusters, each cluster is a list of taxon indices
    """
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform

    # Convert to condensed form (upper triangle, row-major)
    n = len(taxon_names)
    condensed = squareform(D, checks=False)

    # Hierarchical clustering
    Z = linkage(condensed, method=method)

    # Determine cut threshold
    if target_groups <= 1:
        return [list(range(n))]

    # Use fcluster with target number of groups
    clusters_idx = fcluster(Z, t=target_groups, criterion="maxclust")

    # Group by cluster label
    cluster_map: Dict[int, List[int]] = {}
    for i, c in enumerate(clusters_idx):
        cluster_map.setdefault(c, []).append(i)

    clusters = list(cluster_map.values())

    # Post-processing: split oversized clusters
    final_clusters = []
    for cluster in clusters:
        if len(cluster) <= max_group_size:
            final_clusters.append(cluster)
        else:
            # Recursively split oversized cluster
            sub_D = D[np.ix_(cluster, cluster)]
            sub_names = [taxon_names[i] for i in cluster]
            sub_target = max(2, len(cluster) // max_group_size + 1)
            sub_clusters = cluster_taxa(
                sub_D, sub_names,
                max_group_size=max_group_size,
                target_groups=sub_target,
                method=method,
            )
            for sc in sub_clusters:
                final_clusters.append([cluster[i] for i in sc])

    return final_clusters


def build_backbone_nj(
    D: np.ndarray,
    taxon_names: List[str],
    cluster_centroids: List[int],
) -> Optional[str]:
    """Build backbone NJ tree from cluster centroids.

    Args:
        D: Full n x n distance matrix
        taxon_names: List of taxon names
        cluster_centroids: Index of centroid taxon for each cluster

    Returns:
        Newick string of backbone tree, or None on error
    """
    from fusang_v4_dahp_v1 import build_nj

    # Extract centroid-to-centroid distance matrix
    centroid_D = D[np.ix_(cluster_centroids, cluster_centroids)]
    centroid_names = [taxon_names[i] for i in cluster_centroids]

    nwk = build_nj(centroid_D, centroid_names)
    return nwk


def select_centroids(
    D: np.ndarray,
    cluster_indices: List[List[int]],
) -> List[int]:
    """Select centroid taxon for each cluster (minimize avg distance to others).

    Args:
        D: Full n x n distance matrix
        cluster_indices: List of clusters, each containing taxon indices

    Returns:
        List of centroid indices (one per cluster)
    """
    centroids = []
    for cluster in cluster_indices:
        if len(cluster) <= 1:
            centroids.append(cluster[0])
            continue
        sub_D = D[np.ix_(cluster, cluster)]
        # Centroid = taxon with minimum average distance to all others
        avg_dists = sub_D.mean(axis=1)
        centroid_local = np.argmin(avg_dists)
        centroids.append(cluster[centroid_local])
    return centroids


def compute_cluster_stats(
    D: np.ndarray,
    cluster_indices: List[List[int]],
    taxon_names: List[str],
) -> List[dict]:
    """Compute summary statistics for each cluster.

    Args:
        D: Full n x n distance matrix
        cluster_indices: List of clusters (taxon indices)
        taxon_names: List of all taxon names

    Returns:
        List of dicts with cluster statistics
    """
    stats = []
    for i, cluster in enumerate(cluster_indices):
        sub_D = D[np.ix_(cluster, cluster)]
        n_c = len(cluster)
        upper_idx = np.triu_indices(n_c, k=1)
        pairwise = sub_D[upper_idx]

        stats.append({
            "cluster_id": i,
            "size": n_c,
            "mean_dist": float(np.mean(pairwise)) if len(pairwise) > 0 else 0.0,
            "std_dist": float(np.std(pairwise)) if len(pairwise) > 0 else 0.0,
            "max_dist": float(np.max(pairwise)) if len(pairwise) > 0 else 0.0,
            "min_dist": float(np.min(pairwise)) if len(pairwise) > 0 else 0.0,
            "median_dist": float(np.median(pairwise)) if len(pairwise) > 0 else 0.0,
            "taxa": [taxon_names[j] for j in cluster],
        })
    return stats


def run_level0(
    sequences: List[str],
    taxon_names: List[str],
    k: int = L0_K,
    gap: str = L0_GAP,
    metric: str = L0_DISTANCE_METHOD,
    n_threads: int = 4,
    verbose: bool = True,
) -> Tuple[
    np.ndarray,                    # Full distance matrix D
    List[List[int]],               # Cluster assignments (indices)
    Optional[str],                  # Backbone NJ tree (Newick) or None
    List[int],                     # Centroid indices
    List[dict],                    # Cluster statistics
]:
    """Run Level 0: k-mer distance clustering.

    Args:
        sequences: List of DNA sequences
        taxon_names: List of taxon names
        k: k-mer size
        gap: Gap pattern
        metric: Distance metric
        n_threads: Threads for k-mer computation
        verbose: Print progress

    Returns:
        Tuple of (D, clusters, backbone_nwk, centroids, cluster_stats)
    """
    n = len(taxon_names)
    if verbose:
        print(f"[L0] Computing k-mer distance matrix (n={n}, k={k}, gap={gap})",
              file=sys.stderr)

    with Timer("L0: k-mer distance", verbose=verbose):
        D = compute_l0_distance(sequences, taxon_names, k=k, gap=gap,
                                metric=metric, n_threads=n_threads)

    # Get adaptive config
    l0_config = get_l0_config(n)
    if not l0_config["should_split"]:
        if verbose:
            print(f"[L0] n={n} <= 200, skipping L0 splitting", file=sys.stderr)
        clusters = [list(range(n))]
        centroids = [0]  # single centroid
        backbone_nwk = None
        stats = compute_cluster_stats(D, clusters, taxon_names)
        return D, clusters, backbone_nwk, centroids, stats

    if verbose:
        print(f"[L0] Clustering into ~{l0_config['target_groups']} groups "
              f"(max_size={l0_config['max_group_size']})", file=sys.stderr)

    with Timer("L0: clustering", verbose=verbose):
        clusters = cluster_taxa(
            D, taxon_names,
            max_group_size=l0_config["max_group_size"],
            target_groups=l0_config["target_groups"],
        )

    with Timer("L0: centroid selection + backbone", verbose=verbose):
        centroids = select_centroids(D, clusters)
        backbone_nwk = build_backbone_nj(D, taxon_names, centroids)

    stats = compute_cluster_stats(D, clusters, taxon_names)

    if verbose:
        print(f"[L0] {len(clusters)} clusters: "
              + ", ".join(f"c{i}={s['size']}(d={s['mean_dist']:.3f})"
                         for i, s in enumerate(stats)),
              file=sys.stderr)

    return D, clusters, backbone_nwk, centroids, stats
