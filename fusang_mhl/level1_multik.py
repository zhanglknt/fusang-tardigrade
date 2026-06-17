"""
Level 1: Multi-k ensemble + boundary classifier feature extraction.

For clusters output from Level 0, refine within-cluster resolution
using multi-k distance fusion and extract features for the boundary
classifier to decide whether to continue splitting.
"""

import sys
import os
import numpy as np
from typing import Dict, List, Tuple, Set, Optional, Any

from .config import L1_DEFAULTS
from .mlh_utils import Timer


def compute_l1_distance(
    sequences: List[str],
    taxon_names: List[str],
    ks: Tuple[int, ...] = (5, 7, 9),
    gap: Optional[str] = "contiguous",
    fusion_method: str = "average",
    n_threads: int = 4,
) -> np.ndarray:
    gap_pattern = None if gap == "contiguous" else gap
    from fusang_v4_dahp_v1 import compute_multik_distance_matrix
    D = compute_multik_distance_matrix(
        sequences, taxon_names,
        ks=ks,
        gap_pattern=gap_pattern,
        n_threads=n_threads,
        fusion_method=fusion_method,
    )
    return np.array(D, dtype=np.float64)


def build_l1_nj_tree(
    D: np.ndarray,
    taxon_names: List[str],
) -> Optional[str]:
    from fusang_v4_dahp_v1 import build_nj
    return build_nj(D, taxon_names)


def _pairwise_dists(D: np.ndarray) -> np.ndarray:
    n = D.shape[0]
    if n <= 1:
        return np.array([])
    iu = np.triu_indices(n, k=1)
    return D[iu]


def extract_group_a(
    cluster_indices: List[int],
    parent_size: int,
    sibling_sizes: List[int],
) -> List[float]:
    n_c = len(cluster_indices)
    log_size = np.log10(n_c) if n_c > 0 else 0.0
    size_ratio = n_c / parent_size if parent_size > 0 else 1.0
    sizes = sibling_sizes + [n_c]
    sorted_sizes = sorted(sizes)
    if len(sorted_sizes) > 1:
        rank = sorted_sizes.index(n_c) / (len(sorted_sizes) - 1)
    else:
        rank = 0.5
    return [float(log_size), float(size_ratio), float(rank)]


def extract_group_b(
    D: np.ndarray,
    cluster_indices: List[int],
) -> List[float]:
    n_c = len(cluster_indices)
    if n_c <= 1:
        return [0.0] * 10
    sub_D = D[np.ix_(cluster_indices, cluster_indices)]
    pd = _pairwise_dists(sub_D)
    if len(pd) == 0:
        return [0.0] * 10
    from scipy import stats as scipy_stats
    mean_d = float(np.mean(pd))
    std_d = float(np.std(pd))
    min_d = float(np.min(pd))
    max_d = float(np.max(pd))
    median_d = float(np.median(pd))
    skew_d = float(scipy_stats.skew(pd))
    kurt_d = float(scipy_stats.kurtosis(pd))
    q25 = float(np.percentile(pd, 25))
    q75 = float(np.percentile(pd, 75))
    iqr_d = q75 - q25
    return [mean_d, std_d, min_d, max_d, median_d,
            skew_d, kurt_d, q25, q75, iqr_d]


def extract_group_c(
    D_k5: np.ndarray,
    D_k7: np.ndarray,
    D_k9: np.ndarray,
    cluster_indices: List[int],
) -> List[float]:
    n_c = len(cluster_indices)
    if n_c <= 1:
        return [0.0] * 6
    feats = []
    for D_src in [D_k5, D_k7, D_k9]:
        sub = D_src[np.ix_(cluster_indices, cluster_indices)]
        pd = _pairwise_dists(sub)
        feats.append(float(np.mean(pd)) if len(pd) > 0 else 0.0)
    if n_c > 1:
        diff_val = float(np.max(np.abs(D_k5 - D_k9)))
        eps = 1e-9
        pd5 = _pairwise_dists(D_k5[np.ix_(cluster_indices, cluster_indices)])
        pd9 = _pairwise_dists(D_k9[np.ix_(cluster_indices, cluster_indices)])
        m5 = float(np.mean(pd5)) if len(pd5) > 0 else eps
        m9 = float(np.mean(pd9)) if len(pd9) > 0 else eps
        cv_k5 = float(np.std(pd5)) / max(m5, eps)
        cv_k9 = float(np.std(pd9)) / max(m9, eps)
    else:
        diff_val = 0.0
        cv_k5 = 0.0
        cv_k9 = 0.0
    feats.extend([diff_val, cv_k5, cv_k9])
    return feats


def extract_group_d(
    feature_matrix: np.ndarray,
    centroid_idx: int,
    cluster_indices: List[int],
) -> List[float]:
    if len(cluster_indices) <= 1:
        return [0.0] * 8
    eps = 1e-12
    entropies = -np.sum(
        feature_matrix[cluster_indices] * np.log(feature_matrix[cluster_indices] + eps),
        axis=1,
    )
    mean_ent = float(np.mean(entropies))
    std_ent = float(np.std(entropies))
    sparsity = np.mean(feature_matrix[cluster_indices] < 1e-6, axis=1)
    mean_sparse = float(np.mean(sparsity))
    std_sparse = float(np.std(sparsity))
    centroid_vec = feature_matrix[centroid_idx]
    member_vecs = feature_matrix[cluster_indices]
    from sklearn.metrics.pairwise import cosine_distances
    c_dists = cosine_distances(
        centroid_vec.reshape(1, -1), member_vecs
    )[0]
    mean_cd = float(np.mean(c_dists))
    if len(cluster_indices) > 1:
        m_dists = cosine_distances(member_vecs)
        triu = m_dists[np.triu_indices(len(m_dists), k=1)]
        mean_md = float(np.mean(triu)) if len(triu) > 0 else 0.0
    else:
        mean_md = 0.0
    from sklearn.metrics.pairwise import euclidean_distances
    e_dists = euclidean_distances(
        centroid_vec.reshape(1, -1), member_vecs
    )[0]
    mean_ed = float(np.mean(e_dists))
    total_var = float(np.sum(np.var(feature_matrix[cluster_indices], axis=0)))
    return [mean_ent, std_ent, mean_sparse, std_sparse,
            mean_cd, mean_md, mean_ed, total_var]


def extract_group_e(sequences: List[str]) -> List[float]:
    gc_contents = []
    lengths = []
    for seq in sequences:
        if not seq:
            continue
        gc = (seq.count("G") + seq.count("C")) / len(seq)
        gc_contents.append(gc)
        lengths.append(len(seq))
    if not gc_contents:
        return [0.0] * 4
    return [float(np.mean(gc_contents)), float(np.std(gc_contents)),
            float(np.mean(lengths)), float(np.std(lengths))]


def _colless_index(node) -> float:
    children = node.clades if hasattr(node, "clades") else []
    if len(children) <= 1:
        return 0.0
    sizes = []
    for ch in children:
        tips = list(ch.get_terminals()) if hasattr(ch, "get_terminals") else []
        sizes.append(len(tips))
    if len(sizes) < 2:
        return 0.0
    val = abs(sizes[0] - sizes[1])
    for ch in children:
        val += _colless_index(ch)
    return float(val)


def extract_group_f(
    D: np.ndarray,
    cluster_indices: List[int],
) -> List[float]:
    if len(cluster_indices) <= 2:
        return [0.0] * 6
    try:
        from fusang_v4_dahp_v1 import build_nj
        from io import StringIO
        from Bio import Phylo
        sub_D = D[np.ix_(cluster_indices, cluster_indices)]
        sub_names = ["t" + str(i) for i in range(len(cluster_indices))]
        nwk = build_nj(sub_D, sub_names)
        if nwk is None:
            return [0.0] * 6
        tree = Phylo.read(StringIO(nwk), "newick")
        heights = [tree.distance(tip) for tip in tree.get_terminals()]
        tree_height = float(max(heights)) if heights else 0.0
        bls = [c.branch_length or 0.0 for c in tree.find_clades()]
        mean_bl = float(np.mean(bls)) if bls else 0.0
        std_bl = float(np.std(bls)) if bls else 0.0
        max_internal = 0.0
        n_internal = 0
        for c in tree.get_nonterminals():
            if not c.is_terminal():
                bl = c.branch_length or 0.0
                if bl > max_internal:
                    max_internal = bl
                n_internal += 1
        colless = _colless_index(tree.root) if tree.root else 0.0
        return [tree_height, mean_bl, std_bl, max_internal, float(n_internal), colless]
    except Exception:
        return [0.0] * 6


def extract_group_g(
    n_total: int,
    current_level: int,
    ancestor_sizes: List[int],
) -> List[float]:
    log_n = np.log10(n_total) if n_total > 0 else 0.0
    mean_anc = float(np.mean(ancestor_sizes)) if ancestor_sizes else 0.0
    return [float(log_n), float(current_level), mean_anc]


def extract_cluster_features(
    D: np.ndarray,
    D_k5: Optional[np.ndarray],
    D_k7: Optional[np.ndarray],
    D_k9: Optional[np.ndarray],
    seqs_for_cluster: List[str],
    cluster_indices: List[int],
    centroid_idx: int,
    feature_matrix: Optional[np.ndarray],
    n_total: int,
    parent_size: int,
    sibling_sizes: List[int],
    current_level: int,
    ancestor_sizes: List[int],
    verbose: bool = False,
) -> np.ndarray:
    """Extract full 50-dimension feature vector for a cluster."""
    features = []

    # Group A: Size (3)
    features.extend(
        extract_group_a(cluster_indices, parent_size, sibling_sizes)
    )

    # Group B: Distance distribution (10)
    features.extend(extract_group_b(D, cluster_indices))

    # Group C: Multi-k consistency (6)
    if D_k5 is not None and D_k7 is not None and D_k9 is not None:
        features.extend(
            extract_group_c(D_k5, D_k7, D_k9, cluster_indices)
        )
    else:
        features.extend([0.0] * 6)

    # Group D: K-mer frequency (8)
    if feature_matrix is not None:
        features.extend(
            extract_group_d(feature_matrix, centroid_idx, cluster_indices)
        )
    else:
        features.extend([0.0] * 8)

    # Group E: Sequence conservation (4)
    sub_seqs = [seqs_for_cluster[i] for i in range(len(cluster_indices))]
    features.extend(extract_group_e(sub_seqs))

    # Group F: NJ topology (6)
    features.extend(extract_group_f(D, cluster_indices))

    # Group G: Global context (3)
    features.extend(
        extract_group_g(n_total, current_level, ancestor_sizes)
    )

    while len(features) < 50:
        features.append(0.0)

    return np.array(features[:50], dtype=np.float32)


def run_level1(
    sequences: List[str],
    taxon_names: List[str],
    cluster_indices: List[List[int]],
    D_full: np.ndarray,
    ks: Tuple[int, ...] = (5, 7, 9),
    n_threads: int = 4,
    verbose: bool = True,
) -> Tuple[List[np.ndarray], List[str], List[np.ndarray]]:
    n_clusters = len(cluster_indices)
    if verbose:
        print(
            f"[L1] Processing {n_clusters} clusters with multi-k ensemble (k={ks})",
            file=sys.stderr,
        )

    cluster_Ds = []
    cluster_nwks = []
    cluster_features = []

    for ci, indices in enumerate(cluster_indices):
        if verbose:
            print(
                f"  [L1] Cluster {ci}/{n_clusters} (n={len(indices)})",
                file=sys.stderr,
            )

        sub_seqs = [sequences[i] for i in indices]
        sub_names = [taxon_names[i] for i in indices]

        with Timer(f"L1: multi-k cluster {ci}", verbose=verbose):
            D_sub = compute_l1_distance(
                sub_seqs, sub_names,
                ks=ks, gap="contiguous",
                n_threads=max(1, n_threads // max(n_clusters, 1) + 1),
            )

        with Timer(f"L1: NJ cluster {ci}", verbose=verbose):
            nwk = build_l1_nj_tree(D_sub, sub_names)

        feat = extract_cluster_features(
            D_sub, None, None, None,
            sub_seqs, list(range(len(indices))), 0,
            None,
            n_total=len(sequences),
            parent_size=max(len(c) for c in cluster_indices),
            sibling_sizes=[len(c) for c in cluster_indices if c != indices],
            current_level=1,
            ancestor_sizes=[len(c) for c in cluster_indices],
            verbose=False,
        )

        cluster_Ds.append(D_sub)
        cluster_nwks.append(nwk)
        cluster_features.append(feat)

        if verbose:
            pd = _pairwise_dists(D_sub)
            md = float(np.mean(pd)) if len(pd) > 0 else 0.0
            print(
                f"    mean_d={md:.4f}, "
                f"nwk={'OK' if nwk else 'FAIL'}",
                file=sys.stderr,
            )

    return cluster_Ds, cluster_nwks, cluster_features
