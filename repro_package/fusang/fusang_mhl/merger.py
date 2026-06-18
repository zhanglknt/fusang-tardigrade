"""
Bridge Taxa constrained supertree merger for Fusang MHL.

Implements bottom-up merging of multi-level phylogenetic trees using
representative taxa + NJ backbone, then subtree expansion.

Strategy:
  1. Each child subtree selects ONE representative taxon (centroid)
  2. NJ backbone is built on all representative taxa
  3. Each representative leaf is replaced by its full subtree
  4. Non-bridge taxa from parent that aren't in any child remain as leaves
"""

import sys
import os
from typing import Dict, List, Tuple, Set, Optional, Any

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from .config import MERGE_DEFAULTS
from .mlh_utils import Timer, ensure_dir


# ---------------------------------------------------------------------------
# Representative Taxon Selection
# ---------------------------------------------------------------------------

def select_representative_taxa(
    D: np.ndarray,
    cluster_indices: List[int],
    n_reps: int = 1,
) -> List[int]:
    """Select representative taxa from a cluster (centroid + diverse).

    Args:
        D: Full n x n distance matrix
        cluster_indices: Indices of taxa in this cluster
        n_reps: Number of representatives to select

    Returns:
        List of global indices of selected representative taxa
    """
    if len(cluster_indices) <= 0:
        return []

    n_reps = min(n_reps, len(cluster_indices))

    # Centroid = taxon with minimum average distance to others
    sub_D = D[np.ix_(cluster_indices, cluster_indices)]
    avg_dists = sub_D.mean(axis=1)
    sorted_local = list(np.argsort(avg_dists))

    return [cluster_indices[i] for i in sorted_local[:n_reps]]


def select_bridge_taxa(
    sequences: List[str],
    taxon_names: List[str],
    D: np.ndarray,
    cluster_indices: List[int],
    n_bridge: int = 3,
    strategy: str = "centroid+random",
    centroid_idx: Optional[int] = None,
    random_state: int = 42,
) -> List[int]:
    """Select bridge taxa from a cluster.

    Args:
        sequences: All sequences (full dataset)
        taxon_names: All taxon names (full dataset)
        D: Full n x n distance matrix
        cluster_indices: Indices of taxa in this cluster
        n_bridge: Number of bridge taxa to select
        strategy: "centroid+random" | "centroid+diverse" | "random"
        centroid_idx: Pre-computed centroid index (global index)
        random_state: RNG seed

    Returns:
        List of global indices of selected bridge taxa
    """
    if len(cluster_indices) <= 1:
        return list(cluster_indices)

    n_bridge = min(n_bridge, len(cluster_indices))

    if strategy == "centroid+random":
        return _select_bridge_centroid_random(
            D, cluster_indices, n_bridge, centroid_idx, random_state,
        )
    elif strategy == "centroid+diverse":
        return _select_bridge_centroid_diverse(
            D, cluster_indices, n_bridge, centroid_idx,
        )
    elif strategy == "random":
        rng = np.random.RandomState(random_state)
        return list(rng.choice(cluster_indices, size=n_bridge, replace=False))
    else:
        raise ValueError(f"Unknown bridge strategy: {strategy}")


def _select_bridge_centroid_random(
    D: np.ndarray,
    cluster_indices: List[int],
    n_bridge: int,
    centroid_idx: Optional[int],
    random_state: int,
) -> List[int]:
    """Select centroid + random taxa as bridge."""
    bridge = []
    if centroid_idx is not None and centroid_idx in cluster_indices:
        bridge.append(centroid_idx)
    else:
        sub_D = D[np.ix_(cluster_indices, cluster_indices)]
        avg_dists = sub_D.mean(axis=1)
        centroid_local = int(np.argmin(avg_dists))
        bridge.append(cluster_indices[centroid_local])

    remaining = [i for i in cluster_indices if i not in bridge]
    if remaining:
        rng = np.random.RandomState(random_state)
        n_extra = min(n_bridge - len(bridge), len(remaining))
        extra = list(rng.choice(remaining, size=n_extra, replace=False))
        bridge.extend(extra)

    return bridge


def _select_bridge_centroid_diverse(
    D: np.ndarray,
    cluster_indices: List[int],
    n_bridge: int,
    centroid_idx: Optional[int],
) -> List[int]:
    """Select centroid + taxa farthest from centroid (diverse)."""
    bridge = []

    if centroid_idx is not None and centroid_idx in cluster_indices:
        cidx = centroid_idx
        bridge.append(cidx)
    else:
        sub_D = D[np.ix_(cluster_indices, cluster_indices)]
        avg_dists = sub_D.mean(axis=1)
        centroid_local = int(np.argmin(avg_dists))
        cidx = cluster_indices[centroid_local]
        bridge.append(cidx)

    centroid_global_idx = bridge[0]
    dists_to_centroid = D[centroid_global_idx, cluster_indices]
    sorted_indices = sorted(
        [(d, idx) for d, idx in zip(dists_to_centroid, cluster_indices) if idx != centroid_global_idx],
        key=lambda x: -x[0],
    )
    for _, idx in sorted_indices[:n_bridge - 1]:
        bridge.append(idx)

    return bridge


# ---------------------------------------------------------------------------
# Subtree Expansion
# ---------------------------------------------------------------------------

def _expand_representative_taxa(
    node: Any,
    rep_map: Dict[str, str],
) -> None:
    """Recursively replace representative leaf nodes with their full subtrees.

    Args:
        node: TreeNode (from fusang_v2)
        rep_map: {representative_taxon_name: full_subtree_newick_string}
    """
    from fusang_v2 import TreeNode

    if not hasattr(node, 'children') or not node.children:
        # Leaf node: check if it's a representative to expand
        name = node.name
        if name in rep_map:
            subtree_nwk = rep_map.pop(name)  # Remove to prevent re-expansion
            try:
                subtree_root = TreeNode.from_newick(subtree_nwk)
                # Replace this leaf with the subtree root
                node.name = subtree_root.name
                node.children = subtree_root.children
                node.dist = subtree_root.dist
                # Fix parent references
                for ch in node.children:
                    ch.up = node
                # Do NOT recursively expand — the subtree is complete
                # (recursive expansion would cause infinite loop if the
                # subtree Newick contains the same representative name)
            except Exception as e:
                print(f"[merger] Warning: could not expand {name}: {e}",
                      file=sys.stderr)
        return

    # Internal node: recurse
    for ch in node.children:
        _expand_representative_taxa(ch, rep_map)


# ---------------------------------------------------------------------------
# Merge Level (Main Entry)
# ---------------------------------------------------------------------------

def merge_level(
    parent_taxa: List[int],
    children_data: List[Tuple[List[int], str]],
    D: np.ndarray,
    taxon_names: List[str],
    n_bridge_taxa: int = MERGE_DEFAULTS["n_bridge_taxa"],
    bridge_strategy: str = MERGE_DEFAULTS["bridge_selection"],
    verbose: bool = True,
) -> Optional[str]:
    """Merge children subtrees into parent level.

    Strategy: Each child subtree contributes ONE representative (centroid)
    to the NJ backbone. After NJ is built, each representative leaf
    is expanded into its full subtree.

    Args:
        parent_taxa: Global indices of all taxa in parent cluster
        children_data: List of (child_cluster_indices, child_newick_string)
        D: Full n x n distance matrix
        taxon_names: All taxon names
        n_bridge_taxa: Number of bridge taxa (unused in new approach, kept for API compat)
        bridge_strategy: Bridge selection strategy (unused in new approach)
        verbose: Print progress

    Returns:
        Newick string of merged tree, or None on failure
    """
    if verbose:
        print(f"[merge] Merging {len(children_data)} subtrees "
              f"({len(parent_taxa)} total taxa)", file=sys.stderr)

    if len(children_data) <= 1:
        if children_data:
            return children_data[0][1]
        return None

    # Collect all taxa that belong to children
    child_all = set()
    for child_indices, child_nwk in children_data:
        child_all.update(child_indices)

    # Select ONE representative per child (centroid)
    rep_map = {}  # representative_name -> child_subtree_newick
    rep_indices = []  # global indices of representatives
    all_child_taxa = set()

    for child_indices, child_nwk in children_data:
        all_child_taxa.update(child_indices)

        if len(child_indices) == 0:
            continue

        # Select centroid as representative
        sub_D = D[np.ix_(child_indices, child_indices)]
        avg_dists = sub_D.mean(axis=1)
        centroid_local = int(np.argmin(avg_dists))
        centroid_global = child_indices[centroid_local]
        centroid_name = taxon_names[centroid_global]

        rep_map[centroid_name] = child_nwk
        rep_indices.append(centroid_global)

    # Add orphan taxa (in parent but not in any child) as additional leaves
    orphan_indices = [pt for pt in parent_taxa if pt not in all_child_taxa]

    # All taxa in the merge NJ backbone: reps + orphans
    merge_taxa = list(set(rep_indices)) + orphan_indices

    if len(merge_taxa) < 2:
        if verbose:
            print("[merge] <2 taxa for backbone, returning first child",
                  file=sys.stderr)
        return children_data[0][1] if children_data else None

    # Build distance matrix for backbone
    m = len(merge_taxa)
    merge_D = np.zeros((m, m))
    merge_names = [taxon_names[i] for i in merge_taxa]
    for i in range(m):
        for j in range(i + 1, m):
            d = D[merge_taxa[i], merge_taxa[j]]
            merge_D[i, j] = d
            merge_D[j, i] = d

    # Build NJ backbone on representative taxa
    from fusang_v4_dahp_v1 import build_nj
    from fusang_v2 import TreeNode

    with Timer("merge: NJ backbone", verbose=verbose):
        nwk_merge = build_nj(merge_D, merge_names)

    if nwk_merge is None:
        if verbose:
            print("[merge] NJ failed, returning first child", file=sys.stderr)
        return children_data[0][1] if children_data else None

    # Expand representative taxa leaves into full subtrees
    try:
        root = TreeNode.from_newick(nwk_merge)
        _expand_representative_taxa(root, rep_map)
        result_nwk = root.to_newick()

        # Count leaves
        leaves = []
        def collect_leaves(n):
            if not n.children:
                leaves.append(n.name)
            else:
                for c in n.children:
                    collect_leaves(c)
        collect_leaves(root)

        if verbose:
            print(f"[merge] Merge OK: {len(leaves)} leaves "
                  f"(expected: {len(parent_taxa)})", file=sys.stderr)
        return result_nwk
    except Exception as e:
        if verbose:
            print(f"[merge] Expansion failed: {e}", file=sys.stderr)
            import traceback
            traceback.print_exc(file=sys.stderr)
        return nwk_merge


def constrained_nj_merge(
    bridge_D: np.ndarray,
    bridge_names: List[str],
    subtree_map: Dict[str, Any],
    constraints: Optional[List[List[str]]] = None,
) -> Optional['TreeNode']:
    """Build NJ tree with monophyletic constraints for bridge taxa.

    Args:
        bridge_D: Distance matrix for bridge taxa
        bridge_names: Names of bridge taxa
        subtree_map: {bridge_taxon_name: full_subtree_newick_string}
        constraints: Optional monophyly constraints

    Returns:
        TreeNode root of merged tree, or None on failure
    """
    from fusang_v4_dahp_v1 import build_nj
    from fusang_v2 import TreeNode

    nwk_bridge = build_nj(bridge_D, bridge_names)
    if nwk_bridge is None:
        return None

    bridge_root = TreeNode.from_newick(nwk_bridge)
    _expand_representative_taxa(bridge_root, subtree_map)
    return bridge_root


def build_bridge_distance_matrix(
    D: np.ndarray,
    subtree_bridge_map: List[Tuple[str, List[int], Any]],
    all_bridge_indices: List[int],
) -> np.ndarray:
    """Build distance matrix for bridge taxa NJ merge."""
    m = len(all_bridge_indices)
    bm = np.zeros((m, m))
    for i, gi in enumerate(all_bridge_indices):
        for j, gj in enumerate(all_bridge_indices):
            if i != j:
                bm[i, j] = D[gi, gj]
                bm[j, i] = bm[i, j]
    return bm


# ---------------------------------------------------------------------------
# Bottom-up Merge (Main Entry)
# ---------------------------------------------------------------------------

def bottom_up_merge(
    hierarchy: List[List[Tuple[List[int], Optional[str]]]],
    D: np.ndarray,
    taxon_names: List[str],
    n_bridge_taxa: int = MERGE_DEFAULTS["n_bridge_taxa"],
    bridge_strategy: str = MERGE_DEFAULTS["bridge_selection"],
    verbose: bool = True,
) -> Optional[str]:
    """Bottom-up merge of all levels.

    Args:
        hierarchy: List of levels, each level is a list of
                   (cluster_indices, newick_string_or_None)
                   Level 0 = top (root), last level = bottom (leaves)
        D: Full distance matrix
        taxon_names: All taxon names
        n_bridge_taxa: Bridge taxa count
        bridge_strategy: Bridge selection strategy
        verbose: Print progress

    Returns:
        Newick string of final merged tree
    """
    n_levels = len(hierarchy)
    if n_levels == 0:
        return None

    current_trees = hierarchy[-1]

    for level_idx in range(n_levels - 2, -1, -1):
        parent_level = hierarchy[level_idx]
        children_level = current_trees

        if verbose:
            print(f"[merge] Level {level_idx} -> {level_idx+1}: "
                  f"{len(children_level)} subtrees", file=sys.stderr)

        new_trees = []
        for parent_indices, parent_nwk in parent_level:
            matching_children = []
            for child_indices, child_nwk in children_level:
                if child_nwk is None:
                    continue
                child_set = set(child_indices)
                parent_set = set(parent_indices)
                if child_set.issubset(parent_set) and len(child_set) < len(parent_set):
                    matching_children.append((child_indices, child_nwk))

            if not matching_children:
                if parent_nwk is not None:
                    new_trees.append((parent_indices, parent_nwk))
                else:
                    if children_level:
                        new_trees.append(children_level[0])
                continue

            if len(matching_children) == 1:
                new_trees.append(matching_children[0])
                continue

            merged_nwk = merge_level(
                parent_indices, matching_children,
                D, taxon_names,
                verbose=verbose,
            )
            if merged_nwk is not None:
                new_trees.append((parent_indices, merged_nwk))
            else:
                new_trees.append(matching_children[0])

        current_trees = new_trees

    if current_trees:
        return current_trees[0][1]
    return None
