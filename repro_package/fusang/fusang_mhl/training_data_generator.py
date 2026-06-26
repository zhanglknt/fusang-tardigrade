"""
Training data generator for the boundary classifier.

Generates training samples from simulated trees:
  - Label 1: cluster contains >1 natural subclade (should split)
  - Label 0: cluster is monophyletic and close (should stop)

V4 FIXES (2026-06-26):
  - Added generate_structured_tree() to create trees with 2-3 clear clades
  - Whole-tree labels are now GROUND TRUTH (1 for structured, 0 for coalescent)
  - Removed reliance on data-driven _label_whole_tree() (still kept as fallback)
  - Limit sub-cluster samples to 15 per tree (was unlimited → imbalanced)
  - Remove n_samples_target early break — generate ALL configs
  - generate_training_data() uses 50% structured + 50% coalescent trees
"""

import sys
import os
import json
import tempfile
import pickle
from typing import List, Tuple, Dict, Any, Optional

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from .config import L1_DEFAULTS
from .mlh_utils import Timer, ensure_dir, setup_logger

import getpass
USER = getpass.getuser()
LOG_DIR = os.path.join(tempfile.gettempdir(), f"fusang_mhl_{USER}")
ensure_dir(LOG_DIR)
logger = setup_logger("train_gen", log_file=os.path.join(LOG_DIR, "train_gen.log"))


# ---------------------------------------------------------------------------
# Simulation Parameters
# ---------------------------------------------------------------------------

DEFAULT_SIM_CONFIGS = [
    # (n, L, sub, indel, n_trees)
    # n=50
    (50,   500,  0.01, 0.000, 50),
    (50,   500,  0.05, 0.01,  50),
    (50,   500,  0.08, 0.03,  40),
    # n=100
    (100,  500,  0.02, 0.000, 45),
    (100,  500,  0.05, 0.01,  45),
    (100,  500,  0.08, 0.03,  35),
    (100,  2000, 0.03, 0.01,  25),
    # n=200
    (200,  500,  0.03, 0.000, 40),
    (200,  500,  0.05, 0.02,  40),
    (200,  500,  0.10, 0.08,  30),
    (200,  2000, 0.05, 0.02,  25),
    # n=300
    (300,  1000, 0.02, 0.000, 30),
    (300,  1000, 0.04, 0.02,  30),
    (300,  1000, 0.06, 0.04,  25),
    # n=500
    (500,  1000, 0.02, 0.000, 20),
    (500,  1000, 0.03, 0.02,  20),
    (500,  1000, 0.05, 0.03,  20),
    (500,  2000, 0.03, 0.01,  15),
]


# ---------------------------------------------------------------------------
# Core Generation Logic
# ---------------------------------------------------------------------------

_CODE_MAP = {0: 'A', 1: 'T', 2: 'C', 3: 'G'}


def _array_to_fasta_dict(seqs_array: np.ndarray, n: int) -> Dict[str, str]:
    """Convert (n, L) int8 array to {taxon_name: sequence_string} dict."""
    result = {}
    for i in range(n):
        seq_str = ''.join(_CODE_MAP[b] for b in seqs_array[i])
        result[f't{i+1:04d}'] = seq_str
    return result


def generate_simulated_tree(
    n: int,
    L: int,
    sub: float,
    indel: float,
    seed: int,
    verbose: bool = False,
) -> Tuple[str, List[str], List[str]]:
    """Generate one simulated coalescent tree + aligned sequences.

    Returns:
        newick_str: str (true tree in Newick)
        sequences: List[str]
        taxon_names: List[str]
    """
    from tree_simulation import make_coalescent_tree, simulate_seqs
    root_node, leaves = make_coalescent_tree(n, seed=seed)

    def _node_to_nwk(node):
        if node.left is None and node.right is None:
            return node.name or f't{node.idx}'
        left_str = _node_to_nwk(node.left)
        right_str = _node_to_nwk(node.right)
        bl = f':{node.branch_length:.6f}' if hasattr(node, 'branch_length') else ''
        return f'({left_str},{right_str}){node.name or ""}{bl}'
    true_nwk = _node_to_nwk(root_node) + ';'

    leaf_seqs = simulate_seqs(root_node, n, L, sub, seed, indel_rate=indel)
    seqs_dict = _array_to_fasta_dict(leaf_seqs, n)
    taxon_names = list(seqs_dict.keys())
    sequences = [seqs_dict[n] for n in taxon_names]
    return true_nwk, sequences, taxon_names


def generate_structured_tree(
    n: int,
    L: int,
    sub: float,
    indel: float,
    seed: int,
    n_clades: int = 2,
) -> Tuple[str, List[str], List[str]]:
    """Generate a tree with 2-3 clear clades.

    Each clade is evolved independently under a coalescent process
    with slightly different substitution parameters, creating detectable
    clusters in the k-mer distance matrix.

    The true Newick string encodes clade structure so that
    get_true_clades() can extract the ground-truth clades.

    Returns:
        newick_str: str (true tree with clade structure)
        sequences: List[str]
        taxon_names: List[str]
    """
    from tree_simulation import make_coalescent_tree, simulate_seqs

    rng = np.random.RandomState(seed)

    # Determine clade sizes (each clade has >= 4 taxa)
    base = max(4, n // n_clades)
    sizes = [base] * n_clades
    rem = n - base * n_clades
    for i in range(rem):
        sizes[rng.randint(0, n_clades - 1)] += 1

    clade_newick_parts = []
    all_seqs_arrays = []
    all_names = []
    offset = 0

    for ci, size in enumerate(sizes):
        clade_seed = seed + ci * 10000 + 7

        # Different substitution/indel rates per clade
        # This creates detectable clustering in distance space
        clade_sub = max(0.001, sub * rng.uniform(0.4, 2.5))
        clade_indel = max(0.0, indel * rng.uniform(0.4, 2.0))

        # Generate coalescent tree for this clade
        root, leaves = make_coalescent_tree(size, seed=clade_seed)

        # Rename leaves to be globally unique
        for i, leaf in enumerate(leaves):
            leaf.name = f't{offset + i + 1:04d}'

        # Get Newick for this clade's tree
        def _to_nwk(node, _ci=ci):
            if node.left is None and node.right is None:
                return node.name
            left = _to_nwk(node.left)
            right = _to_nwk(node.right)
            bl = f':{node.branch_length:.6f}' if hasattr(node, 'branch_length') else ''
            return f'({left},{right}){bl}'

        clade_nwk = _to_nwk(root)
        clade_newick_parts.append(clade_nwk)

        # Simulate sequences for this clade
        seqs_arr = simulate_seqs(
            root, size, L, clade_sub, clade_seed, indel_rate=clade_indel
        )
        all_seqs_arrays.append(seqs_arr)

        # Record names
        for i in range(size):
            all_names.append(f't{offset + i + 1:04d}')

        offset += size

    # Build combined Newick: (clade0_nwk,clade1_nwk,...);
    combined_nwk = '(' + ','.join(clade_newick_parts) + ');'

    # Combine sequence arrays
    combined_seqs = np.concatenate(all_seqs_arrays, axis=0)

    # Convert to sequence strings
    sequences = []
    for i in range(n):
        seq_str = ''.join(_CODE_MAP[b] for b in combined_seqs[i])
        sequences.append(seq_str)

    return combined_nwk, sequences, all_names


def get_true_clades(
    newick_str: str,
    taxon_names: List[str],
    threshold: float = 0.01,
) -> List[set]:
    """Extract 'true' clade groups from the true tree.

    Uses a simple approach: cut the true tree at `threshold` branch length
    to get natural clades. Returns list of sets of taxon names.
    """
    try:
        from fusang_v2 import TreeNode
        root = TreeNode.from_newick(newick_str)
        clades = []
        _extract_clades_at_threshold(root, threshold, clades)
        if not clades:
            return [set(taxon_names)]
        return clades
    except Exception:
        return [set(taxon_names)]


def _extract_clades_at_threshold(node, threshold, clades):
    """Recursively extract clades with root-ward branch <= threshold."""
    if not hasattr(node, 'children') or not node.children:
        return
    if hasattr(node, 'dist') and node.dist <= threshold:
        leaves = _get_leaf_names(node)
        if len(leaves) >= 4:
            clades.append(set(leaves))
    for ch in node.children:
        _extract_clades_at_threshold(ch, threshold, clades)


def _get_leaf_names(node):
    if not hasattr(node, 'children') or not node.children:
        return {node.name}
    names = set()
    for ch in node.children:
        names.update(_get_leaf_names(ch))
    return names


def label_cluster(
    cluster_taxon_names: List[str],
    true_clades: List[set],
) -> int:
    """Label a cluster based on true clades.

    Returns:
        1: cluster contains >1 true clade (should split)
        0: cluster is entirely within one true clade (stop)
    """
    cluster_set = set(cluster_taxon_names)
    n_clades_touched = sum(1 for c in true_clades if len(c & cluster_set) > 0)
    if n_clades_touched > 1 and len(cluster_set) >= 4:
        return 1  # Should split
    return 0  # Monophyletic, stop


def _label_whole_tree(D: np.ndarray, true_clades: List[set], n: int) -> int:
    """Label whether the ENTIRE dataset should be split at the top level (L0).

    V4 FALLBACK: used only when is_structured flag is not available.
    Uses hierarchical clustering on D to detect 2+ natural clusters.

    Returns:
        1: tree has 2+ substantial clusters → should split
        0: tree is star-like or has 1 dominant clade → stop
    """
    if n < 40:
        return 0

    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform

    condensed = squareform(D, checks=False)
    Z = linkage(condensed, method="average")

    for pct in [30, 40, 50, 60, 70]:
        t = Z[:, 2].max() * pct / 100.0
        if t <= 0.0:
            continue
        labels_arr = fcluster(Z, t=t, criterion="distance")
        counts = np.bincount(labels_arr)
        min_size = max(4, int(n * 0.10))
        substantial = sum(1 for c in counts if c >= min_size)
        if substantial >= 2:
            return 1

    # Fallback: check true_clades from tree topology
    if len(true_clades) >= 2:
        clade_sizes = sorted([len(c) for c in true_clades], reverse=True)
        max_pct = clade_sizes[0] / n
        if max_pct <= 0.75:
            substantial = sum(1 for s in clade_sizes if s >= max(4, n * 0.10))
            if substantial >= 2:
                return 1

    return 0


def _extract_features_for_indices(
    D: np.ndarray,
    sequences: List[str],
    taxon_names: List[str],
    indices: List[int],
    n_total: int,
) -> Optional[np.ndarray]:
    """Extract 50-dim features for a given set of indices."""
    from .level1_multik import extract_cluster_features
    try:
        features = extract_cluster_features(
            D=D, D_k5=None, D_k7=None, D_k9=None,
            seqs_for_cluster=[sequences[i] for i in indices],
            cluster_indices=list(range(len(indices))),
            centroid_idx=0,
            feature_matrix=None,
            n_total=n_total,
            parent_size=n_total,
            sibling_sizes=[],
            current_level=0,
            ancestor_sizes=[],
        )
        if features is not None and len(features) == 50:
            return features
    except Exception:
        pass
    return None


def generate_training_sample(
    newick_str: str,
    sequences: List[str],
    taxon_names: List[str],
    n_samples_per_tree: int = 20,
    verbose: bool = False,
    is_structured: bool = False,
) -> List[Dict]:
    """Generate training samples from one simulated tree.

    Produces TWO types of samples:
    1. Sub-cluster samples: features for clusters at various distance
       thresholds, labeled by whether they contain >1 true clade.
    2. Whole-tree sample: features for the ENTIRE dataset, labeled by
       whether the tree has structure warranting top-level split.

    Args:
        newick_str: True tree Newick
        sequences: List of sequences
        taxon_names: List of taxon names
        n_samples_per_tree: Max number of sub-cluster samples (V4: capped at 15)
        is_structured: Whether the tree has known clade structure.
                      If True, whole-tree label = 1 (ground truth).
                      If False, whole-tree label = 0 (coalescent = no structure).

    Returns:
        List of dicts: {features: List[float], label: int, meta: dict}
    """
    n = len(taxon_names)
    if n < 4:
        return []

    # Compute distance matrix
    from .level0_kmer import compute_l0_distance
    with Timer("  [train] k-mer distance", verbose=verbose):
        D = compute_l0_distance(sequences, taxon_names)

    # Get true clades
    true_clades = get_true_clades(newick_str, taxon_names)

    # Run hierarchical clustering at multiple cut points
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    condensed = squareform(D, checks=False)
    Z = linkage(condensed, method="average")

    samples = []
    rng = np.random.RandomState(hash(tuple(taxon_names)) % (2**31))

    # ── Type 1: Sub-cluster samples ──────────────────────────────────
    # Sample clusters at multiple distance thresholds
    z_max = Z[:, 2].max()
    cut_thresholds = []
    for pct in [10, 20, 30, 40, 50, 60, 70, 80, 90]:
        t = z_max * pct / 100.0
        if t > 0:
            cut_thresholds.append(t)

    seen_cluster_sets = set()
    all_cluster_samples = []

    for t in cut_thresholds:
        try:
            labels_arr = fcluster(Z, t=t, criterion="distance")
            n_clusters = labels_arr.max()
            for cid in range(1, n_clusters + 1):
                indices = [i for i, l in enumerate(labels_arr) if l == cid]
                if len(indices) < 4 or len(indices) > n * 0.95:
                    continue
                cluster_key = frozenset(indices)
                if cluster_key in seen_cluster_sets:
                    continue
                seen_cluster_sets.add(cluster_key)

                cluster_names = [taxon_names[i] for i in indices]
                features = _extract_features_for_indices(
                    D, sequences, taxon_names, indices, n,
                )
                if features is None:
                    continue

                label = label_cluster(cluster_names, true_clades)
                all_cluster_samples.append({
                    "features": features,
                    "label": label,
                    "meta": {
                        "sample_type": "sub_cluster",
                        "n_taxa": n,
                        "cluster_size": len(indices),
                        "dist_threshold": round(t, 6),
                        "n_total_clusters": n_clusters,
                    }
                })
        except Exception as e:
            if verbose:
                print(f"  [train] Cluster sampling error: {e}", file=sys.stderr)
            continue

    # V4 FIX: Limit sub-cluster samples to 15 per tree (was unlimited)
    # This makes whole-tree samples ~6% → ~30% of total
    if len(all_cluster_samples) > 15:
        idx = rng.choice(len(all_cluster_samples), size=15, replace=False)
        all_cluster_samples = [all_cluster_samples[i] for i in idx]

    samples.extend(all_cluster_samples)

    # ── Type 2: Whole-tree sample ──────────────────────────────────
    # V4 FIX: Use GROUND TRUTH instead of data-driven _label_whole_tree()
    #   is_structured=True  → label=1 (tree has clear clades → should SPLIT)
    #   is_structured=False → label=0 (coalescent/star-like → should STOP)
    all_indices = list(range(n))
    whole_features = _extract_features_for_indices(
        D, sequences, taxon_names, all_indices, n,
    )
    if whole_features is not None:
        # GROUND TRUTH labeling
        whole_label = 1 if is_structured else 0

        samples.append({
            "features": whole_features,
            "label": whole_label,
            "meta": {
                "sample_type": "whole_tree",
                "n_taxa": n,
                "cluster_size": n,
                "n_true_clades": len(true_clades),
                "is_structured": is_structured,
            }
        })

    return samples


# ---------------------------------------------------------------------------
# Batch Generation
# ---------------------------------------------------------------------------

def generate_training_data(
    configs: List[Tuple] = DEFAULT_SIM_CONFIGS,
    output_pkl: str = "training_data.pkl",
    n_samples_target: int = 20000,
    verbose: bool = True,
    structured_ratio: float = 0.5,
) -> List[Dict]:
    """Generate training data for boundary classifier.

    V4 FIXES:
      - Generates structured trees (with clear clades) for 50% of samples
      -Uses ground-truth whole-tree labels (not data-driven)
      -Removes early break: processes ALL configs
      -Subsamples to n_samples_target AFTER generating all data

    Args:
        configs: List of (n, L, sub, indel, n_trees) tuples
        output_pkl: Output pickle file path
        n_samples_target: Target number of samples (subsampled after generation)
        structured_ratio: Fraction of trees that are structured (default 0.5)
        verbose: Print progress

    Returns:
        List of sample dicts
    """
    all_samples = []
    total_trees = sum(c[4] for c in configs)

    if verbose:
        print(f"[train_gen] Generating data from {total_trees} trees "
              f"({structured_ratio*100:.0f}% structured)...", file=sys.stderr)

    for ci, (n, L, sub, indel, n_trees) in enumerate(configs):
        if verbose:
            print(f"\n[train_gen] Config {ci+1}/{len(configs)}: "
                  f"n={n}, L={L}, sub={sub}, indel={indel}, trees={n_trees}",
                  file=sys.stderr)

        for ti in range(n_trees):
            seed = ci * 1000 + ti

            # V4: 50% chance of structured tree
            rng_local = np.random.RandomState(seed)
            is_structured = rng_local.uniform() < structured_ratio

            try:
                if is_structured and n >= 8:
                    true_nwk, seqs, names = generate_structured_tree(
                        n, L, sub, indel, seed,
                        n_clades=rng_local.randint(2, 4),  # 2 or 3 clades
                    )
                else:
                    true_nwk, seqs, names = generate_simulated_tree(
                        n, L, sub, indel, seed, verbose=False,
                    )

                samples = generate_training_sample(
                    true_nwk, seqs, names,
                    n_samples_per_tree=20,
                    verbose=False,
                    is_structured=is_structured,
                )
                all_samples.extend(samples)

            except Exception as e:
                if verbose:
                    print(f"  [train_gen] Tree {ti+1} FAILED: {e}",
                          file=sys.stderr)
                continue

        if verbose:
            n_pos = sum(1 for s in all_samples if s["label"] == 1)
            n_neg = len(all_samples) - n_pos
            print(f"  [train_gen] Collected {len(all_samples)} samples "
                  f"({n_pos} pos, {n_neg} neg)",
                  file=sys.stderr)

    # Subsample to n_samples_target if we have too many
    if len(all_samples) > n_samples_target:
        rng_final = np.random.RandomState(42)
        idx = rng_final.choice(len(all_samples), size=n_samples_target, replace=False)
        all_samples = [all_samples[i] for i in idx]
        if verbose:
            print(f"\n[train_gen] Subsampled to {n_samples_target} samples",
                  file=sys.stderr)

    # Statistics
    n_pos = sum(1 for s in all_samples if s["label"] == 1)
    n_neg = len(all_samples) - n_pos
    if verbose:
        print(f"\n[train_gen] Collected {len(all_samples)} samples: "
              f"{n_pos} positive, {n_neg} negative", file=sys.stderr)

    # Per-type statistics
    from collections import Counter
    type_counts = Counter(s.get("meta", {}).get("sample_type", "?") for s in all_samples)
    if verbose:
        print(f"[train_gen] Sample types: {dict(type_counts)}", file=sys.stderr)

    # Save (optional)
    if output_pkl is not None:
        output_path = os.path.abspath(output_pkl)
        ensure_dir(os.path.dirname(output_path))
        with open(output_path, "wb") as f:
            pickle.dump(all_samples, f, protocol=pickle.HIGHEST_PROTOCOL)
        if verbose:
            print(f"[train_gen] Saved to {output_path}", file=sys.stderr)

    if verbose:
        print(f"[train_gen] Total samples: {len(all_samples)}", file=sys.stderr)

    return all_samples


def load_training_data(pkl_path: str) -> List[Dict]:
    """Load training data from pickle file."""
    with open(pkl_path, "rb") as f:
        data = pickle.load(f)
    return data


def split_training_data(
    data: List[Dict],
    test_size: float = 0.2,
    random_state: int = 42,
) -> Tuple[List[Dict], List[Dict]]:
    """Split training data into train/test sets."""
    rng = np.random.RandomState(random_state)
    indices = rng.permutation(len(data))
    n_test = int(len(data) * test_size)
    test_idx = set(indices[:n_test])
    train = [data[i] for i in range(len(data)) if i not in test_idx]
    test = [data[i] for i in test_idx]
    return train, test


if __name__ == "__main__":
    output = generate_training_data(verbose=True)
    print(f"\nDone. Output: {output}")
