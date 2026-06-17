"""
Training data generator for the boundary classifier.

Generates ~50,000 labeled samples from simulated trees:
  - Label 1: cluster contains >1 natural subclade (should split)
  - Label 0: cluster is monophyletic and close (should stop)

Uses tree_simulation.py to generate coalescent trees + simulated sequences,
then runs hierarchical clustering and computes 50-dim feature vectors.
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
    (50,   500, 0.01, 0.001, 20),
    (100,  500, 0.02, 0.005, 20),
    (200,  500, 0.05, 0.02,  30),
    (300,  1000, 0.03, 0.01, 15),
    (500,  1000, 0.02, 0.01, 10),
    (50,   500, 0.08, 0.05,  15),  # high substitution
    (100,  500, 0.08, 0.05,  15),
    (200,  500, 0.08, 0.05,  20),
    (50,   1000, 0.01, 0.001, 15),  # long sequences
    (200,  1000, 0.05, 0.02,  20),
]


# ---------------------------------------------------------------------------
# Core Generation Logic
# ---------------------------------------------------------------------------

_CODE_MAP = {0: 'A', 1: 'T', 2: 'C', 3: 'G'}


def _array_to_fasta_dict(seqs_array: np.ndarray, n: int) -> Dict[str, str]:
    """Convert (n, L) int8 array to {taxon_name: sequence_string} dict.

    Names match tree_simulation.py convention: t0001, t0002, ...
    """
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
):
    """Generate one simulated tree + aligned sequences.

    Note: `indel` param is accepted for API compatibility but
    tree_simulation.py currently only simulates substitutions (JC69).
    For indel simulation, consider post-hoc gap insertion (TODO).

    Returns:
        newick_str: str (true tree in Newick)
        sequences: List[str]
        taxon_names: List[str]
    """
    from tree_simulation import make_coalescent_tree, simulate_seqs
    # True tree — returns (root_node, leaves) tuple
    root_node, leaves = make_coalescent_tree(n, seed=seed)

    # Convert root to Newick string (true tree)
    def _node_to_nwk(node):
        if node.left is None and node.right is None:
            return node.name or f't{node.idx}'
        left_str = _node_to_nwk(node.left)
        right_str = _node_to_nwk(node.right)
        bl = f':{node.branch_length:.6f}' if hasattr(node, 'branch_length') else ''
        return f'({left_str},{right_str}){node.name or ""}{bl}'
    true_nwk = _node_to_nwk(root_node) + ';'

    # Simulate sequences — returns (n, L) int8 array
    leaf_seqs = simulate_seqs(root_node, n, L, sub, seed)

    # Convert to {name: str} dict
    seqs_dict = _array_to_fasta_dict(leaf_seqs, n)
    taxon_names = list(seqs_dict.keys())
    sequences = [seqs_dict[n] for n in taxon_names]
    return true_nwk, sequences, taxon_names


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
        # Collect all leaves under each internal node at threshold
        clades = []
        _extract_clades_at_threshold(root, threshold, clades)
        if not clades:
            # Entire tree is one clade
            return [set(taxon_names)]
        return clades
    except Exception:
        # Fallback: entire tree is one clade
        return [set(taxon_names)]


def _extract_clades_at_threshold(node, threshold, clades):
    """Recursively extract clades with root-ward branch <= threshold."""
    if not hasattr(node, 'children') or not node.children:
        return
    # If this node's branch length is small, its subtree is a clade
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


def generate_training_sample(
    newick_str: str,
    sequences: List[str],
    taxon_names: List[str],
    n_samples_per_tree: int = 20,
    verbose: bool = False,
) -> List[Dict]:
    """Generate training samples from one simulated tree.

    Args:
        newick_str: True tree Newick
        sequences: List of sequences
        taxon_names: List of taxon names
        n_samples_per_tree: Number of clusters to sample

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

    # Sample clusters at different granularity levels
    for n_clusters in [2, 3, 4, 5, min(8, n // 10)]:
        if n_clusters < 2:
            continue
        try:
            labels = fcluster(Z, t=n_clusters, criterion="maxclust")
            for cid in range(1, n_clusters + 1):
                indices = [i for i, l in enumerate(labels) if l == cid]
                if len(indices) < 3 or len(indices) > n * 0.8:
                    continue
                cluster_names = [taxon_names[i] for i in indices]

                # Extract features
                from .level1_multik import extract_cluster_features
                features = extract_cluster_features(
                    D=D, D_k5=None, D_k7=None, D_k9=None,
                    seqs_for_cluster=[sequences[i] for i in indices],
                    cluster_indices=list(range(len(indices))),
                    centroid_idx=0,
                    feature_matrix=None,
                    n_total=n,
                    parent_size=n,
                    sibling_sizes=[],
                    current_level=0,
                    ancestor_sizes=[],
                )
                if features is None or len(features) != 50:
                    continue

                label = label_cluster(cluster_names, true_clades)
                samples.append({
                    "features": features,
                    "label": label,
                    "meta": {
                        "n_taxa": n,
                        "cluster_size": len(indices),
                        "n_clusters": n_clusters,
                    }
                })
        except Exception as e:
            if verbose:
                print(f"  [train] Cluster sampling error: {e}", file=sys.stderr)
            continue

    return samples


# ---------------------------------------------------------------------------
# Batch Generation
# ---------------------------------------------------------------------------

def generate_training_data(
    configs: List[Tuple] = DEFAULT_SIM_CONFIGS,
    output_pkl: str = "training_data.pkl",
    n_samples_target: int = 50000,
    verbose: bool = True,
) -> str:
    """Generate training data for boundary classifier.

    Args:
        configs: List of (n, L, sub, indel, n_trees) tuples
        output_pkl: Output pickle file path
        n_samples_target: Target number of samples
        verbose: Print progress

    Returns:
        Path to output pickle file
    """
    all_samples = []
    total_trees = sum(c[4] for c in configs)

    if verbose:
        print(f"[train_gen] Generating ~{n_samples_target} samples "
              f"from {total_trees} trees...", file=sys.stderr)

    for ci, (n, L, sub, indel, n_trees) in enumerate(configs):
        if verbose:
            print(f"\n[train_gen] Config {ci+1}/{len(configs)}: "
                  f"n={n}, L={L}, sub={sub}, indel={indel}, trees={n_trees}",
                  file=sys.stderr)

        for ti in range(n_trees):
            if len(all_samples) >= n_samples_target:
                break
            seed = ci * 1000 + ti
            try:
                true_nwk, seqs, names = generate_simulated_tree(
                    n, L, sub, indel, seed, verbose=False,
                )
                samples = generate_training_sample(
                    true_nwk, seqs, names,
                    n_samples_per_tree=20, verbose=False,
                )
                all_samples.extend(samples)

                if verbose and (ti + 1) % 5 == 0:
                    print(f"  Tree {ti+1}/{n_trees}: "
                          f"{len(samples)} samples (total: {len(all_samples)})",
                          file=sys.stderr)

            except Exception as e:
                if verbose:
                    print(f"  [train_gen] Tree {ti+1} FAILED: {e}",
                          file=sys.stderr)
                continue

        if len(all_samples) >= n_samples_target:
            break

    # Balance dataset
    n_pos = sum(1 for s in all_samples if s["label"] == 1)
    n_neg = sum(1 for s in all_samples if s["label"] == 0)
    if verbose:
        print(f"\n[train_gen] Collected {len(all_samples)} samples: "
              f"{n_pos} positive, {n_neg} negative", file=sys.stderr)

    # Save
    output_path = os.path.abspath(output_pkl)
    ensure_dir(os.path.dirname(output_path))
    with open(output_path, "wb") as f:
        pickle.dump(all_samples, f, protocol=pickle.HIGHEST_PROTOCOL)

    if verbose:
        print(f"[train_gen] Saved to {output_path}", file=sys.stderr)
        print(f"[train_gen] Total samples: {len(all_samples)}", file=sys.stderr)

    return output_path


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
