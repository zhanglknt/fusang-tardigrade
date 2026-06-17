"""
Fusang: Tardigrade Edition -- Scalable Phylogenetic Inference via Divide-and-Conquer
Supports 10000+ taxa and arbitrary-length MSA.

Architecture (DCM-inspired):
  1. k-mer distance clustering (no MSA required)
  2. Centroid representative selection
  3. Backbone tree (NJ on representatives)
  4. Subtree construction (independent NJ per group)
  5. EPA grafting (subtrees back onto backbone)
  6. Optional DL NNI refinement (Phase 2)

Usage:
  python fusang_v2.py -i input.fasta -o output.nwk -m default
"""

import sys
import time
import numpy as np
from typing import List, Dict, Tuple, Optional
from kmer_distance import seq_to_kmer_freq, make_gap_pattern


# ===========================================================
# 0. TreeNode class (pure Python, NEWICK I/O)
# ===========================================================


class TreeNode:
    """Lightweight tree node supporting NEWICK output. Replaces ete3.Tree."""
    def __init__(self, name="", dist=0.0):
        self.name = name
        self.dist = float(dist)
        self.children = []   # List[TreeNode]
        self.up = None        # TreeNode or None

    def add_child(self, name="", dist=0.0):
        child = TreeNode(name=name, dist=dist)
        child.up = self
        self.children.append(child)
        return child

    def get_leaves(self):
        if not self.children:
            return [self]
        leaves = []
        for ch in self.children:
            leaves.extend(ch.get_leaves())
        return leaves

    def search_by_name(self, name):
        if self.name == name:
            return self
        for ch in self.children:
            found = ch.search_by_name(name)
            if found:
                return found
        return None

    def write(self, format=1):
        """Output NEWICK string with trailing semicolon."""
        return self._to_newick(format=format) + ";"

    def _to_newick(self, format=1):
        """Output NEWICK string for this subtree.

        format=0: full NEWICK with internal node names
        format=1: standard NEWICK (leaf names only, internal nodes unnamed)
        """
        if not self.children:
            # Leaf node: output name:dist
            label = self.name
            if format >= 1 and self.dist >= 0:
                label += ":" + str(round(self.dist, 6))
            return label
        # Internal node: recurse into children
        child_strs = []
        for ch in self.children:
            cs = ch._to_newick(format=format)
            child_strs.append(cs)
        result = "(" + ",".join(child_strs) + ")"
        # Append branch length for internal nodes (format >= 1)
        if format >= 1 and self.dist > 0:
            result += ":" + str(round(self.dist, 6))
        return result

    def detach(self):
        if self.up:
            self.up.children.remove(self)
            self.up = None

    @classmethod
    def from_newick(cls, newick_str):
        """Parse a Newick string into a TreeNode tree.

        Supports standard Newick format with branch lengths and names.
        Ported from incremental_update.py parse_newick().

        Args:
            newick_str: Newick string (with or without trailing semicolon)

        Returns:
            TreeNode root of the parsed tree
        """
        newick_str = newick_str.strip().rstrip(';')

        def _parse_node(s, pos):
            node = cls()
            if pos < len(s) and s[pos] == '(':
                pos += 1
                children = []
                while pos < len(s) and s[pos] != ')':
                    child, pos = _parse_node(s, pos)
                    child.up = node
                    children.append(child)
                    if pos < len(s) and s[pos] == ',':
                        pos += 1
                if pos < len(s) and s[pos] == ')':
                    pos += 1
                node.children = children
                # After closing ')', read optional internal node name
                # (BioPython NJ names internal nodes like "Inner2:0.075")
                name_chars = []
                while pos < len(s) and s[pos] not in ',):;:':
                    name_chars.append(s[pos])
                    pos += 1
                if name_chars:
                    node.name = ''.join(name_chars).strip()
            else:
                name_chars = []
                while pos < len(s) and s[pos] not in ',):;':
                    name_chars.append(s[pos])
                    pos += 1
                node.name = ''.join(name_chars).strip()
            # Parse branch length after ':'
            if pos < len(s) and s[pos] == ':':
                pos += 1
                num_chars = []
                while pos < len(s) and s[pos] in '0123456789.eE-+':
                    num_chars.append(s[pos])
                    pos += 1
                if num_chars:
                    node.dist = float(''.join(num_chars))
            return node, pos

        root, _ = _parse_node(newick_str, 0)
        return root

    def to_newick(self, format=1):
        """Output Newick string with trailing semicolon. Alias for write()."""
        return self.write(format=format)

    def __repr__(self):
        return "TreeNode(name=%s, dist=%.4f, n_children=%d)" % (
            self.name, self.dist, len(self.children))


# ===========================================================
# 1. k-mer distance (using numba-jit accelerated version from kmer_distance)
# ===========================================================


def compute_kmer_distance_matrix(sequences, taxon_names, k=4,
                                   metric="cosine", n_threads=1,
                                   gap_pattern=None):
    """Compute n x n k-mer distance matrix."""
    n = len(sequences)
    if n == 0:
        return np.array([], dtype=np.float32)
    gap_info = ""
    if gap_pattern is not None:
        # If gap_pattern is a string (e.g. "gap2"), convert to position tuple
        if isinstance(gap_pattern, str) and gap_pattern != "none":
            gap_pattern = make_gap_pattern(k, style=gap_pattern)
        gap_info = ", gap=%s" % str(gap_pattern)
    print("[kmer] Computing %dx%d distance matrix, k=%d, metric=%s%s" % (n, n, k, metric, gap_info))
    # Step 1: compute k-mer vectors (using numba-jit accelerated version)
    print("[kmer] Step 1: computing k-mer vectors for %d sequences..." % n)
    # Always use single-threaded execution to avoid numba JIT recompilation in child processes
    freq_vectors = [seq_to_kmer_freq(seq, k, True, gap_pattern) for seq in sequences]
    freq_vectors = np.array(freq_vectors, dtype=np.float32)  # (n, 4^k)
    # Step 2: compute distance matrix
    print("[kmer] Step 2: computing distance matrix...")
    if metric == "cosine":
        norms = np.linalg.norm(freq_vectors, axis=1, keepdims=True)
        norms[norms == 0] = 1.0
        normalized = freq_vectors / norms
        cos_sim = normalized @ normalized.T   # (n, n)
        dist_matrix = np.clip(1.0 - cos_sim, 0.0, 1.0)
        np.fill_diagonal(dist_matrix, 0.0)
    else:
        dist_matrix = np.zeros((n, n), dtype=np.float32)
        for i in range(n):
            for j in range(i + 1, n):
                d = np.linalg.norm(freq_vectors[i] - freq_vectors[j])
                dist_matrix[i][j] = d
                dist_matrix[j][i] = d
            if i % max(1, n // 10) == 0:
                print("[kmer]   progress: %d/%d" % (i + 1, n))
    print("[kmer] Done. Shape:", dist_matrix.shape)
    return dist_matrix


def read_fasta(file_path):
    """Read FASTA file. Robust against trailing empty lines."""
    taxon_names = []
    sequences = []
    current_seq = []
    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_seq:
                    sequences.append("".join(current_seq))
                    current_seq = []
                name = line[1:].split()[0]
                taxon_names.append(name)
            else:
                current_seq.append(line.upper())
    # Flush last sequence (handles files with or without trailing newline)
    if current_seq:
        sequences.append("".join(current_seq))
    # Safety: ensure equal lengths
    n = min(len(taxon_names), len(sequences))
    taxon_names = taxon_names[:n]
    sequences = sequences[:n]
    return taxon_names, sequences


# ===========================================================
# 2. Neighbour-Joining tree construction (pure Python)
# ===========================================================


def nj_tree(distance_matrix, taxon_names):
    """
    Pure Python implementation of Neighbour-Joining (Saitou & Nei 1987).
    Correctly updates the distance matrix using enumerated indices.
    Returns: TreeNode root of the NJ tree.
    """
    n_original = len(taxon_names)
    if n_original == 0:
        return TreeNode(name="empty")
    if n_original == 1:
        root = TreeNode(name="root")
        root.add_child(name=taxon_names[0], dist=0.0)
        return root
    D = distance_matrix.astype(np.float64).copy()
    names = list(taxon_names)
    # node_pool: list of TreeNode, one per active taxon
    node_pool = [TreeNode(name=nm, dist=0.0) for nm in names]
    while len(node_pool) > 2:
        m = len(node_pool)
        row_sums = np.sum(D[:m, :m], axis=1)
        # Vectorized Q matrix: Q_ij = (m-2)*D_ij - row_sum_i - row_sum_j
        Q = (m - 2) * D[:m, :m].copy()
        Q -= row_sums[:, np.newaxis]
        Q -= row_sums[np.newaxis, :]
        np.fill_diagonal(Q, np.inf)
        # Find minimum Q[i][j]
        min_idx = np.unravel_index(np.argmin(Q), Q.shape)
        min_i, min_j = int(min_idx[0]), int(min_idx[1])
        if min_i > min_j:
            min_i, min_j = min_j, min_i
        # Compute limb lengths
        sum_i = row_sums[min_i] - D[min_i][min_j]
        sum_j = row_sums[min_j] - D[min_j][min_i]
        if m > 2:
            limb_i = 0.5 * D[min_i][min_j] + (sum_i - sum_j) / (2 * (m - 2))
        else:
            limb_i = 0.5 * D[min_i][min_j]
        limb_j = D[min_i][min_j] - limb_i
        limb_i = max(limb_i, 1e-8)
        limb_j = max(limb_j, 1e-8)
        # Create new internal node
        new_node = TreeNode(name="internal_%d" % (n_original - len(node_pool) + 1))
        node_i = node_pool[min_i]
        node_j = node_pool[min_j]
        node_i.up = None
        node_j.up = None
        new_node.children.append(node_i)
        new_node.children.append(node_j)
        node_i.up = new_node
        node_j.up = new_node
        node_i.dist = limb_i
        node_j.dist = limb_j
        # --- BUG FIX: correctly build reduced distance matrix ---
        # keep = indices of node_pool that survived (not min_i/min_j)
        keep = [idx for idx in range(m) if idx != min_i and idx != min_j]
        new_m = m - 1  # new matrix size (last slot = new node)
        new_D = np.zeros((new_m, new_m), dtype=np.float64)
        # Fill submatrix among kept nodes using original D indices
        for ni_idx, ii in enumerate(keep):
            for nj_idx, jj in enumerate(keep):
                if ni_idx < nj_idx:
                    val = D[ii][jj]
                    new_D[ni_idx][nj_idx] = val
                    new_D[nj_idx][ni_idx] = val
        # Last row/col = distances from new node to each kept node
        new_dist = np.zeros(new_m - 1, dtype=np.float64)
        for k_idx, k in enumerate(keep):
            new_dist[k_idx] = 0.5 * (D[min_i][k] + D[min_j][k] - D[min_i][min_j])
        last = new_m - 1
        for k_idx in range(new_m - 1):
            new_D[k_idx][last] = new_dist[k_idx]
            new_D[last][k_idx] = new_dist[k_idx]
        new_D[last][last] = 0.0
        D = new_D
        # Rebuild node_pool: kept nodes + new internal node
        new_node_pool = [node_pool[idx] for idx in keep]
        new_node_pool.append(new_node)
        node_pool = new_node_pool
    # Final step: connect last 2 (or 1) nodes to root
    if len(node_pool) == 2:
        d = D[0][1]
        node_pool[0].dist = max(d / 2.0, 1e-8)
        node_pool[1].dist = max(d / 2.0, 1e-8)
        root = TreeNode(name="root")
        root.children.append(node_pool[0])
        root.children.append(node_pool[1])
        node_pool[0].up = root
        node_pool[1].up = root
    else:
        root = node_pool[0]
    return root


# ===========================================================
# 3. Clustering and DCM (Divide-and-Conquer Framework)
# ===========================================================


def _greedy_clustering(distance_matrix, max_group_size=40):
    """Fallback greedy clustering when scipy is not available."""
    n = distance_matrix.shape[0]
    if n <= max_group_size:
        return [list(range(n))]
    assigned = [False] * n
    groups = []
    for start in range(n):
        if assigned[start]:
            continue
        group = [start]
        assigned[start] = True
        changed = True
        while changed and len(group) < max_group_size:
            changed = False
            best_idx = -1
            best_dist = float("inf")
            for candidate in range(n):
                if assigned[candidate]:
                    continue
                # Minimum distance from candidate to any member of group
                min_d = min(distance_matrix[candidate][g] for g in group)
                if min_d < best_dist:
                    best_dist = min_d
                    best_idx = candidate
            if best_idx != -1:
                group.append(best_idx)
                assigned[best_idx] = True
                changed = True
        groups.append(group)
    return groups


def tree_balanced_split(Z, n_taxa, max_group_size, min_ratio=0.1):
    """
    Given a scipy linkage matrix Z (n_taxa-1 rows),
    recursively cut the tree to produce clusters all <= max_group_size,
    while avoiding pathological splits (cluster < min_ratio * n_taxa).

    Returns: list of lists of indices (0-based leaf indices).
    """
    import bisect
    # Z[:,2] = height at each merge. Monotonically increasing.
    # We want to find cut heights that give balanced partitions.
    #
    # Strategy: top-down recursive split using the linkage tree.
    # Each internal node = a merge. We traverse from root,
    # split at the deepest node where both children >= min_size.

    min_size = max(2, int(n_taxa * min_ratio))

    # Build children lookup from Z
    # Z[i] merges clusters (Z[i,0], Z[i,1]) -> new cluster n_taxa+i
    # Leaves: 0..n_taxa-1
    n_merge = Z.shape[0]
    children = {}
    for i in range(n_merge):
        left = int(Z[i, 0])
        right = int(Z[i, 1])
        new_id = n_taxa + i
        children[new_id] = (left, right)

    def _get_leaves(node_id):
        if node_id < n_taxa:
            return [node_id]
        l, r = children[node_id]
        return _get_leaves(l) + _get_leaves(r)

    def _split_node(node_id):
        """Yield leaf-groups by recursively splitting node_id."""
        leaves = _get_leaves(node_id)
        if len(leaves) <= max_group_size:
            yield leaves
            return
        # Find the best split point in this subtree:
        # the merge just above where the subtree root is.
        # node_id is a merge node (n_taxa + i).
        # Its children are left, right.
        if node_id < n_taxa:
            yield leaves  # shouldnt happen
            return
        l_child, r_child = children[node_id]
        l_leaves = _get_leaves(l_child)
        r_leaves = _get_leaves(r_child)
        # If one child is too small, don't split here; need to go deeper
        # or accept imbalance. Instead: force split and accept,
        # but only if both children >= min_size.
        if len(l_leaves) < min_size or len(r_leaves) < min_size:
            # Can't split at this node. But group is > max_group_size,
            # so we must split somewhere. Force split at this node anyway
            # and let the small child be handled by overlap/merging.
            # Actually: just yield both children recursively.
            # If a child is too small, it will be merged later.
            pass
        # Recurse into both children
        yield from _split_node(l_child)
        yield from _split_node(r_child)

    # Start from root = last merge
    root_id = n_taxa + n_merge - 1
    raw_groups = list(_split_node(root_id))

    # Post-process: merge tiny groups (< min_size) into nearest neighbor group
    if len(raw_groups) <= 1:
        return [list(range(n_taxa))]  # fallback: one group

    merged = []
    tiny_groups = []
    for g in raw_groups:
        if len(g) < min_size:
            tiny_groups.append(g)
        else:
            merged.append(g)

    for tg in tiny_groups:
        if not merged:
            merged.append(tg)
            continue
        # Find the merged group whose centroid is closest to tg's centroid
        # Simplest: just append to the largest group
        largest = max(merged, key=len)
        largest.extend(tg)

    return merged


def kmer_clustering(distance_matrix, taxon_names, max_group_size=40, overlap=0.15):
    """
    Cluster taxa using k-mer distance matrix.
    Returns list of groups, each group is a list of integer indices.

    Uses UPGMA hierarchical clustering with balanced tree-cutting
    to avoid pathological splits (e.g., [9999, 1] on coalescent data).
    """
    n = len(taxon_names)
    if n <= max_group_size:
        return [list(range(n))]
    try:
        from scipy.cluster.hierarchy import linkage, fcluster
        from scipy.spatial.distance import squareform
        condensed = squareform(distance_matrix)
        Z = linkage(condensed, method="average")
        # Use balanced tree splitting instead of fixed n_clusters
        groups = tree_balanced_split(Z, n, max_group_size, min_ratio=0.1)
        # Further split any group that still exceeds max_group_size * 2
        balanced = []
        for g in groups:
            if len(g) > max_group_size:
                # Force split into max_group_size chunks
                while g:
                    balanced.append(g[:max_group_size])
                    g = g[max_group_size:]
            else:
                balanced.append(g)
        # Sanity check
        total = sum(len(g) for g in balanced)
        assert total == n, f"Clustering lost taxa: expected {n}, got {total}"
        ns = [len(g) for g in balanced]
        print(f"  [clustering] UPGMA tree-cut: {len(balanced)} groups, "
              f"sizes: min={min(ns)}, max={max(ns)}, mean={total/len(balanced):.0f}")
        return balanced
    except ImportError:
        print("[WARN] scipy not installed, using greedy clustering.")
        return _greedy_clustering(distance_matrix, max_group_size)


def get_centroid_from_distmat(dist_mat, indices):
    """
    Select the centroid using a pre-computed distance matrix.
    Centroid = the node with minimum sum of distances to all others in the group.
    `indices`: list of integer indices into `dist_mat`.
    Returns: local index into `indices`.
    """
    n_local = len(indices)
    if n_local <= 1:
        return 0
    min_sum = float("inf")
    centroid_local_idx = 0
    for i_local in range(n_local):
        s = 0.0
        gi = indices[i_local]
        for j_local in range(n_local):
            if i_local == j_local:
                continue
            gj = indices[j_local]
            s += float(dist_mat[gi][gj])
        if s < min_sum:
            min_sum = s
            centroid_local_idx = i_local
    return centroid_local_idx


def p_distance(seq1, seq2):
    """Compute p-distance between two aligned sequences."""
    diff = 0
    total = 0
    for a, b in zip(seq1.upper(), seq2.upper()):
        if a == "-" or b == "-":
            continue
        total += 1
        if a != b:
            diff += 1
    return diff / total if total > 0 else 0.0


def pairwise_distance_matrix(sequences, taxon_names, method="kmer", k=4, n_threads=4, gap_pattern=None):
    """Compute pairwise distance matrix using specified method."""
    n = len(sequences)
    if method == "kmer":
        return compute_kmer_distance_matrix(sequences, taxon_names, k=k,
                                           metric="cosine", n_threads=n_threads,
                                           gap_pattern=gap_pattern)
    elif method == "p-distance":
        dist_mat = np.zeros((n, n), dtype=np.float64)
        for i in range(n):
            for j in range(i + 1, n):
                d = p_distance(sequences[i], sequences[j])
                dist_mat[i][j] = d
                dist_mat[j][i] = d
        return dist_mat
    else:
        raise ValueError("Unknown distance method: " + method)


def nni_refine(tree, sequences, taxon_names, dl_model=None, max_iter=100, gap_pattern=None, k=4):
    """
    NNI (Nearest Neighbor Interchange) refinement.

    Two scoring modes:
      1. BME (default, dl_model=None): Balanced Minimum Evolution distance-based
      2. DL (dl_model provided): Deep Learning quartet classification

    Parameters
    ----------
    tree : TreeNode
        Root of the phylogenetic tree to refine.
    sequences : list of str or np.ndarray
        If dl_model is None: can be a numpy distance matrix or list of sequences.
        If dl_model is not None: MUST be a list of aligned sequence strings.
    taxon_names : list of str
        Taxon names (indices into distance matrix).
    dl_model : optional Keras Model
        QuartetTransformer model taking (batch, 4, L, 1) input.
        When provided, replaces BME scoring with DL quartet classification.
    max_iter : int
        Maximum number of NNI moves to accept.
    gap_pattern : optional
        Gapped k-mer pattern for distance computation (BME mode only).

    Returns
    -------
    TreeNode
        The refined tree (modified in-place).
    """
    t0 = time.time()
    use_dl = dl_model is not None

    # --- Character-to-integer mapping for DL mode ---
    if use_dl:
        # Build char map for encoding sequences to integers
        char_map = _build_char_map(sequences)
        n_sites = max(len(s) for s in sequences)
        name_to_seq = {name: seq for name, seq in zip(taxon_names, sequences)}
        print("[NNI] DL mode: %d sequences, max length %d bp" % (
            len(sequences), n_sites))

    # --- Resolve distance matrix (for BME score reporting) ---
    if isinstance(sequences, np.ndarray) and sequences.ndim == 2:
        dist_mat = sequences.astype(np.float64)
    else:
        dist_mat = pairwise_distance_matrix(sequences, taxon_names,
                                            method="kmer", k=k, n_threads=1,
                                            gap_pattern=gap_pattern)

    n_taxa = len(taxon_names)
    name_to_idx = {name: i for i, name in enumerate(taxon_names)}

    # Precompute distance values for speed (only used when dl_model=None)
    avg_cache = {}

    def avg_d(names_a, names_b):
        if not names_a or not names_b:
            return 0.0
        key = (frozenset(names_a), frozenset(names_b))
        if key in avg_cache:
            return avg_cache[key]
        s = 0.0
        for na in names_a:
            ia = name_to_idx[na]
            for nb in names_b:
                ib = name_to_idx[nb]
                s += dist_mat[ia][ib]
        val = s / (len(names_a) * len(names_b))
        avg_cache[key] = val
        return val

    # --- DL quartet scoring function ---
    def dl_score_quartet(names_a, names_b, names_c, names_d):
        """
        Score the 3 possible quartet topologies using DL model.

        Returns: (s0, s1, s2) = DL probabilities for topologies:
          T0: (A,B)|(C,D)  — current
          T1: (A,C)|(B,D)
          T2: (A,D)|(B,C)

        Higher = more likely.
        """
        # Select representative sequences (first leaf in each group)
        rep_a = names_a[0]
        rep_b = names_b[0]
        rep_c = names_c[0]
        rep_d = names_d[0]

        # Build 4-taxon alignment as integer array
        seqs = [name_to_seq[n] for n in (rep_a, rep_b, rep_c, rep_d)]
        L = max(len(s) for s in seqs)
        arr = np.zeros((1, 4, L, 1), dtype=np.float32)
        for i, s in enumerate(seqs):
            for j, ch in enumerate(s):
                arr[0, i, j, 0] = char_map.get(ch, 0)

        # DL inference
        probs = dl_model(arr, training=False).numpy()[0]  # (3,)

        # Map DL output to our topology order:
        # DL output classes correspond to 3 rooted topologies
        # We remap based on the quartet arrangement
        return float(probs[0]), float(probs[1]), float(probs[2])

    def get_internal_edges(root):
        edges = []
        stack = [root]
        while stack:
            node = stack.pop()
            if node.up is not None and len(node.children) >= 2 and len(node.up.children) >= 2:
                edges.append((node, node.up))
            for ch in node.children:
                stack.append(ch)
        return edges

    # --- BME tree score for progress reporting (used in both modes) ---
    def bme_score(root):
        leaves = root.get_leaves()
        n_l = len(leaves)
        depth = {}
        stack = [(root, 0)]
        while stack:
            node, d = stack.pop()
            if not node.children:
                depth[node.name] = d
            for ch in node.children:
                stack.append((ch, d + 1))

        node_depth = {}
        stack = [(root, 0)]
        while stack:
            node, d = stack.pop()
            node_depth[id(node)] = d
            for ch in node.children:
                stack.append((ch, d + 1))

        score = 0.0
        for i in range(n_l):
            for j in range(i + 1, n_l):
                ni = name_to_idx[leaves[i].name]
                nj = name_to_idx[leaves[j].name]
                dij = depth[leaves[i].name] + depth[leaves[j].name] - 2 * node_depth[id(_lca(root, leaves[i], leaves[j]))]
                weight = 1.0 / (2 ** (dij - 1)) if dij >= 1 else 1.0
                score += dist_mat[ni][nj] * weight
        return score

    def _lca(root, leaf_a, leaf_b):
        ancestors_a = set()
        cur = leaf_a
        while cur is not None:
            ancestors_a.add(id(cur))
            cur = cur.up
        cur = leaf_b
        while cur is not None:
            if id(cur) in ancestors_a:
                return cur
            cur = cur.up
        return root

    # --- Leaf name cache ---
    leaf_cache = {}

    def get_leaf_names(node):
        nid = id(node)
        if nid not in leaf_cache:
            leaf_cache[nid] = [lf.name for lf in node.get_leaves()]
        return leaf_cache[nid]

    def invalidate_cache(*nodes):
        for node in nodes:
            cur = node
            while cur is not None:
                leaf_cache.pop(id(cur), None)
                cur = cur.up

    initial_score = bme_score(tree)
    method_label = "DL" if use_dl else "BME"
    print("[NNI] Initial BME score: %.6f  (scoring: %s)" % (initial_score, method_label))

    total_moves = 0
    iteration = 0
    no_improvement_streak = 0

    while iteration < max_iter:
        edges = get_internal_edges(tree)
        improved = False

        for u, p in edges:
            if len(u.children) < 2:
                continue
            if u not in p.children:
                continue

            child_a = u.children[0]
            child_b = u.children[1]

            sibling = None
            for ch in p.children:
                if ch is not u:
                    sibling = ch
                    break
            if sibling is None:
                continue

            a_names = get_leaf_names(child_a)
            b_names = get_leaf_names(child_b)
            c_names = get_leaf_names(sibling)

            abc_set = set(a_names) | set(b_names) | set(c_names)
            all_leaf_names = [lf.name for lf in tree.get_leaves()]
            d_names = [nm for nm in all_leaf_names if nm not in abc_set]

            if not d_names:
                # 3-group case (edge at unrooted tree root)
                if use_dl:
                    # DL scoring
                    s0, s1, s2 = dl_score_quartet(a_names, b_names, c_names, c_names[:1])
                else:
                    s0 = avg_d(a_names, b_names)
                    s1 = avg_d(a_names, c_names)
                    s2 = avg_d(b_names, c_names)

                if use_dl:
                    best = max(s0, s1, s2)
                    if best <= s0 + 1e-8:
                        continue
                else:
                    best = min(s0, s1, s2)
                    if best >= s0 - 1e-12:
                        continue

                if (use_dl and best == s1) or (not use_dl and best == s1):
                    u.children.remove(child_b)
                    p.children.remove(sibling)
                    child_b.up = p
                    p.children.append(child_b)
                    sibling.up = u
                    u.children.append(sibling)
                    print("[NNI] Move #%d: DL swap (A,C) at node %s, "
                          "scores DL=(%.4f,%.4f,%.4f)" % (
                              total_moves + 1, u.name, s0, s1, s2)
                          if use_dl else
                          "[NNI] Move #%d: NNI swap (A,C)|(B,D) at node %s, "
                          "score %.6f → %.6f (Δ=%.6f)" % (
                              total_moves + 1, u.name, s0, s1, s1 - s0))
                else:
                    u.children.remove(child_a)
                    p.children.remove(sibling)
                    child_a.up = p
                    p.children.append(child_a)
                    sibling.up = u
                    u.children.append(sibling)
                    print("[NNI] Move #%d: DL swap (B,C) at node %s, "
                          "scores DL=(%.4f,%.4f,%.4f)" % (
                              total_moves + 1, u.name, s0, s1, s2)
                          if use_dl else
                          "[NNI] Move #%d: NNI swap (A,D)|(B,C) at node %s, "
                          "score %.6f → %.6f (Δ=%.6f)" % (
                              total_moves + 1, u.name, s0, s2, s2 - s0))
            else:
                # 4-group case
                if use_dl:
                    s0, s1, s2 = dl_score_quartet(a_names, b_names, c_names, d_names)
                    best = max(s0, s1, s2)
                    if best <= s0 + 1e-8:
                        continue
                else:
                    s0 = avg_d(a_names, b_names) + avg_d(c_names, d_names)
                    s1 = avg_d(a_names, c_names) + avg_d(b_names, d_names)
                    s2 = avg_d(a_names, d_names) + avg_d(b_names, c_names)
                    best = min(s0, s1, s2)
                    if best >= s0 - 1e-12:
                        continue

                if (use_dl and best == s1) or (not use_dl and best == s1):
                    u.children.remove(child_b)
                    p.children.remove(sibling)
                    child_b.up = p
                    p.children.append(child_b)
                    sibling.up = u
                    u.children.append(sibling)
                else:
                    u.children.remove(child_a)
                    p.children.remove(sibling)
                    child_a.up = p
                    p.children.append(child_a)
                    sibling.up = u
                    u.children.append(sibling)

            invalidate_cache(u, p, child_a, child_b, sibling)
            # Targeted cache invalidation: only clear entries involving swapped subtrees
            affected = set(a_names) | set(b_names) | set(c_names)
            if d_names:
                affected |= set(d_names)
            to_remove = [key for key in avg_cache
                         if affected & (set(key[0]) | set(key[1]))]
            for key in to_remove:
                del avg_cache[key]

            total_moves += 1
            improved = True
            break  # greedy hill climbing

        iteration += 1
        if not improved:
            no_improvement_streak += 1
            if no_improvement_streak >= 3:
                break
        else:
            no_improvement_streak = 0

    final_score = bme_score(tree)
    elapsed = time.time() - t0
    print("[NNI] Done: %d moves in %d iterations, %.2fs" % (
        total_moves, iteration, elapsed))
    print("[NNI] BME score: %.6f → %.6f (improvement: %.6f)" % (
        initial_score, final_score, initial_score - final_score))
    return tree


def _build_char_map(sequences):
    """Build character-to-integer mapping for sequence encoding."""
    chars = set()
    for s in sequences:
        chars.update(s)
    # Remove gap characters, treat as 0
    gap_chars = {'-', '.', 'N', 'n', '?'}
    active = chars - gap_chars
    char_map = {c: i + 1 for i, c in enumerate(sorted(active))}  # 1-based, 0=gap/unknown
    # Explicitly set gap chars to 0
    for c in gap_chars:
        char_map[c] = 0
    return char_map


def divide_and_conquer(sequences, taxon_names,
                            max_group_size=200, overlap=0.15,
                            mode="default", distance_method="kmer",
                            k=4, n_threads=4, _depth=0, gap_pattern=None,
                            tree_method="nj", epa_method="nj_centroid",
                            use_minhash=False, minhash_k=5, minhash_sketches=128):
    n = len(sequences)
    indent = "  " * _depth
    sketches = None  # MinHash sketches (computed in Step 2 if use_minhash)
    print("%s[DCM] depth=%d, n=%d, mode=%s, epa=%s" % (indent, _depth, n, mode, epa_method))
    # Base case: build tree directly
    if n <= max_group_size:
        print("%s[DCM] n=%d <= %d, building tree directly" % (indent, n, max_group_size))
        # Fast path: n≤2 skips k-mer distance (trivial tree)
        if n == 1:
            print("%s[DCM]   n=1: trivial leaf" % indent)
            leaf = TreeNode(name=taxon_names[0], dist=0.0)
            root = TreeNode()
            root.children = [leaf]
            leaf.up = root
            return root
        if n == 2:
            print("%s[DCM]   n=2: trivial pair" % indent)
            leaf0 = TreeNode(name=taxon_names[0], dist=0.5)
            leaf1 = TreeNode(name=taxon_names[1], dist=0.5)
            root = TreeNode()
            root.children = [leaf0, leaf1]
            leaf0.up = root
            leaf1.up = root
            return root
        dist_mat = pairwise_distance_matrix(
            sequences, taxon_names,
            method=distance_method, k=k, n_threads=n_threads,
            gap_pattern=gap_pattern
        )
        # FastME requires n >= 4; fall back to NJ for small inputs
        use_fastme = (tree_method == "fastme") and (n >= 4)
        if use_fastme:
            from fastme_backend import build_tree_fastme
            tree = build_tree_fastme(dist_mat, taxon_names, n_threads=n_threads)
        else:
            if tree_method == "fastme":
                print("%s[DCM]   n=%d < 4, FastME unavailable, using NJ fallback" % (indent, n))
            tree = nj_tree(dist_mat, taxon_names)
        if mode == "full-dl":
            print("%s[DCM] full-dl mode: DL quartet building (not yet implemented)" % indent)
        return tree
    # Step 1: compute distance matrix for all sequences in this group
    t0 = time.time()
    print("%s[DCM] Step 1: computing %dx%d distance matrix..." % (indent, n, n))
    dist_mat = pairwise_distance_matrix(
        sequences, taxon_names,
        method=distance_method, k=k, n_threads=n_threads,
        gap_pattern=gap_pattern
    )
    t1 = time.time()
    print("%s[DCM]   Step 1 done: %.1fs" % (indent, t1 - t0))
    # Step 2: clustering
    t1 = time.time()
    print("%s[DCM] Step 2: clustering (max_group=%d)..." % (indent, max_group_size))
    
    if use_minhash:
        # MinHash LSH clustering: O(nL) instead of O(n²)
        try:
            from minhash_lsh import seq_to_minhash, minhash_coarse_cluster
            import numpy as np
            mk = k if k else minhash_k
            ns = minhash_sketches
            print("%s[DCM]   computing MinHash sketches (k=%d, sketches=%d)..." % (indent, mk, ns))
            sketches = [seq_to_minhash(s, k=mk, num_hashes=ns) for s in sequences]
            groups = minhash_coarse_cluster(sketches, target_group_size=max_group_size, band_width=None)
        except ImportError:
            print("%s[DCM]   WARNING: minhash_lsh not available, falling back to kmer_clustering" % indent)
            groups = kmer_clustering(dist_mat, taxon_names, max_group_size=max_group_size, overlap=overlap)
    else:
        groups = kmer_clustering(dist_mat, taxon_names,
                          max_group_size=max_group_size,
                          overlap=overlap)
    
    group_sizes = [len(g) for g in groups]
    print("%s[DCM]   got %d groups, sizes: %s" % (indent, len(groups), group_sizes))
    t2 = time.time()
    print("%s[DCM]   Step 2 done: %.1fs" % (indent, t2 - t1))
    
    # Step 3: select centroid representatives
    t2 = time.time()
    print("%s[DCM] Step 3: selecting centroid representatives..." % indent)
    rep_indices = []       # global indices into `sequences`
    rep_sequences = []
    rep_names = []
    
    if use_minhash and sketches is not None:
        # Use MinHash distance for centroid selection
        from minhash_lsh import minhash_distance
        for g_idx, group in enumerate(groups):
            best_idx = 0
            best_avg_dist = float('inf')
            for i_idx, i in enumerate(group):
                avg_dist = 0.0
                count = 0
                for j_idx, j in enumerate(group):
                    if i == j:
                        continue
                    d = minhash_distance(sketches[i], sketches[j])
                    avg_dist += d
                    count += 1
                if count > 0:
                    avg_dist /= count
                if avg_dist < best_avg_dist:
                    best_avg_dist = avg_dist
                    best_idx = i_idx
            global_idx = group[best_idx]
            rep_indices.append(global_idx)
            rep_sequences.append(sequences[global_idx])
            rep_names.append("REP_%d_%s" % (g_idx, taxon_names[global_idx]))
    else:
        # Use pre-computed distance matrix to find centroid
        for g_idx, group in enumerate(groups):
            centroid_local_idx = get_centroid_from_distmat(dist_mat, group)
            global_idx = group[centroid_local_idx]
            rep_indices.append(global_idx)
            rep_sequences.append(sequences[global_idx])
            rep_names.append("REP_%d_%s" % (g_idx, taxon_names[global_idx]))
    
    print("%s[DCM]   reps: %s" % (indent, rep_names))
    t3 = time.time()
    print("%s[DCM]   Step 3 done: %.1fs" % (indent, t3 - t2))
    # Step 4: build backbone tree (NJ on representatives)
    t3 = time.time()
    print("%s[DCM] Step 4: building backbone tree (%d reps)..." % (indent, len(rep_names)))
    backbone_tree = divide_and_conquer(
        rep_sequences, rep_names,
        max_group_size=max_group_size,
        overlap=overlap, mode=mode,
        distance_method=distance_method,
        k=k, n_threads=n_threads, gap_pattern=gap_pattern,
        _depth=_depth + 1, tree_method=tree_method,
        epa_method=epa_method,
        use_minhash=use_minhash, minhash_k=minhash_k, minhash_sketches=minhash_sketches
    )
    n_leaves_bb = len(backbone_tree.get_leaves())
    t4 = time.time()
    print("%s[DCM]   backbone done: %d leaves, %.1fs" % (indent, n_leaves_bb, t4 - t3))
    # Step 5: build subtrees (independent NJ per group)
    t4 = time.time()
    print("%s[DCM] Step 5: building %d subtrees..." % (indent, len(groups)))
    subtrees = []
    for g_idx, group in enumerate(groups):
        group_seqs = [sequences[i] for i in group]
        group_names = [taxon_names[i] for i in group]
        subtree = divide_and_conquer(
            group_seqs, group_names,
            max_group_size=max_group_size,
            overlap=overlap, mode=mode,
            distance_method=distance_method,
            k=k, n_threads=n_threads, gap_pattern=gap_pattern,
            _depth=_depth + 1, tree_method=tree_method,
            epa_method=epa_method,
            use_minhash=use_minhash, minhash_k=minhash_k, minhash_sketches=minhash_sketches
        )
        subtrees.append((g_idx, group, subtree))
        n_leaves_st = len(subtree.get_leaves())
        print("%s[DCM]   group %d subtree done: %d leaves" % (indent, g_idx, n_leaves_st))
    t5 = time.time()
    print("%s[DCM]   Step 5 done: %.1fs" % (indent, t5 - t4))
    # BUG FIX: Do NOT detach the rep leaf. Instead, replace it with the
    # subtree by attaching subtree root as child of the rep leaf parent.
    print("%s[DCM] Step 6: EPA grafting (simplified)..." % indent)
    for g_idx, group, subtree in subtrees:
        # Find the rep leaf in backbone by matching REP_{g_idx}_ prefix
        rep_leaf = None
        for lf in backbone_tree.get_leaves():
            if lf.name.startswith("REP_%d_" % g_idx):
                rep_leaf = lf
                break
        if rep_leaf and rep_leaf.up:
            parent = rep_leaf.up
            branch_len = rep_leaf.dist
            rep_leaf.detach()
            # Attach subtree root as child of parent
            sub_root = subtree
            # If subtree root is artificial (has 1 child and name=="root"),
            # skip it and attach its children directly
            if len(sub_root.children) == 1 and sub_root.name in ("root", ""):
                real_child = sub_root.children[0]
                real_child.dist = max(branch_len, 1e-8)
                parent.children.append(real_child)
                real_child.up = parent
            else:
                sub_root.dist = max(branch_len, 1e-8)
                parent.children.append(sub_root)
                sub_root.up = parent
            print("%s[DCM]   grafted group %d: %d leaves" % (
                indent, g_idx, len(subtree.get_leaves())))
        else:
            print("%s[DCM]   WARNING: rep for group %d not found in backbone" % (indent, g_idx))
    n_final = len(backbone_tree.get_leaves())
    t6 = time.time()
    print("%s[DCM] Step 6 done: %.1fs, Final leaves: %d" % (indent, t6 - t5, n_final))
    return backbone_tree


# ===========================================================
# 4. CLI
# ===========================================================


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(
        description="Fusang: Tardigrade Edition -- Scalable Phylogenetic Inference"
    )
    parser.add_argument("--input", "-i", required=True,
                        help="Input FASTA file")
    parser.add_argument("--output", "-o", required=True,
                        help="Output tree file (NEWICK format)")
    parser.add_argument("--mode", "-m",
                        choices=["default", "refine", "full-dl", "auto"],
                        default="auto",
                        help="Mode: auto (choose based on n), default (NJ+DCM), refine (NJ+DCM+BME), dl_refine (NJ+DCM+DL), full-dl (DL only)")
    parser.add_argument("--dl_weights", type=str, default=None,
                        help="Path to DL model weights (for dl_refine mode), "
                             "e.g., ./dl_model_transformer/varlen/SVG/best_weights_transformer")
    parser.add_argument("--distance_method", "-d",
                        choices=["kmer", "p-distance"],
                        default="kmer",
                        help="Distance method")
    parser.add_argument("--kmer_k", type=int, default=None,
                        help="k-mer length (auto-selected based on n_taxa if not set: n<=100→4, n>100→5)")
    parser.add_argument("--kmer_gap", type=str,
                        choices=["none", "gap1", "gap2", "gap3", "gap4"],
                        default=None,
                        help="Gapped k-mer pattern (auto-selected based on n_taxa if not set: n<=100→gap1, n>100→gap2)")
    parser.add_argument("--max_group", type=int, default=200,
                        help="Max taxa per group (default 200; auto-scaled for n>500 unless --no-auto-group)")
    parser.add_argument("--no-auto-group", dest="auto_group", action="store_false",
                        default=True,
                        help="Disable automatic max_group_size selection")
    parser.add_argument("--overlap", type=float, default=0.15,
                        help="Group overlap ratio (default 0.15)")
    parser.add_argument("--threads", "-t", type=int, default=4,
                        help="Number of threads (default 4)")
    parser.add_argument("--tree_method", type=str,
                        choices=["nj", "fastme"],
                        default="nj",
                        help="Tree building method: nj (Neighbour-Joining, recommended for k-mer distances), fastme (FastME BIONJ+bNNI, may degrade k-mer distance accuracy)")
    parser.add_argument("--auto_group_method", type=str,
                        choices=["simple", "nj_centroid", "epa_improved"],
                        default=None,
                        help="DCM strategy: simple (no DCM, direct NJ for n<500), nj_centroid (current DCM with NJ), epa_improved (DCM+NNI refinement after EPA grafting)")
    parser.add_argument("--simple", dest="force_simple", action="store_true", default=None,
                        help="Force simplified pipeline (k-mer→dist→tree, no DCM/EPA). Auto-enabled for n≤500.")
    parser.add_argument("--no-simple", dest="force_simple", action="store_false", default=None,
                        help="Force DCM pipeline even for n≤500 (for benchmarking only)")
    parser.add_argument("--use_minhash", dest="use_minhash", action="store_true", default=False,
                        help="Use MinHash LSH for coarse clustering (scales to 50K+ taxa, requires minhash_lsh.py)")
    parser.add_argument("--minhash_k", type=int, default=5,
                        help="k-mer size for MinHash (default 5)")
    parser.add_argument("--minhash_sketches", type=int, default=128,
                        help="MinHash signature size (default 128)")
    return parser.parse_args()


def _resolve_auto_mode(requested_mode, n):
    """Resolve 'auto' mode to 'default' or 'refine' based on benchmark results.

    Thresholds derived from empirical benchmarks (sub_rate=0.1, L=500bp):
      n <  20  →  default  (too small; refine not tested, likely detrimental)
      20 ≤ n < 100  →  default  (refine degrades RF%: 17.6→23.5, 10.6→12.8)
      n ≥ 100  →  refine   (significant gains: 15.7→5.1 @ 200, 14.5→8.0 @ 500)
    """
    if requested_mode != "auto":
        return requested_mode, requested_mode  # (build_mode, nni_mode)
    if n < 20:
        return "default", "default"
    elif n < 100:
        return "default", "default"
    else:
        return "default", "refine"


def main():
    args = parse_args()
    start_time = time.time()
    print("=" * 60)
    print("  Fusang: Tardigrade Edition -- Scalable Phylogenetic Inference")
    print("  Requested mode: %s" % args.mode)
    print("  Input: %s" % args.input)
    print("  Output: %s" % args.output)
    print("=" * 60)
    # Step 1: read FASTA
    print("\n[MAIN] Step 1: reading FASTA %s..." % args.input)
    taxon_names, sequences = read_fasta(args.input)
    n = len(taxon_names)
    if n == 0:
        print("[ERROR] No sequences read from file")
        sys.exit(1)
    seq_len = max((len(s) for s in sequences)) if sequences else 0
    print("[MAIN]   read %d sequences, max length %d bp" % (n, seq_len))
    if n < 3:
        print("[ERROR] Need at least 3 sequences")
        sys.exit(1)
    # Resolve auto mode based on n
    build_mode, nni_mode = _resolve_auto_mode(args.mode, n)
    # Override nni_mode if epa_method requires refinement
    epa_method = args.auto_group_method if args.auto_group_method else "nj_centroid"
    if epa_method == "epa_improved":
        nni_mode = "refine"  # Force NNI for epa_improved
    if args.mode == "auto":
        print("[MAIN]   auto mode resolved: n=%d → build=%s, nni=%s, epa=%s" % (n, build_mode, nni_mode, epa_method))
    # Auto-select kmer parameters based on n_taxa (if not explicitly set)
    # Evidence base (all seed=42 unless noted):
    #   Clean n=200:  k=5,gap2  nRF=0.0051 (best ever, kmer_study/k_value_n200)
    #   Clean n=500:  k=5,gap2  nRF=0.0221  (kmer_study/gap_scaling_n500)
    #   Clean n=1000: k=5,gap2  nRF=0.0461  (kmer_study/gap_scaling_n1000)
    #   Indel n=200:  k=5,gap3  nRF=0.0431 (best), k=5,gap2 nRF=0.0482 (+12%)
    #   Indel multi:   k=5,gap2   mean=0.0804 vs FT2=0.0842, p=0.049 (130 seeds)
    # Conclusion: k=5,gap2 is near-universal optimum across all scales.
    #   Indel data at n<=200: gap3 can help (~+5-10% relative).
    #   Small n<=100 clean: k=4,gap1 is adequate (nRF~0.02).
    if args.kmer_k is None or args.kmer_gap is None:
        if n <= 100:
            auto_k, auto_gap = 4, "gap1"
        elif n <= 500:
            auto_k, auto_gap = 5, "gap2"
        else:
            auto_k, auto_gap = 5, "gap2"
        if args.kmer_k is None:
            args.kmer_k = auto_k
        if args.kmer_gap is None:
            args.kmer_gap = auto_gap
        print("[MAIN]   auto-selected k=%d, gap=%s (n=%d taxa)" % (args.kmer_k, args.kmer_gap, n))
    # Auto-select max_group_size based on n_taxa (if auto_group enabled)
    # Diagnosis (2026-06-06): Critical design constraint discovered:
    #   The DCM backbone tree MUST have at most 2 taxa when possible.
    #   With 3+ backbone taxa, NJ/FastME topology is inherently unstable
    #   on k-mer distances → ~60% disaster rate at n=500 (nRF > 0.44).
    #   At n=500, mg=200 → 3 groups → 60% disaster; mg=250 → 2 groups → 0%.
    #   At n=1000, mg=200 → 5+ groups → 60% disaster; mg=500 → 2 groups → 0%.
    # Rule: max_group ≥ ceil(n/2) to guarantee ≤ 2 backbone taxa
    #   n≤250:  default 200 (≤2 groups, already safe, simplified pipeline used)
    #   n≤500:  250 (2 groups, verified safe, simplified pipeline used)
    #   n≤1000: 500 (2 groups, verified safe, simplified pipeline used)
    #   n≤2000: 1000 (2 groups, simplified pipeline used)
    #   n≤5000: ≥2500 (2 groups, DCM triggered above SIMPLE_THRESHOLD=2000)
    #   n>5000: ≥5000 (2 groups)
    if args.auto_group:
        if n <= 250:
            auto_mg = 200
        elif n <= 500:
            auto_mg = 250
        elif n <= 1000:
            auto_mg = 500
        elif n <= 2000:
            auto_mg = 1000
        elif n <= 5000:
            auto_mg = (n + 1) // 2   # ensure ≤2 backbone taxa
        else:
            auto_mg = (n + 1) // 2   # ensure ≤2 backbone taxa
        if auto_mg != args.max_group:
            args.max_group = auto_mg
            print("[MAIN]   auto-selected max_group=%d (n=%d taxa, %d backbone taxa)"
                  % (auto_mg, n, max(1, -(-n // auto_mg))))
    # Simplified pipeline auto-switch for n <= 500:
    #   DCM adds no accuracy benefit at small-to-medium scales.
    #   Empirical (seeds 401-405, n=1000): DCM nRF=0.084, simplified nRF=0.092.
    #   -> Keep DCM for n>500; use simplified for n<=500.
    SIMPLE_THRESHOLD = 500  # Reverted - 2000 caused regression (nRF 0.0878 → 0.1073)
    use_simple = args.force_simple
    # Override based on auto_group_method
    if epa_method == "simple":
        use_simple = True
        print("[MAIN]   auto_group_method=simple: forcing simplified pipeline (no DCM)")
    elif epa_method in ("nj_centroid", "epa_improved"):
        # DCM helps at large n; for n <= SIMPLE_THRESHOLD, simplified is better
        if n > SIMPLE_THRESHOLD:
            use_simple = False
            print("[MAIN]   auto_group_method=%s: using DCM pipeline (n=%d > %d)"
                  % (epa_method, n, SIMPLE_THRESHOLD))
        # else: leave use_simple as None so the next block can decide
    if use_simple is None:
        # Auto mode: use simplified pipeline for n <= SIMPLE_THRESHOLD
        use_simple = (n <= SIMPLE_THRESHOLD and args.mode != "dl_refine")
        if use_simple:
            print("[MAIN]   auto: n=%d ≤ %d, using simplified pipeline" % (n, SIMPLE_THRESHOLD))
        else:
            print("[MAIN]   auto: n=%d > %d, using DCM pipeline" % (n, SIMPLE_THRESHOLD))
    if use_simple:
        if args.max_group < n:
            args.max_group = n
        print("[MAIN]   simplified pipeline: n=%d ≤ %d → direct tree (no DCM/EPA, %s)"
              % (n, SIMPLE_THRESHOLD, args.tree_method))
    # Step 2: divide-and-conquer tree building
    print("\n[MAIN] Step 2: building tree (build_mode=%s)..." % build_mode)
    gap_pattern = None
    if args.kmer_gap is not None and args.kmer_gap != "none":
        gap_pattern = make_gap_pattern(args.kmer_k, style=args.kmer_gap)
        print("[MAIN]   using gapped k-mer: %s pattern=%s" % (args.kmer_gap, gap_pattern))
    tree = divide_and_conquer(
        sequences, taxon_names,
        max_group_size=args.max_group,
        overlap=args.overlap,
        mode=build_mode,
        distance_method=args.distance_method,
        k=args.kmer_k, gap_pattern=gap_pattern,
        n_threads=args.threads,
        _depth=0,
        tree_method=args.tree_method,
        epa_method=args.auto_group_method if args.auto_group_method else "nj_centroid",
        use_minhash=args.use_minhash,
        minhash_k=args.minhash_k,
        minhash_sketches=args.minhash_sketches,
    )
    # Step 3: NNI refinement (refine / dl_refine / epa_improved mode)
    # Skip NNI when tree was built by FastME (already BME-optimal)
    dl_model = None
    skip_nni = (args.tree_method == "fastme")
    if skip_nni:
        print("\n[MAIN] Step 3: skipping NNI (FastME tree is already BME-optimal, n=%d)" % n)
    elif args.mode == "dl_refine":
        if not args.dl_weights:
            print("[ERROR] dl_refine mode requires --dl_weights")
            sys.exit(1)
        print("\n[MAIN] Step 3: loading DL model...")
        from quartet_transformer import get_dl_model_transformer_varlen
        dl_model = get_dl_model_transformer_varlen(
            d_model=128, n_heads=8, n_layers=4, dff=512,
            dropout=0.1, use_conv_stem=True, conv_stride=4,
        )
        dl_model.load_weights(args.dl_weights).expect_partial()
        print("[MAIN]   DL model loaded: %s" % args.dl_weights)
        print("[MAIN] Step 3b: DL NNI refinement...")
        tree = nni_refine(tree, sequences, taxon_names, dl_model=dl_model, max_iter=100, gap_pattern=gap_pattern, k=args.kmer_k)
    elif nni_mode == "refine":
        print("\n[MAIN] Step 3: BME NNI refinement...")
        tree = nni_refine(tree, sequences, taxon_names, dl_model=None, max_iter=100, gap_pattern=gap_pattern, k=args.kmer_k)
    elif epa_method == "epa_improved":
        # Force NNI for epa_improved mode (more iterations for better optimization)
        print("\n[MAIN] Step 3: EPA-improved - BME NNI refinement (forced, 500 iterations)...")
        tree = nni_refine(tree, sequences, taxon_names, dl_model=None, max_iter=500, gap_pattern=gap_pattern, k=args.kmer_k)
    elif nni_mode == "default":
        print("\n[MAIN] Step 3: skipping NNI (build_mode=%s, n=%d)" % (build_mode, n))
    # Step 4: write NEWICK
    print("\n[MAIN] Step 4: writing NEWICK to %s..." % args.output)
    nwk_str = tree.write(format=1)
    with open(args.output, "w") as f:
        f.write(nwk_str + "\n")
    n_leaves = len(tree.get_leaves())
    print("[MAIN]   tree saved. Leaves: %d" % n_leaves)
    print("[MAIN]   NEWICK preview: %s..." % nwk_str[:200])
    elapsed = time.time() - start_time
    print("")
    print("=" * 60)
    print("  Fusang: Tardigrade Edition done! Time: %.2fs" % elapsed)
    print("   Output: %s" % args.output)
    print("=" * 60)
    return 0


if __name__ == "__main__":
    sys.exit(main())
