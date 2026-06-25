"""
tree_simulation.py
==================
轻量级系统发育树模拟模块（无 TensorFlow 依赖）。

从 validate_distill.py 提取的核心函数：
  - Node: 树节点类
  - make_coalescent_tree(): Kingman coalescent 随机树生成
  - simulate_seqs(): 沿树 JC69 模型序列进化模拟
  - get_quartet_topology(): 从树拓扑判断四元组标签

编码约定：
  0=A, 1=T, 2=C, 3=G (与 simulate_seqs 输出一致)

作者：Fusang v3 项目
日期：2026-06-13
"""

import random
import numpy as np


class Node:
    """树节点"""
    __slots__ = ('idx', 'is_leaf', 'name', 'left', 'right', 'parent', 'branch_length')

    def __init__(self, idx, is_leaf, name=None, left=None, right=None):
        self.idx = idx
        self.is_leaf = is_leaf
        self.name = name
        self.left = left
        self.right = right
        self.parent = None
        self.branch_length = 0.0

    def to_newick(self):
        """递归生成 Newick 格式字符串（不含末尾分号）。"""
        if self.is_leaf:
            name = self.name if self.name else f"t{self.idx:04d}"
            if self.branch_length > 0:
                return f"{name}:{self.branch_length:.6f}"
            return name
        left_str = self.left.to_newick() if self.left else ""
        right_str = self.right.to_newick() if self.right else ""
        inner = f"({left_str},{right_str})"
        if self.name:
            inner += self.name
        if self.branch_length > 0:
            inner += f":{self.branch_length:.6f}"
        return inner


def make_coalescent_tree(n_taxa: int, seed: int):
    """
    生成 Kingman coalescent 随机系统发育树。

    参数：
        n_taxa: 叶节点（taxa）数量
        seed: 随机种子

    返回：
        (root, nodes) — 根节点和所有节点的列表
    """
    rng = random.Random(seed)
    nodes = []
    lineages = []

    for i in range(n_taxa):
        node = Node(idx=i, is_leaf=True, name=f"t{i+1:04d}")
        nodes.append(node)
        lineages.append(node)

    next_idx = n_taxa
    while len(lineages) > 1:
        i, j = rng.sample(range(len(lineages)), 2)
        a, b = lineages[i], lineages[j]
        for k in sorted([i, j], reverse=True):
            lineages.pop(k)
        parent = Node(idx=next_idx, is_leaf=False, left=a, right=b)
        a.parent = parent
        b.parent = parent
        a.branch_length = rng.expovariate(10.0)  # mean=0.1
        b.branch_length = rng.expovariate(10.0)
        nodes.append(parent)
        lineages.append(parent)
        next_idx += 1

    return nodes[-1], nodes


def simulate_seqs(root, n_taxa: int, L: int, mu: float, seed: int,
                  indel_rate: float = 0.0, indel_seed_offset: int = 100000) -> np.ndarray:
    """
    沿系统发育树模拟 JC69 序列进化，支持可选的 indel 模拟。

    参数：
        root: 树的根节点
        n_taxa: taxa 数量
        L: 序列长度
        mu: 突变率（JC69 模型）
        seed: 随机种子
        indel_rate: indel 发生率（每序列每碱基，默认 0 不模拟 indel）
        indel_seed_offset: indel 模拟的种子偏移量

    返回：
        leaf_seqs: (n_taxa, L) int8 数组，编码 0=A, 1=T, 2=C, 3=G
    """
    rng = np.random.RandomState(seed)
    n_nodes = 2 * n_taxa - 1
    seqs = -np.ones((n_nodes, L), dtype=np.int8)

    # Root sequence
    seqs[root.idx] = rng.randint(0, 4, size=L)

    # Post-order: propagate from root to leaves
    def dfs(node):
        if node is None:
            return
        if node.left:
            br_len = node.left.branch_length
            prob = 0.25 * (1 - np.exp(-4.0 / 3.0 * mu * br_len))
            parent_seq = seqs[node.idx]
            child_seq = parent_seq.copy()
            mask = rng.random(L) < prob
            child_seq[mask] = rng.randint(0, 4, size=mask.sum())
            seqs[node.left.idx] = child_seq
            dfs(node.left)
        if node.right:
            br_len = node.right.branch_length
            prob = 0.25 * (1 - np.exp(-4.0 / 3.0 * mu * br_len))
            parent_seq = seqs[node.idx]
            child_seq = parent_seq.copy()
            mask = rng.random(L) < prob
            child_seq[mask] = rng.randint(0, 4, size=mask.sum())
            seqs[node.right.idx] = child_seq
            dfs(node.right)

    dfs(root)
    leaf_seqs = np.array([seqs[i] for i in range(n_taxa)], dtype=np.int8)

    # --- Indel simulation (post-hoc, per-leaf) ---
    if indel_rate > 0.0:
        indel_rng = np.random.RandomState(seed + indel_seed_offset)
        for i in range(n_taxa):
            leaf_seqs[i] = _apply_indels_to_sequence(
                leaf_seqs[i], L, indel_rate, indel_rng)

    return leaf_seqs


def _apply_indels_to_sequence(seq: np.ndarray, target_L: int,
                               indel_rate: float, rng: np.random.RandomState) -> np.ndarray:
    """对单条序列施加独立 indel，然后裁剪/填充到 target_L。"""
    seq_list = list(seq)

    n_events = rng.poisson(indel_rate * target_L)
    for _ in range(n_events):
        if not seq_list:
            break
        pos = rng.randint(0, len(seq_list))
        event_len = rng.randint(1, 6)
        if rng.random() < 0.5:
            new_bases = rng.randint(0, 4, size=event_len).tolist()
            for j, b in enumerate(new_bases):
                seq_list.insert(pos + j, b)
        else:
            end = min(pos + event_len, len(seq_list))
            del seq_list[pos:end]

    if len(seq_list) > target_L:
        seq_list = seq_list[:target_L]
    elif len(seq_list) < target_L:
        pad = rng.randint(0, 4, size=target_L - len(seq_list)).tolist()
        seq_list.extend(pad)

    return np.array(seq_list, dtype=np.int8)


def get_quartet_topology(root, a: int, b: int, c: int, d: int) -> int:
    """
    从树拓扑判断四元组的真实拓扑标签。

    参数：
        root: 树的根节点
        a, b, c, d: 4 个叶节点的索引

    返回：
        0 = (a,b|c,d), 1 = (a,c|b,d), 2 = (a,d|b,c)
    """
    def path_to_root(node, target_idx):
        def dfs(n, target, path):
            if n is None:
                return None
            path.append(n.idx)
            if n.idx == target:
                return list(path)
            if n.left:
                r = dfs(n.left, target, list(path))
                if r:
                    return r
            if n.right:
                r = dfs(n.right, target, list(path))
                if r:
                    return r
            return None
        return dfs(node, target_idx, [])

    def tree_distance(n1, n2):
        path1 = path_to_root(root, n1)
        path2 = path_to_root(root, n2)
        i = 0
        while i < min(len(path1), len(path2)) and path1[i] == path2[i]:
            i += 1
        return len(path1) + len(path2) - 2 * i

    leaves = [a, b, c, d]
    dm = np.zeros((4, 4))
    for i in range(4):
        for j in range(i + 1, 4):
            d_val = tree_distance(leaves[i], leaves[j])
            dm[i][j] = dm[j][i] = d_val

    pairs = [(0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2)]
    best = 0
    best_sum = float('inf')
    for idx, (i1, j1, i2, j2) in enumerate(pairs):
        s = dm[i1][j1] + dm[i2][j2]
        if s < best_sum:
            best_sum = s
            best = idx
    return best
