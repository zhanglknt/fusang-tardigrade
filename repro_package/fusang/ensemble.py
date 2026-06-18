#!/usr/bin/env python3
"""
fusang-tardigrade 多 k/gap 集成模块
=====================================
通过拼接多个 (k, gap) 配置的 k-mer 频率向量，
构建「多尺度距离矩阵」，实现鲁棒的免比对建树。

两种集成策略：
  A. Feature Concatenation: 拼接 k-mer 向量 → 单棵 NJ 树
  B. Distance Averaging: 各配置独立距离矩阵 → 加权平均 → NJ 树

参考: Fusang: Tardigrade Edition
"""

import numpy as np
from scipy.spatial.distance import pdist, squareform
from typing import List, Dict, Tuple, Optional
import time

# 本地模块
from kmer_distance import (
    make_gap_pattern,
    seq_to_kmer_freq,
    compute_kmer_distance_matrix,
    kmer_distance,
)


# ============================================================
# 策略 A: 特征拼接（Feature Concatenation）
# ============================================================

def ensemble_feature_concat(
    sequences: List[str],
    taxon_names: List[str],
    configs: List[Tuple[int, str]],  # [(k, gap_style), ...]
    normalize_per_config: bool = True,
    n_threads: int = 4,
) -> np.ndarray:
    """
    特征拼接集成：
    1. 对每个 (k, gap) 配置，独立计算 k-mer 频率向量
    2. 将所有配置的向量水平拼接
    3. 对拼接向量做 L2 归一化
    4. 计算 cosine 距离矩阵

    参数：
        sequences: 序列列表
        taxon_names: 分类单元名
        configs: [(k, gap_style), ...] 如 [(4,"gap1"), (5,"gap2"), (6,"gap2")]
        normalize_per_config: 是否先对各配置独立归一化再拼接
        n_threads: 并行线程数

    返回：
        dist_matrix: (n, n) 距离矩阵
    """
    n = len(sequences)
    print(f"\n[Ensemble::FeatureConcat] {len(configs)} configs: {configs}")

    all_features = []  # List[np.ndarray], each shape (n, vocab_size)

    for k, gap_style in configs:
        pattern = make_gap_pattern(k, gap_style)
        pk = len(pattern)
        vocab_size = 4 ** pk
        t0 = time.time()

        # 逐序列计算 k-mer 频率
        feat = np.zeros((n, vocab_size), dtype=np.float32)
        for i, seq in enumerate(sequences):
            freq = seq_to_kmer_freq(seq, k=k, normalize=normalize_per_config,
                                    gap_pattern=pattern)
            feat[i, :] = freq

        elapsed = time.time() - t0
        all_features.append(feat)
        print(f"  k={k}, {gap_style}: vocab={vocab_size}, time={elapsed:.1f}s")

    # 拼接所有特征 → shape (n, total_vocab)
    concat = np.hstack(all_features)
    total_dim = concat.shape[1]
    print(f"  Concatenated dim: {total_dim}")

    # L2 归一化拼接向量（保证 cosine 距离有效）
    norms = np.linalg.norm(concat, axis=1, keepdims=True)
    norms[norms == 0] = 1.0  # 防止除零
    concat_normalized = concat / norms

    # 计算 cosine 距离（用 pdist 加速，cosine 距离 = 1 - u·v）
    # scipy pdist 'cosine' metric 就是 1 - u·v/(|u||v|)
    dist_condensed = pdist(concat_normalized, metric='cosine')
    dist_matrix = squareform(dist_condensed)

    return dist_matrix


# ============================================================
# 策略 B: 距离矩阵加权平均（Distance Averaging）
# ============================================================

def ensemble_distance_avg(
    sequences: List[str],
    taxon_names: List[str],
    configs: List[Tuple[int, str]],
    weights: Optional[List[float]] = None,
    metric: str = 'cosine',
    n_threads: int = 4,
) -> np.ndarray:
    """
    距离矩阵加权平均集成：
    1. 对每个 (k, gap) 配置，独立计算距离矩阵
    2. 加权平均所有距离矩阵
    3. 返回平均距离矩阵

    参数：
        weights: 各配置的权重，None 则等权
    """
    n = len(sequences)
    n_configs = len(configs)
    if weights is None:
        weights = [1.0 / n_configs] * n_configs

    print(f"\n[Ensemble::DistAvg] {n_configs} configs, weights={[f'{w:.3f}' for w in weights]}")

    avg_dist = np.zeros((n, n), dtype=np.float64)

    for idx, (k, gap_style) in enumerate(configs):
        pattern = make_gap_pattern(k, gap_style)
        t0 = time.time()

        dist_mat = compute_kmer_distance_matrix(
            sequences, taxon_names,
            k=k, metric=metric, n_threads=n_threads,
            gap_pattern=pattern,
        )
        elapsed = time.time() - t0
        avg_dist += weights[idx] * dist_mat
        print(f"  k={k}, {gap_style}: time={elapsed:.1f}s")

    return avg_dist


# ============================================================
# 快速 NJ 建树（从距离矩阵）
# ============================================================

def nj_from_distance(dist_matrix: np.ndarray,
                     names: List[str]) -> str:
    """
    从距离矩阵用 NJ 算法建树，返回 Newick 字符串。
    使用 scipy 的 hierarchical clustering 实现。
    简单实现：http://evolution.gs.washington.edu/phylip/newicktree.html
    """
    from scipy.cluster.hierarchy import (
        linkage, to_tree, ClusterNode, _convert_to_newick
    )
    # 注意: scipy 的 NJ 实际上用的是 UPGMA 的 linkage
    # 我们需要真正的 NJ。让我用 biopython（如果可用）或手写简化版。
    # 这里使用一个已知正确的 NJ 实现。
    pass


# ============================================================
# 工具函数
# ============================================================

def read_fasta(path: str) -> Tuple[List[str], List[str]]:
    """读取 FASTA，返回 (sequences, taxon_names)"""
    names = []
    seqs = []
    current_name = None
    current_seq = []
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_name is not None:
                    names.append(current_name)
                    seqs.append(''.join(current_seq))
                current_name = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_name is not None:
            names.append(current_name)
            seqs.append(''.join(current_seq))
    return seqs, names


def read_newick(path: str) -> str:
    with open(path, 'r') as f:
        return f.read().strip()


# ============================================================
# NJ 实现 (Neighbour-Joining, 从距离矩阵)
# ============================================================

def nj_tree(dist_matrix: np.ndarray, names: List[str]) -> str:
    """
    标准 Neighbour-Joining 算法 (Saitou & Nei 1987)。
    输入：(n,n) 距离矩阵 + n 个名称
    输出：Newick 字符串
    """
    n = dist_matrix.shape[0]
    if n == 1:
        return f"{names[0]};"
    if n == 2:
        d = dist_matrix[0, 1] / 2.0
        return f"({names[0]}:{d},{names[1]}:{d});"

    D = dist_matrix.copy().astype(np.float64)
    # 活跃节点 ID 列表
    active_ids = list(range(n))
    children = {}  # node_id -> str (leaf) or (child_a, child_b, edge_a, edge_b)
    for i in range(n):
        children[i] = names[i]

    next_id = n

    while len(active_ids) > 2:
        m = len(active_ids)
        # 节点发散度 r[i]
        r = D.sum(axis=1) / (m - 2)
        r_col = r.reshape(-1, 1)

        # Q 矩阵: Q(i,j) = D(i,j) - r[i] - r[j]
        Q = D - r_col - r_col.T

        # 找最小 Q（跳过对角线）
        np.fill_diagonal(Q, np.inf)
        min_flat = np.argmin(Q)
        i, j = min_flat // m, min_flat % m

        # 分支长度
        d_ij = D[i, j]
        edge_i = max(0.0, 0.5 * d_ij + 0.5 * (r[i] - r[j]))
        edge_j = max(0.0, d_ij - edge_i)
        if edge_j < 0:
            edge_j = d_ij / 2.0
            edge_i = d_ij / 2.0

        # 新节点
        new_id = next_id
        next_id += 1
        children[new_id] = (active_ids[i], active_ids[j], edge_i, edge_j)

        # 新节点到其他节点的距离
        new_dists = 0.5 * (D[i, :] + D[j, :] - d_ij)

        # 重建距离矩阵（去掉 i, j，添加新节点）
        keep = [k for k in range(m) if k != i and k != j]
        new_m = m - 1
        D_new = np.zeros((new_m, new_m), dtype=np.float64)

        # 保留的节点间距离
        for a in range(new_m - 1):
            for b in range(new_m - 1):
                D_new[a, b] = D[keep[a], keep[b]]

        # 新节点（最后一行/列）
        for a in range(new_m - 1):
            D_new[a, new_m - 1] = new_dists[keep[a]]
            D_new[new_m - 1, a] = new_dists[keep[a]]

        D = D_new
        new_active = [active_ids[k] for k in keep] + [new_id]
        active_ids = new_active

    # 最后两个节点
    d = D[0, 1] / 2.0
    root_id = next_id
    children[root_id] = (active_ids[0], active_ids[1], max(0.0, d), max(0.0, d))

    return _to_newick(root_id, children) + ";"


def _to_newick(node_id: int, children: dict) -> str:
    """递归生成 Newick 字符串"""
    val = children[node_id]
    if isinstance(val, str):
        return val  # leaf
    child_a, child_b, edge_a, edge_b = val
    left = _to_newick(child_a, children)
    right = _to_newick(child_b, children)
    return f"({left}:{edge_a:.6f},{right}:{edge_b:.6f})"


# ============================================================
# nRF (normalized Robinson-Foulds distance)
# ============================================================

def calc_nrf(tree1_newick: str, tree2_newick: str) -> float:
    """计算两棵新格式树之间的 normalized RF 距离。"""
    from io import StringIO
    try:
        # 尝试使用 dendropy（如果有）
        import dendropy
        t1 = dendropy.Tree.get(data=tree1_newick, schema="newick")
        t2 = dendropy.Tree.get(data=tree2_newick, schema="newick")
        t1.encode_bipartitions()
        t2.encode_bipartitions()
        rf = dendropy.calculate.treecompare.symmetric_difference(t1, t2)
        n = len(t1.taxon_namespace)
        max_rf = 2 * (n - 3) if n > 3 else 2 * (n - 2)
        if max_rf > 0:
            return rf / max_rf
        return 0.0
    except ImportError:
        pass

    # 备选：使用 ete3
    try:
        from ete3 import Tree
        t1 = Tree(tree1_newick)
        t2 = Tree(tree2_newick)
        rf = t1.robinson_foulds(t2)
        rf_abs = rf[0]
        n = len(list(t1.iter_leaves()))
        max_rf = 2 * (n - 3) if n > 3 else 2 * (n - 2)
        if max_rf > 0:
            return rf_abs / max_rf
        return 0.0
    except ImportError:
        pass

    # 终极备选：用 calc_nrf_simple
    pass


# ============================================================
# 主程序（诊断用）
# ============================================================

if __name__ == "__main__":
    import argparse
    import os
    import sys

    # 添加父目录到 path 以便导入 calc_nrf_simple
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from calc_nrf_simple import calc_nrf, get_bipartitions_from_newick

    def calc_nrf_simple(tree1_nwk: str, tree2_nwk: str) -> float:
        """计算两棵 Newick 字符串之间的 nRF。"""
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.nwk', delete=False) as f1, \
             tempfile.NamedTemporaryFile(mode='w', suffix='.nwk', delete=False) as f2:
            f1.write(tree1_nwk)
            f1.close()
            f2.write(tree2_nwk)
            f2.close()
            result = calc_nrf(f1.name, f2.name)
            os.unlink(f1.name)
            os.unlink(f2.name)
        return result

    parser = argparse.ArgumentParser(
        description="Fusang: Tardigrade Edition - Ensemble Benchmark"
    )
    parser.add_argument("--fasta", default="test_n200.fasta",
                        help="FASTA 文件")
    parser.add_argument("--true-tree", default="test_n200_true.nwk",
                        help="标准树")
    parser.add_argument("--configs", default="4,gap1,5,gap2,6,gap2,7,gap2",
                        help="配置列表: k1,gap1,k2,gap2,...")
    parser.add_argument("--n-threads", type=int, default=4)
    args = parser.parse_args()

    DATA_DIR = os.path.dirname(os.path.abspath(__file__))
    fasta_path = os.path.join(DATA_DIR, args.fasta)
    true_path = os.path.join(DATA_DIR, args.true_tree)

    print(f"Loading: {fasta_path}")
    sequences, taxon_names = read_fasta(fasta_path)
    n = len(sequences)
    print(f"  {n} sequences")

    true_newick = read_newick(true_path)
    print(f"  True tree loaded")

    # 解析配置
    parts = args.configs.split(",")
    configs = []
    for i in range(0, len(parts), 2):
        k = int(parts[i])
        gap = parts[i + 1]
        configs.append((k, gap))

    # ============================================================
    # 1. Baseline: 各单配置
    # ============================================================
    print("\n" + "=" * 60)
    print("BASELINE: Single config performance")
    print("=" * 60)

    baseline_results = {}
    for k, gap_style in configs:
        pattern = make_gap_pattern(k, gap_style)
        t0 = time.time()
        dist_mat = compute_kmer_distance_matrix(
            sequences, taxon_names, k=k, metric='cosine',
            n_threads=args.n_threads, gap_pattern=pattern
        )
        tree = nj_tree(dist_mat, taxon_names)
        nrf = calc_nrf_simple(tree, true_newick)
        elapsed = time.time() - t0
        print(f"  k={k},{gap_style}: nRF={nrf:.6f}, time={elapsed:.1f}s")
        baseline_results[(k, gap_style)] = {"nrf": nrf, "tree": tree}

    best_nrf = min(r["nrf"] for r in baseline_results.values())
    print(f"  Best single: nRF={best_nrf:.6f}")

    # ============================================================
    # 2. Ensemble A: Feature Concatenation
    # ============================================================
    print("\n" + "=" * 60)
    print("ENSEMBLE A: Feature Concatenation")
    print("=" * 60)

    t0 = time.time()
    concat_dist = ensemble_feature_concat(
        sequences, taxon_names, configs,
        normalize_per_config=True, n_threads=args.n_threads
    )
    concat_tree = nj_tree(concat_dist, taxon_names)
    concat_nrf = calc_nrf_simple(concat_tree, true_newick)
    elapsed = time.time() - t0
    print(f"  nRF={concat_nrf:.6f}, time={elapsed:.1f}s")

    # ============================================================
    # 3. Ensemble B: Distance Averaging
    # ============================================================
    print("\n" + "=" * 60)
    print("ENSEMBLE B: Distance Averaging")
    print("=" * 60)

    t0 = time.time()
    avg_dist = ensemble_distance_avg(
        sequences, taxon_names, configs,
        weights=None, n_threads=args.n_threads
    )
    avg_tree = nj_tree(avg_dist, taxon_names)
    avg_nrf = calc_nrf_simple(avg_tree, true_newick)
    elapsed = time.time() - t0
    print(f"  nRF={avg_nrf:.6f}, time={elapsed:.1f}s")

    # ============================================================
    # 4. Summary
    # ============================================================
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  {'Method':<30s} {'nRF':>10s}")
    print(f"  {'-'*40}")
    for (k, gap), r in baseline_results.items():
        print(f"  Single k={k},{gap:<5s}  {r['nrf']:10.6f}")
    print(f"  {'-'*40}")
    print(f"  Ensemble:FeatureConcat   {concat_nrf:10.6f}")
    print(f"  Ensemble:DistAvg         {avg_nrf:10.6f}")
    print(f"  {'-'*40}")

    best_ensemble = min(concat_nrf, avg_nrf)
    improvement = best_nrf - best_ensemble
    if improvement > 0:
        print(f"  Ensemble IMPROVES over best single by {improvement:.6f}")
    else:
        print(f"  Ensemble does NOT improve over best single (diff={improvement:.6f})")
