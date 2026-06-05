"""
k-mer 距离模块
用于 Fusang v2 分级建树架构

功能：
1. 对每条序列计算 k-mer 频率向量（Numba JIT 加速）
2. 计算序列间余弦距离矩阵
3. 支持 FASTA（未比对）和 MSA（已比对）两种输入
4. 支持并行加速（ProcessPoolExecutor）
"""

import numpy as np
from typing import List, Dict, Tuple, Optional
from functools import lru_cache
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import numba
from numba import njit, types


# =========================================================
# 标准 gapped k-mer pattern
# =========================================================

def make_gap_pattern(k: int, style: str = 'gap1') -> tuple:
    """
    生成标准 gapped k-mer pattern

    参数：
        k: 有效碱基数（pattern 长度）
        style: 'gap1' — 每 2 个有效碱基后跳 1 位
               'gap2' — 在中间跳 2 位

    返回：
        tuple of int, 相对位置列表

    示例 (k=5):
        gap1 → (0,1,2,4,5)    : 每 3 个位置取 2 个，窗口 = 2*ceil(5/2)+1=7
        gap2 → (0,1,2,5,6)    : 前半 + 后半各取部分，中间跳
    """
    if style == 'gap1':
        # 交替取-取-跳: 位置 0,1, skip 2, 3,4, skip 5, ...
        # 每组取2个跳1个, 需要的窗口长度 = k + (k-1)//2
        positions = []
        pos = 0
        count = 0
        while count < k:
            positions.append(pos)
            pos += 1
            count += 1
            if count < k:
                positions.append(pos)
                pos += 1
                count += 1
            if count < k:
                pos += 1  # skip one
        return tuple(positions)

    elif style == 'gap2':
        # 前 ceil(k/2) 个连续 + 跳 2 + 后 floor(k/2) 个连续
        first_half = k - k // 2
        second_half = k // 2
        return tuple(list(range(first_half)) + list(range(first_half + 2, first_half + 2 + second_half)))

    elif style == 'gap3':
        # 前 ceil(k/2) 个连续 + 跳 3 + 后 floor(k/2) 个连续
        first_half = k - k // 2
        second_half = k // 2
        return tuple(list(range(first_half)) + list(range(first_half + 3, first_half + 3 + second_half)))

    elif style == 'gap4':
        # 前 ceil(k/2) 个连续 + 跳 4 + 后 floor(k/2) 个连续
        first_half = k - k // 2
        second_half = k // 2
        return tuple(list(range(first_half)) + list(range(first_half + 4, first_half + 4 + second_half)))

    else:
        raise ValueError(f"Unknown gap style: {style}")


# =========================================================
# Numba-JIT 加速的 k-mer 频率计算（核心热点）
# =========================================================

@njit(cache=True, fastmath=True)
def _kmer_freq_numba(seq_int, k, vec_size):
    """
    Numba-JIT 编译的内部函数：计算单条序列的 k-mer 频率向量。
    seq_int: int8 数组，A=0,C=1,G=2,T=3，无效位置=-1
    返回：(freq 数组, valid_len)
    """
    n = len(seq_int)
    freq = np.zeros(vec_size, dtype=np.float32)
    valid_len = 0

    for i in range(n - k + 1):
        # 检查窗口内是否全有效
        valid = True
        for j in range(k):
            if seq_int[i + j] < 0:
                valid = False
                break
        if not valid:
            continue

        # 计算 k-mer code（base-4 编码）
        code = 0
        for j in range(k):
            code += np.int32(seq_int[i + j]) * np.int32(4 ** (k - j - 1))
        freq[code] += 1.0
        valid_len += 1

    return freq, valid_len


@njit(cache=True, fastmath=True)
def _kmer_freq_gapped_numba(seq_int, pattern, vec_size):
    """
    Numba-JIT 编译的内部函数：计算 gapped k-mer 频率向量。
    seq_int: int8 数组，A=0,C=1,G=2,T=3，无效位置=-1
    pattern: int 数组，gapped 位置列表（如 [0,1,3,4]）
    vec_size: 向量大小 = 4^len(pattern)
    返回：(freq 数组, valid_len)
    """
    n = len(seq_int)
    pk = len(pattern)
    max_pos = pattern[pk - 1]
    freq = np.zeros(vec_size, dtype=np.float32)
    valid_len = 0

    for i in range(n - max_pos):
        # 检查 pattern 指定位置是否全有效
        valid = True
        for j in range(pk):
            if seq_int[i + pattern[j]] < 0:
                valid = False
                break
        if not valid:
            continue

        # 计算 gapped k-mer code（base-4 编码）
        code = 0
        for j in range(pk):
            code += np.int32(seq_int[i + pattern[j]]) * np.int32(4 ** (pk - j - 1))
        freq[code] += 1.0
        valid_len += 1

    return freq, valid_len


def _seq_to_int_array(seq_upper):
    """把序列字符串转成 int8 数组（A=0,C=1,G=2,T=3，无效=-1）"""
    arr = np.full(len(seq_upper), -1, dtype=np.int8)
    for i, c in enumerate(seq_upper):
        if c == 'A':
            arr[i] = 0
        elif c == 'C':
            arr[i] = 1
        elif c == 'G':
            arr[i] = 2
        elif c == 'T':
            arr[i] = 3
    return arr


def seq_to_kmer_freq(seq: str, k: int = 5, normalize: bool = True, gap_pattern: Optional[tuple] = None):
    """
    将一条序列转换为 k-mer 频率向量（4^k 维）
    使用 Numba JIT 加速，首次调用会编译（~1.5s），后续调用快 10-20x。

    参数：
        seq: 序列字符串（ACGT，忽略大小写，'-' 和 'N' 跳过）
        k: k-mer 长度（默认 5，4^5=1024 维）。支持 k=5/6/7。
            当 gap_pattern 提供时，k 表示 pattern 中的碱基数（非窗口长度）。
        normalize: 是否归一化为概率分布（默认 True）
        gap_pattern: 可选的 gapped k-mer pattern（tuple of int，相对位置列表）。
            例如 (0,1,3,4) 表示跳过位置 2 的 4-mer。
            None 时使用默认连续模式 (0,1,...,k-1)。

    返回：
        numpy array，shape=(4^k,)，dtype=float32
    """
    seq_upper = seq.upper()
    seq_int = _seq_to_int_array(seq_upper)

    if gap_pattern is None:
        # 默认连续 k-mer
        vec_size = 4 ** k
        freq, valid_len = _kmer_freq_numba(seq_int, k, vec_size)
    else:
        # gapped k-mer
        pattern_arr = np.array(gap_pattern, dtype=np.int64)
        pk = len(gap_pattern)
        vec_size = 4 ** pk
        freq, valid_len = _kmer_freq_gapped_numba(seq_int, pattern_arr, vec_size)

    if normalize and valid_len > 0:
        freq = freq / valid_len

    return freq


# =========================================================
# 距离计算
# ==========================================================

def kmer_distance(vec1, vec2, metric: str = 'cosine'):
    """
    计算两个 k-mer 频率向量之间的距离

    参数：
        vec1, vec2: k-mer 频率向量（numpy array）
        metric: 距离度量方式
            'cosine' : 余弦距离 (1 - cos_sim)，范围 [0, 2]
            'euclidean': 欧氏距离
            'jaccard' : Jaccard 距离

    返回：
        float: 距离值，0 表示完全相同
    """
    if metric == 'cosine':
        dot = np.dot(vec1, vec2)
        norm1 = np.linalg.norm(vec1)
        norm2 = np.linalg.norm(vec2)
        if norm1 == 0 or norm2 == 0:
            return 1.0
        cos_sim = dot / (norm1 * norm2)
        return max(0.0, min(1.0, 1.0 - cos_sim))

    elif metric == 'euclidean':
        return np.linalg.norm(vec1 - vec2)

    elif metric == 'jaccard':
        b1 = (vec1 > 0).astype(np.float32)
        b2 = (vec2 > 0).astype(np.float32)
        intersection = np.dot(b1, b2)
        union = np.dot(b1, b1) + np.dot(b2, b2) - intersection
        if union == 0:
            return 0.0
        return 1.0 - intersection / union

    else:
        raise ValueError(f"Unsupported metric: {metric}")


def compute_kmer_distance_matrix(sequences: List[str],
                                 taxon_names: List[str],
                                 k: int = 5,
                                 metric: str = 'cosine',
                                 n_threads: int = 4,
                                 gap_pattern: Optional[tuple] = None):
    """
    计算所有序列对之间的 k-mer 距离矩阵

    参数：
        sequences: 序列字符串列表（FASTA 或 MSA）
        taxon_names: 分类单元名称列表（与 sequences 一一对应）
        k: k-mer 长度（默认 5，支持 5/6/7）
        metric: 距离度量方式（默认 'cosine'）
        n_threads: 并行线程数（默认 4）
        gap_pattern: 可选 gapped k-mer pattern（tuple of int）

    返回：
        numpy array, shape=(n, n), dtype=float32
        distance_matrix[i][j] = distance(taxon_i, taxon_j)
    """
    n = len(sequences)
    if n == 0:
        return np.array([], dtype=np.float32)

    pattern_desc = f"gap_pattern={gap_pattern}" if gap_pattern else "contiguous"
    print(f"[kmer_distance] 计算 {n}×{n} 距离矩阵，k={k}，{pattern_desc}，度量={metric}，线程数={n_threads}")

    # Step 1: 并行计算所有序列的 k-mer 频率向量
    print(f"[kmer_distance] Step 1: 计算 {n} 条序列的 k-mer 频率向量...")

    if n_threads > 1 and n > 1:
        with ProcessPoolExecutor(max_workers=min(n_threads, n)) as executor:
            futures = {
                executor.submit(seq_to_kmer_freq, seq, k, True, gap_pattern): i
                for i, seq in enumerate(sequences)
            }
            freq_vectors = [None] * n
            done = 0
            for future in as_completed(futures):
                idx = futures[future]
                freq_vectors[idx] = future.result()
                done += 1
                if done % max(1, n // 10) == 0:
                    print(f"[kmer_distance]   频率向量计算进度: {done}/{n}")
    else:
        freq_vectors = [seq_to_kmer_freq(seq, k, True, gap_pattern) for seq in sequences]

    freq_vectors = np.array(freq_vectors, dtype=np.float32)  # (n, 4^k)

    # Step 2: 计算距离矩阵（利用矩阵运算加速）
    print(f"[kmer_distance] Step 2: 计算距离矩阵...")

    if metric == 'cosine':
        norms = np.linalg.norm(freq_vectors, axis=1, keepdims=True)
        norms[norms == 0] = 1.0
        normalized = freq_vectors / norms

        cos_sim = normalized @ normalized.T  # (n, n)

        dist_matrix = np.clip(1.0 - cos_sim, 0.0, 1.0)

        np.fill_diagonal(dist_matrix, 0.0)

    else:
        dist_matrix = np.zeros((n, n), dtype=np.float32)
        for i in range(n):
            for j in range(i + 1, n):
                d = kmer_distance(freq_vectors[i], freq_vectors[j], metric)
                dist_matrix[i][j] = d
                dist_matrix[j][i] = d
            if i % max(1, n // 10) == 0:
                print(f"[kmer_distance]   距离矩阵计算进度: {i+1}/{n}")

    print(f"[kmer_distance] 距离矩阵计算完成。形状: {dist_matrix.shape}")
    return dist_matrix


# =========================================================
# FASTA / 距离矩阵文件 I/O
# =========================================================

def read_fasta(file_path: str):
    """
    读取 FASTA 文件

    参数：
        file_path: FASTA 文件路径

    返回：
        (taxon_names, sequences)
    """
    taxon_names = []
    sequences = []
    current_seq = []

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_seq:
                    sequences.append(''.join(current_seq))
                    current_seq = []
                name = line[1:].split()[0]
                taxon_names.append(name)
            else:
                current_seq.append(line.upper())

        if current_seq:
            sequences.append(''.join(current_seq))

    return taxon_names, sequences


def write_distance_matrix(dist_matrix, taxon_names: List[str],
                         output_path: str):
    """
    将距离矩阵写入文件（PHYLIP 格式）
    """
    n = len(taxon_names)
    with open(output_path, 'w') as f:
        f.write(f"{n}\n")
        for i in range(n):
            name_padded = f"{taxon_names[i]:<10}"[:10]
            dist_str = ' '.join(f"{dist_matrix[i][j]:.6f}" for j in range(n))
            f.write(f"{name_padded} {dist_str}\n")


if __name__ == '__main__':
    print("=== kmer_distance 模块测试 ===")

    seq1 = "ACGT" * 50
    seq2 = "ACGT" * 50
    seq3 = "TGCA" * 50

    test_seqs = [seq1, seq2, seq3]
    test_names = ["seq1", "seq2", "seq3"]

    dist_mat = compute_kmer_distance_matrix(test_seqs, test_names, k=3, n_threads=1)
    print(f"测试距离矩阵:\n{dist_mat}")
    print(f"seq1-seq2 距离（应接近 0）: {dist_mat[0][1]:.4f}")
    print(f"seq1-seq3 距离（应较大）: {dist_mat[0][2]:.4f}")
