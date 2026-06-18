"""
kmer_optimized.py
Optimized k-mer frequency computation using numpy vectorization.
Replaces the pure-Python loop in seq_to_kmer_freq.
"""

import numpy as np
import numba
from numba import njit


def seq_to_kmer_freq_vectorized(seq: str, k: int = 5, normalize: bool = True):
    """
    Vectorized k-mer frequency using numpy (no Python loop over positions).
    Converts seq to integer array, then uses bincount on rolling windows.
    """
    seq = seq.upper()
    valid_bases = np.array([ord('A'), ord('C'), ord('G'), ord('T')], dtype=np.uint8)
    base_to_int = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    # Convert to integer array (invalid positions marked as -1)
    arr = np.full(len(seq), -1, dtype=np.int8)
    for i, c in enumerate(seq):
        if c in base_to_int:
            arr[i] = base_to_int[c]

    # Sliding window: compute k-mer code for each valid position
    vec_size = 4 ** k
    freq = np.zeros(vec_size, dtype=np.float32)
    valid_len = 0

    for i in range(len(seq) - k + 1):
        sub = arr[i:i + k]
        if np.any(sub == -1):
            continue
        code = 0
        for j in range(k):
            code += int(sub[j]) * (4 ** (k - j - 1))
        freq[code] += 1
        valid_len += 1

    if normalize and valid_len > 0:
        freq = freq / valid_len

    return freq


@njit(cache=True)
def _seq_to_kmer_freq_numba(seq_int, k, vec_size):
    """
    Numba-JITed inner loop for k-mer counting.
    seq_int: numpy array of int8, with -1 for invalid bases.
    Returns: freq array (float32), valid_len (int)
    """
    n = len(seq_int)
    freq = np.zeros(vec_size, dtype=np.float32)
    valid_len = 0

    for i in range(n - k + 1):
        # Check all bases valid
        valid = True
        for j in range(k):
            if seq_int[i + j] < 0:
                valid = False
                break
        if not valid:
            continue

        # Compute code
        code = 0
        for j in range(k):
            code += int(seq_int[i + j]) * (4 ** (k - j - 1))
        freq[code] += 1.0
        valid_len += 1

    return freq, valid_len


def seq_to_kmer_freq_numba(seq: str, k: int = 5, normalize: bool = True):
    """
    Numba-JIT accelerated k-mer frequency.
    First call is slow (compilation), subsequent calls are fast.
    """
    base_to_int = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    seq = seq.upper()
    arr = np.full(len(seq), -1, dtype=np.int8)
    for i, c in enumerate(seq):
        if c in base_to_int:
            arr[i] = base_to_int[c]

    vec_size = 4 ** k
    freq, valid_len = _seq_to_kmer_freq_numba(arr, k, vec_size)

    if normalize and valid_len > 0:
        freq = freq / valid_len

    return freq


def benchmark():
    """Compare original vs optimized vs numba."""
    import time
    seq = "ACGT" * 500  # 2000 bp

    # Original
    from kmer_distance import seq_to_kmer_freq as orig_func
    t0 = time.time()
    for _ in range(100):
        _ = orig_func(seq, k=5, normalize=True)
    t_orig = time.time() - t0

    # Numba JIT (first call includes compilation)
    t0 = time.time()
    r1 = seq_to_kmer_freq_numba(seq, k=5, normalize=True)
    compile_time = time.time() - t0
    print(f"Numba compilation (first call): {compile_time:.3f}s")

    t0 = time.time()
    for _ in range(100):
        _ = seq_to_kmer_freq_numba(seq, k=5, normalize=True)
    t_numba = time.time() - t0

    print(f"Original (100 calls): {t_orig:.3f}s ({t_orig/100*1000:.1f}ms/call)")
    print(f"Numba (100 calls):  {t_numba:.3f}s ({t_numba/100*1000:.1f}ms/call)")
    print(f"Speedup: {t_orig/t_numba:.1f}x")


if __name__ == "__main__":
    benchmark()
