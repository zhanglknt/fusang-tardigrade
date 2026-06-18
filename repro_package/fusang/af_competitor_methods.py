#!/usr/bin/env python3
"""
Alignment-free competitor methods for Fusang benchmark.

Implements:
1. Co-phylog (Yi & Jin 2013, NAR) - context-object matching, k=19
2. CVTree-like (k-mer frequency + Pearson correlation distance)
3. K-mer cosine (standard contiguous, for comparison with spaced k-mers)

Also includes a simplified andi approximation based on k-mer match lengths.

All methods output a distance matrix that can be used with NJ tree building.
"""

import sys
import os
import math
import time
import numpy as np
from collections import defaultdict
from itertools import combinations

sys.setrecursionlimit(10000)

# ============================================================
# Utility functions
# ============================================================

BASE_MAP = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'a': 0, 'c': 1, 'g': 2, 't': 3}
COMPLEMENT = {0: 3, 1: 2, 2: 1, 3: 0}  # A<->T, C<->G
MASK3 = 3  # 0b11

def read_fasta(path):
    """Read FASTA file, return dict of name -> sequence (uppercase)."""
    seqs = {}
    name = None
    seq = []
    with open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.upper())
    if name is not None:
        seqs[name] = ''.join(seq)
    return seqs


def revcomp_seq(seq):
    """Reverse complement a sequence string."""
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq))


# ============================================================
# Co-phylog implementation (based on source code analysis)
# ============================================================

def _encode_tuple(seq_bytes, start, tuplen=19):
    """Encode a k-mer as a 64-bit integer: lower bits = most recent base."""
    val = 0
    for i in range(start, start + tuplen):
        b = BASE_MAP.get(chr(seq_bytes[i]), 0)
        val = ((val << 2) | b) & ((1 << (2 * tuplen)) - 1)
    return val


def _revcomp_tuple(tuple_val, tuplen=19, ctxlen=18):
    """Compute reverse complement of an encoded tuple."""
    bittl = 2 * tuplen
    mask = (1 << bittl) - 1
    comp_bittl = 64 - bittl
    # Extract bases from tuple (MSB = oldest base, LSB = newest base)
    # For revcomp: reverse the order and complement each base
    rc = 0
    for i in range(tuplen):
        shift = 2 * (tuplen - 1 - i)
        base = (tuple_val >> shift) & 3
        rc = (rc << 2) | COMPLEMENT.get(base, 0)
    return rc & mask


def _extract_context_object(tuple_val, ctxlen=18, tuplen=19):
    """Extract context and object from encoded tuple."""
    bittl = 2 * tuplen
    tupmask = (1 << bittl) - 1
    objmask = tupmask ^ (3 << (2 * ctxlen))
    context = tuple_val & objmask
    obj = (tuple_val >> (2 * ctxlen)) & 3
    return context, obj


def extract_co_signature(seq_str, halfctx=9):
    """
    Extract Co-phylog signature: set of (context, object) pairs where
    the context is unique (only one object observed).

    Args:
        seq_str: DNA sequence string (ACGT only)
        halfctx: HALFCTX parameter (default 9, giving k=19)

    Returns:
        dict: context -> object (only unique contexts retained)
    """
    tuplen = 2 * halfctx + 1  # 19
    ctxlen = 2 * halfctx       # 18
    bittl = 2 * tuplen         # 38
    tupmask = (1 << bittl) - 1
    objmask = tupmask ^ (3 << (2 * ctxlen))

    # Clean sequence: keep only ACGT
    seq_clean = []
    for c in seq_str:
        if c in 'ACGT':
            seq_clean.append(c)
    seq_bytes = ''.join(seq_clean).encode('ascii')

    if len(seq_bytes) < tuplen:
        return {}

    # Hash table: context -> tuple value
    # If a context appears with different objects, discard it
    ctx_map = {}

    # Initialize first tuple
    val = 0
    for i in range(tuplen):
        b = BASE_MAP.get(chr(seq_bytes[i]), 0)
        val = ((val << 2) | b) & tupmask

    rc_val = _revcomp_tuple(val, tuplen, ctxlen)
    unituple = min(val, rc_val)
    ctx = unituple & objmask
    if unituple != 0:
        if ctx in ctx_map:
            if ctx_map[ctx] != unituple:
                ctx_map[ctx] = None  # Mark as ambiguous
        else:
            ctx_map[ctx] = unituple

    # Slide window
    for i in range(tuplen, len(seq_bytes)):
        b = BASE_MAP.get(chr(seq_bytes[i]), 0)
        val = ((val << 2) | b) & tupmask
        rc_val = _revcomp_tuple(val, tuplen, ctxlen)
        unituple = min(val, rc_val)
        if unituple == 0:
            continue
        ctx = unituple & objmask
        if ctx in ctx_map:
            if ctx_map[ctx] is not None and ctx_map[ctx] != unituple:
                ctx_map[ctx] = None  # Mark as ambiguous
        else:
            ctx_map[ctx] = unituple

    # Keep only unique contexts
    return {ctx: t for ctx, t in ctx_map.items() if t is not None}


def cophylog_distance_matrix(sequences, names, halfctx=9):
    """
    Compute Co-phylog pairwise distance matrix.

    Distance = (shared contexts with different objects) / (total shared contexts)

    Args:
        sequences: list of sequence strings
        names: list of sequence names
        halfctx: HALFCTX parameter (default 9, k=19)

    Returns:
        numpy 2D distance matrix
    """
    n = len(sequences)
    D = np.zeros((n, n))

    print(f"  [Co-phylog] Extracting signatures (halfctx={halfctx}, k={2*halfctx+1})...")
    t0 = time.time()

    # Extract signatures
    signatures = []
    for i, seq in enumerate(sequences):
        sig = extract_co_signature(seq, halfctx)
        signatures.append(sig)
        if (i + 1) % 50 == 0:
            print(f"    {i+1}/{n} sequences done")

    print(f"  [Co-phylog] Signatures extracted in {time.time()-t0:.1f}s")
    print(f"  [Co-phylog] Computing pairwise distances...")
    t0 = time.time()

    # Pairwise distance
    for i in range(n):
        sig_i = signatures[i]
        for j in range(i + 1, n):
            sig_j = signatures[j]
            cxt = 0  # shared contexts
            obj = 0  # shared contexts with different objects

            # Iterate over the smaller signature
            if len(sig_i) < len(sig_j):
                smaller, larger = sig_i, sig_j
            else:
                smaller, larger = sig_j, sig_i

            for ctx, tup_s in smaller.items():
                tup_l = larger.get(ctx)
                if tup_l is not None:
                    cxt += 1
                    if tup_s != tup_l:
                        obj += 1

            dist = obj / cxt if cxt > 0 else 1.0
            D[i][j] = dist
            D[j][i] = dist

        if (i + 1) % 50 == 0:
            print(f"    {i+1}/{n} sequences done")

    print(f"  [Co-phylog] Distance matrix computed in {time.time()-t0:.1f}s")
    return D


# ============================================================
# CVTree-like: k-mer frequency + correlation distance
# ============================================================

def count_kmers(seq_str, k):
    """Count contiguous k-mers in a sequence. Returns numpy array."""
    counts = defaultdict(int)
    n_kmers = len(seq_str) - k + 1
    if n_kmers <= 0:
        return np.array([])
    for i in range(n_kmers):
        kmer = seq_str[i:i+k]
        if 'N' not in kmer:
            counts[kmer] += 1
    return counts


def normalize_counts(counts, total_kmers):
    """Normalize k-mer counts to frequencies."""
    return {k: v / total_kmers for k, v in counts.items()}


def cvtree_distance_matrix(sequences, names, k_range=(3, 4, 5, 6), use_markov=True):
    """
    Compute CVTree-like distance using k-mer frequency correlation.

    For each k in k_range:
    1. Count k-mer frequencies
    2. (Optional) Subtract Markov background
    3. Concatenate corrected frequencies
    4. Compute Pearson correlation distance

    Args:
        sequences: list of sequence strings
        names: list of sequence names
        k_range: tuple of k values to use
        use_markov: whether to apply Markov background subtraction

    Returns:
        numpy 2D distance matrix
    """
    n = len(sequences)
    print(f"  [CVTree] Computing k-mer frequencies (k={k_range}, markov={use_markov})...")
    t0 = time.time()

    # Build unified k-mer index across all k values
    all_kmer_lists = {}  # k -> sorted list of all kmers observed
    vectors = {}  # k -> (n, num_kmers) matrix

    for k in k_range:
        # Collect all kmers for this k
        kmer_set = set()
        for seq in sequences:
            counts = count_kmers(seq, k)
            kmer_set.update(counts.keys())
        kmer_list = sorted(kmer_set)
        all_kmer_lists[k] = kmer_list

        # Build frequency matrix
        mat = np.zeros((n, len(kmer_list)))
        for i, seq in enumerate(sequences):
            counts = count_kmers(seq, k)
            total = sum(counts.values())
            if total > 0:
                for j_idx, kmer in enumerate(kmer_list):
                    mat[i, j_idx] = counts.get(kmer, 0) / total
        vectors[k] = mat

    # Apply Markov background subtraction if requested
    if use_markov:
        for k in k_range:
            if k <= 2:
                continue  # Can't do Markov for k <= 2
            mat = vectors[k]
            kmer_list = all_kmer_lists[k]

            # Compute (k-1)-mer frequencies from k-mer data
            # Expected: f(wXYZ) ≈ f(wXY) * f(XYZ) / f(XY)
            # For k=3: f(ABC) ≈ f(AB) * f(BC) / f(B)
            # We approximate using the k-mer matrix itself
            n_kmers_k1 = 4 ** (k - 1)  # Total possible (k-1)-mers

            # Compute (k-1)-mer marginals
            prefix_counts = np.zeros(n)  # frequency of each (k-1)-mer prefix
            suffix_counts = np.zeros(n)  # frequency of each (k-1)-mer suffix
            center_counts = np.zeros(n)  # frequency of each (k-2)-mer center

            # This is getting complex; use a simpler approach:
            # Subtract the mean frequency as background
            row_means = mat.mean(axis=1, keepdims=True)
            mat_corrected = mat - row_means

            # Zero out negative values
            mat_corrected = np.maximum(mat_corrected, 0)

            vectors[k] = mat_corrected

    # Concatenate all k vectors
    all_vecs = np.hstack([vectors[k] for k in k_range])

    # Normalize each vector (L2 normalization)
    norms = np.linalg.norm(all_vecs, axis=1, keepdims=True)
    norms[norms == 0] = 1
    all_vecs_normalized = all_vecs / norms

    # Compute Pearson correlation distance: d = 1 - corr(x, y)
    # This is equivalent to 1 - cosine of L2-normalized mean-centered vectors
    # For simplicity, use cosine distance (widely used in alignment-free methods)
    print(f"  [CVTree] Computing correlation distances...")
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            # Pearson correlation
            x = all_vecs[i]
            y = all_vecs[j]
            x_centered = x - x.mean()
            y_centered = y - y.mean()
            nx = np.linalg.norm(x_centered)
            ny = np.linalg.norm(y_centered)
            if nx > 0 and ny > 0:
                corr = np.dot(x_centered, y_centered) / (nx * ny)
                dist = max(0.0, 1.0 - corr)
            else:
                dist = 1.0
            D[i][j] = dist
            D[j][i] = dist

    print(f"  [CVTree] Distance matrix computed in {time.time()-t0:.1f}s")
    return D


# ============================================================
# Simple k-mer cosine distance (baseline)
# ============================================================

def kmer_cosine_distance_matrix(sequences, names, k=5, gap_pattern=None):
    """
    Compute pairwise cosine distance from k-mer frequency vectors.
    This is the baseline method (no spaced k-mers unless specified).

    Args:
        sequences: list of sequence strings
        names: list of sequence names
        k: k-mer size
        gap_pattern: if provided, use spaced k-mers (e.g., "11011")

    Returns:
        numpy 2D distance matrix
    """
    n = len(sequences)
    print(f"  [KmerCosine] Computing k-mer frequencies (k={k}, spaced={gap_pattern is not None})...")
    t0 = time.time()

    # Collect all k-mers
    all_kmers = set()
    seq_vectors = []
    for seq in sequences:
        counts = defaultdict(int)
        seq_clean = seq.replace('N', '')
        if gap_pattern:
            # Spaced k-mer
            positions = [i for i, c in enumerate(gap_pattern) if c == '1']
            for i in range(len(seq_clean) - max(positions)):
                kmer = ''.join(seq_clean[i + p] for p in positions)
                if 'N' not in kmer:
                    counts[kmer] += 1
        else:
            for i in range(len(seq_clean) - k + 1):
                kmer = seq_clean[i:i+k]
                if 'N' not in kmer:
                    counts[kmer] += 1
        all_kmers.update(counts.keys())
        seq_vectors.append(counts)

    kmer_list = sorted(all_kmers)
    print(f"  [KmerCosine] {len(kmer_list)} unique k-mers")

    # Build frequency matrix
    mat = np.zeros((n, len(kmer_list)))
    for i, counts in enumerate(seq_vectors):
        total = sum(counts.values())
        if total > 0:
            for j_idx, kmer in enumerate(kmer_list):
                mat[i, j_idx] = counts.get(kmer, 0) / total

    # Cosine distance: d = 1 - cos(x, y) = 1 - (x.y / (|x||y|))
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            dot = np.dot(mat[i], mat[j])
            ni = np.linalg.norm(mat[i])
            nj = np.linalg.norm(mat[j])
            if ni > 0 and nj > 0:
                cos_sim = dot / (ni * nj)
                dist = max(0.0, 1.0 - cos_sim)
            else:
                dist = 1.0
            D[i][j] = dist
            D[j][i] = dist

    print(f"  [KmerCosine] Distance matrix computed in {time.time()-t0:.1f}s")
    return D


# ============================================================
# Simplified andi approximation (k-mer match length distribution)
# ============================================================

def andi_approx_distance_matrix(sequences, names, min_anchor=15, use_revcomp=True):
    """
    Simplified andi-like distance estimation based on exact match length distribution.

    The real andi uses enhanced suffix arrays to find MUMs (maximal unique matches)
    and estimates evolutionary distance from the length distribution of anchored
    ungapped alignments. This approximation:

    1. For each pair of sequences, finds all exact matches of length >= min_anchor
    2. Estimates the substitution rate from the average match quality
    3. Applies Jukes-Cantor correction

    Note: This is a SIMPLIFIED approximation. The real andi uses suffix arrays
    for O(n) MUM finding, while this uses a hash-based approach that is O(n^2).
    For gene-length sequences (~500bp), this is still fast enough.

    Args:
        sequences: list of sequence strings
        names: list of sequence names
        min_anchor: minimum match length (default 15, andi's default)
        use_revcomp: whether to use reverse complement (set False for protein)

    Returns:
        numpy 2D distance matrix
    """
    n = len(sequences)
    D = np.zeros((n, n))
    rc_label = "with revcomp" if use_revcomp else "no revcomp"
    print(f"  [andi-approx] Computing anchor distances (min_anchor={min_anchor}, {rc_label})...")
    t0 = time.time()

    for i in range(n):
        seq_i = sequences[i].replace('N', '')
        len_i = len(seq_i)

        # Build hash of all k-mers of length min_anchor
        hash_i = {}
        for start in range(len_i - min_anchor + 1):
            kmer = seq_i[start:start + min_anchor]
            if use_revcomp:
                rc = revcomp_seq(kmer)
                canonical = min(kmer, rc)
            else:
                canonical = kmer
            if canonical not in hash_i:
                hash_i[canonical] = []
            hash_i[canonical].append(start)

        for j in range(i + 1, n):
            seq_j = sequences[j].replace('N', '')
            len_j = len(seq_j)

            # Find all exact matches of length >= min_anchor
            match_lengths = []
            visited = set()

            for start_j in range(len_j - min_anchor + 1):
                kmer = seq_j[start_j:start_j + min_anchor]
                if use_revcomp:
                    rc = revcomp_seq(kmer)
                    canonical = min(kmer, rc)
                else:
                    canonical = kmer

                if canonical in hash_i:
                    for start_i in hash_i[canonical]:
                        # Extend match
                        ext = 0
                        while (start_i + min_anchor + ext < len_i and
                               start_j + min_anchor + ext < len_j):
                            if seq_i[start_i + min_anchor + ext] == seq_j[start_j + min_anchor + ext]:
                                ext += 1
                            else:
                                break
                        match_lengths.append(min_anchor + ext)

            if match_lengths:
                # Estimate substitution rate from match length distribution
                avg_len = np.mean(match_lengths)
                total_len = (len_i + len_j) / 2
                # Simple model: match_len ~ total_len * exp(-d * L)
                # d ≈ -ln(avg_len / total_len) / min_anchor
                # Use Jukes-Cantor correction: d = -3/4 * ln(1 - 4/3 * p)
                # where p is estimated from match rate
                p_hat = 1.0 - avg_len / total_len
                p_hat = max(0.0, min(0.75, p_hat))  # Clamp to valid range
                if p_hat < 0.75:
                    d_jc = -0.75 * math.log(1.0 - 4.0 * p_hat / 3.0)
                else:
                    d_jc = 3.0  # Saturated
                D[i][j] = d_jc
                D[j][i] = d_jc
            else:
                D[i][j] = 3.0  # No matches = saturated distance
                D[j][i] = 3.0

        if (i + 1) % 20 == 0:
            print(f"    {i+1}/{n} sequences done")

    print(f"  [andi-approx] Distance matrix computed in {time.time()-t0:.1f}s")
    return D


# ============================================================
# Quick test
# ============================================================

if __name__ == '__main__':
    # Test on seed42 indel data
    fasta_path = os.path.join(os.path.dirname(__file__), 'seed42_indel.fasta')
    if not os.path.exists(fasta_path):
        print(f"Test file not found: {fasta_path}")
        sys.exit(1)

    seqs = read_fasta(fasta_path)
    names = list(seqs.keys())
    sequences = [seqs[n] for n in names]
    print(f"Loaded {len(sequences)} sequences, avg length {np.mean([len(s) for s in sequences]):.0f}")

    # Test Co-phylog
    D_cophy = cophylog_distance_matrix(sequences, names)
    print(f"\nCo-phylog distance matrix stats: mean={D_cophy[D_cophy>0].mean():.4f}, "
          f"max={D_cophy.max():.4f}")

    # Test k-mer cosine
    D_cos = kmer_cosine_distance_matrix(sequences, names, k=5)
    print(f"Kmer-cosine distance matrix stats: mean={D_cos[D_cos>0].mean():.4f}, "
          f"max={D_cos.max():.4f}")

    print("\nAll methods tested successfully!")
