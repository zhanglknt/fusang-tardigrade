#!/usr/bin/env python3
"""Compute nRF values for benchmarked tools."""
import sys
import re
from io import StringIO

try:
    from Bio import Phylo
except ImportError:
    print("Biopython not available!", file=sys.stderr)
    sys.exit(1)

def clean_newick(s):
    return re.sub(r':-[\d.]+', ':0.00001', s)

def get_splits(tree):
    tips = [l.name for l in tree.get_terminals()]
    all_tip_names = frozenset(tips)
    n_tips = len(tips)
    if n_tips < 3:
        return set()
    first_tip = list(tree.find_clades(name=tips[0]))[0]
    tree.root_with_outgroup(first_tip)
    splits = set()
    for clade in tree.get_nonterminals():
        leaf_names = frozenset(l.name for l in clade.get_terminals())
        m = len(leaf_names)
        if 1 <= m < n_tips and leaf_names != all_tip_names:
            splits.add(min(leaf_names, all_tip_names - leaf_names, key=sorted))
    return splits

def calc_nrf(t1, t2):
    s1, s2 = get_splits(t1), get_splits(t2)
    if not s1 or not s2:
        return 1.0
    rf = len(s1 - s2) + len(s2 - s1)
    n = len(t1.get_terminals())
    return rf / (2.0 * (n - 3)) if n > 3 else 1.0

def compute_nrf(est_file, true_file):
    with open(est_file, 'r') as f:
        content = f.read().strip()
        # Mashtree embeds log lines before the Newick; extract last line
        lines = content.split('\n')
        est = clean_newick(lines[-1].strip())
    with open(true_file, 'r') as f:
        true = clean_newick(f.read().strip())
    t1 = Phylo.read(StringIO(est), 'newick')
    t2 = Phylo.read(StringIO(true), 'newick')
    return calc_nrf(t1, t2)

if __name__ == '__main__':
    base = 'd:/系统发育树项目/Fusang/Fusang-main'

    pairs = [
        ('mashtree_clean.nwk', 'test_indel_n200_true.nwk', 'MashTree Clean'),
        ('mash_fastme_clean.nwk', 'test_indel_n200_true.nwk', 'Mash+FastME Clean'),
        ('mashtree_indel02.nwk', 'test_indel_n200_i002_true.nwk', 'MashTree Indel02'),
        ('mash_fastme_indel02.nwk', 'test_indel_n200_i002_true.nwk', 'Mash+FastME Indel02'),
    ]

    results = {}
    for est_name, true_name, label in pairs:
        est_path = f'{base}/{est_name}'
        true_path = f'{base}/{true_name}'
        try:
            nrf_val = compute_nrf(est_path, true_path)
            results[label] = nrf_val
            print(f'{label} nRF: {nrf_val:.6f}')
        except Exception as e:
            results[label] = None
            print(f'{label} FAILED: {e}')

    print('\nSUMMARY:')
    for k, v in results.items():
        if v is not None:
            print(f'  {k}: {v:.6f}')
        else:
            print(f'  {k}: FAILED')
