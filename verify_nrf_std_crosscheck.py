#!/usr/bin/env python3
"""Independent cross-check of nrf_standard.py using Bio.Phylo.

Splits are extracted with a completely different parser (Biopython) and a
different traversal (clade-based instead of recursive-tuple based).  The two
implementations must agree on (a) the number of non-trivial splits for every
fully resolved tree (n-3) and (b) the exact RF distance / nRF for every pair.
"""
import sys
import glob
import os
import random

from Bio import Phylo
from io import StringIO

import nrf_standard as ns

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "indelible_n200_indel")


def splits_biopython(nwk):
    tree = Phylo.read(StringIO(nwk), "newick")
    all_leaves = {t.name for t in tree.get_terminals()}
    ref = min(all_leaves)
    out = set()
    for clade in tree.find_clades():
        if clade.is_terminal():
            continue
        ls = frozenset(t.name for t in clade.get_terminals())
        if len(ls) < 2 or len(ls) > len(all_leaves) - 2:
            continue
        comp = all_leaves - ls
        out.add(ls if ref not in ls else frozenset(comp))
    return out, len(all_leaves)


def rf_bio(a_nwk, b_nwk):
    sa, n = splits_biopython(a_nwk)
    sb, _ = splits_biopython(b_nwk)
    return len(sa - sb) + len(sb - sa), n


def main():
    random.seed(7)
    trues = sorted(glob.glob(os.path.join(DATA, "seed*_indel_true.nwk")))
    ft2s = sorted(glob.glob(os.path.join(DATA, "seed*_indel_ft2.nwk")))
    fuses = sorted(glob.glob(os.path.join(DATA, "seed*_indel_fusang.nwk")))
    # 12 random seeds, three pair types each
    checks = 0
    mism = 0
    for i in random.sample(range(len(trues)), 12):
        with open(trues[i]) as fh:
            t = fh.read()
        with open(ft2s[i]) as fh:
            f = fh.read()
        with open(fuses[i]) as fh:
            u = fh.read()
        for a_name, a, b_name, b in (("true", t, "ft2", f), ("true", t, "fusang", u),
                                     ("ft2", f, "fusang", u)):
            rf_ref, n_ref = ns.rf_distance(a, b), len(ns.leaves(ns.parse_newick(a)))
            rf_bio_v, n_bio = rf_bio(a, b)
            ok = (rf_ref == rf_bio_v) and (n_ref == n_bio)
            checks += 1
            if not ok:
                mism += 1
                print("MISMATCH", trues[i], a_name, b_name, rf_ref, rf_bio_v, n_ref, n_bio)
        # split-count contract on the resolved trees
        for nm, s in (("true", t), ("fusang", u)):
            sb, n = splits_biopython(s)
            ss = ns.splits(ns.parse_newick(s))
            checks += 1
            if not (sb == ss and len(ss) == n - 3):
                mism += 1
                print("SPLIT MISMATCH", trues[i], nm, len(ss), len(sb), n)
    print("checked %d cases, %d mismatches" % (checks, mism))
    return 1 if mism else 0


if __name__ == "__main__":
    sys.exit(main())
