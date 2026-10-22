"""
Standard root-independent (unrooted) Robinson-Foulds normalised distance.

nRF = RF / (2 * (n - 3))

where RF = |S1 - S2| + |S2 - S1| is the symmetric difference over the sets of
NON-TRIVIAL splits (bipartitions) of the two unrooted trees, and n is the number
of leaves (shared leaf set).  Denominator 2(n-3) is the maximum possible RF for
two fully resolved trees on n leaves (each has n-3 internal edges).

Split extraction is root-independent: for every internal edge we take the leaf set
of one side, reduced to its canonical form (the side NOT containing the
lexicographically smallest leaf).  This is the standard definition and differs
from calc_nrf_simple.py, whose per-node "smaller child clade" rule drops the
complementary side of each split and therefore under-counts.

Validation contract: for a fully resolved unrooted tree with n leaves the number
of non-trivial splits returned must be exactly n - 3.
"""
from __future__ import annotations

import re
from typing import Iterable, Set, Tuple, FrozenSet


def read_newick(path_or_str: str) -> str:
    if "\n" not in path_or_str and len(path_or_str) < 4096 and not path_or_str.strip().startswith("("):
        with open(path_or_str) as fh:
            return fh.read().strip()
    return path_or_str.strip()


def _tokenize(nwk: str):
    return re.findall(r"\(|\)|,|;|[^(),;\s]+", nwk)


def parse_newick(nwk: str):
    """Return nested tuple structure: ('leaf', name) or ('node', [children])."""
    tokens = _tokenize(nwk)
    pos = 0

    def parse_node():
        nonlocal pos
        if tokens[pos] == "(":
            pos += 1
            children = []
            while True:
                children.append(parse_node())
                if tokens[pos] == ",":
                    pos += 1
                    continue
                if tokens[pos] == ")":
                    pos += 1
                    break
                raise ValueError(f"unexpected token {tokens[pos]!r}")
            # optional label / branch length on internal node
            if pos < len(tokens) and tokens[pos] not in ("(", ")", ",", ";"):
                pos += 1
            return ("node", children)
        else:
            label = tokens[pos]
            if ":" in label:
                label = label.split(":")[0]
            pos += 1
            if pos < len(tokens) and tokens[pos] == ":":
                pos += 1
                if pos < len(tokens):
                    pos += 1
            return ("leaf", label)

    tree = parse_node()
    return tree


def leaves(tree) -> list:
    kind, payload = tree
    if kind == "leaf":
        return [payload]
    out = []
    for ch in payload:
        out.extend(leaves(ch))
    return out


def splits(nwk_or_tree) -> Set[FrozenSet[str]]:
    """Non-trivial splits, canonicalised so the side lacking the smallest leaf is kept."""
    tree = parse_newick(nwk_or_tree) if isinstance(nwk_or_tree, str) else nwk_or_tree
    all_leaves = sorted(set(leaves(tree)))
    ref = all_leaves[0]
    result: Set[FrozenSet[str]] = set()

    def walk(node):
        kind, payload = node
        if kind == "leaf":
            return {payload}
        acc = set()
        for ch in payload:
            acc |= walk(ch)
        # for every child subtree, the split is (child_leaves | complement)
        if len(acc) == len(all_leaves):
            # root: each child side gives a split
            for ch in payload:
                side = set(leaves(ch))
                result.add(_canon(side, all_leaves, ref))
        else:
            side = set(acc)
            result.add(_canon(side, all_leaves, ref))
        return acc

    walk(tree)
    # drop trivial splits: size 1, size n-1 (pendant edge, which surfaces when the
    # Newick root has a leaf child, e.g. "(seq0,(...));") and the full set
    return {s for s in result if 1 < len(s) < len(all_leaves) - 1}


def _canon(side: Set[str], all_leaves, ref: str) -> FrozenSet[str]:
    comp = set(all_leaves) - side
    if ref in side:
        return frozenset(comp)
    return frozenset(side)


def rf_distance(t1, t2) -> int:
    s1, s2 = splits(t1), splits(t2)
    return len(s1 - s2) + len(s2 - s1)


def nrf(t1, t2) -> float:
    s1, s2 = splits(t1), splits(t2)
    n = len(sorted(set(leaves(parse_newick(t1) if isinstance(t1, str) else t1))))
    return (len(s1 - s2) + len(s2 - s1)) / (2 * (n - 3))


def nrf_files(path1: str, path2: str) -> float:
    return nrf(read_newick(path1), read_newick(path2))


if __name__ == "__main__":
    import sys
    print(nrf_files(sys.argv[1], sys.argv[2]))
