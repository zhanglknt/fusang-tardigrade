#!/usr/bin/env python3
"""Generate FASTA with real phylogenetic signal.
Generates a random tree, then simulates sequences along branches
using a simple substitution model (Jukes-Cantor like).

Usage: python gen_test_data.py <n> <length> <sub_rate> <seed>
Output: test_n{N}.fasta + test_n{N}_true.nwk
"""

import sys
import random
import os

def generate_tree_newick(n, seed=42):
    """Generate random tree in NEWICK format using coalescent-like method."""
    random.seed(seed)
    tips = [f"seq{i}" for i in range(n)]
    active = [(tips[i], tips[i]) for i in range(n)]
    ic = 0
    while len(active) > 1:
        idxs = sorted(random.sample(range(len(active)), 2), reverse=True)
        j, i = idxs[0], idxs[1]
        left, right = active[i], active[j]
        inode = f"I{ic}"
        ic += 1
        node_str = f"({left[0]}:0.1000,{right[0]}:0.1000)"
        active.pop(j)
        active.pop(i)
        active.append((node_str, inode))
    return active[0][0] + ";"


def simulate_sequences(true_nwk_path, fasta_path, length=1000, sub_rate=0.3, seed=123):
    """Simulate sequences along the true tree.
    Each branch accumulates mutations ~Poisson(sub_rate * branch_len).
    """
    import re
    from io import StringIO
    from Bio import Phylo
    from Bio.Phylo.BaseTree import Clade

    random.seed(seed)

    with open(true_nwk_path) as f:
        nwk = f.read().strip()

    tree = Phylo.read(StringIO(nwk), "newick")

    bases = ['A', 'T', 'C', 'G']

    # DFS: assign sequence to each node
    # Root gets random sequence
    root_seq = [random.choice(bases) for _ in range(length)]
    seqs = {tree.root: root_seq}

    def dfs(node):
        for child in node.clades:
            parent_seq = seqs[node]
            bl = child.branch_length if child.branch_length else 0.1
            # Number of mutations on this branch
            n_mut = max(0, int(bl * sub_rate * length))
            child_seq = parent_seq[:]
            if n_mut > 0:
                positions = random.sample(range(length), min(n_mut, length))
                for pos in positions:
                    cur = child_seq[pos]
                    choices = [b for b in bases if b != cur]
                    child_seq[pos] = random.choice(choices)
            seqs[child] = child_seq
            dfs(child)

    dfs(tree.root)

    # Write FASTA
    terminals = {leaf.name: leaf for leaf in tree.get_terminals()}
    with open(fasta_path, 'w') as f:
        for name in sorted(terminals.keys()):
            leaf = terminals[name]
            seq = "".join(seqs[leaf])
            f.write(f">{name}\n")
            for k in range(0, len(seq), 80):
                f.write(seq[k:k+80] + "\n")

    print(f"  Generated {fasta_path}: {len(terminals)} seqs, length={length}")


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 50
    length = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
    sub_rate = float(sys.argv[3]) if len(sys.argv) > 3 else 0.3
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 42

    base = f"test_n{n}"
    fasta = f"{base}.fasta"
    true_nwk = f"{base}_true.nwk"

    print(f"Generating test data: n={n}, L={length}, sub_rate={sub_rate}, seed={seed}")

    nwk_str = generate_tree_newick(n, seed=seed)
    with open(true_nwk, 'w') as f:
        f.write(nwk_str)
    print(f"  True tree: {true_nwk}")
    print(f"  NEWICK: {nwk_str[:100]}...")

    simulate_sequences(true_nwk, fasta, length=length, sub_rate=sub_rate, seed=seed+1)

    print(f"Done: {fasta}, {true_nwk}")
