"""
INDELible data source adapter for the boundary-classifier pipeline.

Replaces the old JC69 simulator (``tree_simulation.make_coalescent_tree`` +
``simulate_seqs``) with datasets produced by INDELible 1.03 via
``gen_indelible_sim.py`` (which must be run inside WSL, since INDELible is not
available on Windows).

Conventions
-----------
``gen_indelible_sim.py`` writes, for every (seed, tag) job::

    <pool>/seed{seed}_{tag}.fasta        # unaligned, variable-length sequences
    <pool>/seed{seed}_{tag}_true.nwk     # true tree, leaves named seq0..seqN-1

The rest of the pipeline expects leaves named ``t0001`` .. ``tNNNN``, so this
module renames leaves on load.  Renaming MUST happen before
``training_data_generator.get_true_clades()`` is called, otherwise the
sub-cluster labels silently collapse to 0.

Public API mirrors the old generator functions so call sites stay unchanged::

    load_indelible_tree(pool, tag, seed)            -> (newick, seqs, names)
    load_indelible_structured(pool, [(tag, seed)])  -> (newick, seqs, names)
"""

import os
import re
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_POOL = os.path.join(_ROOT, "indelible_bc", "pool")

_SEQ_RE = re.compile(r"seq(\d+)(?![0-9])")


def indelible_paths(pool, tag, seed):
    """Return (fasta_path, newick_path) for one INDELible job."""
    return (
        os.path.join(pool, "seed%d_%s.fasta" % (seed, tag)),
        os.path.join(pool, "seed%d_%s_true.nwk" % (seed, tag)),
    )


def read_fasta(path):
    """Read a FASTA file into an ordered {name: sequence} dict."""
    records = {}
    order = []
    cur = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                cur = line[1:].split()[0]
                records[cur] = []
                order.append(cur)
            elif cur is not None:
                records[cur].append(line)
    return [(n, "".join(records[n])) for n in order]


def rename_newick(nwk, mapping):
    """Rewrite leaf names in a Newick string using ``mapping`` (old -> new)."""
    def _sub(m):
        old = "seq%s" % m.group(1)
        return mapping.get(old, old)
    return _SEQ_RE.sub(_sub, nwk)


def _load_one(pool, tag, seed, offset):
    """Load one INDELible dataset; rename leaves to t{offset+1:04d}..."""
    fas, nwk_path = indelible_paths(pool, tag, seed)
    if not os.path.exists(fas):
        raise FileNotFoundError("INDELible FASTA missing: %s" % fas)

    records = read_fasta(fas)
    mapping = {}
    names = []
    for i, (old, _seq) in enumerate(records):
        new = "t%04d" % (offset + i + 1)
        mapping[old] = new
        names.append(new)
    sequences = [seq for _old, seq in records]

    if os.path.exists(nwk_path):
        with open(nwk_path) as fh:
            nwk = fh.read().strip()
    else:
        raise FileNotFoundError("INDELible true tree missing: %s" % nwk_path)

    nwk = rename_newick(nwk, mapping)
    return nwk, sequences, names


def load_indelible_tree(pool, tag, seed):
    """Single-tree dataset (the coalescent / homogeneous case)."""
    nwk, seqs, names = _load_one(pool, tag, seed, 0)
    return nwk, seqs, names


def load_indelible_structured(pool, jobs):
    """Concatenate several independently simulated clades into one tree.

    ``jobs`` is a list of (tag, seed) pairs, one per clade.  The clade Newick
    strings are joined exactly like ``training_data_generator``
    ``generate_structured_tree()`` does, i.e. ``(clade0,clade1,...);`` with a
    zero-length root stem on each clade so that ``get_true_clades()`` still
    recognises them.
    """
    parts = []
    all_seqs = []
    all_names = []
    offset = 0
    for tag, seed in jobs:
        nwk, seqs, names = _load_one(pool, tag, seed, offset)
        nwk = nwk.rstrip().rstrip(";").rstrip()
        if nwk.endswith(")"):
            nwk += ":0.000000"
        parts.append(nwk)
        all_seqs.extend(seqs)
        all_names.extend(names)
        offset += len(names)
    combined = "(" + ",".join(parts) + ");"
    return combined, all_seqs, all_names


def tag_for(config_idx, tree_idx, clade_idx=None):
    """Deterministic job tag for a training-tree job."""
    if clade_idx is None:
        return "tr_c%d_t%d" % (config_idx, tree_idx)
    return "tr_c%d_t%d_k%d" % (config_idx, tree_idx, clade_idx)


def tag_for_e2e(tree_type, n, label, seed_i, clade_idx=None):
    """Deterministic job tag for an E2E scenario job."""
    if clade_idx is None:
        return "e2e_%s_n%d_%s_s%d" % (tree_type, n, label, seed_i)
    return "e2e_%s_n%d_%s_s%d_k%d" % (tree_type, n, label, seed_i, clade_idx)


if __name__ == "__main__":
    pool = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_POOL
    for fn in sorted(os.listdir(pool))[:5]:
        print(fn)
