"""
P0-5: Mash vs IMMI benchmark results computation.
Computes IMMI (k-mer NJ) nRF vs TRUE trees for 30 seeds (clean + indel).
Compares against pre-computed Mash trees.
"""
import json
import os
import sys
import numpy as np
from Bio import Phylo, SeqIO
from io import StringIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from fusang_mhl.level0_kmer import compute_l0_distance
from scipy.stats import wilcoxon


def compute_nrf(true_nwk: str, inferred_nwk: str):
    """Compute normalized Robinson-Foulds distance."""
    try:
        true_tree = Phylo.read(StringIO(true_nwk), 'newick')
        inf_tree = Phylo.read(StringIO(inferred_nwk), 'newick')
        from Bio.Phylo.Consensus import _BitString

        def bitstrs(tree):
            terms = sorted(t.name for t in tree.get_terminals())
            result = {}
            def visit(clade):
                if clade.is_terminal():
                    return
                for c in clade:
                    visit(c)
                term_set = {t.name for t in clade.get_terminals()}
                bs = _BitString(''.join('1' if t in term_set else '0' for t in terms))
                result[bs] = True
            for c in true_tree.root:
                visit(c)
            return set(result.keys())

        # Use both trees for bitstrings
        def get_bitstrs(tree):
            terms = sorted(t.name for t in tree.get_terminals())
            result = set()
            def visit(clade):
                if clade.is_terminal():
                    return
                for c in clade:
                    visit(c)
                term_set = {t.name for t in clade.get_terminals()}
                bs = _BitString(''.join('1' if t in term_set else '0' for t in terms))
                result.add(bs)
            for c in tree.root:
                visit(c)
            return result

        # Use consistent term ordering from true tree
        terms = sorted(t.name for t in true_tree.get_terminals())

        def get_bitstrs_with_terms(tree, terms):
            result = set()
            def visit(clade):
                if clade.is_terminal():
                    return
                for c in clade:
                    visit(c)
                term_set = {t.name for t in clade.get_terminals()}
                bs = _BitString(''.join('1' if t in term_set else '0' for t in terms))
                result.add(bs)
            for c in tree.root:
                visit(c)
            return result

        bs1 = get_bitstrs_with_terms(true_tree, terms)
        bs2 = get_bitstrs_with_terms(inf_tree, terms)
        rf = len(bs1.symmetric_difference(bs2))
        n = len(terms)
        return rf / (2 * (n - 3))
    except Exception as e:
        print(f'  nRF error: {e}', flush=True)
        return None


def build_immi_tree(fasta_path: str, names: list) -> str:
    """Build IMMI (k-mer NJ) tree, return newick string."""
    seqs_dict = {}
    for rec in SeqIO.parse(fasta_path, 'fasta'):
        seqs_dict[rec.id] = str(rec.seq)
    seqs_list = [seqs_dict[n] for n in names]
    D = compute_l0_distance(seqs_list, names)
    n = len(names)
    lower_tri = [[float(D[i][j]) for j in range(i + 1)] for i in range(n)]
    dm = DistanceMatrix(list(names), lower_tri)
    tree = DistanceTreeConstructor().nj(dm)
    buf = StringIO()
    Phylo.write(tree, buf, 'newick')
    return buf.getvalue().strip()


def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    benchmark_dir = os.path.join(base_dir, 'benchmark_mash')

    print("=" * 60)
    print("P0-5: Computing IMMI vs Mash nRF (30 seeds)")
    print("=" * 60)

    clean_nrfs = []
    indel_nrfs = []

    for seed in range(1, 31):
        for cond, nrfs in [('clean', clean_nrfs), ('indel', indel_nrfs)]:
            fasta = os.path.join(benchmark_dir, cond, f'seed{seed:03d}.fasta')
            true_nwk_path = os.path.join(benchmark_dir, cond, f'seed{seed:03d}_true.nwk')
            if not os.path.exists(fasta) or not os.path.exists(true_nwk_path):
                print(f'  MISSING: {fasta}', flush=True)
                continue
            with open(true_nwk_path) as f:
                true_nwk = f.read().strip()

            # Get names from FASTA
            names = [rec.id for rec in SeqIO.parse(fasta, 'fasta')]
            try:
                immi_nwk = build_immi_tree(fasta, names)
                nrf = compute_nrf(true_nwk, immi_nwk)
                if nrf is not None:
                    nrfs.append(nrf)
                    if seed <= 3 or seed == 30:
                        print(f'  {cond} seed{seed:02d}: IMMI nRF = {nrf:.4f}', flush=True)
            except Exception as e:
                print(f'  Error {cond} seed{seed}: {e}', flush=True)

    print(f'\nIMMI Clean: {np.mean(clean_nrfs):.4f} +/- {np.std(clean_nrfs, ddof=1):.4f} (n={len(clean_nrfs)})')
    print(f'IMMI Indel: {np.mean(indel_nrfs):.4f} +/- {np.std(indel_nrfs, ddof=1):.4f} (n={len(indel_nrfs)})')

    # Mash comparison (single tree, per previously computed files)
    mash_clean_nrf = None
    mash_indel_nrf = None
    for cond, mash_file, var_name in [
        ('clean', 'mashtree_clean.nwk', 'mash_clean_nrf'),
        ('indel', 'mashtree_indel02.nwk', 'mash_indel_nrf')
    ]:
        mfp = os.path.join(base_dir, mash_file)
        if os.path.exists(mfp):
            with open(mfp) as f:
                lines = f.read().strip().split('\n')
                mash_nwk = lines[-1] if len(lines) > 1 else lines[0]
            # Compare vs seed001 true tree
            true_path = os.path.join(benchmark_dir, cond, 'seed001_true.nwk')
            if os.path.exists(true_path):
                with open(true_path) as f:
                    true_nwk = f.read().strip()
                nrf = compute_nrf(true_nwk, mash_nwk)
                if cond == 'clean':
                    mash_clean_nrf = nrf
                else:
                    mash_indel_nrf = nrf
                print(f'Mash {cond} (vs seed001): nRF = {nrf:.4f}', flush=True)
        else:
            print(f'Mash file not found: {mfp}', flush=True)

    # Wilcoxon test: clean vs indel
    if len(clean_nrfs) >= 10 and len(indel_nrfs) >= 10:
        min_len = min(len(clean_nrfs), len(indel_nrfs))
        stat, p = wilcoxon(indel_nrfs[:min_len], clean_nrfs[:min_len])
        print(f'\nWilcoxon (indel vs clean): stat={stat:.1f}, p={p:.4e}')

    # Save results
    results = {
        'immi_clean': {
            'mean': float(np.mean(clean_nrfs)),
            'std': float(np.std(clean_nrfs, ddof=1)),
            'n': len(clean_nrfs)
        },
        'immi_indel': {
            'mean': float(np.mean(indel_nrfs)),
            'std': float(np.std(indel_nrfs, ddof=1)),
            'n': len(indel_nrfs)
        },
        'mash_clean_vs_seed1': mash_clean_nrf,
        'mash_indel_vs_seed1': mash_indel_nrf,
        'note': 'IMMI=k-mer NJ k=5 gap2 cosine, vs TRUE coalescent trees. Mash=single tree vs seed001 true tree.',
        'immi_clean_nrfs': [round(x, 4) for x in clean_nrfs],
        'immi_indel_nrfs': [round(x, 4) for x in indel_nrfs],
    }
    out_path = os.path.join(base_dir, 'mash_vs_immi_results.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)

    print(f'\n=== FINAL SUMMARY (30 seeds, vs TRUE tree) ===')
    print(f'IMMI (Fusang) Clean: {results["immi_clean"]["mean"]:.4f} +/- {results["immi_clean"]["std"]:.4f}')
    print(f'IMMI (Fusang) Indel: {results["immi_indel"]["mean"]:.4f} +/- {results["immi_indel"]["std"]:.4f}')
    print(f'Mash Clean (vs seed001): {mash_clean_nrf}')
    print(f'Mash Indel (vs seed001): {mash_indel_nrf}')
    if mash_indel_nrf and results["immi_indel"]["mean"] > 0:
        fold = mash_indel_nrf / results["immi_indel"]["mean"]
        print(f'Fold improvement (Mash/IMMI on indel): {fold:.2f}x')
    print(f'Results saved to: {out_path}')


if __name__ == '__main__':
    main()
