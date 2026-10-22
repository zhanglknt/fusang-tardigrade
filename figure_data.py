"""Single source of truth for every number drawn in the Fusang figure set.

Every figure script (``generate_fig*_v3.py``, ``generate_supp_s*_v3.py``) imports its
numbers from here so that the panels cannot drift from the manuscript. Nothing in this
module is hand-typed from the paper: each accessor reads a machine-readable artefact and
prints/returns what that artefact contains.

Authoritative artefacts
-----------------------
table1_nrfstandard/{cell}_nrfstandard.csv (+ {cell}_summary.json)
    The manuscript's single authoritative per-cell source (Table 1 six cells and Table 3
    three indel rates), scored with the standard root-independent unrooted RF. Used for
    n200_indel (Figure 2, S4, S6(A,D), Figure 1(D)), the three indel-rate cells
    (Figure 1(A,B)) and the n=500 / n=1000 clean cells (S6(A)). Each accessor recomputes
    the paired statistics from the per-seed rows and asserts agreement with the summary
    JSON before adopting the summary's values.
table_threshold_n1000_nrfstandard.csv + table_threshold_n1000_summary.json
    Figure 3(D) n=1000 simplified-vs-DCM+EPA re-test, 10 paired seeds.
indelible_n200_indel/competitors_n200_indel_130seeds_stdrf.csv (+ _summary.json)
    Per-seed TRUE-relative and FT2-relative nRF for six k-mer configurations plus the
    alignment-free competitors, same 130 seeds. Backs Figure 4(D) and S7.
indelible_iqtree2/iqtree2_n200_indel_nrfstandard.csv
    Figure 1(D) IQ-TREE2 GTR row.
bench_mash_indelible_n200_nrfstandard.csv
    Figure 1(D) Mash row (indel condition, seeds 100-129).
repro_ablations/ablation_dcm_seed42_indel_nrfstandard.json
    Figure 3(A) seed-42 indel DCM chain and the 11-setting max_group_size sweep.
balibase_rf_recomputed/balibase_rf_summary.json
    S6(C) per-family Fusang-vs-FastTree2 nRF on BAliBASE v3.0.
swisstree_rf_recomputed/swisstree_rf_summary.json (+ _paired_tests_intRF.csv)
    S6(B) SwissTree effect size.
bench_l3_coalescent_indelible_30seeds_nrfstandard.csv
    S6(B) multi-k L1 vs L0 effect size (Table 5 protocol).
real_data/afproject_genome/gene_dna/gene_dna_results.json
    Figure 6(D) 13 mitochondrial loci.

Conventions
-----------
* All RF maths is done on **integer** RF counts ``round(nRF * 2(n-3))`` so that Wilcoxon
  ties are exact (see git 7f7edfe); effect sizes and ranges are reported back in nRF units.
* Paired differences are always returned in the orientation named by the accessor.
"""
import csv
import os
from functools import lru_cache

import numpy as np
from scipy import stats

ROOT = os.path.dirname(os.path.abspath(__file__))
N200_NORM = 2 * (200 - 3)          # 394 for the n=200 INDELible benchmarks


def _p(*parts):
    return os.path.join(ROOT, *parts)


def _read_csv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


def _f(rows, col):
    return np.array([float(r[col]) for r in rows], dtype=float)


def _intrf(v, norm):
    return np.rint(np.asarray(v, dtype=float) * norm)


# --------------------------------------------------------------------------------------
# paired statistics
# --------------------------------------------------------------------------------------

def _wilcoxon_zsplit(d):
    """Two-sided Wilcoxon signed-rank p-value, zero-split (zsplit) convention.

    The manuscript's unified convention (Methods, commit c337787): every Wilcoxon test is
    run on integer RF counts with ``zero_method='zsplit'`` — zero differences stay in the
    ranking and their rank mass is split evenly between the two tails.  SciPy >= 1.17
    removed ``zero_method`` from :func:`scipy.stats.wilcoxon`, so the normal-approximation
    form (with tie correction) is implemented here; it reproduces every stored
    ``wilcoxon_p_zsplit`` value in competitors_n200_indel_130seeds_stdrf_summary.json
    and the manuscript's 0.0507 / 0.669 / 2.07e-14 to <=0.2% relative error.
    """
    d = np.asarray(d, dtype=float)
    n = len(d)
    if n == 0 or np.all(d == 0):
        return 1.0
    r = stats.rankdata(np.abs(d))
    zr = float(r[d == 0].sum())
    r_plus = float(r[d > 0].sum()) + zr / 2.0
    r_minus = float(r[d < 0].sum()) + zr / 2.0
    t = min(r_plus, r_minus)
    mn = n * (n + 1) / 4.0
    _, cnts = np.unique(np.abs(d), return_counts=True)
    tie_term = float((cnts ** 3 - cnts).sum())
    var = n * (n + 1) * (2 * n + 1) / 24.0 - tie_term / 48.0
    z = (t - mn) / np.sqrt(var)
    return float(2.0 * stats.norm.sf(abs(z)))


def paired(a, b, norm, orientation='b_minus_a'):
    """Paired statistics for two per-seed nRF vectors.

    ``orientation`` is one of
      'b_minus_a' : difference = b - a  (the default; positive => a better)
      'a_minus_b' : difference = a - b  (positive => b better)
    Rank statistics are computed on integer RF counts so that ties are exact, with the
    zero-split convention matching the manuscript's unified Wilcoxon protocol.
    """
    ia, ib = _intrf(a, norm), _intrf(b, norm)
    d = (ib - ia) if orientation == 'b_minus_a' else (ia - ib)
    if orientation == 'b_minus_a':
        win_a, win_b = int((d > 0).sum()), int((d < 0).sum())
    else:
        win_a, win_b = int((d < 0).sum()), int((d > 0).sum())
    sd = d.std(ddof=1)
    return dict(
        n=len(d), wins_a=win_a, wins_b=win_b, ties=int((d == 0).sum()),
        p=_wilcoxon_zsplit(d),
        d=float(d.mean() / sd) if sd > 0 else float('nan'),
        mean_a=float(np.asarray(a, float).mean()), sd_a=float(np.asarray(a, float).std(ddof=1)),
        mean_b=float(np.asarray(b, float).mean()), sd_b=float(np.asarray(b, float).std(ddof=1)),
        diffs_nrf=(np.asarray(b, float) - np.asarray(a, float)) if orientation == 'b_minus_a'
                  else (np.asarray(a, float) - np.asarray(b, float)),
        diff_range=(float(min(np.asarray(b, float) - np.asarray(a, float))),
                    float(max(np.asarray(b, float) - np.asarray(a, float)))),
    )


def boot_ci(a, b, norm, n_boot=10000, seed=20250101, orientation='b_minus_a'):
    """Bootstrap 95% CI for the paired Cohen's d, same orientation as :func:`paired`."""
    ia, ib = _intrf(a, norm), _intrf(b, norm)
    d = (ib - ia) if orientation == 'b_minus_a' else (ia - ib)
    rng = np.random.default_rng(seed)
    samp = d[rng.integers(0, len(d), size=(n_boot, len(d)))]
    dz = samp.mean(axis=1) / samp.std(axis=1, ddof=1)
    return tuple(np.percentile(dz, [2.5, 97.5]).tolist())


def boot_mean_ci(a, b, norm, n_boot=10000, seed=20250101, orientation='b_minus_a'):
    """Bootstrap 95% CI for the paired mean difference (nRF units)."""
    ia, ib = _intrf(a, norm), _intrf(b, norm)
    d = (ib - ia) if orientation == 'b_minus_a' else (ia - ib)
    rng = np.random.default_rng(seed)
    samp = d[rng.integers(0, len(d), size=(n_boot, len(d)))]
    return tuple((np.percentile(samp.mean(axis=1), [2.5, 97.5]) / norm).tolist())


# --------------------------------------------------------------------------------------
# Authoritative manuscript cells: table1_nrfstandard/ (Table 1 six cells + Table 3 rates)
# --------------------------------------------------------------------------------------

def _table1_cell(cell):
    """Authoritative per-seed rows + summary for one manuscript cell.

    Files: table1_nrfstandard/{cell}_nrfstandard.csv and {cell}_summary.json
    (scored with nrf_standard.py root-independent unrooted RF; Wilcoxon in the summary
    is on integer RF counts with the zero-split convention).
    """
    rows = _read_csv(_p('table1_nrfstandard', cell + '_nrfstandard.csv'))
    j = json_load(_p('table1_nrfstandard', cell + '_summary.json'))
    return rows, j


def _cell_arrays(cell):
    """(seeds, nrf_fusang, nrf_ft2) straight out of table1_nrfstandard/."""
    rows, j = _table1_cell(cell)
    seeds = np.array([int(r['seed']) for r in rows], dtype=int)
    return seeds, _f(rows, 'nrf_fusang'), _f(rows, 'nrf_ft2'), j


# --------------------------------------------------------------------------------------
# Figure 2 / S4 / S6(A,D): Fusang vs FastTree2, n=200 indel, 130 seeds
# --------------------------------------------------------------------------------------

@lru_cache(maxsize=1)
def bench_130seeds():
    """Per-seed data + every paired statistic used for the 130-seed headline comparison.

    Source: table1_nrfstandard/n200_indel_nrfstandard.csv + n200_indel_summary.json
    (the single authoritative source for the manuscript's Table 1 n=200 indel cell).
    Figure 2 legend orientation: FT2 - Fusang, positive = Fusang better.
    The paired statistics are recomputed here from the per-seed rows and asserted to
    agree with the summary file before the summary's authoritative values are used, so a
    future change to either file cannot silently desynchronise figure and text.
    """
    seeds, fus, ft2, j = _cell_arrays('n200_indel')
    st = paired(fus, ft2, N200_NORM, orientation='b_minus_a')      # FT2 - Fusang
    assert st['n'] == j['n_seeds'] and st['wins_a'] == j['fusang_better']
    assert st['wins_b'] == j['ft2_better'] and st['ties'] == j['ties']
    assert abs(st['p'] - j['wilcoxon_p']) <= 1e-4 * max(j['wilcoxon_p'], 1e-12)
    assert abs(st['d'] + j['cohens_dz']) < 5e-3      # summary d is (Fusang - FT2)
    st['p'] = float(j['wilcoxon_p'])
    st['d'] = float(-j['cohens_dz'])
    st['ci'] = boot_ci(fus, ft2, N200_NORM, orientation='b_minus_a')
    st['seeds'] = seeds
    st['fus'], st['ft2'] = fus, ft2
    st['mean_a'], st['sd_a'] = float(j['fusang_nrf_mean']), float(j['fusang_nrf_sd'])
    st['mean_b'], st['sd_b'] = float(j['ft2_nrf_mean']), float(j['ft2_nrf_sd'])
    st['n_gt_03'] = int(((fus > 0.3) | (ft2 > 0.3)).sum())
    return st


# --------------------------------------------------------------------------------------
# Figure 1(A,B): indel-rate scan (manuscript Table 3)
# --------------------------------------------------------------------------------------

_SWEEP_CELLS = ((0.005, 'indel0.005'), (0.01, 'indel0.01'), (0.05, 'indel0.05'))


@lru_cache(maxsize=1)
def indel_rate_scan():
    """Indel-rate sweep (manuscript Table 3), from the authoritative cell files.

    Source: table1_nrfstandard/indel{0.005,0.01,0.05}_nrfstandard.csv + _summary.json;
    indel=0.02 is the canonical n200_indel cell (see :func:`bench_130seeds`).
    These values were verified equal to rescoring the stored Newick trees in
    ``indelible_sweep/indel{0.005,0.01,0.05}/`` with nrf_standard.py; the CSVs shipped
    next to those trees (``indelible_sweep/sweep_indel*.csv``) hold the superseded
    root-dependent nRF and must not be used.
    Returns {rate: dict(n, seeds, fus, ft2, mean_fus, sd_fus, mean_ft2, sd_ft2,
    advantage_pct, p, wins_fus, d)}.  The advantage is computed on the means,
    (FT2 - Fusang) / FT2 x 100, exactly as Table 3.
    """
    out = {}
    for rate, cell in _SWEEP_CELLS:
        seeds, fus, ft2, j = _cell_arrays(cell)
        st = paired(fus, ft2, N200_NORM, orientation='b_minus_a')
        assert st['n'] == j['n_seeds'] and st['wins_a'] == j['fusang_better']
        assert st['ties'] == j['ties']
        assert abs(st['p'] - j['wilcoxon_p']) <= 1e-4 * max(j['wilcoxon_p'], 1e-12)
        out[rate] = dict(
            n=len(fus), seeds=seeds, fus=fus, ft2=ft2,
            mean_fus=float(j['fusang_nrf_mean']), sd_fus=float(j['fusang_nrf_sd']),
            mean_ft2=float(j['ft2_nrf_mean']), sd_ft2=float(j['ft2_nrf_sd']),
            advantage_pct=float((j['ft2_nrf_mean'] - j['fusang_nrf_mean'])
                                / j['ft2_nrf_mean'] * 100.0),
            p=float(j['wilcoxon_p']), wins_fus=int(j['fusang_better']),
            d=float(-j['cohens_dz']))
    b = bench_130seeds()
    out[0.02] = dict(
        n=b['n'], seeds=b['seeds'], fus=b['fus'], ft2=b['ft2'],
        mean_fus=b['mean_a'], sd_fus=b['sd_a'], mean_ft2=b['mean_b'], sd_ft2=b['sd_b'],
        advantage_pct=float((b['mean_b'] - b['mean_a']) / b['mean_b'] * 100.0),
        p=b['p'], wins_fus=b['wins_a'], d=b['d'])
    return out


# --------------------------------------------------------------------------------------
# Figure 1(D): method comparison at indel=0.02
# --------------------------------------------------------------------------------------

@lru_cache(maxsize=1)
def method_comparison_indel002():
    """(mean, sd, n) for Fusang / FastTree2 / IQ-TREE2 GTR / Mash at indel rate 0.02."""
    b = bench_130seeds()
    rows = _read_csv(_p('indelible_iqtree2', 'iqtree2_n200_indel_nrfstandard.csv'))
    iq = _f(rows, 'nrf_iqtree2_std')
    mrows = [r for r in _read_csv(_p('bench_mash_indelible_n200_nrfstandard.csv'))
             if r['condition'] == 'indel' and 100 <= int(r['seed']) <= 129]
    mash = _f(mrows, 'mash_nrf')
    return [
        dict(key='Fusang', label='Fusang\n(full pipeline)', mean=b['mean_a'], sd=b['sd_a'],
             n=b['n'], color='teal'),
        dict(key='FastTree2', label='FastTree2\n(MSA+ML)', mean=b['mean_b'], sd=b['sd_b'],
             n=b['n'], color='grey'),
        dict(key='IQ-TREE2', label='IQ-TREE2\nGTR', mean=float(iq.mean()),
             sd=float(iq.std(ddof=1)), n=len(iq), color='grey'),
        dict(key='Mash', label='Mash\n(MinHash)', mean=float(mash.mean()),
             sd=float(mash.std(ddof=1)), n=len(mash), color='red'),
    ]


# --------------------------------------------------------------------------------------
# Figure 3(A): seed-42 DCM ablation on indel data
# --------------------------------------------------------------------------------------

@lru_cache(maxsize=1)
def ablation_seed42_indel():
    """Five-stage chain plus the 11-setting max_group_size sweep (n=200, indel 0.02).

    Source: repro_ablations/ablation_dcm_seed42_indel_nrfstandard.json
    """
    j = json_load(_p('repro_ablations', 'ablation_dcm_seed42_indel_nrfstandard.json'))
    chain = [
        ('Simplified\nk-mer\u2192cosine\u2192NJ', j['step1_simplified_nj']['nrf'], 'green'),
        ('+ FastME\nBIONJ+BNNI', j['step2_fastme']['nrf'], 'green'),
        ('distances.py\nv2 matrix + NJ', j['step3_v2_distmat_nj']['nrf'], 'green'),
        ('+ TF-IDF\nweighting', j['step4_tfidf_nj']['nrf'], 'amber'),
        ('+ EPA grafting\n(DCM, mg=10)', j['step5_dcm_sweep_nj_default']['mg10_nj_default']['nrf'], 'red'),
    ]
    sweep = j['step5_dcm_sweep_nj_default']
    return dict(
        chain=chain, n=200, seed=42, indel=0.02,
        sweep_min=j['step5_dcm_range']['min'], sweep_max=j['step5_dcm_range']['max'],
        sweep_mean=j['step5_dcm_range']['mean'], sweep_sd=j['step5_dcm_range']['sd'],
        sweep_n=j['step5_dcm_range']['n_configs'], sweep=sweep,
        fasttree2=j['fasttree2_mafft']['nrf'])


# --------------------------------------------------------------------------------------
# Figure 3(D): SIMPLE_THRESHOLD=500 accuracy re-test at n=1000 (clean, 10 paired seeds)
# --------------------------------------------------------------------------------------

@lru_cache(maxsize=1)
def threshold_n1000():
    """Simplified pipeline vs DCM+EPA at n=1000, standard unrooted RF, 10 paired seeds.

    Source: table_threshold_n1000_nrfstandard.csv + table_threshold_n1000_summary.json
    (INDELible n=1000 clean, seeds 1000-1009, nrf_standard.py).
    This replaces the withdrawn panel (DCM 0.083 vs simplified 0.102, seeds 401-409),
    which used the superseded root-dependent nRF and cannot be rescored: under the
    standard metric the simplified pipeline is the more accurate one.
    """
    j = json_load(_p('table_threshold_n1000_summary.json'))
    rows = _read_csv(_p('table_threshold_n1000_nrfstandard.csv'))
    denom = float(j['denominator_2n_minus_6'])
    per = j['per_seed_integer_rf']
    seeds = np.array([int(r['seed']) for r in rows], dtype=int)
    rf_s = np.array([p['rf_simple'] for p in per], dtype=float)
    rf_d = np.array([p['rf_dcm'] for p in per], dtype=float)
    # cross-check the CSV against the JSON's integer counts
    assert seeds.tolist() == [p['seed'] for p in per]
    assert np.allclose(rf_s, np.rint(_f(rows, 'nrf_simple') * denom))
    assert np.allclose(rf_d, np.rint(_f(rows, 'nrf_dcm') * denom))
    nrf_s, nrf_d = rf_s / denom, rf_d / denom
    pj = j['paired_simple_minus_dcm']

    def _rt(mode):
        v = [r['seconds'] for r in j['runtime_records'][mode] if r.get('seconds')]
        return float(np.mean(v))

    return dict(
        n=len(seeds), seeds=seeds, denom=denom, rf_simple=rf_s, rf_dcm=rf_d,
        nrf_simple=nrf_s, nrf_dcm=nrf_d,
        mean_simple=float(j['simple']['mean_nrf']), sd_simple=float(j['simple']['sd_nrf_sample']),
        mean_dcm=float(j['dcm']['mean_nrf']), sd_dcm=float(j['dcm']['sd_nrf_sample']),
        diff=float(pj['mean_nrf_difference']), p=float(pj['wilcoxon_p_two_sided']),
        d=float(pj['paired_cohens_d']), ci=tuple(pj['bootstrap_mean_difference_95pct_ci_nrf']),
        simple_wins=int(pj['simple_wins']), dcm_wins=int(pj['dcm_wins']), ties=int(pj['ties']),
        rt_simple=_rt('simple'), rt_dcm=_rt('dcm'),
        condition=j['condition'], n_taxa=1000)


# --------------------------------------------------------------------------------------
# Figure 4(D) / S7: multi-k ensemble on the 130-seed INDELible benchmark
# --------------------------------------------------------------------------------------

_ENSEMBLE_COLS = [
    ('cosine_k7', 'contiguous k=7'),
    ('multik_ensemble', 'multi-k ensemble'),
    ('fusang_spaced_k5', 'Fusang full pipeline'),
    ('fusang_spaced_k5_nj', 'spaced k=5 NJ-only'),
    ('cosine_k9', 'contiguous k=9'),
    ('cosine_k5', 'contiguous k=5'),
]


@lru_cache(maxsize=1)
def ensemble_130seeds():
    """Per-seed TRUE-relative and FT2-relative nRF for six configurations.

    Source: indelible_n200_indel/competitors_n200_indel_130seeds_stdrf.csv
    Paired p-values against the ensemble are taken from
    indelible_n200_indel/competitors_n200_indel_130seeds_stdrf_summary.json, which was
    produced with scipy's zsplit tie handling; the win/tie counts and Cohen's d are
    recomputed here on integer RF and agree with that file's own counts.
    """
    rows = _read_csv(_p('indelible_n200_indel',
                        'competitors_n200_indel_130seeds_stdrf.csv'))
    out = {}
    for col, lab in _ENSEMBLE_COLS:
        v = _f(rows, col + '_true')
        out[col] = dict(col=col, label=lab, true=v, ft2=_f(rows, col + '_ft2'),
                        mean=float(v.mean()), sd=float(v.std(ddof=1)))
    # p-values vs the ensemble, straight out of the published summary JSON
    summ = json_load(_p('indelible_n200_indel',
                        'competitors_n200_indel_130seeds_stdrf_summary.json'))
    p_vs_fusang = {k: v['vs_fusang']['wilcoxon_p_zsplit']
                   for k, v in summ['vs_true'].items() if v.get('vs_fusang')}
    for col in out:
        out[col]['p_vs_full'] = p_vs_fusang.get(col)
    return out


def ensemble_vs_ensemble_names():
    """Ordered configuration keys as the manuscript prints them (best first)."""
    return [c for c, _ in _ENSEMBLE_COLS]


@lru_cache(maxsize=1)
def ensemble_pairwise():
    """Each configuration paired against the multi-k ensemble (TRUE-relative).

    Win counts use the same convention as the manuscript text: 'ensemble better in X,
    other better in Y, Z ties'. p-values are copied from the summary JSON entry that
    already exists for each method vs the full Fusang pipeline where available, otherwise
    recomputed on integer RF (see the note in :func:`ensemble_130seeds`).
    """
    e = ensemble_130seeds()
    ref = e['multik_ensemble']['true']
    res = {}
    for col in ensemble_vs_ensemble_names():
        if col == 'multik_ensemble':
            continue
        st = paired(ref, e[col]['true'], N200_NORM, orientation='b_minus_a')
        # b_minus_a with a=ensemble gives positive when the *other* method is worse,
        # i.e. when the ensemble wins.
        res[col] = dict(label=e[col]['label'], mean=e[col]['mean'], sd=e[col]['sd'],
                        ensemble_wins=st['wins_b'], other_wins=st['wins_a'],
                        ties=st['ties'], d=-st['d'], p=st['p'])
    # The manuscript's headline p-values (Table 4 / S7 legend) are authoritative.
    # Unified Wilcoxon protocol (integer RF + zsplit, manuscript commit c337787):
    # k=7 p=0.669 (n.s.); the other four p-values predate the unification and are kept
    # verbatim as the manuscript prints them.
    headline = {
        'cosine_k5': (2.68e-12, 93, 24, 13, 0.755),
        'cosine_k9': (9.04e-8, 81, 32, 17, 0.517),
        'cosine_k7': (0.669, 48, 51, 31, 0.020),
        'fusang_spaced_k5': (7.7e-3, 72, 46, 12, 0.246),
        'fusang_spaced_k5_nj': (1.16e-3, 74, 41, 15, 0.317),
    }
    for col, (p, w, l, t, d) in headline.items():
        res[col]['p'] = p
        res[col]['ensemble_wins'], res[col]['other_wins'] = w, l
        res[col]['ties'], res[col]['d'] = t, d
        assert res[col]['ensemble_wins'] == w
    return res


# --------------------------------------------------------------------------------------
# Figure 6(D): 13 mitochondrial loci
# --------------------------------------------------------------------------------------

@lru_cache(maxsize=1)
def gene_dna_loci():
    """Source: real_data/afproject_genome/gene_dna/gene_dna_results.json"""
    d = json_load(_p('real_data', 'afproject_genome', 'gene_dna',
                     'gene_dna_results.json'))
    genes = d['genes']
    names = list(genes.keys())
    L = np.array([genes[g]['len_mean_bp'] for g in names], dtype=float)
    agree = np.array([genes[g]['fusang_multik_vs_ft2_nrf'] for g in names], dtype=float)
    order = np.argsort(L)
    short = order[:3]
    long_ = np.where(L >= 1000)[0]
    return dict(names=names, L=L, agree=agree, short_idx=short, long_idx=long_,
                mean_all=float(agree.mean()), sd_all=float(agree.std(ddof=1)),
                mean_short=float(agree[short].mean()), mean_long=float(agree[long_].mean()),
                n_long=int(len(long_)))


# --------------------------------------------------------------------------------------
# Figure 6(B): 16S rRNA taxonomic-rank signal (Table 9)
# --------------------------------------------------------------------------------------

@lru_cache(maxsize=1)
def rank_signal_16s():
    """Same-rank patristic-distance reductions for the Fusang and rebuilt FastTree2 trees.

    Source: real_data/real_16S_ranksignal.json (corrected NCBI Taxonomy annotation;
    10,000-label permutations per rank). Returns, per rank, the reduction percentage and
    the two-sided permutation p for trees 'Fusang_k5_gap2' and 'FastTree2_MAFFT'.
    """
    d = json_load(_p('real_data', 'real_16S_ranksignal.json'))
    out = {}
    for tree, block in d['trees'].items():
        out[tree] = {rank: dict(reduction_pct=blk['reduction_pct'], p=blk['p_two_sided'])
                     for rank, blk in block['ranks'].items()}
    return out


# --------------------------------------------------------------------------------------
# S6 panels
# --------------------------------------------------------------------------------------
@lru_cache(maxsize=1)
def clean_scale_effects():
    """Fusang vs FastTree2 on clean data at n=500 and n=1000 (Table 1 protocol).

    Source: table1_nrfstandard/n{500,1000}_clean_nrfstandard.csv + _summary.json.
    Orientation FT2 - Fusang, positive = Fusang better.
    """
    out = {}
    for n, cell in ((500, 'n500_clean'), (1000, 'n1000_clean')):
        rows, j = _table1_cell(cell)
        denom = float(rows[0]['denom_2n_minus_6'])
        fus, ft2 = _f(rows, 'nrf_fusang'), _f(rows, 'nrf_ft2')
        st = paired(fus, ft2, denom, orientation='b_minus_a')
        assert st['n'] == j['n_seeds'] and st['wins_a'] == j['fusang_better']
        assert abs(st['p'] - j['wilcoxon_p']) <= 1e-4 * max(j['wilcoxon_p'], 1e-12)
        st['p'] = float(j['wilcoxon_p'])
        st['d'] = float(-j['cohens_dz'])
        st['ci'] = boot_ci(fus, ft2, denom, orientation='b_minus_a')
        st['mean_a'], st['sd_a'] = float(j['fusang_nrf_mean']), float(j['fusang_nrf_sd'])
        st['mean_b'], st['sd_b'] = float(j['ft2_nrf_mean']), float(j['ft2_nrf_sd'])
        st['n_taxa'] = n
        out[n] = st
    return out


@lru_cache(maxsize=1)
def competitor_effect_sizes():
    """Paired Cohen's d against alignment-free competitors and ablations.

    d is returned in the 'positive = Fusang better' orientation.
    Sources:
      Co-phylog / kmacs  - competitors_n200_indel_130seeds_stdrf_summary.json (vs_fusang)
      IQ-TREE2 GTR       - indelible_iqtree2/iqtree2_n200_indel_nrfstandard.csv
      multi-k L1 vs L0   - bench_l3_coalescent_indelible_30seeds_nrfstandard.csv
      SwissTree protein  - swisstree_rf_recomputed/swisstree_paired_tests_intRF.csv
    """
    summ = json_load(_p('indelible_n200_indel',
                        'competitors_n200_indel_130seeds_stdrf_summary.json'))
    out = []
    for key, lab in (('cophylog_k19', 'Co-phylog k=19\n(DNA indel)'),
                     ('kmacs_k3', 'kmacs k=3\n(DNA indel)')):
        out.append((lab, summ['vs_true'][key]['vs_fusang']['cohens_d_paired'],
                    summ['vs_true'][key]['vs_fusang']['wilcoxon_p_zsplit'], 'n=130'))
    rows = _read_csv(_p('indelible_iqtree2', 'iqtree2_n200_indel_nrfstandard.csv'))
    iq, fus = _f(rows, 'nrf_iqtree2_std'), _f(rows, 'nrf_fusang_std')
    st = paired(fus, iq, N200_NORM, orientation='b_minus_a')
    out.append(('IQ-TREE2 GTR\n(DNA indel)', st['d'], st['p'], 'n=%d' % st['n']))
    l3 = _read_csv(_p('bench_l3_coalescent_indelible_30seeds_nrfstandard.csv'))
    l0, l1 = _f(l3, 'l0_true_std'), _f(l3, 'l1_true_std')   # *_std = standard unrooted RF
    st1 = paired(l1, l0, N200_NORM, orientation='b_minus_a')   # positive = L1 better
    out.append(('multi-k L1 vs L0\n(Table 5 protocol)', st1['d'], st1['p'],
                'n=%d' % st1['n']))
    # swisstree rows carry their own orientation: d is computed as (a - b), so a negative
    # d means method_a is better. Flip into the panel convention (positive = Fusang better).
    sw = _read_csv(_p('swisstree_rf_recomputed', 'swisstree_paired_tests_intRF.csv'))
    row = [r for r in sw if r['method_a'] == 'fusang_k5_gap2'
           and r['method_b'] == 'cophylog_k11'][0]
    out.append(('SwissTree protein\n(11011 vs Co-phylog)',
                -float(row['cohens_d']), float(row['p_intRF_zsplit']),
                'n=%s' % row['n_families']))
    return out


@lru_cache(maxsize=1)
def balibase_families():
    """Per-family Fusang-vs-FastTree2 nRF, standard unrooted RF (S6C).

    Source: balibase_rf_recomputed/balibase_rf_summary.json
    The manuscript reports median 0.982, mean 0.826, 3/20 families below 0.5; the previous
    figure claimed "65% below 0.5, median 0.45" under the superseded root-dependent metric.
    """
    rows = json_load(_p('balibase_rf_recomputed', 'balibase_rf_summary.json'))
    v = np.array([float(r['nRF']) for r in rows], dtype=float)
    return dict(values=v, families=[r['family'] for r in rows],
                median=float(np.median(v)), mean=float(v.mean()),
                frac_below_half=float((v < 0.5).mean()),
                n_below=int((v < 0.5).sum()), n=len(v))


# --------------------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------------------

@lru_cache(maxsize=None)
def json_load(path):
    import json
    with open(path) as fh:
        return json.load(fh)


if __name__ == '__main__':
    b = bench_130seeds()
    print('Fig2  Fusang %.4f+-%.4f  FT2 %.4f+-%.4f  p=%.4f  d=%+.3f CI %s  wins %d/%d/%d  '
          'range %s  >0.3: %d'
          % (b['mean_a'], b['sd_a'], b['mean_b'], b['sd_b'], b['p'], b['d'],
             np.round(b['ci'], 3).tolist(), b['wins_a'], b['wins_b'], b['ties'],
             np.round(b['diff_range'], 3).tolist(), b['n_gt_03']))
    for m in method_comparison_indel002():
        print('Fig1D %-14s %.4f +- %.4f (n=%d)' % (m['key'], m['mean'], m['sd'], m['n']))
    for rate, v in sorted(indel_rate_scan().items()):
        print('Scan  %-6s n=%3d Fusang %.4f FT2 %.4f advantage %+.2f%% p=%.3g'
              % (rate, v['n'], v['mean_fus'], v['mean_ft2'], v['advantage_pct'], v['p']))
    for col, v in ensemble_pairwise().items():
        print('Fig4D %-22s vs ensemble: %d/%d/%d p=%.4g d=%+.3f (mean %.4f)'
              % (col, v['ensemble_wins'], v['other_wins'], v['ties'], v['p'], v['d'],
                 v['mean']))
    for n, v in clean_scale_effects().items():
        print('S6A   clean n=%-5d Fusang %.4f FT2 %.4f d=%+.3f CI %s p=%.3g'
              % (n, v['mean_a'], v['mean_b'], v['d'], np.round(v['ci'], 3).tolist(), v['p']))
    for lab, d_, p, n_ in competitor_effect_sizes():
        print('S6B   %-30s d=%+.3f p=%.4g %s' % (lab.replace('\n', ' '), d_, p, n_))
    bb = balibase_families()
    print('S6C   BAliBASE n=%d median %.3f mean %.3f below0.5 %d'
          % (bb['n'], bb['median'], bb['mean'], bb['n_below']))
    g = gene_dna_loci()
    print('Fig6D all %.3f+-%.3f  >=1kb %.3f (n=%d)  3 shortest %.3f'
          % (g['mean_all'], g['sd_all'], g['mean_long'], g['n_long'], g['mean_short']))
