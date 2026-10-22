#!/usr/bin/env python3
"""Authoritative standard-RF rescoring for manuscript Table 1.

Table 1 reports Fusang vs FastTree2 nRF at n=200/500/1000 under clean and
indel conditions. The published values were produced with the standard
root-independent unrooted RF (nrf_standard.py) but no per-seed CSV was ever
written to disk; the summary_*.json files sitting beside the trees still hold
values from the SUPERSEDED root-dependent implementation (calc_nrf_simple.py)
and must not be quoted.

This script recomputes every Table 1 cell from the archived trees, writes a
per-seed CSV plus a summary JSON per cell, and prints a comparison against the
values printed in the manuscript so any drift is caught immediately.

Scoring: nrf_standard.nrf (|S1\\S2|+|S2\\S1| over non-trivial splits, divided
by 2(n-3); each split canonicalised to the side not containing the
lexicographically smallest leaf).

Statistics: paired Wilcoxon is run on INTEGER RF counts (nRF x 2(n-3)) with
zero_method='zsplit'; feeding float nRF breaks rank ties via ~1e-17 rounding
noise and shifts p-values.
"""
import glob
import json
import os
import re
import statistics as st

import numpy as np
from scipy.stats import wilcoxon

import nrf_standard as N

# label -> (directory, condition)
CELLS = [
    ('n200_clean', 'indelible_n200_clean', 'clean'),
    ('n200_indel', 'indelible_n200_indel', 'indel'),
    ('n500_clean', 'indelible_n500_clean', 'clean'),
    ('n500_indel', 'indelible_n500_indel', 'indel'),
    ('n1000_clean', 'indelible_n1000_clean', 'clean'),
    ('n1000_indel', 'indelible_n1000_indel', 'indel'),
]

# Table 3 indel-rate scan (n=200, L=500 bp, 130 seeds each).
SWEEP = [
    ('indel0.005', 'indelible_sweep/indel0.005', 'indel'),
    ('indel0.01', 'indelible_sweep/indel0.01', 'indel'),
    ('indel0.05', 'indelible_sweep/indel0.05', 'indel'),
]

# Values currently printed in NAR_MANUSCRIPT_REVISED.md Table 1, for drift check.
MANUSCRIPT = {
    'n200_clean': (0.234, 0.210, -11.6),
    'n200_indel': (0.199, 0.205, 2.8),
    'n500_clean': (0.244, 0.213, -14.7),
    'n500_indel': (0.208, 0.209, 0.6),
    'n1000_clean': (0.257, 0.212, -21.4),
    'n1000_indel': (0.226, 0.212, -6.9),
    # Table 3: (fusang, ft2, relative advantage %)
    'indel0.005': (0.224, 0.206, -8.9),
    'indel0.01': (0.215, 0.205, -5.1),
    'indel0.05': (0.180, 0.202, 10.9),
}

CELLS = CELLS + SWEEP

OUT_DIR = 'table1_nrfstandard'


def collect(directory, cond):
    """Return list of (seed, n_leaves, rf_fusang, rf_ft2, denom)."""
    rows = []
    for tp in sorted(glob.glob(os.path.join(directory, 'seed*_%s_true.nwk' % cond))):
        m = re.search(r'seed(\d+)_', os.path.basename(tp))
        if not m:
            continue
        s = m.group(1)
        fp = os.path.join(directory, 'seed%s_%s_fusang.nwk' % (s, cond))
        kp = os.path.join(directory, 'seed%s_%s_ft2.nwk' % (s, cond))
        if not (os.path.exists(fp) and os.path.exists(kp)):
            continue
        true = open(tp).read().strip()
        fus = open(fp).read().strip()
        ft2 = open(kp).read().strip()
        n = len(N.leaves(N.parse_newick(true)))
        denom = 2 * (n - 3)
        rf_f, rf_t = N.rf_distance(fus, true), N.rf_distance(ft2, true)
        rows.append((int(s), n, rf_f, rf_t, denom))
    return sorted(rows)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    report = {}
    for label, d, cond in CELLS:
        rows = collect(d, cond)
        if not rows:
            print('%-12s NO DATA (dir=%s)' % (label, d))
            continue
        nf = np.array([r[2] / r[4] for r in rows])
        nt = np.array([r[3] / r[4] for r in rows])
        rf_f = np.array([r[2] for r in rows])
        rf_t = np.array([r[3] for r in rows])

        # paired Wilcoxon on INTEGER RF counts
        diff_int = rf_f - rf_t
        nz = diff_int[diff_int != 0]
        if len(nz) >= 5:
            W, p = wilcoxon(diff_int, zero_method='zsplit')
        else:
            W, p = float('nan'), float('nan')

        # Exact sign-flip permutation p, for the degenerate case where every
        # seed points the same way.  There the signed-rank statistic is pinned
        # at its floor (W = 0, or W = max), only one of the 2**k sign
        # assignments realises it, and the large-sample normal approximation
        # badly overstates the p-value (it returns ~1.73e-6 where the exact
        # value is 2/2**30 = 1.86e-9).  The |d| tie structure does not affect
        # the count: W = 0 requires every sign to be negative, which is one
        # assignment out of 2**k whatever the ranks are.
        unanimous = bool(((nz > 0).all() or (nz < 0).all())) if len(nz) else False
        perm_p = float(2.0 / (2 ** len(nz))) if unanimous else None
        # The p-value the manuscript reports for this cell.
        reported_p = perm_p if perm_p is not None else float(p)
        dz = float(np.mean(diff_int) / np.std(diff_int, ddof=1)) \
            if np.std(diff_int, ddof=1) > 0 else float('nan')

        # bootstrap CI on the paired mean nRF difference (Fusang - FT2)
        d_nrf = nf - nt
        rng = np.random.default_rng(20261022)
        boot = np.array([rng.choice(d_nrf, size=len(d_nrf), replace=True).mean()
                         for _ in range(10000)])
        ci = [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))]

        summary = {
            'label': label,
            'directory': d,
            'condition': cond,
            'n_seeds': len(rows),
            'seed_range': [rows[0][0], rows[-1][0]],
            'n_taxa': rows[0][1],
            'fusang_nrf_mean': float(nf.mean()),
            'fusang_nrf_sd': float(nf.std(ddof=1)),
            'ft2_nrf_mean': float(nt.mean()),
            'ft2_nrf_sd': float(nt.std(ddof=1)),
            'mean_diff_fusang_minus_ft2': float(d_nrf.mean()),
            'bootstrap_ci95_mean_diff': ci,
            'wilcoxon_on_integer_rf': True,
            'wilcoxon_W': float(W),
            'wilcoxon_p_normal_approx': float(p),
            'permutation_p_exact': perm_p,
            'unanimous_all_seeds': unanimous,
            'wilcoxon_p': reported_p,
            'cohens_dz': dz,
            'fusang_better': int((diff_int < 0).sum()),
            'ft2_better': int((diff_int > 0).sum()),
            'ties': int((diff_int == 0).sum()),
            'metric': 'nrf_standard (root-independent unrooted RF)',
        }
        report[label] = summary

        # computed here (not after the JSON write) so it lands in the per-cell
        # summary file as well as the aggregate one
        rel = (nt.mean() - nf.mean()) / nt.mean() * 100.0
        summary['relative_advantage_pct'] = float(rel)

        csv_path = os.path.join(OUT_DIR, '%s_nrfstandard.csv' % label)
        with open(csv_path, 'w', encoding='utf-8') as fh:
            fh.write('seed,n_taxa,rf_fusang,rf_ft2,denom_2n_minus_6,'
                     'nrf_fusang,nrf_ft2,rf_diff_fusang_minus_ft2\n')
            for s, n, rf_f_, rf_t_, den in rows:
                fh.write('%d,%d,%d,%d,%d,%.6f,%.6f,%d\n'
                         % (s, n, rf_f_, rf_t_, den,
                            rf_f_ / den, rf_t_ / den, rf_f_ - rf_t_))
        with open(os.path.join(OUT_DIR, '%s_summary.json' % label), 'w',
                  encoding='utf-8') as fh:
            json.dump(summary, fh, indent=1)

        ms = MANUSCRIPT.get(label)
        flag = ''
        if ms:
            ok = (abs(nf.mean() - ms[0]) < 0.001 and abs(nt.mean() - ms[1]) < 0.001
                  and abs(rel - ms[2]) < 0.15)
            flag = ('  [MATCH]' if ok else
                    '  [DRIFT vs manuscript %.3f/%.3f/%.1f%%]' % (ms[0], ms[1], ms[2]))
        # print the value the MANUSCRIPT reports, not the raw normal
        # approximation, so a reviewer running this script sees the numbers
        # they are trying to check
        tag = ' (exact perm, W=0)' if unanimous else ''
        print('%-12s n=%3d  Fusang %.4f+-%.4f  FT2 %.4f+-%.4f  p=%.3g%s  '
              'd=%+.2f  wins %d/%d/%d%s'
              % (label, len(rows), nf.mean(), nf.std(ddof=1),
                 nt.mean(), nt.std(ddof=1), reported_p, tag, dz,
                 summary['fusang_better'], summary['ft2_better'],
                 summary['ties'], flag))

    with open(os.path.join(OUT_DIR, 'TABLE1_all_summary.json'), 'w',
              encoding='utf-8') as fh:
        json.dump(report, fh, indent=1)
    print('\nWrote per-cell CSV + summary JSON to %s/' % OUT_DIR)


if __name__ == '__main__':
    main()
