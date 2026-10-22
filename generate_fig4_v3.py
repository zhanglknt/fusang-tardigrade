"""Figure 4 (v3): Spaced k-mer parameter optimization — 4 panels.

A: k x gap heatmap at n=200 (data/S1_grid_search_results.csv)
B: adaptive parameter selection logic (schematic)
C: 5-repeat stability validation summary (Supplementary Table S5)
D: multi-k ensemble vs single-k configurations and the full pipeline (130 seeds)

Panel A note: ``data/S1_grid_search_results.csv`` was produced by the superseded
root-dependent nRF scorer and cannot be rescored (the per-configuration trees are not
retained); the panel is therefore drawn as a shape illustration and explicitly labelled.
Panel D reads the standard unrooted RF 130-seed benchmark via ``figure_data``:
``indelible_n200_indel/competitors_n200_indel_130seeds_stdrf.csv`` plus the paired
p-values/wins recorded in ``competitors_n200_indel_130seeds_stdrf_summary.json`` and
printed in manuscript Table 4.
"""
import os
os.environ['MPLCONFIGDIR'] = 'matplotlib_cache'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import csv
from matplotlib.patches import FancyBboxPatch

import figure_data as FD

INK, GREY = '#1a2332', '#64748b'
BLUE, TEAL, AMBER, RED, GREEN = '#2563eb', '#0d9488', '#d97706', '#dc2626', '#16a34a'

fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.6))
fig.patch.set_facecolor('white')

def panel_label(ax, s):
    ax.text(-0.14, 1.06, s, transform=ax.transAxes, fontsize=14,
            fontweight='bold', color=INK, va='top')

# ---------- A: heatmap k x gap at n=200 ----------
ax = axes[0][0]
rows = [r for r in csv.DictReader(open('data/S1_grid_search_results.csv'))
        if int(r['n_taxa']) == 200]
ks = sorted(set(int(r['k']) for r in rows))
gs = ['none', 'gap1', 'gap2', 'gap3']
M = np.full((len(gs), len(ks)), np.nan)
for r in rows:
    M[gs.index(r['gap_pattern']), ks.index(int(r['k']))] = float(r['nRF'])
im = ax.imshow(M, cmap='viridis_r', aspect='auto')
for i in range(len(gs)):
    for j in range(len(ks)):
        if not np.isnan(M[i, j]):
            ax.text(j, i, f'{M[i,j]:.3f}', ha='center', va='center', fontsize=8,
                    color='white' if M[i, j] > np.nanmedian(M) else INK)
ax.set_xticks(range(len(ks))); ax.set_xticklabels([f'k={k}' for k in ks], fontsize=9)
ax.set_yticks(range(len(gs)))
ax.set_yticklabels(['gap0\n(contig.)', 'gap1', 'gap2', 'gap3'], fontsize=9)
fig.colorbar(im, ax=ax, label='nRF (lower = better)', shrink=0.85)
# mark default
if 5 in ks:
    j = ks.index(5); i = gs.index('gap2')
    ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec=RED, lw=2.4))
ax.set_title('Grid search at n=200 (red box = default k=5,gap2)', fontsize=10, color=INK, pad=6)
ax.text(0.5, -0.16, 'superseded root-dependent nRF implementation\n(trees not retained; shown for shape only)',
        transform=ax.transAxes, fontsize=6.4, color=GREY, style='italic', ha='center')
panel_label(ax, 'A')

# ---------- B: adaptive logic schematic ----------
ax = axes[0][1]
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis('off')
ax.set_title('Adaptive k,gap selection by dataset size', fontsize=10.5, color=INK, pad=6)
panel_label(ax, 'B')
ax.add_patch(FancyBboxPatch((0.5, 5.6), 4.0, 2.6, boxstyle='round,pad=0.1,rounding_size=0.25',
                            fc=TEAL, ec='none'))
ax.text(2.5, 7.65, 'n \u2264 100', fontsize=11, color='white', ha='center',
        va='center', fontweight='bold')
ax.text(2.5, 6.35, 'k=4, gap1\nsmaller pattern for\nlimited taxon sampling', fontsize=9,
        color='white', ha='center', va='center')
ax.add_patch(FancyBboxPatch((5.5, 5.6), 4.0, 2.6, boxstyle='round,pad=0.1,rounding_size=0.25',
                            fc=BLUE, ec='none'))
ax.text(7.5, 7.65, 'n > 100', fontsize=11, color='white', ha='center',
        va='center', fontweight='bold')
ax.text(7.5, 6.35, 'k=5, gap2\nlonger pattern for\nricher k-mer statistics', fontsize=9,
        color='white', ha='center', va='center')
ax.add_patch(FancyBboxPatch((2.2, 1.6), 5.6, 2.6, boxstyle='round,pad=0.1,rounding_size=0.25',
                            fc='#f4f7fa', ec=AMBER, lw=2))
ax.text(5.0, 3.55, 'genome scale (>>1 kb): multi-k ensemble', fontsize=9.5, color=AMBER,
        ha='center', va='center', fontweight='bold')
ax.text(5.0, 2.45, 'k=5 saturates (1,024 possible 5-mers);\nensemble of k=5,7,9 rescues signal — no manual k choice',
        fontsize=8.5, color=INK, ha='center', va='center')

# ---------- C: adaptive selection step plot (5-repeat stability) ----------
ax = axes[1][0]
scales = [20, 50, 100, 200, 500]
sel_k   = [4, 4, 4, 5, 5]   # 5/5 repeats identical at every scale
sel_gap = [1, 1, 1, 2, 2]
ax.step(scales, sel_k, where='mid', color=TEAL, lw=2.4, label='selected k')
ax.plot(scales, sel_k, 'o', color=TEAL, ms=7, mfc='white', mew=2)
ax.step(scales, sel_gap, where='mid', color=AMBER, lw=2.4, ls='--', label='selected gap')
ax.plot(scales, sel_gap, 's', color=AMBER, ms=6.5, mfc='white', mew=2)
for n, k, g in zip(scales, sel_k, sel_gap):
    ax.annotate(f'k={k},gap{g}', (n, k), textcoords='offset points', xytext=(0, 8),
                ha='center', fontsize=7.5, color=INK)
ax.set_xscale('log'); ax.set_xticks(scales)
ax.set_xticklabels([str(s) for s in scales], fontsize=9)
ax.set_xlabel('Dataset size n (log scale)', fontsize=9.5, color=INK)
ax.set_ylabel('Selected parameter value', fontsize=9.5, color=INK)
ax.set_yticks([1, 2, 3, 4, 5]); ax.set_ylim(0.4, 5.9)
ax.legend(fontsize=8.5, frameon=False, loc='center right')
ax.set_title('Adaptive k,gap selection vs n\n5/5 repeats identical; nRF deviation < 0.005',
             fontsize=10, color=INK, pad=6)
ax.grid(alpha=0.25)
for s in ('top', 'right'): ax.spines[s].set_visible(False)
panel_label(ax, 'C')

# ---------- D: multi-k ensemble vs configurations (130 seeds) ----------
ax = axes[1][1]
pw = FD.ensemble_pairwise()
order = FD.ensemble_vs_ensemble_names()          # best-first ordering
e = FD.ensemble_130seeds()
means = [e[c]['mean'] for c in order]
sds = [e[c]['sd'] for c in order]
names = [e[c]['label'] for c in order]
n_seeds = len(e['multik_ensemble']['true'])

def _fmt_p(p):
    # three significant figures, matching how the manuscript prints these p-values
    return 'p=%.3g' % p

bar_cols = []
for c in order:
    bar_cols.append(GREEN if c == 'multik_ensemble' else
                    (TEAL if c == 'cosine_k7' else GREY))
bars = ax.bar(range(len(order)), means, yerr=sds, capsize=3, width=0.62,
              color=bar_cols, error_kw=dict(lw=1.0, ecolor=INK))
for i, c in enumerate(order):
    ax.text(i, means[i] + sds[i] + 0.010, f'{means[i]:.4f}', ha='center',
            fontsize=8.2, fontweight='bold', color=INK)
    if c != 'multik_ensemble':
        v = pw[c]
        # in-bar label: p (prefix stripped so the string fits inside the narrow bar)
        # over ensemble wins / 130
        note = f'{_fmt_p(v["p"])[2:]}\n{v["ensemble_wins"]}/130'
        ax.text(i, 0.010, note, ha='center', va='bottom', fontsize=6.8, color='white')
ax.text(order.index('multik_ensemble'), 0.010, 'reference\nmethod', ha='center',
        va='bottom', fontsize=7.0, color='white', style='italic')
ax.text(0.015, 0.985, 'labels in bars: Wilcoxon p (integer RF, zsplit)\n'
                      'over ensemble wins / 130 seeds',
        transform=ax.transAxes, fontsize=6.4, color=GREY, style='italic', va='top')
ax.set_xticks(range(len(order)))
ax.set_xticklabels(names, fontsize=7.4, rotation=24, ha='right')
ax.set_ylabel('nRF vs true tree (mean ± SD)', fontsize=9.5, color=INK)
ax.set_ylim(0, 0.30)
ax.set_title(f'Multi-k ensemble vs alternatives ({n_seeds} seeds, TRUE-relative)\n'
             'beats k=5, k=9 and the full pipeline; ties k=7 (p=0.669)',
             fontsize=9, color=INK, pad=6)
ax.spines[['top', 'right']].set_visible(False)
panel_label(ax, 'D')

fig.suptitle('Figure 4. Spaced k-mer parameter optimization', fontsize=14,
             fontweight='bold', color=INK, y=0.985)
fig.tight_layout(rect=[0, 0, 1, 0.955])
fig.savefig('Figure4.pdf')
fig.savefig('Figure4.png', dpi=300)
print('Figure4 done |', ' '.join('%s %.4f+-%.4f' % (c, e[c]['mean'], e[c]['sd'])
                                 for c in order))
