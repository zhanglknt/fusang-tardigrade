"""Supplementary Figure S4 (v3): 130-seed indel benchmark detail — 4 panels.

All data from ``figure_data.bench_130seeds()``
(``indelible_n200_indel/bench_n200_indel_130seeds_nrfstandard.csv``; standard root-
independent unrooted RF). All 130 seeds are valid; none exceeds the pre-specified
nRF>0.3 exclusion threshold. The previous ``table3_corrected.csv`` input (120 seeds,
superseded root-dependent metric, 8 catastrophic exclusions) is no longer used.

A: violins for all 130 paired seeds
B: per-seed scatter Fusang vs FT2 with the pre-specified exclusion threshold
C: sorted paired differences (Fusang − FT2; negative = Fusang better)
D: histogram of paired differences with Wilcoxon p, paired Cohen's d and bootstrap CI
"""
import os
os.environ['MPLCONFIGDIR'] = 'matplotlib_cache'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import figure_data as FD


def trunc2(x):
    """Truncate toward zero at 2 decimals, matching the manuscript's CI rounding."""
    import math
    return (math.floor(x * 100) if x >= 0 else math.ceil(x * 100)) / 100

INK, GREY = '#1a2332', '#64748b'
TEAL, AMBER, RED = '#0d9488', '#d97706', '#dc2626'

st = FD.bench_130seeds()
fus, ft2, seeds = st['fus'], st['ft2'], st['seeds']
diffs = fus - ft2  # S4 orientation: Fusang − FastTree2; negative = Fusang better
n = len(diffs)
d_fs = -st['d']                       # paired Cohen's d in the S4 orientation
ci_fs = (-st['ci'][1], -st['ci'][0])  # same CI, flipped orientation

fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.6))
fig.patch.set_facecolor('white')

def panel_label(ax, s):
    ax.text(-0.14, 1.06, s, transform=ax.transAxes, fontsize=14,
            fontweight='bold', color=INK, va='top')

# ---------- A: violins ----------
ax = axes[0][0]
vp = ax.violinplot([fus, ft2], positions=[1, 2], widths=0.62, showmeans=True)
for b, c in zip(vp['bodies'], [TEAL, GREY]):
    b.set_facecolor(c); b.set_alpha(0.65)
for key in ('cbars', 'cmins', 'cmaxes', 'cmeans'):
    vp[key].set_color(INK); vp[key].set_linewidth(1.1)
ax.set_xticks([1, 2]); ax.set_xticklabels(['Fusang', 'FastTree2'], fontsize=10)
ax.set_ylabel('nRF vs true tree', fontsize=9.5, color=INK)
ax.set_title(f'Distributions, all {n} paired seeds\n'
             f'(none exceeds the nRF>0.3 threshold: {st["n_gt_03"]} excluded)',
             fontsize=9.5, color=INK, pad=6)
ax.grid(alpha=0.25, axis='y')
ax.spines[['top', 'right']].set_visible(False)
panel_label(ax, 'A')

# ---------- B: per-seed scatter with exclusion rule ----------
ax = axes[0][1]
ax.scatter(fus, ft2, s=13, color=TEAL, alpha=0.6, edgecolor='none', zorder=3)
lims = (0, max(fus.max(), ft2.max()) * 1.12)
ax.plot(lims, lims, color=GREY, lw=1, ls='--', zorder=1)
ax.axvline(0.3, color=RED, lw=1.0, ls=':', zorder=2)
ax.axhline(0.3, color=RED, lw=1.0, ls=':', zorder=2)
ax.text(0.298, lims[1] * 0.02, 'nRF=0.3\nexclusion rule', fontsize=7.5, color=RED,
        va='bottom', ha='right')
ax.text(lims[1] * 0.03, lims[1] * 0.86,
        f'all {n} seeds below the threshold\n(no outlier exclusion applied)',
        fontsize=8, color=INK)
ax.set_xlim(lims); ax.set_ylim(lims)
ax.set_xlabel('Fusang nRF', fontsize=9.5, color=INK)
ax.set_ylabel('FastTree2 nRF', fontsize=9.5, color=INK)
ax.set_title('Per-seed Fusang vs FastTree2 (dashed: parity)', fontsize=10, color=INK, pad=6)
ax.grid(alpha=0.25)
ax.spines[['top', 'right']].set_visible(False)
panel_label(ax, 'B')

# ---------- C: sorted paired differences ----------
ax = axes[1][0]
order = np.argsort(diffs)
ax.bar(np.arange(1, n + 1), diffs[order], width=1.0,
       color=[TEAL if d < 0 else RED for d in diffs[order]], edgecolor='none')
ax.axhline(0, color=INK, lw=0.8)
ax.set_xlabel('Seed (sorted by paired difference)', fontsize=9.5, color=INK)
ax.set_ylabel('nRF difference (Fusang − FT2)', fontsize=9.5, color=INK)
ax.set_title(f'Paired differences, range [{diffs.min():.3f}, {diffs.max():.3f}]\n'
             f'negative = Fusang better ({st["wins_a"]}/{n} seeds)',
             fontsize=9.5, color=INK, pad=6)
ax.grid(alpha=0.25, axis='y')
ax.spines[['top', 'right']].set_visible(False)
panel_label(ax, 'C')

# ---------- D: histogram ----------
ax = axes[1][1]
ax.hist(diffs, bins=22, color=TEAL, alpha=0.75, edgecolor='white')
ax.axvline(0, color=INK, lw=1.0)
ax.axvline(diffs.mean(), color=AMBER, ls='--', lw=1.6,
           label=f'mean = {diffs.mean():+.4f}')
ax.set_xlabel('nRF difference (Fusang − FT2)', fontsize=9.5, color=INK)
ax.set_ylabel('Seeds', fontsize=9.5, color=INK)
ax.set_title(f'130 seeds; Wilcoxon p={st["p"]:.4f} (zsplit, integer RF), '
             f'paired d={d_fs:.2f}\nbootstrap 95% CI for d [{trunc2(ci_fs[0]):+.2f}, '
             f'{trunc2(ci_fs[1]):+.2f}]; Fusang better {st["wins_a"]}/{n} '
             f'(FT2 {st["wins_b"]}, ties {st["ties"]})',
             fontsize=9.2, color=INK, pad=6)
ax.legend(fontsize=8.5, frameon=False)
ax.grid(alpha=0.25, axis='y')
ax.spines[['top', 'right']].set_visible(False)
panel_label(ax, 'D')

fig.suptitle('Supplementary Figure S4. 130-seed indel benchmark: distributions and outlier analysis',
             fontsize=12.5, fontweight='bold', color=INK, y=0.985)
fig.tight_layout(rect=[0, 0, 1, 0.955])
fig.savefig('Supplementary_Figure_S4.pdf'); fig.savefig('Supplementary_Figure_S4.png', dpi=300)
print('S4 done | n=%d | Fusang %.4f+-%.4f FT2 %.4f+-%.4f | p=%.4f d=%+.3f CI %s'
      % (n, st['mean_a'], st['sd_a'], st['mean_b'], st['sd_b'], st['p'], d_fs,
         np.round(ci_fs, 3).tolist()))
