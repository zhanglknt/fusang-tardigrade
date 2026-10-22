"""Figure 2 (v3): 130-seed statistical benchmark — 4 panels.

A: violin Fusang vs FT2
B: sorted per-seed paired differences (FT2 - Fusang; positive = Fusang better)
C: cumulative wins across the seed set (chronological seed order)
D: paired-difference histogram with Wilcoxon p and Cohen's d

All numbers come from ``figure_data.bench_130seeds()``, which reads
``indelible_n200_indel/bench_n200_indel_130seeds_nrfstandard.csv`` — the canonical
130-seed per-seed file produced with the standard root-independent unrooted RF. The
previous ``table3_corrected.csv`` input was written by the superseded root-dependent
implementation and is no longer used.
"""
import os
os.environ['MPLCONFIGDIR'] = 'matplotlib_cache'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import figure_data as FD

INK, GREY = '#1a2332', '#64748b'
BLUE, TEAL, AMBER, RED = '#2563eb', '#0d9488', '#d97706', '#dc2626'

st = FD.bench_130seeds()
fus, ft2, seeds = st['fus'], st['ft2'], st['seeds']
diffs = ft2 - fus  # positive = Fusang better (Figure 2 legend orientation)
n = len(diffs)

fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.6))
fig.patch.set_facecolor('white')


def panel_label(ax, s):
    ax.text(-0.14, 1.06, s, transform=ax.transAxes, fontsize=14,
            fontweight='bold', color=INK, va='top')


# ---------- A: violins ----------
ax = axes[0][0]
vp = ax.violinplot([fus, ft2], showmedians=True, showextrema=True)
for body, c in zip(vp['bodies'], [TEAL, GREY]):
    body.set_facecolor(c); body.set_alpha(0.75); body.set_edgecolor(INK); body.set_linewidth(1)
for part in ('cbars', 'cmins', 'cmaxes', 'cmedians'):
    vp[part].set_color(INK); vp[part].set_linewidth(1.2)
rng = np.random.default_rng(7)
for i, arr in enumerate([fus, ft2]):
    ax.scatter(rng.normal(i + 1, 0.035, len(arr)), arr, s=5, color=INK, alpha=0.25, zorder=3)
ax.set_xticks([1, 2]); ax.set_xticklabels(['Fusang', 'FastTree2'], fontsize=10)
ax.set_ylabel('nRF vs true tree', fontsize=9.5, color=INK)
ax.set_title(f'Distributions (n={n} seeds)\nFusang {fus.mean():.3f}±{fus.std(ddof=1):.3f}  vs  '
             f'FT2 {ft2.mean():.3f}±{ft2.std(ddof=1):.3f}', fontsize=9.5, color=INK, pad=6)
ax.spines[['top', 'right']].set_visible(False)
panel_label(ax, 'A')

# ---------- B: sorted per-seed differences ----------
ax = axes[0][1]
order = np.argsort(diffs)
x = np.arange(1, n + 1)
ax.bar(x, diffs[order], width=1.0,
       color=[TEAL if d > 0 else RED for d in diffs[order]], edgecolor='none')
ax.axhline(0, color=INK, lw=0.9)
ax.set_xlabel('Seed (sorted by paired difference)', fontsize=10, color=INK)
ax.set_ylabel('nRF difference (FT2 − Fusang)', fontsize=9.5, color=INK)
ax.set_title(f'Per-seed differences, range [{diffs.min():.3f}, {diffs.max():.3f}]\n'
             f'positive = Fusang better ({st["wins_a"]}/{n} seeds)',
             fontsize=9.5, color=INK, pad=6)
ax.spines[['top', 'right']].set_visible(False)
panel_label(ax, 'B')

# ---------- C: cumulative wins (chronological seed order) ----------
ax = axes[1][0]
seed_order = np.argsort(seeds)
cum_win = np.cumsum(diffs[seed_order] > 0) / (np.arange(n) + 1) * 100
ax.plot(np.arange(1, n + 1), cum_win, color=TEAL, lw=2.2)
ax.axhline(50, color=GREY, lw=1, ls='--')
ax.text(n * 0.98, 51, '50% (tie)', fontsize=8, color=GREY, ha='right')
ax.set_xlabel('Seeds evaluated (chronological order)', fontsize=10, color=INK)
ax.set_ylabel('Cumulative Fusang win rate (%)', fontsize=9.5, color=INK)
ax.set_ylim(35, 75)
ax.set_title(f'Win rate stabilises at {st["wins_a"]}/{n} '
             f'({100.0 * st["wins_a"] / n:.1f}%)',
             fontsize=9.5, color=INK, pad=6)
ax.spines[['top', 'right']].set_visible(False)
panel_label(ax, 'C')

# ---------- D: paired-difference histogram ----------
ax = axes[1][1]
ax.hist(diffs, bins=22, color=TEAL, edgecolor='white', lw=0.6, alpha=0.9)
ax.axvline(0, color=INK, lw=1)
ax.axvline(diffs.mean(), color=AMBER, lw=2, ls='--')
ax.text(diffs.mean(), ax.get_ylim()[1] * 0.92, f'  mean {diffs.mean():+.4f}',
        color=AMBER, fontsize=9, fontweight='bold', va='top')
ax.set_xlabel('Paired nRF difference (FT2 − Fusang)', fontsize=10, color=INK)
ax.set_ylabel('Seeds', fontsize=9.5, color=INK)
# manuscript unified Wilcoxon protocol (integer RF + zsplit): p=0.0507, borderline n.s.
ax.set_title(f'Paired Wilcoxon p={st["p"]:.4f} (borderline, n.s. at α=0.05)\n'
             f"paired Cohen's d={st['d']:+.2f} (FT2 − Fusang; positive = Fusang better)",
             fontsize=9.5, color=INK, pad=6)
ax.spines[['top', 'right']].set_visible(False)
panel_label(ax, 'D')

fig.suptitle('Figure 2. 130-seed statistical benchmark (n=200, indel rate 0.02)',
             fontsize=14, fontweight='bold', color=INK, y=0.985)
fig.tight_layout(rect=[0, 0, 1, 0.955])
fig.savefig('Figure2.pdf')
fig.savefig('Figure2.png', dpi=300)
print('Figure2 done | n=%d | Fusang %.4f+-%.4f FT2 %.4f+-%.4f | p=%.4f d=%+.3f CI %s |'
      ' wins %d/%d/%d | diff range %.3f..%.3f | >0.3: %d'
      % (n, st['mean_a'], st['sd_a'], st['mean_b'], st['sd_b'], st['p'], st['d'],
         np.round(st['ci'], 3).tolist(), st['wins_a'], st['wins_b'], st['ties'],
         diffs.min(), diffs.max(), st['n_gt_03']))
