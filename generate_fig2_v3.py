"""Figure 2 (v3): 130-seed statistical benchmark — 4 panels.
Data: table3_corrected.csv (120 paired seeds; 8 catastrophic nRF>0.3 exclusions -> 112).
A: violin Fusang vs FT2 (112 seeds)
B: sorted per-seed paired differences (FT2 - Fusang; positive = Fusang better)
C: cumulative wins across seed set
D: paired-difference histogram with Wilcoxon p and Cohen's d
"""
import os
os.environ['MPLCONFIGDIR'] = 'matplotlib_cache'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import csv
from scipy.stats import wilcoxon

INK, GREY = '#1a2332', '#64748b'
BLUE, TEAL, AMBER, RED = '#2563eb', '#0d9488', '#d97706', '#dc2626'

rows = list(csv.DictReader(open('table3_corrected.csv')))
pairs = [(int(r['seed']), float(r['fusang_nrf']), float(r['ft2_nrf'])) for r in rows]
cat = {s for s, f, t in pairs if f > 0.3 or t > 0.3}
valid = [(s, f, t) for s, f, t in pairs if s not in cat]
fus = np.array([f for _, f, _ in valid])
ft2 = np.array([t for _, _, t in valid])
diffs = ft2 - fus  # positive = Fusang better
w112 = wilcoxon(fus, ft2)
dz = diffs.mean() / diffs.std(ddof=1)  # paired Cohen's d (FT2-Fusang orientation)
wins_all = sum(1 for _, f, t in pairs if f < t)
wins_112 = int((diffs > 0).sum())

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
ax.set_title(f'Distributions (n=112 seeds)\nFusang {fus.mean():.3f}±{fus.std(ddof=1):.3f}  vs  '
             f'FT2 {ft2.mean():.3f}±{ft2.std(ddof=1):.3f}', fontsize=9.5, color=INK, pad=6)
ax.spines[['top', 'right']].set_visible(False)
panel_label(ax, 'A')

# ---------- B: sorted per-seed differences ----------
ax = axes[0][1]
order = np.argsort(diffs)
x = np.arange(1, len(diffs) + 1)
ax.bar(x, diffs[order], width=1.0,
       color=[TEAL if d > 0 else RED for d in diffs[order]], edgecolor='none')
ax.axhline(0, color=INK, lw=0.9)
ax.set_xlabel('Seed (sorted by paired difference)', fontsize=10, color=INK)
ax.set_ylabel('nRF difference (FT2 − Fusang)', fontsize=9.5, color=INK)
ax.set_title(f'Per-seed differences, range [{diffs.min():.3f}, {diffs.max():.3f}]\n'
             f'positive = Fusang better ({wins_112}/112 seeds)',
             fontsize=9.5, color=INK, pad=6)
ax.spines[['top', 'right']].set_visible(False)
panel_label(ax, 'B')

# ---------- C: cumulative wins ----------
ax = axes[1][0]
seed_order = np.argsort([s for s, _, _ in valid])
cum_win = np.cumsum(diffs[seed_order] > 0) / (np.arange(len(diffs)) + 1) * 100
ax.plot(np.arange(1, len(diffs) + 1), cum_win, color=TEAL, lw=2.2)
ax.axhline(50, color=GREY, lw=1, ls='--')
ax.text(len(diffs) * 0.98, 51, '50% (tie)', fontsize=8, color=GREY, ha='right')
ax.set_xlabel('Seeds evaluated (chronological order)', fontsize=10, color=INK)
ax.set_ylabel('Cumulative Fusang win rate (%)', fontsize=9.5, color=INK)
ax.set_ylim(35, 75)
ax.set_title(f'Win rate: {wins_all}/120 valid seeds (50.8%) →\n'
             f'{wins_112}/112 after outlier exclusion (53.6%)',
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
ax.set_title(f'Paired Wilcoxon p={w112.pvalue:.3f} (borderline)\n'
             f"paired Cohen's d={abs(dz):.2f} (small, direction: Fusang better)",
             fontsize=9.5, color=INK, pad=6)
ax.spines[['top', 'right']].set_visible(False)
panel_label(ax, 'D')

fig.suptitle('Figure 2. 130-seed statistical benchmark (n=200, indel rate 0.02)',
             fontsize=14, fontweight='bold', color=INK, y=0.985)
fig.tight_layout(rect=[0, 0, 1, 0.955])
fig.savefig('Figure2.pdf')
fig.savefig('Figure2.png', dpi=300)
print('Figure2 done | 120-valid wins:', wins_all, '| 112 wins:', wins_112,
      '| p=%.4f' % w112.pvalue, '| d=%.3f' % dz,
      '| diff range %.3f..%.3f' % (diffs.min(), diffs.max()))
