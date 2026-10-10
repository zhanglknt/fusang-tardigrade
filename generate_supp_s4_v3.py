"""Supplementary Figure S4 (v3): 130-seed benchmark detail — 4 panels.
All data from table3_corrected.csv (120 paired seeds, 8 with nRF>0.3).
A: violins for all 120 paired seeds (outliers visible)
B: per-seed scatter Fusang vs FT2 with exclusion threshold
C: sorted paired diffs, all 120 (excluded seeds in red) — exclusion rationale
D: histogram of paired diffs after exclusion (112 seeds; p=0.052, d=+0.20)
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
TEAL, AMBER, RED = '#0d9488', '#d97706', '#dc2626'

rows = list(csv.DictReader(open('table3_corrected.csv')))
seed = np.array([int(r['seed']) for r in rows])
fus = np.array([float(r['fusang_nrf']) for r in rows])
ft2 = np.array([float(r['ft2_nrf']) for r in rows])
out = (fus > 0.3) | (ft2 > 0.3)
print(f'120 paired; excluded (nRF>0.3): {out.sum()} seeds {sorted(seed[out])}')
keep = ~out
diffs = ft2[keep] - fus[keep]  # positive = Fusang better
d_paired = diffs.mean() / diffs.std(ddof=1)
p = wilcoxon(ft2[keep], fus[keep]).pvalue
print(f'112 kept: fus {fus[keep].mean():.4f}+-{fus[keep].std(ddof=1):.4f} '
      f'ft2 {ft2[keep].mean():.4f}+-{ft2[keep].std(ddof=1):.4f} p={p:.4f} d={d_paired:.3f}')

fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.6))
fig.patch.set_facecolor('white')

def panel_label(ax, s):
    ax.text(-0.14, 1.06, s, transform=ax.transAxes, fontsize=14,
            fontweight='bold', color=INK, va='top')

# ---------- A: violins all 120 ----------
ax = axes[0][0]
vp = ax.violinplot([fus, ft2], positions=[1, 2], widths=0.62, showmeans=True)
for b, c in zip(vp['bodies'], [TEAL, GREY]):
    b.set_facecolor(c); b.set_alpha(0.65)
for key in ('cbars', 'cmins', 'cmaxes', 'cmeans'):
    vp[key].set_color(INK); vp[key].set_linewidth(1.1)
ax.set_xticks([1, 2]); ax.set_xticklabels(['Fusang', 'FastTree2'], fontsize=9.5)
ax.set_ylabel('nRF vs true tree', fontsize=9.5, color=INK)
ax.set_title('All 120 paired seeds (outliers visible)', fontsize=10.5, color=INK, pad=6)
ax.grid(alpha=0.25, axis='y')
for s in ('top', 'right'): ax.spines[s].set_visible(False)
panel_label(ax, 'A')

# ---------- B: scatter with exclusion threshold ----------
ax = axes[0][1]
ax.scatter(ft2[keep], fus[keep], s=18, color=TEAL, alpha=0.65, edgecolor='none',
           label=f'included (n={keep.sum()})')
ax.scatter(ft2[out], fus[out], s=34, color=RED, marker='x', lw=1.6,
           label=f'excluded nRF>0.3 (n={out.sum()})')
lim = [0, max(ft2.max(), fus.max()) * 1.05]
ax.plot(lim, lim, ':', color=GREY, lw=1.2)
ax.axhline(0.3, color=RED, ls='--', lw=0.9, alpha=0.6)
ax.axvline(0.3, color=RED, ls='--', lw=0.9, alpha=0.6)
ax.set_xlabel('FastTree2 nRF', fontsize=9.5, color=INK)
ax.set_ylabel('Fusang nRF', fontsize=9.5, color=INK)
ax.set_title('Per-seed agreement and exclusion rule', fontsize=10.5, color=INK, pad=6)
ax.legend(fontsize=8, frameon=False, loc='upper left')
ax.grid(alpha=0.25)
for s in ('top', 'right'): ax.spines[s].set_visible(False)
panel_label(ax, 'B')

# ---------- C: sorted diffs all 120 ----------
ax = axes[1][0]
d_all = ft2 - fus
order = np.argsort(d_all)
x = np.arange(1, len(d_all) + 1)
cols = [RED if out[i] else (TEAL if d_all[i] > 0 else GREY) for i in order]
ax.bar(x, d_all[order], width=1.0, color=cols, edgecolor='none')
ax.axhline(0, color=INK, lw=0.8)
ax.set_xlabel('Seed (sorted by paired difference)', fontsize=9.5, color=INK)
ax.set_ylabel('nRF difference (FT2 \u2212 Fusang)', fontsize=9.5, color=INK)
ax.set_title('All 120 seeds: 8 catastrophic outliers (red)\ndrive the long tails — excluded from headline stats',
             fontsize=9.5, color=INK, pad=6)
ax.grid(alpha=0.25, axis='y')
for s in ('top', 'right'): ax.spines[s].set_visible(False)
panel_label(ax, 'C')

# ---------- D: histogram 112 ----------
ax = axes[1][1]
ax.hist(diffs, bins=22, color=TEAL, alpha=0.75, edgecolor='white')
ax.axvline(0, color=INK, lw=1.0)
ax.axvline(diffs.mean(), color=AMBER, ls='--', lw=1.6,
           label=f'mean = {diffs.mean():+.4f}')
ax.set_xlabel('nRF difference (FT2 \u2212 Fusang)', fontsize=9.5, color=INK)
ax.set_ylabel('Seeds', fontsize=9.5, color=INK)
ax.set_title(f'112 seeds after exclusion\nWilcoxon p={p:.3f}, paired d={d_paired:+.2f}, '
             f'Fusang better {(diffs>0).sum()}/112', fontsize=10, color=INK, pad=6)
ax.legend(fontsize=8.5, frameon=False)
ax.grid(alpha=0.25, axis='y')
for s in ('top', 'right'): ax.spines[s].set_visible(False)
panel_label(ax, 'D')

fig.suptitle('Supplementary Figure S4. 130-seed indel benchmark: distributions and outlier analysis',
             fontsize=12.5, fontweight='bold', color=INK, y=0.985)
fig.tight_layout(rect=[0, 0, 1, 0.955])
fig.savefig('Supplementary_Figure_S4.pdf'); fig.savefig('Supplementary_Figure_S4.png', dpi=300)
print('S4 done')
