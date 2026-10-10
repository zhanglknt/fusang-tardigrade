"""Supplementary Figure S2 (v3): dimensionality vs accuracy — 6 panels.
A-E: nRF vs feature dimension (4^k), one line per gap pattern, per n
F:   nRF vs n at k=5 (1024-dim) for each gap pattern
"""
import os
os.environ['MPLCONFIGDIR'] = 'matplotlib_cache'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import csv

INK, GREY = '#1a2332', '#64748b'
BLUE, TEAL, AMBER, RED = '#2563eb', '#0d9488', '#d97706', '#dc2626'
GCOL = {'none': BLUE, 'gap1': AMBER, 'gap2': TEAL, 'gap3': RED}

rows = list(csv.DictReader(open('data/S1_grid_search_results.csv')))
NS = [50, 100, 200, 500, 1000]
KS = [4, 5, 6, 7, 8]
GS = ['none', 'gap1', 'gap2', 'gap3']
val = {}
for r in rows:
    val[(int(r['n_taxa']), int(r['k']), r['gap_pattern'])] = float(r['nRF'])

fig, axes = plt.subplots(2, 3, figsize=(12.2, 7.2))
fig.patch.set_facecolor('white')

def panel_label(ax, s):
    ax.text(-0.14, 1.05, s, transform=ax.transAxes, fontsize=13,
            fontweight='bold', color=INK, va='top')

labels = 'ABCDEF'
for i, n in enumerate(NS):
    ax = axes[i // 3][i % 3]
    for g in GS:
        xs, ys = [], []
        for k in KS:
            v = val.get((n, k, g))
            if v is not None:
                xs.append(4 ** k); ys.append(v)
        if xs:
            ax.plot(xs, ys, 'o-', color=GCOL[g], lw=1.6, ms=4.5,
                    label='contiguous' if g == 'none' else g)
    ax.axhline(0.02, color=GREY, ls=':', lw=0.9, alpha=0.6)
    ax.set_xscale('log')
    ax.set_xlabel('Feature dimension (4$^k$)', fontsize=9, color=INK)
    ax.set_ylabel('nRF', fontsize=9, color=INK)
    ax.set_ylim(0, 0.21)
    ax.set_title(f'n = {n}', fontsize=10.5, color=INK, pad=5)
    ax.legend(fontsize=7, frameon=False)
    ax.grid(alpha=0.25)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
    panel_label(ax, labels[i])

# ---------- F: nRF vs n at k=5 ----------
ax = axes[1][2]
for g in GS:
    xs, ys = [], []
    for n in NS:
        v = val.get((n, 5, g))
        if v is not None:
            xs.append(n); ys.append(v)
    ax.plot(xs, ys, 'o-', color=GCOL[g], lw=1.8, ms=5,
            label='contiguous (k=5)' if g == 'none' else f'k=5,{g}')
ax.set_xscale('log'); ax.set_xticks(NS); ax.set_xticklabels([str(n) for n in NS], fontsize=8.5)
ax.set_xlabel('Dataset size n (log scale)', fontsize=9.5, color=INK)
ax.set_ylabel('nRF', fontsize=9.5, color=INK)
ax.set_title('k=5 across dataset sizes\n(1024-dim: gap2 best at n>=100)', fontsize=10, color=INK, pad=5)
ax.legend(fontsize=7.5, frameon=False)
ax.grid(alpha=0.25)
for s in ('top', 'right'): ax.spines[s].set_visible(False)
panel_label(ax, 'F')

fig.suptitle('Supplementary Figure S2. Dimensionality vs phylogenetic accuracy',
             fontsize=12.5, fontweight='bold', color=INK, y=0.99)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig('Supplementary_Figure_S2.pdf'); fig.savefig('Supplementary_Figure_S2.png', dpi=300)
print('S2 done')
