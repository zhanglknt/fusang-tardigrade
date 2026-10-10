"""Supplementary Figure S1 (v3): k-mer grid search — 6 panels.
A-E: nRF heatmaps (k x gap) at n=50/100/200/500/1000 (data/S1_grid_search_results.csv)
F:   nRF vs n for default k=5,gap2 vs best-per-n configuration
"""
import os
os.environ['MPLCONFIGDIR'] = 'matplotlib_cache'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import csv

INK, GREY = '#1a2332', '#64748b'
TEAL, AMBER, RED = '#0d9488', '#d97706', '#dc2626'

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
    M = np.full((len(KS), len(GS)), np.nan)
    for a, k in enumerate(KS):
        for b, g in enumerate(GS):
            M[a, b] = val.get((n, k, g), np.nan)
    im = ax.imshow(M, cmap='RdYlGn_r', vmin=0, vmax=0.20, aspect='auto')
    best = np.nanargmin(M)
    for a in range(len(KS)):
        for b in range(len(GS)):
            if not np.isnan(M[a, b]):
                star = '*' if (a, b) == np.unravel_index(best, M.shape) else ''
                ax.text(b, a, f'{M[a,b]:.3f}{star}', ha='center', va='center',
                        fontsize=6.8, color='white' if M[a, b] > 0.10 else INK,
                        fontweight='bold' if star else 'normal')
    ax.set_xticks(range(len(GS)))
    ax.set_xticklabels(['gap0', 'gap1', 'gap2', 'gap3'], fontsize=8)
    ax.set_yticks(range(len(KS)))
    ax.set_yticklabels([f'k={k}' for k in KS], fontsize=8)
    ax.set_title(f'n = {n}', fontsize=10.5, color=INK, pad=5)
    fig.colorbar(im, ax=ax, shrink=0.8)
    panel_label(ax, labels[i])

# ---------- F: nRF vs n, default vs best ----------
ax = axes[1][2]
def_nrf, best_nrf, best_cfg = [], [], []
for n in NS:
    def_nrf.append(val.get((n, 5, 'gap2'), np.nan))
    cands = {(k, g): v for (nn, k, g), v in val.items() if nn == n}
    bk = min(cands, key=cands.get)
    best_nrf.append(cands[bk])
    best_cfg.append(f'k={bk[0]},{bk[1]}')
ax.plot(NS, def_nrf, 's-', color=TEAL, lw=2.2, ms=6.5, label='default k=5,gap2')
ax.plot(NS, best_nrf, 'o--', color=AMBER, lw=1.8, ms=6, label='best per-n config')
for n, b, c in zip(NS, best_nrf, best_cfg):
    ax.annotate(c, (n, b), textcoords='offset points', xytext=(0, -13),
                ha='center', fontsize=6.8, color=AMBER)
ax.set_xscale('log'); ax.set_xticks(NS); ax.set_xticklabels([str(n) for n in NS], fontsize=8.5)
ax.set_xlabel('Dataset size n (log scale)', fontsize=9.5, color=INK)
ax.set_ylabel('nRF', fontsize=9.5, color=INK)
ax.set_title('Default vs best configuration', fontsize=10.5, color=INK, pad=5)
ax.legend(fontsize=8.5, frameon=False, loc='upper right')
ax.grid(alpha=0.25)
for s in ('top', 'right'): ax.spines[s].set_visible(False)
panel_label(ax, 'F')

fig.suptitle('Supplementary Figure S1. Full k-mer parameter grid search (* = best per panel)',
             fontsize=12.5, fontweight='bold', color=INK, y=0.99)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig('Supplementary_Figure_S1.pdf'); fig.savefig('Supplementary_Figure_S1.png', dpi=300)
print('S1 done')
