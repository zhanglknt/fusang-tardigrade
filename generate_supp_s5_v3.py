"""Supplementary Figure S5 (v3): 16S rRNA validation detail — 4 panels.
A: taxonomic signal, original analysis vs archived-tree rerun (grouped bars)
B: Fusang 73-taxon cladogram (phylum-coloured tips)
C: FastTree2 73-taxon cladogram (same colouring) — direct comparison
D: external-reference recovery metrics (monophyly, nRF vs NCBI taxonomy)
"""
import os
os.environ['MPLCONFIGDIR'] = 'matplotlib_cache'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import json
from Bio import Phylo

INK, GREY = '#1a2332', '#64748b'
BLUE, TEAL, AMBER, RED, GREEN = '#2563eb', '#0d9488', '#d97706', '#dc2626', '#16a34a'
PURPLE, PINK = '#7c3aed', '#db2777'

meta = json.load(open('real_data/real_16S_metadata.json'))
TOP_PHYLA = ['Pseudomonadota', 'Bacillota', 'Actinomycetota',
             'Bacteroidota', 'Cyanobacteriota', 'Euryarchaeota']
PCOL = dict(zip(TOP_PHYLA, [BLUE, TEAL, AMBER, PURPLE, PINK, GREEN]))

fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.2))
fig.patch.set_facecolor('white')

def panel_label(ax, s):
    ax.text(-0.12, 1.05, s, transform=ax.transAxes, fontsize=14,
            fontweight='bold', color=INK, va='top')

def draw_cladogram(ax, nwk_path, title):
    tree = Phylo.read(nwk_path, 'newick')
    tip_y = {t: i for i, t in enumerate(tree.get_terminals())}
    def rec(clade, x):
        if clade.is_terminal():
            return tip_y[clade]
        ys = [rec(c, x + 1) for c in clade.clades]
        yc = (ys[0] + ys[-1]) / 2.0
        for c, y in zip(clade.clades, ys):
            ax.plot([x, x + 1], [y, y], color=INK, lw=0.7)
            ax.plot([x, x], [min(y, yc), max(y, yc)], color=INK, lw=0.7)
        return yc
    rec(tree.root, 0)
    maxd = max(len(tree.get_path(t)) for t in tree.get_terminals())
    for tip, y in tip_y.items():
        ph = meta.get(tip.name, {}).get('phylum')
        ax.plot(maxd + 0.4, y, 'o', ms=2.9, color=PCOL.get(ph, GREY), markeredgecolor='none')
    ax.set_xlim(-0.5, maxd + 1.6)
    ax.set_ylim(-1, len(tip_y))
    ax.axis('off')
    ax.set_title(title, fontsize=10, color=INK, pad=6)

# ---------- A: original vs rerun signal ----------
ax = axes[0][0]
ranks = ['Order', 'Phylum']
orig = [12.6, 4.6]
rerun = [7.3, 3.8]
x = np.arange(2); w = 0.34
b1 = ax.bar(x - w/2, orig, w, color=TEAL, label='original analysis (74 taxa)')
b2 = ax.bar(x + w/2, rerun, w, color=BLUE, label='verification rerun (archived 73 taxa)')
for bars, ps in ((b1, ['p<0.01', 'p<0.05']), (b2, ['p=0.027', 'p=0.003'])):
    for b, pv in zip(bars, ps):
        ax.text(b.get_x() + b.get_width()/2, b.get_height() + 0.3,
                f'{b.get_height():.1f}%\n{pv}', ha='center', fontsize=7.6, color=INK)
ax.set_xticks(x); ax.set_xticklabels(ranks, fontsize=9.5)
ax.set_ylabel('Same-rank distance reduction (%)', fontsize=9, color=INK)
ax.set_ylim(0, 16.5)
ax.set_title('Taxonomic signal, two independent runs', fontsize=10.5, color=INK, pad=6)
ax.legend(fontsize=7.8, frameon=False, loc='upper right')
ax.grid(alpha=0.25, axis='y')
for s in ('top', 'right'): ax.spines[s].set_visible(False)
panel_label(ax, 'A')

# ---------- B/C: cladograms ----------
draw_cladogram(axes[0][1], 'real_data/real_16S_full_fusang_k5_gap2.nwk',
               'Fusang tree (k=5,gap2, no alignment, 1.2 s)')
panel_label(axes[0][1], 'B')
draw_cladogram(axes[1][0], 'real_data/real_16S_full_ft2.nwk',
               'FastTree2 tree (MAFFT + ML)')
panel_label(axes[1][0], 'C')
handles = [plt.Line2D([0], [0], marker='o', color='none', markerfacecolor=PCOL[p],
                      markersize=5, label=p) for p in TOP_PHYLA]
handles.append(plt.Line2D([0], [0], marker='o', color='none', markerfacecolor=GREY,
                          markersize=5, label='other'))
axes[0][1].legend(handles=handles, fontsize=5.6, frameon=False, loc='lower right',
                  handletextpad=0.1, borderaxespad=0.1, labelspacing=0.25)

# ---------- D: recovery metrics ----------
ax = axes[1][1]
x = np.arange(2); w = 0.34
fus_v = [8 / 12 * 100, 0.68 * 100]
ft2_v = [10 / 12 * 100, 0.45 * 100]
b1 = ax.bar(x - w/2, fus_v, w, color=TEAL, label='Fusang')
b2 = ax.bar(x + w/2, ft2_v, w, color=GREY, label='FastTree2')
ax.text(x[0] - w/2, fus_v[0] + 2, '8/12\n(66.7%)', ha='center', fontsize=8, color=INK, fontweight='bold')
ax.text(x[0] + w/2, ft2_v[0] + 2, '10/12\n(83.3%)', ha='center', fontsize=8, color=INK, fontweight='bold')
ax.text(x[1] - w/2, fus_v[1] + 2, 'nRF=0.68', ha='center', fontsize=8, color=INK, fontweight='bold')
ax.text(x[1] + w/2, ft2_v[1] + 2, 'nRF=0.45', ha='center', fontsize=8, color=INK, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(['Monophyletic groups\nrecovered (% of 12)', 'nRF vs NCBI taxonomy\n(\u00d7100, \u2193 better)'], fontsize=8.5)
ax.set_ylim(0, 105)
ax.set_ylabel('value', fontsize=9, color=INK)
ax.set_title('Recovery against NCBI taxonomy (descriptive)', fontsize=10.5, color=INK, pad=6)
ax.legend(fontsize=8.5, frameon=False, loc='upper left')
ax.grid(alpha=0.25, axis='y')
for s in ('top', 'right'): ax.spines[s].set_visible(False)
panel_label(ax, 'D')

fig.suptitle('Supplementary Figure S5. Real 16S rRNA validation (73\u201374 taxa)',
             fontsize=13, fontweight='bold', color=INK, y=0.985)
fig.tight_layout(rect=[0, 0, 1, 0.955])
fig.savefig('Supplementary_Figure_S5.pdf'); fig.savefig('Supplementary_Figure_S5.png', dpi=300)
print('S5 done')
