"""Supplementary Figure S5 (v3): 16S rRNA validation detail — 4 panels.
A: taxonomic signal (real_data/real_16S_ranksignal.json; corrected NCBI annotation)
B: Fusang 73-taxon cladogram (phylum-coloured tips)
C: rebuilt FastTree2 73-taxon cladogram (MAFFT 7.526 --localpair + FastTree2 -nt -gtr)
D: external-reference recovery metrics (monophyly, nRF vs NCBI taxonomy)
All values follow the corrected NCBI Taxonomy annotation (Supplementary Note S17);
earlier evaluations against the erroneous shipped labels are superseded.
"""
import os
os.environ['MPLCONFIGDIR'] = 'matplotlib_cache'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import json
from Bio import Phylo

import figure_data as FD

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
    return len(tip_y)

# ---------- A: taxonomic signal (Table 9 values, corrected annotation) ----------
ax = axes[0][0]
rs = FD.rank_signal_16s()
ranks = ['phylum', 'order', 'family', 'genus']
fus = [rs['Fusang_k5_gap2'][r]['reduction_pct'] for r in ranks]
ft2 = [rs['FastTree2_MAFFT'][r]['reduction_pct'] for r in ranks]
x = np.arange(len(ranks)); w = 0.36
b1 = ax.bar(x - w/2, fus, w, color=TEAL, label='Fusang (k=5,gap2)')
b2 = ax.bar(x + w/2, ft2, w, color=GREY, label='rebuilt FastTree2 (MAFFT)')
for b, v in zip(b1, fus):
    ax.text(b.get_x() + b.get_width()/2, v + 1.5, f'{v:.1f}', ha='center',
            fontsize=8.3, fontweight='bold', color=TEAL)
for b, v in zip(b2, ft2):
    ax.text(b.get_x() + b.get_width()/2, v + 1.5, f'{v:.1f}', ha='center',
            fontsize=8.3, fontweight='bold', color=INK)
ax.set_xticks(x)
ax.set_xticklabels([r.capitalize() for r in ranks], fontsize=9.5)
ax.set_ylabel('Same-rank distance reduction (%)', fontsize=9, color=INK)
ax.set_ylim(0, 112)
ax.set_title('Taxonomic signal (Table 9; all p < 1\u00d710\u207b\u2074)', fontsize=10.5,
             color=INK, pad=6)
ax.legend(fontsize=7.8, frameon=False, loc='upper left')
ax.grid(alpha=0.25, axis='y')
for s in ('top', 'right'): ax.spines[s].set_visible(False)
panel_label(ax, 'A')

# ---------- B/C: cladograms ----------
nB = draw_cladogram(axes[0][1], 'real_data/real_16S_full_fusang_k5_gap2.nwk',
                    'Fusang tree (k=5,gap2, no alignment, 0.11 s)')
panel_label(axes[0][1], 'B')
nC = draw_cladogram(axes[1][0], 'real_data/real_16S_full_ft2_mafft.nwk',
                    'Rebuilt FastTree2 tree (MAFFT 7.526 --localpair + FastTree2 -nt -gtr)')
panel_label(axes[1][0], 'C')
handles = [plt.Line2D([0], [0], marker='o', color='none', markerfacecolor=PCOL[p],
                      markersize=5, label=p) for p in TOP_PHYLA]
handles.append(plt.Line2D([0], [0], marker='o', color='none', markerfacecolor=GREY,
                          markersize=5, label='other'))
axes[0][1].legend(handles=handles, fontsize=5.6, frameon=False, loc='lower right',
                  handletextpad=0.1, borderaxespad=0.1, labelspacing=0.25)

# ---------- D: recovery against NCBI taxonomy (corrected annotation) ----------
# Fusang recovers 3/7 assessable phyla and 10/14 assessable orders as monophyletic;
# the rebuilt FastTree2 tree recovers 6/7 phyla and 11/14 orders. nRF vs the
# taxonomy-derived reference: 0.414 vs 0.343 (RF/(2n-6)); 0.547 vs 0.453 (RF/(s1+s2)).
ax = axes[1][1]
groups = ['Phyla\nmonophyletic', 'Orders\nmonophyletic',
          'nRF vs taxonomy\nRF/(2n\u22126)', 'nRF vs taxonomy\nRF/(s\u2081+s\u2082)']
fus_v = [3/7 * 100, 10/14 * 100, 0.414 * 100, 0.547 * 100]
ft2_v = [6/7 * 100, 11/14 * 100, 0.343 * 100, 0.453 * 100]
fus_lab = ['3/7', '10/14', '0.414', '0.547']
ft2_lab = ['6/7', '11/14', '0.343', '0.453']
x = np.arange(len(groups)); w = 0.36
b1 = ax.bar(x - w/2, fus_v, w, color=TEAL, label='Fusang')
b2 = ax.bar(x + w/2, ft2_v, w, color=GREY, label='rebuilt FastTree2')
for b, lab, v in zip(b1, fus_lab, fus_v):
    ax.text(b.get_x() + b.get_width()/2, v + 2, lab, ha='center', fontsize=8,
            color=TEAL, fontweight='bold')
for b, lab, v in zip(b2, ft2_lab, ft2_v):
    ax.text(b.get_x() + b.get_width()/2, v + 2, lab, ha='center', fontsize=8,
            color=INK, fontweight='bold')
ax.set_xticks(x); ax.set_xticklabels(groups, fontsize=8)
ax.set_ylim(0, 112)
ax.set_ylabel('share (%; nRF shown \u00d7100)', fontsize=9, color=INK)
ax.set_title('Recovery against NCBI taxonomy (13 phyla / 40 orders)', fontsize=10,
             color=INK, pad=6)
ax.legend(fontsize=8.5, frameon=False, loc='upper right')
ax.grid(alpha=0.25, axis='y')
for s in ('top', 'right'): ax.spines[s].set_visible(False)
panel_label(ax, 'D')

fig.suptitle('Supplementary Figure S5. Real 16S rRNA validation (73 taxa, corrected NCBI annotation)',
             fontsize=13, fontweight='bold', color=INK, y=0.985)
fig.tight_layout(rect=[0, 0, 1, 0.955])
fig.savefig('Supplementary_Figure_S5.pdf'); fig.savefig('Supplementary_Figure_S5.png', dpi=300)
print('S5 done | B %d tips, C %d tips' % (nB, nC))
