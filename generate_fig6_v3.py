"""Figure 6 (v3): Real-data validation — 4 panels.
A: 16S rRNA tree (73 taxa, archived k=5,gap2 tree), tip labels coloured by phylum
B: taxonomic-rank signal (Table 9 displayed values + archived-rerun note)
C: Fusang v1 vs Tardigrade Edition comparison
D: gene-length DNA benchmark: locus length vs Fusang-FT2 agreement (gene_dna_results.json)
"""
import os
os.environ['MPLCONFIGDIR'] = 'matplotlib_cache'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
import json

INK, GREY = '#1a2332', '#64748b'
BLUE, TEAL, AMBER, RED, GREEN = '#2563eb', '#0d9488', '#d97706', '#dc2626', '#16a34a'
PURPLE, PINK = '#7c3aed', '#db2777'

fig = plt.figure(figsize=(10.4, 8.4))
fig.patch.set_facecolor('white')
gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1.0])

def panel_label(ax, s):
    ax.text(-0.12, 1.04, s, transform=ax.transAxes, fontsize=14,
            fontweight='bold', color=INK, va='top')

# ---------- A: 16S tree (custom cladogram, phylum-coloured tip dots) ----------
from Bio import Phylo
ax = fig.add_subplot(gs[0, 0])
tree = Phylo.read('real_data/real_16S_full_fusang_k5_gap2.nwk', 'newick')
meta = json.load(open('real_data/real_16S_metadata.json'))
TOP_PHYLA = ['Pseudomonadota', 'Bacillota', 'Actinomycetota',
             'Bacteroidota', 'Cyanobacteriota', 'Euryarchaeota']
PCOL = dict(zip(TOP_PHYLA, [BLUE, TEAL, AMBER, PURPLE, PINK, GREEN]))

tip_y = {}
for i, tip in enumerate(tree.get_terminals()):
    tip_y[tip] = i

def draw_clade(clade, x):
    """Recursive cladogram: x = depth, tips evenly spaced."""
    if clade.is_terminal():
        return tip_y[clade]
    ys = [draw_clade(c, x + 1) for c in clade.clades]
    yc = (ys[0] + ys[-1]) / 2.0
    for c, y in zip(clade.clades, ys):
        ax.plot([x, x + 1], [y, y], color=INK, lw=0.7, solid_capstyle='round')
        ax.plot([x, x], [min(y, yc), max(y, yc)], color=INK, lw=0.7)
    return yc

draw_clade(tree.root, 0)
maxd = max(len(tree.get_path(t)) for t in tree.get_terminals())  # edge-count depth
for tip, y in tip_y.items():
    ph = meta.get(tip.name, {}).get('phylum')
    ax.plot(maxd + 0.4, y, 'o', ms=3.2, color=PCOL.get(ph, GREY),
            markeredgecolor='none')
# highlight sister pairs cited in the manuscript
pairs = [('Escherichia_coli', 'Salmonella_enterica', 'E. coli / S. enterica'),
         ('Bacillus_subtilis', 'Geobacillus_kaustophilus', 'B. subtilis / G. kaustophilus')]
name2tip = {t.name: t for t in tip_y}
for a, b, lab in pairs:
    if a in name2tip and b in name2tip:
        y0, y1 = tip_y[name2tip[a]], tip_y[name2tip[b]]
        ym = (y0 + y1) / 2.0
        ax.plot([maxd + 0.9, maxd + 0.9], [y0, y1], color=RED, lw=1.1)
        ax.text(maxd + 1.2, ym, lab, fontsize=5.2, style='italic', color=RED, va='center')
ax.set_xlim(-0.5, maxd + 6.5)
ax.set_ylim(-1, len(tip_y))
ax.axis('off')
ax.set_title('16S rRNA tree, 73 taxa (k=5,gap2; 1.2 s, no alignment)',
             fontsize=10, color=INK, pad=6)
handles = [plt.Line2D([0], [0], marker='o', color='none', markerfacecolor=PCOL[p],
                      markersize=5, label=p) for p in TOP_PHYLA]
handles.append(plt.Line2D([0], [0], marker='o', color='none', markerfacecolor=GREY,
                          markersize=5, label='other (7 phyla)'))
ax.legend(handles=handles, fontsize=5.8, frameon=False, loc='lower right',
          bbox_to_anchor=(0.99, 0.01), handletextpad=0.15, borderaxespad=0.1,
          labelspacing=0.28)
panel_label(ax, 'A')

# ---------- B: taxonomic-rank signal ----------
ax = fig.add_subplot(gs[0, 1])
ranks  = ['Order', 'Phylum', 'Family']
vals   = [12.6, 4.6, -14.2]
notes  = ['p < 0.01', 'p < 0.05', 'n.s.']
cols   = [TEAL, TEAL, GREY]
bars = ax.bar(ranks, vals, color=cols, width=0.55)
for b, v, nt in zip(bars, vals, notes):
    yy = v + 1.2 if v >= 0 else v - 3.4
    ax.text(b.get_x() + b.get_width()/2, yy, f'{v:+.1f}%\n{nt}', ha='center',
            fontsize=8.5, color=INK, fontweight='bold')
ax.axhline(0, color=INK, lw=0.8)
ax.set_ylabel('Same-rank distance reduction (%)', fontsize=9.5, color=INK)
ax.set_ylim(-22, 22)
ax.set_title('Phylogenetic signal by taxonomic rank', fontsize=10.5, color=INK, pad=6)
ax.grid(alpha=0.25, axis='y')
for s in ('top', 'right'): ax.spines[s].set_visible(False)
ax.text(0.02, 0.03, 'Verification rerun on archived tree (73/74 taxa): order 7.3% (p=0.027),\n'
        'phylum 3.8% (p=0.003); family labels not retained in archive (original analysis).',
        transform=ax.transAxes, fontsize=6.8, color=GREY, style='italic')
panel_label(ax, 'B')

# ---------- C: v1 vs Tardigrade ----------
ax = fig.add_subplot(gs[1, 0])
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis('off')
ax.set_title('Fusang v1 (2023) vs Tardigrade Edition', fontsize=10.5, color=INK, pad=6)
panel_label(ax, 'C')
ax.add_patch(FancyBboxPatch((3.55, 8.6), 2.9, 0.95, boxstyle='round,pad=0.06,rounding_size=0.12',
                            fc='#e5e7eb', ec='none'))
ax.text(5.0, 9.07, 'Fusang v1', fontsize=9.5, fontweight='bold', color=GREY, ha='center', va='center')
ax.add_patch(FancyBboxPatch((6.65, 8.6), 3.2, 0.95, boxstyle='round,pad=0.06,rounding_size=0.12',
                            fc=TEAL, ec='none'))
ax.text(8.25, 9.07, 'Tardigrade', fontsize=9.5, fontweight='bold', color='white', ha='center', va='center')
rows = [
    ('Representation', 'neural net features', 'k-mer frequency vectors'),
    ('Scalability', '4-40 taxa', '10,000+ taxa'),
    ('Speed', 'GPU-dependent', 'CPU-only; 10k taxa < 1 min'),
    ('Deployment', 'Docker / cloud', 'single script + bundled binary'),
    ('Indel handling', 'not modelled', 'spaced k-mer sampling'),
]
y = 7.7
for name, v1, v2 in rows:
    ax.text(0.15, y, name, fontsize=8.5, fontweight='bold', color=INK, va='center')
    ax.text(3.7, y, v1, fontsize=8.2, color=GREY, va='center')
    ax.text(6.8, y, v2, fontsize=8.2, color=INK, va='center', fontweight='bold')
    ax.plot([0.1, 9.9], [y - 0.75, y - 0.75], color='#e5e7eb', lw=0.8)
    y -= 1.45

# ---------- D: gene-length scatter ----------
ax = fig.add_subplot(gs[1, 1])
d = json.load(open('real_data/afproject_genome/gene_dna/gene_dna_results.json'))
genes = d['genes']
names = list(genes.keys())
L = np.array([genes[g]['len_mean_bp'] for g in names])
agree = np.array([genes[g]['fusang_multik_vs_ft2_nrf'] for g in names])
ax.scatter(L, agree, s=42, color=TEAL, edgecolor=INK, lw=0.6, zorder=3)
for g, x, yv in zip(names, L, agree):
    ax.annotate(g, (x, yv), textcoords='offset points', xytext=(4, 3),
                fontsize=6.2, color=GREY)
ax.axvline(1000, color=AMBER, ls=':', lw=1.2)
short = L.argsort()[:3]
long_ = L >= 1000
m_short, m_long = agree[short].mean(), agree[long_].mean()
print(f'D: shortest-3 mean agreement = {m_short:.3f} (expect ~0.591); '
      f'>=1kb mean = {m_long:.3f} (expect ~0.373); overall {agree.mean():.3f} (expect 0.451)')
ax.text(0.03, 0.06, f'>=1 kb loci: mean nRF = {m_long:.3f}\n3 shortest loci: mean nRF = {m_short:.3f}',
        transform=ax.transAxes, fontsize=8, color=INK,
        bbox=dict(boxstyle='round,pad=0.3', fc='#f4f7fa', ec=GREY, lw=0.6))
ax.set_xlabel('Locus length (bp)', fontsize=9.5, color=INK)
ax.set_ylabel('nRF: Fusang multi-k vs FastTree2', fontsize=9, color=INK)
ax.set_title('Gene-length DNA benchmark (13 mtDNA loci, 25 fish)',
             fontsize=10, color=INK, pad=6)
ax.set_ylim(0.15, 0.85)
ax.grid(alpha=0.25)
for s in ('top', 'right'): ax.spines[s].set_visible(False)
panel_label(ax, 'D')

fig.suptitle('Figure 6. Real-data validation', fontsize=13.5, fontweight='bold',
             color=INK, y=0.985)
fig.tight_layout(rect=[0, 0, 1, 0.955])
fig.savefig('Figure6.pdf'); fig.savefig('Figure6.png', dpi=300)
print('Figure6 done')
