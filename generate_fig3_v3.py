"""Figure 3 (v3): Pipeline ablation and DCM degradation — 4 panels.
A: step-by-step nRF through pipeline stages (manuscript ablation values)
B: schematic — simplified vs DCM+EPA architectures
C: schematic — why EPA grafting fails for subtrees
D: adaptive threshold empirics (SIMPLE_THRESHOLD=500 basis)
"""
import os
os.environ['MPLCONFIGDIR'] = 'matplotlib_cache'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

INK, GREY = '#1a2332', '#64748b'
BLUE, TEAL, AMBER, RED, PURPLE = '#2563eb', '#0d9488', '#d97706', '#dc2626', '#7c3aed'

fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.6))
fig.patch.set_facecolor('white')

def panel_label(ax, s):
    ax.text(-0.14, 1.06, s, transform=ax.transAxes, fontsize=14,
            fontweight='bold', color=INK, va='top')

# ---------- A: ablation steps ----------
ax = axes[0][0]
steps = ['Simplified\nk-mer\u2192cosine\u2192NJ', '+TF-IDF\ntransform', '+DCM\n(NJ subtrees)', '+EPA\ngrafting']
vals  = [0.005, 0.030, 0.005, 0.388]
cols  = [TEAL, AMBER, TEAL, RED]
bars = ax.bar(steps, vals, color=cols, width=0.58)
for b, v in zip(bars, vals):
    ax.text(b.get_x() + b.get_width()/2, v + 0.012, f'{v:.3f}', ha='center',
            fontsize=10, fontweight='bold', color=INK)
ax.set_ylabel('nRF vs true tree (n=200, clean)', fontsize=9.5, color=INK)
ax.set_ylim(0, 0.45)
ax.set_title('Ablation: EPA grafting is the error source\n(illustrative single seed 42; multi-seed mean 0.014±0.003)',
             fontsize=9.5, color=INK, pad=6)
ax.tick_params(axis='x', labelsize=8.5)
ax.spines[['top', 'right']].set_visible(False)
panel_label(ax, 'A')

# ---------- B: architecture schematic ----------
ax = axes[0][1]
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis('off')
ax.set_title('Two pipelines, one adaptive switch', fontsize=10.5, color=INK, pad=6)
panel_label(ax, 'B')

def box(ax, x, y, w, h, fc, label, fs=9, tc='white'):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.08,rounding_size=0.18',
                                fc=fc, ec='none'))
    ax.text(x + w/2, y + h/2, label, ha='center', va='center', fontsize=fs,
            color=tc, fontweight='bold')

# simplified (top)
box(ax, 0.4, 7.6, 2.6, 1.1, TEAL, 'k-mer\nvectors', 8.5)
box(ax, 3.7, 7.6, 2.6, 1.1, TEAL, 'cosine\ndistance', 8.5)
box(ax, 7.0, 7.6, 2.6, 1.1, TEAL, 'NJ / FastME\n\u2192 tree', 8.5)
for x0 in (3.0, 6.3):
    ax.add_patch(FancyArrowPatch((x0 + 0.1, 8.15), (x0 + 0.65, 8.15),
                                 arrowstyle='-|>', mutation_scale=13, color=GREY, lw=1.8))
ax.text(5.0, 9.25, 'Simplified (default, n\u2264500): direct, no grafting', fontsize=9,
        color=TEAL, ha='center', fontweight='bold')
# DCM (bottom)
box(ax, 0.4, 3.6, 2.0, 1.1, BLUE, 'cluster\n(DCM)', 8.5)
box(ax, 2.9, 3.6, 2.0, 1.1, BLUE, 'backbone\ntree', 8.5)
box(ax, 5.4, 3.6, 2.0, 1.1, BLUE, 'subtrees', 8.5)
box(ax, 7.9, 3.6, 2.0, 1.1, RED, 'EPA\ngrafting', 8.5)
for x0 in (2.4, 4.9, 7.4):
    ax.add_patch(FancyArrowPatch((x0 + 0.1, 4.15), (x0 + 0.42, 4.15),
                                 arrowstyle='-|>', mutation_scale=13, color=GREY, lw=1.8))
ax.text(5.0, 5.35, 'DCM+EPA (n>500): scalable, grafting error contained', fontsize=9,
        color=BLUE, ha='center', fontweight='bold')
ax.text(5.0, 1.9, 'adaptive switch: SIMPLE_THRESHOLD = 500\n'
                  '(n\u2264500 simplified \u00b7 n>500 DCM)', fontsize=9.5, color=INK,
        ha='center', style='italic')

# ---------- C: EPA failure schematic ----------
ax = axes[1][0]
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis('off')
ax.set_title('Why EPA grafting fails for subtrees', fontsize=10.5, color=INK, pad=6)
panel_label(ax, 'C')
# backbone
ax.plot([1.2, 4.2, 7.2], [2.2, 5.4, 8.4], color=INK, lw=2.4)
ax.plot([4.2, 6.0], [5.4, 7.6], color=INK, lw=2.4)
ax.text(7.5, 8.6, 'backbone', fontsize=9, color=INK)
# correct attachment
ax.plot([4.2, 3.4], [5.4, 6.9], color=TEAL, lw=2.2)
ax.plot([3.4, 2.6], [6.9, 7.9], color=TEAL, lw=2.2)
ax.plot([3.4, 4.2], [6.9, 7.9], color=TEAL, lw=2.2)
ax.text(2.2, 8.6, 'true subtree\n(internal structure)', fontsize=8.5, color=TEAL)
# EPA misplacement
ax.plot([1.9, 2.55], [3.0, 4.15], color=RED, lw=2.2, ls='--')
ax.plot([2.55, 1.9], [4.15, 5.2], color=RED, lw=2.2, ls='--')
ax.plot([2.55, 3.2], [4.15, 5.2], color=RED, lw=2.2, ls='--')
ax.text(0.7, 5.6, 'EPA places subtree\nat wrong edge\n(designed for single\nshort reads)', fontsize=8.5, color=RED)
ax.text(5.0, 0.9, 'single-seed illustration: grafting entire structured subtrees\n'
                  'introduces topological errors that dominate the k-mer signal',
        fontsize=8.5, color=GREY, ha='center', va='top')

# ---------- D: threshold empirics ----------
ax = axes[1][1]
x = np.arange(2)
dcm  = [0.084, 0.088]
simp = [0.092, 0.107]
w = 0.34
b1 = ax.bar(x - w/2, dcm, w, color=BLUE, label='DCM+EPA pipeline')
b2 = ax.bar(x + w/2, simp, w, color=TEAL, label='simplified pipeline')
for bars in (b1, b2):
    for b in bars:
        ax.text(b.get_x() + b.get_width()/2, b.get_height() + 0.003,
                f'{b.get_height():.3f}', ha='center', fontsize=8.5, color=INK)
ax.set_xticks(x)
ax.set_xticklabels(['n=1000, 5 seeds 401\u2013405\n(DCM better \u2192 DCM for n>500)',
                    'threshold=2000 trial\n(regression \u2192 reverted)'], fontsize=8.5)
ax.set_xlim(-0.6, 1.6)
ax.set_ylim(0, 0.13)
ax.set_ylabel('nRF vs true tree', fontsize=9.5, color=INK)
ax.legend(fontsize=8.5, frameon=False, loc='upper left')
ax.set_title('Threshold evidence: DCM wins at n=1000;\nthreshold 2000 caused regression \u2192 reverted to 500',
             fontsize=9.5, color=INK, pad=6)
ax.spines[['top', 'right']].set_visible(False)
panel_label(ax, 'D')

fig.suptitle('Figure 3. Pipeline ablation and DCM degradation', fontsize=14,
             fontweight='bold', color=INK, y=0.985)
fig.tight_layout(rect=[0, 0, 1, 0.955])
fig.savefig('Figure3.pdf')
fig.savefig('Figure3.png', dpi=300)
print('Figure3 done')
