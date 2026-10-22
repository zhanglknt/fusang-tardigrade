"""Figure 3 (v3): Pipeline ablation and DCM degradation — 4 panels.

A: step-by-step nRF through pipeline stages (seed 42, n=200, indel rate=0.02)
B: schematic — simplified vs DCM+EPA architectures
C: schematic — why EPA grafting fails for subtrees
D: adaptive threshold empirics (SIMPLE_THRESHOLD=500 basis)

Panel A values come from ``figure_data.ablation_seed42_indel()``
(repro_ablations/ablation_dcm_seed42_indel_nrfstandard.json, standard unrooted RF).
Panel D quotes the superseded-metric n=1000 comparison (0.083±0.005 vs 0.102); the true
trees for seeds 401-409 are no longer in the working tree, so the panel carries the
disclosure required by the manuscript's provenance note (see Results and Figure 3 legend).
"""
import os
os.environ['MPLCONFIGDIR'] = 'matplotlib_cache'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

import figure_data as FD

INK, GREY = '#1a2332', '#64748b'
BLUE, TEAL, AMBER, RED, PURPLE = '#2563eb', '#0d9488', '#d97706', '#dc2626', '#7c3aed'

fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.6))
fig.patch.set_facecolor('white')

def panel_label(ax, s):
    ax.text(-0.14, 1.06, s, transform=ax.transAxes, fontsize=14,
            fontweight='bold', color=INK, va='top')

# ---------- A: ablation steps (n=200, indel rate=0.02, seed 42) ----------
abl = FD.ablation_seed42_indel()
ax = axes[0][0]
steps = [c[0] for c in abl['chain']]
vals = [c[1] for c in abl['chain']]
cols = [dict(green=TEAL, amber=AMBER, red=RED)[c[2]] for c in abl['chain']]
bars = ax.bar(steps, vals, color=cols, width=0.58)
for b, v in zip(bars, vals):
    ax.text(b.get_x() + b.get_width()/2, v + 0.014, f'{v:.3f}', ha='center',
            fontsize=9.5, fontweight='bold', color=INK)
ax.set_ylabel('nRF vs true tree (n=200, indel rate=0.02)', fontsize=9, color=INK)
ax.set_ylim(0, 0.62)
ax.set_title('Ablation chain, seed 42: EPA grafting is the error source\n'
             f'(11 DCM configs mg=10–150 span {abl["sweep_min"]:.3f}–{abl["sweep_max"]:.3f})',
             fontsize=9.5, color=INK, pad=6)
ax.tick_params(axis='x', labelsize=7.6)
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
# n=1000 (seeds 401-409): DCM 0.083 +/- 0.005 vs simplified 0.102 on seed 401.
# These values were computed under the superseded root-dependent nRF implementation and
# cannot be rescored (the seed 401-409 true trees are no longer in the working tree);
# they are shown here only to document the scalability decision, per the manuscript's
# provenance note (Results, "simplified pipeline" section, and Figure 3 legend).
x = [0]
w = 0.32
b1 = ax.bar([x[0] - w/2], [0.083], w, yerr=[0.005], capsize=3, color=BLUE,
            label='DCM+EPA (9 seeds 401\u2013409)', error_kw=dict(lw=1.1, ecolor=INK))
b2 = ax.bar([x[0] + w/2], [0.102], w, color=TEAL,
            label='simplified (seed 401 only)')
ax.text(x[0] - w/2, 0.083 + 0.008, '0.083\n±0.005', ha='center', fontsize=8.5,
        color=INK, fontweight='bold')
ax.text(x[0] + w/2, 0.102 + 0.006, '0.102', ha='center', fontsize=8.5,
        color=INK, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(['n=1000, indel (seed set 401\u2013409)'], fontsize=8.5)
ax.set_xlim(-0.7, 0.7)
ax.set_ylim(0, 0.20)
ax.set_ylabel('nRF vs true tree', fontsize=9.5, color=INK)
ax.legend(fontsize=8.2, frameon=False, loc='upper left', handlelength=1.4)
ax.set_title('Empirical basis of SIMPLE_THRESHOLD=500:\nDCM+EPA adopted at n=1000 scale',
             fontsize=9.5, color=INK, pad=6)
ax.text(0.5, 0.82, 'superseded nRF implementation (seed 401\u2013409 true trees\n'
                   'no longer in the working tree); shown to document scalability,\n'
                   'not as an accuracy result. threshold=2000 trial degraded\n'
                   'accuracy at n=1000 and was reverted to 500.',
        transform=ax.transAxes, fontsize=6.6, color=GREY, style='italic',
        va='top', ha='center')
ax.spines[['top', 'right']].set_visible(False)
panel_label(ax, 'D')

fig.suptitle('Figure 3. Pipeline ablation and DCM degradation', fontsize=14,
             fontweight='bold', color=INK, y=0.985)
fig.tight_layout(rect=[0, 0, 1, 0.955])
fig.savefig('Figure3.pdf')
fig.savefig('Figure3.png', dpi=300)
print('Figure3 done | chain:',
      ' -> '.join('%.3f' % v for v in vals),
      '| sweep %.3f-%.3f (mean %.3f +- %.3f, n=%d configs)'
      % (abl['sweep_min'], abl['sweep_max'], abl['sweep_mean'], abl['sweep_sd'],
         abl['sweep_n']))
