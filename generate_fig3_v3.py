"""Figure 3 (v3): Pipeline ablation and DCM degradation — 4 panels.

A: step-by-step nRF through pipeline stages (seed 42, n=200, indel rate=0.02)
B: schematic — simplified vs DCM+EPA architectures
C: schematic — why EPA grafting fails for subtrees
D: adaptive threshold empirics (SIMPLE_THRESHOLD=500) — reversed under the standard RF

Panel A values come from ``figure_data.ablation_seed42_indel()``
(repro_ablations/ablation_dcm_seed42_indel_nrfstandard.json, standard unrooted RF).
Panel D values come from ``figure_data.threshold_n1000()``
(table_threshold_n1000_nrfstandard.csv + _summary.json): 10 paired seeds 1000-1009 at
n=1000 clean, standard unrooted RF. It supersedes the withdrawn panel (DCM 0.083 vs
simplified 0.102, seeds 401-409), which used the root-dependent nRF and cannot be
rescored; under the standard metric the simplified pipeline is significantly better, so
the threshold is retained on runtime grounds only.
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

# ---------- D: threshold empirics (reversed under the standard RF) ----------
ax = axes[1][1]
th = FD.threshold_n1000()
ax.set_xlim(0.55, 3.25)
ax.set_ylim(0.205, 0.315)
for s, d in zip(th['nrf_simple'], th['nrf_dcm']):
    ax.plot([1, 2], [d, s], '-o', color=GREY, lw=0.9, ms=3.2, alpha=0.7, zorder=2)
ax.errorbar([1], [th['mean_dcm']], yerr=[th['sd_dcm']], fmt='s', ms=8, color=BLUE,
            capsize=4, lw=1.6, ecolor=BLUE, zorder=4)
ax.errorbar([2], [th['mean_simple']], yerr=[th['sd_simple']], fmt='s', ms=8, color=TEAL,
            capsize=4, lw=1.6, ecolor=TEAL, zorder=4)
ax.text(1, 0.240, f"DCM+EPA\n{th['mean_dcm']:.4f} \u00b1 {th['sd_dcm']:.4f}",
        ha='center', va='top', fontsize=8.0, color=BLUE, fontweight='bold')
ax.text(2.07, th['mean_simple'],
        f"simplified\n{th['mean_simple']:.4f} \u00b1 {th['sd_simple']:.4f}",
        ha='left', va='center', fontsize=8.0, color=TEAL, fontweight='bold')
ax.set_xticks([1, 2])
ax.set_xticklabels(['DCM+EPA', 'simplified'], fontsize=9)
ax.set_ylabel('nRF vs true tree', fontsize=9.5, color=INK)
ax.set_title('SIMPLE_THRESHOLD=500: accuracy basis withdrawn and reversed\n'
             f"n=1,000 clean, {th['n']} paired seeds \u2014 simplified "
             f"{th['simple_wins']}/{th['n']} wins",
             fontsize=9.2, color=INK, pad=6)
ax.text(0.97, 0.97,
        f"simplified {th['mean_simple']:.4f} \u00b1 {th['sd_simple']:.4f}\n"
        f"DCM+EPA  {th['mean_dcm']:.4f} \u00b1 {th['sd_dcm']:.4f}\n"
        f"\u0394 = {th['diff']:+.4f} nRF (simple \u2212 DCM)\n"
        f"Wilcoxon p={th['p']:.4f} (integer RF, zsplit)\n"
        f"paired d={th['d']:+.2f}, 95% CI [{th['ci'][0]:+.4f}, {th['ci'][1]:+.4f}]\n"
        f"runtime: DCM {th['rt_dcm']:.0f} s vs simplified {th['rt_simple']:.0f} s / seed\n"
        f"\u2192 threshold kept for ~{th['rt_simple'] / th['rt_dcm']:.1f}\u00d7 speed, "
        "not accuracy",
        transform=ax.transAxes, fontsize=6.9, color=INK, va='top', ha='right',
        bbox=dict(boxstyle='round,pad=0.32', fc='#f4f7fa', ec=GREY, lw=0.6))
ax.text(0.02, 0.015,
        'supersedes the earlier panel (DCM 0.083 vs simplified 0.102, seeds 401\u2013409), '
        'which\nused the root-dependent nRF and cannot be rescored',
        transform=ax.transAxes, fontsize=6.3, color=GREY, style='italic', va='bottom')
ax.grid(alpha=0.2, axis='y')
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
print('         D: n=1000 simplified %.4f+-%.4f vs DCM+EPA %.4f+-%.4f | %d/%d wins '
      'p=%.4f d=%+.2f CI [%+.4f, %+.4f] | runtime DCM %.0f s vs simple %.0f s'
      % (th['mean_simple'], th['sd_simple'], th['mean_dcm'], th['sd_dcm'],
         th['simple_wins'], th['n'], th['p'], th['d'], th['ci'][0], th['ci'][1],
         th['rt_dcm'], th['rt_simple']))
