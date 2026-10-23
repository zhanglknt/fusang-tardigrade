"""Supplementary Figure S3 (v3): DCM degradation — 4 panels.
A: single-seed pipeline stage trace (seed=42, n=200, indel rate=0.02;
   machine-readable artefact: repro_ablations/ablation_dcm_seed42_indel_nrfstandard.json)
B: DCM max_group_size sweep on the same seed (11 settings) vs simplified pipeline
   and MAFFT+FastTree2 reference lines (same JSON artefact)
C: DCM pipeline schematic with EPA grafting failure point
D: nRF stability across scales incl. DCM regime (scalability_results.json, 5 seeds)
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

fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.6))
fig.patch.set_facecolor('white')

def panel_label(ax, s):
    ax.text(-0.14, 1.06, s, transform=ax.transAxes, fontsize=14,
            fontweight='bold', color=INK, va='top')

# ---------- A: single-seed stage trace ----------
# Standard root-independent unrooted RF, n=200 indel rate=0.02 seed 42. All values are
# read directly from the machine-readable artefact
# repro_ablations/ablation_dcm_seed42_indel_nrfstandard.json (no transcription).
ABL = json.load(open('repro_ablations/ablation_dcm_seed42_indel_nrfstandard.json'))
ax = axes[0][0]
stages = ['Simplified\npipeline (NJ)', 'Simplified\n(FastME)', 'distances.py\nmatrix (NJ)',
          'TF-IDF\n(diagnostic)', 'DCM+EPA\ngrafting']
vals = [ABL['step1_simplified_nj']['nrf'], ABL['step2_fastme']['nrf'],
        ABL['step3_v2_distmat_nj']['nrf'], ABL['step4_tfidf_nj']['nrf'],
        ABL['step5_dcm_sweep_nj_default']['mg10_nj_default']['nrf']]
cols = [GREEN, GREEN, GREEN, AMBER, RED]
bars = ax.bar(stages, vals, color=cols, width=0.6)
for b, v in zip(bars, vals):
    ax.text(b.get_x() + b.get_width()/2, v + 0.012, f'{v:.3f}', ha='center',
            fontsize=9, fontweight='bold', color=INK)
ax.set_ylabel('nRF distance', fontsize=9.5, color=INK)
ax.set_ylim(0, 0.62)
ax.set_title('Stage-by-stage trace (n=200, indel 0.02, seed=42;\nstandard root-independent unrooted RF)',
             fontsize=9.8, color=INK, pad=6)
ax.tick_params(axis='x', labelsize=7.5)
ax.grid(alpha=0.25, axis='y')
for s in ('top', 'right'): ax.spines[s].set_visible(False)
panel_label(ax, 'A')

# ---------- B: DCM max_group_size sweep on the same seed ----------
# All values from the same JSON artefact (step5_dcm_sweep_nj_default, 11 settings,
# NJ subtrees); reference lines: simplified pipeline and MAFFT+FastTree2 on this seed.
ax = axes[0][1]
sweep = ABL['step5_dcm_sweep_nj_default']
mg_labels = ['10', '20', '25', '30', '40', '50', '60', '80', '100', '120', '150']
mg_vals = [sweep[f'mg{m}_nj_default']['nrf'] for m in mg_labels]
cols = [RED if i == 0 else AMBER for i in range(len(mg_vals))]
bars = ax.bar(mg_labels, mg_vals, color=cols, width=0.62)
for b, vv in zip(bars, mg_vals):
    ax.text(b.get_x() + b.get_width()/2, vv + 0.010, f'{vv:.2f}', ha='center',
            fontsize=7.6, fontweight='bold', color=INK)
ax.axhline(ABL['step1_simplified_nj']['nrf'], color=GREEN, ls='--', lw=1.6,
           label='simplified pipeline (0.081)')
ax.axhline(ABL['fasttree2_mafft']['nrf'], color=BLUE, ls=':', lw=1.6,
           label='MAFFT+FastTree2 (0.096)')
ax.legend(loc='upper right', fontsize=7.2, frameon=False, handlelength=2.2)
ax.set_xlabel('max_group_size (NJ subtrees)', fontsize=9, color=INK)
ax.set_ylabel('nRF distance', fontsize=9.5, color=INK)
ax.set_ylim(0, 0.62)
ax.set_title('DCM max_group_size sweep (same seed, 11 settings):\nevery setting loses to the simplified pipeline',
             fontsize=9.4, color=INK, pad=6)
ax.tick_params(axis='x', labelsize=8)
ax.grid(alpha=0.25, axis='y')
for s in ('top', 'right'): ax.spines[s].set_visible(False)
panel_label(ax, 'B')

# ---------- C: DCM schematic with failure point ----------
ax = axes[1][0]
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis('off')
ax.set_title('Where DCM loses accuracy', fontsize=10.5, color=INK, pad=6)
panel_label(ax, 'C')

def node(x, y, w, h, text, fc, tc='white', fs=8):
    ax.add_patch(FancyBboxPatch((x - w/2, y - h/2), w, h,
                                boxstyle='round,pad=0.08,rounding_size=0.15', fc=fc, ec='none'))
    ax.text(x, y, text, fontsize=fs, color=tc, ha='center', va='center', fontweight='bold')

steps = [('pairwise\nk-mer cosine\nmatrix', TEAL), ('UPGMA\nbalanced split\n(\u2264200/group)', TEAL),
         ('centroid\nbackbone\ntree', TEAL), ('per-group\nNJ\nsubtrees', TEAL),
         ('EPA\ngrafting', RED)]
x = 1.0
for i, (t, c) in enumerate(steps):
    node(x, 6.6, 1.55, 1.5, t, c, fs=6.8)
    if i < len(steps) - 1:
        ax.annotate('', xy=(x + 1.22, 6.6), xytext=(x + 0.80, 6.6),
                    arrowprops=dict(arrowstyle='-|>', color=INK, lw=1.4))
    x += 2.0
ax.text(9.0, 5.1, '\u2717 primary error source: grafting\nsubtrees with internal structure\n(nRF 0.081 \u2192 0.523)',
        fontsize=7.5, color=RED, ha='center', va='top', fontweight='bold')
# recovery path
ax.annotate('', xy=(5.0, 3.1), xytext=(8.6, 5.7),
            arrowprops=dict(arrowstyle='-|>', color=GREEN, lw=1.6, linestyle='--'))
node(5.0, 2.2, 6.6, 1.35, 'recovery: skip EPA \u2014 direct NJ on full matrix\n(simplified pipeline, n \u2264 500)', GREEN, fs=7.8)
ax.text(5.0, 0.8, 'DCM retained only where it is needed: n > 500 (scalability + small accuracy gain at scale)',
        fontsize=8, color=GREY, ha='center', style='italic')

# ---------- D: nRF across scales incl. DCM regime ----------
ax = axes[1][1]
d = json.load(open('scalability_results.json'))
by_n = {}
for r in d:
    if isinstance(r, dict) and 'nrf_mean' in r:
        key = r['n']
        if key not in by_n or r.get('pipeline'):
            by_n[key] = r
pts = [(r['n'], r['nrf_mean'], r['nrf_std'], r.get('pipeline', '')) for r in by_n.values()]
pts.sort()
ns = [p[0] for p in pts]; m = [p[1] for p in pts]; sd = [p[2] for p in pts]
cols = [TEAL if 'DCM' not in p[3] else AMBER for p in pts]
ax.errorbar(ns, m, yerr=sd, fmt='o', color=INK, ecolor=GREY, elinewidth=1.1,
            capsize=3.5, ms=6.5, mfc='white', mew=1.6, zorder=3)
for x, y, c in zip(ns, m, cols):
    ax.plot(x, y, 'o', color=c, ms=5, zorder=4)
ax.axvspan(1500, 2500, color=AMBER, alpha=0.08)
ax.text(1900, 0.318, 'DCM regime\n(n>500)', fontsize=7.5, color=AMBER, ha='center')
ax.set_xscale('log'); ax.set_xticks(ns)
ax.set_xticklabels([str(n) for n in ns], fontsize=8.5)
ax.set_xlabel('Dataset size n (log scale)', fontsize=9.5, color=INK)
ax.set_ylabel('nRF vs true tree (mean \u00b1 SD, 5 seeds)', fontsize=9, color=INK)
ax.set_ylim(0.26, 0.33)
ax.set_title('DCM holds accuracy at scale\n(clean coalescent, vs true tree)', fontsize=10, color=INK, pad=6)
ax.text(0.02, 0.02, 'legacy JC69 coalescent data + superseded nRF\nimplementation — shown for shape only',
        transform=ax.transAxes, fontsize=6.6, color=GREY, style='italic', va='bottom')
ax.grid(alpha=0.25)
for s in ('top', 'right'): ax.spines[s].set_visible(False)
panel_label(ax, 'D')

fig.suptitle('Supplementary Figure S3. DCM degradation trace and analysis',
             fontsize=13, fontweight='bold', color=INK, y=0.985)
fig.tight_layout(rect=[0, 0, 1, 0.955])
fig.savefig('Supplementary_Figure_S3.pdf'); fig.savefig('Supplementary_Figure_S3.png', dpi=300)
print('S3 done', pts)
