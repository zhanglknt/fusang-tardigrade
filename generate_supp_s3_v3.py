"""Supplementary Figure S3 (v3): DCM degradation — 4 panels.
A: single-seed pipeline stage trace (seed=42, n=200, clean)
B: simplified pipeline multi-seed validation vs full DCM
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
ax = axes[0][0]
stages = ['Simplified\npipeline', 'TF-IDF\nweighting', 'DCM+NJ\nrecovery', 'DCM+FastME\nrecovery', 'DCM+EPA\ngrafting']
vals = [0.005, 0.030, 0.005, 0.013, 0.388]
cols = [GREEN, AMBER, GREEN, GREEN, RED]
bars = ax.bar(stages, vals, color=cols, width=0.6)
for b, v in zip(bars, vals):
    ax.text(b.get_x() + b.get_width()/2, v + 0.012, f'{v:.3f}', ha='center',
            fontsize=9, fontweight='bold', color=INK)
ax.set_ylabel('nRF distance', fontsize=9.5, color=INK)
ax.set_ylim(0, 0.46)
ax.set_title('Stage-by-stage trace (n=200, seed=42, clean)', fontsize=10.5, color=INK, pad=6)
ax.tick_params(axis='x', labelsize=7.5)
ax.grid(alpha=0.25, axis='y')
for s in ('top', 'right'): ax.spines[s].set_visible(False)
panel_label(ax, 'A')

# ---------- B: simplified multi-seed vs DCM ----------
ax = axes[0][1]
names = ['Simplified\nsingle seed', 'Simplified\n10-seed mean', 'Full DCM\n(single seed)']
v = [0.005, 0.014, 0.388]
e = [0, 0.003, 0]
cols = [GREEN, TEAL, RED]
bars = ax.bar(names, v, yerr=e, capsize=5, color=cols, width=0.52,
              error_kw=dict(lw=1.2, ecolor=INK))
for b, vv in zip(bars, v):
    ax.text(b.get_x() + b.get_width()/2, vv + 0.014, f'{vv:.3f}', ha='center',
            fontsize=9.5, fontweight='bold', color=INK)
ax.set_ylabel('nRF distance', fontsize=9.5, color=INK)
ax.set_ylim(0, 0.46)
ax.set_title('Simplified pipeline: multi-seed validation\n(seeds 42-51, n=200, clean)', fontsize=10, color=INK, pad=6)
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
ax.text(9.0, 5.1, '\u2717 primary error source: grafting\nsubtrees with internal structure\n(nRF 0.005 \u2192 0.388)',
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
ax.grid(alpha=0.25)
for s in ('top', 'right'): ax.spines[s].set_visible(False)
panel_label(ax, 'D')

fig.suptitle('Supplementary Figure S3. DCM degradation trace and analysis',
             fontsize=13, fontweight='bold', color=INK, y=0.985)
fig.tight_layout(rect=[0, 0, 1, 0.955])
fig.savefig('Supplementary_Figure_S3.pdf'); fig.savefig('Supplementary_Figure_S3.png', dpi=300)
print('S3 done', pts)
