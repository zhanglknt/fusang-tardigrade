"""Figure 1 (v3): Indel robustness advantage — 4 panels.

A: nRF vs indel rate (manuscript Table 3)
B: relative Fusang advantage over FT2 vs indel rate
C: conceptual schematic — spaced vs contiguous k-mers under an indel
D: method comparison at indel=0.02 (Fusang / FT2 / IQ-TREE2 / Mash)

Numbers never appear in this file: everything comes from ``figure_data.py``, which reads
  - indelible_sweep/indel{0.005,0.01,0.05}/seed*_indel_{true,fusang,ft2}.nwk  (A, B)
    rescored on the fly with nrf_standard.py (the stored sweep_indel*.csv are superseded);
  - indelible_n200_indel/bench_n200_indel_130seeds_nrfstandard.csv           (D, Fusang + FT2)
  - indelible_iqtree2/iqtree2_n200_indel_nrfstandard.csv                     (D, IQ-TREE2)
  - bench_mash_indelible_n200_nrfstandard.csv                                (D, Mash)
"""
import os
os.environ['MPLCONFIGDIR'] = 'matplotlib_cache'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import figure_data as FD

INK, GREY = '#1a2332', '#64748b'
BLUE, TEAL, AMBER, RED, GREEN = '#2563eb', '#0d9488', '#d97706', '#dc2626', '#16a34a'

fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.6))
fig.patch.set_facecolor('white')


def panel_label(ax, s):
    ax.text(-0.14, 1.06, s, transform=ax.transAxes, fontsize=14,
            fontweight='bold', color=INK, va='top')


# ---------- A: nRF vs indel rate (Table 3) ----------
scan = FD.indel_rate_scan()
rates = sorted(scan)
ax = axes[0][0]
fus = [scan[r]['mean_fus'] for r in rates]
ft2 = [scan[r]['mean_ft2'] for r in rates]
ax.plot(rates, ft2, 'o--', color=GREY, lw=1.8, ms=6, label='FastTree2 (MSA+ML)')
ax.plot(rates, fus, 's-', color=TEAL, lw=2.2, ms=6.5, label='Fusang (simplified pipeline, k=5,gap2)')
for i, (x, y) in enumerate(zip(rates, fus)):
    off = (12, 8) if i == 0 else (-9, -15)
    ha = 'left' if i == 0 else 'right'
    ax.annotate(f'{y:.3f}', (x, y), textcoords='offset points', xytext=off,
                ha=ha, fontsize=8, color=TEAL)
for x, y in zip(rates, ft2):
    ax.annotate(f'{y:.3f}', (x, y), textcoords='offset points', xytext=(9, 7),
                ha='left', fontsize=8, color=GREY)
ax.set_xscale('log')
ax.set_xticks(rates)
ax.set_xticklabels([str(r) for r in rates], fontsize=9)
ax.set_ylim(0.165, 0.235)
ax.set_xlabel('Indel rate', fontsize=10, color=INK)
ax.set_ylabel(f'nRF vs true tree (mean, {scan[rates[0]]["n"]} seeds)', fontsize=9.5, color=INK)
ax.legend(fontsize=8.5, frameon=False, loc='upper right')
ax.set_title('Accuracy ranking flips under indels', fontsize=10.5, color=INK, pad=6)
ax.spines[['top', 'right']].set_visible(False)
panel_label(ax, 'A')

# ---------- B: relative advantage ----------
ax = axes[0][1]
adv = [scan[r]['advantage_pct'] for r in rates]
bars = ax.bar([str(r) for r in rates], adv, width=0.58,
              color=[GREY if a < 0 else TEAL for a in adv], edgecolor='none')
for b, a, r in zip(bars, adv, rates):
    off = 0.45 if a >= 0 else -1.35
    ax.text(b.get_x() + b.get_width() / 2, a + off, f'{a:+.1f}%', ha='center',
            fontsize=9.5, fontweight='bold', color=INK)
sig = {0.005: 'FT2 better', 0.01: 'FT2 better', 0.02: 'borderline', 0.05: 'Fusang better'}
for b, a, r in zip(bars, adv, rates):
    off = -2.0 if a >= 0 else 0.5
    ax.text(b.get_x() + b.get_width() / 2, off, sig[r], ha='center',
            fontsize=7.2, color=GREY, style='italic')
ax.axhline(0, color=GREY, lw=0.8)
ax.set_xlabel('Indel rate', fontsize=10, color=INK)
ax.set_ylabel('Relative Fusang advantage over FastTree2 (%)', fontsize=9, color=INK)
ax.set_ylim(-13, 15.5)
ax.set_title('Advantage grows monotonically with indel rate', fontsize=10.5, color=INK, pad=6)
ax.spines[['top', 'right']].set_visible(False)
panel_label(ax, 'B')

# ---------- C: conceptual schematic ----------
ax = axes[1][0]
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis('off')
ax.set_title('Why spaced k-mers tolerate indels', fontsize=10.5, color=INK, pad=6)
panel_label(ax, 'C')
# reference sequence
ax.add_patch(plt.Rectangle((0.6, 7.6), 8.8, 1.0, fc='#e2e8f0', ec=GREY, lw=1))
ax.text(0.6, 9.0, 'Reference sequence', fontsize=9, color=INK, fontweight='bold')
# indel sequence (with insertion)
ax.add_patch(plt.Rectangle((0.6, 4.6), 4.6, 1.0, fc='#e2e8f0', ec=GREY, lw=1))
ax.add_patch(plt.Rectangle((5.2, 4.6), 1.2, 1.0, fc='#fecaca', ec=RED, lw=1.4))
ax.add_patch(plt.Rectangle((6.4, 4.6), 3.0, 1.0, fc='#e2e8f0', ec=GREY, lw=1))
ax.text(0.6, 6.0, 'Sequence with 4-bp insertion', fontsize=9, color=INK, fontweight='bold')
ax.text(5.8, 5.1, 'indel', fontsize=8.5, color=RED, ha='center', fontweight='bold')
# contiguous k-mer (broken)
ax.add_patch(plt.Rectangle((3.4, 3.3), 3.6, 0.8, fc='none', ec=RED, lw=2.2))
ax.text(5.2, 2.7, 'contiguous 5-mer: frameshifted, counts disrupted', fontsize=8.5,
        color=RED, ha='center')
# spaced k-mer (skips)
for x0 in (1.0, 3.4, 6.8):
    ax.add_patch(plt.Rectangle((x0, 7.0), 1.0, 0.5, fc=GREEN, ec='none', alpha=0.85))
ax.text(5.0, 1.6, 'spaced 5-mer (1 0 0 1 0 0 1): samples anchor positions on both\n'
                  'sides of the indel — most spaced patterns survive length shifts',
        fontsize=8.5, color=GREEN, ha='center', va='top')
ax.annotate('', xy=(3.9, 7.5), xytext=(3.9, 4.15),
            arrowprops=dict(arrowstyle='-|>', color=GREEN, lw=2))

# ---------- D: methods at indel=0.02 ----------
ax = axes[1][1]
meth = FD.method_comparison_indel002()
cmap = {'teal': TEAL, 'grey': GREY, 'red': RED}
names = [m['label'] for m in meth]
vals = [m['mean'] for m in meth]
sds = [m['sd'] for m in meth]
cols = [cmap[m['color']] for m in meth]
bars = ax.bar(names, vals, yerr=sds, capsize=4, color=cols, width=0.6,
              error_kw=dict(lw=1.2, ecolor=INK))
top = max(v + s for v, s in zip(vals, sds))
for b, m in zip(bars, meth):
    ax.text(b.get_x() + b.get_width() / 2, m['mean'] + m['sd'] + 0.022,
            f'{m["mean"]:.3f}', ha='center', fontsize=9, fontweight='bold', color=INK)
    ax.text(b.get_x() + b.get_width() / 2, 0.014, f'n={m["n"]}', ha='center',
            fontsize=8, color='white', fontweight='bold')
ax.set_ylabel('nRF vs true tree (mean ± SD)', fontsize=9.5, color=INK)
ax.set_title('Indel rate 0.02: method comparison', fontsize=10.5, color=INK, pad=6)
ax.set_ylim(0, top * 1.18)
ax.spines[['top', 'right']].set_visible(False)
panel_label(ax, 'D')

fig.suptitle('Figure 1. Indel robustness advantage', fontsize=14, fontweight='bold',
             color=INK, y=0.985)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig('Figure1.pdf')
fig.savefig('Figure1.png', dpi=300)
print('Figure1 done |',
      ' '.join('%s %.4f+-%.4f(n=%d)' % (m['key'], m['mean'], m['sd'], m['n']) for m in meth))
print('            | advantage', ' '.join('%s %+.1f%%' % (r, scan[r]['advantage_pct'])
                                          for r in rates))
