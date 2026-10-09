"""Figure 1 (v3): Indel robustness advantage — 4 panels.
A: nRF vs indel rate (Table 3 data, 112 seeds)
B: relative Fusang advantage over FT2 vs indel rate
C: conceptual schematic — spaced vs contiguous k-mers under an indel
D: method comparison at indel=0.02 (Fusang / FT2 / IQ-TREE2 / Mash)
"""
import os
os.environ['MPLCONFIGDIR'] = 'matplotlib_cache'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

INK, GREY = '#1a2332', '#64748b'
BLUE, TEAL, AMBER, RED, GREEN = '#2563eb', '#0d9488', '#d97706', '#dc2626', '#16a34a'

fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.6))
fig.patch.set_facecolor('white')

def panel_label(ax, s):
    ax.text(-0.14, 1.06, s, transform=ax.transAxes, fontsize=14,
            fontweight='bold', color=INK, va='top')

# ---------- A: nRF vs indel rate ----------
ax = axes[0][0]
rates = [0.005, 0.01, 0.02, 0.05]
fus   = [0.137, 0.107, 0.080, 0.066]
ft2   = [0.137, 0.112, 0.085, 0.076]
ax.plot(rates, ft2, 'o--', color=GREY, lw=1.8, ms=6, label='FastTree2 (MSA+ML)')
ax.plot(rates, fus, 's-', color=TEAL, lw=2.2, ms=6.5, label='Fusang (alignment-free)')
for x, y in zip(rates[1:], fus[1:]):
    ax.annotate(f'{y:.3f}', (x, y), textcoords='offset points', xytext=(0, -15),
                ha='center', fontsize=8, color=TEAL)
for x, y in zip(rates[1:], ft2[1:]):
    ax.annotate(f'{y:.3f}', (x, y), textcoords='offset points', xytext=(0, 8),
                ha='center', fontsize=8, color=GREY)
ax.annotate('0.137 (tie)', (rates[0], fus[0]), textcoords='offset points',
            xytext=(10, 6), ha='left', fontsize=8, color=INK)
ax.set_xlabel('Indel rate', fontsize=10, color=INK)
ax.set_ylabel('nRF vs true tree (mean, 112 seeds)', fontsize=9.5, color=INK)
ax.set_xticks(rates); ax.set_ylim(0.04, 0.16)
ax.legend(fontsize=8.5, frameon=False, loc='upper right')
ax.set_title('Accuracy ranking flips under indels', fontsize=10.5, color=INK, pad=6)
ax.spines[['top', 'right']].set_visible(False)
panel_label(ax, 'A')

# ---------- B: relative advantage ----------
ax = axes[0][1]
adv = [(t - f) / t * 100 for f, t in zip(fus, ft2)]
bars = ax.bar([str(r) for r in rates], adv, width=0.58,
              color=[GREY if a == 0 else TEAL for a in adv], edgecolor='none')
for b, a in zip(bars, adv):
    ax.text(b.get_x() + b.get_width()/2, a + 0.35, f'{a:.1f}%', ha='center',
            fontsize=9.5, fontweight='bold', color=INK)
ax.axhline(0, color=GREY, lw=0.8)
ax.set_xlabel('Indel rate', fontsize=10, color=INK)
ax.set_ylabel('Relative advantage over FastTree2 (%)', fontsize=9.5, color=INK)
ax.set_ylim(0, 15.5)
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
names = ['Fusang', 'FastTree2', 'IQ-TREE2\nGTR', 'Mash\n(MinHash)']
vals  = [0.080, 0.085, 0.147, 0.762]
sds   = [0.016, 0.025, 0.027, 0.036]
ns    = ['n=112', 'n=112', 'n=121', 'n=30']
cols  = [TEAL, GREY, GREY, RED]
bars = ax.bar(names, vals, yerr=sds, capsize=4, color=cols, width=0.6,
              error_kw=dict(lw=1.2, ecolor=INK))
for b, v, n in zip(bars, vals, ns):
    ax.text(b.get_x() + b.get_width()/2, v + 0.055, f'{v:.3f}', ha='center',
            fontsize=9, fontweight='bold', color=INK)
    ax.text(b.get_x() + b.get_width()/2, 0.02, n, ha='center', fontsize=8, color='white',
            fontweight='bold')
ax.set_ylabel('nRF vs true tree (mean ± SD)', fontsize=9.5, color=INK)
ax.set_title('Indel rate 0.02: method comparison', fontsize=10.5, color=INK, pad=6)
ax.set_ylim(0, 0.88)
ax.spines[['top', 'right']].set_visible(False)
panel_label(ax, 'D')

fig.suptitle('Figure 1. Indel robustness advantage', fontsize=14, fontweight='bold',
             color=INK, y=0.985)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig('Figure1.pdf')
fig.savefig('Figure1.png', dpi=300)
print('Figure1 done')
