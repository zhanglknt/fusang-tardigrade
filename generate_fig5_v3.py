"""Figure 5 (v3): Scalability and deployment — 4 panels.
A: wall-clock vs n (Table 11) vs RAxML-NG (raxml_bench_results.csv), log-log
B: FastME vs NJ tree-building time (FASTME_BENCHMARK_REPORT.md)
C: peak RAM vs n (scalability_results.json + DCM n=10000)
D: deployment / cross-platform reproducibility summary
"""
import os
os.environ['MPLCONFIGDIR'] = 'matplotlib_cache'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np

INK, GREY = '#1a2332', '#64748b'
BLUE, TEAL, AMBER, RED, GREEN = '#2563eb', '#0d9488', '#d97706', '#dc2626', '#16a34a'

fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.6))
fig.patch.set_facecolor('white')

def panel_label(ax, s):
    ax.text(-0.14, 1.06, s, transform=ax.transAxes, fontsize=14,
            fontweight='bold', color=INK, va='top')

# ---------- A: wall-clock vs n ----------
ax = axes[0][0]
# Fusang, NJ tree builder (Table 11)
n_nj  = [20, 50, 100, 200]; t_nj = [4.9, 13.9, 24.5, 46.1]
# Fusang, FastME simplified + DCM (Table 11)
n_fm  = [500, 1000, 10000]; t_fm = [3.8, 5.2, 54.4]
# RAxML-NG (raxml_bench_results.csv; includes MAFFT)
n_rx  = [20, 50, 100, 200]; t_rx = [4.8, 39.8, 56.8, 244.1]
ax.plot(n_rx, t_rx, 'o--', color=GREY, lw=1.8, ms=6, label='MAFFT + RAxML-NG')
ax.plot(n_nj, t_nj, 's-', color=BLUE, lw=2, ms=6, label='Fusang (NJ builder)')
ax.plot(n_fm, t_fm, 's-', color=TEAL, lw=2.2, ms=6.5, label='Fusang (FastME / DCM)')
# O(n^2 log n) reference anchored at (10000, 54.4)
nn = np.array(n_fm, dtype=float)
ref = 54.4 * (nn**2 * np.log(nn)) / (10000.0**2 * np.log(10000.0))
ax.plot(nn, ref, ':', color=INK, lw=1.4, label='O(n$^2$ log n) reference')
ax.annotate('10,000 taxa\nin 54.4 s', xy=(10000, 54.4), xytext=(1400, 120),
            fontsize=8.5, color=TEAL, fontweight='bold',
            arrowprops=dict(arrowstyle='->', color=TEAL, lw=1.2))
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel('Number of taxa (n)', fontsize=9.5, color=INK)
ax.set_ylabel('Wall-clock time (s)', fontsize=9.5, color=INK)
ax.set_title('Wall-clock scaling, 4-core workstation', fontsize=10.5, color=INK, pad=6)
ax.legend(fontsize=8, frameon=False, loc='lower right')
ax.grid(alpha=0.25, which='both')
for s in ('top', 'right'): ax.spines[s].set_visible(False)
panel_label(ax, 'A')

# ---------- B: FastME vs NJ ----------
ax = axes[0][1]
n_b   = [200, 500, 1000]
nj_t  = [1.348, 13.199, 104.339]
fm_t  = [0.388, 0.872, 3.624]
spdup = [3.5, 15.1, 28.8]
x = np.arange(len(n_b)); w = 0.36
b1 = ax.bar(x - w/2, nj_t, w, color=GREY, label='NJ (BioPython), O(n$^3$)')
b2 = ax.bar(x + w/2, fm_t, w, color=TEAL, label='FastME BIONJ+BNNI, O(n$^2$ log n)')
for i, (xx, v) in enumerate(zip(x, spdup)):
    ax.text(xx, max(nj_t[i], fm_t[i]) * 1.6, f'{v}x', ha='center', fontsize=10,
            fontweight='bold', color=AMBER)
ax.set_yscale('log'); ax.set_ylim(0.2, 900)
ax.set_xticks(x); ax.set_xticklabels([f'n={n}' for n in n_b], fontsize=9.5)
ax.set_ylabel('Tree-building time (s)', fontsize=9.5, color=INK)
ax.set_title('FastME vs NJ tree building', fontsize=10.5, color=INK, pad=6)
ax.legend(fontsize=8, frameon=False, loc='upper left')
ax.grid(alpha=0.25, axis='y')
for s in ('top', 'right'): ax.spines[s].set_visible(False)
panel_label(ax, 'B')

# ---------- C: RAM vs n ----------
ax = axes[1][0]
n_r  = [200, 500, 1000, 2000, 5000, 10000]
ram  = [45, 78, 156, 312, 780, 609]
cols = [BLUE]*5 + [AMBER]
ax.bar(range(len(n_r)), ram, color=cols, width=0.62)
for i, v in enumerate(ram):
    ax.text(i, v + 18, f'{v}', ha='center', fontsize=9, color=INK)
ax.text(5.0, 720, 'DCM: 50 subtrees\nof 200 taxa each', ha='center', fontsize=7.5, color=AMBER)
ax.set_xticks(range(len(n_r)))
ax.set_xticklabels(['200', '500', '1k', '2k', '5k', '10k'], fontsize=9)
ax.set_ylim(0, 950)
ax.set_xlabel('Number of taxa (n)', fontsize=9.5, color=INK)
ax.set_ylabel('Peak RAM (MB)', fontsize=9.5, color=INK)
ax.set_title('Memory footprint (< 1 GB at n=10,000)', fontsize=10.5, color=INK, pad=6)
ax.grid(alpha=0.25, axis='y')
for s in ('top', 'right'): ax.spines[s].set_visible(False)
panel_label(ax, 'C')

# ---------- D: deployment panel ----------
ax = axes[1][1]
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis('off')
ax.set_title('Cross-platform deployment', fontsize=10.5, color=INK, pad=6)
panel_label(ax, 'D')
rows = [
    ('Bundled binaries', 'fastme.exe (Windows, 725 KB) + fastme_linux ship with the repo — no compilation'),
    ('Zero WSL dependency', 'Windows-native FastME (PE32+ x86-64) gives single-command operation on native Windows'),
    ('Bit-identical results', 'Windows-native binary produces bit-identical trees to the Linux build (FastME v2.1.6.4)'),
    ('Auto discovery', 'binary search cascade: bundled native exe > WSL install > system PATH'),
    ('Light stack', 'Python 3.9+ with NumPy / SciPy / Biopython / scikit-learn; MIT license + Zenodo DOI'),
]
y = 8.9
for head, body in rows:
    ax.text(0.35, y, '\u2713', fontsize=11, color=GREEN, fontweight='bold', va='center')
    ax.text(0.95, y, head, fontsize=9.5, color=INK, fontweight='bold', va='center')
    ax.text(0.95, y - 0.62, body, fontsize=8, color=GREY, va='center')
    y -= 1.75

fig.suptitle('Figure 5. Scalability and cross-platform deployment',
             fontsize=13.5, fontweight='bold', color=INK, y=0.985)
fig.tight_layout(rect=[0, 0, 1, 0.955])
fig.savefig('Figure5.pdf'); fig.savefig('Figure5.png', dpi=300)
print('Figure5 done')
