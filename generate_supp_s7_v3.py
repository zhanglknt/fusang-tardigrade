"""Supplementary Figure S7 (v3): multi-k ensemble detail — 4 panels.
A: per-seed scatter, default spaced vs ensemble (seed set B, 30 seeds)
B: accuracy comparison across configurations (seed set B, mean +/- SD)
C: sorted per-seed differences, set B (d=+0.53, p=0.006)
D: sorted per-seed differences, set A (27 seeds, d=-0.22, p=0.30) — instability
"""
import os
os.environ['MPLCONFIGDIR'] = 'matplotlib_cache'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import csv, json
from scipy.stats import wilcoxon

INK, GREY = '#1a2332', '#64748b'
BLUE, TEAL, AMBER, RED = '#2563eb', '#0d9488', '#d97706', '#dc2626'

rows = list(csv.DictReader(open('benchmark_multik_ensemble_n200_indel.csv')))
defB = np.array([float(r['nrf_fusang_original']) for r in rows])
mulB = np.array([float(r['nrf_multik_ensemble']) for r in rows])
k5 = np.array([float(r['nrf_k5_contig']) for r in rows])
k7 = np.array([float(r['nrf_k7_contig']) for r in rows])
k9 = np.array([float(r['nrf_k9_contig']) for r in rows])

d7 = json.load(open('table7_recomputed_recovered.json'))
resA = d7['results'] if isinstance(d7, dict) and 'results' in d7 else d7
if isinstance(resA, dict):
    seedsA = list(resA.keys())
    first = resA[seedsA[0]]
    fk = [k for k in first if 'fusang' in k.lower() or 'default' in k.lower() or 'spaced' in k.lower()][0]
    mk = [k for k in first if 'multi' in k.lower() or 'ensemble' in k.lower()][0]
    defA = np.array([resA[s][fk] for s in seedsA], dtype=float)
    mulA = np.array([resA[s][mk] for s in seedsA], dtype=float)
else:
    defA = np.array([r['fusang'] for r in resA], dtype=float)
    mulA = np.array([r['multik'] for r in resA], dtype=float)

def paired_d(a, b):
    dd = a - b
    return dd.mean() / dd.std(ddof=1)

dB, pB = paired_d(defB, mulB), wilcoxon(defB, mulB).pvalue
dA, pA = paired_d(defA, mulA), wilcoxon(defA, mulA).pvalue
impr = (defB > mulB).sum()
print(f'B: n=30 d={dB:+.3f} p={pB:.4f} improved {impr}/30 | '
      f'means orig {defB.mean():.3f} ens {mulB.mean():.3f} '
      f'k5 {k5.mean():.3f} k7 {k7.mean():.3f} k9 {k9.mean():.3f}')
print(f'A: n=27 d={dA:+.3f} p={pA:.4f}')

fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.6))
fig.patch.set_facecolor('white')

def panel_label(ax, s):
    ax.text(-0.14, 1.06, s, transform=ax.transAxes, fontsize=14,
            fontweight='bold', color=INK, va='top')

# ---------- A: scatter ----------
ax = axes[0][0]
lim = [0.045, 0.155]
ax.plot(lim, lim, '--', color=RED, lw=1.2, alpha=0.7, label='y = x (no change)')
ax.scatter(defB, mulB, s=30, color=TEAL, edgecolor=INK, lw=0.5, zorder=3)
below = (mulB < defB).sum()
ax.text(0.052, 0.135, f'ensemble better: {below}/30 seeds\nmean improvement '
        f'{np.mean(defB - mulB):+.3f}', fontsize=8.5, color=INK,
        bbox=dict(boxstyle='round,pad=0.3', fc='#f4f7fa', ec=GREY, lw=0.6))
ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel('Default spaced k=5,gap2 nRF', fontsize=9.5, color=INK)
ax.set_ylabel('Multi-k ensemble (k=5,7,9) nRF', fontsize=9, color=INK)
ax.set_title('Per-seed comparison (seed set B, n=30)', fontsize=10.5, color=INK, pad=6)
ax.legend(fontsize=8, frameon=False, loc='lower right')
ax.grid(alpha=0.25)
for s in ('top', 'right'): ax.spines[s].set_visible(False)
panel_label(ax, 'A')

# ---------- B: config bars ----------
ax = axes[0][1]
names = ['default\nspaced\nk=5,gap2', 'ensemble\n(k=5,7,9)', 'k=5\ncontig', 'k=7\ncontig', 'k=9\ncontig']
vals = [defB, mulB, k5, k7, k9]
means = [v.mean() for v in vals]; sds = [v.std(ddof=1) for v in vals]
cols = [GREY, TEAL, BLUE, BLUE, BLUE]
bars = ax.bar(names, means, yerr=sds, capsize=4, color=cols, width=0.62,
              error_kw=dict(lw=1.1, ecolor=INK))
for b, m in zip(bars, means):
    ax.text(b.get_x() + b.get_width()/2, b.get_height() + 0.012, f'{m:.3f}',
            ha='center', fontsize=8.5, fontweight='bold', color=INK)
ax.set_ylabel('nRF vs true tree (mean \u00b1 SD)', fontsize=9.5, color=INK)
ax.set_ylim(0, 0.16)
ax.set_title('Configuration comparison (seed set B)', fontsize=10.5, color=INK, pad=6)
ax.tick_params(axis='x', labelsize=7.5)
ax.grid(alpha=0.25, axis='y')
for s in ('top', 'right'): ax.spines[s].set_visible(False)
panel_label(ax, 'B')

# ---------- C/D: sorted diffs ----------
for ax, dd, d, p, lab, n in ((axes[1][0], defB - mulB, dB, pB, 'seed set B', 30),
                             (axes[1][1], defA - mulA, dA, pA, 'seed set A', 27)):
    order = np.argsort(dd)
    x = np.arange(1, len(dd) + 1)
    cols = [TEAL if v > 0 else RED for v in dd[order]]
    ax.bar(x, dd[order], width=1.0, color=cols, edgecolor='none')
    ax.axhline(0, color=INK, lw=0.8)
    wins = (dd > 0).sum()
    ax.set_xlabel('Seed (sorted)', fontsize=9.5, color=INK)
    ax.set_ylabel('nRF difference (default \u2212 ensemble)', fontsize=9, color=INK)
    ax.set_title(f'{lab} (n={n}): d={d:+.2f}, Wilcoxon p={p:.3f}\n'
                 f'ensemble better in {wins}/{n} seeds', fontsize=9.5, color=INK, pad=6)
    ax.grid(alpha=0.25, axis='y')
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
panel_label(axes[1][0], 'C')
panel_label(axes[1][1], 'D')

fig.suptitle('Supplementary Figure S7. Multi-k distance ensemble: per-seed analysis',
             fontsize=13, fontweight='bold', color=INK, y=0.985)
fig.tight_layout(rect=[0, 0, 1, 0.955])
fig.savefig('Supplementary_Figure_S7.pdf'); fig.savefig('Supplementary_Figure_S7.png', dpi=300)
print('S7 done')
