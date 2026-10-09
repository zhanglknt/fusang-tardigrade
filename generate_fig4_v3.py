"""Figure 4 (v3): Spaced k-mer parameter optimization — 4 panels.
A: k x gap heatmap at n=200 (S1_grid_search_results.csv)
B: adaptive parameter selection logic (schematic)
C: 5-repeat stability validation summary (Supplementary Table S5)
D: multi-k ensemble vs default across two independent seed sets
"""
import os
os.environ['MPLCONFIGDIR'] = 'matplotlib_cache'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import csv, json
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from scipy.stats import wilcoxon

INK, GREY = '#1a2332', '#64748b'
BLUE, TEAL, AMBER, RED, GREEN = '#2563eb', '#0d9488', '#d97706', '#dc2626', '#16a34a'

fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.6))
fig.patch.set_facecolor('white')

def panel_label(ax, s):
    ax.text(-0.14, 1.06, s, transform=ax.transAxes, fontsize=14,
            fontweight='bold', color=INK, va='top')

# ---------- A: heatmap k x gap at n=200 ----------
ax = axes[0][0]
rows = [r for r in csv.DictReader(open('data/S1_grid_search_results.csv'))
        if int(r['n_taxa']) == 200]
ks = sorted(set(int(r['k']) for r in rows))
gs = ['none', 'gap1', 'gap2', 'gap3']
M = np.full((len(gs), len(ks)), np.nan)
for r in rows:
    M[gs.index(r['gap_pattern']), ks.index(int(r['k']))] = float(r['nRF'])
im = ax.imshow(M, cmap='viridis_r', aspect='auto')
for i in range(len(gs)):
    for j in range(len(ks)):
        if not np.isnan(M[i, j]):
            ax.text(j, i, f'{M[i,j]:.3f}', ha='center', va='center', fontsize=8,
                    color='white' if M[i, j] > np.nanmedian(M) else INK)
ax.set_xticks(range(len(ks))); ax.set_xticklabels([f'k={k}' for k in ks], fontsize=9)
ax.set_yticks(range(len(gs)))
ax.set_yticklabels(['gap0\n(contig.)', 'gap1', 'gap2', 'gap3'], fontsize=9)
fig.colorbar(im, ax=ax, label='nRF (lower = better)', shrink=0.85)
# mark default
if 5 in ks:
    j = ks.index(5); i = gs.index('gap2')
    ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec=RED, lw=2.4))
ax.set_title('Grid search at n=200 (red box = default k=5,gap2)',
             fontsize=10, color=INK, pad=6)
panel_label(ax, 'A')

# ---------- B: adaptive logic schematic ----------
ax = axes[0][1]
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis('off')
ax.set_title('Adaptive k,gap selection by dataset size', fontsize=10.5, color=INK, pad=6)
panel_label(ax, 'B')
ax.add_patch(FancyBboxPatch((0.5, 5.6), 4.0, 2.6, boxstyle='round,pad=0.1,rounding_size=0.25',
                            fc=TEAL, ec='none'))
ax.text(2.5, 7.65, 'n \u2264 100', fontsize=11, color='white', ha='center',
        va='center', fontweight='bold')
ax.text(2.5, 6.35, 'k=4, gap1\nsmaller pattern for\nlimited taxon sampling', fontsize=9,
        color='white', ha='center', va='center')
ax.add_patch(FancyBboxPatch((5.5, 5.6), 4.0, 2.6, boxstyle='round,pad=0.1,rounding_size=0.25',
                            fc=BLUE, ec='none'))
ax.text(7.5, 7.65, 'n > 100', fontsize=11, color='white', ha='center',
        va='center', fontweight='bold')
ax.text(7.5, 6.35, 'k=5, gap2\nlonger pattern for\nricher k-mer statistics', fontsize=9,
        color='white', ha='center', va='center')
ax.add_patch(FancyBboxPatch((2.2, 1.6), 5.6, 2.6, boxstyle='round,pad=0.1,rounding_size=0.25',
                            fc='#f4f7fa', ec=AMBER, lw=2))
ax.text(5.0, 3.55, 'genome scale (>>1 kb): multi-k ensemble', fontsize=9.5, color=AMBER,
        ha='center', va='center', fontweight='bold')
ax.text(5.0, 2.45, 'k=5 saturates (1,024 possible 5-mers);\nensemble of k=5,7,9 rescues signal — no manual k choice',
        fontsize=8.5, color=INK, ha='center', va='center')

# ---------- C: stability summary ----------
ax = axes[1][0]
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis('off')
ax.set_title('5-repeat stability validation (Suppl. Table S5)', fontsize=10.5, color=INK, pad=6)
panel_label(ax, 'C')
scales = [20, 50, 100, 200, 500]
for i, n in enumerate(scales):
    y = 8.2 - i * 1.45
    ax.text(1.0, y, f'n={n}', fontsize=10, color=INK, fontweight='bold', va='center')
    ax.text(4.4, y, '5/5 repeats: identical k,gap selection', fontsize=9, color=INK, va='center')
    ax.text(9.3, y, '\u2713', fontsize=13, color=GREEN, fontweight='bold', ha='center', va='center')
ax.text(5.0, 0.6, 'nRF deviation across repeats < 0.005 at every scale (100% reproducible)',
        fontsize=9, color=GREY, ha='center', style='italic')

# ---------- D: multi-k across seed sets ----------
ax = axes[1][1]
d7 = json.load(open('table7_recomputed_recovered.json'))
resA = d7['results'] if isinstance(d7, dict) and 'results' in d7 else d7
if isinstance(resA, dict):
    seedsA = list(resA.keys())
    def _getA(entry, key):
        e = resA[entry]
        return e.get(key) or e.get(key.replace('nrf_', ''))
    # fall back: inspect structure
    first = resA[seedsA[0]]
    keysA = list(first.keys())
    fk = [k for k in keysA if 'fusang' in k.lower() or 'default' in k.lower() or 'spaced' in k.lower()][0]
    mk = [k for k in keysA if 'multi' in k.lower() or 'ensemble' in k.lower()][0]
    defA = np.array([resA[s][fk] for s in seedsA], dtype=float)
    mulA = np.array([resA[s][mk] for s in seedsA], dtype=float)
else:
    defA = np.array([r['fusang'] for r in resA], dtype=float)
    mulA = np.array([r['multik'] for r in resA], dtype=float)
rowsB = list(csv.DictReader(open('benchmark_multik_ensemble_n200_indel.csv')))
defB = np.array([float(r['nrf_fusang_original']) for r in rowsB])
mulB = np.array([float(r['nrf_multik_ensemble']) for r in rowsB])

def paired_d(a, b):  # d of (default - multik): positive = ensemble better
    dd = a - b
    return dd.mean() / dd.std(ddof=1)

pA = wilcoxon(defA, mulA).pvalue
pB = wilcoxon(defB, mulB).pvalue
x = np.arange(2)
w = 0.34
b1 = ax.bar(x - w/2, [defA.mean(), defB.mean()], w,
            yerr=[defA.std(ddof=1), defB.std(ddof=1)], capsize=4,
            color=GREY, label='default spaced k=5,gap2', error_kw=dict(lw=1.1, ecolor=INK))
b2 = ax.bar(x + w/2, [mulA.mean(), mulB.mean()], w,
            yerr=[mulA.std(ddof=1), mulB.std(ddof=1)], capsize=4,
            color=TEAL, label='multi-k ensemble (k=5,7,9)', error_kw=dict(lw=1.1, ecolor=INK))
for bars in (b1, b2):
    for b in bars:
        ax.text(b.get_x() + b.get_width()/2, b.get_height() + 0.012,
                f'{b.get_height():.3f}', ha='center', fontsize=8.5, color=INK)
ax.set_xticks(x)
ax.set_xticklabels([f'seed set A (n=27)\nd={paired_d(defA, mulA):+.2f}, p={pA:.2f}',
                    f'seed set B (n=30)\nd={paired_d(defB, mulB):+.2f}, p={pB:.3f}'],
                   fontsize=8.5)
ax.set_ylabel('nRF vs true tree (mean ± SD)', fontsize=9.5, color=INK)
ax.set_ylim(0, 0.16)
ax.legend(fontsize=8, frameon=False, loc='upper right')
ax.set_title('Ensemble effect is seed-set-dependent\n(pooled n=57: d=\u22120.17, p=0.19 — no general advantage at L=500 bp)',
             fontsize=9, color=INK, pad=6)
ax.spines[['top', 'right']].set_visible(False)
panel_label(ax, 'D')

fig.suptitle('Figure 4. Spaced k-mer parameter optimization', fontsize=14,
             fontweight='bold', color=INK, y=0.985)
fig.tight_layout(rect=[0, 0, 1, 0.955])
fig.savefig('Figure4.pdf')
fig.savefig('Figure4.png', dpi=300)
print('Figure4 done | A: %d/%d seeds d=%.3f p=%.3f | B: %d/%d seeds d=%.3f p=%.4f'
      % (len(defA), len(defA), paired_d(defA, mulA), pA,
         len(defB), len(defB), paired_d(defB, mulB), pB))
