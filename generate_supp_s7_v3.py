"""Supplementary Figure S7 (v3): multi-k ensemble per-seed analysis — 4 panels.

All panels read the 130-seed INDELible benchmark (n=200, indel=0.02, TRUE-relative,
standard root-independent unrooted RF) via ``figure_data.ensemble_130seeds()`` and
``figure_data.ensemble_pairwise()``; the earlier two seed-set inputs (27/30 seeds, FT2-
relative, superseded scorer) are no longer used.

A: per-seed scatter, NJ-only spaced k=5,gap2 vs the multi-k ensemble
B: configuration comparison (k=7, ensemble, full pipeline, NJ-only, k=9, k=5)
C: sorted per-seed differences, ensemble vs contiguous k=5
D: sorted per-seed differences, ensemble vs the full pipeline
"""
import os
os.environ['MPLCONFIGDIR'] = 'matplotlib_cache'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import figure_data as FD

INK, GREY = '#1a2332', '#64748b'
BLUE, TEAL, AMBER, RED = '#2563eb', '#0d9488', '#d97706', '#dc2626'
GREEN = '#16a34a'

e = FD.ensemble_130seeds()
pw = FD.ensemble_pairwise()
ens = e['multik_ensemble']['true']
njonly = e['fusang_spaced_k5_nj']['true']
full = e['fusang_spaced_k5']['true']
k5 = e['cosine_k5']['true']
n = len(ens)

fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.6))
fig.patch.set_facecolor('white')

def panel_label(ax, s):
    ax.text(-0.14, 1.06, s, transform=ax.transAxes, fontsize=14,
            fontweight='bold', color=INK, va='top')

def fmt_p(p):
    # three significant figures, matching the manuscript; small p-values keep a trailing
    # zero in the mantissa (9.1046e-8 -> 9.10e-8, not 9.1e-8)
    if p < 1e-3:
        ex = int(np.floor(np.log10(p)))
        return '%.2fe%d' % (p / 10 ** ex, ex)
    return '%.3g' % p

# ---------- A: scatter, NJ-only spaced vs ensemble ----------
ax = axes[0][0]
lo = min(njonly.min(), ens.min()) * 0.94
hi = max(njonly.max(), ens.max()) * 1.05
ax.plot([lo, hi], [lo, hi], '--', color=RED, lw=1.2, alpha=0.7, label='y = x (no change)')
ax.scatter(njonly, ens, s=20, color=TEAL, edgecolor=INK, lw=0.4, alpha=0.75, zorder=3)
v = pw['fusang_spaced_k5_nj']
ax.text(lo + 0.006, hi * 0.93,
        f'ensemble better: {v["ensemble_wins"]}/{n} seeds\n'
        f'p={fmt_p(v["p"])}, paired d={v["d"]:+.3f}', fontsize=8.5, color=INK,
        bbox=dict(boxstyle='round,pad=0.3', fc='#f4f7fa', ec=GREY, lw=0.6))
ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
ax.set_xlabel('Spaced k=5,gap2 NJ-only nRF', fontsize=9.5, color=INK)
ax.set_ylabel('Multi-k ensemble (k=5,7,9) nRF', fontsize=9, color=INK)
ax.set_title(f'Per-seed comparison (n={n} seeds)', fontsize=10.5, color=INK, pad=6)
ax.legend(fontsize=8, frameon=False, loc='lower right')
ax.grid(alpha=0.25)
for s in ('top', 'right'): ax.spines[s].set_visible(False)
panel_label(ax, 'A')

# ---------- B: configuration bars ----------
ax = axes[0][1]
order = FD.ensemble_vs_ensemble_names()
names = [e[c]['label'] for c in order]
means = [e[c]['mean'] for c in order]
sds = [e[c]['sd'] for c in order]
cols = [GREEN if c == 'multik_ensemble' else (TEAL if c == 'cosine_k7' else GREY)
        for c in order]
bars = ax.bar(range(len(order)), means, yerr=sds, capsize=4, color=cols, width=0.62,
              error_kw=dict(lw=1.1, ecolor=INK))
for i, c in enumerate(order):
    ax.text(i, means[i] + sds[i] + 0.008, f'{means[i]:.4f}', ha='center',
            fontsize=7.6, fontweight='bold', color=INK)
ax.set_xticks(range(len(order)))
ax.set_xticklabels(names, fontsize=7.2, rotation=24, ha='right')
ax.set_ylabel('nRF vs true tree (mean \u00b1 SD)', fontsize=9.5, color=INK)
ax.set_ylim(0, 0.27)
ax.set_title(f'Configuration comparison ({n} seeds, TRUE-relative)', fontsize=10.5,
             color=INK, pad=6)
ax.grid(alpha=0.25, axis='y')
for s in ('top', 'right'): ax.spines[s].set_visible(False)
panel_label(ax, 'B')

# ---------- C/D: sorted per-seed differences ----------
def diff_panel(ax, other_col):
    dd = e[other_col]['true'] - ens      # positive = ensemble better
    v = pw[other_col]
    order_d = np.argsort(dd)
    x = np.arange(1, len(dd) + 1)
    cols = [TEAL if val > 0 else RED for val in dd[order_d]]
    ax.bar(x, dd[order_d], width=1.0, color=cols, edgecolor='none')
    ax.axhline(0, color=INK, lw=0.8)
    ax.set_xlabel('Seed (sorted)', fontsize=9.5, color=INK)
    ax.set_ylabel('nRF difference (other − ensemble)', fontsize=9, color=INK)
    ax.set_title(f'ensemble vs {v["label"]}\n'
                 f'p={fmt_p(v["p"])}, d={v["d"]:+.3f}; ensemble better '
                 f'{v["ensemble_wins"]}/{n} (other {v["other_wins"]}, ties {v["ties"]})',
                 fontsize=8.8, color=INK, pad=6)
    ax.grid(alpha=0.25, axis='y')
    for s in ('top', 'right'): ax.spines[s].set_visible(False)

diff_panel(axes[1][0], 'cosine_k5')
panel_label(axes[1][0], 'C')
diff_panel(axes[1][1], 'fusang_spaced_k5')
axes[1][1].text(0.03, 0.80,
                'under the standard RF the ensemble now *does*\nexceed the default pipeline '
                '(earlier ordering reversed)', transform=axes[1][1].transAxes,
                fontsize=7.2, color=GREY, style='italic', va='top')
panel_label(axes[1][1], 'D')

fig.suptitle('Supplementary Figure S7. Multi-k distance ensemble: per-seed analysis (130 seeds)',
             fontsize=13, fontweight='bold', color=INK, y=0.985)
fig.tight_layout(rect=[0, 0, 1, 0.955])
fig.savefig('Supplementary_Figure_S7.pdf'); fig.savefig('Supplementary_Figure_S7.png', dpi=300)
print('S7 done |', ' '.join('%s %.4f' % (c, e[c]['mean']) for c in order))
for c in ('cosine_k5', 'fusang_spaced_k5'):
    v = pw[c]
    print('         ensemble vs %-20s wins %d/%d ties %d p=%s d=%+.3f'
          % (c, v['ensemble_wins'], n, v['ties'], fmt_p(v['p']), v['d']))
