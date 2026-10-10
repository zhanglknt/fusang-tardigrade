"""Supplementary Figure S6 (v3): effect size analysis — 4 panels.
A: Cohen's d forest plot, Fusang vs FastTree2 across benchmarks
   (orientation: FT2 - Fusang; positive = Fusang better)
B: effect sizes vs alignment-free competitors (from manuscript-reported d)
C: BAliBASE per-family nRF distribution (FT2-relative, 20 families)
D: bootstrap distribution of the 112-seed paired mean difference
"""
import os
os.environ['MPLCONFIGDIR'] = 'matplotlib_cache'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import csv, json

INK, GREY = '#1a2332', '#64748b'
TEAL, AMBER, RED = '#0d9488', '#d97706', '#dc2626'

fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.6))
fig.patch.set_facecolor('white')

def panel_label(ax, s):
    ax.text(-0.14, 1.06, s, transform=ax.transAxes, fontsize=14,
            fontweight='bold', color=INK, va='top')

# ---------- A: forest plot Fusang vs FT2 ----------
ax = axes[0][0]
rows = list(csv.DictReader(open('table3_corrected.csv')))
fus = np.array([float(r['fusang_nrf']) for r in rows])
ft2 = np.array([float(r['ft2_nrf']) for r in rows])
keep = ~((fus > 0.3) | (ft2 > 0.3))
diffs = ft2[keep] - fus[keep]
rng = np.random.default_rng(7)

es = json.load(open('figures/effect_sizes.json'))
ov = es['overall_130seed']
n500 = es['scale_clean']['n500']; n1000 = es['scale_clean']['n1000']
items = [
    (f'indel n=200 (130 seeds, all valid)', ov['cohens_d'], ov['ci_95_low'], ov['ci_95_high'], TEAL),
    (f'clean n=500 (30 seeds)', n500['cohens_d'], n500['ci_95_low'], n500['ci_95_high'], RED),
    (f'clean n=1000 (30 seeds)', n1000['cohens_d'], n1000['ci_95_low'], n1000['ci_95_high'], RED),
]
for i, (lab, d, lo, hi, c) in enumerate(items):
    y = len(items) - i
    ax.plot([lo, hi], [y, y], color=c, lw=2.4, solid_capstyle='round')
    ax.plot(d, y, 'o', color=c, ms=9, mec=INK, mew=0.8)
    ax.text(hi + 0.25, y, f'd={d:+.2f} [{lo:+.2f}, {hi:+.2f}]', fontsize=7.8,
            color=INK, va='center')
ax.axvline(0, color=INK, lw=1.0)
ax.set_yticks([len(items) - i for i in range(len(items))])
ax.set_yticklabels([it[0] for it in items], fontsize=9)
ax.set_xlim(-8, 4.5)
ax.set_xlabel("Cohen's d (FT2 \u2212 Fusang; + = Fusang better)", fontsize=9, color=INK)
ax.text(-7.6, 0.35, '\u2190 FT2 better', fontsize=8, color=RED)
ax.text(4.2, 0.35, 'Fusang better \u2192', fontsize=8, color=TEAL, ha='right')
ax.set_title('Fusang vs FastTree2 across benchmarks', fontsize=10.5, color=INK, pad=6)
ax.grid(alpha=0.25, axis='x')
for s in ('top', 'right'): ax.spines[s].set_visible(False)
panel_label(ax, 'A')

# ---------- B: effect sizes vs AF competitors ----------
ax = axes[0][1]
comp = [
    ('Co-phylog k=19\n(DNA indel)', 11.91),
    ('kmacs k=3\n(DNA indel)', 2.54),
    ('IQ-TREE2 GTR\n(DNA indel)', 3.10),
    ('SwissTree\n(protein)', 1.32),
    ('L1 vs L0\n(coalescent)', 3.82),
]
labs = [c[0] for c in comp]; ds = [c[1] for c in comp]
bars = ax.bar(range(len(ds)), ds, color=TEAL, width=0.58)
for b, d in zip(bars, ds):
    ax.text(b.get_x() + b.get_width()/2, d + 0.25, f'{d:.2f}', ha='center',
            fontsize=8.5, fontweight='bold', color=INK)
ax.set_xticks(range(len(ds))); ax.set_xticklabels(labs, fontsize=7.2)
ax.set_ylabel("paired Cohen's d (Fusang better)", fontsize=9, color=INK)
ax.set_ylim(0, 13.5)
ax.set_title('Effect sizes vs competitors / ablations', fontsize=10.5, color=INK, pad=6)
ax.grid(alpha=0.25, axis='y')
for s in ('top', 'right'): ax.spines[s].set_visible(False)
panel_label(ax, 'B')

# ---------- C: BAliBASE per-family ----------
ax = axes[1][0]
bb = json.load(open('balibase/bench_results/balibase_results.json'))
nrf = np.array([f['nRF_vs_FT2'] for f in bb])
nrf_c = np.clip(nrf, 0, 2.0)
jit = rng.uniform(-0.16, 0.16, len(nrf_c))
ax.scatter(nrf_c, jit, s=42, color=TEAL, alpha=0.8, edgecolor=INK, lw=0.5, zorder=3)
ax.axvline(np.median(nrf_c), color=AMBER, ls='--', lw=1.6,
           label=f'median = {np.median(nrf_c):.3f}')
ax.axvline(0.5, color=GREY, ls=':', lw=1.2)
good = (nrf_c < 0.5).mean() * 100
ax.text(0.52, 0.20, f'{good:.0f}% of families < 0.5', fontsize=8.5, color=INK)
ax.set_yticks([])
ax.set_xlabel('nRF vs FastTree2 tree (FT2-relative, clipped at 2.0)', fontsize=9, color=INK)
ax.set_xlim(-0.08, 2.1)
ax.set_title(f'BAliBASE v3.0 protein families (n={len(nrf_c)})', fontsize=10.5, color=INK, pad=6)
ax.legend(fontsize=8.5, frameon=False, loc='upper right')
ax.grid(alpha=0.25, axis='x')
for s in ('top', 'right', 'left'): ax.spines[s].set_visible(False)
panel_label(ax, 'C')

# ---------- D: bootstrap distribution ----------
ax = axes[1][1]
boot_mean = np.array([rng.choice(diffs, len(diffs), replace=True).mean()
                      for _ in range(10000)])
cim = np.percentile(boot_mean, [2.5, 97.5])
ax.hist(boot_mean, bins=40, color=TEAL, alpha=0.75, edgecolor='white')
ax.axvline(0, color=INK, lw=1.0)
ax.axvline(boot_mean.mean(), color=AMBER, ls='--', lw=1.6,
           label=f'mean = {boot_mean.mean():+.4f}')
ax.axvspan(cim[0], cim[1], color=AMBER, alpha=0.12)
ax.text(cim[0], ax.get_ylim()[1] * 0.02, f'95% CI [{cim[0]:+.4f}, {cim[1]:+.4f}]',
        fontsize=8, color=INK, va='bottom')
ax.set_xlabel('Bootstrap mean nRF difference (FT2 \u2212 Fusang)', fontsize=9, color=INK)
ax.set_ylabel('Frequency', fontsize=9.5, color=INK)
ax.set_title('112-seed paired difference, 10,000 bootstrap resamples',
             fontsize=10, color=INK, pad=6)
ax.legend(fontsize=8.5, frameon=False)
ax.grid(alpha=0.25, axis='y')
for s in ('top', 'right'): ax.spines[s].set_visible(False)
panel_label(ax, 'D')

fig.suptitle('Supplementary Figure S6. Effect size analysis across benchmarks',
             fontsize=13, fontweight='bold', color=INK, y=0.985)
fig.tight_layout(rect=[0, 0, 1, 0.955])
fig.savefig('Supplementary_Figure_S6.pdf'); fig.savefig('Supplementary_Figure_S6.png', dpi=300)
print('S6 done')
