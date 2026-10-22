"""Supplementary Figure S6 (v3): effect size analysis across benchmarks — 4 panels.

A: Cohen's d forest plot, Fusang vs FastTree2 across benchmarks
   (orientation: FT2 - Fusang; positive = Fusang better)
B: effect sizes vs alignment-free competitors and ablations
C: BAliBASE per-family Fusang-vs-FastTree2 nRF distribution (20 families)
D: bootstrap distribution of the 130-seed paired mean difference (10,000 resamples)

All values are computed by ``figure_data.py`` from the standard unrooted RF artefacts;
the earlier ``figures/effect_sizes.json`` input was written under the superseded
root-dependent metric (and claims like "65% of BAliBASE families below 0.5, median 0.45"
are withdrawn).
"""
import os
os.environ['MPLCONFIGDIR'] = 'matplotlib_cache'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import figure_data as FD

INK, GREY = '#1a2332', '#64748b'
TEAL, AMBER, RED = '#0d9488', '#d97706', '#dc2626'

fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.6))
fig.patch.set_facecolor('white')

def panel_label(ax, s):
    ax.text(-0.14, 1.06, s, transform=ax.transAxes, fontsize=14,
            fontweight='bold', color=INK, va='top')

# ---------- A: forest plot Fusang vs FT2 ----------
ax = axes[0][0]
bench = FD.bench_130seeds()
cs = FD.clean_scale_effects()
items = [
    ('indel n=200 (130 seeds)', bench['d'], bench['ci'][0], bench['ci'][1], TEAL),
    ('clean n=500 (30 seeds)', cs[500]['d'], cs[500]['ci'][0], cs[500]['ci'][1], RED),
    ('clean n=1000 (30 seeds)', cs[1000]['d'], cs[1000]['ci'][0], cs[1000]['ci'][1], RED),
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
lo_all = min(it[2] for it in items)
ax.set_xlim(lo_all - 2.0, 4.5)
ax.set_xlabel("Cohen's d (FT2 \u2212 Fusang; + = Fusang better)", fontsize=9, color=INK)
ax.text(lo_all - 1.8, 0.35, '\u2190 FT2 better', fontsize=8, color=RED)
ax.text(4.2, 0.35, 'Fusang better \u2192', fontsize=8, color=TEAL, ha='right')
ax.set_title('Fusang vs FastTree2 across benchmarks', fontsize=10.5, color=INK, pad=6)
ax.grid(alpha=0.25, axis='x')
for s in ('top', 'right'): ax.spines[s].set_visible(False)
panel_label(ax, 'A')

# ---------- B: effect sizes vs AF competitors ----------
ax = axes[0][1]
comp = FD.competitor_effect_sizes()
labs = [c[0] for c in comp]
ds = [c[1] for c in comp]
ps = [c[2] for c in comp]
ns = [c[3] for c in comp]
bars = ax.bar(range(len(ds)), ds, color=TEAL, width=0.58)
for b, d, p in zip(bars, ds, ps):
    ax.text(b.get_x() + b.get_width()/2, d + 0.25, f'{d:+.2f}', ha='center',
            fontsize=8.5, fontweight='bold', color=INK)
ax.set_xticks(range(len(ds)))
ax.set_xticklabels(labs, fontsize=7.0, rotation=22, ha='right')
ax.set_ylabel("paired Cohen's d (Fusang better)", fontsize=9, color=INK)
ax.set_ylim(0, max(ds) * 1.15)
ax.set_title('Effect sizes vs competitors / ablations\n(standard unrooted RF; integer-RF paired tests)',
             fontsize=9.6, color=INK, pad=6)
ax.grid(alpha=0.25, axis='y')
for s in ('top', 'right'): ax.spines[s].set_visible(False)
panel_label(ax, 'B')

# ---------- C: BAliBASE per-family ----------
ax = axes[1][0]
bb = FD.balibase_families()
nrf = bb['values']
rng = np.random.default_rng(7)
jit = rng.uniform(-0.16, 0.16, len(nrf))
ax.scatter(nrf, jit, s=42, color=TEAL, alpha=0.8, edgecolor=INK, lw=0.5, zorder=3)
ax.axvline(bb['median'], color=AMBER, ls='--', lw=1.6,
           label=f"median = {bb['median']:.3f}")
ax.axvline(0.5, color=GREY, ls=':', lw=1.2)
ax.text(0.52, 0.20, f"only {bb['n_below']}/{bb['n']} families < 0.5", fontsize=8.5,
        color=INK)
ax.text(0.52, 0.10, f"mean = {bb['mean']:.3f}", fontsize=8.5, color=INK)
ax.text(0.02, 0.97, 'the earlier "65% below 0.5, median 0.45" claim\n'
                    '(superseded root-dependent metric) is withdrawn',
        transform=ax.transAxes, fontsize=7.0, color=GREY, style='italic', va='top')
ax.set_yticks([])
ax.set_xlabel('nRF: Fusang vs FastTree2 tree (standard unrooted RF)', fontsize=9, color=INK)
ax.set_xlim(-0.08, 1.12)
ax.set_title(f'BAliBASE v3.0 protein families (n={bb["n"]})',
             fontsize=10.5, color=INK, pad=6)
ax.legend(fontsize=8.5, frameon=False, loc='upper right')
ax.grid(alpha=0.25, axis='x')
for s in ('top', 'right', 'left'): ax.spines[s].set_visible(False)
panel_label(ax, 'C')

# ---------- D: bootstrap distribution ----------
ax = axes[1][1]
rng = np.random.default_rng(7)
diffs = bench['ft2'] - bench['fus']   # FT2 − Fusang; positive = Fusang better
boot_mean = np.array([rng.choice(diffs, len(diffs), replace=True).mean()
                      for _ in range(10000)])
cim = np.percentile(boot_mean, [2.5, 97.5])
ax.hist(boot_mean, bins=40, color=TEAL, alpha=0.75, edgecolor='white')
ax.axvline(0, color=INK, lw=1.0)
ax.axvline(boot_mean.mean(), color=AMBER, ls='--', lw=1.6,
           label=f'mean = {boot_mean.mean():+.4f}')
ax.axvspan(cim[0], cim[1], color=AMBER, alpha=0.12)
ax.text(cim[1], ax.get_ylim()[1] * 0.02, f' 95% CI [{cim[0]:+.4f}, {cim[1]:+.4f}]',
        fontsize=8, color=INK, va='bottom')
ax.set_xlabel('Bootstrap mean nRF difference (FT2 \u2212 Fusang)', fontsize=9, color=INK)
ax.set_ylabel('Frequency', fontsize=9.5, color=INK)
ax.set_xticks(np.arange(-0.002, 0.016, 0.004))
ax.set_xticklabels([f'{x:+.3f}' for x in np.arange(-0.002, 0.016, 0.004)], fontsize=8)
ax.set_title('130-seed paired difference, 10,000 bootstrap resamples',
             fontsize=10, color=INK, pad=6)
ax.legend(fontsize=8.5, frameon=False)
ax.grid(alpha=0.25, axis='y')
for s in ('top', 'right'): ax.spines[s].set_visible(False)
panel_label(ax, 'D')

fig.suptitle('Supplementary Figure S6. Effect size analysis across benchmarks',
             fontsize=13, fontweight='bold', color=INK, y=0.985)
fig.tight_layout(rect=[0, 0, 1, 0.955])
fig.savefig('Supplementary_Figure_S6.pdf'); fig.savefig('Supplementary_Figure_S6.png', dpi=300)
print('S6 done | A:', ['%+.2f [%.2f, %.2f]' % (d, lo, hi) for _, d, lo, hi, _ in items])
print('       | B:', ['%s %+.2f' % (l.replace('\n', ' '), d) for l, d, _, _ in comp])
print('       | C: median %.3f mean %.3f below0.5 %d/%d'
      % (bb['median'], bb['mean'], bb['n_below'], bb['n']))
