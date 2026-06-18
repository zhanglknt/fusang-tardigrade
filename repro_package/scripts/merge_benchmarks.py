"""Merge NJ, MHL, FT2, and old Fusang benchmark results."""
import csv, statistics as st
from scipy import stats

# Load new NJ+MHL data
new_data = {}
with open('benchmark_mhl_n200_30seeds.csv') as f:
    for r in csv.DictReader(f):
        s = int(r['seed'])
        n = int(r['n'])
        nj = float(r['nj_nrf']) if r['nj_nrf'] != 'nan' else None
        mhl = float(r['mhl_nrf']) if r['mhl_nrf'] != 'nan' else None
        new_data[s] = {'n': n, 'nj': nj, 'mhl': mhl}

# Load old FT2+Fusang data
old_data = {}
with open('indel_benchmark_seeds100_129.csv') as f:
    for r in csv.DictReader(f):
        s = int(r['seed'])
        old_data[s] = {
            'ft2': float(r['ft2_nrf']),
            'fusang_old': float(r['fusang_nrf']),
            'fusang_better': r['fusang_better'] == 'yes',
        }

# Merge: n=200 only, all non-NaN
nj_vals = []
mhl_vals = []
ft2_vals = []
fusold_vals = []

for s in range(100, 130):
    n = new_data.get(s, {})
    o = old_data.get(s, {})
    if n.get('n') != 200:
        continue
    if n.get('nj') is None or n.get('mhl') is None:
        continue  # bad ref
    nj_vals.append(n['nj'])
    mhl_vals.append(n['mhl'])
    ft2_vals.append(o['ft2'])
    fusold_vals.append(o['fusang_old'])

N = len(nj_vals)
print(f"\n=== n=200 valid seeds: {N} ===")
print(f"{'Method':18s}  {'nRF mean':>8s}  {'nRF std':>8s}")
print("-" * 40)

methods = [
    ("NJ (baseline)", nj_vals),
    ("Fusang old (NJ)", fusold_vals),
    ("MHL (new)", mhl_vals),
    ("FT2 (MAFFT+ML)", ft2_vals),
]
for name, vals in methods:
    print(f"{name:18s}  {st.mean(vals):8.4f}  {st.stdev(vals):8.4f}")

print("\nPaired t-tests vs NJ baseline:")
for name, vals in methods:
    if name == "NJ (baseline)":
        continue
    t, p = stats.ttest_rel(vals, nj_vals)
    diff = st.mean(vals) - st.mean(nj_vals)
    d = diff / max(st.stdev(nj_vals), st.stdev(vals), 1e-6)
    better = sum(1 for i in range(N) if vals[i] < nj_vals[i])
    print(f"  {name:18s}: delta={diff:+.4f}, t={t:+.2f}, p={p:.4f}, d={d:+.2f}, better={better}/{N}")

print("\nKey pairwise comparisons:")
pairs = [
    ("FT2", ft2_vals, "MHL", mhl_vals),
    ("FT2", ft2_vals, "FusangOld", fusold_vals),
    ("MHL", mhl_vals, "FusangOld", fusold_vals),
]
for a_name, a_vals, b_name, b_vals in pairs:
    t, p = stats.ttest_rel(a_vals, b_vals)
    diff = st.mean(a_vals) - st.mean(b_vals)
    d = diff / max(st.stdev(a_vals), st.stdev(b_vals), 1e-6)
    print(f"  {a_name} vs {b_name}: delta={diff:+.4f}, t={t:+.2f}, p={p:.4f}, d={d:+.2f}")

# Save merged CSV
with open('benchmark_all_methods_n200.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['seed', 'nj_nrf', 'mhl_nrf', 'ft2_nrf', 'fusang_old_nrf'])
    for s in range(100, 130):
        n = new_data.get(s, {})
        o = old_data.get(s, {})
        nj = n.get('nj')
        mhl = n.get('mhl')
        ft2 = o.get('ft2', '')
        fuo = o.get('fusang_old', '')
        if nj is not None:
            w.writerow([s, nj, mhl if mhl is not None else '', ft2, fuo])

print("\nMerged CSV saved: benchmark_all_methods_n200.csv")
