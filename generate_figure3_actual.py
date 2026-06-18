"""
Generate Figure 3: Multi-k Ensemble Results (Actual Data)

Uses real data from benchmark_multik_ensemble_n200_indel.csv
"""

import os
os.environ['MPLCONFIGDIR'] = 'D:/系统发育树项目/Fusang/Fusang-main/matplotlib_cache'

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

BASE_DIR = Path("D:/系统发育树项目/Fusang/Fusang-main")

# Load actual data
df = pd.read_csv(BASE_DIR / "benchmark_multik_ensemble_n200_indel.csv")
print(f"Loaded {len(df)} seeds from multi-k ensemble benchmark")

# Extract data for each method
k5_data = df['nrf_k5_contig'].dropna().tolist()
k7_data = df['nrf_k7_contig'].dropna().tolist()
k9_data = df['nrf_k9_contig'].dropna().tolist()
ensemble_data = df['nrf_multik_ensemble'].dropna().tolist()
original_data = df['nrf_fusang_original'].dropna().tolist()

print(f"\nData summary:")
print(f"  k=5 contiguous: {len(k5_data)} seeds, mean={np.mean(k5_data):.4f}")
print(f"  k=7 contiguous: {len(k7_data)} seeds, mean={np.mean(k7_data):.4f}")
print(f"  k=9 contiguous: {len(k9_data)} seeds, mean={np.mean(k9_data):.4f}")
print(f"  ensemble (avg): {len(ensemble_data)} seeds, mean={np.mean(ensemble_data):.4f}")
print(f"  original (k=5,gap2): {len(original_data)} seeds, mean={np.mean(original_data):.4f}")

# Create figure
fig, ax = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor('white')

# Data for boxplot
data = [k5_data, k7_data, k9_data, ensemble_data]
labels = ['k=5\ncontiguous', 'k=7\ncontiguous', 'k=9\ncontiguous', 'avg(k=5,7,9)\nensemble']

# Colors
colors = ['#4C72B0', '#55A868', '#C44E52', '#8172B3']
box_colors = ['lightblue', 'lightgreen', 'lightcoral', 'lavender']

# Create boxplot
bp = ax.boxplot(data, tick_labels=labels, patch_artist=True, showmeans=True)

# Color boxes
for patch, color in zip(bp['boxes'], box_colors):
    patch.set_facecolor(color)
    patch.set_edgecolor('black')
    patch.set_linewidth(1.5)

# Style medians
for median in bp['medians']:
    median.set_color('red')
    median.set_linewidth(2)

# Add significance marker for ensemble vs k=5 (p=0.006)
# Find the y position for the marker
ensemble_mean = np.mean(ensemble_data)
k5_mean = np.mean(k5_data)
y_max = max(ensemble_mean, k5_mean) + 0.02
ax.annotate('* p=0.006', xy=(1, k5_mean), xytext=(4, y_max + 0.01),
            arrowprops=dict(arrowstyle='-', lw=1.5, color='red'),
            fontsize=11, fontweight='bold', color='red', ha='center')

# Formatting
ax.set_ylabel('Normalized Robinson-Foulds distance (nRF)', fontsize=12, fontweight='bold')
ax.set_title('Figure 3. Multi-k Ensemble Improves Accuracy\n(n=200 Indel-Rich Data, 20 seeds)', 
             fontsize=14, fontweight='bold', pad=20)
ax.set_ylim([-0.02, 0.18])
ax.grid(True, axis='y', alpha=0.3, linestyle='--')

# Add n values
for i, d in enumerate(data):
    n = len(d)
    ax.text(i+1, -0.01, f'n={n}', ha='center', fontsize=9, style='italic')

# Add mean±SD as text
for i, d in enumerate(data):
    mean = np.mean(d)
    sd = np.std(d, ddof=1)
    ax.text(i+1, 0.15, f'{mean:.3f}±{sd:.3f}', ha='center', 
            fontsize=9, bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray'))

plt.tight_layout()
output_path = BASE_DIR / "Figure3_multik_ensemble_actual.png"
plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"\nSaved: {output_path}")

pdf_path = BASE_DIR / "Figure3_multik_ensemble_actual.pdf"
plt.savefig(pdf_path, bbox_inches='tight', facecolor='white')
print(f"Saved: {pdf_path}")

plt.close()
print("\nFigure 3 generation complete!")
