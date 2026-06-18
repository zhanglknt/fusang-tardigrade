"""
Generate Figure 4: Scalability Curves

Shows execution time and RAM usage vs number of taxa (n).
Compares: NJ, FastME, DCM methods.

Uses data from Table 7 in the manuscript.
"""

import os
os.environ['MPLCONFIGDIR'] = 'D:/系统发育树项目/Fusang/Fusang-main/matplotlib_cache'

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

BASE_DIR = Path("D:/系统发育树项目/Fusang/Fusang-main")

# Data from Table 7 (manuscript)
# n | NJ time (s) | FastME time (s) | RAM (MB) | nRF vs true
data = {
    'n': [200, 500, 1000, 2000, 5000, 10000],
    'nj_time': [2.8, 18.9, 27.0, 184, None, None],  # NJ O(n³) - omitted at n=5000+
    'fastme_time': [None, None, None, 55, 399, 70],  # FastME O(n²)
    'ram_mb': [45, 78, 156, 312, 780, 609],
    'nrf': [0.078, 0.081, 0.083, 0.085, 0.088, None]
}

n = np.array(data['n'])
nj_time = np.array([2.8, 18.9, 27.0, 184, np.nan, np.nan])
fastme_time = np.array([np.nan, np.nan, np.nan, 55, 399, 70])
ram = np.array(data['ram_mb'])

# Create figure with dual y-axes
fig, ax1 = plt.subplots(figsize=(12, 6))
fig.patch.set_facecolor('white')

# Plot execution time (left y-axis)
ax1.set_xlabel('Number of Taxa (n)', fontsize=12, fontweight='bold')
ax1.set_ylabel('Execution Time (seconds)', fontsize=12, fontweight='bold', color='blue')
ax1.tick_params(axis='y', labelcolor='blue')

# NJ time (O(n³))
nj_mask = ~np.isnan(nj_time)
ax1.plot(n[nj_mask], nj_time[nj_mask], 'o-', color='blue', linewidth=2, 
          markersize=8, label='NJ (O(n³))')

# FastME time (O(n²))
fastme_mask = ~np.isnan(fastme_time)
ax1.plot(n[fastme_mask], fastme_time[fastme_mask], 's-', color='green', linewidth=2, 
          markersize=8, label='FastME (O(n²))')

ax1.set_yscale('log')
ax1.grid(True, alpha=0.3, linestyle='--')
ax1.set_xscale('log')

# Plot RAM usage (right y-axis)
ax2 = ax1.twinx()
ax2.set_ylabel('RAM Usage (MB)', fontsize=12, fontweight='bold', color='red')
ax2.tick_params(axis='y', labelcolor='red')

ax2.plot(n, ram, 'd-', color='red', linewidth=2, markersize=8, label='RAM')

# Formatting
plt.title('Figure 4. IMMI L0-1 Scalability\n(Execution Time and RAM Usage vs Number of Taxa)', 
         fontsize=14, fontweight='bold', pad=20)

# Combined legend
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left', fontsize=10)

ax1.set_ylim([1, 10000])
ax2.set_ylim([10, 2000])

plt.tight_layout()
output_path = BASE_DIR / "Figure4_scalability.png"
plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"Saved: {output_path}")

pdf_path = BASE_DIR / "Figure4_scalability.pdf"
plt.savefig(pdf_path, bbox_inches='tight', facecolor='white')
print(f"Saved: {pdf_path}")

plt.close()
print("\nFigure 4 generation complete!")
