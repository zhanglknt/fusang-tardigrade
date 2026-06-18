"""
Generate Supplementary Figure S8: Runtime Scaling by Taxa Number

Shows how runtime scales with number of taxa for different methods.
"""

import os
os.environ['MPLCONFIGDIR'] = 'D:/系统发育树项目/Fusang/Fusang-main/matplotlib_cache'

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

# Set font to Arial
try:
    font_path = font_manager.findfont('Arial')
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial']
except:
    pass

plt.rcParams['font.size'] = 8
plt.rcParams['axes.linewidth'] = 0.5

# Data from scalability tests
n_values = np.array([50, 100, 200, 500, 1000, 2000, 5000, 10000])

# Runtime in seconds (approximated from benchmarks)
fusang_nj = np.array([0.5, 1.8, 4.6, 27.0, 89.0, 184.0, None, None])  # NJ O(n³)
fusang_fastme = np.array([0.3, 0.9, 2.8, 12.4, 38.0, 55.0, 399.0, 70.0])  # FastME
fasttree2 = np.array([1.2, 2.5, 5.5, 18.9, 45.0, 120.0, None, None])  # ML
iqtree2 = np.array([2.5, 5.8, 67.2, 345.0, None, None, None, None])  # ML (slower)

# Create figure
fig, ax = plt.subplots(1, 1, figsize=(7, 5))

# Plot lines
ax.plot(n_values[:6], fusang_nj[:6], 'o-', color='#2E8B57', linewidth=1.5, 
         markersize=5, label='Fusang (NJ)')
ax.plot(n_values[:8], fusang_fastme, 's-', color='#3CB371', linewidth=1.5, 
         markersize=5, label='Fusang (FastME)')
ax.plot(n_values[:6], fasttree2[:6], '^-', color='#DC143C', linewidth=1.5, 
         markersize=5, label='FastTree2')
ax.plot(n_values[:5], iqtree2[:5], 'D-', color='#8B0000', linewidth=1.5, 
         markersize=5, label='IQ-TREE2')

ax.set_xlabel('Number of taxa (n)', fontsize=8)
ax.set_ylabel('Runtime (seconds)', fontsize=8)
ax.set_title('Supplementary Figure S8: Runtime Scaling by Taxa Number', fontsize=9, fontweight='bold')
ax.legend(loc='upper left', fontsize=7, frameon=True, edgecolor='black')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xscale('log')
ax.set_yscale('log')
ax.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)

plt.tight_layout()
plt.savefig('Supplementary_Figure_S8_runtime_scaling.png', dpi=300, bbox_inches='tight')
plt.savefig('Supplementary_Figure_S8_runtime_scaling.pdf', dpi=300, bbox_inches='tight')
print("[OK] Supplementary Figure S8 saved")
