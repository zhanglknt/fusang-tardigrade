"""
Generate Supplementary Figure S9: Memory Usage Scaling

Shows how memory usage scales with number of taxa for different methods.
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

# Memory in MB (approximated from benchmarks)
fusang_nj = np.array([12, 18, 45, 156, 312, 780, None, None])  # NJ
fusang_fastme = np.array([15, 22, 45, 134, 267, 445, 612, 609])  # FastME
fasttree2 = np.array([34, 45, 156, 890, 2100, None, None, None])  # ML (needs alignment)
iqtree2 = np.array([45, 78, 312, 1800, 4500, None, None, None])  # ML (slower)

# Create figure
fig, ax = plt.subplots(1, 1, figsize=(7, 5))

# Plot lines
ax.plot(n_values[:8], fusang_nj[:8], 'o-', color='#2E8B57', linewidth=1.5, 
         markersize=5, label='Fusang (NJ)')
ax.plot(n_values[:8], fusang_fastme, 's-', color='#3CB371', linewidth=1.5, 
         markersize=5, label='Fusang (FastME)')
ax.plot(n_values[:6], fasttree2[:6], '^-', color='#DC143C', linewidth=1.5, 
         markersize=5, label='FastTree2')
ax.plot(n_values[:5], iqtree2[:5], 'D-', color='#8B0000', linewidth=1.5, 
         markersize=5, label='IQ-TREE2')

ax.set_xlabel('Number of taxa (n)', fontsize=8)
ax.set_ylabel('Memory usage (MB)', fontsize=8)
ax.set_title('Supplementary Figure S9: Memory Usage Scaling', fontsize=9, fontweight='bold')
ax.legend(loc='upper left', fontsize=7, frameon=True, edgecolor='black')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xscale('log')
ax.set_yscale('log')
ax.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)

# Add annotations for n=10000
ax.annotate('n=10,000\n609 MB\n70 seconds', 
            xy=(10000, 609), xytext=(8000, 800),
            arrowprops=dict(arrowstyle='->', color='black', lw=0.8),
            fontsize=6, ha='center')

plt.tight_layout()
plt.savefig('Supplementary_Figure_S9_memory_scaling.png', dpi=300, bbox_inches='tight')
plt.savefig('Supplementary_Figure_S9_memory_scaling.pdf', dpi=300, bbox_inches='tight')
print("[OK] Supplementary Figure S9 saved")
