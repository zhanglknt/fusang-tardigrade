"""
Generate Supplementary Figure S5: Computational Cost Breakdown

Shows the computational cost (time and memory) of each IMMI level
and comparison methods.
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
plt.rcParams['xtick.major.width'] = 0.5
plt.rcParams['ytick.major.width'] = 0.5

# Data from benchmarks
methods = ['Fusang\nL0-1', 'Fusang\nmulti-k', 'FastTree2', 'IQ-TREE2\nGTR', 'MAFFT\n(alignment)']
time_seconds = [4.6, 8.3, 5.5, 67.2, 12.4]  # n=200
memory_mb = [45, 78, 156, 312, 234]  # n=200

# Create figure with two subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

# Time comparison
colors = ['#2E8B57', '#3CB371', '#DC143C', '#8B0000', '#4682B4']
bars1 = ax1.barh(methods, time_seconds, color=colors, alpha=0.8)
ax1.set_xlabel('Time (seconds)', fontsize=8)
ax1.set_title('(A) Computational Time (n=200)', fontsize=9, fontweight='bold')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

# Add value labels
for i, (method, time) in enumerate(zip(methods, time_seconds)):
    ax1.text(time + 1, i, '{:.1f}s'.format(time), va='center', fontsize=6)

# Memory comparison
bars2 = ax2.barh(methods, memory_mb, color=colors, alpha=0.8)
ax2.set_xlabel('Memory (MB)', fontsize=8)
ax2.set_title('(B) Memory Usage (n=200)', fontsize=9, fontweight='bold')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

# Add value labels
for i, (method, mem) in enumerate(zip(methods, memory_mb)):
    ax2.text(mem + 10, i, '{} MB'.format(mem), va='center', fontsize=6)

plt.suptitle('Supplementary Figure S5: Computational Cost Breakdown', fontsize=10, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig('Supplementary_Figure_S5_cost_breakdown.png', dpi=300, bbox_inches='tight')
plt.savefig('Supplementary_Figure_S5_cost_breakdown.pdf', dpi=300, bbox_inches='tight')
print("[OK] Supplementary Figure S5 saved")
