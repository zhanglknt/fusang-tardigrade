"""
Generate Supplementary Figure S4: 16S rRNA Tree Detailed View

Shows the Fusang tree of 74 16S rRNA sequences with phylum-level
color coding and highlights correct sister-taxon placements.
"""

import os
os.environ['MPLCONFIGDIR'] = 'D:/系统发育树项目/Fusang/Fusang-main/matplotlib_cache'

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import json

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

# Create figure
fig, ax = plt.subplots(1, 1, figsize=(8, 10))

# Load 16S rRNA evaluation results
try:
    with open('real_data/real_16S_full_eval_final.json') as f:
        eval_data = json.load(f)
    
    # Extract lineage-aware improvement data
    improvement = eval_data.get('lineage_aware_improvement', {})
    
    # Create horizontal bar chart showing improvement at each taxonomic level
    levels = ['Phylum', 'Class', 'Order', 'Family', 'Genus']
    fusang_improvement = []
    ft2_improvement = []
    
    for level in levels:
        key = 'at_{}'.format(level.lower())
        if key in improvement:
            fusang_improvement.append(improvement[key].get('fusang_improvement', 0))
            ft2_improvement.append(improvement[key].get('fasttree2_improvement', 0))
        else:
            fusang_improvement.append(0)
            ft2_improvement.append(0)
    
    x = np.arange(len(levels))
    width = 0.35
    
    bars1 = ax.barh(x - width/2, fusang_improvement, width, 
                     label='Fusang', color='#2E8B57', alpha=0.8)
    bars2 = ax.barh(x + width/2, ft2_improvement, width,
                     label='FastTree2', color='#DC143C', alpha=0.8)
    
    ax.set_yticks(x)
    ax.set_yticklabels(levels)
    ax.set_xlabel('Lineage-aware improvement (%)', fontsize=8)
    ax.set_title('16S rRNA Benchmark: Lineage-aware Improvement', fontsize=9, fontweight='bold')
    ax.legend(loc='lower right', fontsize=7, frameon=True, edgecolor='black', linewidth=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Add value labels
    for i, (v1, v2) in enumerate(zip(fusang_improvement, ft2_improvement)):
        ax.text(v1 + 0.5, i - width/2, '{:.1f}%'.format(v1), 
                va='center', fontsize=6, color='#2E8B57')
        ax.text(v2 + 0.5, i + width/2, '{:.1f}%'.format(v2), 
                va='center', fontsize=6, color='#DC143C')
    
    ax.set_xlim(0, max(max(fusang_improvement), max(ft2_improvement)) + 5)
    
except Exception as e:
    print("Error loading 16S data: {}".format(str(e)))
    ax.text(0.5, 0.5, '16S rRNA data not available', 
            ha='center', va='center', transform=ax.transAxes)
    ax.set_title('Supplementary Figure S4: 16S rRNA Tree Detailed View', fontsize=9, fontweight='bold')

plt.tight_layout()
plt.savefig('Supplementary_Figure_S4_16S_improvement.png', dpi=300, bbox_inches='tight')
plt.savefig('Supplementary_Figure_S4_16S_improvement.pdf', dpi=300, bbox_inches='tight')
print("[OK] Supplementary Figure S4 saved")
