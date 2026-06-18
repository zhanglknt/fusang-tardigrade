"""
Generate Figure 1: IMMI Framework Architecture Diagram

Schematic showing the 4-level architecture:
L0: k-mer feature extraction (O(nL))
L1: cosine distance + NJ backbone (O(n²))
L2: random forest boundary classifier
L3: MAFFT+ML subtree refinement (O(m²L²) per cluster)

Shows information flow and computational cost gradient.
"""

import os
os.environ['MPLCONFIGDIR'] = 'D:/系统发育树项目/Fusang/Fusang-main/matplotlib_cache'

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

BASE_DIR = Path("D:/系统发育树项目/Fusang/Fusang-main")

# Create figure
fig, ax = plt.subplots(figsize=(14, 10))
fig.patch.set_facecolor('white')
ax.set_facecolor('white')

# Define box positions and sizes
box_width = 2.5
box_height = 1.5
x_start = 1
y_levels = [8, 6, 4, 2]  # L0-L3
level_names = ['Level 0\nFeature Extraction', 
              'Level 1\nGlobal Distance Inference',
              'Level 2\nInformation-Aware Partitioning',
              'Level 3\nHigh-Resolution Refinement']

# Colors for each level
colors = ['#4C72B0', '#55A868', '#C44E52', '#8172B3']
info_labels = ['O(nL)\n4^k features', 
              'O(n²)\nDistance matrix',
              'O(n²) + classifier\n50 features',
              'O(m²L²) per cluster\nMSA + ML']

# Draw boxes for each level
for i, (y, name, color, info) in enumerate(zip(y_levels, level_names, colors, info_labels)):
    # Main box
    box = mpatches.FancyBboxPatch(
        (x_start, y - box_height/2), box_width, box_height,
        boxstyle="round,pad=0.1", 
        facecolor=color, edgecolor='black', linewidth=2, alpha=0.7
    )
    ax.add_patch(box)
    
    # Text inside box
    ax.text(x_start + box_width/2, y, name, 
            ha='center', va='center', fontsize=11, fontweight='bold', 
            color='white', zorder=10)
    
    # Info label below box
    ax.text(x_start + box_width/2, y - box_height/2 - 0.3, info,
            ha='center', va='top', fontsize=8, style='italic',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray', alpha=0.8))

# Draw arrows between levels
arrow_props = dict(arrowstyle='->', lw=2.5, color='black')
arrow_props_conditional = dict(arrowstyle='->', lw=2.5, color='red', linestyle='--')

# L0 -> L1
ax.annotate('', xy=(x_start, y_levels[1] + box_height/2), 
            xytext=(x_start + box_width, y_levels[0] - box_height/2),
            arrowprops=arrow_props)

# L1 -> L2
ax.annotate('', xy=(x_start, y_levels[2] + box_height/2), 
            xytext=(x_start + box_width, y_levels[1] - box_height/2),
            arrowprops=arrow_props)

# L2 -> L3 (conditional, dashed red)
ax.annotate('Conditional\n(escalation)', 
            xy=(x_start, y_levels[3] + box_height/2), 
            xytext=(x_start + box_width, y_levels[2] - box_height/2),
            arrowprops=arrow_props_conditional,
            fontsize=10, color='red', fontweight='bold')

# Add information gradient bar on the right
gradient_x = 5
gradient_width = 0.5
for i, (y, color) in enumerate(zip(y_levels, colors)):
    rect = plt.Rectangle((gradient_x, y - box_height/2), gradient_width, box_height,
                          facecolor=color, edgecolor='black', linewidth=1)
    ax.add_patch(rect)
    # Label
    info_text = ['Low', 'Medium-Low', 'Medium-High', 'High'][i]
    ax.text(gradient_x + gradient_width + 0.2, y, f'Info: {info_text}',
            va='center', fontsize=9, fontweight='bold')

# Add title and formatting
ax.set_xlim([0, 7])
ax.set_ylim([0.5, 10])
ax.set_title('Figure 1. IMMI Framework Architecture\n(Information-Matched Multi-Level Inference)', 
             fontsize=16, fontweight='bold', pad=30)

# Add legend for notation
legend_elements = [
    mpatches.Patch(color=colors[0], label='L0: k-mer features (O(nL))'),
    mpatches.Patch(color=colors[1], label='L1: Distance + NJ (O(n²))'),
    mpatches.Patch(color=colors[2], label='L2: Classifier (learned)'),
    mpatches.Patch(color=colors[3], label='L3: MSA+ML (O(m²L²))'),
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=10, 
          framealpha=0.9, edgecolor='black')

ax.axis('off')

plt.tight_layout()
output_path = BASE_DIR / "Figure1_IMMI_architecture.png"
plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"Saved: {output_path}")

pdf_path = BASE_DIR / "Figure1_IMMI_architecture.pdf"
plt.savefig(pdf_path, bbox_inches='tight', facecolor='white')
print(f"Saved: {pdf_path}")

plt.close()
print("\nFigure 1 generation complete!")
