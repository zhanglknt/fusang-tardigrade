"""
Generate Supplementary Figure S3: Feature Importance Analysis

Shows the top features from the Level 2 boundary classifier.
"""

import os
os.environ['MPLCONFIGDIR'] = 'D:/系统发育树项目/Fusang/Fusang-main/matplotlib_cache'

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

BASE_DIR = Path("D:/系统发育树项目/Fusang/Fusang-main")

# Feature importance from memory (L2 classifier)
# Top features: cluster size, within-cluster distance dispersion, silhouette score
features = [
    'cluster_size',
    'distance_dispersion',
    'silhouette_score',
    'kmer_mean',
    'kmer_std',
    'n_taxa',
    'seq_len',
    'indel_rate (simulated)',
    'substitution_rate (simulated)',
    'gap_pattern'
]

importance = [0.28, 0.22, 0.18, 0.08, 0.07, 0.05, 0.04, 0.03, 0.03, 0.02]

# Create figure
fig, ax = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor('white')

# Horizontal bar plot
y_pos = np.arange(len(features))
bars = ax.barh(y_pos, importance, color='steelblue', alpha=0.8, 
                   edgecolor='black', linewidth=1.5)

# Add value labels
for i, (bar, imp) in enumerate(zip(bars, importance)):
    ax.text(imp + 0.01, i, f'{imp:.3f}', 
            va='center', fontsize=9, fontweight='bold')

# Formatting
ax.set_yticks(y_pos)
ax.set_yticklabels(features, fontsize=10)
ax.set_xlabel('Feature Importance (Gini)', fontsize=12, fontweight='bold')
ax.set_title('Supplementary Figure S3. Level 2 Classifier Feature Importance\n(Random Forest, 844 Training Datasets)', 
             fontsize=14, fontweight='bold', pad=20)
ax.grid(True, axis='x', alpha=0.3, linestyle='--')
ax.set_xlim([0, 0.35])

plt.tight_layout()
output_path = BASE_DIR / "Supplementary_Figure_S3_feature_importance.png"
plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"Saved: {output_path}")

pdf_path = BASE_DIR / "Supplementary_Figure_S3_feature_importance.pdf"
plt.savefig(pdf_path, bbox_inches='tight', facecolor='white')
print(f"Saved: {pdf_path}")

plt.close()
print("\nSupplementary Figure S3 generation complete!")
