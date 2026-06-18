"""
Generate Figure 6: Confusion Matrix (Classifier Performance)

Shows the Level 2 boundary classifier performance.
"""

import os
os.environ['MPLCONFIGDIR'] = 'D:/系统发育树项目/Fusang/Fusang-main/matplotlib_cache'

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

BASE_DIR = Path("D:/系统发育树项目/Fusang/Fusang-main")

# Confusion matrix from memory (Table 6)
# Test set: 168 samples (92 pos, 76 neg)
# TP=89, TN=71, FP=5, FN=3 (from corrected Figure 6)
confusion_matrix = np.array([[89, 3],   # TP, FN
                              [5, 71]])  # FP, TN

# Create figure
fig, ax = plt.subplots(figsize=(8, 6))
fig.patch.set_facecolor('white')

# Plot confusion matrix as heatmap
im = ax.imshow(confusion_matrix, cmap='Blues', interpolation='nearest', aspect='auto')

# Add colorbar
cbar = plt.colorbar(im, ax=ax)
cbar.set_label('Count', rotation=270, labelpad=20, fontsize=12)

# Add text annotations
for i in range(2):
    for j in range(2):
        ax.text(j, i, str(confusion_matrix[i, j]),
                ha='center', va='center', fontsize=20, fontweight='bold',
                color='white' if confusion_matrix[i, j] > 50 else 'black')

# Set labels
ax.set_xticks([0, 1])
ax.set_xticklabels(['Predicted: Split', 'Predicted: No Split'], fontsize=11)
ax.set_yticks([0, 1])
ax.set_yticklabels(['Actual: Split', 'Actual: No Split'], fontsize=11)

# Add metrics as text
accuracy = (89 + 71) / (89 + 71 + 5 + 3)
precision = 89 / (89 + 5)
recall = 89 / (89 + 3)
f1 = 2 * precision * recall / (precision + recall)

metrics_text = f'Accuracy: {accuracy:.3f}\nPrecision: {precision:.3f}\nRecall: {recall:.3f}\nF1: {f1:.3f}'
ax.text(1.5, 0.5, metrics_text, ha='left', va='center', fontsize=11,
         bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor='gray'))

# Formatting
ax.set_title('Figure 6. Level 2 Boundary Classifier Performance\n(Confusion Matrix, 168 Test Samples)', 
             fontsize=14, fontweight='bold', pad=20)
ax.grid(False)

plt.tight_layout()
output_path = BASE_DIR / "Figure6_confusion_matrix.png"
plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"Saved: {output_path}")

pdf_path = BASE_DIR / "Figure6_confusion_matrix.pdf"
plt.savefig(pdf_path, bbox_inches='tight', facecolor='white')
print(f"Saved: {pdf_path}")

plt.close()
print("\nFigure 6 generation complete!")
