"""
Fusang MHL (Multi-Level Hierarchical) Phylogenetic Tree Construction Framework

Part of Fusang: Tardigrade Edition.

Architecture:
  Level 0: k-mer cosine distance + hierarchical clustering (pure AF, O(n^2))
  Level 1: multi-k ensemble + boundary classifier decision
  Level 2: DAHP-V2 centroid backbone + optional NNI
  Level 3: MSA+ML within MRC (MAFFT + FastTree2)
  Merge: Bridge Taxa constrained supertree assembly

MHL = Multi-level Hierarchical, information-matched inference.
"""

__version__ = "0.1.0"
__author__ = "Fusang Team"

from .config import (
    MAFFT_BAT, FASTTREE_EXE, MAFFT_DIR, MAFFT_TMP,
    BOUNDARY_MODEL_DIR, TEMP_DIR,
    L0_THRESHOLDS, L1_DEFAULTS, L2_DEFAULTS, L3_DEFAULTS,
    MERGE_DEFAULTS,
)
from .mlh_utils import Timer, ProgressReporter, setup_logger

__all__ = [
    "MAFFT_BAT", "FASTTREE_EXE", "MAFFT_DIR", "MAFFT_TMP",
    "BOUNDARY_MODEL_DIR", "TEMP_DIR",
    "L0_THRESHOLDS", "L1_DEFAULTS", "L2_DEFAULTS", "L3_DEFAULTS",
    "MERGE_DEFAULTS",
    "Timer", "ProgressReporter", "setup_logger",
]
