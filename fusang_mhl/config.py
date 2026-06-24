"""
Configuration constants for Fusang MHL framework.
"""

import os
import tempfile
from pathlib import Path

# ============================================================
# External Tool Paths (Cross-Platform)
# ============================================================
# Auto-detect FUSANG_ROOT relative to this config file
FUSANG_ROOT = Path(__file__).resolve().parent.parent  # fusang_mhl -> project root

# External tools: use environment variables, with sensible defaults
MAFFT_DIR = Path(os.environ.get("MAFFT_DIR", str(FUSANG_ROOT.parent / "bench_tools" / "mafft-win" / "mafft-win")))
MAFFT_BAT = os.environ.get("MAFFT_BAT", os.environ.get("MAFFT", str(MAFFT_DIR / "mafft.bat")))
MAFFT_TMP = MAFFT_DIR / "tmp"
FASTTREE_EXE = os.environ.get("FASTTREE_EXE", os.environ.get("FASTTREE", str(FUSANG_ROOT.parent / "bench_tools" / "FastTree.exe")))

# ============================================================
# Model and Temp Directories
# ============================================================
BOUNDARY_MODEL_DIR = str(FUSANG_ROOT / "dl_model_boundary")
TEMP_DIR = os.path.join(tempfile.gettempdir(), "fusang_mhl")

# ============================================================
# Level 0: k-mer Clustering Configuration
# ============================================================
# Adaptive group sizing based on total taxa count
# Format: (min_n, max_group_size, target_groups)
L0_THRESHOLDS = [
    (0,    200, 1),     # n<=200: no L0 splitting (single group)
    (200,  150, 3),     # n=200-500: max_group=150, ~2-3 groups
    (500,  200, 5),     # n=500-1000: max_group=200, ~3-5 groups
    (1000, 250, 8),     # n=1000-2000: max_group=250, ~4-8 groups
    (2000, 300, 15),    # n>2000: max_group=300, ~7-15 groups
]

# Default k-mer parameters for L0
L0_K = 5
L0_GAP = "gap2"  # spaced k-mer pattern (0,1,2,5,6)
L0_DISTANCE_METHOD = "cosine"

# ============================================================
# Level 1: Multi-k Ensemble + Boundary Classifier
# ============================================================
L1_DEFAULTS = {
    "ks": (5, 7, 9),                    # multi-k values
    "gap": "contiguous",                 # contiguous for ensemble
    "fusion_method": "average",          # average distance fusion
    "split_threshold": 0.5,             # P(split) > threshold -> split
    "min_cluster_size": 10,             # minimum cluster size to consider splitting
    "distance_method": "cosine",
}

# ============================================================
# Level 2: DAHP-V2 Centroid Backbone
# ============================================================
L2_DEFAULTS = {
    "clade_threshold": 0.03,            # cosine distance threshold for clade cutting
    "clade_threshold_large": 0.02,      # stricter for clusters >100 taxa
    "large_cluster_threshold": 100,
    "min_clade_size": 4,                 # minimum taxa for DAHP-V2
    "k": 5,
    "gap": "gap2",
    "enable_nni": False,                 # optional NNI refinement
    "nni_max_rounds": 10,
}

# ============================================================
# Level 3: MSA+ML within MRC
# ============================================================
L3_DEFAULTS = {
    "max_pairwise_dist": 0.01,          # cosine distance threshold for MRC
    "min_taxa": 4,                       # minimum taxa for MSA+ML
    "max_taxa": 200,                     # maximum taxa per MRC (computational limit)
    "mafft_mode": "--auto",              # MAFFT mode for small alignments
    "fasttree_model": "-gtr",            # GTR model
}

# ============================================================
# Merger Configuration
# ============================================================
MERGE_DEFAULTS = {
    "n_bridge_taxa": 3,                  # bridge taxa per subtree
    "bridge_selection": "centroid+random", # how to select bridge taxa
    "constrain_monophyly": True,         # enforce monophyly constraints
}

# ============================================================
# Benchmark Configuration
# ============================================================
BENCHMARK_SCALES = {
    "small":  {"n": 200,  "seeds": 30, "L": 500, "sub": 0.05, "indel": 0.02},
    "medium": {"n": 500,  "seeds": 30, "L": 500, "sub": 0.05, "indel": 0.02},
    "large":  {"n": 1000, "seeds": 20, "L": 500, "sub": 0.05, "indel": 0.01},
    "xl":     {"n": 2000, "seeds": 10, "L": 500, "sub": 0.03, "indel": 0.01},
    "xxl":    {"n": 5000, "seeds": 5,  "L": 500, "sub": 0.02, "indel": 0.01},
}

# ============================================================
# Normalization
# ============================================================
NRF_NORMALIZATION = "2*(n-3)"  # nRF = RF / (2*(n-3))

# ============================================================
# Logging
# ============================================================
LOG_FORMAT = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
LOG_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"
