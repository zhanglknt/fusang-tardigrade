"""
ML-based adaptive splitting for the MHL pipeline.

Replaces hardcoded L0_THRESHOLDS with a trained boundary classifier
that decides whether to split a cluster based on 50-dim features
extracted from the k-mer distance matrix and sequence data.

The model (`boundary_rf.pkl`) is trained on simulated trees with
diverse substitution rates, indel rates, and taxon counts.
"""

import sys
import os
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional

_MODEL_DIR = Path(__file__).parent / "models"
_MODEL_PATH = _MODEL_DIR / "boundary_rf.pkl"
_V2_PATH = _MODEL_DIR / "boundary_rf_v2.pkl"   # V2 with indel training data
_V3_PATH = _MODEL_DIR / "boundary_rf_v3.pkl"    # V3 with whole-tree samples (L0-aligned)
_V4b_PATH = _MODEL_DIR / "boundary_rf_v4b.pkl"  # V4b: structured trees, ground-truth whole-tree labels, 6707 samples

_model = None
_model_path_used = None


def _resolve_model_path() -> Path:
    """Return the best available model path (prefer V4b > V3 > V2 > V1)."""
    if _V4b_PATH.exists():
        return _V4b_PATH
    if _V3_PATH.exists():
        return _V3_PATH
    if _V2_PATH.exists():
        return _V2_PATH
    if _MODEL_PATH.exists():
        return _MODEL_PATH
    return None


def _load_model():
    """Lazy-load the boundary classifier via BoundaryClassifier.load()."""
    global _model, _model_path_used
    model_path = _resolve_model_path()
    if model_path is None:
        return None
    if _model is None or _model_path_used != str(model_path):
        from .boundary_classifier import BoundaryClassifier
        _model = BoundaryClassifier().load(str(model_path))
        _model_path_used = str(model_path)
    return _model


def model_available() -> bool:
    """Check if a boundary classifier model file exists."""
    return _resolve_model_path() is not None


def _extract_l0_features(
    D: np.ndarray,
    sequences: List[str],
    n: int,
) -> Optional[np.ndarray]:
    """Extract 50-dim features treating the full dataset as a single cluster.

    Uses the same feature extraction as training_data_generator to ensure
    distribution consistency between training and inference.
    """
    try:
        from .level1_multik import extract_cluster_features
        features = extract_cluster_features(
            D=D, D_k5=None, D_k7=None, D_k9=None,
            seqs_for_cluster=sequences,
            cluster_indices=list(range(n)),
            centroid_idx=0,
            feature_matrix=None,
            n_total=n,
            parent_size=n,
            sibling_sizes=[],
            current_level=0,
            ancestor_sizes=[],
        )
        if features is not None and len(features) == 50:
            return features
    except Exception:
        pass
    return None


def _estimate_split_parameters(n: int, p_split: float) -> Dict:
    """Estimate target_groups and max_group_size from prediction confidence.

    Higher p_split → more groups (dataset is more heterogeneous).
    """
    # Scale target_groups with both n and confidence
    if p_split > 0.9:
        # Very confident: aggressive splitting
        n_groups_factor = 0.18
    elif p_split > 0.7:
        n_groups_factor = 0.12
    else:
        # Borderline: conservative
        n_groups_factor = 0.07

    target_groups = max(2, int(n * n_groups_factor / 50))
    target_groups = min(target_groups, 25)  # Cap at 25 groups

    # max_group_size: inverse of target_groups but with floor
    max_group_size = max(50, n // target_groups)

    return {
        "target_groups": target_groups,
        "max_group_size": max_group_size,
    }


def ml_split_decision(
    D: np.ndarray,
    sequences: List[str],
    taxon_names: List[str],
    verbose: bool = False,
) -> Dict:
    """Use boundary classifier to decide L0 split strategy.

    Args:
        D: n x n k-mer distance matrix
        sequences: List of DNA sequences
        taxon_names: List of taxon names
        verbose: Print decision details

    Returns:
        dict with keys:
            should_split: bool
            target_groups: int
            max_group_size: int
            p_split: float (classifier confidence, None if model unavailable)
            model_used: bool
    """
    n = len(taxon_names)

    # Hard floor: too few taxa to split meaningfully
    if n <= 30:
        return {
            "should_split": False, "target_groups": 1,
            "max_group_size": n, "p_split": None, "model_used": False,
            "reason": "n_too_small",
        }

    # Try ML model
    model = _load_model()
    if model is None:
        # Fallback to heuristic
        if verbose:
            print("[ML-split] No model found, using heuristic fallback",
                  file=sys.stderr)
        return _heuristic_fallback(n)

    features = _extract_l0_features(D, sequences, n)
    if features is None:
        if verbose:
            print("[ML-split] Feature extraction failed, using heuristic",
                  file=sys.stderr)
        return _heuristic_fallback(n)

    # Predict
    try:
        proba = model.predict_proba([features])[0]
        p_split = float(proba[1])
    except Exception as e:
        if verbose:
            print(f"[ML-split] Prediction failed: {e}", file=sys.stderr)
        return _heuristic_fallback(n)

    should_split = p_split > 0.5

    if should_split:
        params = _estimate_split_parameters(n, p_split)
    else:
        params = {"target_groups": 1, "max_group_size": n}

    if verbose:
        decision = "SPLIT" if should_split else "STOP"
        print(f"[ML-split] n={n} → {decision} "
              f"(p_split={p_split:.3f}, groups={params['target_groups']}, "
              f"max_size={params['max_group_size']})",
              file=sys.stderr)

    return {
        "should_split": should_split,
        "target_groups": params["target_groups"],
        "max_group_size": params["max_group_size"],
        "p_split": p_split,
        "model_used": True,
        "reason": "ml_prediction",
    }


def _heuristic_fallback(n: int) -> Dict:
    """Heuristic split decision when ML model is unavailable."""
    if n <= 200:
        return {
            "should_split": False, "target_groups": 1,
            "max_group_size": n, "p_split": None, "model_used": False,
            "reason": "heuristic_n_small",
        }
    elif n <= 500:
        tg = max(2, n // 150)
        return {
            "should_split": True, "target_groups": tg,
            "max_group_size": 150, "p_split": None, "model_used": False,
            "reason": "heuristic_n_medium",
        }
    elif n <= 1000:
        tg = max(3, n // 200)
        return {
            "should_split": True, "target_groups": tg,
            "max_group_size": 200, "p_split": None, "model_used": False,
            "reason": "heuristic_n_large",
        }
    else:
        tg = max(5, n // 250)
        return {
            "should_split": True, "target_groups": min(tg, 25),
            "max_group_size": 250, "p_split": None, "model_used": False,
            "reason": "heuristic_n_xl",
        }
