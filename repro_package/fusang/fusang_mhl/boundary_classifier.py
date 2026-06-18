"""
Boundary Classifier: sklearn Random Forest for split/stop decision.

Trains on simulated data to predict whether a cluster should be
further subdivided (split) or is homogeneous enough (stop).
"""

import sys
import os
import pickle
import numpy as np
from typing import Dict, List, Tuple, Optional, Any

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from .config import BOUNDARY_MODEL_DIR
from .mlh_utils import Timer


class BoundaryClassifier:
    """Random Forest classifier for cluster split/stop decision."""

    def __init__(
        self,
        n_estimators: int = 200,
        max_depth: Optional[int] = 15,
        class_weight: str = "balanced",
        random_state: int = 42,
    ):
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.class_weight = class_weight
        self.random_state = random_state
        self.model = None
        self.scaler = None
        self.feature_names: Optional[List[str]] = None
        self.is_fitted = False

    def _get_model(self):
        """Lazily create the RF model."""
        from sklearn.ensemble import RandomForestClassifier
        return RandomForestClassifier(
            n_estimators=self.n_estimators,
            max_depth=self.max_depth,
            class_weight=self.class_weight,
            random_state=self.random_state,
            n_jobs=-1,
        )

    def fit(
            self,
            X: np.ndarray,
            y: np.ndarray,
            feature_names: Optional[List[str]] = None,
    ) -> "BoundaryClassifier":
        """Fit the classifier.

        Args:
            X: (n_samples, n_features) feature matrix
            y: (n_samples,) labels (0=stop, 1=split)
            feature_names: Optional feature name list
        """
        from sklearn.preprocessing import StandardScaler

        self.scaler = StandardScaler()
        X_scaled = self.scaler.fit_transform(X)

        self.model = self._get_model()
        self.model.fit(X_scaled, y)

        self.feature_names = feature_names
        self.is_fitted = True
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        """Predict split/stop for each sample."""
        if not self.is_fitted:
            raise RuntimeError("Model not fitted yet. Call fit() first.")
        X_s = self.scaler.transform(X)
        return self.model.predict(X_s)

    def predict_proba(self, X: np.ndarray) -> np.ndarray:
        """Predict class probabilities."""
        if not self.is_fitted:
            raise RuntimeError("Model not fitted yet.")
        X_s = self.scaler.transform(X)
        return self.model.predict_proba(X_s)

    def predict_split_decision(
            self,
            X: np.ndarray,
            threshold: float = 0.5,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Predict split decision with probability threshold.

        Returns:
            (decisions, probabilities)
            decisions: bool array, True=split
            probabilities: P(split) array
        """
        proba = self.predict_proba(X)
        p_split = proba[:, 1]  # probability of class 1 (split)
        decisions = p_split > threshold
        return decisions, p_split

    def evaluate(
            self,
            X_test: np.ndarray,
            y_test: np.ndarray,
    ) -> Dict[str, float]:
        """Evaluate on test set."""
        from sklearn.metrics import (
            accuracy_score, precision_score, recall_score,
            f1_score, roc_auc_score,
        )
        if not self.is_fitted:
            raise RuntimeError("Model not fitted yet.")

        y_pred = self.predict(X_test)
        proba = self.predict_proba(X_test)

        return {
            "accuracy": float(accuracy_score(y_test, y_pred)),
            "precision": float(precision_score(y_test, y_pred, zero_division=0)),
            "recall": float(recall_score(y_test, y_pred, zero_division=0)),
            "f1": float(f1_score(y_test, y_pred, zero_division=0)),
            "roc_auc": float(roc_auc_score(y_test, proba[:, 1])),
        }

    def cross_validate(
            self,
            X: np.ndarray,
            y: np.ndarray,
            n_splits: int = 10,
            stratified: bool = True,
    ) -> Dict[str, Any]:
        """Run cross-validation."""
        from sklearn.model_selection import (
            cross_validate as sk_cv,
            StratifiedKFold, KFold,
        )
        if stratified:
            cv = StratifiedKFold(
                n_splits=n_splits,
                shuffle=True,
                random_state=self.random_state,
            )
        else:
            cv = KFold(
                n_splits=n_splits,
                shuffle=True,
                random_state=self.random_state,
            )

        self.model = self._get_model()
        scoring = ["accuracy", "precision", "recall", "f1", "roc_auc"]
        cv_results = sk_cv(
            self.model, self.scaler.fit_transform(X) if self.scaler else X,
            y,
            cv=cv,
            scoring=scoring,
            n_jobs=1,
        )
        # Fit on full data
        self.fit(X, y)

        return {
            "accuracy": {
                "mean": float(np.mean(cv_results["test_accuracy"])),
                "std": float(np.std(cv_results["test_accuracy"])),
            },
            "precision": {
                "mean": float(np.mean(cv_results["test_precision"])),
                "std": float(np.std(cv_results["test_precision"])),
            },
            "recall": {
                "mean": float(np.mean(cv_results["test_recall"])),
                "std": float(np.std(cv_results["test_recall"])),
            },
            "f1": {
                "mean": float(np.mean(cv_results["test_f1"])),
                "std": float(np.std(cv_results["test_f1"])),
            },
            "roc_auc": {
                "mean": float(np.mean(cv_results["test_roc_auc"])),
                "std": float(np.std(cv_results["test_roc_auc"])),
            },
        }

    def feature_importance(self) -> Optional[np.ndarray]:
        """Get feature importances."""
        if not self.is_fitted or self.model is None:
            return None
        return self.model.feature_importances_

    def print_feature_importance(self, top_n: int = 20):
        """Print top feature importances."""
        if not self.is_fitted:
            print("Model not fitted.")
            return
        imp = self.feature_importance()
        if imp is None:
            return
        indices = np.argsort(imp)[::-1]
        names = self.feature_names or [f"f{i}" for i in range(len(imp))]
        print(f"{'Rank':<6} {'Feature':<40} {'Importance':>10}")
        for i in range(min(top_n, len(imp))):
            idx = indices[i]
            print(f"{i+1:<6} {names[idx]:<40} {imp[idx]:>10.4f}")

    def save(self, path: Optional[str] = None):
        """Save model to file."""
        if not self.is_fitted:
            raise RuntimeError("Model not fitted yet.")
        path = path or os.path.join(BOUNDAY_MODEL_DIR, "boundary_rf.pkl")
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "wb") as f:
            pickle.dump({
                "model": self.model,
                "scaler": self.scaler,
                "feature_names": self.feature_names,
                "n_estimators": self.n_estimators,
                "max_depth": self.max_depth,
                "class_weight": self.class_weight,
            }, f)
        print(f"[BoundaryClassifier] Saved to {path}")

    def load(self, path: Optional[str] = None):
        """Load model from file."""
        path = path or os.path.join(BOUNDAY_MODEL_DIR, "boundary_rf.pkl")
        with open(path, "rb") as f:
            data = pickle.load(f)
        self.model = data["model"]
        self.scaler = data["scaler"]
        self.feature_names = data.get("feature_names")
        self.n_estimators = data.get("n_estimators", 200)
        self.max_depth = data.get("max_depth", 15)
        self.class_weight = data.get("class_weight", "balanced")
        self.is_fitted = True
        print(f"[BoundaryClassifier] Loaded from {path}")
        return self


def train_boundary_classifier(
        X: np.ndarray,
        y: np.ndarray,
        test_size: float = 0.2,
        random_state: int = 42,
        verbose: bool = True,
) -> Tuple[BoundaryClassifier, Dict[str, Any]]:
    """Train a boundary classifier with train/test split and CV.

    Args:
        X: Feature matrix
        y: Labels
        test_size: Test set fraction
        random_state: Random seed
        verbose: Print progress

    Returns:
        (trained classifier, CV results dict)
    """
    from sklearn.model_selection import train_test_split

    if verbose:
        print(f"[Train] X shape: {X.shape}, y dist: {np.bincount(y.astype(int))}")

    X_train, X_test, y_train, y_test = train_test_split(
        X, y,
        test_size=test_size,
        stratify=y,
        random_state=random_state,
    )

    clf = BoundaryClassifier()
    cv_results = clf.cross_validate(X_train, y_train, n_splits=10, stratified=True)

    if verbose:
        print("[Train] CV Results:")
        for metric, vals in cv_results.items():
            print(f"  {metric}: {vals['mean']:.4f} +/- {vals['std']:.4f}")

    # Final evaluation on held-out test set
    test_metrics = clf.evaluate(X_test, y_test)
    if verbose:
        print("[Train] Test set results:")
        for metric, val in test_metrics.items():
            print(f"  {metric}: {val:.4f}")

    return clf, {"cv": cv_results, "test": test_metrics}
