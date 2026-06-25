#!/usr/bin/env python3
"""
Train the boundary classifier for Fusang MHL.

Usage:
    python train_boundary_classifier.py [--data training_data.pkl] [--output model_dir]
"""

import sys
import os
import argparse

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV, StratifiedKFold

from fusang_mhl.boundary_classifier import (
    BoundaryClassifier, train_boundary_classifier,
)
from fusang_mhl.training_data_generator import (
    generate_training_data, load_training_data, split_training_data,
)
from fusang_mhl.config import BOUNDARY_MODEL_DIR
from fusang_mhl.mlh_utils import Timer, ensure_dir, setup_logger


def main():
    parser = argparse.ArgumentParser(
        description="Train boundary classifier for Fusang MHL"
    )
    parser.add_argument(
        "--data", type=str, default="fusang_mhl/models/training_data_v2.pkl",
        help="Pre-generated training data (pkl). Default: V2 with indel support.",
    )
    parser.add_argument(
        "--output", type=str, default="fusang_mhl/models",
        help="Output directory for model files",
    )
    parser.add_argument(
        "--n-samples", type=int, default=10000,
        help="Target number of training samples (if generating)",
    )
    parser.add_argument(
        "--model-name", type=str, default="boundary_rf_v2.pkl",
        help="Output model filename",
    )
    parser.add_argument(
        "--grid-search", action="store_true",
        help="Run GridSearchCV for hyperparameter tuning",
    )
    parser.add_argument(
        "--verbose", action="store_true", default=True,
    )
    args = parser.parse_args()

    ensure_dir(args.output)
    logger = setup_logger("train_bc", log_file=os.path.join(args.output, "train.log"))

    # Load or generate training data
    if args.data and os.path.exists(args.data):
        logger.info(f"Loading training data from {args.data}")
        data = load_training_data(args.data)
    else:
        logger.info("Generating training data...")
        pkl_path = os.path.join(args.output, "training_data.pkl")
        data = generate_training_data(
            output_pkl=pkl_path,
            n_samples_target=args.n_samples,
            verbose=args.verbose,
        )

    n_pos = sum(1 for s in data if s["label"] == 1)
    n_neg = sum(1 for s in data if s["label"] == 0)
    logger.info(f"Total samples: {len(data)} (pos={n_pos}, neg={n_neg})")

    # Split (used for grid-search path; default path does its own split)
    train_data, test_data = split_training_data(data, test_size=0.2)
    logger.info(f"Train: {len(train_data)}, Test: {len(test_data)}")

    X_train = np.array([s["features"] for s in train_data])
    y_train = np.array([s["label"] for s in train_data])
    X_test = np.array([s["features"] for s in test_data])
    y_test = np.array([s["label"] for s in test_data])

    logger.info(f"Feature dim: {X_train.shape[1]}")

    if args.grid_search:
        logger.info("Running GridSearchCV...")
        param_grid = {
            "n_estimators": [100, 200, 500],
            "max_depth": [10, 15, 20, None],
            "min_samples_leaf": [1, 5, 10],
            "class_weight": ["balanced", None],
        }
        clf = BoundaryClassifier()
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X_train)

        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        grid = GridSearchCV(
            clf.model, param_grid, cv=cv,
            scoring="roc_auc", n_jobs=-1, verbose=2,
        )
        grid.fit(X_scaled, y_train)

        logger.info(f"Best params: {grid.best_params_}")
        logger.info(f"Best CV ROC-AUC: {grid.best_score_:.4f}")

        # Use best params to train final model
        best_clf = BoundaryClassifier(**grid.best_params_)
        with Timer("Training (best params)", verbose=True):
            best_clf.fit(X_train, y_train)
        test_acc = best_clf.evaluate(X_test, y_test)
        logger.info(f"Test accuracy: {test_acc}")

        # Save
        model_path = os.path.join(args.output, "boundary_rf_best.pkl")
        scaler_path = os.path.join(args.output, "scaler_best.pkl")
        best_clf.save(model_path)
        with open(scaler_path, "wb") as f:
            import pickle
            pickle.dump(scaler, f)
        logger.info(f"Saved model to {model_path}")
    else:
        # Default training
        logger.info("Training with default parameters...")
        X_all = np.array([s["features"] for s in data])
        y_all = np.array([s["label"] for s in data])
        clf, metrics = train_boundary_classifier(
            X_all, y_all,
            test_size=0.2,
            verbose=args.verbose,
        )
        # Save model
        model_path = os.path.join(args.output, args.model_name)
        clf.save(model_path)
        logger.info(f"Model: {model_path}")
        logger.info(f"Metrics: {metrics}")

    logger.info("Done.")


if __name__ == "__main__":
    main()
