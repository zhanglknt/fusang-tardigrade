#!/usr/bin/env python
"""
Fusang Web Server - Scalable Phylogenetic Inference
====================================================

Setup script for pip/conda installation.
After installation, run with: fusang-web
"""

import os
import sys
from pathlib import Path
from setuptools import setup, find_packages

BASE_DIR = Path(__file__).parent

# Read requirements
with open(BASE_DIR / "requirements.txt", encoding="utf-8") as f:
    requirements = [line.strip() for line in f if line.strip() and not line.startswith("#")]

# Find all root-level .py modules (exclude test/benchmark scripts)
root_py_modules = []
exclude_patterns = {"setup", "apply_fix", "apply_clean_fix", "apply_line_fix",
                    "benchmark_competitors", "batch_iqtree2_n200",
                    "train_boundary_classifier", "af_competitor_methods"}

for f in BASE_DIR.glob("*.py"):
    name = f.stem
    if name not in exclude_patterns and not name.startswith("_") and not name.startswith("."):
        root_py_modules.append(name)

setup(
    name="fusang-web",
    version="1.0.0",
    author="Lei Kong, Li Zhang",
    author_email="zhangli@cibr.ac.cn",
    description="Fusang: Scalable Phylogenetic Inference via Divide-and-Conquer — Web Server",
    long_description=(BASE_DIR / "NAR_MANUSCRIPT_IMMI_NAR_FORMAT.md").read_text(encoding="utf-8")[:5000] if (BASE_DIR / "NAR_MANUSCRIPT_IMMI_NAR_FORMAT.md").exists() else "Fusang Web Server",
    long_description_content_type="text/plain",
    url="https://github.com/Fusang/Fusang",
    packages=find_packages(),
    py_modules=root_py_modules,
    include_package_data=True,
    install_requires=requirements,
    python_requires=">=3.8",
    entry_points={
        "console_scripts": [
            "fusang-web=fusang_webapp:main",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    keywords="phylogenetics bioinformatics tree-inference k-mer divide-and-conquer",
)
