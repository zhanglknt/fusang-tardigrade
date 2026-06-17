"""
Utility functions for Fusang MHL framework.
"""

import sys
import time
import logging
import os
from contextlib import contextmanager
from typing import Optional


class Timer:
    """Context manager for timing code blocks."""

    def __init__(self, name: str = "", verbose: bool = True):
        self.name = name
        self.verbose = verbose
        self.elapsed = 0.0
        self._start = 0.0

    def __enter__(self):
        self._start = time.perf_counter()
        return self

    def __exit__(self, *args):
        self.elapsed = time.perf_counter() - self._start
        if self.verbose:
            print(f"  [{self.name}] {self.elapsed:.2f}s", file=sys.stderr)

    @staticmethod
    def format_time(seconds: float) -> str:
        """Format seconds into human-readable string."""
        if seconds < 60:
            return f"{seconds:.1f}s"
        elif seconds < 3600:
            m, s = divmod(seconds, 60)
            return f"{int(m)}m {s:.0f}s"
        else:
            h, rem = divmod(seconds, 3600)
            m, s = divmod(rem, 60)
            return f"{int(h)}h {int(m)}m {s:.0f}s"


class ProgressReporter:
    """Simple progress reporter for iterations."""

    def __init__(self, total: int, desc: str = "", bar_width: int = 40):
        self.total = total
        self.desc = desc
        self.bar_width = bar_width
        self._count = 0
        self._start = time.perf_counter()

    def update(self, n: int = 1):
        self._count += n
        self._print_bar()

    def _print_bar(self):
        frac = self._count / self.total if self.total > 0 else 0
        filled = int(self.bar_width * frac)
        bar = "#" * filled + "-" * (self.bar_width - filled)
        elapsed = time.perf_counter() - self._start
        rate = self._count / elapsed if elapsed > 0 else 0
        eta = (self.total - self._count) / rate if rate > 0 else 0
        sys.stderr.write(
            f"\r  {self.desc}: [{bar}] {self._count}/{self.total} "
            f"({frac:.0%}) {Timer.format_time(elapsed)} elapsed, "
            f"ETA {Timer.format_time(eta)}   "
        )
        sys.stderr.flush()

    def finish(self):
        self._count = self.total
        self._print_bar()
        sys.stderr.write("\n")
        sys.stderr.flush()


def setup_logger(
    name: str = "fusang_mhl",
    level: int = logging.INFO,
    log_file: Optional[str] = None,
) -> logging.Logger:
    """Create and configure a logger."""
    logger = logging.getLogger(name)
    logger.setLevel(level)

    # Console handler
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(logging.Formatter(
        "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S"
    ))
    logger.addHandler(handler)

    # Optional file handler
    if log_file:
        os.makedirs(os.path.dirname(log_file) or ".", exist_ok=True)
        fh = logging.FileHandler(log_file, encoding="utf-8")
        fh.setFormatter(logging.Formatter(
            "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S"
        ))
        logger.addHandler(fh)

    return logger


@contextmanager
def temp_directory(base: str):
    """Create and clean up a temporary directory."""
    import tempfile
    import shutil
    d = tempfile.mkdtemp(dir=base)
    try:
        yield d
    finally:
        try:
            shutil.rmtree(d, ignore_errors=True)
        except Exception:
            pass


def ensure_dir(path: str):
    """Ensure a directory exists."""
    os.makedirs(path, exist_ok=True)


def read_fasta_simple(path: str, encoding: str = "utf-8") -> dict:
    """Simple FASTA reader. Returns {name: sequence}."""
    seqs = {}
    name = None
    VALID = set("ACGTURYSWKMBDHVN.-*")
    with open(path, "r", encoding=encoding, errors="replace") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(
                    "".join(c for c in line.upper() if c in VALID)
                )
    return {k: "".join(v) for k, v in seqs.items()}


def write_fasta(seqs: dict, path: str, encoding: str = "utf-8"):
    """Write sequences to FASTA file."""
    with open(path, "w", encoding=encoding) as f:
        for name, seq in seqs.items():
            f.write(f">{name}\n{seq}\n")


def write_newick(tree_str: str, path: str, encoding: str = "utf-8"):
    """Write Newick string to file."""
    with open(path, "w", encoding=encoding) as f:
        if not tree_str.endswith(";"):
            tree_str += ";"
        f.write(tree_str + "\n")


def check_leaf_completeness(tree_nwk: str, expected_leaves: set) -> bool:
    """Check that a Newick tree contains all expected leaf names."""
    import re
    leaves = set(re.findall(r"[A-Za-z0-9_\.]+(?=:)", tree_nwk))
    # Remove branch length numbers that might look like leaves
    leaves = {l for l in leaves if not l.replace(".", "").isdigit()}
    return expected_leaves == leaves
