"""
fastme_backend.py — FastME (Balanced Minimum Evolution) backend for Fusang.

Replaces O(n³) NJ with O(n² log n) FastME.
Self-contained: no dependency on fusang_v2 (avoid numba import chain).
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import numpy as np


# ── Lightweight TreeNode (avoid importing fusang_v2 → numba) ──────────

class _TreeNode:
    """Lightweight Newick tree node, compatible with fusang_v2.TreeNode."""
    def __init__(self, name: Optional[str] = None, dist: float = 0.0):
        self.name = name
        self.dist = dist          # fusang_v2.TreeNode uses .dist
        self.children: List["_TreeNode"] = []
        self.up: Optional["_TreeNode"] = None

    def get_leaves(self) -> List["_TreeNode"]:
        """Return all leaf nodes."""
        if not self.children:
            return [self]
        result = []
        for c in self.children:
            result.extend(c.get_leaves())
        return result

    def write(self, format=1):
        """Output NEWICK string with trailing semicolon."""
        return self._to_newick(format=format) + ";"

    def _to_newick(self, format=1):
        """Output NEWICK string for this subtree.

        format=0: full NEWICK with internal node names
        format=1: standard NEWICK (leaf names only, internal nodes unnamed)
        """
        if not self.children:
            # Leaf node: output name:dist
            label = self.name or ""
            if format >= 1 and self.dist >= 0:
                label += ":" + str(round(self.dist, 6))
            return label
        # Internal node: recurse into children
        child_strs = []
        for ch in self.children:
            cs = ch._to_newick(format=format)
            child_strs.append(cs)
        result = "(" + ",".join(child_strs) + ")"
        # Append branch length for internal nodes (format >= 1)
        if format >= 1 and self.dist > 0:
            result += ":" + str(round(self.dist, 6))
        return result

    def __repr__(self):
        return "TreeNode(name=%s, dist=%.4f, n_children=%d)" % (
            self.name, self.dist, len(self.children))

    def detach(self):
        """Remove this node from its parent."""
        if self.up:
            self.up.children.remove(self)
            self.up = None


# ── WSL Detection (environment variable override supported) ──────────────

_DEFUALT_WSL_DISTROS = ["Ubuntu-24.04", "Ubuntu-22.04", "Ubuntu"]


def _find_wsl_distro() -> Optional[str]:
    """
    Find a usable WSL distribution name.
    Priority:
      1. Environment variable FASTME_WSL_DISTRO (comma-separated, e.g. "Ubuntu-24.04,Ubuntu-22.04")
      2. Try known distro names (fast, reliable)
      3. Fallback: parse `wsl -l -q` (slower, encoding-sensitive)
    Returns distro name like 'Ubuntu-24.04', or None.
    """
    if os.name != "nt":
        return None

    # 1. Environment variable override (most reliable)
    env_distros = os.environ.get("FASTME_WSL_DISTRO", "").strip()
    if env_distros:
        for d in env_distros.split(","):
            d = d.strip()
            if not d:
                continue
            try:
                r = subprocess.run(
                    ["wsl", "-d", d, "--", "echo", "WSL_OK"],
                    capture_output=True,
                    timeout=30,   # WSL cold start can be slow
                )
                if r.returncode == 0 and b"WSL_OK" in r.stdout:
                    return d
            except Exception:
                pass

    # 2. Try known distro names
    for distro in _DEFUALT_WSL_DISTROS:
        try:
            r = subprocess.run(
                ["wsl", "-d", distro, "--", "echo", "WSL_OK"],
                capture_output=True,
                timeout=30,
            )
            if r.returncode == 0 and b"WSL_OK" in r.stdout:
                return distro
        except Exception:
            pass

    # 3. Fallback: parse `wsl -l -q`
    try:
        r = subprocess.run(
            ["wsl", "-l", "-q"],
            capture_output=True,
            timeout=10,
        )
        if r.returncode != 0 or not r.stdout:
            return None
        # WSL outputs UTF-16 LE on Windows
        text = r.stdout.decode("utf-16-le", errors="replace")
        for line in text.splitlines():
            d = line.strip().strip("\x00")
            if not d or d.startswith("\ufffd") or "docker" in d.lower():
                continue
            # Validate this distro works
            try:
                vr = subprocess.run(
                    ["wsl", "-d", d, "--", "echo", "OK"],
                    capture_output=True,
                    timeout=30,
                )
                if vr.returncode == 0:
                    return d
            except Exception:
                pass
    except Exception:
        pass

    return None


# ── FastME Binary Discovery ────────────────────────────────────────────

def _find_fastme_binary() -> Tuple[Optional[List[str]], bool]:
    """
    Find FastME binary. Returns (command: list[str], is_native_windows: bool).
    Priority:
      0. Running inside Linux/WSL  → use 'fastme' from PATH
      1. fastme_bin/fastme.exe      (Windows native, fastest)
      2. WSL-installed /usr/local/bin/fastme
      3. fastme_bin/fastme_linux    (project-bundled, WSL fallback)
    """
    # 0. Already inside Linux/WSL?
    if os.name == "posix":
        try:
            r = subprocess.run(["which", "fastme"], capture_output=True, timeout=5)
            if r.returncode == 0 and r.stdout.strip():
                return (["fastme"], False)
        except Exception:
            pass

    script_dir = Path(__file__).parent

    # 1. Windows native .exe (preferred — no cross-system call overhead)
    fastme_exe = script_dir / "fastme_bin" / "fastme.exe"
    if fastme_exe.exists():
        return ([str(fastme_exe)], True)

    # 2. WSL-installed binary
    wsl_distro = _find_wsl_distro()
    if wsl_distro:
        try:
            r = subprocess.run(
                ["wsl", "-d", wsl_distro, "--",
                 "test", "-x", "/usr/local/bin/fastme"],
                capture_output=True,
                timeout=15,
            )
            if r.returncode == 0:
                return (["wsl", "-d", wsl_distro, "--", "/usr/local/bin/fastme"], False)
        except Exception:
            pass

    # 3. Project-bundled WSL binary
    fastme_linux = script_dir / "fastme_bin" / "fastme_linux"
    if wsl_distro and fastme_linux.exists() and os.access(str(fastme_linux), os.R_OK):
        wsl_path = _win_to_wsl_path(str(fastme_linux))
        return (["wsl", "-d", wsl_distro, "--", wsl_path], False)

    return (None, False)


# ── Windows → WSL Path Conversion ──────────────────────────────────────

def _win_to_wsl_path(win_path: str) -> str:
    """Convert Windows absolute path to WSL path. D:\\foo\\bar.fa → /mnt/d/foo/bar.fa"""
    if not win_path or os.name != "nt":
        return win_path
    p = win_path.replace("\\", "/")
    m = re.match(r"^([A-Za-z]):/(.*)", p)
    if m:
        return f"/mnt/{m.group(1).lower()}/{m.group(2)}"
    return p


# ── PHYLIP Distance Matrix Writer ─────────────────────────────────────

def _write_phylip_matrix(D: np.ndarray, names: Sequence[str], out_path: str) -> None:
    """Write pairwise distance matrix in PHYLIP format for FastME."""
    n = D.shape[0]
    assert D.shape == (n, n)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(f"{n}\n")
        for i in range(n):
            row = " ".join(f"{D[i,j]:.4f}" for j in range(n))
            f.write(f"{names[i]:<10s} {row}\n")


# ── Newick Parser ───────────────────────────────────────────────────────

def parse_newick_to_treenode(nwk: str) -> _TreeNode:
    """Parse a Newick string into a _TreeNode tree."""
    i = 0
    nwk = nwk.strip().rstrip(";")

    def _parse_node():
        nonlocal i
        if i >= len(nwk):
            return None
        if nwk[i] == "(":
            node = _TreeNode()
            i += 1  # skip '('
            children = []
            while i < len(nwk) and nwk[i] != ")":
                child = _parse_node()
                if child:
                    children.append(child)
                if i < len(nwk) and nwk[i] == ",":
                    i += 1
            if i < len(nwk) and nwk[i] == ")":
                i += 1
            node.children = children
            for ch in children:
                ch.up = node
            _read_branch_length(node)
            return node
        else:
            name_chars = []
            while i < len(nwk) and nwk[i] not in (",", ")", ":"):
                name_chars.append(nwk[i])
                i += 1
            name = "".join(name_chars).strip()
            node = _TreeNode(name=name or None)
            _read_branch_length(node)
            return node

    def _read_branch_length(node: _TreeNode):
        nonlocal i
        # Skip bootstrap/value before colon (e.g., "0.95:0.123")
        if i < len(nwk) and nwk[i] not in (",", ")", ":"):
            while i < len(nwk) and nwk[i] != ":":
                i += 1
        if i < len(nwk) and nwk[i] == ":":
            i += 1  # skip ':'
            num_chars = []
            while i < len(nwk) and (nwk[i].isdigit() or nwk[i] in (".", "e", "-", "+")):
                num_chars.append(nwk[i])
                i += 1
            if num_chars:
                try:
                    node.dist = float("".join(num_chars))
                except ValueError:
                    pass

    root = _parse_node()
    if root is None:
        raise ValueError(f"Failed to parse Newick: {nwk[:80]}")
    return root


# ── Main API: build_tree_fastme() ────────────────────────────────────

def build_tree_fastme(distance_matrix: np.ndarray,
                      taxon_names: List[str],
                      n_threads: int = 4,
                      method: str = "I",
                      nni_type: str = "B",
                      verbose: int = 0) -> _TreeNode:
    """
    Build a phylogenetic tree using FastME.
    Returns _TreeNode (compatible with fusang_v2.TreeNode).
    """
    n = len(taxon_names)
    if n < 4:
        # Fallback to NJ for small inputs (n=3 for backbone trees)
        import sys as _sys
        _sys.path.insert(0, str(Path(__file__).parent))
        from fusang_v2 import neighbor_joining
        return neighbor_joining(distance_matrix, taxon_names)

    work_dir = tempfile.mkdtemp(prefix="fastme_")
    mat_path = os.path.join(work_dir, "dist_matrix.phy")
    tree_path = os.path.join(work_dir, "tree.nwk")

    # Debug log file (bypass Windows console encoding issues)
    _log_path = os.path.join(work_dir, "debug.log")
    def _log(msg: str):
        try:
            with open(_log_path, "a", encoding="utf-8", errors="replace") as f:
                f.write(str(msg) + "\n")
                f.flush()
        except Exception:
            pass

    try:
        _log(f"build_tree_fastme started, n={n}, work_dir={work_dir}")
        _log(f"Python: {sys.executable}, os.name={os.name}")

        # Write PHYLIP distance matrix
        _log("Writing PHYLIP matrix...")
        _write_phylip_matrix(distance_matrix, taxon_names, mat_path)
        _log(f"PHYLIP written: {mat_path}")

        # Find FastME binary
        _log("Finding FastME binary...")
        fastme_cmd, is_native = _find_fastme_binary()
        _log(f"FastME cmd: {fastme_cmd}, native: {is_native}")
        if fastme_cmd is None:
            raise FileNotFoundError(
                "FastME binary not found.\n"
                "Install via: https://gitlab.lirmm.fr/atgc/FastME/\n"
                "Or set env FASTME_WSL_DISTRO=Ubuntu-24.04"
            )

        # If running via WSL, convert paths; native exe uses Windows paths directly
        is_wsl = not is_native
        if is_wsl:
            mat_input = _win_to_wsl_path(mat_path)
            tree_output = _win_to_wsl_path(tree_path)
            _log(f"WSL mat_input: {mat_input}")
            _log(f"WSL tree_output: {tree_output}")
        else:
            mat_input = mat_path
            tree_output = tree_path

        cmd = fastme_cmd + [
            "-i", mat_input,
            "-o", tree_output,
            "-m", method,
            "-n", nni_type,
            "-T", str(n_threads),
            "-v", str(verbose),
        ]
        _log(f"Full cmd: {cmd}")

        # Run FastME — bytes mode (text=False) to avoid Windows encoding crash
        _log("Running FastME subprocess...")
        result = subprocess.run(cmd, capture_output=True, timeout=600)
        _log(f"FastME returncode: {result.returncode}")
        stdout_str = result.stdout.decode("utf-8", errors="replace") if result.stdout else "(empty)"
        stderr_str = result.stderr.decode("utf-8", errors="replace") if result.stderr else "(empty)"
        _log(f"FastME stdout (first 500): {stdout_str[:500]}")
        _log(f"FastME stderr (first 500): {stderr_str[:500]}")

        if result.returncode != 0:
            raise RuntimeError(
                f"FastME failed (rc={result.returncode})\n"
                f"STDOUT: {stdout_str}\nSTDERR: {stderr_str}"
            )

        # Read and parse output
        _log(f"Checking tree_path exists: {os.path.exists(tree_path)}")
        if not os.path.exists(tree_path):
            if is_wsl:
                # Read from WSL directly
                _log("Tree not found at Windows path, reading from WSL...")
                wsl_distro = fastme_cmd[2]  # e.g. 'Ubuntu-24.04'
                r = subprocess.run(
                    ["wsl", "-d", wsl_distro, "--", "cat", tree_output],
                    capture_output=True,
                    timeout=10,
                )
                if r.returncode == 0 and r.stdout:
                    nwk_str = r.stdout.decode("utf-8", errors="replace").strip()
                    _log(f"Read tree from WSL, length={len(nwk_str)}")
                else:
                    _log(f"WSL cat failed, rc={r.returncode}")
                    raise FileNotFoundError(
                        f"FastME output not found: {tree_path}\n"
                        f"(WSL wrote to: {tree_output})"
                    )
            else:
                _log(f"work_dir contents: {os.listdir(work_dir)}")
                raise FileNotFoundError(f"FastME output not found: {tree_path}")
        else:
            with open(tree_path, "r", encoding="utf-8", errors="replace") as f:
                nwk_str = f.read().strip()
            _log(f"Tree read from file, length={len(nwk_str)}")

        lines = [l.strip() for l in nwk_str.split("\n") if l.strip()]
        if not lines:
            raise ValueError(f"FastME output empty: {tree_path}")

        _log("Parsing Newick...")
        tree = parse_newick_to_treenode(lines[-1])
        _log("Parse complete.")
        return tree

    except Exception as e:
        _log(f"EXCEPTION: {type(e).__name__}: {e}")
        import traceback
        _log(traceback.format_exc())
        # Print log path to stderr (bypass console encoding)
        try:
            sys.stderr.write(f"[FastME DEBUG log: {_log_path}]\n")
            sys.stderr.flush()
        except Exception:
            pass
        raise


# ── Standalone Test ───────────────────────────────────────────────────

if __name__ == "__main__":
    # ALL output goes to a file to avoid Windows console encoding crash
    out_file = "C:/Users/admin/AppData/Local/Temp/fastme_test_output.txt"
    log_file = "C:/Users/admin/AppData/Local/Temp/fastme_test_debug.txt"

    def out(msg: str):
        with open(out_file, "a", encoding="utf-8", errors="replace") as f:
            f.write(str(msg) + "\n")
            f.flush()

    out("=" * 60)
    out("Testing FastME backend...")
    out(f"Python: {sys.executable}")
    out(f"os.name: {os.name}")

    n = 5
    names = [f"T{i+1}" for i in range(n)]
    D = np.array([
        [0.0, 0.1, 0.2, 0.3, 0.4],
        [0.1, 0.0, 0.15, 0.25, 0.35],
        [0.2, 0.15, 0.0, 0.2, 0.3],
        [0.3, 0.25, 0.2, 0.0, 0.1],
        [0.4, 0.35, 0.3, 0.1, 0.0],
    ], dtype=np.float64)

    out(f"Taxa: {n}")
    out(f"Distance matrix shape: {D.shape}")

    # Check FastME availability
    distro = _find_wsl_distro()
    out(f"WSL distro: {distro}")

    cmd, _ = _find_fastme_binary()
    out(f"FastME command: {cmd}")

    if cmd is None:
        out("ERROR: FastME binary not found!")
        out("Please install FastME or set FASTME_WSL_DISTRO env var.")
        # Also write to stderr for visibility
        sys.stderr.write("FastME binary not found. Check " + out_file + "\n")
        sys.exit(1)

    # Build tree
    out("Calling build_tree_fastme...")
    tree = build_tree_fastme(D, names, n_threads=2, verbose=1)
    out("Tree built successfully!")

    def count_leaves(node: _TreeNode) -> int:
        if not node.children:
            return 1
        return sum(count_leaves(c) for c in node.children)

    n_leaves = count_leaves(tree)
    out(f"Number of leaves: {n_leaves} (expected: {n})")

    nwk = str(tree)
    out(f"Newick (first 200 chars): {nwk[:200]}")
    out("Done!")

    # Also print to stderr for visibility (short message, ASCII-only)
    sys.stderr.write(f"Test complete. Output: {out_file}\n")
    sys.stderr.write(f"Debug log: {log_file}\n")
