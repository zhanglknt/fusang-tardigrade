#!/usr/bin/env python3
"""Competitor comparison on INDELible n=200 indel data, scored with the STANDARD
normalised Robinson-Foulds distance (``nrf_standard.py``).

This is a drop-in replacement for ``run_competitors_indelible.py``.  The pipelines
(distance matrices, NJ implementation, Fusang CLI configuration, inputs) are
byte-for-byte identical to the ones that produced
``competitors_n200_indel_130seeds_fusangbench.csv``; the ONLY change is the RF
convention:

  old  (calc_nrf_simple.py)  per internal node only the SMALLER child clade is
                             kept, so bifurcating trees yield ~64 of the 197
                             non-trivial splits -> RF is systematically
                             underestimated and rooting-dependent.
  new  (nrf_standard.py)     root-independent standard RF: every non-trivial
                             split canonicalised to the side NOT containing the
                             lexicographically smallest leaf,
                             RF = |S1-S2| + |S2-S1|, nRF = RF / (2(n-3)),
                             n = 200 -> denominator 394.

Every inferred tree is written to disk (<data dir>/trees_stdrf/seed{N}_{method}.nwk)
so that any future change of RF convention can be recomputed without rerunning
the distance/NJ pipelines.

Usage (WSL, fusang_bench):
  python run_competitors_stdrf.py --seeds 100-229 --workers 2 \
      --out-prefix competitors_n200_indel_130seeds_stdrf
"""

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.setrecursionlimit(10000)

from af_competitor_methods import (  # noqa: E402
    cophylog_distance_matrix,
    kmer_cosine_distance_matrix,
    read_fasta,
)
from fusang_v4_dahp_v1 import build_nj  # noqa: E402
import nrf_standard as nrfstd  # noqa: E402

REPO_PY = os.environ.get("FUSANG_PYTHON", "python")
FUSANG_PY = os.path.join(HERE, "fusang_v2.py")

INPROC_METHODS = ["cophylog_k19", "cosine_k5", "cosine_k7", "cosine_k9",
                  "multik_ensemble", "fusang_spaced_k5_nj"]
KMACS_METHODS = ["kmacs_k3", "kmacs_k5", "kmacs_k10"]
ALL_METHODS = ["cophylog_k19", "cosine_k5", "cosine_k7", "cosine_k9",
               "multik_ensemble", "fusang_spaced_k5_nj",
               "kmacs_k3", "kmacs_k5", "kmacs_k10", "fusang_spaced_k5"]

GAP2_MASK = "1110011"

METHOD_LABELS = {
    "fusang_spaced_k5": "Fusang spaced k=5 gap2 (CLI: NJ + BME NNI)",
    "fusang_spaced_k5_nj": "Fusang spaced k=5 gap2 (plain NJ, no NNI)",
    "cophylog_k19": "Co-phylog (k=19)",
    "kmacs_k3": "kmacs (k=3)",
    "kmacs_k5": "kmacs (k=5)",
    "kmacs_k10": "kmacs (k=10)",
    "cosine_k5": "K-mer cosine (k=5)",
    "cosine_k7": "K-mer cosine (k=7)",
    "cosine_k9": "K-mer cosine (k=9)",
    "multik_ensemble": "Multi-k ensemble (k=5,7,9)",
}

TREES_SUBDIR = "trees_stdrf"


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def read_text(path):
    with open(path, encoding="utf-8", errors="replace") as fh:
        return fh.read().strip()


def _splits_of(nwk):
    """(frozenset splits, n_leaves) with std RF convention."""
    tree = nrfstd.parse_newick(nwk)
    lv = nrfstd.leaves(tree)
    return nrfstd.splits(tree), len(set(sorted(lv)))


def nrf_between(sa, sb, n):
    return (len(sa - sb) + len(sb - sa)) / (2 * (n - 3))


def parse_phylip_dmat(path):
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = [ln for ln in fh.read().split("\n") if ln.strip()]
    n = int(lines[0].split()[0])
    names, rows = [], []
    for ln in lines[1:]:
        parts = ln.split()
        names.append(parts[0].strip())
        rows.append([float(x) for x in parts[1:]])
    if len(names) != n or any(len(r) != n for r in rows):
        raise ValueError("malformed matrix %s" % path)
    M = np.array(rows, dtype=np.float64)
    M = (M + M.T) / 2.0
    np.fill_diagonal(M, 0.0)
    return names, M


# --------------------------------------------------------------------------
# per-seed worker
# --------------------------------------------------------------------------
def process_seed(task):
    seed, data_dir, kmacs_dir, trees_dir, run_fusang = task
    tag = "indel"
    fasta = os.path.join(data_dir, "seed%d_%s.fasta" % (seed, tag))
    ft2_path = os.path.join(data_dir, "seed%d_%s_ft2.nwk" % (seed, tag))
    true_path = os.path.join(data_dir, "seed%d_%s_true.nwk" % (seed, tag))
    fusang_path = os.path.join(data_dir, "seed%d_%s_fusang.nwk" % (seed, tag))

    row = {"seed": seed}
    nsplit = {}
    t0 = time.time()
    for m in ALL_METHODS:
        row[m + "_ft2"] = ""
        row[m + "_true"] = ""

    if not (os.path.exists(fasta) and os.path.exists(ft2_path) and os.path.exists(true_path)):
        row["error"] = "missing input"
        return row

    try:
        seqs_dict = read_fasta(fasta)
        names = list(seqs_dict.keys())
        seqs = [seqs_dict[n] for n in names]
        sa_true, n_true = _splits_of(read_text(true_path))
        sa_ft2, n_ft2 = _splits_of(read_text(ft2_path))
    except Exception as exc:
        row["error"] = "load failed: %s" % exc
        return row
    nsplit["TRUE"] = len(sa_true)
    nsplit["FT2"] = len(sa_ft2)
    n = n_true

    def _save_tree(method, nwk):
        p = os.path.join(trees_dir, "seed%d_%s.nwk" % (seed, method))
        with open(p, "w", encoding="utf-8") as fh:
            fh.write(nwk if nwk.endswith(";") else nwk + "\n")
        return p

    def _score(method, nwk):
        """save tree, compute both nRFs, record split count."""
        _save_tree(method, nwk)
        sa, k = _splits_of(nwk)
        nsplit[method] = len(sa)
        if k != n:
            row["error"] = row.get("error", "") + " %s:ntaxa %d!=%d" % (method, k, n)
        row[method + "_ft2"] = nrf_between(sa, sa_ft2, n)
        row[method + "_true"] = nrf_between(sa, sa_true, n)

    # ---------------- distance matrices + NJ ----------------
    mats = {}
    try:
        mats["cophylog_k19"] = cophylog_distance_matrix(seqs, names, halfctx=9)
    except Exception as exc:
        row["error"] = row.get("error", "") + " cophylog:%s" % exc
    for k in (5, 7, 9):
        try:
            mats["cosine_k%d" % k] = kmer_cosine_distance_matrix(seqs, names, k=k)
        except Exception as exc:
            row["error"] = row.get("error", "") + " cosine_k%d:%s" % (k, exc)
    if all(("cosine_k%d" % k) in mats for k in (5, 7, 9)):
        mats["multik_ensemble"] = (
            np.asarray(mats["cosine_k5"], dtype=np.float64)
            + np.asarray(mats["cosine_k7"], dtype=np.float64)
            + np.asarray(mats["cosine_k9"], dtype=np.float64)
        ) / 3.0
    try:
        mats["fusang_spaced_k5_nj"] = kmer_cosine_distance_matrix(
            seqs, names, k=5, gap_pattern=GAP2_MASK)
    except Exception as exc:
        row["error"] = row.get("error", "") + " fusang_gap2_matrix:%s" % exc

    for m in INPROC_METHODS:
        if m not in mats:
            continue
        try:
            nwk = build_nj(np.asarray(mats[m], dtype=np.float64), names)
        except Exception as exc:
            row["error"] = row.get("error", "") + " nj_%s:%s" % (m, exc)
            continue
        if not nwk:
            continue
        _score(m, nwk)

    # ---------------- kmacs ----------------
    for m in KMACS_METHODS:
        k = int(m.split("_k")[1])
        dmat = os.path.join(kmacs_dir, "seed%d_kmacs_k%d.dmat" % (seed, k))
        if not os.path.exists(dmat):
            row["error"] = row.get("error", "") + " kmacs_k%d:missing" % k
            continue
        try:
            knames, M = parse_phylip_dmat(dmat)
            if set(knames) != set(names):
                row["error"] = row.get("error", "") + " kmacs_k%d:name-mismatch" % k
                continue
            nwk = build_nj(M, knames)
            if not nwk:
                continue
            _score(m, nwk)
        except Exception as exc:
            row["error"] = row.get("error", "") + " kmacs_k%d:%s" % (k, exc)

    # ---------------- Fusang spaced k=5 gap2 (CLI) ----------------
    if run_fusang or not os.path.exists(fusang_path):
        try:
            out_dir = os.path.join(data_dir, "competitors_out")
            os.makedirs(out_dir, exist_ok=True)
            tmp_out = os.path.join(out_dir, "seed%d_fusang_cli.nwk" % seed)
            cmd = [REPO_PY, FUSANG_PY, "--input", fasta, "--output", tmp_out,
                   "--tree_method", "nj", "--kmer_k", "5", "--kmer_gap", "gap2"]
            subprocess.run(cmd, capture_output=True, text=True, cwd=HERE, timeout=1800)
            if os.path.exists(tmp_out) and os.path.getsize(tmp_out) > 0:
                fusang_path = tmp_out
        except Exception as exc:
            row["error"] = row.get("error", "") + " fusang_cli:%s" % exc
    if os.path.exists(fusang_path) and os.path.getsize(fusang_path) > 0:
        try:
            _score("fusang_spaced_k5", read_text(fusang_path))
        except Exception as exc:
            row["error"] = row.get("error", "") + " fusang_nrf:%s" % exc
    else:
        row["error"] = row.get("error", "") + " fusang_tree:missing"

    row["nsplits"] = nsplit
    row["seconds"] = round(time.time() - t0, 1)
    return row


def parse_seed_spec(spec):
    out = []
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            lo, hi = part.split("-")
            out.extend(range(int(lo), int(hi) + 1))
        else:
            out.append(int(part))
    return sorted(set(out))


# --------------------------------------------------------------------------
# summary statistics
# --------------------------------------------------------------------------
def summarize(rows, frame):
    """frame in {'ft2','true'}; per-method stats + paired comparison vs Fusang CLI."""
    from scipy import stats

    def vals(m):
        out = []
        for r in rows:
            v = r.get(m + "_" + frame)
            if v is None or v == "":
                continue
            try:
                f = float(v)
            except (TypeError, ValueError):
                continue
            if np.isnan(f):
                continue
            out.append((r["seed"], f))
        return out

    base = dict(vals("fusang_spaced_k5"))
    summary = {}
    for m in ALL_METHODS:
        v = vals(m)
        seeds = [s for s, _ in v]
        arr = np.array([x for _, x in v], dtype=np.float64)
        entry = {
            "label": METHOD_LABELS[m],
            "n": len(arr),
            "mean": float(np.mean(arr)) if len(arr) else None,
            "sd": float(np.std(arr, ddof=1)) if len(arr) > 1 else None,
            "median": float(np.median(arr)) if len(arr) else None,
            "min": float(np.min(arr)) if len(arr) else None,
            "max": float(np.max(arr)) if len(arr) else None,
            "vs_fusang": None,
        }
        if len(arr) == 0:
            entry["status"] = "NO DATA"
            summary[m] = entry
            continue
        if m == "fusang_spaced_k5":
            summary[m] = entry
            continue
        pairs = [(base[s], x) for s, x in zip(seeds, arr) if s in base]
        if len(pairs) >= 5:
            a = np.array([p[0] for p in pairs], dtype=np.float64)   # Fusang CLI
            b = np.array([p[1] for p in pairs], dtype=np.float64)   # method
            diff = b - a
            wins = int(np.sum(b < a))
            losses = int(np.sum(b > a))
            ties = int(np.sum(b == a))
            try:
                _, p = stats.wilcoxon(a, b, zero_method="zsplit")
                p = float(p)
            except Exception:
                p = None
            sd_diff = float(np.std(diff, ddof=1))
            d = float(np.mean(diff) / sd_diff) if sd_diff > 0 else 0.0
            entry["vs_fusang"] = {
                "n_paired": len(pairs),
                "fusang_mean": float(np.mean(a)),
                "method_mean": float(np.mean(b)),
                "delta_method_minus_fusang": float(np.mean(diff)),
                "wins_vs_fusang": wins,
                "losses_vs_fusang": losses,
                "ties": ties,
                "wilcoxon_p_zsplit": p,
                "cohens_d_paired": d,
            }
        summary[m] = entry
    return summary


def rescore_rows(rows, data_dir, trees_dir):
    """Recompute every nRF in the cached rows from the saved Newick files.

    Uses the CURRENT nrf_standard implementation (import-time), so a fix to the
    split extraction can be applied without rerunning any distance/NJ pipeline.
    """
    fixed = 0
    for r in rows:
        seed = r["seed"]
        try:
            sa_true, n = _splits_of(read_text(
                os.path.join(data_dir, "seed%d_indel_true.nwk" % seed)))
            sa_ft2, _ = _splits_of(read_text(
                os.path.join(data_dir, "seed%d_indel_ft2.nwk" % seed)))
        except Exception as exc:
            r["error"] = (r.get("error", "") + " rescore_ref:%s" % exc).strip()
            continue
        nsplit = {"TRUE": len(sa_true), "FT2": len(sa_ft2)}
        for m in ALL_METHODS:
            p = os.path.join(trees_dir, "seed%d_%s.nwk" % (seed, m))
            if not os.path.exists(p):
                continue
            try:
                sa, k = _splits_of(read_text(p))
            except Exception as exc:
                r["error"] = (r.get("error", "") + " rescore_%s:%s" % (m, exc)).strip()
                continue
            nsplit[m] = len(sa)
            if k != n:
                r["error"] = (r.get("error", "") + " %s:ntaxa %d!=%d" % (m, k, n)).strip()
            r[m + "_ft2"] = nrf_between(sa, sa_ft2, n)
            r[m + "_true"] = nrf_between(sa, sa_true, n)
            fixed += 1
        r["nsplits"] = nsplit
    print("  rescored %d (seed, method) entries from %s" % (fixed, trees_dir))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seeds", default="100-229")
    ap.add_argument("--dir", default="indelible_n200_indel")
    ap.add_argument("--kmacs-dir", default=None)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--out-prefix", default=None)
    ap.add_argument("--cache-dir", default=None)
    ap.add_argument("--trees-dir", default=None,
                    help="directory for per-seed inferred trees (default: <dir>/trees_stdrf)")
    ap.add_argument("--run-fusang", action="store_true",
                    help="always (re)run the Fusang CLI instead of reusing "
                         "seed{N}_indel_fusang.nwk")
    ap.add_argument("--rescore", action="store_true",
                    help="do not run any pipeline: reload the per-seed checkpoints, "
                         "recompute every nRF from the Newick files in --trees-dir "
                         "with the CURRENT nrf_standard implementation, and rewrite "
                         "the CSV/JSON outputs (used when only the RF convention "
                         "changed, e.g. after a bugfix)")
    args = ap.parse_args()

    data_dir = args.dir if os.path.isabs(args.dir) else os.path.join(HERE, args.dir)
    kmacs_dir = args.kmacs_dir or os.path.join(data_dir, "kmacs_out_indelible")
    trees_dir = args.trees_dir or os.path.join(data_dir, TREES_SUBDIR)
    os.makedirs(trees_dir, exist_ok=True)
    seeds = parse_seed_spec(args.seeds)
    prefix = args.out_prefix or "competitors_n200_indel_%dseeds_stdrf" % len(seeds)

    print("=" * 78)
    print("  Competitor comparison on INDELible data - STANDARD unrooted RF")
    print("  seeds: %d (%d..%d)   data: %s" % (len(seeds), seeds[0], seeds[-1], data_dir))
    print("  kmacs matrices : %s" % kmacs_dir)
    print("  trees saved to : %s" % trees_dir)
    print("  nRF convention : standard root-independent RF, nRF = RF/(2(n-3)), n=200")
    print("  methods        : %s" % ", ".join(ALL_METHODS))
    print("=" * 78, flush=True)

    cache_dir = args.cache_dir or os.path.join(data_dir, ".cache_" + prefix)
    os.makedirs(cache_dir, exist_ok=True)

    def _cache_path(seed):
        return os.path.join(cache_dir, "seed%d.json" % seed)

    def _save_row(row):
        tmp = _cache_path(row["seed"]) + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(row, fh)
        os.replace(tmp, _cache_path(row["seed"]))

    rows = []
    for fn in sorted(os.listdir(cache_dir)):
        if fn.startswith("seed") and fn.endswith(".json"):
            try:
                with open(os.path.join(cache_dir, fn)) as fh:
                    cached_row = json.load(fh)
                if cached_row.get("seed") in seeds:
                    rows.append(cached_row)
            except Exception:
                pass
    done_seeds = {r["seed"] for r in rows}
    if done_seeds:
        print("  resume: %d seeds already cached in %s" % (len(done_seeds), cache_dir))

    if args.rescore:
        rows = rescore_rows(rows, data_dir, trees_dir)
        for r in rows:
            _save_row(r)
        missing = sorted(set(seeds) - done_seeds)
        if missing:
            print("  WARNING: %d requested seeds have no checkpoint: %s%s"
                  % (len(missing), missing[:10], "..." if len(missing) > 10 else ""))
        tasks = []
    else:
        tasks = [(s, data_dir, kmacs_dir, trees_dir, args.run_fusang)
                 for s in seeds if s not in done_seeds]
    n_total = len(seeds)
    t0 = time.time()
    if args.workers > 1:
        with Pool(processes=args.workers) as pool:
            for i, row in enumerate(pool.imap_unordered(process_seed, tasks)):
                rows.append(row)
                _save_row(row)
                print("[%d/%d] seed %d  %s  (%.0fs elapsed)"
                      % (i + 1 + len(done_seeds), n_total, row["seed"],
                         " ".join("%s=%.4f" % (m.split('_')[0][:6], row[m + "_ft2"])
                                  for m in ALL_METHODS
                                  if isinstance(row.get(m + "_ft2"), float)),
                         time.time() - t0), flush=True)
    else:
        for i, t in enumerate(tasks):
            row = process_seed(t)
            rows.append(row)
            _save_row(row)
            print("[%d/%d] seed %d done" % (i + 1 + len(done_seeds), n_total, row["seed"]),
                  flush=True)

    rows.sort(key=lambda r: r["seed"])

    # ---------------- CSV ----------------
    cols = ["seed"]
    for m in ALL_METHODS:
        cols += [m + "_ft2", m + "_true"]
    extra = [c for c in ("seconds", "error") if any(c in r for r in rows)]
    cols += extra
    csv_path = os.path.join(data_dir, prefix + ".csv")
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print("\nwrote %s" % csv_path)

    # ---------------- trees manifest ----------------
    manifest = {}
    for m in ALL_METHODS:
        got = sum(1 for s in seeds
                  if os.path.exists(os.path.join(trees_dir, "seed%d_%s.nwk" % (s, m))))
        manifest[m] = got
    with open(os.path.join(trees_dir, "_manifest.json"), "w") as fh:
        json.dump({"n_seeds": len(seeds), "trees_per_method": manifest}, fh, indent=2)

    # ---------------- split-resolution diagnostics ----------------
    diag = {}
    for key in ALL_METHODS + ["TRUE", "FT2"]:
        xs = [r.get("nsplits", {}).get(key) for r in rows]
        xs = [x for x in xs if isinstance(x, int)]
        if xs:
            diag[key] = {"n_trees": len(xs), "min_splits": min(xs), "max_splits": max(xs),
                         "mean_splits": round(sum(xs) / len(xs), 2),
                         "fully_resolved_197": sum(1 for x in xs if x == 197)}
    print("\n  split counts (fully resolved n=200 unrooted tree = 197):")
    for k, v in diag.items():
        print("    %-22s n=%3d  min=%3d max=%3d mean=%6.2f  resolved=%d"
              % (k, v["n_trees"], v["min_splits"], v["max_splits"], v["mean_splits"],
                 v["fully_resolved_197"]))

    # ---------------- JSON ----------------
    summary = {
        "data": {
            "simulator": "INDELible v1.03 (GTR+G, indel rate 0.02, root length L=500)",
            "directory": os.path.relpath(data_dir, HERE),
            "n_taxa": 200,
            "seeds": [seeds[0], seeds[-1]],
            "n_seeds_requested": len(seeds),
            "n_seeds_scored": sum(1 for r in rows
                                  if isinstance(r.get("fusang_spaced_k5_ft2"), float)),
            "nrf_convention": ("standard root-independent (unrooted) RF via nrf_standard.py: "
                               "non-trivial splits canonicalised to the side lacking the "
                               "lexicographically smallest leaf; "
                               "RF=|S1-S2|+|S2-S1|; nRF=RF/(2(n-3)), n=200 -> /394"),
            "previous_convention": ("calc_nrf_simple.py kept only the smaller child clade per "
                                    "internal node (~64/197 splits captured) and is therefore "
                                    "biased low and rooting-dependent"),
            "nj_implementation": ("BioPython DistanceTreeConstructor().nj "
                                  "(fusang_v4_dahp_v1.build_nj) for every distance-matrix "
                                  "method; fusang_spaced_k5 additionally uses the official CLI "
                                  "(fusang_v2.py --tree_method nj --kmer_k 5 --kmer_gap gap2), "
                                  "which applies a BME/NNI refinement pass that competitor "
                                  "methods do not receive"),
            "fusang_config": "k=5, gap2 (position mask 1110011 = positions 0,1,2,5,6)",
            "trees_dir": os.path.relpath(trees_dir, HERE),
            "wilcoxon": "scipy.stats.wilcoxon(a, b, zero_method='zsplit')",
        },
        "split_resolution_diagnostics": diag,
        "vs_ft2": summarize(rows, "ft2"),
        "vs_true": summarize(rows, "true"),
        "rows": rows,
    }
    json_path = os.path.join(data_dir, prefix + "_summary.json")
    with open(json_path, "w") as fh:
        json.dump(summary, fh, indent=2)
    print("wrote %s" % json_path)

    # ---------------- console table ----------------
    for frame, title in (("ft2", "FT2-RELATIVE (Tables 4/7 frame)"),
                         ("true", "TRUE-RELATIVE (INDELible true tree)")):
        print("\n" + "=" * 78)
        print("  %s   -- %d seeds" % (title, summary["data"]["n_seeds_scored"]))
        print("=" * 78)
        print("  %-42s %5s %-18s %10s %7s %s"
              % ("method", "n", "mean +/- SD", "p(zsplit)", "d", "wins/ties"))
        s = summary["vs_" + frame]
        for m in ALL_METHODS:
            e = s[m]
            if e["n"] == 0:
                print("  %-42s %5d  %s" % (METHOD_LABELS[m], 0, "NO DATA"))
                continue
            v = e["vs_fusang"]
            if v is None:
                print("  %-42s %5d  %.4f +/- %.4f   %10s %7s"
                      % (METHOD_LABELS[m], e["n"], e["mean"], e["sd"] or 0.0, "--", "--"))
            else:
                print("  %-42s %5d  %.4f +/- %.4f   %10.3g %7.2f  %d/%d"
                      % (METHOD_LABELS[m], e["n"], e["mean"], e["sd"] or 0.0,
                         v["wilcoxon_p_zsplit"] if v["wilcoxon_p_zsplit"] is not None
                         else float("nan"),
                         v["cohens_d_paired"], v["wins_vs_fusang"], v["ties"]))

    errs = [(r["seed"], r.get("error")) for r in rows if r.get("error")]
    if errs:
        print("\n  ERRORS (%d seeds):" % len(errs))
        for s, e in errs[:20]:
            print("    seed %d: %s" % (s, e))


if __name__ == "__main__":
    main()
