#!/usr/bin/env python3
"""
Summarize all benchmark results into a readable markdown report.

Run after all generate_*.py scripts have completed.
"""
import sys, os, json, csv, statistics as st
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
DATA_DIR = Path(__file__).resolve().parent.parent / 'data'


def load_csv(path):
    """Load CSV and return list of dicts."""
    with open(path, newline='') as f:
        return list(csv.DictReader(f))


def format_nrf(val):
    if val is None or val == '':
        return 'N/A'
    return f"{float(val):.4f}"


def main():
    lines = []
    lines.append("# Reproducibility Benchmark Summary")
    lines.append(f"*Generated: 2026-06-18*")
    lines.append("")

    # ============ D1: Clean benchmark seeds 130-229 ============
    clean_csv = DATA_DIR / 'benchmark_clean_seeds130_229.csv'
    if clean_csv.exists():
        data = load_csv(clean_csv)
        valid = [r for r in data if r.get('status') == 'OK' and r.get('fusang_nrf')]
        nrf_vals = [float(r['fusang_nrf']) for r in valid]
        times = [float(r['fusang_time']) for r in valid]
        lines.append(f"## D1: Clean Benchmark (JC69, n=200)")
        lines.append(f"- Seeds: {len(data)} total, {len(valid)} valid")
        if nrf_vals:
            lines.append(f"- nRF (Fusang k=5,gap2): **{st.mean(nrf_vals):.4f} ± {st.stdev(nrf_vals):.4f}**")
        if times:
            lines.append(f"- Avg time: {st.mean(times):.1f}s per seed")
        lines.append(f"- Note: JC69 simulator (not INDELible). nRF values differ from paper.")
        lines.append("")
    else:
        lines.append("## D1: Clean Benchmark — NOT YET GENERATED")
        lines.append("")

    # ============ D4: Multi-k ensemble ============
    multik_csv = DATA_DIR / 'benchmark_multik_seeds200_229.csv'
    if multik_csv.exists():
        data = load_csv(multik_csv)
        valid = [r for r in data if r.get('status') == 'OK']
        lines.append(f"## D4: Multi-k Ensemble (JC69, n=200)")
        lines.append(f"- Seeds: {len(data)} total, {len(valid)} valid")
        for key in ['k5', 'k7', 'k9', 'ensemble']:
            vals = [float(r[f'{key}_nrf']) for r in valid if r.get(f'{key}_nrf') and r[f'{key}_nrf'] != '']
            if vals:
                label = {'k5': 'k=5,gap2', 'k7': 'k=7,gap2', 'k9': 'k=9,contiguous', 'ensemble': 'Ensemble (avg)'}[key]
                lines.append(f"- {label}: nRF={st.mean(vals):.4f} ± {st.stdev(vals):.4f}")
        best_key = None
        best_mean = float('inf')
        for key in ['k5', 'k7', 'k9', 'ensemble']:
            vals = [float(r[f'{key}_nrf']) for r in valid if r.get(f'{key}_nrf') and r[f'{key}_nrf'] != '']
            if vals and st.mean(vals) < best_mean:
                best_mean = st.mean(vals)
                best_key = key
        if best_key:
            lines.append(f"- **Best**: {best_key} (nRF={best_mean:.4f})")
        lines.append("")
    else:
        lines.append("## D4: Multi-k Ensemble — NOT YET GENERATED")
        lines.append("")

    # ============ D11: Scalability nRF ============
    snrf_json = DATA_DIR / 'scalability_nrf.json'
    if snrf_json.exists():
        with open(snrf_json) as f:
            snrf_data = json.load(f)
        lines.append(f"## D11: Scalability nRF (JC69)")
        lines.append(f"- Note: {snrf_data.get('note', '')}")
        summary = snrf_data.get('summary', {})
        for n_key in sorted(summary.keys(), key=int):
            s = summary[n_key]
            lines.append(f"- n={n_key}: nRF={s['nrf_mean']:.4f} ± {s['nrf_std']:.4f} "
                        f"({s['n_seeds']} seeds, {s['time_mean_s']}s avg)")
        lines.append("")
    else:
        lines.append("## D11: Scalability nRF — NOT YET GENERATED")
        lines.append("")

    # ============ Existing data summary ============
    lines.append("## Pre-existing Data (from paper)")
    lines.append("")

    # N200 indel
    indel_csv = DATA_DIR / 'indel_benchmark_seeds100_129.csv'
    if indel_csv.exists():
        data = load_csv(indel_csv)
        valid = [r for r in data if r.get('fusang_nrf') and r['fusang_nrf'] != '']
        if valid:
            fusang_vals = [float(r['fusang_nrf']) for r in valid]
            ft2_vals = [float(r['ft2_nrf']) for r in valid if r.get('ft2_nrf') and r['ft2_nrf'] != '']
            lines.append(f"### n=200 Indel (seeds 100-129, 30 seeds)")
            lines.append(f"- Fusang (k=5,gap2): nRF={st.mean(fusang_vals):.4f} ± {st.stdev(fusang_vals):.4f}")
            if ft2_vals:
                lines.append(f"- FT2: nRF={st.mean(ft2_vals):.4f} ± {st.stdev(ft2_vals):.4f}")
            lines.append("")

    # Scalability
    scal_json = DATA_DIR / 'scalability_results.json'
    if scal_json.exists():
        with open(scal_json) as f:
            scal = json.load(f)
        lines.append("### Scalability (runtime/memory)")
        lines.append("| n | Time (s) | RAM (MB) | Pipeline |")
        lines.append("|---|----------|----------|----------|")
        for entry in scal:
            lines.append(f"| {entry['n']} | {entry['total_s']} | {entry['ram_mb']} | {entry['pipeline']} |")
        lines.append("")

    # ============ External tools needed ============
    lines.append("## Remaining D-issues (Require External Tools)")
    lines.append("")
    lines.append("| # | Issue | Required |")
    lines.append("|---|-------|----------|")
    lines.append("| D6 | IQ-TREE2 benchmark | IQ-TREE2 binary (iqtree.org) |")
    lines.append("| D7 | Indel rate scan | INDELible (abacus.gene.ucl.ac.uk) |")
    lines.append("| D8 | AF method comparison | INDELible data + competitor implementations |")
    lines.append("| D9 | SwissTree full | SwissTree data (AFproject) |")
    lines.append("| D12 | 16S rRNA validation | 16S sequences (SILVA/Greengenes) |")
    lines.append("")

    report = '\n'.join(lines)
    output_path = DATA_DIR / 'benchmark_summary.md'
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(report)
    print(report)
    print(f"\nSaved to: {output_path}")


if __name__ == '__main__':
    main()
