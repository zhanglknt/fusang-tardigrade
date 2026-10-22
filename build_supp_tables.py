"""Build Supplementary_Tables.docx (Supplementary Tables S1-S17) for NAR submission v3.1.

Pattern follows build_v31_tables.py: assemble one markdown document (one caption line
`**Supplementary Table SN.** ...` + one markdown table per supplementary table) and
convert with pandoc. Read-only on all data sources; writes only
Supplementary_Tables.docx in the repository root.
"""
import csv
import json
import os
import re
import subprocess
import tempfile

SRC = r"D:\系统发育树项目\Fusang\Fusang-main"
PANDOC = r"C:\Users\admin\miniconda3\envs\human_selection\Library\bin\pandoc.exe"
OUT = os.path.join(SRC, "Supplementary_Tables.docx")


def p(*a):
    return os.path.join(SRC, *a)


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def fmt(x, nd=4):
    return f"{float(x):.{nd}f}"


_SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")


def fmt_p(x):
    x = float(x)
    if x >= 0.01:
        return f"{x:.4f}".rstrip("0").rstrip(".")
    if x >= 0.001:
        return f"{x:.4f}"
    ex = int(f"{x:.2e}".split("e")[1])
    mant = x / 10 ** ex
    sign = "⁻" if ex < 0 else ""
    return f"{mant:.2f}×10{sign}{str(abs(ex)).translate(_SUP)}"


def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |",
           "|" + "|".join("---" for _ in header) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


sections = []

# ---------------------------------------------------------------- S1
# 74 curated accessions from real_data/download_16s_robust.py (CURATED_16S),
# cross-referenced with real_data/real_16S_metadata.json (73 analysed).
src = open(p("real_data", "download_16s_robust.py"), encoding="utf-8").read()
m = re.search(r"CURATED_16S\s*=\s*(\[.*?\])\s*\n", src, re.DOTALL)
entries = re.findall(r'\(\s*"([^"]+)",\s*"([^"]+)",\s*"([^"]+)",\s*"([^"]+)"\s*\)', m.group(1))
assert len(entries) == 74, len(entries)
meta = json.load(open(p("real_data", "real_16S_metadata.json"), encoding="utf-8"))
meta_by_acc = {v["accession"]: v for v in meta.values()}
rows = []
for acc, org, phylum, order in entries:
    length = meta_by_acc.get(acc, {}).get("length", "—")
    rows.append([acc, org.replace("_", " "), phylum, order, length])
sections.append((
    "**Supplementary Table S1.** NCBI GenBank accessions for 74 16S rRNA sequences. "
    "Organism, phylum and order give the corrected NCBI Taxonomy annotation "
    "(Supplementary Note S17); length is in bp (73 sequences were successfully "
    "downloaded and analysed; one curated accession was not retrieved).",
    md_table(["Accession", "Organism", "Phylum", "Order", "Length (bp)"], rows)))

# ---------------------------------------------------------------- S2
rows = []
for r in csv.DictReader(open(p("indel_k_gap_scan_results.tsv"), encoding="utf-8"), delimiter="\t"):
    if r["status"] == "ok":
        rows.append([r["n_taxa"], r["k"], r["gap"], fmt(r["nrf"], 6), fmt(r["time_s"], 2)])
assert len(rows) == 50, len(rows)
sections.append((
    "**Supplementary Table S2.** Complete k/gap parameter scan on indel data "
    "(indel rate=0.02): nRF for k=4–8 × gap=none/gap1–gap4 at n=100/200. "
    "Single dataset per scale (n=100: `test_indel_n100_i002.fasta`; n=200: "
    "`test_indel_n200_i0.02.fasta`), Fusang pipeline (k-mer → FastME tree), "
    "nRF vs the true tree. Source: `indel_k_gap_scan_results.tsv`.",
    md_table(["n", "k", "gap", "nRF", "time (s)"], rows)))

# ---------------------------------------------------------------- S3
d = json.load(open(p("repro_ablations", "ablation_dcm_seed42_indel_nrfstandard.json"), encoding="utf-8"))
rows = [
    ["1. Simplified pipeline (k-mer cosine → NJ)", d["step1_simplified_nj"]["rf"], fmt(d["step1_simplified_nj"]["nrf"])],
    ["2. FastME on the same distance matrix", d["step2_fastme"]["rf"], fmt(d["step2_fastme"]["nrf"])],
    ["3. v2 distance matrix → NJ", d["step3_v2_distmat_nj"]["rf"], fmt(d["step3_v2_distmat_nj"]["nrf"])],
    ["4. TF-IDF weighting → NJ", d["step4_tfidf_nj"]["rf"], fmt(d["step4_tfidf_nj"]["nrf"])],
]
for name, v in d["step5_dcm_sweep_nj_default"].items():
    mg = name.replace("mg", "").replace("_nj_default", "")
    rows.append([f"5. DCM+EPA (max_group_size={mg}, NJ subtrees)", v["rf"], fmt(v["nrf"])])
rows.append(["MAFFT+FastTree2 (MSA+ML reference)", d["fasttree2_mafft"]["rf"], fmt(d["fasttree2_mafft"]["nrf"])])
sw = d["step5_dcm_range"]
assert abs(sw["mean"] - 0.2330) < 1e-3 and abs(sw["sd"] - 0.1260) < 1e-3
sections.append((
    "**Supplementary Table S3.** DCM degradation step-by-step: nRF at each pipeline stage "
    "(n=200, seed=42, indel rate=0.02), standard root-independent unrooted RF "
    "(max RF = 2(n−3) = 394). The step-5 rows give the full DCM sweep "
    "(max_group_size=10–150, NJ subtrees, 11 configurations; range "
    f"{fmt(sw['min'], 3)}–{fmt(sw['max'], 3)}, mean {fmt(sw['mean'], 3)} ± {fmt(sw['sd'], 3)}), "
    "so the conclusion that EPA grafting dominates the error does not depend on that choice. "
    "Source: `repro_ablations/ablation_dcm_seed42_indel_nrfstandard.json`. "
    "The stage values printed in earlier versions are withdrawn (see Results).",
    md_table(["Pipeline stage / configuration", "RF (splits)", "nRF"], rows)))

# ---------------------------------------------------------------- S4
rows = [[r["seed"], fmt(r["nrf_fusang_std"]), fmt(r["nrf_ft2_std"])]
        for r in read_csv(p("indelible_n200_indel", "bench_n200_indel_130seeds_nrfstandard.csv"))]
assert len(rows) == 130, len(rows)
fus = [float(r[1]) for r in rows]; ft2 = [float(r[2]) for r in rows]
import statistics as st
assert abs(st.mean(fus) - 0.1990) < 5e-4 and abs(st.mean(ft2) - 0.2047) < 5e-4, (st.mean(fus), st.mean(ft2))
sections.append((
    "**Supplementary Table S4.** 130-seed benchmark raw data: per-seed nRF for Fusang and "
    "FastTree2 at n=200, indel rate=0.02 (seeds 100–229, INDELible; standard "
    "root-independent unrooted RF vs the true tree; Fusang "
    f"{st.mean(fus):.4f} ± {st.stdev(fus):.4f}, FastTree2 {st.mean(ft2):.4f} ± {st.stdev(ft2):.4f}). "
    "Source: `indelible_n200_indel/bench_n200_indel_130seeds_nrfstandard.csv`.",
    md_table(["seed", "nRF Fusang", "nRF FastTree2"], rows)))

# ---------------------------------------------------------------- S5
# 5-repeat stability data as stated in the manuscript (Results, parameter stability)
# and encoded in generate_fig4_v3.py panel C: 5/5 repeats identical at every scale.
rows = [
    ["20", "4", "gap1", "5/5", "< 0.005"],
    ["50", "4", "gap1", "5/5", "< 0.005"],
    ["100", "4", "gap1", "5/5", "< 0.005"],
    ["200", "5", "gap2", "5/5", "< 0.005"],
    ["500", "5", "gap2", "5/5", "< 0.005"],
]
sections.append((
    "**Supplementary Table S5.** Adaptive parameter stability: 5-repeat validation of "
    "auto-selected k,gap across n=20/50/100/200/500. All 5 stability repeats at every "
    "scale returned identical k,gap selections and nRF within 0.005 units of each other "
    "(see Results, 'Parameter stability and reproducibility'; selection rule: k=4,gap1 "
    "for n≤100, k=5,gap2 for n>100). These five rows are the auto-selection outcomes "
    "recorded in the figure generator (`generate_fig4_v3.py`, panel C); no per-repeat "
    "raw records are archived for this check.",
    md_table(["n", "selected k", "selected gap", "repeats with identical selection",
              "max nRF deviation across repeats"], rows)))

# ---------------------------------------------------------------- S6
cond_files = [("n=500, clean", "n500_clean_nrfstandard.csv"),
              ("n=500, indel", "n500_indel_nrfstandard.csv"),
              ("n=1000, clean", "n1000_clean_nrfstandard.csv"),
              ("n=1000, indel", "n1000_indel_nrfstandard.csv")]
rows = []
for cond, fn in cond_files:
    for r in read_csv(p("table1_nrfstandard", fn)):
        rows.append([cond, r["seed"], r["n_taxa"], fmt(r["nrf_fusang"]), fmt(r["nrf_ft2"])])
assert len(rows) == 120, len(rows)
summ = json.load(open(p("table1_nrfstandard", "TABLE1_all_summary.json"), encoding="utf-8"))
srows = []
for key, cond in [("n500_clean", "n=500, clean"), ("n500_indel", "n=500, indel"),
                  ("n1000_clean", "n=1000, clean"), ("n1000_indel", "n=1000, indel")]:
    s = summ[key]
    srows.append([cond, f"{s['fusang_nrf_mean']:.4f} ± {s['fusang_nrf_sd']:.4f}",
                  f"{s['ft2_nrf_mean']:.4f} ± {s['ft2_nrf_sd']:.4f}",
                  fmt_p(s["wilcoxon_p"]), f"{s['cohens_dz']:.2f}",
                  f"{s['fusang_better']}/{s['ft2_better']}/{s['ties']}"])
sections.append((
    "**Supplementary Table S6.** Multi-seed benchmark at n=500 (30 seeds, seeds 500–529) "
    "and n=1000 (30 seeds, seeds 1000–1029): all conditions complete including n=1000 "
    "indel data (standard root-independent unrooted RF vs the true tree). All raw nRF "
    "values and statistical tests are provided for reproducibility. "
    "Sources: `table1_nrfstandard/{n500_clean,n500_indel,n1000_clean,n1000_indel}_nrfstandard.csv`.",
    md_table(["condition", "seed", "n", "nRF Fusang", "nRF FastTree2"], rows)
    + "\n\nSummary statistics (paired Wilcoxon on integer RF counts; wins = "
      "Fusang better / FastTree2 better / ties):\n\n"
    + md_table(["condition", "Fusang nRF", "FastTree2 nRF", "p", "Cohen's d", "wins (F/FT2/tie)"], srows)))

# ---------------------------------------------------------------- S7
rf = {r["family"]: r for r in read_csv(p("balibase_rf_recomputed", "balibase_rf_standard.csv"))}
diag = {d_["family"]: d_ for d_ in json.load(open(p("balibase", "diagnosis", "balibase_diagnosis.json"), encoding="utf-8"))}
rows = []
for fam, r in rf.items():
    dg = diag.get(fam, {})
    rows.append([fam, r["n"], dg.get("mean_length", "—") if not dg else f"{dg['mean_length']:.0f}",
                 f"{dg.get('mean_identity', float('nan')):.4f}" if dg else "—",
                 f"{dg.get('mean_gap_fraction', float('nan')):.4f}" if dg else "—",
                 dg.get("reference", "—") if dg else "—",
                 r["RF"], r["s1"], r["s2"], fmt(r["nRF_standard_new"])])
vals = sorted(float(r["nRF_standard_new"]) for r in rf.values())
med = (vals[9] + vals[10]) / 2
mean = sum(vals) / len(vals)
below = sum(1 for v in vals if v < 0.5)
assert abs(med - 0.982) < 5e-3 and abs(mean - 0.826) < 5e-3 and below == 3, (med, mean, below)
assert all(r["dendro_agrees"] == "True" for r in rf.values())
sections.append((
    "**Supplementary Table S7.** BAliBASE v3.0 benchmark: per-family Fusang-vs-FastTree2 "
    "nRF (standard root-independent unrooted RF; median "
    f"{med:.3f}, mean {mean:.3f}, {below}/20 families below 0.5), sequence statistics, "
    "and gap analysis for 20 protein families. No reference trees are available for this "
    "dataset, so these values measure between-method agreement, not accuracy. RF = "
    "differing splits between the Fusang and FastTree2 trees (s₁, s₂ = internal splits "
    "per tree); all values verified against DendroPy 5.0.8. Sources: "
    "`balibase_rf_recomputed/balibase_rf_standard.csv`, `balibase/diagnosis/balibase_diagnosis.json`.",
    md_table(["family", "n taxa", "mean length (aa)", "mean pairwise identity",
              "mean gap fraction", "BAliBASE ref", "RF", "s₁", "s₂",
              "nRF (Fusang vs FT2)"], rows)))

# ---------------------------------------------------------------- S8
order6 = ["n200_clean", "n200_indel", "n500_clean", "n500_indel", "n1000_clean", "n1000_indel"]
raw = [(k, float(summ[k]["wilcoxon_p"])) for k in order6]
m = len(raw)
# Benjamini-Hochberg adjusted p-values (step-up, monotone)
ranked = sorted(enumerate(raw), key=lambda t: t[1][1])
bh = [0.0] * m
runmin = 1.0
for rank, (idx, (_, pv)) in reversed(list(enumerate(ranked))):
    adj = min(pv * m / (rank + 1), 1.0)
    runmin = min(runmin, adj)
    bh[idx] = runmin
rows = []
print("S8 raw p-values:")
for i, (k, pv) in enumerate(raw):
    s = summ[k]
    direction = "FastTree2 better" if s["mean_diff_fusang_minus_ft2"] > 0 else "Fusang better"
    bonf = min(pv * m, 1.0)
    rows.append([k.replace("_", " "), f"{s['n_taxa']}", s["condition"], direction,
                 fmt_p(pv), fmt_p(bonf), fmt_p(bh[i]),
                 "yes" if bonf < 0.008 else "no"])
    print(f"  {k}: {pv:.6g}")
sections.append((
    "**Supplementary Table S8.** Multiple comparison correction: raw and adjusted p-values "
    "for 6 ground truth dataset comparisons (Fusang vs FastTree2; paired Wilcoxon "
    "signed-rank on integer RF counts, exact permutation value when the comparison is "
    "unanimous). Both Bonferroni (α≈0.008) and Benjamini-Hochberg FDR corrections applied. "
    "'Survives' = significant after Bonferroni correction at α=0.05/6≈0.008. "
    "Source: `table1_nrfstandard/TABLE1_all_summary.json`.",
    md_table(["dataset", "n", "condition", "direction", "raw p", "Bonferroni-adjusted p",
              "BH-adjusted p", "survives α≈0.008"], rows)))

# ---------------------------------------------------------------- S9 / S10
comp = read_csv(p("indelible_n200_indel", "competitors_n200_indel_130seeds_stdrf.csv"))
assert len(comp) == 130, len(comp)


def col(name):
    return [float(r[name]) for r in comp]


for c, want in [("cosine_k7_true", 0.1938), ("multik_ensemble_true", 0.1940),
                ("fusang_spaced_k5_true", 0.1990), ("fusang_spaced_k5_nj_true", 0.2003),
                ("cosine_k9_true", 0.2005), ("cosine_k5_true", 0.2078),
                ("cophylog_k19_true", 0.8812), ("kmacs_k3_true", 0.3656),
                ("kmacs_k5_true", 0.4125), ("kmacs_k10_true", 0.4517)]:
    got = st.mean(col(c))
    assert abs(got - want) < 5e-4, (c, got, want)

s9_cfg = [("Fusang full pipeline (k=5,gap2)", "fusang_spaced_k5"),
          ("Fusang NJ-only (k=5,gap2)", "fusang_spaced_k5_nj"),
          ("contiguous k=5", "cosine_k5"), ("contiguous k=7", "cosine_k7"),
          ("contiguous k=9", "cosine_k9"), ("multi-k ensemble (k=5,7,9)", "multik_ensemble")]
rows = []
for i, r in enumerate(comp):
    row = [r["seed"]]
    for _, base in s9_cfg:
        row.append(fmt(r[base + "_true"]))
    for _, base in s9_cfg:
        row.append(fmt(r[base + "_ft2"]))
    rows.append(row)
sections.append((
    "**Supplementary Table S9.** Multi-k ensemble benchmark: per-seed nRF for Fusang "
    "(k=5,gap2; NJ-only and full pipeline), individual contiguous k-mers (k=5,7,9), and "
    "ensemble average (n=200, indel rate=0.02, 130 seeds 100–229, INDELible; TRUE-relative "
    "and FT2-relative), scored with the standard root-independent unrooted RF. "
    "TRUE-relative means: contiguous k=7 0.1938 ± 0.0295, multi-k ensemble 0.1940 ± 0.0296, "
    "Fusang full pipeline 0.1990 ± 0.0296, spaced NJ-only 0.2003 ± 0.0296, contiguous k=9 "
    "0.2005 ± 0.0299, contiguous k=5 0.2078 ± 0.0314 (Table 4). Source: "
    "`indelible_n200_indel/competitors_n200_indel_130seeds_stdrf.csv`.",
    md_table(["seed"] + [f"{lbl} (vs TRUE)" for lbl, _ in s9_cfg]
             + [f"{lbl} (vs FT2)" for lbl, _ in s9_cfg], rows)))

s10_cfg = [("Co-phylog k=19", "cophylog_k19"), ("kmacs k=3", "kmacs_k3"),
           ("kmacs k=5", "kmacs_k5"), ("kmacs k=10", "kmacs_k10"),
           ("cosine k=5", "cosine_k5"), ("cosine k=7", "cosine_k7"),
           ("Fusang spaced NJ-only", "fusang_spaced_k5_nj"),
           ("Fusang spaced full pipeline", "fusang_spaced_k5"),
           ("multi-k ensemble", "multik_ensemble")]
rows = []
for r in comp:
    row = [r["seed"]]
    for _, base in s10_cfg:
        row.append(fmt(r[base + "_true"]))
    for _, base in s10_cfg:
        row.append(fmt(r[base + "_ft2"]))
    rows.append(row)
sections.append((
    "**Supplementary Table S10.** Alignment-free competitor comparison: per-seed nRF for "
    "Co-phylog (k=19), kmacs (k=3,5,10), k-mer cosine (k=5, k=7 contiguous), Fusang "
    "(k=5,gap2 spaced; NJ-only and full pipeline), and the multi-k ensemble on n=200 indel "
    "data (130 seeds 100–229, INDELible; TRUE-relative and FT2-relative), standard "
    "root-independent unrooted RF. TRUE-relative means: Co-phylog 0.8812 ± 0.0353, kmacs "
    "k=3 0.3656 ± 0.0389, k=5 0.4125, k=10 0.4517, Fusang full pipeline 0.1990 (Table 7). "
    "Source: `indelible_n200_indel/competitors_n200_indel_130seeds_stdrf.csv`.",
    md_table(["seed"] + [f"{lbl} (vs TRUE)" for lbl, _ in s10_cfg]
             + [f"{lbl} (vs FT2)" for lbl, _ in s10_cfg], rows)))

# ---------------------------------------------------------------- S11 / S12
sw_rows = read_csv(p("swisstree_rf_recomputed", "swisstree_rf_standard.csv"))
assert len(sw_rows) == 66, len(sw_rows)
corrected = {r["family"]: r for r in read_csv(p("benchmark_swisstree_corrected.csv"))}
gap_frac = {"ST001": ("98", "524", "0.362"), "ST002": ("54", "1127", "0.489"),
            "ST003": ("49", "563", "0.112"), "ST004": ("115", "930", "0.517"),
            "ST005": ("29", "488", "0.324"), "ST007": ("60", "687", "0.809"),
            "ST008": ("42", "860", "0.678"), "ST009": ("39", "259", "0.572"),
            "ST010": ("34", "428", "0.502"), "ST011": ("159", "1707", "0.682"),
            "ST012": ("21", "732", "0.451")}
methods = ["fusang_k5_gap2", "fusang_k4_gap1", "kmer_k5", "kmer_k4", "cophylog_k5", "cophylog_k11"]
m_label = {"fusang_k5_gap2": "Fusang spaced `11011`", "fusang_k4_gap1": "Fusang spaced `1011`",
           "kmer_k5": "k-mer cosine k=5", "kmer_k4": "k-mer cosine k=4",
           "cophylog_k5": "Co-phylog halfctx=5", "cophylog_k11": "Co-phylog halfctx=11"}
fam_nrf = {}
for r in sw_rows:
    fam_nrf.setdefault(r["family"], {})[r["method"]] = float(r["nRF_standard_new"])
# sanity vs Table 10 means
for meth, want in [("fusang_k5_gap2", 0.329), ("kmer_k5", 0.343), ("kmer_k4", 0.346),
                   ("fusang_k4_gap1", 0.354), ("cophylog_k11", 0.727), ("cophylog_k5", 0.786)]:
    got = st.mean(fam_nrf[f][meth] for f in fam_nrf)
    assert abs(got - want) < 5e-3, (meth, got, want)
rows = []
for fam in sorted(fam_nrf):
    c = corrected.get(fam, {})
    g = gap_frac.get(fam, ("—", "—", "—"))
    ref_splits = next(r["s2"] for r in sw_rows if r["family"] == fam and r["method"] == "fusang_k5_gap2")
    rows.append([fam, c.get("n_taxa", g[0]), c.get("avg_seq_len", "—"), g[2], ref_splits]
                + [fmt(fam_nrf[fam][mth]) for mth in methods])
sections.append((
    "**Supplementary Table S11.** SwissTree gene tree benchmark: per-family nRF for all "
    "methods, sequence statistics, and reference tree properties (AFproject standard; "
    "standard root-independent unrooted RF, max RF = 2(n−3)). Aligned-length gap "
    "fractions from MAFFT --auto alignments (Table 12 supplement); ref splits = internal "
    "branches of the SwissTree reference tree. Sources: "
    "`swisstree_rf_recomputed/swisstree_rf_standard.csv`, `benchmark_swisstree_corrected.csv`.",
    md_table(["family", "n taxa", "mean seq length (aa)", "alignment gap fraction",
              "ref tree splits"] + [m_label[m] for m in methods], rows)))

rows = []
for r in sw_rows:
    rows.append([r["family"], m_label.get(r["method"], r["method"]), r["n"], r["s1"], r["s2"],
                 r["RF"], r["shared"], fmt(r["nRF_standard_new"]), r["dendro_agrees"]])
sections.append((
    "**Supplementary Table S12.** SwissTree per-family nRF results for all methods "
    "(k-mer cosine, spaced k-mers, Co-phylog configurations): complete long-format table "
    "(11 families × 6 methods). s₁/s₂ = internal splits of the inferred/reference tree, "
    "RF = differing splits, shared = common splits; nRF = RF/(2n−6), standard "
    "root-independent unrooted RF, all values verified against DendroPy 5.0.8 "
    "(`dendro_agrees`). Source: `swisstree_rf_recomputed/swisstree_rf_standard.csv`.",
    md_table(["family", "method", "n", "s₁", "s₂", "RF", "shared", "nRF", "DendroPy agrees"], rows)))

# ---------------------------------------------------------------- S13
l3 = read_csv(p("bench_l3_coalescent_indelible_30seeds_nrfstandard.csv"))
assert len(l3) == 30, len(l3)
# NB: the *_true columns are the superseded root-dependent nRF; the standard
# root-independent unrooted RF lives in the *_true_std columns.
for c, want in [("l0_true_std", 0.150), ("l1_true_std", 0.137), ("l3_true_std", 0.124)]:
    got = st.mean(float(r[c]) for r in l3)
    assert abs(got - want) < 5e-3, (c, got, want)
rows = [[r["seed"], fmt(r["l0_true_std"]), fmt(r["l1_true_std"]), fmt(r["l3_true_std"])]
        for r in l3]
m0 = st.mean(float(r["l0_true_std"]) for r in l3)
m1 = st.mean(float(r["l1_true_std"]) for r in l3)
m3 = st.mean(float(r["l3_true_std"]) for r in l3)
s0 = st.stdev(float(r["l0_true_std"]) for r in l3)
s1d = st.stdev(float(r["l1_true_std"]) for r in l3)
s3d = st.stdev(float(r["l3_true_std"]) for r in l3)
sections.append((
    "**Supplementary Table S13.** L3 pipeline-level validation: per-seed nRF for L0 "
    "(single-k NJ), L1 (multi-k ensemble NJ), and L3 (MAFFT+FastTree2) against the true "
    f"simulated tree (n=200, L=1000 bp, indel=0.02, 30 seeds, seeds 300–329, all levels "
    f"complete; means L0 {m0:.3f} ± {s0:.3f}, L1 {m1:.3f} ± {s1d:.3f}, L3 {m3:.3f} ± {s3d:.3f}). "
    "Values are standard root-independent unrooted RF from the INDELible v1.03 "
    "regeneration; the superseded hand-written-simulator dataset is not reported (see "
    "Table 5 provenance note). Source: `bench_l3_coalescent_indelible_30seeds_nrfstandard.csv`. "
    "L3 executed on Linux (WSL2 Ubuntu-24.04, MAFFT --auto + FastTree2 -nt -gtr "
    "-nosupport); supersedes the earlier 5-seed Windows run (Supplementary Note S8).",
    md_table(["seed", "L0 (single-k NJ) nRF", "L1 (multi-k ensemble NJ) nRF",
              "L3 (MAFFT+FastTree2) nRF"], rows)))

# ---------------------------------------------------------------- S14
mash = [r for r in read_csv(p("bench_mash_indelible_n200_nrfstandard.csv"))
        if 100 <= int(r["seed"]) <= 129]
assert len(mash) == 60, len(mash)
by_seed = {}
for r in mash:
    by_seed.setdefault(int(r["seed"]), {})[r["condition"]] = r
assert all(set(v) == {"clean", "indel"} for v in by_seed.values())
rows = []
mc = [float(by_seed[s]["clean"]["mash_nrf"]) for s in sorted(by_seed)]
fc = [float(by_seed[s]["clean"]["fusang_nrf"]) for s in sorted(by_seed)]
mi = [float(by_seed[s]["indel"]["mash_nrf"]) for s in sorted(by_seed)]
fi = [float(by_seed[s]["indel"]["fusang_nrf"]) for s in sorted(by_seed)]
assert abs(st.mean(mc) - 0.3486) < 5e-4 and abs(st.mean(fc) - 0.2364) < 5e-4
assert abs(st.mean(mi) - 0.4034) < 5e-4 and abs(st.mean(fi) - 0.1926) < 5e-4
for s in sorted(by_seed):
    c_, i_ = by_seed[s]["clean"], by_seed[s]["indel"]
    rows.append([s, fmt(c_["mash_nrf"]), fmt(c_["fusang_nrf"]),
                 fmt(i_["mash_nrf"]), fmt(i_["fusang_nrf"]),
                 f"{float(i_['mash_nrf']) / float(c_['mash_nrf']):.3f}",
                 f"{float(i_['fusang_nrf']) / float(c_['fusang_nrf']):.3f}"])
sections.append((
    "**Supplementary Table S14.** Mash vs spaced k-mer benchmark (Table 2): per-seed nRF "
    "for Fusang (k=5, gap2 → `1110011`, cosine, NJ) and Mash (k=21, MinHash Jaccard, NJ) "
    "on clean and indel-rich data — n=200, L=500 bp, INDELible v1.03, 30 seeds 100–129, "
    "scored against the true tree with the standard root-independent unrooted RF "
    "(`nrf_standard.py`), together with the per-seed indel/clean ratios "
    f"(clean: Mash {st.mean(mc):.4f} ± {st.stdev(mc):.4f} vs Fusang {st.mean(fc):.4f} ± {st.stdev(fc):.4f}; "
    f"indel: Mash {st.mean(mi):.4f} ± {st.stdev(mi):.4f} vs Fusang {st.mean(fi):.4f} ± {st.stdev(fi):.4f}). "
    "Mash v2.3, `mash triangle -i -k 21 -s 10000`. **Provenance**: the values in the "
    "earlier version of this table came from `benchmark_mash/{clean,indel}/seed001–030.fasta` "
    "(defective hand-written coalescent simulator) and from Mash runs that omitted the "
    "required `-i` flag; those values are withdrawn. See the Table 2 provenance note in "
    "the main text. Source: `bench_mash_indelible_n200_nrfstandard.csv`.",
    md_table(["seed", "Mash nRF (clean)", "Fusang nRF (clean)", "Mash nRF (indel)",
              "Fusang nRF (indel)", "Mash indel/clean", "Fusang indel/clean"], rows)))

# ---------------------------------------------------------------- S15
e2e = json.load(open(p("e2e_indelible_results.json"), encoding="utf-8"))
assert e2e["total_scenarios"] == 150 and len(e2e["individual_results"]) == 150
rows = []
for r in e2e["individual_results"]:
    rows.append([r["tree_type"], r["n"], r["param"], r["seed"], r["expected"],
                 fmt(r["p_split"], 4), "yes" if r["correct"] else "no",
                 "yes" if r["complete"] else "no"])
acc = e2e["overall_accuracy"]
sections.append((
    "**Supplementary Table S15.** Boundary classifier E2E validation: all 150 individual "
    "scenario results including tree type, expected decision, predicted probability, and "
    f"correctness (overall accuracy {acc:.3f}, 95% CI "
    f"[{e2e['overall_ci_lower']:.3f}, {e2e['overall_ci_upper']:.3f}]; tree completeness "
    f"{e2e['tree_completeness']:.3f}). p_split = predicted probability of the SPLIT "
    "decision; scenarios are INDELible-simulated coalescent and structured trees across "
    "n=50–500 and five parameter regimes (default, low_indel, vlow_indel, short_seq, "
    "short_lowindel; 3 seeds per group). Source: `e2e_indelible_results.json`.",
    md_table(["tree type", "n", "regime", "seed", "expected decision", "p_split",
              "correct", "tree complete"], rows)))

# ---------------------------------------------------------------- S16
s16_md = open(p("real_data", "afproject_genome", "Supplementary_Table_S16.md"), encoding="utf-8").read()
s16_md = re.sub(r"^# Supplementary Table S16[^\n]*\n", "", s16_md).strip()
sections.append((
    "**Supplementary Table S16.** Real genome benchmarks (AFproject community datasets): "
    "complete per-method results for fish mtDNA (25 mitochondrial genomes, Fischer et al. "
    "2013 reference tree) and E. coli/Shigella (27 genomes, Skippington & Ragan 2011 "
    "reference tree), including NCBI accessions, k-scan results, multi-k variants, "
    "kmacs/Mash/Co-phylog runs at working and AFproject-published configurations, "
    "MAFFT+FastTree2 (fish), alignment gap-content quantification, and SwissTree "
    "per-family gap fractions; ETE3 unrooted nRF scoring.",
    s16_md))

# ---------------------------------------------------------------- S17
gd = json.load(open(p("real_data", "afproject_genome", "gene_dna", "gene_dna_results.json"), encoding="utf-8"))
rows = []
for gene, g in gd["genes"].items():
    rows.append([gene, f"{g['len_range_bp'][0]}–{g['len_range_bp'][1]}", g["aln_len"],
                 fmt(g["indel_gap_fraction"], 4),
                 fmt(g["fusang_multik_nrf_vs_ref"]), fmt(g["fusang_k5_nrf_vs_ref"]),
                 fmt(g["fusang_k7_nrf_vs_ref"]), fmt(g["ft2_nrf_vs_ref"]),
                 fmt(g["fusang_multik_vs_ft2_nrf"])])
assert len(rows) == 13, len(rows)
sections.append((
    "**Supplementary Table S17.** Gene-length real DNA benchmark: per-locus results for "
    "the 13 mitochondrial protein-coding genes (168–1,866 bp) from the 25 fish genomes — "
    "locus lengths, per-locus alignment gap fractions, ETE3 nRF vs the Fischer et al. 2013 "
    "whole-mitogenome reference tree, and Fusang-vs-FastTree2 agreement nRF, for the "
    "multi-k ensemble (mean of k=5,7,9 cosine distances), single k=5, single k=7, and the "
    "MAFFT+FastTree2 anchor. Source: `real_data/afproject_genome/gene_dna/gene_dna_results.json`.",
    md_table(["locus", "length range (bp)", "aligned length", "gap fraction",
              "multi-k nRF vs ref", "k=5 nRF vs ref", "k=7 nRF vs ref",
              "MAFFT+FT2 nRF vs ref", "multi-k vs FT2 nRF"], rows)))

# ---------------------------------------------------------------- emit
assert len(sections) == 17, len(sections)
all_md = "\n\n\\newpage\n\n".join(f"{cap}\n\n{tbl}" for cap, tbl in sections)
tmp = os.path.join(tempfile.gettempdir(), "supp_tables.md")
open(tmp, "w", encoding="utf-8").write(all_md)
r = subprocess.run([PANDOC, tmp, "-o", OUT, "--standalone"], capture_output=True, text=True)
print(f"pandoc rc={r.returncode} {r.stderr[:300]}")
print(f"{OUT}: {os.path.getsize(OUT):,} bytes")
os.remove(tmp)
