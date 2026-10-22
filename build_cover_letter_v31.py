"""Build cover_letter_v3.1.docx for the NAR submission.

Single source of truth: NAR_MANUSCRIPT_REVISED.md (HEAD).
Every number, seed count and method name written into the letter was checked
against that file (see PROVENANCE below). Values deliberately NOT used because
the manuscript declares them withdrawn / superseded:

    0.147 (IQ-TREE2, old root-dependent nRF)
    0.419 (Co-phylog, chimaeric-era value)
    112 seeds / p=0.052 (legacy JC69 simulator)
    88 scenarios (superseded by 150 E2E scenarios)
    "six phyla" / "74 taxa" (16S; now 73 analysed taxa, 13 phyla)
    n=5 seeds (L3 Windows run; now 30 seeds on Linux)
    154 seeds / 24 seeds (IQ-TREE2 data-availability counts)
    609 MB (memory is now reported as "below 1 GB")
    1.005 / 0.592 / d=0.54 / 0.080 / 0.085 (all withdrawn)

Usage:
    python build_cover_letter_v31.py
"""
import os
import subprocess
import sys

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Pt, Cm, RGBColor

SRC = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(SRC, "cover_letter_v3.1.docx")

FONT = "Times New Roman"
PT = 10.0  # sized so the letter fits on one page

# ---------------------------------------------------------------------------
# PROVENANCE — every figure below, with the place it occurs in
# NAR_MANUSCRIPT_REVISED.md
# ---------------------------------------------------------------------------
PROVENANCE = [
    ("title", "L1 (manuscript title)"),
    ("130 seeds, n=200, L=500 bp", "Abstract; Table 1 caption; Table 3"),
    ("0.005: -8.9%, p=1.75e-11", "Table 3 + indel-rate scan paragraph"),
    ("0.01: -5.1%, p=5.20e-5", "Table 3 + indel-rate scan paragraph"),
    ("0.02: 0.199+-0.030 vs 0.205+-0.030, p=0.0507", "Table 1 / Table 3 / 130-seed section"),
    ("0.05: 0.180 vs 0.202, +10.9%, p=1.27e-11, 99/130", "Table 3 row 0.05"),
    ("ensemble 0.1940+-0.0296 vs 0.1990+-0.0296, p=7.7e-3, 72/130", "Table 4 + paired comparison"),
    ("contiguous k=7 0.1938, p=0.669", "Table 4 + paired comparison"),
    ("L=1000 bp 30 seeds: 0.137 vs 0.150, 8.7%, p=0.00206", "Table 5"),
    ("MAFFT+FastTree2 0.124, ~11% relative", "Table 5 + finding 1"),
    ("fish mtDNA 25 genomes, 9.1% gaps, nRF=0.045", "Table 12 + finding 1"),
    ("gene-scale default k=5 saturates at 0.545", "Table 12 + finding 2"),
    ("Mash clean 0.3486+-0.0351 vs 0.2364+-0.0336, d=3.12, 30/30", "Table 2 / Discussion finding 3"),
    ("Mash indel 0.4034+-0.0844 vs 0.1926+-0.0272, d=2.61, 30/30", "Table 2 / Discussion finding 2"),
    ("Mash degrades 1.16x, cosine improves 0.83x", "Abstract + Introduction contribution 3"),
    ("E. coli/Shigella plateau nRF=0.500", "Table 12 + finding 3"),
    ("10,000 taxa in 69.8 s", "Results scalability + Table 11"),
    ("peak memory below 1 GB", "Figure 5 legend (C)"),
    ("150 independent test scenarios", "Table 6 + Contributions 4"),
    ("73-taxon 16S rRNA dataset", "Data Availability"),
    ("GitHub / Zenodo 10.5281/zenodo.20746742", "Data Availability"),
    ("NSFC 32370682; 2026ZD01910500", "FUNDING"),
    ("corresponding author Li Zhang, knightz@pumc.edu.cn", "L20, CORRESPONDING AUTHOR"),
]


def stamp(run, size=PT, bold=False, italic=False):
    run.font.name = FONT
    run.font.size = Pt(size)
    run.font.bold = bold
    run.font.italic = italic
    run.font.color.rgb = RGBColor(0, 0, 0)
    rpr = run._element.get_or_add_rPr()
    rfonts = rpr.find(qn("w:rFonts"))
    if rfonts is None:
        rfonts = OxmlElement("w:rFonts")
        rpr.append(rfonts)
    for attr in ("w:ascii", "w:hAnsi", "w:eastAsia", "w:cs"):
        rfonts.set(qn(attr), FONT)


def para(doc, text, size=PT, bold=False, italic=False,
         space_after=4, align=None, first_line=0, left=0):
    p = doc.add_paragraph()
    pf = p.paragraph_format
    pf.space_after = Pt(space_after)
    pf.space_before = Pt(0)
    pf.line_spacing = 1.0
    if first_line:
        pf.first_line_indent = Cm(first_line)
    if left:
        pf.left_indent = Cm(left)
    if align is not None:
        p.alignment = align
    r = p.add_run(text)
    stamp(r, size=size, bold=bold, italic=italic)
    return p


def today():
    """Current date, taken from the system clock."""
    out = subprocess.run(["date", "+%d %B %Y"], capture_output=True, text=True)
    return (out.stdout or "").strip() or "2026"


def build():
    doc = Document()
    sec = doc.sections[0]
    sec.top_margin = Cm(1.6)
    sec.bottom_margin = Cm(1.6)
    sec.left_margin = Cm(2.0)
    sec.right_margin = Cm(2.0)

    normal = doc.styles["Normal"]
    normal.font.name = FONT
    normal.font.size = Pt(PT)

    d = today()

    # ---- header -----------------------------------------------------------
    para(doc, f"Date: {d}", space_after=2)
    para(doc, "To: The Editor, Nucleic Acids Research, Oxford University Press",
         space_after=2)
    para(doc, "Re: Submission of a Methods Article", space_after=8)

    para(doc, "Title: Fusang: Tardigrade Edition — K-mer Frequency Vector "
              "Alignment-Free Phylogenetic Inference Resilient to Indel-Rich "
              "Sequence Evolution",
         bold=True, space_after=10)

    para(doc, "Dear Editor,", space_after=5)

    # ---- nature of the manuscript ----------------------------------------
    para(doc,
         "We submit the above manuscript for consideration as a Methods Article in "
         "Nucleic Acids Research. It is a systematic evaluation of k-mer frequency "
         "vector cosine distances — spaced and contiguous patterns — for "
         "phylogenetic inference under indel-rich sequence evolution, delivered as "
         "open-source software (Fusang: Tardigrade Edition) with automatic parameter "
         "selection, a multi-k ensemble requiring no manual choice of k, and an "
         "empirically delimited domain of applicability. Three findings are of "
         "principal interest:",
         space_after=6)

    # ---- headline 1 -------------------------------------------------------
    para(doc, "1. Indel robustness is real but rate-dependent, and we locate the "
              "crossover.", bold=True, space_after=2)
    para(doc,
         "On INDELible-simulated gene-length data (n=200, L=500 bp, 130 independent "
         "seeds, scored against the true tree), Fusang's accuracy relative to "
         "MAFFT+FastTree2 improves monotonically with indel rate and crosses zero "
         "between 0.01 and 0.02: it is significantly worse at 0.005 (−8.9%, "
         "p=1.75×10⁻¹¹) and at 0.01 (−5.1%, p=5.20×10⁻⁵), numerically but not "
         "significantly better at 0.02 (nRF 0.199 ± 0.030 vs 0.205 ± 0.030, "
         "p=0.0507), and clearly better at 0.05 (+10.9%; 0.180 vs 0.202, p=1.27×10⁻¹¹, "
         "better in 99/130 seeds). The advantage is concentrated at high indel "
         "density, and we report it as such.",
         space_after=6)

    # ---- headline 2 -------------------------------------------------------
    para(doc, "2. The multi-k ensemble removes manual k selection and is the most "
              "accurate configuration we tested at gene length.",
         bold=True, space_after=2)
    para(doc,
         "On the same 130-seed benchmark the ensemble (contiguous k=5/7/9 averaged) "
         "attains nRF=0.1940 against 0.1990 for the default spaced-k plus NNI "
         "pipeline (p=7.7×10⁻³, better in 72/130 seeds), while tying the best single "
         "k (contiguous k=7, 0.1938; p=0.669). Under an independent 30-seed protocol "
         "at L=1000 bp, multi-k fusion cuts error by 8.7% relative to single-k (0.137 "
         "vs 0.150, p=0.00206), reaching within ~11% relative error of "
         "MAFFT+FastTree2 (0.124) with no alignment step. On real indel-rich data — "
         "25 fish mitochondrial genomes, 9.1% gap cells in the alignment — multi-k "
         "cosine reaches nRF=0.045, identical to MAFFT+FastTree2 and to the best "
         "published alignment-free configurations, where the gene-scale default (k=5) "
         "saturates at 0.545.",
         space_after=6)

    # ---- headline 3 -------------------------------------------------------
    para(doc, "3. “Alignment-free” is not one category, and we map its boundaries.",
         bold=True, space_after=2)
    para(doc,
         "At gene length (L=500 bp) MinHash Jaccard (Mash, k=21) is far less accurate "
         "than k-mer cosine in both conditions (clean 0.3486 vs 0.2364, paired "
         "d=3.12; indel 0.4034 vs 0.1926, d=2.61; Fusang better in 30/30 seeds in "
         "both), and the two respond to indels in opposite directions (Mash degrades "
         "1.16× while k-mer cosine improves 0.83×). We also state where k-mer cosine "
         "should not be preferred: strain-level bacterial genomes (E. coli/Shigella, "
         "plateau nRF=0.500) and structurally conserved RNA genes.",
         space_after=6)

    # ---- transparency -----------------------------------------------------
    para(doc, "Transparency.", bold=True, space_after=2)
    para(doc,
         "This version supersedes values carried by earlier iterations of the work — "
         "data from a coalescent simulator we regard as defective, and figures from a "
         "root-dependent RF implementation — and states on its face which claims are "
         "withdrawn or reversed. All accuracy figures use the standard "
         "root-independent unrooted RF (RF/2(n−3)), with Wilcoxon signed-rank tests "
         "on integer RF counts.",
         space_after=6)

    # ---- availability -----------------------------------------------------
    para(doc, "Software and data availability.", bold=True, space_after=2)
    para(doc,
         "Fusang: Tardigrade Edition is released under the MIT license with "
         "pre-compiled FastME binaries (Windows/Linux x86-64), benchmark scripts and "
         "analysis code at https://github.com/zhanglknt/fusang-tardigrade, archived "
         "at https://doi.org/10.5281/zenodo.20746742. All benchmark datasets — the "
         "130-seed n=200 indel results, the n=500 and n=1,000 multi-seed data, the "
         "30-seed Linux MAFFT+FastTree2 comparison, the Mash, kmacs and Skmer runs, "
         "the AFproject real-genome benchmarks, the 150-scenario classifier "
         "validation, the 73-taxon 16S rRNA dataset and the BAliBASE v3.0 results — "
         "are in the repository and as Supplementary Data. Inference on 10,000 taxa "
         "completes in 69.8 s with peak memory below 1 GB.",
         space_after=6)

    para(doc, "No competing interests are declared. Support: National Natural Science "
              "Foundation of China (32370682) and the Prevention and Control of "
              "Emerging and Major Infectious Diseases — National Science and "
              "Technology Major Project (2026ZD01910500).",
         space_after=6)

    para(doc, "We confirm that this manuscript has not been published elsewhere and is "
              "not under consideration by any other journal, and that all authors have "
              "read and approved it.",
         space_after=8)

    para(doc, "Sincerely,", space_after=1)
    para(doc, "Li Zhang, Ph.D. — corresponding author, on behalf of all authors "
              "(first author: Lei Kong).", space_after=1)
    para(doc, "Institute of Blood Transfusion, Chinese Academy of Medical Sciences and "
              "Peking Union Medical College; 26 Huacai Road, Longtan Industrial Park, "
              "Chenghua District, Chengdu 610052, Sichuan, China. "
              "Email: knightz@pumc.edu.cn", space_after=0)

    doc.save(OUT)
    return OUT


if __name__ == "__main__":
    path = build()
    print("provenance of every figure used:")
    for item, where in PROVENANCE:
        print(f"  - {item}: {where}")
    print(f"\nwrote {path}")
    sys.exit(0)
