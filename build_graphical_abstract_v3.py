"""Fusang v3.0 Graphical Abstract - card-style, design-system colors."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.font_manager as fm

# ---- design system ----
INK    = "#1a2332"
WHITE  = "#ffffff"
CARD   = "#f4f7fa"
BLUE   = "#2563eb"
TEAL   = "#0d9488"
AMBER  = "#d97706"
PURPLE = "#7c3aed"
GREY   = "#64748b"
RED    = "#dc2626"

FIG_W, FIG_H = 13.2, 6.6
fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=300)
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 132); ax.set_ylim(0, 66); ax.axis("off")
fig.patch.set_facecolor(WHITE)

def card(x, y, w, h, fc, ec="none", lw=0, r=1.2):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                                fc=fc, ec=ec, lw=lw, mutation_aspect=1))

def txt(x, y, s, size, color=INK, weight="normal", ha="center", va="center", style="normal"):
    ax.text(x, y, s, fontsize=size, color=color, fontweight=weight, ha=ha, va=va,
            fontstyle=style, linespacing=1.35)

def arrow(x1, y1, x2, y2, color=GREY, lw=2.4):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=16,
                                 color=color, lw=lw, shrinkA=0, shrinkB=0))

# ---- title ----
txt(66, 63.2, "Fusang: Tardigrade Edition", 21, INK, "bold")
txt(66, 60.4, "Alignment-free phylogenetic inference with k-mer frequency vector cosine distances",
    11.5, GREY)

# ---- top row: 4-step pipeline ----
py, ph = 44.5, 13.5
steps = [
    (2.5,  BLUE,   "INPUT",          "Unaligned FASTA\n(DNA or protein)",            "no alignment step"),
    (34.5, TEAL,   "K-MER VECTORS",  "Spaced + contiguous\nk-mer frequency vectors", "auto k,gap selection\nmulti-k ensemble: no k choice"),
    (66.5, PURPLE, "COSINE DISTANCE","Pairwise cosine\ndistance matrix",             "boundary classifier:\nhomogeneous vs structured"),
    (98.5, AMBER,  "TREE",           "NJ / FastME\nBIONJ+BNNI \u2192 Newick",        "simplified \u2264500 taxa\nDCM >500 taxa"),
]
pw = 28.5
for x, c, head, body, foot in steps:
    card(x, py, pw, ph, CARD, c, 1.8)
    card(x, py + ph - 3.4, pw, 3.4, c)
    txt(x + pw/2, py + ph - 1.7, head, 11, WHITE, "bold")
    txt(x + pw/2, py + 6.6, body, 10.5, INK, "bold")
    txt(x + pw/2, py + 2.2, foot, 8, GREY, style="italic")
for i in range(3):
    arrow(2.5 + pw + 0.6 + i*32, py + ph/2, 2.5 + 32*(i+1) - 0.6, py + ph/2)

# ---- divider ----
txt(66, 41.2, "Validated across simulation and real data", 11, INK, "bold")

# ---- bottom row: three evidence cards ----
ey, eh, ew = 4.5, 34, 40
gap = (132 - 3*ew) / 4

# (i) indel-rich accuracy
x1 = gap
card(x1, ey, ew, eh, CARD, BLUE, 1.8)
txt(x1 + ew/2, ey + eh - 2.6, "Indel-rich accuracy", 12, BLUE, "bold")
txt(x1 + ew/2, ey + eh - 5.2, "n=200, indel rate 0.02, vs true tree (nRF \u2193)", 8, GREY)
axbar = fig.add_axes([(x1+6.5)/132, (ey+8.5)/66, (ew-13)/132, 16/66])
methods = ["Fusang", "FastTree2", "IQ-TREE2"]
vals    = [0.080, 0.085, 0.147]
colors  = [TEAL, GREY, GREY]
bars = axbar.barh(methods[::-1], vals[::-1], color=colors[::-1], height=0.62)
for b, v in zip(bars, vals[::-1]):
    axbar.text(v + 0.004, b.get_y() + b.get_height()/2, f"{v:.3f}", va="center", fontsize=7.5, color=INK)
axbar.set_xlim(0, 0.175); axbar.tick_params(labelsize=8, length=0, pad=2)
for s in ["top", "right", "bottom"]: axbar.spines[s].set_visible(False)
axbar.spines["left"].set_color("#cbd5e1")
axbar.set_xticks([])
txt(x1 + ew/2, ey + 5.6, "matches FastTree2 (p=0.052, n.s.-borderline)", 8.5, INK)
txt(x1 + ew/2, ey + 3.2, "1.8\u00d7 more accurate than IQ-TREE2 GTR (p<0.001)", 8.5, INK, "bold")

# (ii) real data
x2 = 2*gap + ew
card(x2, ey, ew, eh, CARD, TEAL, 1.8)
txt(x2 + ew/2, ey + eh - 2.6, "Real-data validation", 12, TEAL, "bold")
txt(x2 + ew/2, ey + eh - 5.2, "AFproject community benchmark, ETE3 nRF \u2193", 8, GREY)
rows = [
    ("Fish mitogenomes (n=25, ~17 kb)", "0.045", "community-benchmark ceiling \u2014 ties MSA+ML"),
    ("13 mtDNA genes (168\u20131,866 bp)",  "0.451", "Fusang\u2013FastTree2 agreement (13 loci)"),
    ("SwissTree proteins (11 families)", "0.239", "k-mer cosine vs Co-phylog 0.361 (p=0.006)"),
]
ry = ey + eh - 9.5
for name, v, note in rows:
    txt(x2 + 2.2, ry, name, 9, INK, "bold", ha="left")
    txt(x2 + ew - 2.2, ry, v, 13, TEAL, "bold", ha="right")
    txt(x2 + 2.2, ry - 2.6, note, 8, GREY, ha="left")
    ry -= 7.6
txt(x2 + ew/2, ey + 2.8, "no k selection at any scale: 500 bp genes \u2192 17 kb genomes", 8.5, INK, "bold")

# (iii) scalability + usability
x3 = 3*gap + 2*ew
card(x3, ey, ew, eh, CARD, AMBER, 1.8)
txt(x3 + ew/2, ey + eh - 2.6, "Fast, scalable, ready to use", 12, AMBER, "bold")
stats = [
    ("10,000 taxa", "54.4 s", "single CPU core"),
    ("1,000 taxa",  "5.2 s",  "simplified pipeline"),
    ("boundary classifier", "88/88", "E2E scenario routing"),
]
sy = ey + eh - 9.5
for name, v, note in stats:
    txt(x3 + 2.2, sy, name, 9, INK, "bold", ha="left")
    txt(x3 + ew - 2.2, sy, v, 13, AMBER, "bold", ha="right")
    txt(x3 + 2.2, sy - 2.6, note, 8, GREY, ha="left")
    sy -= 7.6
txt(x3 + ew/2, ey + 2.8, "pre-compiled Windows/Linux binaries \u00b7 MIT license", 8.5, INK, "bold")

# ---- footer ----
txt(66, 1.6, "github.com/zhanglknt/fusang-tardigrade", 9.5, BLUE, "bold")

fig.savefig("GraphicalAbstract_v3.png", dpi=300, facecolor=WHITE)
fig.savefig("GraphicalAbstract_v3.pdf", facecolor=WHITE)
print("saved GraphicalAbstract_v3.png / .pdf")
