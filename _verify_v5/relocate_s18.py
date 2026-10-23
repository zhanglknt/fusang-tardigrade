"""Item H: relocate provenance/withdrawal audit detail into new Note S18.

Moves (body -> Supplementary Note S18), leaving one-sentence pointers:
  a) Methods provenance subsection (3 paragraphs: Table 1 recompute, Table 3
     rescore, n=1000 tree reproducibility)
  b) L221 ablation-chain provenance tail (keeps EPA-style sentence in body)
  c) Mash benchmark provenance block (4 paragraphs)
  d) Table 5 footnote provenance tail
Also fixes a broken internal cross-reference at the end of the L221 tail.
"""
import sys

PATH = r"D:/系统发育树项目/Fusang/Fusang-main/NAR_MANUSCRIPT_REVISED.md"

lines = open(PATH, encoding="utf-8").read().split("\n")

# ---------- locate blocks by content ----------
i_m = next(i for i, l in enumerate(lines)
           if l.strip() == "#### Provenance of the accuracy and indel-rate benchmarks")
assert lines[i_m + 2].startswith("All six cells of Table 1"), lines[i_m + 2][:60]
assert lines[i_m + 4].startswith("The three indel-rate rows"), lines[i_m + 4][:60]
assert lines[i_m + 6].startswith("**Reproducibility of the n=1,000 trees.**")

i221 = next(i for i, l in enumerate(lines)
            if l.startswith("For calibration, FastTree2 on the MAFFT alignment"))
i261 = next(i for i, l in enumerate(lines)
            if l.startswith("**Data provenance and withdrawal.** Two independent"))
assert lines[i261 + 2].startswith("1. **Defective simulator.**")
assert lines[i261 + 4].startswith("2. **Defective Mash execution.**")
assert lines[i261 + 6].startswith("**Metric and implementation stability.**")
i345 = next(i for i, l in enumerate(lines)
            if l.startswith("nRF=0: perfect match; nRF=1.0: random tree; nRF = RF/(2(n\u22123)). L0 and L1"))

# ---------- extract ----------
t_a = [lines[i_m + 2], lines[i_m + 4], lines[i_m + 6]]

t221 = lines[i221]
cut = t221.find("**Provenance**")
head221, tail221 = t221[:cut], t221[cut:]
keep_sent = ("The implementation is a deterministic, EPA-*style* grafting: "
             "each subtree root replaces its group's representative leaf on "
             "the backbone, with no likelihood-based placement.")
broken = ("Multi-seed validation confirms the directional advantage at smaller "
          "magnitude (see Section \"130-seed benchmark confirms directional "
          "advantage\").")
assert keep_sent in tail221 and broken in tail221
move221 = tail221.replace(keep_sent, "").replace(broken, "").strip()

t345 = lines[i345]
cut345 = t345.find("**Data provenance and withdrawal.**")
assert cut345 > 0
head345, tail345 = t345[:cut345], t345[cut345:]

t_c = [lines[i261], lines[i261 + 2], lines[i261 + 4], lines[i261 + 6]]
t_c[0] = t_c[0].replace("**Data provenance and withdrawal.** ", "")

# ---------- build Note S18 ----------
s18 = ("**Supplementary Note S18.** Relocated provenance and withdrawal audit "
       "detail. *(a) Accuracy and indel-rate benchmark recomputation.* "
       + " ".join(t_a)
       + " *(b) DCM ablation-chain provenance.* " + move221
       + " *(c) Mash benchmark provenance.* " + " ".join(t_c)
       + " *(d) Pipeline-level validation provenance.* " + tail345)

# ---------- body replacements ----------
summary_m = ("All six accuracy-benchmark cells and the three indel-rate rows "
             "were recomputed from the archived trees under the standard "
             "root-independent RF; the full audit trail \u2014 including the "
             "n=1,000 archived-tree reproducibility analysis \u2014 is in "
             "Supplementary Note S18.")
lines[i_m + 2:i_m + 7] = [summary_m, ""]

i221 = next(i for i, l in enumerate(lines)
            if l.startswith("For calibration, FastTree2 on the MAFFT alignment"))
lines[i221] = (head221 + keep_sent + " Provenance and configuration details "
               "for this chain (the superseded clean-data trace, the "
               "`max_group_size` specification, and an unreproduced informal "
               "seed-100 trace) are in Supplementary Note S18; the "
               "recomputation protocol is documented in Supplementary Note "
               "S12. Multi-seed validation confirms the directional advantage "
               "at smaller magnitude (see the 130-seed benchmark section "
               "below).")

i261 = next(i for i, l in enumerate(lines)
            if l.startswith("**Data provenance and withdrawal.** Two independent"))
pointer_c = ("**Data provenance and withdrawal.** Two independent provenance "
             "problems affected the earlier version of this table \u2014 a "
             "defective hand-written simulator and a defective Mash driver "
             "invocation; both were corrected by INDELible regeneration and a "
             "fixed Mash driver, and the full audit detail is in "
             "Supplementary Note S18.")
lines[i261:i261 + 7] = [pointer_c, ""]

i345 = next(i for i, l in enumerate(lines)
            if l.startswith("nRF=0: perfect match; nRF=1.0: random tree; nRF = RF/(2(n\u22123)). L0 and L1"))
lines[i345] = head345 + "(Provenance and withdrawal details: Supplementary Note S18.)"

# ---------- insert Note S18 after Note S17 ----------
i17 = next(i for i, l in enumerate(lines)
           if l.startswith("**Supplementary Note S17.**"))
lines[i17 + 1:i17 + 1] = ["", s18]

open(PATH, "w", encoding="utf-8", newline="\n").write("\n".join(lines))
print("H relocation OK; total lines now", len(lines))
