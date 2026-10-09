# Fusang v2.8 第四轮统计审稿意见（v4）

**审稿人**：v3-stats-reviewer（统计）
**稿件**：NAR_MANUSCRIPT_REVISED.md（v2.8）
**评审范围**：统计方法正确性、报告完整性、效应量解释、样本量充分性、选择性报告风险。所有引用数值均经独立复核（配对 d 反推隐含相关系数、p_adj 算术、异质性检验、Wilcoxon 正态近似特征值核对）。

---

## 总体印象

v2.8 是迄今统计质量最高的版本：统计框架首次明确划分 confirmatory（5 个 ground truth 数据集）与 exploratory 家族（L145），Table 7 重建消除了 seed set 嵌合问题（统一 100–126，数值自洽且 d 值反推隐含相关系数全部合理），L1 vs L3 的 30-seed Linux 重跑同时如实报告了未校正 p=0.024 与 Bonferroni 校正后 p_adj=0.071 borderline，且 E. coli 真实数据上的诚实负面结果（0.50 vs Mash 0.208）显著提升了可信度。剩余问题集中在三处：**(1) 摘要的两处 headline claim 均省略了正文已有的关键限定**（"significantly outperforms"未提 p_adj=0.071；multi-k "improves accuracy"未提在 27-seed 第二组上未复现且方向相反）；**(2) multi-k ensemble 效应存在真实的 seed-set 异质性**（Table 4 d=0.54 vs Table 7 d≈−0.22，异质性 z≈2.87、p≈0.004）而稿件未做统计处理；(3) Discussion 与 Table 1 的 n=500/n=1000 SD 存在未调和的矛盾。三者均为可在一轮修改内解决的文字/分析问题，不涉及新实验。

## 逐项打分

| 维度 | 分数 | 理由 |
|---|:---:|---|
| 统计方法正确性 | 8/10 | confirmatory/exploratory 框架明确、配对检验与 CI 正确、p_adj=0.024×3=0.072≈0.071 算术自洽；扣 multi-k 异质性未做统计处理、低功效检验推出等价结论 |
| 报告完整性 | 8/10 | 全文 mean±SD 均标注 seeds 数与 seed set、离群值机制透明；扣 Discussion n=500/1000 SD 与 Table 1 矛盾、置换检验未报次数 |
| 效应量解释 | 8/10 | d=11.91/3.82/2.54 均配对冲解释、L3 协议修正（`-nt` 遗漏）的披露堪称典范；扣 d 符号不一致（0.45 vs −0.45）、摘要限定缺失 |
| 样本量充分性 | 8/10 | 112/30 seeds 功效讨论到位、fish 单数据集以描述性呈现；扣 fish tie-at-floor 的机制叙述、n=27 等价声明未附功效 |
| 选择性报告风险 | 7/10 | 负面结果丰富（E. coli、16S、n≥500 clean 全部如实）；扣摘要只呈现显著 seed set 的 multi-k 结果、Co-phylog 竞争者结论基于自研重实现 |

## 相对上轮评审的变化

**已解决（相对 v3 评审及核销报告中未关闭项）**：
1. Minor 2（JSD 单尾 n=10）→ CLOSED：L73 改为 "small exploratory benchmark... informed the heuristic choice... should not be interpreted as a confirmatory test"，单尾表述与 p=0.031 已删除。
2. Minor 3（27 seeds 来源）→ CLOSED：现明确为 "seed set 100–126"（100–126 含 27 个整数，自洽），且修复了三组 seed 混用的嵌合体问题（我 v3 未识别此问题，v2.8 修复是重要改进）。
3. Minor 4（置换检验）→ PARTIAL：Table 9 注已明确 "two-sided permutation tests... with taxonomic group labels permuted"，但置换次数仍未给出。
4. Minor 6（分类器 100% vs CV AUC 0.84）→ CLOSED：L349 新增解释（E2E 场景占据参数空间中良好分离的区域、边界场景表现可能更低），推理合理。
5. Major 2（校正覆盖）→ CLOSED：L145 明确 "The 5 ground truth dataset comparisons constitute the pre-specified confirmatory family; all other reported p-values... are secondary and should be interpreted as exploratory"，且 pipeline 家族单独做 ×3 Bonferroni，框架自洽。
6. 上轮残留 N1（Table 3 列头 FT2-rel 矛盾）→ CLOSED：L259 列头已改 "(TRUE-rel)"。
7. 上轮残留 Major 1（Fusang σ=0.017）→ CLOSED：L278 统一为 σ=0.016。
8. 上轮 Minor 1（算术）→ CLOSED：13.2%（L268、L758）、4.5–4.8%（L268）、0.007（6.3%）（L298）均已修正。

**v2.7–v2.8 新增内容核查（team-lead 指定五项）**：
1. **L1 vs L3（30-seed Linux）**：核验通过。0.583±0.044 vs 0.601±0.055、差 −0.018、d=−0.45 → SD_diff≈0.040，隐含 r≈0.69 合理；Wilcoxon p=0.024 与 t≈2.47（p≈0.020）相容；p_adj=0.072≈0.071 报告值自洽；19/30 胜、旧 Windows run（漏 `-nt` 标志）supersession 披露完整（L551、S8）。正文 L324 如实并列两版 p 值并定位为 exploratory/borderline——**处理正确**；问题仅在摘要（见 Major 1）。
2. **Table 7 重建**：核验通过。统一 seed set 100–126、统一 FT2-relative 参照；新数值内部自洽：Co-phylog 0.408±0.021（d=11.91，隐含 r≈0.25）、kmacs 0.177±0.024（d=2.54，隐含 r≈0.25）、k5 contiguous 0.104（d=−0.21，隐含 r≈0.50）。两个 p=5.6×10⁻⁶ 相同不是笔误而是指纹特征——n=27 配对 Wilcoxon 全胜时正态近似双侧 p≈5.5–5.8×10⁻⁶，两比较均满足全 27 seeds 胜出，与均值差之大一致。multi-k 不一致的处理：L375 诚实披露 27-seed set 上 n.s.（p=0.30）并将显著结果归于 "pre-specified 30-seed benchmark"——披露合格，但构成实质的效应异质性问题（见 Major 2）。
3. **112-seed indel 基准**：与 v3 复核一致（源码级已验证 TRUE-relative；p=0.052 borderline、d=−0.20 [−0.42, 0.02] CI 跨零、60/112 胜），报告规范无回归。
4. **Table 12 真实数据**：无 p 值、无 SD（单数据集本性），正文 L494 明确 "n=25 taxa (44 informative splits), 0.045 corresponds to a 2-split difference"（2/44=0.045 算术正确），ETE3 与内部 calc_nrf 两种约定加了不可跨表比较的 caveat（L181、L488）——单点结果的限定总体到位；遗留问题见 Minor 3/4。
5. **交叉验证逻辑**：kmacs/Mash 复现已发表值（差 ≤1 split）能支撑的推理是 **benchmark harness（数据、NJ、ETE3 评分）正确**，不能延伸为对 Fusang 方法本身的验证。手稿用词 "cross-validating our pipeline"（L181）尚属恰当，但需防止读者误读（见 Minor 4）。

## 必须修改（Major）

**M1. 摘要两处 headline claim 缺失正文已有的关键限定（L26）**
- "significantly outperforms MAFFT+FastTree2 under indel-rich coalescent simulation (…p=0.024)"——正文 L324 与贡献段 L48 均已明确该比较为 secondary/exploratory 且 Bonferroni 后 p_adj=0.071 borderline，摘要却用无限定的 "significantly outperforms"。NAR 摘要是读者最先（且常常唯一）阅读的部分，此表述构成超claim。建议改为 "…(p=0.024 uncorrected; borderline after Bonferroni correction, p_adj=0.071)" 或等效措辞。
- "A multi-k distance ensemble improves accuracy without manual k selection (…p=0.006)"——该效应在 Table 7 的 27-seed set 上未复现且方向相反（见 M2），摘要应加限定（如 "in a pre-specified 30-seed benchmark"）或在修好 M2 后改写。

**M2. multi-k ensemble 效应的 seed-set 异质性需统计处理**
Table 4（seed set 230–259，n=30）：ensemble 0.105 vs default 0.112，d=0.54，p=0.006（ensemble 优）；Table 7（seed set 100–126，n=27）：ensemble 0.111 vs default 0.108，d=+0.22（即 default 优），p=0.30。同协议、同参照（FT2-relative）、不同 seed set，**效应符号相反**。异质性检验：z=(0.54−(−0.22))/√(1/30+1/27)≈2.87，p≈0.004——这不是抽样噪声能轻易解释的差异。当前处理（L375 把显著结果归于"pre-specified benchmark"）有选择性报告之嫌：pre-specification 界定的是检验而非 seed set。要求：(a) 合并两 seed set 做总分析（n=57，per-seed差值 pooled），或 (b) 报告异质性检验并明确声明"multi-k 对 default 的改进在不同 seed set 间不稳定"，相应调低摘要与贡献段中 multi-k 的 claim 强度。注意 Table 5 的 multi-k 巨大效应（d=3.82）来自不同协议（L=1000 bp coalescent），不能作为 L=500 bp 协议下效应存在的证据。

**M3. Discussion 与 Table 1 的 n=500/n=1000 SD 矛盾**
Discussion L541：n=500 clean "0.119 ± 0.020 vs 0.093 ± 0.015"、n=1000 clean "0.115 ± 0.022 vs 0.091 ± 0.016"；Table 1（L218/220）：同条件 "0.119 ± 0.011 vs 0.093 ± 0.013"、"0.115 ± 0.011 vs 0.091 ± 0.010"。同一基准两组 SD 必有一组错误。核查线索：d=1.47 与 Discussion 的 SD 组自洽（隐含 r≈0.52），与 Table 1 的 SD 组矛盾（隐含 r≈−0.08）——推测 Table 1 的 SD 被改小过或 Discussion 沿用旧值。须回源数据统一，并在 S6 核对全部派生统计量（d、p）与所公布的 SD 一致。

## 建议修改（Minor）

**m1. d 符号不一致**：贡献段 L48 "paired d=0.45" vs Results L324 "paired Cohen's d=−0.45"。按 L1−L3 约定应统一为 −0.45（或全文声明取绝对值）。

**m2. 置换检验补全**：Table 9 注已声明双侧与标签置换，但未给置换次数（建议 ≥10,000）与检验统计量定义，补一句即可。

**m3. fish 数据集的"平局地板"表述与事后机制解释**：(a) fish 上 multi-k、kmacs、Mash、Co-phylog（重实现）、MAFFT+FT2 全部 =0.045——该数据集对方法无区分度，建议明确写 "this dataset does not discriminate among these methods" 以免读者把并列最优读作方法优势；(b) multi-k(5,7,9)=0.045 优于 k=9 单值 0.091、且在含饱和 k=5/k=7 矩阵时反而最优（L496）的机制解释（"informative scale-specific variation dominates"）是单数据集上的事后叙述，应标注为推测（speculative）；(c) 摘要 fish claim 建议加 "single dataset" 限定。

**m4. 交叉验证推理的措辞**：kmacs/Mash 复现发表值验证的是 harness，建议 L181/L488 加半句 "validating the benchmark harness (data, NJ, and ETE3 scoring), not the accuracy of Fusang itself"；"difference ≤ one split, attributable to the NJ implementation" 的归因未经证明，宜软化为 "consistent with NJ implementation differences"。

**m5. Co-phylog DNA 失败结论基于 Python 重实现**：重实现在 fish 上得 0.045（vs 发表 0.09），与原版差异达 2×，说明重实现行为与原版有实质偏差。DNA 上 "3.8× worse, d=11.91" 的结论若基于重实现，虽方向保守（重实现不差于原版仍大败），但严谨做法是用原版二进制复核（如 kmacs 的处理方式），或在正文明确论证方向保守性。

**m6. 等价性声明的功效说明**：Table 7 "statistically indistinguishable"（k5 vs spaced，n=27，d=−0.21）：对 d=0.5 的功效仅约 0.73，检测不到中小效应；"provides no measurable advantage"（L375）宜改为 "no significant difference detected"，并附一句功效备注，避免读者将非显著读作等价。

**m7. 存储损坏恢复的验证材料**：L657 声明受影响文件 "regenerated deterministically... validated byte-identical"，建议在复现包中附原始与再生文件的 SHA-256 校验和对照表，使该声明可审计。

**m8. Table 8 n=200 indel 行**：数值为 30-seed（0.077±0.018）而 Notes 挂 112-seed p=0.052，虽已加说明文字（L396），建议格式上把两个基准分行或加脚注编号，避免读者误配。

## 总体推荐

**Minor Revision**

理由：v2.8 的统计纪律已达到发表标准的基础——预指定框架、双侧检验、配对效应量、诚实的负面结果与协议修正披露均到位，我前三轮提出的 2 Critical + 4 Major + 6 Minor 中除置换次数外已全部实质关闭。剩余 3 项 Major（M1 摘要限定、M2 multi-k 异质性处理、M3 SD 矛盾）均为摘要措辞修正、一次合并分析或数据核对，不需要新实验；其中 M2 是唯一涉及结论强度调整的实质问题，但手稿正文已完整披露相关数据，透明度充分。完成上述修改后，我认为统计方面可支持 Accept。

**投稿就绪度评分：85/100**（v2.5 复审 82 → v2.8 升至 85；如 M1–M3 完成预计可达 92+）
