# v2.9 修复独立交叉验证报告（fresh-eyes verifier）

**验证人**: v4-verifier（未参与修复，独立重查）
**稿件**: NAR_MANUSCRIPT_REVISED.md（v2.9，774 行）
**对照**: review_v4_stats.md / review_v4_general.md / review_v4_phylo.md / review_v4_format.md + REVISION_CHANGELOG.md v2.9 节
**方法**: 全文重读 + scipy 独立复算（pooled n=57、n=500/1000）+ Crossref API 全量 DOI 审计 + 脚本化全局一致性检查
**日期**: 2026-10-09

---

## 逐项核销表

### A. Major 项（7 项）

| # | 项目 | 判定 | 证据 |
|---|------|------|------|
| A1 | Abstract p=0.024 加 borderline/p_adj=0.071 限定 | **CLOSED** | L26: "outperforming MAFFT+FastTree2 on the same benchmark (p=0.024; borderline after Bonferroni correction, p_adj=0.071)"；"significantly" 已删 |
| A2 | Abstract multi-k 声称不再引 Table 4 p=0.006，改用 Table 5 | **CLOSED** | L26: "reducing topological error under indel-rich coalescent simulation (nRF=0.583±0.044 vs 0.743±0.046 single-k, 30 seeds, p<0.0001)"；Abstract 中的 "p=0.006" 仅为 SwissTree 跨域比较（合法保留） |
| A3 | fish 天花板表述四处齐改 | **CLOSED** | Abstract L26 "reaches the community-benchmark ceiling … this single small dataset"；Table 12 finding 1 L496 "community-benchmark ceiling … This dataset does not discriminate among these methods"；Limitations L547 "community-benchmark ceiling … matched by all top methods"；Practical recs L573 "reaches the community-benchmark ceiling" |
| A4 | 新增 "Reproducibility of the ensemble effect across seed sets" 段 | **CLOSED** | L377 独立段落，含 heterogeneity z=−2.71 p=0.007 与 pooled n=57（mean diff −0.003, d=−0.17, CI [−0.45, 0.09], Wilcoxon p=0.19, 32/57），并明确 "We therefore do not claim a general accuracy advantage … at L=500 bp"。独立复算见下节：全部统计量吻合（z 值有舍入级偏差，见新发现问题 #3） |
| A5 | Discussion n=500/1000 SD 与 d 修正 | **CLOSED** | L543: "0.119 ± 0.011 vs 0.093 ± 0.013 (Cohen's d=1.77)"、"0.115 ± 0.011 vs 0.091 ± 0.010 (Cohen's d=2.81)"，p<0.001 保留。与 Table 1（L218/L220）一致。独立复算完全吻合（见下节） |
| A6 | Table 5 L=1000 bp 动机 + S8 协议参数 | **CLOSED** | L307: "the coalescent protocol was fixed during multi-layer pipeline development — before the L1-vs-L3 comparison was designed — with L=1000 bp chosen so that the k=9 member … has non-sparse frequency vectors … identical for the superseded 5-seed run and the full 30-seed rerun"；Supp Note S8 描述 L644 含完整协议参数（n=200, L=1000 bp, sub=0.05, indel=0.02, MAFFT+FastTree2 配置） |
| A7 | Data Availability storage 披露精确化 | **CLOSED** | L659: 受影响 "21 competitor-benchmark input files (5 FASTA and 16 tree files, seeds 109–126)"（21=5+16 自洽）；三条恢复路径；"regenerated files validated byte-identical across the independent recovery paths"（消解了"与何比对"的逻辑矛盾）；"53/54 affected per-seed nRF values match the pre-incident master CSV (one FastTree2 value differs by two bipartitions, 0.0787→0.0838)"；隔离（quarantined）+ 校验和记录（recovery log with checksums in the reproducibility package） |

### B. Minor 项抽查

| 来源 | 项目 | 判定 | 证据 |
|------|------|------|------|
| stats m1 | d=−0.45 统一 | **CLOSED** | L48 "paired d=−0.45" 与 L324 "paired Cohen's d=−0.45" 一致；全文无正号 0.45 残留（脚本验证） |
| stats m2 | Table 9 注 10,000 permutations + rerun 披露 | **CLOSED** | L422: "10,000 random permutations per level" + 检验统计量定义；验证 rerun（73/74 taxa）复现 order/phylum 信号（7.3%/3.8%, p=0.027/0.003）；family 标签未归档、该行披露为原始分析 |
| stats m3 | fish 平局地板三件套 | **CLOSED** | (a) L496 "This dataset does not discriminate among these methods"；(b) L498 机制标注 "A plausible but unverified mechanism (this is a single dataset, and the explanation is post hoc)"；(c) Abstract L26 "this single small dataset" |
| stats m4 | harness 措辞 | **CLOSED** | L181 与 L490 两处均为 "cross-validating the benchmark harness (data, NJ, and ETE3 scoring) rather than of Fusang itself"；归因软化为 "consistent with NJ implementation differences" |
| stats m5 | Co-phylog 重实现方向保守论证 | **CLOSED** | L371: 重实现 fish 上优于发表配置 → "the DNA conclusion here is therefore conservative in direction" |
| stats m6 | no significant difference + 功效注 | **CLOSED** | L375: "no significant difference from spaced k-mers … With n=27, power to detect a medium effect (d=0.5) is approximately 0.73, so equivalence should not be inferred" |
| stats m7 | SHA-256 校验和对照表 | **CLOSED（文本级）** | L659 声明 "recovery log with checksums is documented in the reproducibility package"；repro package 内容本身不在本次核验范围 |
| stats m8 | Table 8 n=200 indel 行 | **CLOSED** | L398: "Fusang better (indel) ‡" + L402 脚注 ‡ 说明 30-seed vs 112-seed 基准区别 |
| general | "88/88 simulated scenarios" | **CLOSED** | L26 Abstract: "(88/88 simulated scenarios)" |
| general | Table 12 "n/a ³" 与 "—" 定义 | **CLOSED** | L483 "n/a ³"（kmacs E. coli 超时，脚注³ L492 详释）；L490 表注定义 '"—" = not tested (see notes); "n/a" = not applicable' |
| general | 基因长度 gap 候选数据类型 | **CLOSED** | L547 末尾: "candidate data types for closing this gap include viral quasispecies collections and rapidly evolving gene families with curated reference trees" |
| phylo | L42 中性化 | **CLOSED** | L42 加粗句改为 "**k-mer frequency vectors (spaced and contiguous) have not been systematically evaluated…**"；贡献点 (1) 改为 "k-mer frequency vectors with cosine distance"；与 Discussion L512 定位一致 |
| phylo | "hypothesized"×2 | **CLOSED** | L527 "hypothesized robustness from spaced sampling"、L529 "spaced sampling is hypothesized to provide additional robustness"；"theoretical motivation"/"inherent robustness" 零残留（脚本验证） |
| phylo | "among the most robust" | **CLOSED** | L563: "among the most robust alignment-free approaches for gene-length phylogenetic inference" |
| phylo Minor 1 | Co-phylog 重实现为何更好的机制 | **PARTIAL** | L371/L490 已加方向保守论证与 Table 7↔Table 12 交叉引用，但未解释重实现优于发表值的具体机制（如 NJ 实现差异）；审稿人主诉求（避免读者怀疑对标不等价）已实质满足 |
| phylo Minor 5 | boundary classifier on/off 或 Practical recs 说明 | **NOT-FIXED** | Practical recommendations（L567–576）未提及 classifier；仍仅有 L559 "feature ablation remains future work"（低严重度遗留项，changelog 未声称修复此条） |
| phylo Minor 6 | FFP/FSWM 未测原因 | **CLOSED** | L563: "(both were outside the scope of the current benchmark harness)" |
| format Major 1 | 字数声明 | **CLOSED** | L774: "approximately 11,200 words"。实测正文（Intro–Practical recs，不含表格行）11,380 词，1.6% 低估，在 "approximately" 容差内（备注见新发现问题 #6） |
| format Minor 3 | Table 12 脚注¹锚点 | **CLOSED** | L473 表题 "ETE3 nRF¹"，L490 脚注 "¹ nRF computed with ETE3…" |
| format Minor 4 | 补充材料三组分块 | **CLOSED** | L584–596 Fig S1–S7 → L598–628 Table S1–S16 → L630–648 Note S1–S10，脚本验证仅两次类型转换（Figure→Table→Note）且各组内严格递增、连续无缺 |
| format Minor 7 | TRUE→true | **CLOSED** | 15 处 bare TRUE 全部为 TRUE-relative/TRUE-rel 参照系标签或 L223 定义式 "ground-truth (TRUE) tree"（按修复方案保留）；正文散文无 "TRUE tree"；NOT×2 零残留（脚本验证） |
| format Minor 8 | borderline 去重 5→3 | **CLOSED** | L1-vs-L3 的 borderline 仅剩 L26（Abstract）+ L324（Results）+ L48 指针（"see Results for multiple-comparison context"，不带词）；L50/L276/L402 的 borderline 属 112-seed p=0.052 不同比较，L349 属分类器决策边界语境，均不在清理范围 |

---

## 独立复算结果

### 1. pooled n=57（multi-k ensemble vs default spaced，L=500 bp）

数据源：`table7_recomputed_recovered.json`（seed set 100–126, n=27, 字段 fusang/multik）+ `benchmark_multik_ensemble_n200_indel.csv`（seed set 230–259, n=30, 列 nrf_fusang_original/nrf_multik_ensemble）。配对差值 = multik − fusang（负 = ensemble 更优）。scipy.stats.wilcoxon（zero_method='wilcox'），bootstrap CI 10,000 次（seed=42）。

| 统计量 | 手稿值（L377） | 独立复算 | 一致性 |
|--------|--------------|----------|--------|
| Set A (n=27) mean | 0.111 vs 0.108 | 0.1105 vs 0.1075 | ✓ |
| Set A paired d | +0.22 | +0.219 | ✓ |
| Set A Wilcoxon p | 0.30 | 0.2996 | ✓ |
| Set B (n=30) mean | 0.105 vs 0.112 | 0.1045 vs 0.1121 | ✓ |
| Set B paired d | −0.53 | −0.529 | ✓ |
| Set B Wilcoxon p | 0.006 | 0.0057 | ✓ |
| Pooled n=57 mean diff | −0.003 | −0.00255 | ✓ |
| Pooled paired d | −0.17 | −0.171 | ✓ |
| Pooled d 95% CI | [−0.45, 0.09] | [−0.451, 0.088] | ✓ |
| Pooled Wilcoxon p | 0.19 | 0.1945 | ✓ |
| Ensemble wins | 32/57 | 32/57 (56.1%) | ✓ |
| Heterogeneity z | −2.71, p=0.007 | **−2.82, p=0.0048** | 方向/显著性一致，数值有舍入级偏差（见新发现 #3） |

**结论**：pooled 分析全部统计量与手稿吻合；手稿的 claim 降级（"do not claim a general accuracy advantage … at L=500 bp"）与数据一致。

### 2. n=500 / n=1000 clean（Discussion L543 vs Table 1）

数据源：`benchmark_n500_clean_30seeds.csv`、`benchmark_n1000_clean_30seeds.csv`（列 nrf_ft2/nrf_fusang）。

| 条件 | 量 | 手稿值 | 独立复算 | 一致性 |
|------|-----|--------|----------|--------|
| n=500 clean | Fusang | 0.119 ± 0.011 | 0.1185 ± 0.0111 | ✓ |
| n=500 clean | FT2 | 0.093 ± 0.013 | 0.0930 ± 0.0125 | ✓ |
| n=500 clean | paired d | 1.77 | 1.770 | ✓ |
| n=500 clean | Wilcoxon p | <0.001 | 3.1×10⁻⁶ | ✓ |
| n=1000 clean | Fusang | 0.115 ± 0.011 | 0.1147 ± 0.0109 | ✓ |
| n=1000 clean | FT2 | 0.091 ± 0.010 | 0.0907 ± 0.0097 | ✓ |
| n=1000 clean | paired d | 2.81 | 2.807 | ✓ |
| n=1000 clean | Wilcoxon p | <0.001 | 1.9×10⁻⁶ | ✓ |

**结论**：v2.9 修正后的 Discussion 数值与 Table 1 及源数据三方一致，v2.8 的 SD 矛盾已消除；旧值（0.119±0.020、0.115±0.022 于 n=1000、d=1.47、d=1.26）全文零残留（脚本验证；L396 的 "0.115 ± 0.022" 为 Table 8 n=100 clean 条件，非残留——已确认）。

---

## Crossref 抽查结果

方法：对 REFERENCES 全部 31 条中的 30 条 DOI 逐条调用 `https://api.crossref.org/works/{DOI}`（带 User-Agent），比对标题、作者、年份、卷页。

| Ref | DOI | 解析 | 标题一致 | 备注 |
|-----|-----|------|---------|------|
| 1 | 10.1093/bioinformatics/bty407 | OK | ✓ | 9 作者全列 ✓ |
| 2 | 10.1038/nmicrobiol.2016.48 | OK | ✓ | 17 作者 → et al. 正确保留 |
| 3 | （无 DOI） | — | — | DIMACS 技术报告，已声明无 DOI，可接受 |
| 4 | 10.1186/gb-2010-11-4-r37 | OK | ✓ | |
| 5 | 10.1126/science.1151532 | OK | ✓ | v2.9 修正后正确（原 1146308 为 404） |
| 6 | 10.1093/nar/gkad805 | OK | ✓ | Wang Z. … Zhang L.（10 作者），51(20), 10909–10923 ✓（原 gkad751 指向 TTD 库论文，已修正） |
| 7 | 10.1186/s13015-015-0032-x | OK | ✓ 标题 | **作者错误**：手稿 "Zielezinski, A. and Karlowski, W.M." → 应为 "Horwege, S. and Leimeister, C.-A."（见新发现问题 #1） |
| 8 | 10.1186/1471-2105-10-56 | OK | ✓ | Diaz 等 5 作者, 2009（原错误条目已修正） |
| 9 | 10.1093/bioinformatics/btg005 | OK | ✓ | |
| 10 | 10.1186/s13059-017-1319-7 | OK | ✓ | |
| 11 | 10.1093/bib/bbx067 | OK | ✓ | Bernard 等 8 作者, 20(2), 426–435（2019 为印刷年，Crossref online 2017，可接受） |
| 12 | 10.1093/bioinformatics/btu331 | OK | ✓ | 30(14), 2000–2008 ✓ |
| 13 | 10.1093/bioinformatics/btu177 | OK | ✓ | |
| 14 | 10.1186/s13059-019-1755-7 | OK | ✓ | 19 作者 → et al. 正确保留 |
| 15 | 10.1089/106652799318337 | OK | ✓ | |
| 16 | 10.1093/sysbio/syr010 | OK | ✓ | |
| 17 | 10.1093/molbev/msv150 | OK | ✓ | |
| 18 | 10.1093/oxfordjournals.molbev.a025808 | OK | ✓ | |
| 19 | 10.1093/molbev/msp098 | OK | ✓ | |
| 20 | 10.1371/journal.pone.0009490 | OK | ✓ | |
| 21 | 10.1093/molbev/mst010 | OK | ✓（截断） | 手稿省略副标题 "Improvements in Performance and Usability"（极轻） |
| 22 | 10.1093/bioinformatics/btz305 | OK | ✓ | |
| 23 | 10.1093/molbev/msaa015 | OK | ✓ | |
| 24 | 10.1186/s13059-016-0997-x | OK | ✓ | |
| 25 | 10.1093/bioinformatics/btu815 | OK | ✓ | Haubold/Klötzl/Pfaffelhuber, 31(8), 1169–1175 ✓（原 btv047 指向 KAPPA，已修正） |
| 26 | 10.1093/nar/gkt003 | OK | ✓ | Yi & Jin, 41(7), e75 ✓（原 gkt165 指向 Adams 杂交论文，已修正） |
| 27 | 10.1186/1471-2105-6-83 | OK | ✓ | Lunter 等 5 作者, 2005, 6:83 ✓（原错误 DOI 指向 Ooi 等，已修正） |
| 28 | 10.1093/molbev/msn275 | OK | ✓ | |
| 29 | 10.1186/s13059-019-1632-4 | OK | ✓ | |
| 30 | 10.1093/bioinformatics/18.3.440 | OK | ✓ | |
| 31 | 10.1073/pnas.0813249106 | OK | ✓ | |

**作者截断规则核验**：refs 2（17 作者）与 14（19 作者）正确保留 "et al."；refs 1/4/5/6/8/9/10/11/12/13/15/16/17/18/19/20/21/22/23/24/25/26/27/28/29/30/31 全列且均 ≤10 作者（与 Crossref 逐一比对通过）。**唯一例外是 ref 7（作者名单本身错误）**。

---

## 新发现问题（按严重度）

1. **【中等·需修复】Ref 7 作者名单错误**：手稿 L706 "Morgenstern, B., Zhu, B., Zielezinski, A. and Karlowski, W.M. (2015)"。Crossref（10.1186/s13015-015-0032-x）实际作者为 **Morgenstern, B., Zhu, B., Horwege, S. and Leimeister, C.-A.**（Zielezinski/Karlowski 是 ref 10 的作者，疑似复制错位）。v2.9 的 7 条参考文献修复未覆盖此条。一行修复。
2. **【轻微】Table 4 的 d 值与新增段落不一致**：L301 "Cohen's d (vs default spaced) = 0.54"（无符号）vs L377 "d=−0.53, p=0.006"。独立复算精确值 −0.529（应作 −0.53）。L301 的 0.54 是旧舍入且缺符号，与 stats m1 符号统一的精神相悖（虽 m1 严格针对 L1–L3 的 −0.45）。建议 L301 改为 −0.53。
3. **【极轻】异质性 z 值舍入偏差**：手稿 z=−2.71, p=0.007；用未取整的两组 d（+0.219 / −0.529）按 z=(d₁−d₂)/√(1/n₁+1/n₂) 精确复算为 **z=−2.82, p=0.0048**。手稿值疑似用取整输入（如 −0.50）计算。方向与显著性结论不变（均为 p<0.01 的真实异质性），可留可改；若改，同步更新 p。
4. **【轻微·遗留】Supplementary Notes S1/S2/S3/S6/S7 正文零引用**：正文仅引用 Note S4/S5/S8/S9/S10。v4 format 审稿只验证了 Figure/Table 的正文引用覆盖，未覆盖 Note。属 v2.8 遗留而非本轮退化，但投稿前宜补 5 处正文引用（如 Methods 仿真参数节引 S1、indel-rate 节引 S2、BAliBASE 节引 S3、统计框架节引 S6、SwissTree 节引 S7）。
5. **【轻微·遗留】phylo Minor 5 未处理**：Practical recommendations 仍无 boundary classifier 仅影响 n>1000 DCM 路径的说明，无 on/off 消融（审稿人允许二选一，两者均未做）。低严重度。
6. **【可选】字数声明略低**：声明 ~11,200，实测 11,380（不含表格行，1.6% 低估；changelog 自称 actual 11,376）。"approximately" 容差内，可接受；如求精确可改 ~11,400。
7. **【可选】Ref 21 标题截断**：MAFFT v7 完整标题含副标题 "Improvements in Performance and Usability"，手稿省略。该截断在文献中常见，非合规硬伤。

---

## 总体判定

**7 项 Major 全部 CLOSED（其中 A4/A5 经独立复算数值确认）；Minor 项除 phylo Minor 5（NOT-FIXED，低严重度遗留）与 phylo Minor 1（PARTIAL）外全部 CLOSED。**

**存在 1 项需修复的新发现问题：ref 7 作者名单错误（一行修复）**；另有 2 项建议顺手修复的舍入不一致（L301 d=0.54→−0.53；z=−2.71→−2.82）与若干可选润色。全局一致性（Abstract 216 词 ≤250、引文 [1]–[31] 首现严格单调且无缺漏、Table 1–12 首引严格单调、补充材料三组分块递增、旧值零残留）全部通过脚本验证。

完成上述一行参考文献修复后，v2.9 可视为已核销 v4 四位审稿人的全部实质意见。
