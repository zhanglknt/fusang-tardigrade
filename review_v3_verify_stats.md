# Fusang v2.5 统计复审核销报告

**复审人**：v3-stats-reviewer
**原评审**：review_v3_stats.md（2 Critical + 4 Major + 6 Minor）
**复审范围**：逐项核销 + 独立验证 + 新发现问题

**独立验证说明**：针对 team-lead 提供的背景①，我未照单全收，直接核查了源脚本 `run_indel_benchmark_v9.py:135-160`——`calc_nrf()` 对 method="ft2" 与 method="fusang" **均以 `seed{seed}_indel_true.nwk` 作为参照树**（第 137、144 行），证实 130-seed 基准两列确为 TRUE-relative，原"FT2-relative"确属标签误标。关键结论的统计有效性由此得到源码级确认。

---

## 逐项核销

### Critical 1：跨参考系比较与配对检验参照树不明 → **CLOSED**

- 源码验证：两方法 nRF 均对 TRUE 树计算（见上），配对 Wilcoxon 前提成立。
- 文本修复到位：Table 1 标题（L207）"All nRF values are TRUE-relative"；Table 3 标题（L252）同；L261 明确"Both columns are TRUE-relative and directly comparable per seed"；L219 明确"paired Wilcoxon signed-rank test on TRUE-relative nRF"；n=1000 indel 行以 † 脚注单独标注 FT2-relative 且声明不可比（L219、L373）。
- 摘要（L26）与贡献段（L50）均改为"against the simulated ground-truth tree"表述。
- "13.3% 优势"公式在两列同框架后数学上成立。
- **残留**：Table 3 列头（L254）仍印 "Fusang nRF ↓ (FT2-rel)"，与自身表题直接矛盾——属修复遗漏的标签错误（计入新发现 Minor-N1，不影响核销结论）。

### Critical 2：Mash "collapse vs modest" 叙事与 Table 2 矛盾 → **CLOSED**

- 摘要（L26）："both MinHash (Mash) and k-mer cosine distances degrade severely... (1.93× vs 1.97×), though cosine distances retain marginally lower absolute error"——对比如实。
- 正文 L244 明确："a formal between-method test on matched seeds was not performed, this small absolute difference should be interpreted descriptively"——正是我要求的"未做检验则声明描述性"。
- Discussion（L471、473）与 Methods（L510）口径一致；"collapse" 措辞全文已删除（grep 无命中）。
- L510 还主动披露并作废了开发期单 seed 错误值（nRF=1.005，解析错误），透明度加分。

### Major 1：112-seed 基准 SD 三处不一致 → **PARTIAL**

- FT2 已全稿统一为 **0.084±0.019**：摘要 L26、L50、L219、L250、L269 五处一致 ✓（原 0.085±0.025 已消除）。
- **残留**：L273 仍为 "Fusang variance (σ=0.017)"，与全稿统一的 ±0.016（L26、L50、L269）矛盾。一行之改，遗留至今。

### Major 2：多重比较校正覆盖不一致 → **PARTIAL**

- 已修复部分：Table 4（原 Table 7）双 p 值问题改为规范表述——"Wilcoxon signed-rank test (vs default spaced, **pre-specified primary test**): p=0.006 (paired t-test as **sensitivity analysis**: p=0.007)"（L295），主检验+敏感性框架正确，选择报告嫌疑消除 ✓。
- 未修复部分：校正家族仍仅限"5 个 ground truth 数据集"（L144-145 未变）；SwissTree 两个 Co-phylog 比较（p=0.005/0.014、p=0.001/0.006）、Table 5 L1 vs L0（p<0.0001）、Table 7 竞争者比较（p<0.001、p=0.0002）等 confirmatory 色彩的检验仍未纳入任何校正，也未统一标注 exploratory。
- 考虑到各未校正检验 p 值与临界值距离较远（校正后结论不会改变），严重度由 Major 降为 Minor 量级，但形式上仍应补一句统一声明。

### Major 3：SwissTree d=1.08 隐含强负相关 → **CLOSED**

- L429 新增逐家族数据披露：Co-phylog halfctx=5 胜出的 3 个家族（ST008/ST009/ST012）及具体 nRF 值，全量数据指向 S11–S12。
- 披露数据证实差值确为**双向**（3 个家族 Co-phylog 优、其余 Fusang 大幅优），SD_diff=0.180 源于真实的跨家族异质性而非数据对齐错误——我要求的"展示逐家族差值以排除转置错误"已满足。
- 旁证：halfctx=5 的 Wilcoxon p 由 v2.4 的 0.006 更正为 0.014（t 检验 p=0.005 不变），表明作者确实重算了该比较；摘要引用的 1.5× 对应 halfctx=11 的 Wilcoxon p=0.006，内部一致 ✓。

### Major 4：功效声明 ~55% 低估 → **CLOSED**

- L147 已改为 "30-seed benchmarks provide moderate power (**~75%** for d=0.5)"，与我的复核值（75–78%）一致 ✓。

### Minor 1：Table 7 改进幅度与 indel 扫描百分比算术 → **NOT-FIXED**

- L293 仍为 "0.008 (6.7% relative reduction)"——实际 0.112−0.105=**0.007**（6.25%），0.008 对应 7.1%，两个数字仍不兼容。
- L250/L261/L263 仍为 "13.3%"（复核 0.010/0.076=**13.2%**）与 "4.6–4.7%"（复核 4.46%/4.76%，应为 4.5–4.8%）。
- 现两列同框架后这些百分比有了合法意义，但算术舍入错误原样保留。

### Minor 2：cosine vs JSD 单尾 n=10 预实验 → **NOT-FIXED**

- L73 原文未动："10 seeds, preliminary... Wilcoxon p=0.031, **one-tailed**"。既未补双侧正式检验，也未改为启发式选择的声明。措辞上"modest but consistent advantage"较前有收敛，但方法论缺陷依旧。

### Minor 3：27 seeds 缺失原因 → **NOT-FIXED**

- L348 仍为 "27 seeds with valid reference trees"，3 个种子的缺失机制（是否 MCAR/MAR）全文无交代。

### Minor 4：16S 置换检验细节 → **NOT-FIXED**

- L176 仍仅 "using permutation tests"；Table 9（L402-408）无置换次数、统计量定义、单/双边任何说明。8/12 vs 10/12 恢复率对比仍未加注不可过度解读的提示。

### Minor 5：Table 2 注跨框架比较 → **CLOSED**

- 原误导句（0.394 vs 0.112 跨框架对比）已删除；Table 2 注（L236）改为同框架说明并主动披露作废旧单 seed 值。L242/L473 的比较均为同框架（0.394 vs 0.376 clean TRUE-relative）✓。

### Minor 6：分类器 100% vs CV AUC=0.84 张力 → **NOT-FIXED**

- L344 保留同一模拟引擎的既有声明，但未新增：测试表现远超 CV 表现的讨论、按场景难度分层的准确率、测试参数范围与训练空间边缘的覆盖关系。

---

## 新发现问题（v2.5 引入或修复过程中暴露）

- **Minor-N1（建议升 Major-级校对优先级）**：Table 3 列头 L254 残留 "(FT2-rel)"，与表题（L252）、正文（L261）直接矛盾。这是 Critical 1 修复的"最后一块拼图"——读者看到列头仍会误以为跨框架。**必须改。**
- **Minor-N2**：Table 8（L386）indel 行数值为 30-seed（0.077±0.018），Notes 却挂"112 seeds: p=0.052"，两个不同基准混在一格，建议拆注。
- **Minor-N3**：SwissTree halfctx=5 的 Wilcoxon p 由 0.006 改为 0.014 但全文无变更说明。虽不强制披露修订历史，但建议作者留存重算记录备查。

---

## 修订后逐项打分

| 维度 | v2.4 | v2.5 | 理由 |
|---|:---:|:---:|---|
| 统计方法正确性 | 6 | **8** | 核心参照系问题解决并源码验证；残余为标签遗留（N1）与未扩展的校正家族 |
| 报告完整性 | 6 | **7** | FT2 SD 五处统一；残留 Fusang σ=0.017 一处、27 seeds 未解释、置换检验细节缺 |
| 样本量与功效 | 7 | **9** | 功效声明已更正为 ~75%；n=5 处理依旧堪称范例 |
| 效应量解释 | 6 | **9** | Mash 叙事彻底修正并声明描述性；SwissTree 逐家族披露到位；配对 d 声明完善 |
| 抗 p-hacking 稳健性 | 6 | **7** | 主检验+敏感性框架规范；JSD 单尾预实验与校正覆盖仍未动 |

## 修订后推荐意见

**Minor Revision**

两个 Critical 均已实质性关闭（其一经源码独立验证），6 项 Major/Minor 修复到位。剩余全部为文字、标签与算术层面：1 处矛盾列头（N1）、1 处 SD 残留（σ=0.017）、3 处百分比舍入、若干方法学细节补注。无需任何新实验或重算，一轮编辑性修改即可达标。

## 修订后投稿就绪度评分

**82 / 100**（v2.4 为 62）

- 统计方法与执行：27/30（核心问题解决；扣校正家族未扩展与 JSD 单尾）
- 报告完整性与一致性：18/25（扣 N1 列头、σ=0.017 残留、27 seeds、置换细节）
- 结论与数据匹配度：19/20（Mash 叙事已匹配；扣百分比舍入）
- 透明度与可复现性：18/25（逐家族披露、作废值声明、主/敏感性框架加分；扣预实验单尾与校正声明缺失）
