# 复审核销报告（系统发育学方法方向）— v3-phylo-reviewer

**稿件版本**: v2.5（NAR_MANUSCRIPT_REVISED.md，725 行）
**对照评审**: review_v3_phylo.md（v2.4）
**核销结论图例**: CLOSED = 已解决 / PARTIAL = 部分解决 / NOT-FIXED = 未解决

---

## 一、Critical 问题核销

### C1. spaced k-mer 卖点与证据矛盾 — **PARTIAL（约 85% 解决）**

已修复的证据：
- Running Title 已改为 "K-mer frequency vector phylogenetics"（L5），标题（L1）不再含 "spaced"
- 摘要（L26）改为"systematically evaluates k-mer frequency vector cosine distances — across spaced and contiguous patterns"，不再以 spaced 为卖点
- Discussion 第 1 点（L469）明确写道"**Cosine distance is the primary accuracy driver; spaced patterns provide domain-dependent secondary benefits**"，并引用具体 regime 证据（蛋白 k=4,gap1 最优 nRF=0.239；indel 下 gap3 优于 gap2 10.5%）——这正是我上次建议的"领域依赖次级收益"定位
- SwissTree 结果段（L435）明确"Spaced k-mers show no significant advantage on protein data (p=0.31)"

残留问题：
- Introduction 贡献段（L42）仍写"The core contribution of this work is therefore a systematic evaluation of **spaced k-mer frequency vector cosine distances**"，且贡献点 (1) 仍表述为"the combination of **spaced** k-mer representation with cosine distance provides effective phylogenetic signal"——与 Discussion 的"cosine 为主驱动"定位仍有轻微张力。建议将 L42 加粗部分改为"k-mer frequency vector cosine distances (spaced and contiguous)"。
- 蛋白最优构型 k=4,gap1（0.239）vs contiguous k=5（0.244）差异极小（0.005），文中以"best protein-family configuration is spaced"作为 spaced 有效区证据（L469），但该差异未单独做显著性检验（p=0.31 是 spaced vs contiguous 的汇总检验）。作为"次级收益"定性可以接受，但措辞从"clearest benefits"降为"slightly better"更稳妥。

### C2. MinHash collapse 叙事与数字不符 — **CLOSED**

- 摘要（L26）已改为诚实版："both MinHash (Mash) and k-mer cosine distances degrade severely relative to their clean baselines (1.93× vs 1.97×), though cosine distances retain marginally lower absolute error"
- 正文（L244）明确承认"because a formal between-method test on matched seeds was not performed, this small absolute difference should be interpreted descriptively"——虽未补做成对检验（我原建议第 2 条），但显式披露检验缺失并降级为描述性表述，符合诚实报告标准
- Discussion 第 2 点（L471）与之一致；L246/L528 的"not recommended / should be avoided"表述基于基因长度场景的绝对错误率（~76% bipartition error），合理
- 旧版"collapse vs modest degradation"框架已彻底移除

---

## 二、Major 问题核销

### M3. 参考框架不统一 — **PARTIAL**

已修复：
- Table 1（L207, L218-219）标题与脚注明确"All nRF values are TRUE-relative except where noted"，n=1000 indel 行以 † 单独标注 FT2-relative 且声明不可互比——原"FT2-rel"误标已修正
- Table 2（L225）与 Table 5（L304）均新增说明：coalescent 模拟协议独立于 Table 1/3 的 guide-tree 协议，"absolute nRF values are higher for all methods under this protocol and should not be compared across benchmarks"——协议差异已显式标注，这正是我上次要求的核心
- Table 7（L350）标注 FT2-relative

**新发现问题（残留矛盾）**：Table 3 存在内部自相矛盾——表头（L252）声明"All nRF values are TRUE-relative"，正文（L250, L261）也说"Both columns are TRUE-relative and directly comparable per seed"，**但表格列标签（L254）仍写着"Fusang nRF ↓ (FT2-rel)"**。team-lead 说明称"Table 1/3 原标签误标已修正"，但 Table 3 的列标签实际未修正。这是 v2.5 遗留的一处直接矛盾，投稿前必须修正（将列标签 "(FT2-rel)" 改为 "(TRUE-rel)"）。

### M4. spaced k-mer indel 鲁棒性缺理论推导 — **NOT-FIXED**

- 主文仍只有直观论证：L42"spaced k-mers provide theoretical robustness at high indel rates"（"theoretical"一词无推导支撑）、L392 gap3 优势归因于"wider spacing better tolerating indel-induced length variation"、L484"inherent robustness from spaced sampling"
- 未新增单 indel 事件下 spaced vs contiguous 被破坏 k-mer 数量的期望分析，也未引用 Leimeister & Morgenstern 的 spaced-word 理论框架做形式化论证（文献 [7][13] 已在引文列表中，具备引用条件）
- Supplementary Note S2 标题为"Gap optimality on indel data: mechanism and robustness analysis"，可能含机制分析，但主文未引用其具体结论；从稿件本身无法确认其是否包含形式化推导
- 缓解因素：卖点已从"spaced 鲁棒"降级为"次级收益"，该问题的严重度随之从 Major 降为 Minor；但既然 L42 仍保留"theoretical robustness"字样，建议要么给出推导、要么改为"hypothesized robustness"

### M5. 竞品缺口（Skmer/kmacs/Alfpy）— **PARTIAL**

- Limitations（L518）现以真实文献明确承认缺口：Skmer [30]（真实引用，Genome Biol. 2019）、FFP [31]（PNAS 2009）、spaced-word 方法族 [13]、kmacs [12]（Bioinformatics 2014），并指明"particularly Skmer and kmacs under indel-rich conditions, represents important future work"
- Introduction（L42）也已将 kmacs、spaced-word frequency methods、Alfpy/AFproject 纳入现状综述，"首次系统评测"类表述已改为"have not been systematically evaluated... under realistic evolutionary conditions (indel-rich)"，更为准确
- 未做的部分：实际 benchmark 仍未补充（Skmer/kmacs 未跑）；Co-phylog 参数扫描仍未做（L518 仅保留承认性文字）。作为修订轮的妥协，"明确承认 + 真实引用 + 降级声称"可以接受，但严重度维持 Minor-Major 边界：若审稿人坚持"与最直接可比方法的对照是方法学论文的必要条件"，此点仍可能成为拒稿理由

### M6. L1 vs L3 仅 n=5 — **CLOSED（框架层面）**

- 摘要（L26）已完全移除 L1≈L3 等价性表述，不再出现在 abstract 层面
- Table 5 引言（L302）以加粗声明"this is a preliminary result... interpret as encouraging preliminary evidence rather than definitive equivalence"
- 结果段（L319, L323）两处重申 n=5 限制、p=0.24 不构成等价检验、"initial evidence warranting larger-scale validation"
- 残留轻微问题：贡献点 1（L48）仍以该 n=5 结果开篇且未在贡献点内标注"preliminary"（正文其他处均有标注），建议在 L48 末尾加"(preliminary, n=5)"。Linux 30 种子验证仍未完成（L508, L520 列为 future work），但披露充分，属可接受的待办而非隐瞒

---

## 三、Minor 问题核销

### M7. SIMPLE_THRESHOLD=1000 缺精度依据 — **NOT-FIXED**
L105/L199 仍断言"For n≤1000, the simplified pipeline matches DCM+EPA accuracy"，但消融证据仍只来自 n=200（seed=42 + 10 种子）；n=500–1000 区间仍无 DCM vs simplified 直接精度对照。建议将措辞改为"based on ablation at n=200, we adopt SIMPLE_THRESHOLD=1000 as an engineering default; direct DCM-vs-simplified accuracy comparison at n=500–1000 remains future work"。

### M8. Boundary classifier 与建树精度无因果消融 — **NOT-FIXED**
仍无"使用/不使用分类器时最终树 nRF 差异"的消融实验；分类器验证（L327-344）仍是独立组件评估。缓解：L344 对同一模拟引擎的局限披露充分，88 场景 Wilson CI 统计表述规范。维持 Minor，建议降级为 Supplementary 工程组件或补一个 on/off 消融。

### M9. 16S 验证 (k,gap) 参数标注自相矛盾 — **NOT-FIXED**
矛盾仍然存在且可精确定位：Methods（L176）写"k=5,**gap2**, simplified pipeline, FastME"，而 Results（L396）与 Table 9 表题（L402）均写"k=5,**gap1**"。两处必有一处错误，投稿前必须统一。

---

## 四、新发现问题（v2.5 引入或首次定位）

1. **[Minor，必须修] Table 3 列标签自相矛盾**（详见 M3）：L254 列头"(FT2-rel)" vs L252/L261"Both columns are TRUE-relative"。
2. **[Minor] Table 1 "Winner" 列将不显著结果标注为"Fusang (n.s.)"**（L212；Table 8 L386 同）：n=200 indel 行 0.077 vs 0.080 标 "Fusang (n.s.)"。虽然带了 n.s. 标注，但在 Winner 列给不显著结果判胜负有轻度误导，建议改为"Tie (n.s.)"或"Fusang numerically lower (n.s.)"。
3. **[Minor] 贡献点 1（L48）未标 preliminary**（详见 M6）。
4. **[观察项，非问题] Table 5 中 L1 vs L0 的 Cohen's d=3.55**（L321）已附机制解释（高基线+低方差导致），表述合理。

---

## 五、修订后逐项打分（1–10）

| 项目 | v2.4 → v2.5 | 理由 |
|---|---|---|
| 新颖性 | 5 → 5 | 重新定位不改变贡献本质：仍是"系统评测 + multi-k ensemble + 自适应管线"的增量式组合贡献。 |
| 方法严谨性 | 7 → 7 | 统计框架不变；Table 3 列标签矛盾与 16S 参数矛盾属执行层面疏漏，未影响方法本身。 |
| 验证充分性 | 6 → 6 | L3 仍 n=5、Skmer/kmacs 仍未实测、阈值 1000 仍无 n=500–1000 精度对照；但协议差异标注后，已有验证的可解释性显著提升。 |
| 声称-证据匹配 | 8 → 9 | 两个 Critical 叙事实质性修复：spaced 降级为次级收益、Mash 改为双严重退化+描述性差异；摘要与正文基本一致。扣 1 分：Table 3 列标签、16S (k,gap)、L42 残留 spaced 中心表述、Winner 列判负。 |
| 与现有文献的关系 | 6 → 7 | Skmer/kmacs/FFP/spaced-word 现均以真实文献引用并在 Introduction 正确定位；"未系统评测"表述更准确。未实测直接对照仍扣 3 分。 |

---

## 六、复审推荐意见

**Minor Revision**

两个 Critical 问题已实质性解决（C2 完全关闭；C1 在 Running Title/摘要/Discussion 层面完成重构，仅 Introduction L42 有残留张力）。剩余必须修复项均为一行级文字修正（Table 3 列标签、16S gap1/gap2、L48 加 preliminary、Table 1 Winner 列），不涉及新实验。坚持性保留意见（不阻塞接收）：Skmer/kmacs 实测对照与 L3 Linux 30 种子验证应在接收前或最终版中完成，若编辑要求"与最直接可比方法的对照"为必要条件，则维持 Major Revision。

## 七、修订后投稿就绪度评分

**76 / 100**（v2.4 为 62）

加分：Mash 叙事修复（+6）、spaced 重新定位（+5）、参考系协议标注（+4）、文献补强（+2）。
剩余扣分：竞品未实测（-8）、L3 仅 n=5（-6）、理论推导缺失（-4）、四处一行级内部矛盾（-6）。
