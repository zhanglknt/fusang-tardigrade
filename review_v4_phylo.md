# 第四轮盲审意见（系统发育学方法方向）— v4-phylo-reviewer

**稿件**: Fusang: Tardigrade Edition — K-mer Frequency Vector Alignment-Free Phylogenetic Inference Resilient to Indel-Rich Sequence Evolution（NAR 方法学论文，v2.8）
**评审历史**: 本审稿人此前完成 v2.4 初审（Major Revision, 62/100）与 v2.5 复审（Minor Revision, 76/100）
**本轮性质**: 完整重读全文盲审（非 diff 审查）

---

## 总体印象

v2.8 相对 v2.4/v2.5 发生了实质性的证据升级而非仅文字修补：L3 对比在 Linux 上完成全量 30 种子并主动披露了此前 Windows 运行的协议错误（漏 `-nt` 标志，蛋白质模式跑 DNA 输入）；kmacs 与 Skmer 补测使经典竞品对照基本齐备，且 Table 7 在统一 seed set 上重算消除了此前三组 seed 混用的嵌合体问题；新增的 AFproject 真实基因组基准（Table 12）既是本文最强的正面证据（fish mtDNA 上 multi-k 与 MSA+ML 及最优已发表 AF 方法打平），也诚实地划出了方法边界（E. coli/Shigella 菌株级数据上 k-mer cosine 停在 0.50）。multi-k ensemble 在 >2,000 倍序列长度范围内自动救回饱和的 k=5 默认参数，是全文最扎实的实用性贡献。剩余问题集中在措辞层面：摘要中出现两处未带限定的"significant/matches"声称，与正文自报的 borderline 校正后 p 值和 n=25 小样本饱和基准不完全对齐；Introduction 的 spaced k-mer 中心化表述仍与 Discussion 的"cosine 为主驱动"定位存在残留张力。

---

## 逐项打分

| 项目 | 分数 (1–10) | 理由 |
|---|---|---|
| 新颖性 | 5 | 核心组件（spaced k-mer、cosine、NJ、DCM、multi-k）仍是既有技术的组合；multi-k 在跨长度尺度上的自动适应（Table 12 发现 2）是本文最具原创性的实用贡献，但整体仍属系统性评测型而非新方法型贡献。 |
| 方法严谨性 | 8 | 统一 seed set 重算、L3 协议错误的主动披露与作废、配对 Cohen's d 的计算说明（L367）、Bonferroni 家族的预指定（L144-145）均达到高标准；扣分项：fish 基准未讨论饱和效应，spaced 鲁棒性仍无理论推导。 |
| 验证充分性 | 8 | L3 30 种子完整、kmacs k 扫描 + Skmer 双 k 值结构不适用测试、两个真实全基因组基准 + SwissTree gap 率追溯性定量（中位 50%）；剩余缺口：真实基因长度（500–1,000 bp）indel-rich DNA 数据仍未验证（作者自己在 L545 承认）、Co-phylog DNA 仍仅默认参数。 |
| 声称-证据匹配 | 8 | E. coli 边界的处理堪称诚实报告的范本（L498、L572）；但摘要两处声称强度超出证据：`significantly outperforms MAFFT+FastTree2 (p=0.024)` 未提校正后 p_adj=0.071 borderline；fish `matches MSA+ML and the best published AF methods` 未提 n=25/2-split 粒度/多方法同分意味着基准饱和。 |
| 竞品覆盖 | 8 | kmacs（k 扫描）、Skmer（结构性不适用已实测论证）、andi、Co-phylog 四家经典方法齐备，AFproject 发表配置的交叉验证（复现值与已发表值差 ≤1 split）是加分项；剩余：FFP 与 spaced-word frequency 家族（FSWM/SpaM）仍未测（L561 已承认），Co-phylog DNA 无参数扫描。 |

---

## 相对上轮评审的变化

**已解决（相对 v2.4 初审与 v2.5 复审）**：
1. **L3 仅 n=5（初审 M6，复审 CLOSED 待补）→ 彻底解决**：Linux 全量 30 种子（Table 5, L307-328），统计完整（p=0.024, d=−0.45, 19/30 wins），且主动披露并作废了旧 Windows 运行的双重协议错误（MAFFT 空比对 + 漏 `-nt`，L551, Supp S8）。这一处理方式符合最佳实践。
2. **竞品缺口（初审 M5，复审 PARTIAL）→ 基本解决**：kmacs 以 k∈{3,5,10} 扫描补测（Table 7, L360, L373）；Skmer 在 k=31/k=21 双配置、27 种子全部 coverage 除零，结构性不适用结论有实测支撑（L379）；Table 7 统一 seed set 100–126 消除了旧版三组 seed 混用问题；Table 12 中 kmacs/Mash 以 AFproject 发表配置复现（差 ≤1 split）交叉验证了管线。
3. **Table 3 参考系标签矛盾（复审新发现）→ 已修复**：现两列均明确标注 TRUE-relative（L257-266）。
4. **16S (k,gap) 参数矛盾（复审 M9）→ 已修复**：Methods（L177）、Results（L406）、Table 9（L412）统一为 k=5,gap2。
5. **Table 1 Winner 列 n.s. 判负（复审新发现 2）→ 已修复**：现为 "Tie (n.s.)"（L216-217）。
6. **SIMPLE_THRESHOLD=1000 无依据（初审 M7）→ 以诚实披露方式解决**：L543 明确承认"n=500–1000 过渡区未做简化 vs DCM 直接精度对照"。仍无数据，但已从无依据断言降级为声明的工程默认。
7. **参考系/协议不可比标注（初审 M3）→ 持续保持**：Table 1/2/4/5 均标注协议独立性；Table 12 进一步标注 ETE3 nRF 与内部 nRF 不可比（L181, L488）。
8. **真实数据验证（初审最大缺口之一）→ 大幅改善**：Table 12 + S16 新增 fish mtDNA 与 E. coli/Shigella 两个带可信参考树的社区基准，正反两面结果均如实报告。

**仍在（遗留）**：
1. Introduction L42 仍将"the combination of **spaced** k-mer representation with cosine distance"表述为有效信号来源、加粗"spaced k-mer frequency vectors have not been systematically evaluated"——与 Discussion L512"cosine is the primary accuracy driver; spaced patterns provide domain-dependent secondary benefits"的定位仍有张力（复审 C1 残留，未动）。
2. spaced k-mer indel 鲁棒性仍无理论推导（初审 M4）：L42"theoretical motivation"、L527"inherent robustness from spaced sampling"仍为断言。随 spaced≈contiguous 的全面实证（Table 7 p=0.26；Table 10 p=0.31），此问题严重度已低，但既然 L527 仍宣称"inherent robustness"，应给出期望值分析或改写为假设。
3. Boundary classifier 仍无 on/off 因果消融（初审 M8）：L557 仅新增"feature ablation remains future work"的承认。
4. 蛋白最优 spaced 构型 k=4,gap1（0.239）vs contiguous k=5（0.244）仍未单独检验显著性（L512 引用为 spaced"clearest benefits"之一，证据强度弱）。

**新增问题（v2.8 引入）**：见下节 Major 1、2 与 Minor 1。

---

## 必须修改（Major）

1. **摘要 L26 "significantly outperforms MAFFT+FastTree2 under indel-rich coalescent simulation (…p=0.024)" 措辞过强**。稿件自身统计框架（L144-145）将 pipeline-level 比较归为 exploratory，正文（L324）与 Limitations（L551）均如实报告 Bonferroni 校正后 p_adj=0.071 borderline，L1 胜 19/30（63.3%）也非压倒性。摘要应改为"outperforms MAFFT+FastTree2 (p=0.024, borderline after multiple-comparison correction)"或等效限定——这是全文最高可见度位置与正文证据强度的不一致，必须对齐。
2. **摘要 L26 与正文 L494 的 fish "matches MSA+ML and the best published alignment-free methods" 声称需加饱和度限定**。n=25（44 informative splits）下 nRF=0.045 仅对应 2-split 差异；multi-k、kmacs、Mash、Co-phylog、MAFFT+FastTree2 五种方法同分 0.045，k=9/k=11 单 k 也达 0.091——基准明显处于多方法同达下限的饱和区，其判别力是粗粒度的。声称字面为真，但在摘要不加 n=25 与"benchmark saturation"限定，读者会误读为精细判别下的打平。建议摘要补 "(25 taxa; multiple methods reach the same floor)"，正文 L494 补一句饱和效应讨论。
3. **Introduction L42 的 spaced 中心化残留表述需与全文定位对齐**（两轮遗留）：将加粗句改为"k-mer frequency vectors (spaced and contiguous) have not been systematically evaluated…"，并把贡献点 (1) 中"the combination of spaced k-mer representation with cosine distance"改为"k-mer frequency vectors with cosine distance"，使 Introduction 与 Discussion L512 的诚实定位一致。

---

## 建议修改（Minor）

1. **Table 12 fish 行的 Co-phylog（Python 重实现）0.045 优于其已发表值 0.09**（L484, L488）：稿件已注明重实现差异，但建议明确说明为何重实现会更好（NJ 实现/参数差异），否则读者可能怀疑对标配置不完全等价；DNA Table 7 中的 Co-phylog 同为重实现，建议在 Table 7 脚注交叉引用该说明。
2. **Co-phylog DNA 仍仅 k=19 默认配置**（L561 已承认）：鉴于 fish 上 Co-phylog AFproject 配置（halfctx=5,k=11）表现良好，DNA 基准至少应补 halfctx=5 一档，或将 L561"most robust alignment-free approach"的总结性措辞弱化为"among the most robust"。
3. **spaced 鲁棒性理论推导（遗留）**：给出单 indel 事件下 spaced vs contiguous k-mer 期望破坏数量的简单分析（可引用 [7][13] 的 spaced-word 理论框架），或将 L42"theoretical motivation"与 L527"inherent robustness"改为"hypothesized"。
4. **Table 5 协议从（旧版未标注 L）变为 L=1000 bp**（L307）：请确认 Supp Note S8 明确记录该参数（以及是否为协议修正的一部分），避免新旧运行参数差异引起审计疑问。
5. **Boundary classifier 因果消融（遗留）**：仍建议补一个使用/不使用分类器时 DCM 管线最终 nRF 的 on/off 对照；若不做，建议在 Practical recommendations 中明确该组件仅影响 n>1000 的 DCM 路径。
6. **FFP 与 spaced-word frequency（FSWM/SpaM）家族仍未测试**（L561 已承认）：FSWM 与本文方法在表示层面（spaced-word 频率）最同源，建议至少在 Discussion 说明未测原因（如实现可得性），为下一版或回应信预留立场。
7. **E. coli 数据集上自家未运行 andi/Co-phylog**（L486 引用已发表值 0.08、kmacs 超时 L490 ³）：建议在 Table 12 脚注说明"published values"与"our run"的区别已尽可能标注，避免与 fish 行混读；目前标注已较好，仅需微调。

---

## 总体推荐

**Minor Revision**

理由：v2.8 在证据层面补齐了我前两轮指出的全部关键实验缺口（L3 30 种子、kmacs/Skmer、真实基因组基准），并以罕见的透明度处理了自身的协议错误与不利结果（E. coli 边界）。剩余问题集中在三处措辞级修改（摘要两处声称限定 + Introduction spaced 表述对齐）与若干建议级补强，无需新实验。完成 Major 1–3 后，本文将达到 NAR 方法学论文的可接收标准。

**投稿就绪度：84/100**（v2.4: 62 → v2.5: 76 → v2.8: 84）
