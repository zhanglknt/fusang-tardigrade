# NAR 格式复审报告（核销单）— v3-format-reviewer

**稿件**: NAR_MANUSCRIPT_REVISED.md (v2.5)
**对照基线**: review_v3_format.md（v2.4 评审，4 Critical + 6 Major + 5 Minor）
**复审日期**: 2026-10-09
**方法**: 全部定量指标脚本化重测（引用序列、首引行号、正文提及计数、词数/字符数），并对关键段落人工核对。

---

## 总体结论

v2.5 修复了 v2.4 评审中**最致命的三类问题**：参考文献引用顺序制、Figure 3–6 悬孤、3 处引文错配，且 Abstract 已改为单段非结构式（实测 193 词）。但**表格编号体系仍未按首引顺序排列（Critical #2 未核销），且 Table 6 正文零引用（新退化）；4 项补充材料（Supp Fig S6、Supp Table S7/S11/S12）仍悬孤**。修订方向正确、完成度约 75%，剩余问题全部为机械性修复。

---

## 逐项核销

### Critical

| # | 原问题 | 裁定 | 证据 |
|---|--------|------|------|
| C1 | 参考文献未按引用顺序编号 | **CLOSED** | 脚本重测：正文首现序列 = 1→2→3→…→31 完全单调（唯一异常 "0" 系 L68 数学式 `[0,1]^(4^k)` 误检）；31 条文献全部被正文引用，最大引号 31 与列表条数一致，无超界 |
| C2 | 表格编号与首引顺序脱节 | **NOT-FIXED** | 重测首引序：Table 7 & 10（L168，Methods 前瞻引用）→ Table 1（L205）→ Table 4（L221）→ Table 3（L225）→ Table 2（L242）→ Table 5（L261）→ Table 8（L377）→ Table 9（L396）→ Table 11（L441）。表确实被重编过号（scalability 现为 Table 11、boundary classifier 现为 Table 6），但首引顺序仍非 1→11 单调 |
| C3 | Figure 3–6 正文零引用 | **CLOSED** | 重测 Figure 1–6 正文提及各 1 次；Figure 3 已补于消融段（L188）、Figure 4/5/6 均已在相应 Results 小节出现 |
| C4 | refs 22/23（Wong 2008、Zielezinski 2017）未被引用 | **CLOSED** | Wong 2008 现为 ref 5（L38 处 [3,5] 区段内被引）、Zielezinski 2017 现为 ref 10（L42 区段被引）；脚本确认 uncited = ∅ |

### Major

| # | 原问题 | 裁定 | 证据 |
|---|--------|------|------|
| M5 | Supplementary 材料大面积悬孤 | **PARTIAL** | 已补：Supp Fig S1–S5、S7；Supp Table S1–S6、S8–S10、S13–S15（正文均有 ≥1 次引用，如 S15→L327、S10→L348、S3→L188）。**仍悬孤 ×4**：Supp Fig S6（仅出现于列表 L547）、Supp Table S7（L563）、S11（L583）、S12（L587）正文零引用 |
| M6 | 引文-文献错配 ×3 | **CLOSED** | ① kmacs→ref 12，条目正确（Leimeister & Morgenstern 2014, *Bioinformatics*, 30, 2000–2008, btu331，L667）；② bioNJ→[18] Gascuel 1997 *MBE* 14:685（L109 引用，L679 条目正确）；③ "MashTree" 已改写为 "Mash v2.3 [24]"（L165），工具-文献一致。新增 Gascuel(18)/spaced-word(13)/Skmer(30)/FFP(31) 四条均已被正文引用 |
| M7 | 无编号 author-year 行内引用 | **CLOSED（残留 2 处，降级为新 Minor）** | SpaMz→[13]（spaced-word 文献，L42 区段）、Skmer→[30]、FFP→[31]（L518 数字化引用）；kr/SpaMz 不可考提及已移除。**残留**：L369 "(Haubold et al. 2015, Bioinformatics)"、L416 "(Zielezinski et al. 2019, *Genome Biology*)" 两处冗余括注，与已编号的 [25]/[14] 重复 |
| M8 | 参考文献格式内部不一致 | **PARTIAL** | 期刊缩写已统一（Genome Biol. 全部统一，L651/663/671/691/703 ✅；DOI 前缀统一 "DOI:"）。**作者截断仍两套并存**：首作者+et al.（refs 6, 8, 11, 15, 16, 17, 20, 22, 23, 25, 27, 29）vs 全列/多作者+et al.（refs 1, 2, 5, 7, 10, 12, 13, 14, 24, 26, 28, 30, 31）。NAR 惯例为列全作者（≤10）后接 et al. |
| M9 | Abstract 结构式三段 | **CLOSED** | 实测：单段、无 Background/Results/Conclusion 标签，193 词（≤250 ✅；与声称 198 词差 5 词，系计数口径差异，可忽略） |

### Minor

| # | 原问题 | 裁定 | 证据 |
|---|--------|------|------|
| m10 | 文末字数声明不实 | **CLOSED** | 声明改为 "approximately 9,800 words"；实测正文不含表格 9,551 词、含表格约 10.3k，声明落在合理区间 |
| m11 | 全大写强调 | **PARTIAL** | NOT 由 4 降至 2（L329 等）；IMPORTANT 已清除；**TRUE 反由 27 升至 31**（"TRUE tree" 仍未按建议改为小写定义式写法） |
| m12 | 防御性措辞密度 | **NOT-FIXED** | preliminary 7 次、borderline 7 次，与 v2.4 持平或略增；重复免责段仍在 |
| m13 | 标题大小写不统一 | **CLOSED** | Running Title 已改为 "K-mer frequency vector phylogenetics"（36 字符 ≤50 ✅），小节标题风格趋于一致 |
| m14 | Zenodo DOI 占位符 | **CLOSED** | 真实 DOI 10.5281/zenodo.20746742 已填入 Data Availability |

---

## 新发现的问题（v2.5 引入或首次检出）

| 严重度 | 问题 | 证据 |
|--------|------|------|
| **Major** | **Table 6（boundary classifier 验证表）正文零引用** —— 系重编号后的新退化：boundary classifier 一节（L325–344）通篇无 "Table 6" 字样，该表仅有 caption（L334） | 脚本：Table 6 在 body 中仅出现于 L334 caption 行 |
| Minor | 冗余 author-year 括注 ×2 | L369 "(Haubold et al. 2015, Bioinformatics)"；L416 "(Zielezinski et al. 2019, *Genome Biology*)" —— 应删去括注仅保留 [25]/[14] |

---

## 修订后逐项打分（1-10）

| 维度 | v2.4 | v2.5 | 理由 |
|------|:---:|:---:|------|
| 硬性合规 | 7 | **9** | Abstract 单段 193 词 ✅、Running Title 36 字符 ✅、字数声明已更正 ✅、Zenodo DOI 已填 ✅；14 板块仍齐备 |
| 参考文献规范 | 3 | **7** | 引用顺序制已建立、31 条全被引、3 处错配全修、新增 4 条可考文献；扣分：作者截断两套并存、2 处冗余 author-year 括注 |
| 图表与交叉引用 | 3 | **5** | Figure 1–6 全部有引 ✅、Supp 材料 18/22 已补引；扣分：表格首引顺序仍乱（Critical C2 未核销）、Table 6 零引用、S6/S7/S11/S12 悬孤 |
| 语言与术语 | 6 | **6** | k-mer 术语维持统一；TRUE×31、preliminary×7、borderline×7 未改善 |
| 结构完整性 | 9 | **9** | 维持 |

---

## 剩余修复清单（按优先级，全部机械性）

1. **表格重排序（最高优先）**：当前首引序 7,10→1→4→3→2→5→8→9→11。二选一：
   - (a) 按首引行号重映射编号：旧7→1、旧10→2、旧1→3、旧4→4、旧3→5、旧2→6、旧5→7、Table 6 补引后→8、旧8→9、旧9→10、旧11→11；或
   - (b) 删除/改写 L168 Methods 中对 Table 7/10 的前瞻引用（改为 "see Results"），使首引从 Results 开始自然单调。推荐 (b)+(a) 组合，工作量最小。
2. **Table 6 补引**：在 L327 句末或 L342 段首加 "(Table 6)"。
3. **补 4 处 Supp 引用**：Supp Fig S6→Discussion 效应量处（约 L498）；Supp Table S7→BAliBASE 句（L502 附近）；Supp Table S11/S12→SwissTree 一节（L416–427）。
4. **统一作者截断**：refs 6/8/11/15/16/17/20/22/23/25/27/29 补全作者列表（≤10 后接 et al.）。
5. **删 2 处冗余 author-year 括注**（L369、L416）。
6.（可选）"TRUE tree" 首次定义后改小写；preliminary/borderline 各减至 ≤2 处。

---

## 推荐意见

**Minor Revision**（由 Major 降级）

Critical 4 项中 3 项已核销，剩余 C2（表格顺序）与新退化（Table 6 零引用）同属一个问题域，一次重编号即可同时解决；Major 剩余项均为 10 分钟级机械修复。修复清单 1–5 完成后可直接投稿。

## 投稿就绪度评分

**74 / 100**（v2.4: 55 → v2.5: 74，+19）

- 章节骨架与硬性指标：37/40
- 参考文献体系：18/25（顺序制 ✅、错配 ✅；作者截断与冗余括注 −7）
- 图表交叉引用：10/20（表格乱序、Table 6 零引用、4 项 Supp 悬孤 −10）
- 语言与文体：9/15（TRUE 大写与防御性措辞未改善 −6）

完成剩余修复清单后预计 90–93 分，即投稿就绪。
