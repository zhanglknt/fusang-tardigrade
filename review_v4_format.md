# NAR 格式合规评审报告（第四轮，v4）— v3-format-reviewer

**稿件**: NAR_MANUSCRIPT_REVISED.md (v2.8)
**对照基线**: review_v3_format.md（v2.4 首评）与 review_v3_verify_format.md（v2.5 复核）
**评审日期**: 2026-10-09
**方法**: 全文通读 + 全部交叉引用指标脚本化重测（引用首现序列、表格首引序、图表/补充材料正文提及计数、词数与字符数实测）。

---

## 总体印象

v2.8 是该稿件四轮修订中格式状态最好的版本：上轮遗留的全部 Critical 与 Major 交叉引用缺陷（表格首引顺序、Table 6 悬孤、S6/S7/S11/S12 悬孤、author-year 括注残留）均已核销，新增的 Table 12 与 Supplementary Table S16 无缝嵌入且编号体系保持 Table 1–12 首引严格单调、参考文献 1–31 首现严格单调、无未引用与超界条目。Abstract 实测 197 词、单段非结构式，Running Title 36 字符，均达标。ETE3 nRF 与内部 nRF 双约定的 caveat 在 Methods（L181）与 Table 12 表注（L488）两处声明，表述清晰。剩余问题均为 Minor 级（文末字数声明再次失实、Table 12 孤儿脚注 ¹、补充材料列表 Note/Table 交错排列），不构成投稿障碍。

---

## 逐项打分（1-10）

| 维度 | 分数 | 理由 |
|------|:---:|------|
| Abstract 合规 | 9 | 实测 197 词（≤250）、单段非结构式、含真实数据验证句；扣 1 分：Abstract 称 "significantly outperforms MAFFT+FastTree2 (p=0.024)" 而正文 L324/L551 明言 Bonferroni 校正后 borderline（p_adj=0.071），措辞强度不一致 |
| 结构完整性 | 10 | 14 个必备板块全齐且顺序正确；NAR Methods Article 的工具可用性要求由独立 Python 工具 + GitHub + Zenodo DOI（10.5281/zenodo.20746742）满足，无需 Web Server |
| 参考文献规范 | 8 | 31 条引用顺序制完美（首现 1→31 严格单调、无未引、无超界、无错配；Skmer=[29]/PatternHunter=[30] 对调正确）；扣 2 分：作者截断仍两套并存（refs 1/2/6/8/11/15/16/17/20/22/23/24/25/27/30 为首作者+et al.，其余全列） |
| 图表编号一致性 | 9 | Table 1–12 首引顺序实测严格单调（L210→230→230→282→307→332→353→387→406→426→451→469），Figure 1–6 各被正文引用 ≥1 次，Supp Fig S1–S7 与 Supp Table S1–S16 全部有正文引用（S11/S12 经 L426/L428 复数形式引用）；扣 1 分：Table 12 存在孤儿脚注 ¹（表注定义了 ¹ 但表头/caption 无 ¹ 锚点） |
| 补充材料组织 | 8 | S16 新条目完整（含 NCBI accessions、AFproject 配置、gap 量化），正文两处引用（L469/L181）；扣 2 分：列表排列将 Note S1–S6 插在 Table S10 与 S11 之间、Note S8–S10 置于 S16 之后，类型交错影响检索；建议按 Figure→Table→Note 三组分块重排 |

---

## 相对上轮评审的变化

**上轮（v2.5 复核）遗留问题的核销状态：**

| 上轮问题 | 状态 | 证据 |
|----------|------|------|
| C2 表格首引顺序乱（NOT-FIXED） | **CLOSED** | 首引序实测 1→12 严格单调；Table 6 已在 L332 补引 |
| Table 6 正文零引用（新退化） | **CLOSED** | L332 "(Table 6; per-scenario results: Supplementary Table S15)" |
| Supp Fig S6 / Table S7/S11/S12 悬孤 | **CLOSED** | S6→L276；S7→L545；S11/S12→L426、L428 |
| author-year 冗余括注 ×2（Haubold/Zielezinski） | **CLOSED** | 脚本扫描 body 无残留（Table 7 单元格内 "Yi & Jin 2013" 为作者署名式标注，非引用，可接受） |
| 作者截断两套并存（PARTIAL） | **仍开放（Minor）** | 与 v2.5 相同的两套风格 |
| TRUE×31 / defensive 措辞 | **部分改善** | preliminary 由 7 降至 0 ✅；borderline 7 次持平；TRUE×31 持平 |
| 文末字数声明 | **回退（新问题）** | 声明 9,800 词，实测正文（不含表格行）11,191 词——新增 AFproject 章节后未同步更新声明，低估约 14% |

**本轮新增内容的核查：**

| 项目 | 结果 |
|------|------|
| Table 12（AFproject 真实基因组基准） | ✅ 编号/位置正确（Scalability 之后，保持 1–12 单调）、正文 5 次引用、表注含 ETE3 不可比声明 |
| Supplementary Table S16 | ✅ 条目完整（L640）、正文引用 ×2（L181、L469）、含 NCBI accessions 声明 |
| Abstract 扩展 | ✅ 实测 197 词 ≤250（与自报 197 一致）、单段、新句数据与正文一致 |
| refs 29/30 对调 | ✅ 实测：Skmer=[29]（L353/L379/L561）、PatternHunter=[30]（L510），首现序列 ...29→30→31 单调，条目-编号匹配无误 |
| ETE3 vs 内部 nRF 双约定 | ✅ 表述清晰：Methods L181 与 Table 12 表注 L488 两处声明 "informative splits only vs counts trivial splits"，且 Table 12 标题即标注 "ETE3 nRF"；读者可明确区分两套数值不可直接比较 |
| Data Availability 扩展 | ✅ NCBI accessions 与 AFproject 数据集已声明；另含 "storage-corruption incident" 披露句（见 Minor 5） |
| Limitations / Practical recommendations 各 +2 bullet | ✅ 新增 boundary-classifier feature ablation、n=500–1000 过渡带未测、深分歧基因组建议、菌株级边界建议，编号引用均正确 |

---

## 必须修改（Major）

1. **文末字数声明失实（回退问题）**：L772 "Main text: approximately 9,800 words"，实测正文（Introduction–Discussion，不含表格行）11,191 词，含表格约 12,300 词。新增 AFproject 章节后未更新声明。改为 "approximately 11,200 words"（或压缩正文）。此为第二轮出现的同类问题，投稿前必须修正。

## 建议修改（Minor）

2. **Abstract 与正文的统计措辞强度不一致**：Abstract "significantly outperforms MAFFT+FastTree2 (nRF=0.583±0.044 vs 0.601±0.055, 30 seeds, p=0.024)" vs 正文 L324/L551 明言 "borderline after Bonferroni correction (p_adj=0.071)"、"interpreted as exploratory"。建议 Abstract 改为 "outperforms MAFFT+FastTree2 on a 30-seed indel-rich benchmark (p=0.024; borderline after multiple-comparison correction)"，与正文的审慎口径对齐。
3. **Table 12 孤儿脚注 ¹**：表注 L488 以 "¹ nRF computed with ETE3..." 开头，但表题与表头均无 ¹ 锚点。将 ¹ 挂到表题 "ETE3 nRF¹" 或表头 "nRF ↓¹" 上；² ³ 锚点正常（L485/L481）。
4. **补充材料列表类型交错**：现顺序为 Fig S1–S7 → Table S1–S10 → **Note S1–S6** → Table S11 → Note S7 → Table S12–S16 → Note S8–S10。建议按 Figure → Table → Note 三组分块重排（Table S11–S16 提到 Note S1 之前，Note S1–S10 连续排列），便于编辑与审稿人检索。
5. **Data Availability 中的 "storage-corruption incident" 句**（L657）：披露精神可取，但在 Data Availability 段落中详述数据恢复过程略显冗长且可能引发编辑不必要追问。建议压缩为一句（"Benchmark inputs affected by a storage incident were regenerated deterministically and validated byte-identical; recovery details in the reproducibility package"）并将完整过程移至 SI/Supplementary Note。
6. **参考文献作者截断统一**（遗留）：refs 1/2/6/8/11/15/16/17/20/22/23/24/25/27/30 补全作者（≤10 全列后接 et al.），与 refs 4/5/7/10/12/13/18/21/26/28/29/31 的全列风格一致。
7. **全大写降级**（遗留）：TRUE×31 建议首次定义（L181 附近已有 "ground-truth (TRUE) tree" 语境）后改小写 "true tree"；NOT×2（L334、L569）改为常规否定句式。
8. **borderline×7 去重**：同一 L1 vs L3 比较的 borderline 声明出现于 Abstract、L48、L324、L551、L568 共 5 处，保留 Methods 统计框架 1 处 + 正文结论 1 处即可。

---

## 总体推荐

**Minor Revision**

四轮修订后，全部 Critical 与 Major 交叉引用缺陷已清零，编号体系（参考文献 1–31、表格 1–12、图 1–6、补充材料 S1–S16）经脚本实测完全一致。剩余 1 项 Major（字数声明）为 5 分钟级修正，其余为 Minor 润饰。完成 Major 1 与 Minor 2–3 后即可达到投稿就绪状态；据格式维度判断，本稿已接近 Accept with minor changes。

**投稿就绪度：91 / 100**（v2.4: 55 → v2.5: 74 → v2.8: 91）
