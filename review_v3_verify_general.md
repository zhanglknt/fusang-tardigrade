# NAR 稿件 v2.5 复审核销报告（v3-general-reviewer）

**稿件**: NAR_MANUSCRIPT_REVISED.md（v2.5）
**复核基准**: review_v3_general.md 中提出的 2 Critical + 4 Major + 3 Minor
**复核日期**: 2026-10-09
**复核方式**: 独立通读全文，逐项对照原意见核查，不依赖修订说明

---

## 一、逐项核销

### Critical 1：参考系混用 → **PARTIAL**（主体修复，残留一处自相矛盾）

已修复部分：
- Table 1 标题与脚注已统一为 TRUE-relative，n=1000 indel 行以 † 明确标注 FT2-relative 且声明"not comparable with TRUE-relative values elsewhere"（L207、L218–219）
- Table 3 表题已改为"All nRF values are TRUE-relative"，表注明确"Both columns are TRUE-relative and directly comparable per seed"，13.3% 优势现于同一参照系内计算（L252、L261）
- n=1000 indel 的表述已从"high accuracy"改为"topological divergence/similarity"框架，并保留 FT2-relative 限定（L373、L500、L516）

**残留问题（必须修复）**：
- **Table 3 列头未改**：L254 仍为 `| Indel Rate | Fusang nRF ↓ (FT2-rel) | FastTree2 nRF ↓ (TRUE-rel) |`——与同表标题（L252"All nRF values are TRUE-relative"）和表注（L261"Both columns are TRUE-relative"）直接矛盾。这正是原 Critical 问题的可视化残留，细心审稿人一眼可见，会质疑"标签统一"是否仅为表面功夫。**一行即可修复**。
- Table 4（L279）完全没有参考系标注；其数值（0.112/0.105/0.099）与明确标注 FT2-relative 的 Table 7（L350）对应行完全一致，读者无法判断 Table 4 的参照系。需补一句标注。
- 次要：Table 1 中 n=1000 indel 的 0.037 仍以加粗呈现（L216），在 TRUE-relative 表格中加粗一个 FT2-relative 数值，视觉上仍暗示"最优"；13.3% 与按表内数字计算的 13.2%（(0.076−0.066)/0.076）存在舍入差异。

### Critical 2：Mash n=1 vs n=30 矛盾 → **CLOSED**

- L225 明确 Fusang 与 Mash"each with 30 seeds"；Table 2 注释（L236）、Results（L240–246）、Discussion（L510）、Practical recommendations（L526）与 Abstract（L26）全部统一为 30-seed 基准，n=1 叙述已彻底消除
- 更重要的是叙事已改为诚实版：原"MinHash collapse"主张降级为"both degrade severely (1.93× vs 1.97×), cosine retains marginally lower absolute error"，并明确声明"a formal between-method test on matched seeds was not performed, this small absolute difference should be interpreted descriptively"（L244）。历史错误值（0.162/1.005）以"development 期间错误解析、已被取代"的方式一次性交代，处理得当

### Major 3：16S 不利结果的叙事张力 → **CLOSED**（框架层面；证据缺口见攻击点预测）

- 已按建议重新框定为适用边界："This result delineates an applicability boundary of the method rather than contradicting the indel-robustness findings"（L398），L412、L502 口径一致
- 底层证据缺口（indel 优势仅有模拟证据）无法靠编辑修复，转入审稿人攻击点预测继续跟踪

### Major 4：绝对精度语境 → **CLOSED**

- 贡献声明 1（L48）保留"both methods remaining substantially distant from the TRUE tree (nRF≈0.58–0.59)"；Table 5 正文（L321）明确"58.3% of bipartitions differ"；读者不会在低水平相当性上被误导

### Major 5：L3 (n=5, p=0.24) 移出 Abstract → **CLOSED**（附残留提示）

- Abstract（L26，单段 193 词）已无任何 L3/n=5 数字，符合要求
- **残留提示（Minor）**：贡献声明第 1 条（L48）仍以 n=5 的 L1 vs L3 数字开头，n=5 初步比较在 Introduction 中仍占据最高曝光位置。虽有完整限定语，建议将贡献 1 改写为以 multi-k ensemble 的 30-seed 结果（L1 vs L0, p<0.0001, d=3.55）为主体、L3 仅作附带说明

### Minor 6：边界分类器 100% 准确率 → **PARTIAL**（残留轻微）

- Results（L344）已补同引擎限定；贡献声明 4（L54）已改为"100% accuracy on simulated scenarios"——明显改善
- 残留：Abstract（L26）"88/88 scenarios, 95% CI [0.958, 1.0]"仍未带"simulated"限定词，建议加一词

### Minor 7：写作冗余与过度对冲 → **PARTIAL**

- 整体叙事已收紧（尤其 Mash 部分），但：Mash 错误树文件故事仍出现 3 次（L236、L510、L526）；DCM 0.005/0.388 数字重复 3 处以上（L188–195、L488、L541）；L3 的 n=5 限定语在 L302、L315、L319、L323、L508、L525 重复 6 次——诚实但臃肿。建议在语言润色阶段对 L3 限定语做"首次完整+后续引用"处理

### Minor 8：可复现性承诺未闭环 → **CLOSED**

- Zenodo DOI 已分配（L608：10.5281/zenodo.20746742）；Docker 镜像（L180）、双平台二进制、原始基准数据、INDELible 配置均在声明中。本人未验证链接可访问性，建议投稿前实测一次

---

## 二、新发现问题

1. **[Minor] Table 3 列头自相矛盾**（L254 vs L252/L261）——见 Critical 1 残留，单列于此因其是当前全文最显眼的一致性硬伤。
2. **[Minor] Methods 中 spaced k-mer 记号混乱**（L62）："k=5, g=2 (gap1 notation: 10101 ... spanning 13 nucleotides)"——10101 仅跨 5 个位置，与"13 nt"不符；"For gap2 (11011011011), 3 positions are skipped"——该模式每位点间跳过 2 个位置而非 3 个；且与 Table 10 的记号（k=4,gap1=1011；k=5,gap2=11011）不自洽。gap 编号体系需要在 Methods 一处定义清楚。
3. **[Minor] Table 4 缺参考系标注**（L279）——见 Critical 1 残留。
4. **[Minor] 摘要词数**：实测 193 词（wc），与声称的 198 略有出入，无碍（NAR 上限内），仅作记录。
5. **[Minor] 13.3% vs 13.2%** 舍入差异（L250/L261/L711）。

---

## 三、修订后逐项打分（1-10）

| 维度 | v2.4 | v2.5 | 理由 |
|------|:---:|:---:|------|
| 写作质量 | 7 | 8 | Abstract 大幅收紧、叙事转向诚实版、16S 重框定流畅；扣分项为 L3 限定语六次重复与 Mash 故事三次重复 |
| 新颖性 | 6 | 6 | 定位未变（系统评估+ensemble+工程），诚实降级反而使贡献边界更清晰；仍为扎实增量 |
| 影响力 | 6 | 6 | 验证结构未变；SwissTree 跨域结果与工具交付仍是主要影响力支点 |
| 诚实性与平衡 | 9 | 10 | MinHash collapse 主张主动降级为描述性表述、未做配对检验主动声明、16S 重框定为适用边界、n=5 结果从 Abstract 撤出——作者对批评的回应方式本身就是方法学诚信的示范 |
| 可复现性 | 6 | 8 | Zenodo DOI 落地、Docker 镜像、双平台二进制、全量原始数据；留 2 分给投稿前的链接实测与一键复现脚本验证 |

---

## 四、审稿人攻击点预测（更新版）

1. **[可立即消除] Table 3 的自我矛盾**：列头"(FT2-rel)"与表题/表注"Both columns are TRUE-relative"直接冲突。审稿人会问："到底哪个标签是对的？如果这次是标签错误，上次是吗？"——这会重新打开已经基本关上的参考系问题。**投稿前必须修复（一行改动）**，连同 Table 4 补标注。
2. **[固有，已缓解但未消除] indel 优势的全部证据来自模拟**：16S 已被诚实地框定为边界证据，但稿件仍没有任何真实 indel-rich 数据集（病毒 quasispecies、快速进化基因家族）上的正面结果。审稿人仍可主张"simulation-only advantage"。缓解已到位（边界声明+L=500bp 局限披露），彻底消除需补实验，可作为修回后的应对预案。
3. **[新浮出，温和] 标题级贡献与数据方向的错位**：作者自己的数据显示 contiguous k=5 优于 spaced（Table 7, p=0.0002）、multi-k ensemble 全部由 contiguous k-mer 构成、SwissTree 上 spaced vs contiguous 无显著差异（p=0.31）——真正的精度驱动是 cosine 距离 + multi-k 融合，而非 spacing。作者已在 Discussion 发现 1（L469）中主动承认这一点，且主标题/Running title 已去掉"spaced"，缓解得当；但关键词、方法名与工具叙事仍以 spaced k-mer 为标识，审稿人可能要求进一步淡化或改写为"k-mer frequency vector + cosine distance"框架。

---

## 五、结论

**推荐意见：Minor Revision（原 Major Revision 升级）**

两个 Critical 中，Mash 矛盾已彻底关闭；参考系问题的实质重构已完成且方向正确，仅剩 Table 3 列头一行未同步——这是唯一可能重新引爆 Critical 问题的残留，但修复成本极低。Major 级问题全部关闭。v2.5 的修订质量高、回应诚实且有的放矢。

**投稿就绪度评分：80 / 100**（v2.4 为 62）

- 剩余扣分：Table 3 列头矛盾（−6）、Table 4 未标注（−2）、Methods gap 记号混乱（−2）、L3 在贡献声明中的曝光位置（−3）、真实 indel-rich 证据缺口（固有，−5）、写作冗余（−2）
- **达到 90+ 的路径**：Table 3 列头改为 TRUE-rel（一行）→ Table 4 补参考系标注（一行）→ Methods gap 记号统一（一段）→ Abstract 给分类器准确率加"simulated"限定（一词）→ 投稿前实测 GitHub/Zenodo 链接。以上全部可在一日内完成，无需任何新实验。
