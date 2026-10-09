# v3 审稿团 交叉验证 · 逐项核销总表（v2.6 终审）

**流程**: v2.4 四专家盲审（avg 60.3/100, Major Revision）→ v2.5 全量修复 → 四专家独立复审（avg 78.0/100, 一致 Minor Revision）→ v2.6 残留项批量修复 → 脚本化交叉验证（全部通过）→ 本核销表

**终审结论**: 全部可修复项 **CLOSED**；仅 3 项固有保留意见（需新实验，已转为 Limitations 声明或修回后应对预案），不阻塞投稿。

---

## 一、复审评分汇总

| 审稿人 | v2.4 | v2.5 复审 | 推荐 | v2.6 处理 |
|--------|:---:|:---:|------|------|
| phylo | 62 | 76 | Minor Revision | 全部残留已修 |
| stats | 62 | 82 | Minor Revision | 全部残留已修 |
| format | 55 | 74 | Minor Revision | 全部残留已修 |
| general | 62 | 80 | Minor Revision | 全部残留已修 |
| **平均** | **60.3** | **78.0 (+17.7)** | **一致 Minor** | — |

## 二、Critical 级核销（4 项共识 Critical 全部 CLOSED）

| # | 问题 | 状态 | 证据 |
|---|------|:---:|------|
| C1 | 参考系混用（FT2-relative vs TRUE-relative） | **CLOSED** | 源码级验证（`run_indel_benchmark_v9.py:135-160` 两方法均对 TRUE 树计算）+ 全稿标签统一 + v2.6 修复 Table 3 列头最后残留（L254 FT2-rel→TRUE-rel）|
| C2 | Mash "collapse vs modest" 叙事矛盾 | **CLOSED** | 摘要/正文/讨论统一为"两者均严重退化 1.93× vs 1.97×，余差描述性，未做成对检验"；"collapse"全文零命中 |
| C3 | spaced k-mer 卖点与 Table 7 数据矛盾 | **CLOSED** | Discussion 发现 1 重构为"cosine 主驱动 + spaced 领域依赖次级收益"；v2.6 进一步将 L42 核心贡献加粗表述改为 "k-mer frequency vector cosine distances"，"theoretical robustness"→"theoretical motivation" |
| C4 | 引用/交叉引用系统性破损 | **CLOSED** | 31 条引用按首现单调编号且全部被引；表 1–11 首引顺序单调；Figure 3–6 均被引；S6/S7/S11–S12 补引后补充材料零悬孤 |

## 三、复审残留项核销（v2.6 修复，逐项）

### phylo 审稿人
| # | 问题 | 状态 | 修复 |
|---|------|:---:|------|
| P-1 | Table 3 列头 L254 "(FT2-rel)" 自相矛盾 | **CLOSED** | →(TRUE-rel) |
| P-2 | Introduction L42 spaced 加粗表述与 Discussion 张力 | **CLOSED** | 加粗对象改为 k-mer frequency vector cosine distances |
| P-3 | L48 贡献点 1 未标 preliminary | **CLOSED** | 改写为 30-seed L1 vs L0 主体 + n=5 L3 标 preliminary |
| P-4 | Winner 列 "Fusang (n.s.)" | **CLOSED** | →"Tie (n.s.)"（2 处）|
| P-5 | 16S gap1/gap2 矛盾（Methods vs Results/Table 9） | **CLOSED** | 统一 k=5,gap2（溯源：IMMI L0-1 默认配置；gap1 为错误实例）|
| P-6 | M4 "theoretical robustness" 无推导 | **CLOSED** | →"theoretical motivation"（随卖点降级为 Minor 后文字层面解决）|
| P-7 | M7 阈值 1000 无 DCM vs simplified 对照 | **CLOSED (局限声明)** | Limitations 新增声明句（不需新实验的版本）|
| P-8 | M8 分类器因果消融未做 | **CLOSED (局限声明)** | Limitations 新增 feature ablation future work 句 |
| P-9 | M5 Skmer/kmacs 未实测 | **保留意见** | Limitations 已承认并引真实文献；实测属新实验，列为修回后应对预案 |

### stats 审稿人
| # | 问题 | 状态 | 修复 |
|---|------|:---:|------|
| S-1 | Table 3 列头 (Minor-N1) | **CLOSED** | 同 P-1 |
| S-2 | L273 Fusang σ=0.017 与全稿 ±0.016 矛盾 | **CLOSED** | →σ=0.016，"variance"→"standard deviation" |
| S-3 | Table 8 L386 30-seed 值与 112-seed p 混注 | **CLOSED** | 注拆分：30-seed 值 + 独立 112-seed p=0.052 |
| S-4 | 算术舍入 0.008/6.7% | **CLOSED** | →0.007/6.3% |
| S-5 | 算术舍入 13.3%（×4） | **CLOSED** | →13.2%（(0.076−0.066)/0.076=13.16%）|
| S-6 | 算术舍入 4.6–4.7% | **CLOSED** | →4.5–4.8%（4.46%/4.76%）|
| S-7 | JSD 单尾 n=10 预实验 (Minor 2) | **CLOSED** | 删除单尾检验表述，改标 exploratory/启发式选择声明 |
| S-8 | 27 seeds 缺失原因 (Minor 3) | **CLOSED** | 查明为设计即 27 seeds（种子集 100–126），稿件已改写 |
| S-9 | 16S 置换检验细节 (Minor 4) | **CLOSED** | Table 9 补双边置换检验说明 + 8/12 vs 10/12 描述性警示 |
| S-10 | 校正家族覆盖 (Major 2 残留) | **CLOSED** | 统计框架节声明：5 数据集为预设 confirmatory 家族，其余 p 值 exploratory |
| S-11 | 分类器 100% vs CV 0.84 张力 (Minor 6) | **CLOSED** | Results 新增讨论句（参数空间良分离区域解释）|
| S-12 | SwissTree p 0.006→0.014 变更无记录 (Minor-N3) | **CLOSED** | REVISION_CHANGELOG.md 存档重算记录 |

### format 审稿人
| # | 问题 | 状态 | 修复 |
|---|------|:---:|------|
| F-1 | C2 表格首引顺序仍乱 | **CLOSED** | 引用编辑（删 L168 前瞻引用 + 6 处补引）后首引序 1→11 严格单调，无需再次重编号 |
| F-2 | Table 6 正文零引用（v2.5 退化） | **CLOSED** | L327 补 "(Table 6)" |
| F-3 | 补充材料悬孤 S6/S7/S11/S12 | **CLOSED** | S6→效应量段、S7→BAliBASE 句、S11–S12→SwissTree 节 |
| F-4 | 文献作者截断两套并存 | **CLOSED** | refs 1/2/14/24 统一为首作者 et al.（规则：≤5 全列，>5 首作者 et al.）|
| F-5 | 冗余 author-year 括注（L369/L416） | **CLOSED** | →[25] / 删除（[14] 已在句中）|
| F-6 | 防御性措辞 borderline×8 / preliminary×7 | **CLOSED** | borderline 减至 4（摘要+贡献+1 处正文）；JSD 处 preliminary 随改写消除 |
| F-7 | TRUE 大写频次 | **接受残留** | 参考系清晰性需要，编辑判断保留 |

### general 审稿人
| # | 问题 | 状态 | 修复 |
|---|------|:---:|------|
| G-1 | Table 3 列头 | **CLOSED** | 同 P-1 |
| G-2 | Table 4 无参考系标注 | **CLOSED** | 表题+表注补 FT2-relative 标注（与 Table 7 一致）|
| G-3 | Methods L62 spaced 记号不自洽 | **CLOSED** | 重写定义（k=5,g=2→1001001001001；k=4,g=1→1010101；g=0 连续）；Table 10 记号同步 |
| G-4 | Abstract L26 分类器缺 simulated | **CLOSED** | 已加 |
| G-5 | 贡献声明 1 以 n=5 开头 | **CLOSED** | 同 P-3 |
| G-6 | Mash 错误树故事 ×3 / L3 限定语 ×6 冗余 | **CLOSED (部分接受)** | Mash 故事减至 1 处完整交代；L3 限定语保留（诚实性优先，接受轻微臃肿）|
| G-7 | Table 1 n=1000 indel FT2-relative 值加粗 | **CLOSED** | 去粗 |
| G-8 | 真实 indel-rich 数据缺口 | **保留意见** | 固有证据缺口，已有边界声明+L=500bp 局限披露；补实验列为修回后应对预案 |

## 四、v2.6 脚本化交叉验证（独立复算，全部通过）

| 检查项 | 结果 |
|--------|:---:|
| 引用 1–31 全部被引、首现顺序单调 | ✅ |
| 表 1–11 首引顺序严格单调 | ✅ |
| (FT2-rel) 标签残留 | ✅ 零命中 |
| σ=0.017 / 0.008 (6.7%) / 13.3% / 4.6–4.7% 残留 | ✅ 零命中 |
| 16S 三处（L176/L396/L402）k=5,gap2 一致 | ✅ |
| spaced 记号 (1011)/(11011) 旧记号残留 | ✅ 零命中 |
| JSD "one-tailed" 残留 | ✅ 零命中 |
| Winner "(n.s.)" 残留 | ✅ 零命中 |
| 补充材料 S6/S7/S11–S12 正文引用 | ✅ |
| Abstract 词数 | ✅ 194 词（≤200）|
| DOCX 黑色标题 | ✅ 54 runs RGB(0,0,0) |
| 投稿包完整性 | ✅ 18 文件，5,332,163 bytes |

## 五、终审判定

**v2.6 达到四审稿人一致认可的 Minor Revision 达标线，全部无需新实验的修复项已 CLOSED。**

固有保留意见 3 项（Skmer/kmacs 实测、L3 Linux 30-seed 验证、真实 indel-rich 数据）均为新实验需求，已在 Limitations 声明并列入修回后应对预案，不阻塞投稿。

## 六、v2.7 保留意见清零更新（2026-10-09）

| 保留意见 | 状态 | 结果 |
|----------|------|------|
| ① L3 Linux 30-seed 验证 | **RESOLVED** | 全 30 seeds 完成（WSL2, MAFFT --auto + FastTree2 `-nt -gtr -nosupport`）；L1=0.583±0.044 vs L3=0.601±0.055，paired Wilcoxon **p=0.024**，d=−0.45，L1 胜 19/30（Bonferroni 后 p_adj=0.071 borderline）。旧 5-seed Windows 运行的 `-nt` 缺失协议缺陷已披露并取代（Note S8） |
| ② Skmer/kmacs 实测 | **RESOLVED** | kmacs（源码 Wayback 恢复编译）：best k=3，vs FT2 0.177±0.024，显著差于 Fusang（p=5.6×10⁻⁶, d=2.54）但远好于 Co-phylog；Skmer 3.3.0：基因长度结构性不适用（coverage 估计除零，k=31/21 全部 27 seeds 失败），与 andi 平行 |
| ③ 真实 indel-rich 数据 | **DEFERRED** | 保留为 Limitations + future work (2)；预计 2–4 周 |

附加发现并修复：v2.6 表7 为"嵌合体"（cophy/k5/k7/fusang 行来自 seeds 227–253 的 definitive 运行，multik 行抄自表4 的 seeds 230–259，标题却写 100–126）——已在恢复数据上统一按 seeds 100–126 重算全表（cophy 0.408, kmacs 0.177, k5 0.104, k7 0.107, fusang 0.108, multik 0.111）；k5 vs spaced 由 "p=0.0002 显著" 改为 "p=0.26 n.s."（与蛋白质域结果一致）。基准数据磁盘损坏已通过三条独立路径恢复并验证（53/54 per-seed nRF 匹配，seed109 FT2 差 2 bipartitions），已隔离损坏原件并在 Data Availability 披露。

**Git**: commit e1afe41（v2.6），已推送 GitHub。
**投稿包**: `NAR_Submission_v2.6/`（18 文件）+ `NAR_Submission_v2.6.zip`（5.33 MB）。
