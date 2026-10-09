# NAR 稿件第四轮整体审稿意见（v4 · v2.8 稿件）

**稿件**: Fusang: Tardigrade Edition — K-mer Frequency Vector Alignment-Free Phylogenetic Inference Resilient to Indel-Rich Sequence Evolution
**评审人角色**: 资深学术编辑（整体学术质量与影响力）
**评审日期**: 2026-10-09
**评审基准**: 完整重读 v2.8 全文（772 行），对照 v2.4 轮（review_v3_general.md）与 v2.5 复核（review_v3_verify_general.md）的历史意见

---

## 总体印象

v2.8 相对 v2.4 发生了实质性而非装饰性的升级：三项此前被判为"固有、无法靠编辑修复"的保留意见——L3 小样本、kmacs/Skmer 缺测、真实 indel-rich 数据缺失——全部以新实验而非新措辞解决，且每一项解决方式都附带了额外的诚实披露（旧 L3 实现曾漏掉 `-nt` 标志、E. coli 上 0.50 vs 0.08 的不利结果如实报告、storage-corruption 事件公开）。真实数据章节是最大增量：fish mtDNA 上 multi-k cosine 达到社区基准最优（0.045）且管线经已发表分数交叉验证，SwissTree gap 率中位数 50% 的量化巧妙地把既有蛋白基准重新定性为真实 indel-rich 验证。稿件当前的科学姿态、统计纪律与证据结构已达到 NAR Methods 的发表标准；剩余问题集中在 Abstract 措辞的显著性表述、fish 数据集区分度的语境呈现、以及 L3 协议变更的说明——全部是文字级修复，无需新实验。

## 逐项打分（1-10）

| 维度 | 分数 | 理由 |
|------|:---:|------|
| 写作质量 | 8.5 | L3 叙事一体化、竞品章节条理清晰、16S/EC 两个边界案例处理成熟；Abstract 数字密度偏高（197 词塞入 8 组统计量），全文词数注释未更新 |
| 新颖性 | 6.5 | 仍为扎实增量，但"k 随序列长度缩放、multi-k ensemble 自动解决"是 v2.8 新增的实质性实用发现，比 v2.4 的纯评估定位更有贡献感 |
| 诚实性与平衡性 | 10 | E. coli 适用边界如实报告、旧 L3 协议错误主动披露、storage-corruption 事件公开、k=5 基因组尺度饱和如实展示——唯一扣分项是 Abstract 对 p=0.024 用"significantly"而未带校正后 borderline 语境（正文与贡献声明均诚实标注） |
| 可复现性 | 8.5 | Zenodo DOI、175 文件 repro package、kmacs/Mash 复现已发表分数交叉验证管线、Skmer 不适用记录归档；storage-corruption 的"byte-identical to originals"验证逻辑表述需精确化 |
| 影响力 | 7 | 真实数据验证 + 明确的适用域地图（深分歧好/浅分歧差/结构化基因差）+ 工具交付，使 NAR 读者可以据此做出实际方法选择；基因长度真实 DNA 数据仍是声明过的 gap |

## 相对上轮评审的变化（实质提升 vs 增量修补）

**实质提升（三项均为新实验，非措辞修补）**：

1. **L3 Linux 全量验证（Table 5）**：从 5-seed Windows 无效比较（且实现有误——漏 `-nt` 标志，即蛋白质模式距离用于 DNA，v2.4 的 L3 数字整体无效）升级为 30-seed Linux 修正协议完整比较，结果方向反转且对作者有利（L1 0.583 vs L3 0.601, p=0.024, 胜 19/30）。更重要的是处理方式：作者没有掩盖旧实现的错误，而是公开 supersession 并在 Limitations 中保留"deliberately harsh protocol, nRF≈0.58–0.60"的绝对精度语境。这是我历轮意见中见到的最佳回应——不仅补了实验，还修正了一个此前无人发现的实现级错误。
2. **竞品补测（Table 7 重建 + Table 12）**：kmacs 以最接近的古典同源方法身份被系统测试（k 扫描 k∈{3,5,10}，best k=3 仍差 1.6×，p=5.6×10⁻⁶）；Skmer 的结构性不适用（500bp 上 division-by-zero，k=31/k=21 均测）被记录而非略过；Table 7 在统一 seed set 100–126 上重建，消除了嵌合体数据问题。竞品覆盖从"两个不完整对比"升级为"四方法系统 head-to-head + 两方法不适用记录"。
3. **真实数据验证（Table 12，全新章节）**：fish mtDNA (n=25, 9.1% gaps) multi-k 达 0.045 与社区最优打平；管线用 kmacs/Mash 已发表配置复现已发表分数（差 ≤1 split）完成外部交叉验证——这是方法论上很聪明的一步，直接预防了"你们自己算的分数可信吗"的攻击；E. coli/Shigella 上 0.50 vs 锚点 0.08 如实报告为适用边界；SwissTree gap 率中位数 50%（11–81%）把既有蛋白基准升格为定量证明的 indel-rich 真实数据。

**上轮遗留问题的修复确认**：Table 3 列头已改 TRUE-rel（L259）✓；Table 4 已标注 FT2-relative + seed set（L284）✓；Methods spaced k-mer 记号已统一为 1001001001001/1010101 且与 Table 10 一致（L62）✓；贡献声明 1 已改写为以 30-seed L1 vs L0（p<0.0001, d=3.82）为主体（L48）✓；13.3%→13.2% 舍入统一（L255/L266/L758）✓；统计框架明确了 confirmatory/exploratory 家庭划分（L145）✓；边界分类器补了 100% vs CV 0.84 落差的解释（L349）✓。上轮 5 条"一日可修"清单中 4.5 条已完成；唯一未完成的是 Abstract 分类器句仍缺"simulated"限定（L26）。

**增量修补**：n=1000 阈值未在过渡区直接验证的披露（L543）、特征消融 future work（L557）、16S 恢复率比较降为 descriptive（L408）——均为诚实性微调，不改变证据结构。

## 当前最弱点（审稿人最可能攻击处，按严重度排序）

1. **fish 数据集的区分度问题（最可能被攻击）**：n=25（44 informative splits，0.045=2 splits）上，multi-k、kmacs、Mash、Co-phylog、MAFFT+FastTree2 **全部**打平在 0.045——这是一个天花板数据集，连 Mash 都能达到最优。审稿人会说："在一个任何合理方法都能拿 0.045 的数据集上打平 MSA+ML 不构成区分性证据；'matches MSA+ML'只能证明你没拖后腿。"作者已有部分防御（2-split 语境、交叉验证管线），但 Table 12 finding 1 的标题措辞（"matches both MSA+ML and the best published alignment-free methods"）把无区分度的打平当作正面证据陈述，且该句进入了 Abstract。正确框架应是"达到社区基准天花板（所有最优方法在此收敛）"而非"匹配 MSA+ML"。这是措辞问题而非数据问题，但正是审稿人最爱抓的那类措辞。
2. **Abstract 显著性表述过强**：Abstract 写"significantly outperforms MAFFT+FastTree2 (nRF=0.583±0.044 vs 0.601±0.055, 30 seeds, p=0.024)"，但正文与贡献声明均诚实标注 Bonferroni 后 p_adj=0.071 borderline、该比较属 exploratory。正文诚实 + Abstract 强硬的组合，会被统计敏感的审稿人定性为选择性表述——与 v2.4 轮"L3 n=5 进 Abstract"同构的问题模式：作者把诚实的限定留在正文，把好看的数字放进 Abstract。修复成本一行。
3. **L3 协议变更未解释**：Table 5 从 v2.4 的（隐含 L=500）变为明确 L=1000 bp，与主基准协议（Tables 1/3, L=500）不同。表注声明了"不可跨基准比较"，但没有任何一句解释**为什么**改用 L=1000。审稿人会怀疑这是为让 L1/L3 比较更有利而调参（长序列下 MSA 成本上升、k-mer 频率估计更稳——两个方向都可能有利）。一句动机说明（或 L=500 下的补充结果）即可消除。
4. **storage-corruption 披露的验证逻辑（次级）**：Data Availability 写"regenerated deterministically... validated byte-identical to originals"——若原文件已损坏，"byte-identical to originals"是与什么比对的？合理推测是校验和或备份，但表述上自相矛盾，会引来追问：哪些表受影响？验证链是什么？建议在 repro package 文档中给出精确的受影响文件清单与校验和比对记录，正文引用之。

## 必须修改（Major）

1. **Abstract 的 p=0.024 表述降级**（L26）：改为"significantly outperforms MAFFT+FastTree2 (p=0.024; borderline after multiple-comparison correction)"或"matches or exceeds MAFFT+FastTree2"，与正文/贡献声明的诚实标注对齐。一行修复。
2. **fish 结果的证据定位修正**（L494 finding 1 标题句 + Abstract L26）：将"matches both MSA+ML and the best published alignment-free methods"重新定性为"reaches the community-benchmark ceiling (all top methods converge at nRF=0.045 on this small deep-divergence dataset)"，并在 finding 1 中明确指出该数据集**不区分**最优方法——证据力在于"免比对、免选 k 达到社区最优"，而非"匹配 MSA"。这反而更准确且同样有力。
3. **L3 协议变更的动机说明**（L309 表注或 Methods）：一句话解释为何 pipeline-level 基准采用 L=1000 bp（如"to ensure sufficient k-mer frequency signal at multi-k resolutions / to align with the coalescent protocol's signal density"），或声明该协议在 v2.7 重建时即为预先设定、与旧 5-seed 运行协议一致/不一致的明确交代。
4. **storage-corruption 披露精确化**（L657）：列明受影响的基准/表格清单，明确"byte-identical"的比对对象（原始校验和/备份），并指向 repro package 中的具体文档。

## 建议修改（Minor）

1. Abstract 分类器句加"simulated"限定（L26："88/88 **simulated** scenarios"）——上轮遗留，仍未完成。
2. 全文词数注释更新（L772 仍写 ~9,800 词；v2.8 实际主文本约 11,000+ 词）。
3. Abstract 数字密度：197 词含 8 组统计量，建议删去 Table 4 的 ensemble p=0.006 细节（Abstract 已有 L3 结果承载 ensemble 主张），保留 4–5 组核心数字。
4. Table 12 中"—"的语义注明（未测试 vs 不适用），例如 kmacs (k=3) 在 E. coli 列的 "—" 是超时（脚注³已解释，但表内符号宜统一为 "n/a (timeout)"）。
5. Limitations 中"remaining gap is real DNA data at gene length (500–1,000 bp)"（L545）建议同时前瞻性地指出 1–2 个候选数据类型（如病毒准种、快速进化基因家族），把 gap 转化为可检验的下一步——目前的表述已诚实但偏被动。

## 总体推荐

**Minor Revision**

三项历史遗留的实质性保留意见已全部以新实验解决，解决过程还额外暴露并修正了一个此前未知的实现级错误（L3 的 `-nt` 遗漏），这种"修一处、亮一处"的处理方式显著提升了稿件可信度。当前的全部 Major 问题均为措辞与文档级修复（Abstract 显著性、fish 证据定位、L3 协议动机、storage 披露精确化），无一项需要新计算或新实验。完成上述 4 项 Major 与 5 项 Minor 后，本稿达到投稿状态；若审稿人按此框架运作，预期审稿结果为 Minor Revision 至 Accept 之间。

**投稿就绪度：88/100**（v2.4 轮 62 → v2.5 轮 80 → 本轮 88）

- 剩余扣分：Abstract 显著性表述（−4）、fish 证据定位（−3）、L3 协议动机缺失（−2）、storage 披露表述（−2）、"simulated"限定残留（−1）
- 若 4 项 Major 全部完成，预计达 93–95，剩余空间属于语言润色与审稿人不可控偏好。
