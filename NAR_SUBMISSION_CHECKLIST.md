# NAR 投稿包清单 — Fusang: Tardigrade Edition

**目标期刊**: Nucleic Acids Research (NAR), Methods Article
**稿件类型**: Computational Biology / Phylogenetics
**截止日期**: 待定

---

## 一、NAR 投稿要求（基于官方指南）

### 稿件格式
- **文件格式**: 接受 Word (.docx) 或 LaTeX (Overleaf)
- **字体**: 推荐 Arial 或 Times New Roman, 12pt
- **行距**: 双倍行距
- **页码**: 底部居中
- **行号**: 连续行号（便于审稿）
- **图/表位置**: 稿件中需嵌入图/表位置，或统一放在文末

### 长度限制（Methods Article）
- **正文**: 约 5000 词（含摘要、参考文献）
- **摘要**: 约 150-200 词
- **参考文献**: 无严格数量限制，但 Methods 文章通常 30-50 篇
- **图表**: 无严格数量限制，但 NAR 彩色在线免费

### 图表要求
- **分辨率**: 最低 300 DPI（推荐 600 DPI）
- **格式**: TIFF 或 PNG（矢量图用 EPS 或 PDF）
- **色彩模式**: RGB（在线版），如打印版需 CMYK
- **字体**: Arial 或 Helvetica, 8-12pt
- **图注**: 单独提供（Figure Legends），不在图中嵌入
- **表注**: 表格上方需有说明

### 补充材料（Supplementary Data）
- **格式**: PDF 或 HTML（NAR Online）
- **图表**: 高分辨率 TIFF/PNG
- **长度**: 无严格限制，但需合理
- **引用**: 正文中用 "Supplementary Figure S1" 格式

### 必需文件清单
| 文件 | 状态 | 备注 |
|---|---|---|
| Manuscript (Word/LaTeX) | ✅ 草稿完成 | NAR_MANUSCRIPT.md → 需转 Word |
| Cover letter | ❌ 待写 | |
| Figures (高分辨率) | ⚠️ 需检查 | 现有 figures 可能不适合 NAR |
| Tables | ✅ 稿件中包含 | 需单独提取为 Word 表格 |
| Supplementary Material | ⚠️ 仅有列表 | 需撰写实际内容 |
| Author contributions | ❌ 待写 | NAR 要求贡献声明 |
| Conflict of interest | ❌ 待写 | 必须声明 |
| Data availability | ✅ 草稿完成 | 需补充 GitHub/Zenodo 链接 |
| Funding statement | ❌ 待填 | [To be added] |
| Reviewer suggestions | 可选 | 3-5 人 |
| Highlights (如有) | 可选 | NAR 不要求 |

---

## 二、当前稿件状态

### ✅ 已完成部分
1. **标题 + 摘要**: 完整
2. **引言**: 完整（含 Fusang v1 → v2 演进）
3. **材料与方法**: 完整（k-mer 提取、距离计算、pipeline、FastME 集成、参数选择、基准设计）
4. **结果**: 完整（简化 pipeline、spaced k-mer 优势、indel 鲁棒性、130-seed 基准、参数稳定性、16S rRNA 验证、可扩展性）
5. **讨论**: 完整（20 年鸿沟、架构演进、pipeline 简洁性、局限与未来工作、实践建议）
6. **数据可用性**: 草稿完成（需补充实际链接）

### ⚠️ 待完成部分
1. **参考文献**: 不完整（仅 17 条，实际引用约 23 条）→ 需补全
2. **经费声明**: 空白 → 需填写
3. **致谢**: 不完整 → 需补充
4. **补充材料**: 仅有列表，无实际内容 → 需撰写
5. **图注**: 稿件中引用了 Figure 1-6，但图注未单独提供 → 需撰写
6. **GitHub URL**: 占位符 `[GitHub URL]` → 需替换为实际链接

### ❌ 缺失材料
1. **Cover letter**: 未写
2. **作者贡献声明**: 未写
3. **利益冲突声明**: 未写
4. **图表文件**: 需检查现有 figures 是否符合 NAR 要求

---

## 三、行动计划

### 阶段 1：稿件完善（1-2 天）
- [ ] 补全参考文献（从正文引用中提取所有文献，格式化）
- [ ] 填写经费声明
- [ ] 完善致谢
- [ ] 替换 GitHub URL 占位符
- [ ] 撰写补充材料实际内容（Supplementary Figures S1-S6, Tables S1-S7, Notes S1-S3）
- [ ] 撰写图注（Figure Legends）

### 阶段 2：图表准备（1 天）
- [ ] 检查现有 figures（paper_figures/ 目录）是否适合 NAR
- [ ] 如需要，重新生成高分辨率 figures（600 DPI, TIFF）
- [ ] 确保图表字体一致（Arial, 8-12pt）

### 阶段 3：投稿材料准备（1 天）
- [ ] 撰写 Cover letter
- [ ] 撰写作者贡献声明
- [ ] 撰写利益冲突声明
- [ ] 整理 Data availability 声明（最终版）
- [ ] 准备 Reviewer suggestions（可选）

### 阶段 4：格式转换与检查（1 天）
- [ ] 将 NAR_MANUSCRIPT.md 转换为 Word (.docx) 或 LaTeX
- [ ] 检查格式（行号、页码、双倍行距）
- [ ] 检查参考文献格式（NAR 要求带括号的数字引用）
- [ ] 最终校对

### 阶段 5：投稿（1 天）
- [ ] 在 NAR 在线投稿系统注册
- [ ] 上传所有文件
- [ ] 填写投稿表格
- [ ] 提交

---

## 四、NAR 参考文献格式示例

NAR 使用 **带括号的数字引用**，如：
> ... multiple sequence alignment (MSA) [1,2].

参考文献列表格式：
```
1. Bernard,G. et al. (2019) Alignment-free inference of hierarchical orthologous groups. Nucleic Acids Res., 47, W202–W208.
```

注意：
- 作者姓名：姓氏全拼 + 名字首字母
- 年份括号
- 标题（句子大小写）
- 期刊名（缩写，斜体）
- 卷号（斜体），页码

---

## 五、优先级

1. **高优先级**（阻塞投稿）：
   - 补全参考文献
   - 填写经费声明
   - 撰写补充材料
   - 准备 Figures（高分辨率）

2. **中优先级**（投稿前完成）：
   - Cover letter
   - 作者贡献声明
   - 利益冲突声明

3. **低优先级**（可选）：
   - Reviewer suggestions
   - Highlights

---

**当前阻塞项**: 补充材料实际内容 + 参考文献补全

**预计完成时间**: 3-5 天（假设每天工作 4-6 小时）
