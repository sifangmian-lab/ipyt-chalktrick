# Plan: 《粉笔虚线现象（Chalk Trick / Stick-Slip）文献综述》

## 背景
源文件：/mnt/agents/upload/2903.pdf —— "The Chalk Trick Pre-selection Report"（IYPT 风格报告）。
主题：粉笔-黑板虚线现象的 stick-slip 机理、摩擦学、弹性梁振动、混沌行为。
目标：生成一份中文学术文献综述（.md + .docx），为该报告的进一步研究铺垫。

## Stage 1 — 文献调研（并行两个 explore 子代理）
- A: stick-slip 基础理论、摩擦机理（Rabinowicz、Bowden-Tabor 接触结理论、Coulomb 模型局限、速度/时间依赖摩擦）
- B: 摩擦诱导振动与应用场景（琴弦、刹车尖叫、粉笔黑板）、弹性梁自由振动建模（Euler-Bernoulli）、stick-slip 混沌（Wiercigroch）、粉笔问题已有研究（IYPT、Lewin 视频、相关论文）
- 子代理可用 web_search / scholar 插件工具（如有）。
- 输出：各自写调研简报 md 到 /mnt/agents/output/research_a.md, research_b.md

## Stage 2 — 写作（paper-writing 技能）
- 加载 /app/.agents/skills/paper-writing/SKILL.md
- 一个写作子代理，基于 2903.pdf 内容 + 两份调研简报，撰写文献综述 md（引言、摩擦学基础、stick-slip 理论与建模、摩擦诱导振动的典型系统、粉笔问题研究现状、研究空白与展望），含规范引用。
- 输出：/mnt/agents/output/literature_review.md

## Stage 3 — 排版（docx 技能）
- 加载 docx 技能，将 final md 转为 /mnt/agents/output/粉笔虚线现象文献综述.docx
