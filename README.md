# ipyt-chalktrick

粉笔虚线（Chalk Trick）现象——用粉笔在黑板上划线时出现虚线/断线现象的实验研究与理论建模项目。

## 项目简介

本项目研究粉笔划线"虚线戏法"现象：粉笔在黑板上运动时，划痕呈现周期性虚线而非连续实线。项目包含文献综述、实验方案与指导手册、理论模型（N 下降段完整模型、划痕相对强度时空依赖模拟）以及数值计算代码。

## 阅读顺序与当前有效版本

接手本项目的 agent 请先阅读 [Agent 交接说明](docs/研究记录/Agent交接说明.md)，再阅读 [研究进度评审](docs/研究记录/研究进度评审_2026-10-24.md)、[参数—证据总表](docs/研究记录/参数与证据总表.md)、[模型更新影响评审](docs/研究记录/模型更新影响评审_2026-10-25.md) 和 [模型 v2 规格](docs/研究记录/模型v2规格_2026-10-25.md)。

- **实验方案当前有效版本：**[虚线戏法 ChalkTrick 实验与理论建模方案（Markdown）](docs/实验指导/虚线戏法ChalkTrick实验与理论建模方案.md)（含 2026-10-25 第十节增补）；同目录 DOCX 为历史排版版。**2026-10-25 增补 DOCX/PDF：**[方案增补 DOCX](docs/实验指导/虚线戏法实验方案_2026-10-25增补.docx) / [PDF](docs/实验指导/虚线戏法实验方案_2026-10-25增补.pdf)；[模型报告更新与勘误 DOCX](docs/模型报告/粉笔_N下降段_模型更新与勘误_2026-10-25.docx) / [PDF](docs/模型报告/粉笔_N下降段_模型更新与勘误_2026-10-25.pdf)。**A2 备注：**两个增补 DOCX 已按预载修正后的数字重生成；对应 PDF 仍为 A2 前版本（本机无 Word/COM/LibreOffice），待补转。
- **模型更新状态（2026-10-25，含 A2 修正）：**原主模型与 S1–S3 是 **v1 退化模型**（固定铰轴、轴向冻结、尖端法向罚弹簧、无真正黏着），已归档为数值基准。当前目标是 **v2：弹性杆 + 手端正交弹簧（$k_x,k_y,k_\theta$）**。A2 修正了 v2 手端 y 弹簧的预载符号，并新增刚性接触（`contact='rigid'`，LCP/事件驱动）：S4 基线 A 得到“**$N$ 平台 + 有限强度截断**”（$N$ 涨落 0.25%、$I_c/I_{\rm plateau}=1.002$）；S5b 的刚性极限给出 $t_{\rm stop}=1.550$ ms、$I_c/I_{\rm plateau}=1.000$、约束残差 $3\times10^{-13}$。罚接触与刚性极限的定量差约 7%–15%，$k_n$ 必须作为物理参数标定。
- **P1/P2/A2 增补：**S5 收敛性、S5b 接触极限（LCP）、S6 参数平面、S7 本构对比；E3 正交柔度模块、E5 痕端测量链、E6 轮廓提取、E6A 截断判据、E7/E8 更新；关键文献待补清单；回归测试 `code/tools/check_regression.py`（含 v2 基线 + 刚性接触）。
- **数值后处理修正：**主模型/ S2 的 $I(\xi)$ 直方图权重曾误用 $q\cdot|d\xi|$ 与 $\mu N\,dt$，已统一为 `code/deposit_profile.py`（权重 $q\,dt$，满足 $\int I\,d\xi=I_{2\mathrm{nd}}$）；S2/S3 强度图与结果说明已更新，主模型报告需配读 [勘误](docs/模型报告/粉笔_N下降段_报告勘误_2026-10-25.md)。
- **内部链接约定：**`code/deposit_profile.py` 是空间强度 $I(\xi)$ 的唯一后处理入口；`code/v2_model.py` 是 v2 仿真共享内核；`docs/研究记录/模型v2规格_2026-10-25.md` 是当前模型规格草案。
- **运行环境：**S1–S7 依赖 NumPy、SciPy、Matplotlib；仓库记录环境为 Anaconda Python 3.14（numpy 2.4.6 / scipy 1.18.0 / matplotlib 3.11.0）。2026-10-25 A2 复算在系统 Python 3.10.11（numpy 2.2.6 / scipy 1.15.3 / matplotlib 3.10.9）完成，回归全部通过；中文字体已加入 Windows 字体回退（msyh/simhei）。
- **数据状态：**当前仓库未见可用于最终标定的 E1/E4 原始数据与标定表，也没有截断的原始图像/力数据；E6A 截断判据目前只有执行清单 [E6A截断判据执行清单](docs/实验指导/E6A截断判据执行清单_2026-10-25.md)。不得把数值结果写成“实验已经证实”，所有 v2 数值结果均为占位参数下的数值预言。
- **符号约定：**接触刚度 **k_n**（旧 S1–S3 的 `k_c`）；结构/手模型刚度 `k_x`、`k_y`、`k_θ`（或 `EI`）、`k_a`；率态临界刚度 **k_crit**。旧文件名和旧结果表保留历史符号，按上述映射阅读。
- **旧版手册提示：**`虚线戏法第一阶段实验详解与搭建指导.pdf` 与 `实验E1指导手册_静摩擦系数及其时间依赖.pdf` 仍是早期单向弹性臂/手刚度口径；涉及手模型时以当前方案 Markdown 的 E3 规程和 v2 规格为准。

## 目录结构

```
├── docs/                    # 文档
│   ├── 文献综述/            # 文献调研与综述（Markdown + Word）
│   ├── 研究记录/            # 交接、评审、参数—证据总表、v2 规格、E7/E8 更新
│   ├── 实验指导/            # 方案 Markdown/DOCX、E3/E5/E6/E6A 规程与增补 PDF
│   └── 模型报告/            # N 下降段报告（PDF）+ 2026-10-25 勘误/增补
├── code/                    # 主模型、v2_model 内核（含刚性接触）、deposit_profile、S1–S7+S5b 源码
│   └── tools/               # 回归测试、文档生成、LaTeX 公式、字体/PDF 工具
├── simulations/             # S1–S7 + S5b：v1 基准、v2 基线/收敛性/接触极限/参数平面/本构对比
│                            # 各目录含 CSV、结果图与结果说明.md
├── figures/                 # 结果图、公式图与过程预览图
└── assets/fonts/            # 报告用字体（思源宋/黑体子集）
```

`docs/参考文献库/` 与部分研究记录目前在 `trunk` 确认，`simulations/` 与最新方案在 `develop`；两分支内容尚未同步，不要因某一路径在当前分支缺失就断定资料不存在。

## 分支结构

- `trunk`：主分支（由原始 `main` 重命名而来），保存整理后的稳定版本
- `develop`：开发分支，日常修改在此进行，稳定后合入 `trunk`

当前两分支有未同步内容。合并前须先列出差异并取得用户明确确认，不自行合并或覆盖。

## 运行与复现说明

- `code/粉笔_N下降段_完整模型.py` 仍含 `/mnt/agents/output/` 绝对路径，运行前请改为本地路径；`code/划痕相对强度_时空依赖.py` 已改为项目相对路径，输出到 `figures/结果图/`。
- S1–S7 与 S5b 依赖 NumPy、SciPy、Matplotlib，输出 CSV/PNG；S2–S7 支持 `--outdir`，S4–S7 另有 `--tmax`；S5b 用 `contact='rigid'` 调刚性接触。`code/deposit_profile.py` 可独立运行自检；`code/v2_model.py` 是 v2 共享内核（罚接触默认，`contact='rigid'` 为 LCP/事件驱动极限）。
- **回归测试：**`python code/tools/check_regression.py` 检查 v1 名义 $t^*$/$I_{2\mathrm{nd}}$、v2 基线（预载修正后 $N$ 涨落 <2%、$I_c/I_p>0.95$、$t_{\rm stop}\approx1.816$ ms）、v2 刚性接触（$t_{\rm stop}\approx1.550$ ms、约束残差 <1e-9、无释放）与 $I(\xi)$ 守恒；当前全部通过。
- **增补 DOCX/PDF 生成：**`python code/tools/build_update_docs.py`（需要 python-docx）；PDF 可由 Word/COM 转换。
- 修改数值模型后，至少检查分离事件方向、接触掩膜、初始内态、积分容差、网格/分箱敏感性和输出单位；空间强度必须满足 $\int I\,d\xi=I_{2\mathrm{nd}}$，不要用 `q*|dξ|` 或 `μN*dt` 代替统一后处理。
- `code/tools/render_latex_eqs.py` 依赖外部 LaTeX；`build_pdf.py` 依赖 ReportLab、Pillow 与 `assets/fonts/` 中的字体。

