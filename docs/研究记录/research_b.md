# 文献调研简报：粉笔在黑板上画虚线（Chalk Trick）的 Stick-Slip 机理

> 调研日期：本轮任务；来源：Web/学术检索（Google Scholar 类、出版商页面、arXiv、大学仓库）。
> 原则：所有条目均为可核查的真实文献/资料；个别著录细节不把握的标注"待核实"。

---

## 1. 粉笔—黑板问题本身（Chalk Trick / Painlevé 悖论的直接研究）

### 1.1 E. Rabinowicz, "Stick and Slip", Scientific American, 194(5): 109–118 (1956)
链接（PDF）：http://www3.geosc.psu.edu/courses/Geosc508/Rabinowicz1956SciAm.pdf
- 经典科普性综述，系统讲解 stick-slip 机理（静摩擦大于动摩擦、摩擦系数随滑动速度下降的负斜率），并**直接用粉笔在黑板上留下的断续痕迹作为 stick-slip 的示范插图**：粉笔轴与运动方向成锐角时画连续线，成钝角时发生"粘—滑—粘"循环留下点划线。
- 指出握得越紧滑移距离越短——正是点间距的调节规律，可作为粉笔虚线问题最早的正式文献锚点。

### 1.2 IPT 2023 赛题 "The Chalk Trick"（International Physicists' Tournament）
- 赛题文本（pre-selection report 范例）：https://ipt.science/media/model%20reports/2903.pdf
- 题目陈述直接引用 Walter Lewin 画虚线视频：改变粉笔与黑板的接触角，连续运动可产生虚线；要求从痕迹反推运动参数与粉笔尺寸。注意这是 IPT（国际物理学家锦标赛），非 IYPT；IYPT 中未见同名正式题目（待核实）。示例解答用 slip-stick（弹簧-滑块+传送带）模型建模，认为角度是关键参数。

### 1.3 V. Sundström (Chalmers), "The Chalk Trick" 研究报告，收录于 I/IPT Conference Proceedings 2023–2024
PDF：https://zenodo.org/records/10819172/files/I_IPT_Conference_Proceedings_2023_2024-6.pdf
- 针对 IPT 2023 题目建立的 1D 弹簧-变摩擦力模型，可拟合实验点间距数据；并建造了机械装置复现画虚线。核心发现：点间距与画线速度呈强线性相关（Pearson 0.99）；弹簧刚度单调影响点间距；法向力影响不明显；粉笔触板角度限制可实现的画线速度。
- 是目前对粉笔虚线现象最系统的实验-模型对照研究，可作为综述中"直接证据"部分的核心。

### 1.4 A. R. Champneys & P. L. Várkonyi, "The Painlevé paradox in contact mechanics", IMA Journal of Applied Mathematics, 81(3): 538–588 (2016)
DOI: 10.1093/imamat/hxw027；arXiv:1601.03545
- 120 年 Painlevé 悖论权威综述，**明确以"推粉笔而非拉粉笔时在黑板上的奇异抖动（sprag-slip 振荡）"和 Walter Lewin 画虚线作为核心实例**，并给出 Lewin 手绘的粉笔痕迹图（经授权首次发表）。
- 提出该现象不只是 stick-slip，而是每个循环含离面（lift-off）与碰撞的 sprag-slip，可由 Painlevé 悖论经"reverse chatter"不稳定产生。为粉笔虚线提供刚性接触力学层面的严格理论框架。

### 1.5 R. I. Leine, B. Brogliato & H. Nijmeijer, "Periodic motion and bifurcations induced by the Painlevé paradox", European Journal of Mechanics A/Solids, 21(5): 869–896 (2002)
作者主页 PDF：https://www.inm.uni-stuttgart.de/institut/mitarbeiter/leine/papers/journal_publications/Leine_x_Brogliato_x_Nijmeijer_-_Periodic_Motion_and_Bifurcations_Induced_by_the_Painleve_Paradox.pdf
- 研究"摩擦碰撞振子"（Frictional Impact Oscillator）的周期运动与分岔，证明 hopping/bouncing（跳动）运动与 Painlevé 悖论（滑动模式解的非唯一/不存在）密切相关，且在物理上现实的摩擦系数下即可出现。
- 即 28 页"Painlevé 分岔"论文（Lewin 虚线现象讨论中常被引用）；给出粉笔跳动类系统的分岔结构，是连接"库仑摩擦刚性模型"与"粉笔弹跳虚线"的关键理论文献。

### 1.6 补充（背景资料，非学术文献）
- Walter Lewin 原始视频："Walter Lewin - Dotted lines"（YouTube, https://www.youtube.com/watch?v=xODkb-i-t2M）；其自述要领：粉笔不能太短、向前"推"而非"拉"、压而不死压即会跳动。
- Wikipedia "Painlevé paradox" 词条的 Physical realizations 一节即以 Lewin 画虚线为图例。Physics Forums 等讨论帖可作社会传播背景。

---

## 2. 摩擦诱导振动/噪声的典型系统

### 2.1 小提琴弓弦（Helmholtz 运动谱系）
- **C. V. Raman, "On the mechanical theory of the vibrations of bowed strings and of musical instruments of the violin family, with experimental verification of the results: Part I", Bull. Indian Assoc. Cultivation of Science, 15: 1–158 (1918)**——弓弦力学的奠基性专著，首次给出 Helmholtz 运动（尖角沿弦往返、每周期一粘一滑）的定量理论。与粉笔问题同源：都是"摩擦驱动弹性体的自激粘滑振荡"。
- **J. C. Schelleng, "The bowed string and the player", J. Acoust. Soc. Am., 53(1): 26–41 (1973)**——给出最小/最大弓压力与弓速、弓-琴码距离的参数图（"Schelleng 图"）；力不足出现多次滑移（表面噪声），力过大出现非周期粗音。类比于粉笔问题中"压力/速度参数平面内的不同运动模式分区"。
- **L. Cremer, "The Physics of the Violin", MIT Press (1984, 英译)**（另见 Fletcher & Rossing 的总结）——系统处理弓弦系统的现代理论；指出 Helmholtz 尖角在传播中变圆、过弓时被摩擦粘/滑切换重新"锐化"的机制。说明摩擦非线性在波传播系统中如何维持周期粘滑，对理解粉笔尖端的周期冲击有方法借鉴。
- **M. E. McIntyre & J. Woodhouse, "On the fundamentals of bowed-string dynamics", Acustica, 43: 93–108 (1979)**（卷期页码待核实）——把弓弦问题表述为反射波逐次映射，确立"摩擦特性曲线+弦的特征线"分析方法；现代弓弦模拟的基准文献。

### 2.2 刹车尖叫（brake squeal）
- **N. M. Kinkaid, O. M. O'Reilly & P. Papadopoulos, "Automotive disc brake squeal", Journal of Sound and Vibration, 267(1): 105–166 (2003)**——刹车尖叫权威综述，梳理 stick-slip、sprag-slip、负摩擦-速度斜率、模态耦合等机制；sprag-slip 与粉笔跳动在几何上完全同构（锐角接触+转动自由度耦合）。
- **A. Papinniemi, J. C. S. Lai, J. Zhao & L. Loader, "Brake squeal: a literature review", Applied Acoustics, 63(4): 391–400 (2002)**——早期机制综述，指出没有单一机制能解释全部尖叫事件。
- **A. Akay, "Acoustics of friction", J. Acoust. Soc. Am., 111(4): 1525–1548 (2002)**——摩擦声学综述，统一讨论各类摩擦发声（含粉笔尖叫、湿手指、琴弦），是从"声"角度组织文献的最佳入口。
- **R. T. Spurr, "A theory of brake squeal", Proc. Automobile Div. Inst. Mech. Eng., 1: 33–40 (1961)**——sprag-slip 理论源头；对粉笔"推跳"的几何不稳定性有直接参考价值。

### 2.3 湿手指沿杯口 / 玻璃竖琴
- **A. P. French, "In Vino Veritas: A study of wineglass acoustics", American Journal of Physics, 51(8): 688–694 (1983)**——湿手指摩擦激发酒杯振动的经典教学研究，证明其为 stick-slip 自激振动并测量频谱。
- **T. D. Rossing, "Acoustics of the glass harmonica", J. Acoust. Soc. Am., 95(2): 1106–1111 (1994)**——玻璃竖琴声学；湿手指-玻璃界面的 stick-slip 是发声源。
- **S. Derler 等, "Stick–slip phenomena in the friction of human skin", Tribology International（约 2013，卷期页码待核实）**——https://www.sciencedirect.com/science/article/abs/pii/S0043164812003729 用高时间分辨率三向测力台测量指腹在湿玻璃上的粘滑摩擦，给出静摩擦、稳态滑动与粘滑三种摩擦系数的系统对比。这是"湿手指"类系统的定量摩擦学数据，可作为粉笔-黑板摩擦系数测量的实验范式参照。

---

## 3. Stick-Slip 系统的混沌/非线性行为

- **K. Popp & P. Stelter, "Stick-slip vibrations and chaos", Philosophical Transactions of the Royal Society A, 332(1624): 89–105 (1990)**, DOI: 10.1098/rsta.1990.0102——领域奠基论文：数值与实验展示摩擦振子的倍周期通往混沌、多种摩擦模型定性行为一致。粉笔虚线的点间距若出现不规则化，可用该文的非线性动力学语言描述。
- **B. F. Feeny, A. Guran, N. Hinrichs & K. Popp, "A historical review on dry friction and stick-slip phenomena", Applied Mechanics Reviews, 51(5): 321–341 (1998)**, DOI: 10.1115/1.3099008——干摩擦与粘滑现象的历史综述，涵盖从达芬奇到现代非线性动力学的全景。
- **U. Galvanetto, "Some discontinuous bifurcations in a two-block stick-slip system", Journal of Sound and Vibration, 248(4): 653–669 (2001)**（页码待核实）——两滑块系统的非光滑分岔（grazing-sliding 等），展示粘滑切换本身即产生非光滑动力系统特有的分岔类型。
- **M. Denny, "Stick-slip motion: an important example of self-excited oscillation", European Journal of Physics, 25(2): 311–322 (2004)**——教学向论文，把 stick-slip 作为自激振荡范例（负阻尼视角），适合综述中概念框架的表述。
- **（补充）"Bifurcation analysis of stick–slip vibration in a 2-DOF nonlinear dynamical system with dry friction", 中文/国际期刊 2022（ScienceDirect, 待核实具体期刊）**——提出"sticking phase plot"工具，归纳 Border-Collision、Grazing-Sliding 等四类分岔及粘滑混沌，展示该方向仍活跃。

---

## 4. 弹性梁/悬臂粘滑测量方法与高速摄影、声学测量

- **C. M. Mate, G. M. McClelland, R. Erlandsson & S. Chiang, "Atomic-scale friction of a tungsten tip on a graphite surface", Physical Review Letters, 59(17): 1942–1945 (1987)**——AFM/FFM 纳米摩擦测量的开山之作：悬臂扭转测横向力，首次观测到原子尺度周期性 stick-slip。粉笔虚线研究中的"弹性悬臂+接触点粘滑"思想在微尺度上的对应范式。
- **M. Evstigneev et al., "Rate description of the stick-slip motion in friction force microscopy experiments", Physical Review E, 71: 056119 (2005)**, DOI: 10.1103/PhysRevE.71.056119——用热激活逃逸率理论描述悬臂弹簧力的锯齿振荡，给出平均摩擦力-拉速关系的解析近似；适合用于粉笔尖端热/速率效应讨论的模型参照。
- **S. M. Rubinstein, G. Cohen & J. Fineberg, "Detachment fronts and the onset of friction", Nature, 430: 1005–1009 (2004)**, DOI: 10.1038/nature02830——用高速摄影实时观测接触界面真实接触面积与脱黏前沿（detachment fronts）的传播，证明宏观滑动由界面波状破裂前缘介导。这是"高速摄影用于摩擦跳动/粘滑界面观测"的标杆方法学文献。
- **H. Shlomai & J. Fineberg, "The structure of slip-pulses and supershear ruptures driving slip in bimaterial friction", Nature Communications, 7: 11787 (2016)**——高速同步测量滑移速度、真实接触面积与应力，直接观测 slip pulse（滑移脉冲）与超剪切破裂；为粉笔尖端周期性脱黏/冲击的界面波动图像提供实验模板。
- **A. Djellouli, K. Bertoldi 等（Harvard/Nottingham/CNRS 合作）, 关于橡胶-玻璃滑动尖叫（squeak）由超声速脱黏脉冲（slip pulses）介导的研究, Nature（约 2025–2026 发表；确切标题/卷期 DOI 待核实）**——参见 Harvard MRSEC 报道：https://www.mrsec.harvard.edu/pages/news-2026-nsf-mrsec-Physics-of-a-squeak.php 。用最高 100 万 fps 全内反射高速成像+同步声学测量，证明球鞋/软界面尖叫不是随机 stick-slip 而是周期性传播的滑移脉冲，频率由脉冲重复率决定。与粉笔尖叫-虚线的对应性极强：点间距≈脉冲重复率×画线速度，是方法论（高速成像+声学同步）上最直接的现代样板。
- **M. R. Reeves et al., "A study of brake disc modal behaviour during squeal generation using high-speed electronic speckle pattern interferometry and near-field sound pressure measurements", Proc. IMechE Part D, 214: 285–296 (2000)**——高速 ESPI+近场声压联合测量刹车尖叫模态，示范"振动场光学测量+声学"联用方案，可移植到粉笔-黑板实验。

---

## 备注与不确定项
1. IPT 2023 题目名为 "The Chalk Trick"；IYPT 正式赛题中未检索到粉笔虚线题目（不排除某年开放题中出现，待核实）。
2. Djellouli/Bertoldi 团队的 squeak 论文据 MRSEC 报道发表于 Nature，具体卷期与 DOI 需进一步核实。
3. McIntyre & Woodhouse (1979)、Galvanetto (2001)、Derler 等的卷期页码标注了"待核实"，引用前请核对。
4. Painlevé 原始文献：P. Painlevé, "Sur les lois du frottement de glissement", C. R. Acad. Sci. 121: 112–115 (1895)——历史源头，可列入。
5. 中文教科书层面已有讨论：如《Mechanics, Waves and Thermodynamics: An Example-based Approach》（Cambridge, 2016）第 7.3 节专门讨论粉笔尖叫的 stick-slip 分析（含 F = μN/(1−μtanθ) 发散条件），可作为教学化表述参照。
