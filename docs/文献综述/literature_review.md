# 粉笔-黑板虚线现象（Chalk Trick）的黏滑机理：文献综述

## 1. 引言

以一支粉笔在黑板上快速画出一条由断续斑点构成的"虚线"，是流传已久的物理演示。麻省理工学院 Walter Lewin 教授的授课视频使这一现象广为人知：当粉笔与运动方向成较大（钝角）倾角、并以"推"而非"拉"的方式施压时，即使手部匀速移动，粉笔也会在黑板上留下等间距的点状痕迹。国际物理学家锦标赛（IPT）2023 年将之列为正式赛题 "The Chalk Trick"，要求参赛者解释虚线的形成机理，并探讨能否从痕迹反推粉笔-黑板之间的相对运动参数乃至粉笔的几何尺寸 [1]。本文综述的源研究 [2]（IPT 参赛报告 2903.pdf）采用黏滑（slip-stick）模型，以弹性橡胶梁模拟人手的柔顺性，将粉笔运动分解为静止（黏）、跳跃与画线三个阶段循环，用 Euler-Bernoulli 梁振动理论计算跳距，并辅以实验测得的静、动摩擦系数（$\mu_s=0.62\pm0.03$、$\mu_k=0.51\pm0.02$）与梁参数标定，对虚线节距给出了定量预测。

尽管该现象演示效果强烈、物理内涵丰富，但令人意外的是，专门针对粉笔虚线的同行评议系统研究至今几近空白：除 Rabinowicz 1956 年的科普性论述 [3] 和 Champneys 与 Várkonyi 2016 年从 Painlevé 悖论角度给出的理论分析 [4] 外，系统工作主要见于 IPT 参赛报告 [2][5]。因此，对这一现象的学术定位必须借助于黏滑与摩擦诱导振动这一更为成熟的研究传统。本综述的目的在于：（i）梳理黏滑现象的经典理论与摩擦本构的精细描述；（ii）总结黏滑系统的建模与非线性动力学成果；（iii）考察琴弦、刹车尖叫、湿手指等典型摩擦发声系统对粉笔问题的类比价值；（iv）评估粉笔-黑板问题的研究现状与不同建模路线之间的张力；（v）归纳可用的测量方法，并指出研究空白与未来方向。

## 2. 黏滑现象的经典理论基础

黏滑（stick-slip）是指两固体接触界面在持续驱动下交替发生"黏着"与"突然滑移"的自激振动。其经典判据由 Rabinowicz 系统确立：当静摩擦系数高于动摩擦系数（或摩擦随滑动速度下降），且驱动系统具备足够的弹性储能能力时，稳态滑动失稳而转为黏滑振荡 [3][6]。Rabinowicz 在 1956 年的著名综述中即以粉笔痕迹为插图：粉笔轴与运动方向成锐角时画出连续线，成钝角时则发生"黏—滑—黏"循环留下点划线，且握得越紧滑移距离越短 [3]。这几乎是粉笔虚线问题最早的正式文献锚点，也直接启发了后续研究以"角度"与"按压力"作为关键控制参数。

黏滑判据的微观基础来自黏着摩擦理论。Bowden 与 Tabor 证明，真实接触仅发生在微凸体顶端并形成塑性接触结（junction），摩擦力等于结的总面积乘以界面剪切强度，由此解释了 Amontons 定律 [7]；他们进一步的接触电阻测量表明，滑动时接触结在数十毫秒内破坏并重建，为"黏"阶段结生长、"滑"阶段结断裂的循环图像提供了直接证据 [8]。Greenwood 与 Williamson 的统计接触模型则证明，即使纯弹性接触，真实接触面积也近似正比于法向载荷，并给出区分弹性/塑性接触的塑性指数 [9]。对粉笔问题而言，这些结果意味着：粉笔（石膏/碳酸钙）尖端与黑板漆面均为粗糙表面，"黏"阶段能储存多少弹性能取决于真实接触面积与结强度，而虚线节距随按压力的变化可由接触点数的统计估算来预测。值得强调的是，Brace 与 Byerlee 将同一"弹性储能—摩擦失稳"机制用于解释地震复发 [10]，表明弹簧-滑块物理跨越十余个数量级的普适性——粉笔虚线节距与地震滑移量在模型层面是同构的，这为借用地震学中成熟的速率-状态摩擦理论（见第 3 节）研究粉笔问题提供了合法性。在接触力学层面，粉笔尖端与黑板的单点接触可近似为 Hertz 接触，Cattaneo-Mindlin 部分滑移理论给出切向加载下黏着区/微滑区的演化，全黏滑移的临界切向力即对应"黏→滑"转换条件 [48]，这为细化源报告的整体 Coulomb 假设提供了标准工具。

## 3. 摩擦本构的精细描述

源报告 [2] 采用双系数 Coulomb 模型（恒定 $\mu_s>\mu_k$），并明确承认其"过于简单"，将静摩擦的时间依赖与动摩擦的速度依赖列为改进方向。这一自评与摩擦学界半个世纪的进展高度吻合。

Dieterich 的岩石摩擦实验发现，静摩擦强度随静止接触时间呈对数增长（老化，aging），其机制是微接触点的蠕变使真实接触面积增大 [11]；他随后提出含内部状态变量的速率-状态摩擦律（rate-and-state friction），即摩擦系数依赖于滑动速度与一个表征接触"年龄"的状态变量 [12]。Ruina 将该本构系统化，并给出稳态滑动的线性稳定性判据：当稳态摩擦呈速度弱化（$a-b<0$）且加载弹簧刚度低于临界刚度 $k_{crit}$ 时发生黏滑失稳 [13]；Rice 与 Ruina 进一步把速率-状态本构与弹性体耦合，给出临界刚度的解析公式 [14]，Marone 则总结了 $k_{crit}=-(a-b)\sigma_n/D_c$ 的估算式及周期黏滑、慢滑移与稳定滑动的参数分区 [15]。对粉笔问题，这套理论可直接回答"何时画虚线、何时画实线"这一核心问题：虚线对应于系统柔度超过临界值的失稳区，而画线速度、按压力（正应力 $\sigma_n$）与粉笔-手臂系统的等效刚度共同决定系统位于哪个参数区。源报告观测到虚线仅出现在 10°–25° 倾角范围内 [2]，正是这类稳定性边界的实验体现。

"黏"阶段的精细图像则由 Heslot 等的纸-玻璃实验给出：即使名义静止阶段，界面也在缓慢蠕变，静摩擦的增强源于蠕变导致的真实接触面积增长 [16]；Berthoud 等定量验证了静摩擦随停留时间的对数老化律，并将其物理起源归于多接触界面的塑性蠕变 [17]。这提示粉笔的"黏点"并非绝对静止，而是存在亚毫米级蠕变——它决定虚线斑点边缘的堆积形貌，也是"从痕迹反推运动参数"时必须考虑的系统误差。Baumberger 与 Caroli 的权威综述统一讨论了从宏观黏滑到微观钉扎与老化的干摩擦物理，明确指出 Coulomb 双系数模型只是低速、低加速度下的粗糙近似，无法解释蠕变、记忆效应与速度依赖 [18]。

与此形成对照的是 Ballard 的分析 [19]：在纯 Coulomb 摩擦定律与线弹性耦合的框架内，稳态滑动解随摩擦系数变化可出现分岔（解失去唯一性），这种分岔本身即可解释某些黏滑事件，而无需引入任何更精细的摩擦律。这一视角对粉笔问题具有方法论意义：源报告观测到的黏-滑切换，究竟源于摩擦系数的双值性（本构假设），还是源于弹性-接触系统的分岔（结构性质），在原则上是可以通过改变梁刚度与摩擦参数加以区分的——这也说明 Coulomb 模型虽然在定量上粗糙，却未必在定性上错误。

## 4. 黏滑系统的建模与非线性动力学

弹簧-滑块-传送带（spring-block-belt）是黏滑建模的范式模型 [20]：滑块经弹簧与恒速驱动端相连，弹簧储能至摩擦力阈值后滑块突跳；其从单体到链式多振子的推广（如 Coulomb 摩擦谐振子链 [49]）同样呈现丰富的黏滑动力学。源报告 [2] 的橡胶梁-粉笔系统正是该范式的连续体版本，其以悬臂梁等效刚度 $k=3EI/L^3$ 代入黏滑条件、并用梁的自由振动方程（特征方程 $\cosh(\lambda\sqrt{\omega})\cos(\lambda\sqrt{\omega})+1=0$）计算跳跃相历时，实质上把"手"的柔顺性归因于一根Euler-Bernoulli 梁——这与 Landau-Lifshitz [21] 和 Timoshenko 等 [22] 的经典弹性体振动理论一脉相承。数值处理上，Wiercigroch 提出的光滑开关函数可替代静/动摩擦的分段切换，消除黏-滑转换时刻的数值奇异性 [23]——源报告亦引用了该文，而含摩擦、碰撞的不连续系统更有系统的非光滑动力学工具箱（极限环、擦边分岔等）可用 [24]。

黏滑系统的节距不必是严格周期的。Popp 与 Stelter 的奠基性工作表明，干摩擦振子可经倍周期分岔走向混沌，且多种摩擦模型给出定性一致的行为 [25]；两滑块系统还呈现 grazing-sliding 等非光滑分岔 [26]。Denny 从负阻尼角度将黏滑表述为自激振荡的教学范例 [27]，Feeny 等的历史综述则给出了从达芬奇到现代非线性动力学的全景 [28]。源报告在实验中确实观测到"所有表观参数相同而行为不同"的异常轨迹，并谨慎地将之归因于文献 [25] 报道的混沌行为 [2]——这为虚线节距的不规则起伏提供了现成的动力学语言，也提示"从痕迹反推参数"在混沌区可能因病态反问题而失效。

更具张力的文献分歧来自刚性接触力学。Painlevé 早在 1895 年即指出 Coulomb 摩擦的刚体动力学可发生解的非唯一或不存在（Painlevé 悖论）[29]。Leine、Brogliato 与 Nijmeijer 证明摩擦碰撞振子的跳动（hopping/bouncing）周期运动与该悖论密切相关，并在现实摩擦系数下即可出现 [30]。Champneys 与 Várkonyi 的权威综述更明确以"推粉笔"为核心实例，主张粉笔虚线并非纯粹的切向黏滑，而是每个循环包含离面（lift-off）与碰撞的 sprag-slip 振荡，可由 Painlevé 悖论经"reverse chatter"失稳产生 [4]。这与源报告的纯黏滑图像存在实质差异：在 sprag-slip 图像中，粉笔与黑板在"滑"阶段完全脱离接触（源报告描述的"粉笔跳起来越过一段距离再接触"其实已隐含此成分），接触几何（推 vs 拉、锐角 vs 钝角）通过法向-切向耦合而非仅通过摩擦系数双值性起作用。两种机制并不互斥——在柔性加载下二者可叠加——但区分它们需要离面位移的直接测量（高速摄影，见第 7 节），这是当前文献中最值得澄清的张力点。

## 5. 典型摩擦诱导振动系统的类比

粉笔问题并非孤例，它与多类摩擦诱导振动系统共享物理机制，可从其成熟理论中借鉴。

**小提琴琴弦。** 弓弦系统是研究最透彻的黏滑自激振动：Raman 最早给出 Helmholtz 运动（尖角沿弦往返、每周期一黏一滑）的定量理论 [31]；Schelleng 给出最小/最大弓压力与弓速的参数图，力不足出现多次滑移噪声、力过大则音粗，这与粉笔问题中"压力-速度平面内不同运动模式分区"完全对应 [32]；Cremer 系统处理了摩擦黏/滑切换在波传播中对 Helmholtz 尖角的"再锐化"机制 [33]；McIntyre 与 Woodhouse 建立的反射波逐次映射方法成为现代弓弦模拟基准 [34]。差异在于：琴弦是分布参数波动系统，黏滑由行波反馈维持；粉笔-梁系统更接近集中参数振子（至少在一阶模态近似下），因此弓弦文献的"参数分区图"思想比其具体方程更具可移植性。

**刹车尖叫。** 刹车尖叫与粉笔"推跳"在几何上同构：Spurr 提出的 sprag-slip 机制 [35]——锐角接触使摩擦力的切向分量耦合出转动/法向自由度，等效于一个负刚度失稳——正是 Champneys 与 Várkonyi 用于粉笔的机制 [4]。Kinkaid 等的权威综述梳理了黏滑、sprag-slip、负摩擦-速度斜率与模态耦合四类机制 [36]，Papinniemi 等则强调没有单一机制能解释全部尖叫事件 [37]。这一教训对粉笔研究尤为重要：仅凭虚线痕迹反推机理，可能面临"多机制产生相似痕迹"的反问题欠定性。

**湿手指与皮肤摩擦。** French 证明湿手指激发酒杯振动是黏滑自激振动并测其频谱 [38]，Rossing 对玻璃竖琴的研究类似 [39]；Derler 等用高时间分辨率三向测力台给出指腹在湿玻璃上静摩擦、稳态滑动与黏滑三区的系统摩擦数据 [40]，为粉笔-黑板摩擦系数的分 regime 测量提供了实验范式。Akay 的摩擦声学综述将粉笔尖叫、湿手指、琴弦统一置于"摩擦发声"框架下 [41]，是从声学角度组织类比的入口。

## 6. 粉笔-黑板问题的研究现状

直接针对粉笔虚线的文献可归纳为四条路线。**其一**，Rabinowicz [3] 最早以粉笔痕迹说明黏滑，给出"钝角→断续、握力→节距"的定性规律，但未建模。**其二**，Champneys 与 Várkonyi [4]（及 Leine 等 [30]）的 sprag-slip/Painlevé 路线，强调离面碰撞与接触几何，属刚性接触力学框架，尚无可与节距数据直接对比的定量预测。**其三**，IPT 2023 参赛报告：源报告 [2] 以橡胶梁模拟人手的三阶段循环模型，用梁振动周期决定跳距、以静摩擦阈值决定黏相时长，预测节距与画线速度、倾角、梁刚度的关系，并与实验定性符合；Sundström 报告 [5] 则采用更简洁的一维弹簧-变摩擦力模型，拟合出节距与画线速度的强线性关系（Pearson 系数 0.99），发现弹簧刚度单调影响节距而法向力影响不显著，并建造机械装置复现虚线。比较可见：源报告的连续梁模型物理内容更丰富（能给出跳跃相历时与振动频率），但 Sundström 的简约模型在数据拟合上更成功；两者都采用 Coulomb 型双系数摩擦，均未纳入时间老化与速度弱化，也都未检验 sprag-slip 假说。三条建模路线的并存，恰恰说明该现象尚缺一个能统一解释角度窗口、节距-速度标度与痕迹形貌的公认理论。值得注意的是，在蠕变控制的摩擦 regime 下，黏滑失稳被证明是 Hopf 型分岔，且黏滑周期在分岔点附近按标度律发散 [42]——这为"画线速度接近临界值时虚线周期急剧变长"的实验现象提供了定量拟合工具。教学文献中已有对粉笔尖叫的简化分析（如含 $F=\mu N/(1-\mu\tan\theta)$ 发散条件的讨论 [50]），其发散条件与 sprag-slip 的几何失稳一脉相承，可作为连接两条路线的中间环节。

## 7. 测量与实验方法

可靠的测量是区分上述机制的瓶颈。**摩擦学测量**方面，源报告用斜面法测得 $\mu_s$、$\mu_k$ [2]，属常规手段；更精细的分 regime 测量可参照皮肤摩擦的三向测力台方案 [40]。**梁参数标定**方面，源报告以附加质量-位移的静力标定获得 $EI$ [2]，其理论基础即经典梁弯曲公式 [21]。**高速摄影**已成为黏滑界面观测的标杆方法：Rubinstein、Cohen 与 Fineberg 用高速摄影实时观测真实接触面积与脱黏前沿（detachment fronts）的传播，证明宏观滑移由界面破裂前缘介导 [43]；Shlomai 与 Fineberg 进一步同步测量滑移速度、接触面积与应力，直接观测滑移脉冲与超剪切破裂 [44]。近期关于软界面尖叫由周期传播的滑移脉冲介导、并以百万帧每秒成像与声学同步测量的工作（文献信息待核实）与粉笔"尖叫-虚线"的对应性极强：点间距近似等于脉冲重复率乘以画线速度。**AFM/摩擦力显微镜**方面，Mate 等首次观测原子尺度周期黏滑 [45]，Evstigneev 等用热激活逃逸率理论描述悬臂锯齿振荡 [46]——"弹性悬臂+接触点黏滑"在微尺度上的对应范式，可为粉笔尖端力学提供标度类比。**声学方法**方面，Akay 奠定了摩擦声学的分析框架 [41]，Reeves 等以高速电子散斑干涉与近场声压联合测量刹车尖叫模态 [47]，示范了"光学振动场测量+声学"联用方案。源报告建议以粉笔敲击黑板的声响作为测量手段 [2]——文献表明这完全可行：黏滑事件的时间戳可从声脉冲序列中提取（黏滑频率达数 kHz 量级，可用梁基频估算 [22]），与高速摄影同步即可分别标定黏相与跳相，是低成本且信息量大增的改进路径。

## 8. 研究空白与展望

综合上述文献，可以明确指出：**目前尚无专门针对粉笔-黑板虚线现象的同行评议系统研究**——最接近的正式出版物是 Rabinowicz 的科普论述 [3] 与 Champneys-Várkonyi 的理论分析 [4]，系统的实验-模型对照仅见于 IPT 参赛报告 [2][5]。据此提出以下未来方向：（i）**本构升级**：将速率-状态摩擦律（含时间老化与速度弱化 [11][12][13]）引入粉笔-梁模型，利用 $k_{crit}$ 判据定量预测虚线/实线的参数边界，解释角度窗口的来源；（ii）**机制判别**：通过离面位移的直接高速摄影观测，判别纯黏滑与 sprag-slip/reverse chatter [4][30] 的相对贡献；（iii）**非线性行为**：系统考察节距的倍周期化与混沌 [25]，明确反演问题的适定参数区；（iv）**声学反演**：发展以声脉冲序列反推黏滑时间结构的方法 [41][47]，并检验"点间距=重复率×速度"的脉冲图像；（v）**反问题的统计学**：将"从痕迹反推运动参数与粉笔尺寸"表述为带模型误差的反问题，量化 Coulomb 近似对反演结果的偏差。粉笔-黑板系统以其低成本、强可视性与丰富的动力学行为，有望成为黏滑与非光滑动力学教学与研究的标准模型系统。

## 参考文献

[1] International Physicists' Tournament. Problem "The Chalk Trick", IPT 2023 [EB/OL]. https://ipt.science/media/model%20reports/2903.pdf.

[2] The Chalk Trick Pre-selection Report [R/OL]. IPT model report 2903, 2023. https://ipt.science/media/model%20reports/2903.pdf.

[3] Rabinowicz E. Stick and slip [J]. Scientific American, 1956, 194(5): 109-118.

[4] Champneys A R, Várkonyi P L. The Painlevé paradox in contact mechanics [J]. IMA Journal of Applied Mathematics, 2016, 81(3): 538-588.

[5] Sundström V. The Chalk Trick [C]// I/IPT Conference Proceedings 2023–2024. Zenodo, 2024. https://zenodo.org/records/10819172.

[6] Rabinowicz E. The intrinsic variables affecting the stick-slip process [J]. Proceedings of the Physical Society, 1958, 71(4): 668-675.

[7] Bowden F P, Tabor D. The Friction and Lubrication of Solids [M]. Oxford: Clarendon Press, 1950.

[8] Bowden F P, Tabor D. The area of contact between stationary and between moving surfaces [J]. Proceedings of the Royal Society A, 1939, 169: 391-413.

[9] Greenwood J A, Williamson J B P. Contact of nominally flat surfaces [J]. Proceedings of the Royal Society A, 1966, 295(1442): 300-319.

[10] Brace W F, Byerlee J D. Stick-slip as a mechanism for earthquakes [J]. Science, 1966, 153(3739): 990-992.

[11] Dieterich J H. Time-dependent friction and the mechanics of stick-slip [J]. Pure and Applied Geophysics, 1978, 116(4-5): 790-806.

[12] Dieterich J H. Modeling of rock friction: 1. Experimental results and constitutive equations [J]. Journal of Geophysical Research, 1979, 84(B5): 2161-2168.

[13] Ruina A. Slip instability and state variable friction laws [J]. Journal of Geophysical Research, 1983, 88(B12): 10359-10370.

[14] Rice J R, Ruina A L. Stability of steady frictional slipping [J]. Journal of Applied Mechanics, 1983, 50(2): 343-349.

[15] Marone C. Laboratory-derived friction laws and their application to seismic faulting [J]. Annual Review of Earth and Planetary Sciences, 1998, 26: 643-696.

[16] Heslot F, Baumberger T, Perrin B, et al. Creep, stick-slip, and dry-friction dynamics: experiments and a heuristic model [J]. Physical Review E, 1994, 49(6): 4973-4988.

[17] Berthoud P, Baumberger T, G'Sell C, et al. Physical analysis of the state- and rate-dependent friction law: static friction [J]. Physical Review B, 1999, 59(22): 14313-14327.

[18] Baumberger T, Caroli C. Solid friction from stick-slip down to pinning and ageing [J]. Advances in Physics, 2006, 55(3-4): 279-348.

[19] Ballard P. Can the 'stick-slip' phenomenon be explained by a bifurcation in the steady sliding frictional contact problem? [J]. Discrete and Continuous Dynamical Systems Series S, 2016, 9(2): 363-381.

[20] Awrejcewicz J, Olejnik P. Occurrence of stick-slip phenomenon [J]. Journal of Theoretical and Applied Mechanics, 2007, 45(1): 33-40.

[21] Landau L D, Lifshitz E M. Theory of Elasticity [M]. 3rd ed. Oxford: Pergamon Press, 1986.

[22] Timoshenko S, Young D H, Weaver W. Vibration Problems in Engineering [M]. 4th ed. New York: Wiley, 1974.

[23] Wiercigroch M. A note on the switch function for the stick-slip phenomenon [J]. Journal of Sound and Vibration, 1994, 175(5): 700-704.

[24] Wiercigroch M, de Kraker B, eds. Applied Nonlinear Dynamics and Chaos of Mechanical Systems with Discontinuities [M]. Singapore: World Scientific, 2000.

[25] Popp K, Stelter P. Stick-slip vibrations and chaos [J]. Philosophical Transactions of the Royal Society A, 1990, 332(1624): 89-105.

[26] Galvanetto U. Some discontinuous bifurcations in a two-block stick-slip system [J]. Journal of Sound and Vibration, 2001, 248(4) [页码待核实].

[27] Denny M. Stick-slip motion: an important example of self-excited oscillation [J]. European Journal of Physics, 2004, 25(2): 311-322.

[28] Feeny B F, Guran A, Hinrichs N, et al. A historical review on dry friction and stick-slip phenomena [J]. Applied Mechanics Reviews, 1998, 51(5): 321-341.

[29] Painlevé P. Sur les lois du frottement de glissement [J]. Comptes Rendus de l'Académie des Sciences, 1895, 121: 112-115.

[30] Leine R I, Brogliato B, Nijmeijer H. Periodic motion and bifurcations induced by the Painlevé paradox [J]. European Journal of Mechanics A/Solids, 2002, 21(5): 869-896.

[31] Raman C V. On the mechanical theory of the vibrations of bowed strings and of musical instruments of the violin family: Part I [J]. Bulletin of the Indian Association for the Cultivation of Science, 1918, 15: 1-158.

[32] Schelleng J C. The bowed string and the player [J]. Journal of the Acoustical Society of America, 1973, 53(1): 26-41.

[33] Cremer L. The Physics of the Violin [M]. Cambridge, MA: MIT Press, 1984.

[34] McIntyre M E, Woodhouse J. On the fundamentals of bowed-string dynamics [J]. Acustica, 1979, 43: 93-108 [卷期页码待核实].

[35] Spurr R T. A theory of brake squeal [J]. Proceedings of the Automobile Division, Institution of Mechanical Engineers, 1961, 1: 33-40.

[36] Kinkaid N M, O'Reilly O M, Papadopoulos P. Automotive disc brake squeal [J]. Journal of Sound and Vibration, 2003, 267(1): 105-166.

[37] Papinniemi A, Lai J C S, Zhao J, et al. Brake squeal: a literature review [J]. Applied Acoustics, 2002, 63(4): 391-400.

[38] French A P. In vino veritas: a study of wineglass acoustics [J]. American Journal of Physics, 1983, 51(8): 688-694.

[39] Rossing T D. Acoustics of the glass harmonica [J]. Journal of the Acoustical Society of America, 1994, 95(2): 1106-1111.

[40] Derler S, et al. Stick-slip phenomena in the friction of human skin [J]. Tribology International [卷期页码待核实].

[41] Akay A. Acoustics of friction [J]. Journal of the Acoustical Society of America, 2002, 111(4): 1525-1548.

[42] Baumberger T, Caroli C, Perrin B, et al. Nonlinear analysis of the stick-slip bifurcation in the creep-controlled dry friction regime [J]. Physical Review E, 1995, 51(5): 4005-4010.

[43] Rubinstein S M, Cohen G, Fineberg J. Detachment fronts and the onset of friction [J]. Nature, 2004, 430: 1005-1009.

[44] Shlomai H, Fineberg J. The structure of slip-pulses and supershear ruptures driving slip in bimaterial friction [J]. Nature Communications, 2016, 7: 11787.

[45] Mate C M, McClelland G M, Erlandsson R, et al. Atomic-scale friction of a tungsten tip on a graphite surface [J]. Physical Review Letters, 1987, 59(17): 1942-1945.

[46] Evstigneev M, et al. Rate description of the stick-slip motion in friction force microscopy experiments [J]. Physical Review E, 2005, 71: 056119.

[47] Reeves M R, et al. A study of brake disc modal behaviour during squeal generation using high-speed electronic speckle pattern interferometry and near-field sound pressure measurements [J]. Proceedings of the IMechE Part D, 2000, 214: 285-296.

[48] Johnson K L. Contact Mechanics [M]. Cambridge: Cambridge University Press, 1985.

[49] Wensrich C. Slip-stick motion in harmonic oscillator chains subject to Coulomb friction [J]. Tribology International, 2006, 39: 490.

[50] Mechanics, Waves and Thermodynamics: An Example-based Approach [M]. Cambridge: Cambridge University Press, 2016 [作者待核实].
