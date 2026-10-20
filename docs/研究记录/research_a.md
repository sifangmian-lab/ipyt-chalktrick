# 文献调研简报：粉笔黑板虚线（Chalk Trick）的 Stick-Slip 机理——理论基础文献

调研日期：2026（以 web_search 核查为准）。以下文献均经网络检索交叉核实（卷期页码/DOI 与多个独立引用源一致）；个别无法完全确认的细节已标注"待核实"。

---

## 主题一：Stick-Slip 经典理论（Rabinowicz、Bowden–Tabor 黏着摩擦、粗糙度与塑性变形）

1. **Rabinowicz, E. "The Intrinsic Variables Affecting the Stick-Slip Process." Proceedings of the Physical Society, 71(4): 668–675, 1958. DOI: 10.1088/0370-1328/71/4/316.**
   - 要点：系统测定黏滑振幅与加载刚度、驱动速度、摩擦副材料的关系，确立了"静摩擦高于动摩擦 + 系统弹性储能"这一黏滑判据的实验基础。粉笔虚线正是"黏"（粉笔尖端静止、弹性弯曲储能）与"滑"（突然释放、快速跳过一段）交替的直接表现，该文提供了用弹簧刚度—速度参数描述虚线节距的定量框架。

2. **Rabinowicz, E. "Stick and Slip." Scientific American, 194(5): 109–119, 1956. DOI: 10.1038/scientificamerican0556-109.**
   - 要点：面向大众的黏滑综述，明确指出小提琴琴弦发声、轴承尖叫、黑板吱吱声（blackboard squeak）等现象均源于黏滑自激振动。这是把"粉笔/黑板"问题明确归入黏滑范式的直接文献依据，适合作为综述引言引用。

3. **Bowden, F. P. & Tabor, D. The Friction and Lubrication of Solids. Oxford: Clarendon Press, 1950（第 1 卷；1986/2001 有 Oxford Classic Texts 再版）。**
   - 要点：经典专著，提出黏着摩擦理论：真实接触只发生在微凸体顶端，接触点塑性变形形成"结"（junction），摩擦力 = 结面积 × 界面剪切强度，解释了 Amontons 定律（μ 与表观面积无关）。对粉笔问题：粉笔（石膏/CaCO₃）尖端与黑板的真实接触面积决定黏着强度，进而决定"黏"阶段能储存多少弹性能。

4. **Bowden, F. P. & Tabor, D. "The Area of Contact between Stationary and between Moving Surfaces." Proceedings of the Royal Society A, 169: 391–413, 1939. DOI: 10.1098/rspa.1939.0005.**
   - 要点：用接触电阻等方法实测静止与滑动状态下的真实接触面积，证明滑动时接触结在数十毫秒内破坏与重建。支撑"黏"阶段接触结随时间生长、"滑"阶段接触结断裂的图像——这正是虚线分段形成（黏→结生长→断裂→跳移）的微观机制。

5. **Greenwood, J. A. & Williamson, J. B. P. "Contact of Nominally Flat Surfaces." Proceedings of the Royal Society A, 295(1442): 300–319, 1966. DOI: 10.1098/rspa.1966.0242.**
   - 要点：GW 统计接触模型：将表面视为高度随机分布的微凸体集合，证明即使纯弹性接触，真实接触面积也近似正比于载荷；提出区分弹性/塑性接触的"塑性指数"。黑板漆面与粉笔均为粗糙表面，该模型可用于估算不同按压力下粉笔—黑板的接触点数与黏着总强度，从而预测虚线节距随压力的变化。

6. **Brace, W. F. & Byerlee, J. D. "Stick-Slip as a Mechanism for Earthquakes." Science, 153(3739): 990–992, 1966. DOI: 10.1126/science.153.3739.990.**
   - 要点：证明岩石断层的黏滑失稳可作为地震机制，把黏滑从工程摩擦问题提升为普适的"弹性储能—摩擦失稳"现象。说明同一套弹簧-滑块物理既解释地震复发间隔也解释粉笔虚线节距，可作为跨尺度普适性的论据。

---

## 主题二：摩擦系数的精细描述（静摩擦的时间老化、动摩擦的速度依赖、Coulomb 模型的局限）

1. **Dieterich, J. H. "Time-Dependent Friction and the Mechanics of Stick-Slip." Pure and Applied Geophysics, 116(4-5): 790–806, 1978. DOI: 10.1007/BF00876539.**
   - 要点：实验发现静摩擦强度随静止接触时间对数增长（aging/creep 机制：微接触点蠕变使真实面积增大），并提出动摩擦的速度依赖本构。对粉笔问题：虚线中每个"黏点"停留时间越长，下一次滑移所需的力越大、跳距越远——该文为节距—频率关系提供本构依据。

2. **Dieterich, J. H. "Modeling of Rock Friction: 1. Experimental Results and Constitutive Equations." Journal of Geophysical Research, 84(B5): 2161–2168, 1979. DOI: 10.1029/JB084iB05p02161.**
   - 要点：提出含状态变量的速率-状态摩擦律（rate-and-state friction）雏形：μ 取决于滑动速度与一个表征接触"年龄"的内部状态变量。这是描述粉笔尖端在快速滑移—静止循环中摩擦演化的现代标准本构。

3. **Ruina, A. "Slip Instability and State Variable Friction Laws." Journal of Geophysical Research, 88(B12): 10359–10370, 1983. DOI: 10.1029/JB088iB12p10359.**
   - 要点：系统化速率-状态摩擦律（Dieterich 律与 Ruina 律），给出稳态滑动的线性稳定性判据：稳态速度弱化（a−b<0）+ 弹簧刚度低于临界刚度 k_crit 时出现黏滑失稳。可直接用于推导粉笔出现虚线（黏滑）而非实线（稳定滑动）的参数条件。

4. **Heslot, F., Baumberger, T., Perrin, B., Caroli, B. & Caroli, C. "Creep, Stick-Slip, and Dry-Friction Dynamics: Experiments and a Heuristic Model." Physical Review E, 49(6): 4973–4988, 1994. DOI: 10.1103/PhysRevE.49.4973.**
   - 要点：纸-玻璃多接触界面实验表明，即使名义"黏"的阶段界面也在缓慢蠕变（creep），静摩擦的增强源于蠕变导致的真实接触面积增长；并提出把接触面视为独立微结集合的启发式模型。对粉笔：说明"黏点"并非绝对静止而是亚毫米级蠕变，这决定虚线边缘的模糊/堆积形貌。

5. **Baumberger, T. & Caroli, C. "Solid Friction from Stick-Slip Down to Pinning and Ageing." Advances in Physics, 55(3-4): 279–348, 2006. DOI: 10.1080/00018730600732186.**
   - 要点：权威综述，统一讨论从宏观黏滑到微观钉扎与老化的干摩擦物理，明确指出 Coulomb 双系数模型（恒定 μs > μk）只是低速、低加速度下的粗糙近似，无法解释蠕变、记忆效应与速度依赖。可作为"Coulomb 模型局限性"一节的支柱文献。

6. **Berthoud, P., Baumberger, T., G'Sell, C. & Hiver, J.-M. "Physical Analysis of the State- and Rate-Dependent Friction Law: Static Friction." Physical Review B, 59(22): 14313–14327, 1999. DOI: 10.1103/PhysRevB.59.14313.**
   - 要点：定量验证静摩擦随停留时间的对数老化律，并把状态变量的物理起源归因于多接触界面的塑性蠕变（活化体积与材料流变相关）。为"粉笔停留时间 → 黏着力 → 虚线间隔"链条提供材料学层面的微观解释。

---

## 主题三：Stick-Slip 建模（弹簧-滑块-传送带模型、摩擦自激振动、非线性/混沌行为）

1. **Wiercigroch, M. "A Note on the Switch Function for the Stick-Slip Phenomenon." Journal of Sound and Vibration, 175(5): 700–704, 1994. DOI: 10.1006/jsvi.1994.1352（待核实，卷期页码已由多个引用源确认）。**
   - 要点：给出用光滑开关函数（switch function）替代静/动摩擦分段切换的建模技巧，使黏滑系统可用连续动力学工具分析。对粉笔虚线模拟（数值积分黏→滑切换时的奇异性处理）有直接方法论价值。

2. **Popp, K. & Stelter, P. "Stick-Slip Vibrations and Chaos." Philosophical Transactions of the Royal Society A, 332: 89–105, 1990. DOI: 待核实。**
   - 要点：对弹簧-滑块类干摩擦振子的黏滑振动做非线性分析，展示倍周期分岔通往混沌的路径，说明虚线节距可以是不规则（混沌）而非严格周期的。可解释实际粉笔画出的虚线间隔不均匀现象。

3. **Wiercigroch, M. & de Kraker, B. (Eds.). Applied Nonlinear Dynamics and Chaos of Mechanical Systems with Discontinuities. World Scientific, Series on Nonlinear Science A, Vol. 28, 2000.**
   - 要点：含摩擦、碰撞等不连续机械系统非线性动力学的论文集，涵盖黏滑极限环、擦边分岔（grazing bifurcation）等。提供将"粉笔—手—黑板"视为不连续振子进行分岔分析的理论工具箱。

4. **Baumberger, T., Caroli, C., Perrin, B. & Ronsin, O. "Nonlinear Analysis of the Stick-Slip Bifurcation in the Creep-Controlled Dry Friction Regime." Physical Review E, 51(5): 4005–4010, 1995. DOI: 10.1103/PhysRevE.51.4005（DOI 后缀按 APS 惯例推断，页码已由引用源确认，DOI 待核实）。**
   - 要点：在蠕变控制 regime 下对黏滑分岔做弱非线性（振幅方程）分析，证明失稳是超临界/亚临界 Hopf 型分岔并给出黏滑周期在分岔点附近发散的标度律。可用于定量拟合"画线速度接近临界值时虚线周期急剧变长"的实验现象。

5. **Awrejcewicz, J. & Olejnik, P. "Occurrence of Stick-Slip Phenomenon." Journal of Theoretical and Applied Mechanics, 45(1): 33–40, 2007.**
   - 要点：简明综述黏滑的发生条件（速度弱化摩擦、足够柔的驱动、低速），并给出弹簧-滑块-传送带基准模型的相图。适合作为综述中"建模方法"一节的入门引文。

6. **Marone, C. "Laboratory-Derived Friction Laws and Their Application to Seismic Faulting." Annual Review of Earth and Planetary Sciences, 26: 643–696, 1998. DOI: 10.1146/annurev.earth.26.1.643.**
   - 要点：总结速率-状态摩擦律下弹簧-滑块模型的全部稳定性结果（k 与 k_crit 比较、临界滑移距离 Dc），并讨论周期性黏滑、慢滑移与稳定滑动的参数区域。粉笔虚线节距类比于地震复发滑移量，该文提供 k_crit = −(a−b)σ_n/D_c 的定量估算式。

---

## 主题四：弹性理论与梁振动在摩擦系统中的应用（简述）

1. **Landau, L. D. & Lifshitz, E. M. Theory of Elasticity（《弹性理论》，Course of Theoretical Physics, Vol. 7）. Pergamon Press, 3rd ed., 1986（俄文原版更早；中译本由高等教育出版社出版）。**
   - 要点：提供弹性体应力-应变、弯曲（梁/板挠曲）与接触问题的基本方程。粉笔在"黏"阶段的弯曲储能可用细长杆弯曲（E-B 梁）近似估算，其等效弹簧刚度 k = 3EI/L³ 直接代入黏滑稳定性判据。

2. **Timoshenko, S., Young, D. H. & Weaver, W. Vibration Problems in Engineering. Wiley, 4th ed., 1974.**
   - 要点：Euler-Bernoulli 梁自由/受迫振动的标准推导与固有频率公式。粉笔尖端受摩擦力扰动的横向振动模态决定黏滑释放时的动态响应与可闻"吱吱"声频（数 kHz 量级），可用梁振动基频估算。

3. **Johnson, K. L. Contact Mechanics. Cambridge University Press, 1985/1987. DOI: 10.1017/CBO9781139171731.**
   - 要点：Hertz 接触及切向加载（Cattaneo–Mindlin 部分滑移理论）的经典教材，给出黏着区/微滑区在切向力作用下的演化。粉笔尖端与黑板的接触可近似为球-面 Hertz 接触，全黏滑移的临界切向力即对应"黏→滑"转换条件。

4. **Rice, J. R. & Ruina, A. L. "Stability of Steady Frictional Slipping." Journal of Applied Mechanics, 50(2): 343–349, 1983. DOI: 10.1115/1.3167042.**
   - 要点：把速率-状态摩擦本构与弹性加载系统耦合，解析给出稳定滑动的临界弹簧刚度公式。这是"弹性元件（粉笔杆身/梁）+ 摩擦界面"联合建模的典范，直接把弹性理论与摩擦失稳联系起来。

---

## 备注与待核实事项
- 未找到专门针对"chalk trick / 黑板粉笔虚线"的同行评议论文（搜索 "chalk blackboard stick-slip American Journal of Physics" 无果）；教学层面多见于物理演示资料（如 Walker《The Flying Circus of Physics》提及黑板尖叫），建议综述中以 Rabinowicz 1956（明确提及黑板吱吱声）作为最接近的一手引文。
- Popp & Stelter (1990) 的卷号 332、页码 89–105 由 arXiv:2405.13142 参考文献列表确认，DOI 未核实。
- Wiercigroch (1994) 期刊卷期页码（JSV 175(5):700–704）已由多个独立引用源确认，DOI 按 Elsevier 旧格式推断，标注待核实。
