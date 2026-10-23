# 新参考文献与 N 下降段方向对照笔记

整理日期：2026-10-23
配套模型：《粉笔书写单周期压力下降段的建模与数值模拟》（第一阶段，自由轴刚体棒 + 板面接触柔度 + 库仑摩擦含方向翻转）

## 0. 定位

第一阶段模型骨架：自由轴刚体棒（绕定点 O 转动）+ 板面接触柔度 $k_c$（罚弹簧，$N=k_c\delta$，$\delta=L\cos\theta-h$）+ 库仑摩擦（tanh 正则化方向翻转）+ $\delta=0$ 事件终止。计算对象是单周期黏滑循环中"滑移事件的后半段"：法向力从交接值 $N_0$ 几何卸载至零（分离），中途经历"启动段 → 摩擦翻转 → 滑行段"两阶段，输出分离时间 $t^*$、本段笔痕强度 $I_{2nd}=\int\mu N|v_{rel}|dt$ 与空间强度分布 $I(\xi)$（含末端回划折叠尖峰）。

新加参考文献恰好落在模型的四条"接缝"上，以下逐块对照，并给出每条对照可直接执行的动作。

## 1. 接缝一：摩擦本构（tanh 翻转 → 状态变量）

| 文献 | 存放位置 | 与模型的对应 | 可直接执行的动作 |
|---|---|---|---|
| Åström & Canudas de Wit 2008, *Revisiting the LuGre friction model*, IEEE Control Syst. Mag. 28(6): 101–114 | references/期刊论文/Astrom_CanudasDeWit_2008_Revisiting_LuGre_IEEE.pdf | 鬃毛偏转内态 $z$：$\dot z=v-\sigma_0\|v\|z/g(v)$，$F=\sigma_0 z+\sigma_1\dot z+f(v)$；$g(v)$ 容纳库仑 + Stribeck 效应。替换 tanh 翻转后，"黏"相获得真实切向预滑柔度，摩擦翻转由内态记忆自然涌现 | ① 做 $\varepsilon$ 灵敏度扫描（翻转时刻 $t_{flip}$ 对 $\varepsilon$ 的依赖），评估预言③的稳健性；② 以 LuGre 重写摩擦项，对比翻转前后的 $N(t)$、$I(\xi)$ 差异 |
| González-Carbajal et al. 2024, *The role of dynamic friction in the appearance of periodic oscillations in mechanical systems*, Nonlinear Dyn. 112: 21587–21603 | references/期刊论文/Gonzalez-Carbajal_etal_2024_Role_dynamic_friction_periodic_oscillations_NonlinDyn.pdf | Brown–McPhee 光滑摩擦律（C¹）下，稳定周期轨道可经超临界 Hopf 或极限环鞍结（SN）产生，参数平面分单稳/双稳三区 | 多周期冲击映射拼好后，"虚线能否自持"用 Hopf/SN 分岔语言描述；该律光滑、可微，适合数值延拓 |
| Elmer 1997, *Nonlinear dynamics of dry friction*, arXiv:chao-dyn/9707004 | references/预印本/Elmer_1997_Nonlinear_dynamics_of_dry_friction_arXiv_chao-dyn9707004.pdf | 弹簧-滑块 + "静摩擦随黏着时间增长、动摩擦速度依赖"的轻量参数化；三种吸引子（黏滑/连续滑动/无黏着振荡）与局部分岔、全局分岔全景 | 第一版上升段模型的摩擦律首选（比完整速率-状态律轻量）；为相图提供分类学参照 |

## 2. 接缝二：分离机制（罚弹簧归零 vs 刚性悖论）

| 文献 | 存放位置 | 与模型的对应 | 可直接执行的动作 |
|---|---|---|---|
| Leine, Brogliato & Nijmeijer 2002, *Periodic motion and bifurcations induced by the Painlevé paradox*, Eur. J. Mech. A/Solids 21(5): 869–896 | references/期刊论文/Leine_Brogliato_Nijmeijer_2002_Periodic_motion_bifurcations_Painleve_EJMAS.pdf | 摩擦碰撞振子的跳跃（hopping）与 Painlevé 悖论（解不存在/不唯一）密切相关；dynamic jam 与 reverse chatter 是失稳出口 | 见下方"刚性极限收敛研究" |
| Champneys & Várkonyi 2016, *The Painlevé paradox in contact mechanics*, IMA J. Appl. Math. 81(3): 538–588 | references/预印本/Champneys_Varkonyi_2016_Painleve_paradox_arXiv1601.03545.pdf | 以"推粉笔"为核心实例：粉笔虚线是含离面跃起与碰撞的 sprag-slip 振荡，可由悖论经 reverse chatter 产生。其图 1(c)：钝角 + 更大法向力 → 点间距更短——与模型预言①"重按短而浓"定性一致，可作外部佐证引用 | ① 同上收敛研究；② 报告/综述中引用其图 1 作为预言①的实验对照 |

**关键观察**：力矩方程中摩擦项的力臂 $\mu h$ 正是产生 Painlevé 悖论的几何耦合（摩擦力矩 $\mu N h$ vs 法向力矩 $NL\sin\theta$ 的竞争）。当前分离是 $\delta$ 光滑归零；刚性极限下应出现"无碰撞跃起"。

**刚性极限收敛研究（建议优先级最高）**：扫描 $k_c=10^3\sim10^7$ N/m，观察 $t^*$、$t_{flip}$、$I_{2nd}$ 是收敛还是发散。
- 若收敛：罚模型自洽，$k_c$ 的物理意义（接触等效刚度）确立，标定路线可行；
- 若在某区间发散：刚性悖论在起作用，分离须改用法向-切向 LCP / 能量恢复碰撞律描述。

这一个扫描同时回应"$k_c$ 标定"与"第二阶段建模"两个问题。

**待补文献**：Sanders J. W. 2018, *Could chalk hopping be caused by reverse chatter?* Proc. ASME DSCC 2018, V003T39A001（doi:10.1115/DSCC2018-8906）——本问题上最对症的一篇，建议经 ASME 或作者索取。

## 3. 接缝三：系统即 2 自由度参数激励振子

| 文献 | 存放位置 | 与模型的对应 | 可直接执行的动作 |
|---|---|---|---|
| Sulollari, van Dalen & Cabboi 2025, *Parametric excitation and friction modulation for a forced 2-DOF system*, Nonlinear Dyn. 113: 12793–12816 | references/期刊论文/Sulollari_etal_2025_Parametric_excitation_friction_modulation_2DOF_NonlinDyn.pdf | 移动带 + Hertz-Damp 法向接触 + 切向自由度：法向-切向经摩擦耦合形成参数激励。本模型中 $\delta$（接触弹簧）即法向自由度，$\theta$ 转动即切向自由度，整个 N 下降段就是一次法向力调制半周期；相对速度换向（翻转）对应其 F–v 曲线上的负阻尼段 | 用其"运动直接分离法 + 线性化稳定性 + Routh–Hurwitz"工具，在不做全周期积分的情况下预言"滑动失稳进入黏滑"的参数边界（即画实线/画虚线的分界线）；多周期阶段的方法模板 |

## 4. 接缝四：双稳态与实验对照

| 文献 | 存放位置 | 与模型的对应 | 可直接执行的动作 |
|---|---|---|---|
| Tonazzi et al. 2020, *Numerical and experimental analysis of the bi-stable state for frictional continuous system*, Nonlinear Dyn. 102: 1361–1374 | references/期刊论文/Tonazzi_etal_2020_Bistable_state_frictional_continuous_system_NonlinDyn.pdf | 连续摩擦振子在某一速度区间内"连续滑动/黏滑"双稳共存，由初始条件或扰动选择；测得双稳速度区间随载荷、倾角变化 | 实验设计：板速 $v$ 先升后降扫描，若存在双稳区则转变速度出现迟滞回线——把"画线/画虚线边界"从一条线变为有面积的区域，是参数扫描的自然延伸 |
| Barsotti, Bennati & Migliaccio 2024, *Non-linear dynamics of simple elastic systems undergoing friction-ruled stick–slip motions*, CivilEng 5(2): 420–434 | references/期刊论文/Barsotti_etal_2024_Nonlinear_dynamics_stick_slip_CivilEng.pdf | 刚块-双支承模型，黏着相/滑动相各自解析求解 | 分段闭式解可作积分器基准测试；其"解析骨架 + 数值"风格与本课题一致 |

## 5. 上升段与 E1 实验的闭环

| 文献 | 存放位置 | 与模型的对应 | 可直接执行的动作 |
|---|---|---|---|
| Putelat, Dawes & Willis 2011, *On the microphysical foundations of rate-and-state friction*, J. Mech. Phys. Solids 59(5): 1062–1075 | references/期刊论文/Putelat_etal_2011_Microphysical_foundations_rate_state_friction_JMPS.pdf | 状态演化律可写成向稳态松弛的一阶方程，且其函数形式可直接由黏滑起始的实验观测估计 | E1（静摩擦系数及其时间依赖）正是这种测量：用 E1 数据定演化律 → 搭上升段 → 交接参数 $(N_0,\omega_0)$ 由输出量导出而非自由参数，单周期首尾一致的延拓条件获得物理内容 |

**Gu et al. 1984**（*Slip motion and stability of a single degree of freedom elastic system with rate and state dependent friction*, JMPS 32(3): 167–196）是速率-状态稳定性分析的基准文献：临界刚度、稳态滑动失稳边界、周期-1 极限环到倍周期的路径。其滑移事件历时/滑移量对加载刚度的标度，可用作 $t^*(N_0)$、$\ell(N_0)$ 扫描的定量对照基准。**注意：库内该 PDF 已损坏（下载截断），Rice 主页公开副本同样截断，需经图书馆/出版商重取。**

## 6. 建议的优先顺序

1. **$k_c$ 刚性极限扫描 + $\varepsilon$ 灵敏度扫描**——纯数值、见效快；同时决定预言③的稳健性与 $k_c$ 标定路线，并把 Painlevé 问题转化为可计算问题。
2. **$(N_0,\omega_0)$ 自洽化**——E1 数据 + Elmer/Putelat 演化律搭上升段，使交接参数成为模型输出。
3. **多周期与稳定性工具**——Sulollari 的线性稳定性分析预测画线/虚线边界；拼自由段 + 恢复系数成冲击映射，用 González-Carbajal 的分岔语言分析自持条件。
4. **实验对照**——Tonazzi 式迟滞扫描；以 Champneys–Várkonyi 图 1 佐证预言①。

## 7. 资料管理备注

- Gu et al. 1984 PDF 损坏待重取（见第 5 节）。
- Sanders 2018（DSCC）建议索取全文（见第 2 节）。
- 本笔记对应模型源码：粉笔_N下降段_完整模型.py、划痕相对强度_时空依赖.py；模型报告：粉笔_N下降段_建模与数值模拟.pdf。
