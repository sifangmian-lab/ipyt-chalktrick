# -*- coding: utf-8 -*-
"""构建《粉笔书写单周期压力下降段的建模与数值模拟》PDF（ReportLab 路线）"""
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import cm
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Image,
                                Table, TableStyle, PageBreak)
from reportlab.platypus import Preformatted
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors

# ---------------- 字体注册（宋体正文 / 黑体标题 / 等宽代码，子集化 TTF） ----------------
FD = '/mnt/agents/output/fonts/'
pdfmetrics.registerFont(TTFont('Song',  FD + 'SongSC.ttf'))
pdfmetrics.registerFont(TTFont('SongB', FD + 'SongSC-Bold.ttf'))
pdfmetrics.registerFont(TTFont('Hei',   FD + 'HeiSC-Bold.ttf'))
pdfmetrics.registerFont(TTFont('HeiR',  FD + 'HeiSC-Bold.ttf'))
pdfmetrics.registerFont(TTFont('Mono',  FD + 'MonoSC.ttf'))

# ---------------- 样式 ----------------
S = {}
S['title']  = ParagraphStyle('title', fontName='Hei', fontSize=17, leading=26,
                             alignment=TA_CENTER, spaceAfter=6)
S['subtitle'] = ParagraphStyle('subtitle', fontName='HeiR', fontSize=12, leading=18,
                               alignment=TA_CENTER, spaceAfter=4)
S['info']   = ParagraphStyle('info', fontName='Song', fontSize=9.5, leading=14,
                             alignment=TA_CENTER, textColor=colors.black, spaceAfter=14)
S['h1']     = ParagraphStyle('h1', fontName='Hei', fontSize=13.5, leading=20,
                             spaceBefore=14, spaceAfter=8)
S['h2']     = ParagraphStyle('h2', fontName='Hei', fontSize=11.5, leading=17,
                             spaceBefore=10, spaceAfter=6)
S['body']   = ParagraphStyle('body', fontName='Song', fontSize=10.5, leading=17,
                             firstLineIndent=21, alignment=TA_LEFT, spaceAfter=4)
S['eq']     = ParagraphStyle('eq', fontName='Song', fontSize=10.5, leading=14,
                             alignment=TA_CENTER, spaceBefore=6, spaceAfter=8)
S['caption']= ParagraphStyle('caption', fontName='Song', fontSize=9.5, leading=14,
                             alignment=TA_CENTER, spaceBefore=4, spaceAfter=10)
S['abstract'] = ParagraphStyle('abstract', fontName='Song', fontSize=10.5, leading=17,
                               firstLineIndent=21, spaceAfter=4)
S['code']   = ParagraphStyle('code', fontName='Mono', fontSize=7.2, leading=9.6,
                             leftIndent=6, spaceBefore=4, spaceAfter=4)
S['tocell'] = ParagraphStyle('tocell', fontName='Song', fontSize=9.5, leading=13,
                             alignment=TA_CENTER)

EQ = '/mnt/agents/output/eq/'

def eq(name):
    from PIL import Image as PILImage
    w, h = PILImage.open(EQ + name + '.png').size
    disp_w = min(w / 300 * 72 * 1.06, 15 * cm)
    return Image(EQ + name + '.png', width=disp_w, height=disp_w * h / w)

def P(text, style='body'):
    return Paragraph(text, S[style])

# ---------------- 页码 ----------------
def on_page(canv, doc):
    canv.saveState()
    canv.setFont('Song', 9)
    canv.drawCentredString(A4[0] / 2, 1.2 * cm, f'— {doc.page} —')
    canv.restoreState()

# ---------------- 文档内容 ----------------
story = []

story.append(P('粉笔书写单周期压力下降段的建模与数值模拟', 'title'))
story.append(P('——自由轴刚体棒模型（第一阶段）', 'subtitle'))
story.append(P('研究日期：2026 年 10 月 20 日', 'info'))

story.append(P('摘要', 'h1'))
story.append(P('针对粉笔书写过程中单个周期内压力从交接值 N<sub>0</sub> 下降到零（笔尖与板面分离）的后半段，'
               '本文建立了"自由轴刚体棒 + 板面接触柔度 + 库仑摩擦（含摩擦方向翻转）"的单自由度模型。'
               '模型给出压力 N(θ) 的精确几何本构与分离角的解析表达式，证明下降段在无外力驱动下即可自启动；'
               '解析分析揭示出"启动段—滑行段"的两阶段图像与摩擦方向翻转的临界角速度。'
               '数值积分（DOP853，事件检测）给出名义工况（N<sub>0</sub> = 4 N）下分离时间 t* = 1.29 ms、'
               '本段笔痕强度 I<sub>2</sub>nd = 1.8×10<super>−4</super> J；对交接压力 N<sub>0</sub> ∈ [1, 8] N 与交接角速度 ω<sub>0</sub> ∈ {0, 3, 6} rad/s '
               '的扫描给出分离时间、笔痕强度的依赖规律及三条可供实验检验的模型预言；'
               '进一步以 Archard 型磨耗律把压力史映射为划痕相对强度的时间分布与空间分布，'
               '揭示末端回划折叠导致的强度尖峰与"重按短痕"效应。', 'abstract'))
story.append(P('<b>关键词</b>：粉笔书写；接触力学；干摩擦；事件驱动数值模拟；笔痕强度', 'body'))
story.append(Spacer(1, 8))

story.append(P('1　问题描述与研究区间', 'h1'))
story.append(P('书写过程可分解为若干压力周期。总计划是先研究单个周期内笔痕强度随压力的变化规律，'
               '再延拓到多周期（只需保证单周期首尾状态一致）。本文只处理单个周期内的压力下降段，'
               '即图 1 中从压力最大值（交接值 N<sub>0</sub>）到笔尖分离（N = 0）的区间 [t<sub>0</sub>, t*]。'
               '区间内恒保持接触，分离时刻即为区间终点，分离后的自由飞行段不在本文范围内。'))
story.append(P('实验装置方面：手握端简化为绕定点 O 在竖直平面内自由转动的轴（无力矩），'
               '板面在机构驱动下沿书写反方向匀速移动，粉笔尖端在板面上拖行并留下粉痕。'))

story.append(P('2　模型建立', 'h1'))
story.append(P('2.1　基本假设', 'h2'))
story.append(P('（1）粉笔视为匀质刚体细棒，长 L，质量 m；弯曲弹性不进入静力学，留作第二阶段的修正；'
               '（2）手握端 O 固定于板面上方 h 处，h &lt; L，轴无力矩；'
               '（3）板面平直，接触柔度以罚弹簧 k_c 集中于板面（k_c 为标定参数）；'
               '（4）库仑动摩擦，摩擦方向由尖端相对板面的滑动速度自动判定，允许翻转；'
               '（5）忽略轴向变形；重力保留但在卸载段为小修正。'))
story.append(P('2.2　几何关系与符号约定', 'h2'))
story.append(P('取 θ 为棒与板面法线的夹角，θ &gt; 0 表示笔尖偏向板面运动来的一侧（trailing 位形）。'
               '板以速度 v 沿 −x 方向运动，粉笔相对板面沿 +x 书写。θ 增大的方向为顺时针。'
               '压入量 δ(θ) = L cosθ − h，接触压力由板面柔度给出。分离角 θ* 由 N = 0 定义。'))
story.append(eq('eq_constit'))
story.append(P('式 (1)', 'caption'))
story.append(eq('eq_tstar'))
story.append(P('式 (2)', 'caption'))
story.append(P('2.3　控制方程与力矩分析', 'h2'))
story.append(P('对 O 点取矩。注意 θ 增大为顺时针方向，故角加速度与力矩 z 分量反号：'))
story.append(eq('eq_torque'))
story.append(P('式 (3)', 'caption'))
story.append(P('合并后得到整个下降段的控制方程（单自由度）：'))
story.append(eq('eq_model'))
story.append(P('式 (4)', 'caption'))
story.append(P('三项力矩的物理作用：法向力矩 NL sinθ 恒指向卸载方向（板上顶力把棒"顶倒"）；'
               '摩擦力矩经 tanh 正则化，v_rel &gt; 0 时卸载、v_rel &lt; 0 时恢复；'
               '重力矩为单摆恢复项，在卸载全程占比不超过 1%，仅在接近分离时显现。'))
story.append(P('2.4　初值交接与分离事件', 'h2'))
story.append(P('下降段初值由前半段（压力上升段）末态交接给出：N<sub>0</sub> 唯一确定 θ<sub>0</sub>，角速度 ω<sub>0</sub> 为交接参数。'))
story.append(eq('eq_ic'))
story.append(P('式 (5)', 'caption'))
story.append(P('数值积分以 δ(θ) = 0（方向为 θ 增大）为终端事件，事件触发即分离，区间终止。'))

story.append(P('3　解析骨架', 'h1'))
story.append(P('3.1　自启动判据', 'h2'))
story.append(P('ω<sub>0</sub> = 0 时下降段也能自行启动的条件为：'))
story.append(eq('eq_start'))
story.append(P('式 (6)', 'caption'))
story.append(P('右端量级约为 0.02 N，任何真实交接压力（≳ 0.1 N）均满足——'
               '<b>下降段无需外力驱动，接触几何本身即完成卸载</b>。'))
story.append(P('3.2　两阶段图像与摩擦翻转', 'h2'))
story.append(P('数值与解析分析一致表明卸载分两阶段：'
               '<b>启动段</b>（0 &lt; t &lt; t_flip）法向力矩与摩擦力矩同向卸载，棒强加速；'
               '当尖端水平速度赶上板速时相对滑动换向，摩擦翻转进入<b>滑行段</b>，'
               'μh ≈ L sinθ 近抵消，棒以近匀速滑行至分离。翻转临界角速度为'))
story.append(eq('eq_flip'))
story.append(P('式 (7)', 'caption'))
story.append(P('接近分离时 N 随角度亏缺线性归零：'))
story.append(eq('eq_linear'))
story.append(P('式 (8)', 'caption'))
story.append(P('3.3　能量关系与积分形式的分离时间', 'h2'))
story.append(P('下降段无耗散机制（摩擦翻转后近似不做净功），可用能量恒等式把 t* 写成积分形式作独立验证：'))
story.append(eq('eq_energy'))
story.append(P('式 (9)', 'caption'))

story.append(P('4　数值方法与结果', 'h1'))
story.append(P('4.1　参数与数值方法', 'h2'))
story.append(P('采用 DOP853 变步长积分（rtol = 10<super>−10</super>），以 δ(θ) 降零为终端事件并二分回插定位 t*；'
               '笔痕强度沿轨迹数值积分获得。参数取值见表 1。'))
tab1 = Table([
    [P('参数', 'tocell'), P('符号', 'tocell'), P('取值', 'tocell'), P('说明', 'tocell')],
    [P('棒长', 'tocell'), P('L', 'tocell'), P('0.08 m', 'tocell'), P('粉笔', 'tocell')],
    [P('轴—板距离', 'tocell'), P('h', 'tocell'), P('0.075 m', 'tocell'), P('h &lt; L', 'tocell')],
    [P('棒质量', 'tocell'), P('m', 'tocell'), P('7.08 g', 'tocell'), P('ρ = 2300 kg/m³，r = 3.5 mm', 'tocell')],
    [P('转动惯量', 'tocell'), P('I_O', 'tocell'), P('1.51×10<super>−5</super> kg·m²', 'tocell'), P('mL²/3', 'tocell')],
    [P('接触刚度', 'tocell'), P('k_c', 'tocell'), P('2×10<super>4</super> N/m', 'tocell'), P('标定参数', 'tocell')],
    [P('摩擦系数', 'tocell'), P('μ', 'tocell'), P('0.35', 'tocell'), P('动摩擦', 'tocell')],
    [P('板速', 'tocell'), P('v', 'tocell'), P('0.5 m/s', 'tocell'), P('反向移动', 'tocell')],
    [P('分离角', 'tocell'), P('θ*', 'tocell'), P('20.36°', 'tocell'), P('arccos(h/L)', 'tocell')],
], colWidths=[3.0*cm, 1.6*cm, 4.2*cm, 5.4*cm])
tab1.setStyle(TableStyle([
    ('FONTNAME', (0,0), (-1,-1), 'Song'),
    ('FONTSIZE', (0,0), (-1,-1), 9.5),
    ('GRID', (0,0), (-1,-1), 0.5, colors.black),
    ('BACKGROUND', (0,0), (-1,0), colors.Color(0.92,0.92,0.92)),
    ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
    ('ALIGN', (0,0), (-1,-1), 'CENTER'),
    ('TOPPADDING', (0,0), (-1,-1), 3), ('BOTTOMPADDING', (0,0), (-1,-1), 3),
]))
story.append(tab1)
story.append(P('表 1　模型参数', 'caption'))

story.append(P('4.2　名义工况结果', 'h2'))
story.append(P('N<sub>0</sub> = 4 N、ω<sub>0</sub> = 0 的名义工况：t* = 1.292 ms；摩擦翻转发生于 0.526 ms；'
               '分离角速度 ω(t*) = 7.12 rad/s；尖端最大水平速度 0.53 m/s ≳ v，翻转确实发生。'
               '由于尖端滑动速度 v_rel 全程显著偏离板速 v（翻转前减小、翻转后反向增大），'
               '笔痕可观测量必须用相对滑动速度自洽定义：'))
story.append(eq('eq_obs'))
story.append(P('式 (10)', 'caption'))
story.append(P('由此得本段笔痕强度 I<sub>2</sub>nd = 1.75×10<super>−4</super> J，痕迹总长 ℓ = 0.145 mm，'
               '平均线密度 μ〈N〉 = 1.21 N。注意板面位移 v t* = 0.646 mm 明显大于痕迹总长——'
               '摩擦翻转后尖端几乎"钉"在板面原处滑行，二者不再相等。'
               'N(t) 全程单调下凹（图 1b），θ(t) 单调上升至 θ*（图 1a）。'
               '收敛性自检：两次独立积分 t* 偏差为 0。'))
img = Image('/mnt/agents/output/粉笔_N下降段_结果.png', width=15.5*cm, height=15.5*cm*1050/1575)
story.append(img)
story.append(P('图 1　数值结果。(a) θ(t)，虚线为摩擦翻转时刻，点线为分离时刻；'
               '(b) 压力 N(t) 下降段全程；(c) 分离时间 t* 对交接压力 N<sub>0</sub> 的依赖；'
               '(d) 本段笔痕强度 I<sub>2</sub>nd 对 N<sub>0</sub> 的依赖。', 'caption'))

story.append(P('4.3　参数扫描与模型预言', 'h2'))
story.append(P('对 N<sub>0</sub> ∈ {1, 2, 4, 8} N 与 ω<sub>0</sub> ∈ {0, 3, 6} rad/s 扫描，结果见表 2。'))
rows = [[P('N<sub>0</sub> / N', 'tocell'), P('ω<sub>0</sub> = 0：t* / ms', 'tocell'), P('ω<sub>0</sub> = 0：I<sub>2</sub>nd / J', 'tocell'),
         P('ω<sub>0</sub> = 0：ℓ / mm', 'tocell'), P('ω<sub>0</sub> = 3：t* / ms', 'tocell'), P('ω<sub>0</sub> = 6：t* / ms', 'tocell')]]
data = [('1.0','1.128','1.02×10<super>−4</super>','0.429','0.500','0.285'),
        ('2.0','1.124','1.56×10<super>−4</super>','0.291','0.740','0.535'),
        ('4.0','1.292','1.75×10<super>−4</super>','0.145','1.115','1.026'),
        ('8.0','2.121','2.87×10<super>−4</super>','0.161','2.034','1.988')]
for r in data:
    rows.append([P(x, 'tocell') for x in r])
tab2 = Table(rows, colWidths=[1.8*cm, 2.9*cm, 3.3*cm, 2.6*cm, 2.9*cm, 2.9*cm])
tab2.setStyle(TableStyle([
    ('FONTNAME', (0,0), (-1,-1), 'Song'),
    ('GRID', (0,0), (-1,-1), 0.5, colors.black),
    ('BACKGROUND', (0,0), (-1,0), colors.Color(0.92,0.92,0.92)),
    ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
    ('ALIGN', (0,0), (-1,-1), 'CENTER'),
    ('TOPPADDING', (0,0), (-1,-1), 3), ('BOTTOMPADDING', (0,0), (-1,-1), 3),
]))
story.append(tab2)
story.append(P('表 2　参数扫描结果（t* 为分离时间，I<sub>2</sub>nd 与 ℓ 采用自洽定义，式 (10)）', 'caption'))
story.append(P('由扫描可提炼三条可供实验检验的模型预言：'
               '（1）t* 与 I<sub>2</sub>nd 总体随 N<sub>0</sub> 增大，且 ω<sub>0</sub> = 0 时 I<sub>2</sub>nd 对 N<sub>0</sub> 呈亚线性幂律（指数约 0.5），'
               '同时痕迹变短、线密度增大——"重按短而浓，轻按长而淡"；'
               '（2）交接角速度 ω<sub>0</sub> 越大，t* 越短、本段笔痕越淡——前半段释放越猛，后半段沉积越少，'
               '与"弹跳式书写痕迹变淡"的直觉一致；'
               '（3）N(t) 呈"先下凹后近线性"的两阶段形态，摩擦翻转时刻可由力传感器实测曲线检验。'))

story.append(P('4.4　划痕相对强度的时空依赖', 'h2'))
story.append(P('把压力史与运动学史转化为可观测的划痕强度，需引入沉积律。'
               '采用 Archard 型磨耗律：单位时间沉积粉量正比于摩擦功耗，即沉积率为'))
story.append(eq('eq_deposit'))
story.append(P('式 (11)', 'caption'))
story.append(P('其中 α 为标定常数（含磨损系数与密度因子，求相对强度时约去），该沉积律为 Archard 型磨耗律。'
               '时间依赖直接由 q(t) 给出（图 2a）；空间依赖需引入板面坐标映射——'
               '以 t = 0 时尖端触点为原点，ξ 为尖端相对板面的累计位移：'))
story.append(eq('eq_ximap'))
story.append(P('式 (12)', 'caption'))
story.append(P('当 ξ(t) 单调时空间强度即 q(t)/|v_rel(t)|；'
               '但摩擦翻转后尖端相对板面回划，ξ(t) 非单调，同一板面位置被多次覆盖，'
               '强度须在回划折叠处按覆盖次数求和：'))
story.append(eq('eq_foldsum'))
story.append(P('式 (13)', 'caption'))
story.append(P('数值实现上按 ξ 分箱直方图累加 q·|dξ| 并作轻微高斯平滑（对应有限光斑测量）。'
               '名义工况下 ξ(t) 先前进 0.123 mm 再回划 0.022 mm（图 2b），'
               '折叠叠加在 ξ ≈ 0.10～0.12 mm 处形成强度尖峰（图 2c 黑线右端）。'
               '扫描结果表明：ω<sub>0</sub> 较大时痕迹整体向负 ξ 延伸（尖端前冲回拖），沉积近乎为零（图 2c）；'
               'ω<sub>0</sub> = 0 时随 N<sub>0</sub> 增大，痕迹显著变短而峰值增高（图 2d），'
               '与表 2 的线密度规律一致——这为"由划痕空间轮廓反演交接状态"提供了直接的模拟工具。'))
from PIL import Image as PILImage
_w2, _h2 = PILImage.open('/mnt/agents/output/划痕相对强度_时空依赖.png').size
img2 = Image('/mnt/agents/output/划痕相对强度_时空依赖.png',
             width=15.5*cm, height=15.5*cm*_h2/_w2)
story.append(img2)
story.append(P('图 2　划痕相对强度的时空依赖。(a) 时间域：归一化压力 N(t) 与沉积率 q(t)；'
               '(b) 空间映射 ξ(t)，末端回划即折叠源；'
               '(c) 空间强度分布 I(ξ)：固定 N<sub>0</sub> = 4 N，变交接角速度 ω<sub>0</sub>；'
               '(d) 空间强度分布：ω<sub>0</sub> = 0，变交接压力 N<sub>0</sub>。', 'caption'))

story.append(P('5　模型边界与后续工作', 'h1'))
story.append(P('（1）当 ω<sub>0</sub> 大到在摩擦翻转前即分离时，系统进入"直接弹起"区，两阶段图像退化；'
               '（2）若实测 N(t) 叠加百赫兹级高频振荡，需引入粉笔一阶弯曲模态（振动方程），进入第二阶段模型；'
               '（3）t* 之后接自由段与恢复系数碰撞律，即可组装多周期冲击映射，实现单周期向多周期的延拓；'
               '（4）k_c 为标定参数，建议由板下力传感器实测 N(t) 反演。'))
story.append(P('数值计算源码见附录 A（主模型）与附录 B（划痕时空分布），可直接运行复现全部图表。'))

story.append(PageBreak())
story.append(P('附录 A　数值计算源码（主模型）', 'h1'))
story.append(P('文件：粉笔_N下降段_完整模型.py。运行环境：Python 3 + NumPy + SciPy + Matplotlib。', 'body'))
story.append(Spacer(1, 4))
with open('/mnt/agents/output/粉笔_N下降段_完整模型.py', encoding='utf-8') as f:
    code = f.read()
story.append(Preformatted(code, S['code']))

story.append(PageBreak())
story.append(P('附录 B　数值计算源码（划痕相对强度的时空依赖）', 'h1'))
story.append(P('文件：划痕相对强度_时空依赖.py。依赖附录 A 的 run_segment 函数，同一环境运行。', 'body'))
story.append(Spacer(1, 4))
with open('/mnt/agents/output/划痕相对强度_时空依赖.py', encoding='utf-8') as f:
    code2 = f.read()
story.append(Preformatted(code2, S['code']))

doc = SimpleDocTemplate('/mnt/agents/output/粉笔_N下降段_建模与数值模拟.pdf',
                        pagesize=A4,
                        leftMargin=2.4*cm, rightMargin=2.4*cm,
                        topMargin=2.2*cm, bottomMargin=2.2*cm,
                        title='粉笔书写单周期压力下降段的建模与数值模拟',
                        author='')
doc.build(story, onFirstPage=on_page, onLaterPages=on_page)
print('PDF built.')
