# -*- coding: utf-8 -*-
"""
粉笔书写 N 下降段（单周期后半段）完整数值模拟
模型：自由轴刚体棒 + 板材接触柔度 + 库仑摩擦（含摩擦方向翻转）
区间：[t0, t*]  接触保持，N: N0 -> 0（分离即停，不接自由段）
"""
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# ---------------- 物理参数（SI） ----------------
L   = 0.08            # 粉笔长 m
h   = 0.075           # 轴 O 到板面距离 m（h < L）
r   = 3.5e-3          # 粉笔半径 m
rho = 2300.0          # 密度 kg/m^3（石膏）
g   = 9.8
A   = np.pi * r**2
m   = rho * A * L
I_O = m * L**2 / 3.0
mu  = 0.35            # 动摩擦系数
kc  = 2.0e4           # 等效接触刚度 N/m（标定参数）
v   = 0.5             # 板速 m/s（反向移动，仅笔痕强度与摩擦翻转用到）
eps = 0.02            # 摩擦方向正则化速度尺度 m/s

th_star = np.arccos(h / L)   # 分离角 N=0
print(f"棒质量 m = {m*1e3:.2f} g,  I_O = {I_O:.3e} kg·m²")
print(f"分离角 θ* = {np.degrees(th_star):.2f}°,  μh = {mu*h:.4f} m·(无量纲臂)")

# ---------------- 模型 ----------------
def delta(th):
    return L * np.cos(th) - h          # 压入量（>0 表示接触）

def rhs(t, y):
    th, w = y
    d = delta(th)
    if d <= 0.0:
        return [w, -(m*g*L/2)*np.sin(th)/I_O]  # 已分离（事件应已截停，仅为保险）
    N = kc * d
    v_rel = v - L*np.cos(th)*w                  # 尖端相对板面的滑动速度
    tau = (m*g*L/2)*np.sin(th) - N*(L*np.sin(th) + mu*h*np.tanh(v_rel/eps))
    return [w, -tau/I_O]   # 注意：θ 增大=顺时针，θ̈ = -τ_z/I_O

def ev_sep(t, y):            # 事件：压入量降到 0（θ 增大 ⇒ δ 减小）
    return delta(y[0])
ev_sep.terminal = True
ev_sep.direction = -1

# ---------------- 单段计算 ----------------
def run_segment(N0, omega0=0.0, t_max=0.05):
    th0 = np.arccos(np.clip((h + N0/kc)/L, -1.0, 1.0))
    sol = solve_ivp(rhs, (0, t_max), [th0, omega0],
                    method='DOP853', rtol=1e-10, atol=1e-12,
                    events=ev_sep, dense_output=True, max_step=5e-5)
    separated = len(sol.t_events[0]) > 0
    t_end = sol.t_events[0][0] if separated else t_max
    tt = np.linspace(0, t_end, 2000)
    th, w = sol.sol(tt)
    N  = kc * np.clip(delta(th), 0.0, None)
    v_rel = v - L*np.cos(th)*w                    # 相对滑动速度
    I2 = mu * np.trapezoid(N*np.abs(v_rel), tt)   # 本段笔痕强度（摩擦功，一致定义）
    ell = np.trapezoid(np.abs(v_rel), tt)         # 痕迹长度
    xdot = L*np.cos(th)*np.abs(w)                 # 尖端水平速度幅值
    return dict(N0=N0, omega0=omega0, th0=th0, t_star=t_end,
                separated=separated, tt=tt, th=th, w=w, N=N, I2=I2, ell=ell,
                xdot_max=xdot.max(), flip=(xdot>v).any())

# ---------------- 名义工况 ----------------
case = run_segment(N0=4.0, omega0=0.0)
th_flip = case['th'][np.argmax(L*np.cos(case['th'])*np.abs(case['w']) > v)]
t_flip  = case['tt'][np.argmax(L*np.cos(case['th'])*np.abs(case['w']) > v)]
print(f"\n名义工况 N0 = 4 N, ω0 = 0：")
print(f"  θ0 = {np.degrees(case['th0']):.3f}°")
print(f"  分离时间 t* = {case['t_star']*1e3:.3f} ms")
print(f"  摩擦翻转时刻 t_flip = {t_flip*1e3:.3f} ms,  ω(t*) = {case['w'][-1]:.2f} rad/s")
print(f"  尖端最大水平速度 = {case['xdot_max']:.2f} m/s（板速 {v} m/s ⇒ 翻转={'是' if case['flip'] else '否'}）")
print(f"  本段笔痕强度 I_2nd = {case['I2']:.3e} J")
print(f"  板面位移 v·t* = {v*case['t_star']*1e3:.3f} mm,  痕迹长 ℓ = ∫|v_rel|dt = {case['ell']*1e3:.3f} mm")
print(f"  平均线密度 μ〈N〉 = I_2nd/ℓ = {case['I2']/case['ell']:.3f} N")

# ---------------- 收敛性自检 ----------------
case_h = run_segment(N0=4.0, omega0=0.0)
print(f"  收敛性：重复积分 t* 偏差 = {abs(case_h['t_star']-case['t_star']):.1e} s")

# ---------------- 参数扫描 ----------------
N0_list = [1.0, 2.0, 4.0, 8.0]
om_list = [0.0, 3.0, 6.0]
rows = []
for N0 in N0_list:
    for om in om_list:
        r_ = run_segment(N0, om)
        rows.append((N0, om, r_['t_star']*1e3, r_['I2'], r_['w'][-1]))
print("\n扫描结果（t* 单位 ms，I_2nd 单位 J，ω* 单位 rad/s）：")
print(f"{'N0/N':>6} {'ω0':>5} {'t*':>9} {'I_2nd':>11} {'ω*':>8}")
for N0, om, ts, I2, ws in rows:
    print(f"{N0:>6.1f} {om:>5.1f} {ts:>9.3f} {I2:>11.3e} {ws:>8.2f}")

# ---------------- 作图 ----------------
fig, ax = plt.subplots(2, 2, figsize=(12, 8))

ax[0,0].plot(case['tt']*1e3, np.degrees(case['th']), 'b-', lw=1.8)
ax[0,0].axvline(t_flip*1e3, color='gray', ls='--', lw=1)
ax[0,0].axvline(case['t_star']*1e3, color='r', ls=':', lw=1.4)
ax[0,0].set_xlabel('t / ms'); ax[0,0].set_ylabel('θ / °')
ax[0,0].set_title('(a) 倾斜角 θ(t)：虚线=摩擦翻转，点线=分离')
ax[0,0].grid(alpha=0.3)

ax[0,1].plot(case['tt']*1e3, case['N'], 'g-', lw=1.8)
ax[0,1].axvline(t_flip*1e3, color='gray', ls='--', lw=1)
ax[0,1].set_xlabel('t / ms'); ax[0,1].set_ylabel('N / N')
ax[0,1].set_title('(b) 压力 N(t)：下降段全程')
ax[0,1].grid(alpha=0.3)

colors = {0.0:'k', 3.0:'b', 6.0:'r'}
for om in om_list:
    sel = [(N0, ts) for (N0, o, ts, I2, ws) in rows if o == om]
    ax[1,0].plot([x[0] for x in sel], [x[1] for x in sel], 'o-', color=colors[om], label=f'ω0={om}')
ax[1,0].set_xlabel('N_0 / N'); ax[1,0].set_ylabel('t* / ms')
ax[1,0].set_title('(c) 分离时间 vs 交接压力'); ax[1,0].legend(); ax[1,0].grid(alpha=0.3)

for om in om_list:
    sel = [(N0, I2) for (N0, o, ts, I2, ws) in rows if o == om]
    ax[1,1].plot([x[0] for x in sel], [x[1] for x in sel], 's-', color=colors[om], label=f'ω0={om}')
ax[1,1].set_xlabel('N_0 / N'); ax[1,1].set_ylabel('I_2nd / J')
ax[1,1].set_title('(d) 本段笔痕强度 vs 交接压力'); ax[1,1].legend(); ax[1,1].grid(alpha=0.3)

plt.tight_layout()
plt.savefig('/mnt/agents/output/粉笔_N下降段_结果.png', dpi=150)
print("\n图已保存：/mnt/agents/output/粉笔_N下降段_结果.png")
