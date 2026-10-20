# -*- coding: utf-8 -*-
"""
划痕相对强度的时空依赖模拟（修正版）
沉积率：q(t) = μ·N(t)·|v_rel(t)|，v_rel = v − L·cosθ·θ̇
空间强度：I(ξ) 经 ξ(t) = ∫v_rel dt 折叠求和（含回划）
"""
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# ---------------- 物理参数（与下降段模型一致） ----------------
L, h = 0.08, 0.075
r, rho, g = 3.5e-3, 2300.0, 9.8
A = np.pi * r**2
m = rho * A * L
I_O = m * L**2 / 3.0
mu, kc, v, eps = 0.35, 2.0e4, 0.5, 0.02
th_star = np.arccos(h / L)

def delta(th): return L * np.cos(th) - h

def rhs(t, y):
    th, w = y
    d = delta(th)
    if d <= 0.0:
        return [w, -(m*g*L/2)*np.sin(th)/I_O]
    N = kc * d
    v_rel = v - L*np.cos(th)*w
    tau = (m*g*L/2)*np.sin(th) - N*(L*np.sin(th) + mu*h*np.tanh(v_rel/eps))
    return [w, -tau/I_O]

def ev_sep(t, y): return delta(y[0])
ev_sep.terminal, ev_sep.direction = True, -1

def run_segment(N0, omega0=0.0, t_max=0.05, npts=20000):
    th0 = np.arccos(np.clip((h + N0/kc)/L, -1.0, 1.0))
    sol = solve_ivp(rhs, (0, t_max), [th0, omega0], method='DOP853',
                    rtol=1e-10, atol=1e-12, events=ev_sep,
                    dense_output=True, max_step=5e-5)
    separated = len(sol.t_events[0]) > 0
    t_end = sol.t_events[0][0] if separated else t_max
    tt = np.linspace(0, t_end, npts)
    th, w = sol.sol(tt)
    N = kc * np.clip(delta(th), 0.0, None)
    v_rel = v - L*np.cos(th)*w
    # 一致的沉积量：μ∫N|v_rel|dt
    I2 = mu * np.trapezoid(N*np.abs(v_rel), tt)
    # 痕迹长度：∫|v_rel|dt
    ell = np.trapezoid(np.abs(v_rel), tt)
    return tt, th, w, N, v_rel, separated, I2, ell

# ---------------- 空间强度分布 ----------------
def trace_profiles(tt, th, w, N, v_rel, nbins=220):
    q = mu * N * np.abs(v_rel)                 # 单位时间沉积率
    xi = np.concatenate([[0.0], np.cumsum(0.5*(v_rel[:-1]+v_rel[1:])*np.diff(tt))])
    dxi = np.abs(np.diff(xi))
    qm = 0.5*(q[:-1] + q[1:])
    xim = 0.5*(xi[:-1] + xi[1:])
    edges = np.linspace(xi.min(), xi.max(), nbins + 1)
    I, _ = np.histogram(xim, bins=edges, weights=qm*dxi)
    I = I / np.gradient(edges)[:-1]
    # 轻度平滑消除分箱伪影（高斯窗宽 ~2% 痕迹长）
    from scipy.ndimage import gaussian_filter1d
    I = gaussian_filter1d(I, sigma=2.0)
    return q, xi, edges[:-1], I

# ---------------- 图 ----------------
fig, ax = plt.subplots(2, 2, figsize=(12, 8))

# (a) 时间域
res = run_segment(4.0, 0.0)
tt, th, w, N, v_rel = res[:5]
q_n, xi_n, xc_n, I_n = trace_profiles(*res[:5])
ax[0,0].plot(tt*1e3, N/np.max(N), 'b-', lw=1.6, label='N(t)/N₀')
ax[0,0].plot(tt*1e3, q_n/np.max(q_n), 'r--', lw=1.6, label='q(t)/max')
ax[0,0].set_xlabel('t / ms'); ax[0,0].set_ylabel('归一化强度')
ax[0,0].set_title('(a) 时间域：压力史与沉积率')
ax[0,0].legend(); ax[0,0].grid(alpha=0.3)

# (b) 空间映射 ξ(t)
ax[0,1].plot(tt*1e3, xi_n*1e3, 'g-', lw=1.6)
ax[0,1].set_xlabel('t / ms'); ax[0,1].set_ylabel('ξ / mm')
ax[0,1].set_title('(b) 空间映射 ξ(t)：末端回划 = 折叠源')
ax[0,1].grid(alpha=0.3)

# (c) 固定 N₀=4，变 ω₀
for om, c in zip([0.0, 3.0, 6.0], ['k', 'b', 'r']):
    r_ = run_segment(4.0, om)
    q, xi, xc, I = trace_profiles(*r_[:5])
    ax[1,0].plot(xc*1e3, I/np.max(I), c+'-', lw=1.6, label=f'ω₀={om:.0f}')
ax[1,0].set_xlabel('ξ / mm'); ax[1,0].set_ylabel('I(ξ)/max')
ax[1,0].set_title('(c) 空间强度分布：N₀=4 N，变交接角速度')
ax[1,0].legend(); ax[1,0].grid(alpha=0.3)

# (d) 固定 ω₀=0，变 N₀
for N0, c in zip([1.0, 4.0, 8.0], ['k', 'b', 'r']):
    r_ = run_segment(N0, 0.0)
    q, xi, xc, I = trace_profiles(*r_[:5])
    ax[1,1].plot(xc*1e3, I/np.max(I), c+'-', lw=1.6, label=f'N₀={N0:.0f} N')
ax[1,1].set_xlabel('ξ / mm'); ax[1,1].set_ylabel('I(ξ)/max')
ax[1,1].set_title('(d) 空间强度分布：ω₀=0，变交接压力')
ax[1,1].legend(); ax[1,1].grid(alpha=0.3)

plt.tight_layout()
plt.savefig('/mnt/agents/output/划痕相对强度_时空依赖.png', dpi=150)

# ---------------- 定量输出 ----------------
print("名义工况 N₀=4, ω₀=0：")
print(f"  板面位移 v·t* = {v*tt[-1]*1e3:.3f} mm，痕迹总长 ∫|v_rel|dt = {res[7]*1e3:.3f} mm")
print(f"  ξ_max = {xi_n.max()*1e3:.3f} mm，末端回划量 = {(xi_n.max()-xi_n[-1])*1e3:.3f} mm")
print(f"  I_2nd(旧近似 μv∫Ndt) 对比 I_2nd(一致定义 μ∫N|v_rel|dt) = {res[6]:.3e} J")
print("\n扫描（一致定义）：")
print(f"{'N₀/N':>6} {'ω₀':>5} {'t*/ms':>8} {'I_2nd/J':>11} {'ℓ/mm':>8}")
for N0 in [1.0, 2.0, 4.0, 8.0]:
    for om in [0.0, 3.0, 6.0]:
        r_ = run_segment(N0, om)
        print(f"{N0:>6.1f} {om:>5.1f} {r_[0][-1]*1e3:>8.3f} {r_[6]:>11.3e} {r_[7]*1e3:>8.3f}")
print("\n图已保存：划痕相对强度_时空依赖.png")
