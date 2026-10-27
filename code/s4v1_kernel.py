# -*- coding: utf-8 -*-
"""
粉笔 N 下降段：纯 NumPy 自适应积分内核（S4v1 用，v1 方程）
=================================================
动机
----
S1/S2/S3 用 scipy.integrate.solve_ivp(DOP853, events, dense_output)。
当前执行环境无法安装 scipy（pip 在沙箱内写临时文件被拒，PyPI 大文件下载超时），
因此在本文件中用 NumPy 重新实现同一模型的积分，控制方程、参数、事件定义与
code/仿真S1_kc刚性极限收敛扫描.py 逐行对齐，并用其已发布的 CSV 做回归校验
（见 validate_s1.py）。

与 S1 的唯一实现差异
--------------------
1. 积分器为 Dormand–Prince 5(4)（S1 为 DOP853 8(7)），步长用 PI 控制；
   两者都只影响"数值误差"层面，双容差复算仍可给出误差带。
2. I_2nd 与 ℓ 作为附加状态量随 ODE 一起积分（dI2/dt = μ·N·|v_rel|，
   dℓ/dt = |v_rel|），而不是在 4001 点稠密输出上做梯形积分；
   这样积分量的精度与 θ、ω 同阶。S1 的做法误差约 1e-6 相对量级。
3. 分离事件用步内二分细化，时刻精度 < 1e-12 s。

模型（与主模型一致）
--------------------
δ(θ) = L cosθ − h,          N = k_c·max(δ, 0)
v_rel = v − L cosθ·θ̇
I_O θ̈ = N(L sinθ + μ h tanh(v_rel/ε)) − (m g L/2) sinθ
"""
from __future__ import annotations

import numpy as np

# ---------------- 物理参数（与 code/粉笔_N下降段_完整模型.py 一致） ----------------
L = 0.08                 # 粉笔长 m
H = 0.075                # 轴 O 到板面距离 m（默认值；S4 中按 h/L 扫描）
RADIUS = 3.5e-3          # 粉笔半径 m
RHO = 2300.0             # 密度 kg/m^3
GRAV = 9.8
M = RHO * np.pi * RADIUS ** 2 * L
I_O = M * L ** 2 / 3.0
V_BOARD = 0.5            # 板速 m/s
EPS = 0.02               # tanh 正则化尺度 m/s

TH_STAR_DEFAULT = np.arccos(H / L)

# ---------------- Dormand–Prince 5(4) 系数 ----------------
_C = (0.0, 1 / 5, 3 / 10, 4 / 5, 8 / 9, 1.0)
_A = (
    (),
    (1 / 5,),
    (3 / 40, 9 / 40),
    (44 / 45, -56 / 15, 32 / 9),
    (19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729),
    (9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656),
)
_B5 = (35 / 384, 0.0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84)
_B4 = (5179 / 57600, 0.0, 7571 / 16695, 393 / 640, -92097 / 339200,
       187 / 2100, 1 / 40)


def delta(th, h=H):
    """压入量 δ = L cosθ − h（>0 表示接触）。"""
    return L * np.cos(th) - h


def contact_omega(kc, th, h=H):
    """接触弹簧-转动振子特征角频率（用于上限步长）。"""
    return float(np.sqrt(max(kc * (L * np.sin(th)) ** 2 / I_O, 1.0)))


def make_rhs(kc, mu, h=H):
    """返回 f(t, y)，y = [θ, ω, I2, ℓ]。"""
    def rhs(t, y):
        th, w = y[0], y[1]
        d = L * np.cos(th) - h
        if d <= 0.0:
            # 已分离（事件应已截停，仅为保险）：无接触，仅重力
            acc = -(M * GRAV * L / 2.0) * np.sin(th) / I_O
            return np.array([w, acc, 0.0, 0.0])
        Nn = kc * d
        v_rel = V_BOARD - L * np.cos(th) * w
        tau = (M * GRAV * L / 2.0) * np.sin(th) \
            - Nn * (L * np.sin(th) + mu * h * np.tanh(v_rel / EPS))
        return np.array([w, -tau / I_O,
                         mu * Nn * abs(v_rel),      # dI2/dt
                         abs(v_rel)])               # dℓ/dt
    return rhs


def _rk_step(f, t, y, hstep):
    """单步 DOPRI5，返回 (y5, y4)。"""
    k = []
    for i in range(6):
        yi = y.copy()
        for j in range(i):
            yi = yi + hstep * _A[i][j] * k[j]
        k.append(f(t + _C[i] * hstep, yi))
    y5 = y + hstep * sum(b * ki for b, ki in zip(_B5, k))
    k6 = f(t + hstep, y5)          # FSAL
    y4 = y + hstep * sum(b * ki for b, ki in zip(_B4, list(k) + [k6]))
    return y5, y4


def integrate_segment(kc, N0, omega0=0.0, mu=0.35, h=H,
                      rtol=1e-10, atol=1e-12, t_max=0.05, record=False):
    """积分单段 [0, t*]，返回结果字典。

    事件：δ 由正变负（分离）。步内二分细化到 1e-13 s。
    record=True 时额外返回时间序列（供物理判读，不参与扫描）。
    """
    d0 = N0 / kc
    cos0 = (h + d0) / L
    out = dict(kc=kc, N0=N0, omega0=omega0, mu=mu, h=h, valid=False,
               separated=False, t_star=np.nan, t_flip=np.nan, flip=False,
               I2=np.nan, ell=np.nan, omega_star=np.nan, xdot_max=np.nan,
               theta0=np.nan, arm_ratio0=np.nan, steps=0, note="")
    if cos0 >= 1.0:
        out["note"] = "几何不可达"
        return out

    th0 = float(np.arccos(np.clip(cos0, -1.0, 1.0)))
    out["theta0"] = th0
    out["arm_ratio0"] = mu * h / (L * np.sin(th0))
    out["valid"] = True

    f = make_rhs(kc, mu, h)
    y = np.array([th0, omega0, 0.0, 0.0], dtype=float)
    t = 0.0
    max_step = min(5e-5, 0.05 * 2 * np.pi / contact_omega(kc, th0, h))
    hstep = min(max_step, 1e-6)
    nsteps = 0
    last_vrel_sign = np.sign(V_BOARD - L * np.cos(th0) * omega0)
    t_flip = np.nan
    xdot_max = L * abs(np.cos(th0)) * abs(omega0)
    rec_t, rec_th, rec_w = [], [], []
    if record:
        rec_t.append(0.0); rec_th.append(th0); rec_w.append(omega0)

    while t < t_max:
        hstep = min(hstep, max_step, t_max - t)
        if hstep <= 1e-15:
            break
        y5, y4 = _rk_step(f, t, y, hstep)
        nsteps += 1

        # 误差控制只针对 θ、ω（积分量不参与步长选择）
        err = np.abs(y5[:2] - y4[:2])
        scale = atol + rtol * np.maximum(np.abs(y[:2]), np.abs(y5[:2]))
        err_norm = float(np.sqrt(np.mean((err / scale) ** 2)))

        # 事件检测：δ 由正到非正
        d_old = L * np.cos(y[0]) - h
        d_new = L * np.cos(y5[0]) - h
        if d_old > 0.0 and d_new <= 0.0:
            lo, hi = 0.0, hstep
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                ym, _ = _rk_step(f, t, y, mid)
                if L * np.cos(ym[0]) - h > 0.0:
                    lo = mid
                else:
                    hi = mid
                if hi - lo < 1e-14:
                    break
            ym, _ = _rk_step(f, t, y, hi)
            t = t + hi
            y = ym
            if record:
                rec_t.append(t); rec_th.append(float(y[0])); rec_w.append(float(y[1]))
            out.update(separated=True, t_star=t, omega_star=float(y[1]),
                       I2=float(y[2]), ell=float(y[3]))
            xdot_max = max(xdot_max, float(L * np.cos(min(y[0], th0)) * abs(y[1])))
            break

        # 摩擦翻转检测（v_rel 过零，步内线性插值）
        vr_old = V_BOARD - L * np.cos(y[0]) * y[1]
        vr_new = V_BOARD - L * np.cos(y5[0]) * y5[1]
        if np.sign(vr_new) != last_vrel_sign and np.isfinite(vr_new):
            if last_vrel_sign > 0 and vr_new < 0:
                frac = vr_old / (vr_old - vr_new) if vr_old != vr_new else 0.0
                t_flip = t + frac * hstep
            last_vrel_sign = np.sign(vr_new)

        t = t + hstep
        y = y5
        if record:
            rec_t.append(t); rec_th.append(float(y[0])); rec_w.append(float(y[1]))
        xdot_max = max(xdot_max, float(L * abs(np.cos(y[0])) * abs(y[1])))

        # PI 步长控制
        if err_norm == 0.0:
            fac = 5.0
        else:
            fac = 0.9 * err_norm ** (-0.2)
        hstep = hstep * min(5.0, max(0.2, fac))

    out["steps"] = nsteps
    out["xdot_max"] = xdot_max
    if not out["separated"]:
        out["note"] = "t_max 内未分离"
        out["t_star"] = t
        out["omega_star"] = float(y[1])
        out["I2"] = float(y[2])
        out["ell"] = float(y[3])
    out["flip"] = bool(np.isfinite(t_flip))
    out["t_flip"] = float(t_flip) if np.isfinite(t_flip) else np.nan
    if record:
        out["traj"] = dict(t=np.array(rec_t), th=np.array(rec_th), w=np.array(rec_w),
                           N=np.array([kc * max(L * np.cos(x) - h, 0.0) for x in rec_th]),
                           vrel=np.array([V_BOARD - L * np.cos(a) * b
                                          for a, b in zip(rec_th, rec_w)]))
    return out


def run_dual(kc, N0, omega0=0.0, mu=0.35, h=H,
             rtol_pair=(1e-9, 1e-12)):
    """双容差复算，返回 (结果, 数值误差 t*)。"""
    a = integrate_segment(kc, N0, omega0, mu, h, rtol=rtol_pair[0])
    b = integrate_segment(kc, N0, omega0, mu, h, rtol=rtol_pair[1])
    err_t = abs(a["t_star"] - b["t_star"]) * 1e3 if np.isfinite(a["t_star"]) else np.nan
    return a, err_t
