# -*- coding: utf-8 -*-
"""v2 最小模型共享内核：弹性杆 + 手端正交弹簧（3-DOF 集总）。

被 S4（基线）、S5（收敛性）、S6（参数平面）、S7（本构对比）共用。
模型定义见 docs/研究记录/模型v2规格_2026-10-25.md。

坐标约定
--------
板面 y=0，板沿 x 反向恒速 v；手端 H=(x_H,y_H)，杆长 L，与向下法线夹角 θ；
尖端 T = H + L(sinθ, -cosθ)。法向压入 δ = L cosθ - y_H ≥ 0；v_rel = ẋ_T + v。
手端弹簧 k_x, k_y（可加阻尼 c_x, c_y）、等效弯曲弹簧 k_θ、单边法向罚接触 k_n。
摩擦模式：
  · "lugre"              LuGre 完整形式（σ0,σ1,σ2）
  · "lugre_springslider" LuGre 取 σ1=σ2=0（等价弹簧-滑块，给出干净黏着）
  · "tanh"               μ_k N tanh(v_rel/ε)（无静摩擦，只作对照）
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

L = 0.08
RADIUS = 3.5e-3
RHO = 2300.0
M = RHO * np.pi * RADIUS**2 * L
I_C = M * L**2 / 12.0
GRAV = 9.8

DEFAULT = dict(
    n0=4.0, theta0_deg=20.0,
    kx=300.0, ky=1.0e6, ktheta=1000.0,
    cx=0.0, cy=0.0,
    kn=2.0e4, cn=1.0,
    friction="lugre_springslider",
    mu_s=0.62, mu_k=0.51, vs=0.01, eps_tanh=0.02,
    sigma0=1.0e5, sigma1=0.0, sigma2=0.0,
    v_board=0.5,
)


def _friction(p, N, v_rel, z, n0):
    """返回 (F_f 沿 +v_rel 方向, ż)。"""
    if N <= 1e-9 * n0:
        return 0.0, 0.0
    if p["friction"] == "tanh":
        return p["mu_k"] * N * np.tanh(v_rel / p["eps_tanh"]), 0.0
    mu_s, mu_k, vs = p["mu_s"], p["mu_k"], p["vs"]
    sig0 = p["sigma0"]
    g = N * (mu_k + (mu_s - mu_k) * np.exp(-(v_rel / vs) ** 2))
    zdot = v_rel - sig0 * abs(v_rel) * z / max(g, 1e-12)
    F = sig0 * z + p["sigma1"] * zdot + p["sigma2"] * v_rel
    return F, zdot


def rhs_factory(p):
    """返回 (rhs, y0, meta)。"""
    kx, ky, kth = p["kx"], p["ky"], p["ktheta"]
    cx, cy, cn = p["cx"], p["cy"], p["cn"]
    kn = p["kn"]
    v_board = p["v_board"]
    th0 = np.deg2rad(p["theta0_deg"])
    n0 = p["n0"]
    dz0 = n0 / kn
    y0 = L * np.cos(th0) - dz0
    x0 = 0.0
    y_anchor = y0 + (n0 + M * GRAV) / ky
    x_anchor = 0.0
    v0 = v_board
    F0, _ = _friction(p, n0, v0, 0.0, n0)
    # LuGre 初值取滑动稳态；tanh 无内态
    if p["friction"] == "tanh":
        z0 = 0.0
    else:
        g0 = n0 * (p["mu_k"] + (p["mu_s"] - p["mu_k"]) *
                   np.exp(-(v0 / p["vs"]) ** 2))
        z0 = g0 / p["sigma0"]

    def rhs(t, s):
        x, y, th, vx, vy, w, z = s
        cth, sth = np.cos(th), np.sin(th)
        delta = L * cth - y
        vx_tip = vx + L * cth * w
        v_rel = vx_tip + v_board
        delta_dot = -L * sth * w - vy
        if delta > 0.0:
            N = max(kn * delta + cn * delta_dot, 0.0)
        else:
            N = 0.0
        F, zdot = _friction(p, N, v_rel, z, n0)
        f = -F
        Qx = -kx * (x - x_anchor) - cx * vx + f
        Qy = -ky * (y - y_anchor) - cy * vy - M * GRAV + N
        Qth = (-kth * (th - th0) - M * GRAV * L / 2.0 * sth
               + f * L * cth + N * L * sth)
        Mmat = np.array([
            [M, 0.0, M * L / 2.0 * cth],
            [0.0, M, M * L / 2.0 * sth],
            [M * L / 2.0 * cth, M * L / 2.0 * sth, I_C + M * L**2 / 4.0],
        ])
        rhs_vec = np.array([
            Qx + M * L / 2.0 * w * w * sth,
            Qy - M * L / 2.0 * w * w * cth,
            Qth,
        ])
        acc = np.linalg.solve(Mmat, rhs_vec)
        return [vx, vy, w, acc[0], acc[1], acc[2], zdot]

    meta = dict(x_anchor=x_anchor, y_anchor=y_anchor, th0=th0,
                n0=n0, v_board=v_board, z0=z0)
    return rhs, [x0, y0, th0, 0.0, 0.0, 0.0, z0], meta


def simulate(p, t_max=0.004, rtol=1e-8, max_step=1e-5, npts=4001,
             stop_on_contact_loss=True, v_stop_tol=1e-3, t_ignore=5e-5):
    """积分并返回 (info, traj)。

    traj 额外包含：
      t_stop / i_stop : 首个“停止点”（v_rel 过零或 |v_rel|<v_stop_tol，t>t_ignore）
      slip            : slice(0, i_stop+1)（用于单次划痕的 I(ξ)）
    """
    rhs, y0, meta = rhs_factory(p)

    def ev_contact(t, s):
        return L * np.cos(s[2]) - s[1]

    ev_contact.terminal = True
    ev_contact.direction = -1
    sol = solve_ivp(rhs, (0.0, t_max), y0, method="DOP853", rtol=rtol,
                    atol=1e-10, events=ev_contact, dense_output=True,
                    max_step=max_step)
    tt = np.linspace(0.0, sol.t[-1], npts)
    Y = sol.sol(tt)
    x, y, th, vx, vy, w, z = Y
    delta = L * np.cos(th) - y
    N = np.maximum(p["kn"] * delta, 0.0)
    vx_tip = vx + L * np.cos(th) * w
    v_rel = vx_tip + p["v_board"]
    F = np.zeros_like(v_rel)
    zdot = np.zeros_like(v_rel)
    n0 = meta["n0"]
    ok = N > 1e-9 * n0
    if p["friction"] == "tanh":
        F[ok] = p["mu_k"] * N[ok] * np.tanh(v_rel[ok] / p["eps_tanh"])
    else:
        g = N[ok] * (p["mu_k"] + (p["mu_s"] - p["mu_k"]) *
                     np.exp(-(v_rel[ok] / p["vs"]) ** 2))
        zdot[ok] = v_rel[ok] - p["sigma0"] * np.abs(v_rel[ok]) * z[ok] / g
        F[ok] = (p["sigma0"] * z[ok] + p["sigma1"] * zdot[ok]
                 + p["sigma2"] * v_rel[ok])
    xi = np.concatenate([[0.0],
                         np.cumsum(0.5 * (v_rel[:-1] + v_rel[1:]) * np.diff(tt))])
    q = np.abs(F * v_rel)

    # 停止点：首个 v_rel 过零 / 接近零
    i_stop = len(tt) - 1
    for i in range(1, len(tt)):
        if tt[i] > t_ignore and (v_rel[i] <= 0.0 or abs(v_rel[i]) < v_stop_tol):
            i_stop = i
            break
    info = dict(
        t_end_ms=float(tt[-1] * 1e3),
        touch_end=bool(sol.t_events[0].size > 0),
        t_stop_ms=float(tt[i_stop] * 1e3),
        i_stop=int(i_stop),
        N_min=float(N[:i_stop + 1].min()),
        N_max=float(N[:i_stop + 1].max()),
        N_mean=float(N[:i_stop + 1].mean()),
        N_rel_std=float(N[:i_stop + 1].std() / max(N[:i_stop + 1].mean(), 1e-12)),
        v_rel_min=float(v_rel[:i_stop + 1].min()),
        v_rel_max=float(v_rel[:i_stop + 1].max()),
        stick_frac=float(np.mean(np.abs(v_rel[:i_stop + 1]) < v_stop_tol)),
        xi_span_mm=float((xi[:i_stop + 1].max() - xi[:i_stop + 1].min()) * 1e3),
        xi_end_mm=float(xi[i_stop] * 1e3),
    )
    traj = dict(tt=tt, N=N, v_rel=v_rel, delta=delta, F_f=F, xi=xi, q=q,
                th=th, x=x, y=y, z=z, i_stop=i_stop, meta=meta)
    return info, traj


def dash_profile(traj, nbins=120, sigma_bins=1.0):
    """取 traj 的滑移段 [0, i_stop]，返回 I(ξ) 剖面与截断指标。"""
    i = traj["i_stop"]
    sl = slice(0, i + 1)
    t = traj["tt"][sl]
    q = traj["q"][sl]
    xi = traj["xi"][sl]
    # 延迟导入，避免依赖顺序
    import os
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from deposit_profile import deposit_profile, profile_integral, total_deposit
    xc, I = deposit_profile(t, q, xi=xi, nbins=nbins, sigma_bins=sigma_bins)
    I2 = total_deposit(t, q)
    metrics = dict(n_bins=int(xc.size), I2=float(I2),
                   I_conservation=float(profile_integral(xc, I) / I2) if I2 > 0 else np.nan)
    if xc.size >= 10:
        n = len(I)
        mid = I[int(0.2 * n):int(0.8 * n)]
        plateau = float(np.median(mid))
        I_c = float(np.mean(I[-max(3, n // 20):]))
        metrics.update(
            I_plateau=plateau,
            I_c=I_c,
            cutoff_ratio=float(I_c / plateau) if plateau > 0 else np.nan,
            xi_c_mm=float(xc[-1] * 1e3),
            I_max=float(np.max(I)),
        )
    else:
        metrics.update(I_plateau=np.nan, I_c=np.nan, cutoff_ratio=np.nan,
                       xi_c_mm=np.nan, I_max=np.nan)
    return xc, I, metrics
