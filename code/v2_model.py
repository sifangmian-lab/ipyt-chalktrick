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
    # 2026-10-25 A2 修正：手端 y 弹簧预载方向。
    # 初态静力平衡应为 F_s - M g + N0 = 0，即 F_s = M g - N0（向下），
    # 故 y_anchor = y0 - (N0 - M g)/ky。原实现取 y0 + (N0 + M g)/ky，
    # 使初始净法向力为 +2N0（向上 8 N），是 S4 平台振铃与 k_n 不收敛的诱因之一。
    y_anchor = y0 - (n0 - M * GRAV) / ky
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


# ---------------------------------------------------------------------------
# 刚性单边接触（LCP / 事件驱动）：k_n -> ∞ 的极限参考
# ---------------------------------------------------------------------------
# 罚接触在 k_n -> ∞ 时，接触力振荡频率 ~ sqrt(k_n/M) 发散而振幅不消失，
# 宏观量没有唯一极限（S5 的"k_n 刚性极限不收敛"即此）。这里改用指标 1 的
# 单边约束：接触相内强制 C = y_T = y - L cosθ = 0，接触力 N 作为约束乘子
# 由约束加速度 C̈ = 0 解出；N < 0 时以事件释放（自由飞行相暂不续接）。
# 本模式只支持 spring-slider 本构（σ1 = σ2 = 0，F_f = σ0 z，与 N 解耦），
# 供 S5b 与完整周期模型作极限参考；tanh/粘性 LuGre 仍走罚接触。

def _mass_matrix(th):
    cth, sth = np.cos(th), np.sin(th)
    return np.array([
        [M, 0.0, M * L / 2.0 * cth],
        [0.0, M, M * L / 2.0 * sth],
        [M * L / 2.0 * cth, M * L / 2.0 * sth, I_C + M * L**2 / 4.0],
    ])


def _coriolis_vec(th, w):
    sth, cth = np.sin(th), np.cos(th)
    return np.array([M * L / 2.0 * w * w * sth,
                     -M * L / 2.0 * w * w * cth, 0.0])


def _rigid_terms(p, s, y_anchor):
    """刚性接触相：状态 s=(x, th, vx, w, z) -> (acc0, accN, N, aux)。

    接触约束 C = y - L cosθ = 0，故 y = L cosθ、vy = -L sinθ w。
    N 由 C̈ = a_y + L sinθ α + L cosθ w² = 0 解出：
        a = M^{-1}(Q0 + N J_N) = acc0 + N accN
        N = -(acc0_y + L sinθ acc0_θ + L cosθ w²)
            / (accN_y + L sinθ accN_θ)
    """
    x, th, vx, w, z = s
    cth, sth = np.cos(th), np.sin(th)
    y = L * cth
    vy = -L * sth * w
    v_board = p["v_board"]
    # spring-slider：F = σ0 z，不依赖 N
    F = p["sigma0"] * z
    f = -F
    v_rel = vx + L * cth * w + v_board
    th_ref = np.deg2rad(p["theta0_deg"])
    Q0 = np.array([
        -p["kx"] * x - p["cx"] * vx + f,
        -p["ky"] * (y - y_anchor) - p["cy"] * vy - M * GRAV,
        -p["ktheta"] * (th - th_ref) - M * GRAV * L / 2.0 * sth
        + f * L * cth,
    ])
    Minv = np.linalg.inv(_mass_matrix(th))
    acc0 = Minv @ (Q0 + _coriolis_vec(th, w))
    JN = np.array([0.0, 1.0, L * sth])
    accN = Minv @ JN
    den = accN[1] + L * sth * accN[2]
    N = -(acc0[1] + L * sth * acc0[2] + L * cth * w * w) / den
    acc = acc0 + N * accN
    residual = acc[1] + L * sth * acc[2] + L * cth * w * w
    # 鬃毛内态；N→0 时 g→0，这里截断避免 zdot 奇异
    if N > 1e-12 * max(p["n0"], 1.0):
        g = N * (p["mu_k"] + (p["mu_s"] - p["mu_k"])
                 * np.exp(-(v_rel / p["vs"]) ** 2))
        zdot = v_rel - p["sigma0"] * abs(v_rel) * z / max(g, 1e-12)
    else:
        zdot = 0.0
    aux = dict(acc=acc, zdot=zdot, v_rel=v_rel, F=F, y=y, vy=vy,
               residual=residual)
    return acc0, accN, N, aux


def _rigid_N(p, s, y_anchor):
    return _rigid_terms(p, s, y_anchor)[2]


def rigid_factory(p):
    """刚性接触相的 (rhs, y0, meta)。初态取 C=0 且 N(0)=N0。

    y_anchor 通过"N_anchor 对 anchor 仿射"一次求出，使初始接触力与罚模型
    的名义 N0 一致（δ0 = N0/k_n -> 0 的匹配口径）。
    """
    if p["friction"] != "lugre_springslider":
        raise ValueError("contact='rigid' 目前只支持 friction='lugre_springslider'")
    th0 = np.deg2rad(p["theta0_deg"])
    n0 = p["n0"]
    v_board = p["v_board"]
    g0 = n0 * (p["mu_k"] + (p["mu_s"] - p["mu_k"])
               * np.exp(-(v_board / p["vs"]) ** 2))
    z0 = g0 / p["sigma0"]
    s0 = (0.0, th0, 0.0, 0.0, z0)
    n_anchor0 = _rigid_N(p, s0, 0.0)
    n_anchor1 = _rigid_N(p, s0, 1.0)
    slope = n_anchor1 - n_anchor0
    if abs(slope) < 1e-30:
        raise ValueError("无法由 y_anchor 调整初始法向力")
    y_anchor = (n0 - n_anchor0) / slope
    n_check = _rigid_N(p, s0, y_anchor)
    meta = dict(x_anchor=0.0, y_anchor=float(y_anchor), th0=th0,
                n0=n0, v_board=v_board, z0=z0, contact="rigid",
                n0_check=float(n_check))

    def rhs(t, s):
        _, _, _, aux = _rigid_terms(p, s, y_anchor)
        acc = aux["acc"]
        return [s[2], s[3], acc[0], acc[2], aux["zdot"]]

    return rhs, list(s0), meta


def simulate_rigid(p, t_max=0.004, rtol=1e-8, max_step=1e-5, npts=4001,
                   v_stop_tol=1e-3, t_ignore=5e-5):
    """刚性接触（LCP 极限）积分；接口与 simulate() 相同。

    traj 的 N 为约束乘子（接触力），delta 恒为 0；另记录约束加速度残差。
    """
    rhs, y0, meta = rigid_factory(p)
    y_anchor = meta["y_anchor"]

    def ev_release(t, s):
        return _rigid_N(p, s, y_anchor)

    ev_release.terminal = True
    ev_release.direction = -1
    sol = solve_ivp(rhs, (0.0, t_max), y0, method="DOP853", rtol=rtol,
                    atol=1e-10, events=ev_release, dense_output=True,
                    max_step=max_step)
    tt = np.linspace(0.0, sol.t[-1], npts)
    S = sol.sol(tt)
    x, th, vx, w, z = S
    cth, sth = np.cos(th), np.sin(th)
    y = L * cth
    vy = -L * sth * w
    delta = np.zeros_like(tt)
    v_rel = vx + L * cth * w + p["v_board"]
    N = np.empty_like(tt)
    residual = np.empty_like(tt)
    zdot = np.empty_like(tt)
    for i in range(len(tt)):
        _, _, Ni, aux = _rigid_terms(p, (x[i], th[i], vx[i], w[i], z[i]),
                                     y_anchor)
        N[i] = Ni
        residual[i] = aux["residual"]
        zdot[i] = aux["zdot"]
    F = p["sigma0"] * z
    q = np.abs(F * v_rel)
    xi = np.concatenate([[0.0],
                         np.cumsum(0.5 * (v_rel[:-1] + v_rel[1:]) * np.diff(tt))])

    def first_stop(arr, thr):
        for i in range(1, len(arr)):
            if tt[i] > t_ignore and arr[i] <= thr:
                return i
        return len(arr) - 1

    i_stop = first_stop(v_rel, v_stop_tol)
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
        contact="rigid",
        release_events=int(sol.t_events[0].size),
        constraint_residual_max=(float(np.max(np.abs(residual)))
                                 if residual.size else 0.0),
        n0_check=float(meta["n0_check"]),
    )
    traj = dict(tt=tt, N=N, v_rel=v_rel, delta=delta, F_f=F, xi=xi, q=q,
                th=th, x=x, y=y, z=z, i_stop=i_stop, meta=meta)
    return info, traj


def simulate(p, t_max=0.004, rtol=1e-8, max_step=1e-5, npts=4001,
             stop_on_contact_loss=True, v_stop_tol=1e-3, t_ignore=5e-5):
    """积分并返回 (info, traj)。

    traj 额外包含：
      t_stop / i_stop : 首个“停止点”（v_rel 过零或 |v_rel|<v_stop_tol，t>t_ignore）
      slip            : slice(0, i_stop+1)（用于单次划痕的 I(ξ)）
    """
    if p.get("contact", "penalty") == "rigid":
        return simulate_rigid(p, t_max=t_max, rtol=rtol, max_step=max_step,
                              npts=npts, v_stop_tol=v_stop_tol,
                              t_ignore=t_ignore)
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
