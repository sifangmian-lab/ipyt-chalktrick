# -*- coding: utf-8 -*-
"""仿真 S4：v2 最小模型（弹性杆 + 手端正交弹簧）基线

模型（3 自由度集总，详见 docs/研究记录/模型v2规格_2026-10-25.md）：
  q = (x_H, y_H, θ)：手端两个正交弹簧 k_x, k_y 连到固定支架，
  杆为刚体 + 等效弯曲弹簧 k_θ = 3EI/L，板面接触用罚法向弹簧 k_n（N≥0），
  切向用 LuGre（含 μ_s 峰），板以恒速 v 沿 -x 运动。

本脚本目的不是拟合实验，而是 P0-3 的“最小实现 + 基线”：
  1) 验证 N(t) 在部分参数区可以呈平台（对比 v1 的几何卸载）；
  2) 给出 I(ξ) 的统一后处理结果；
  3) 暴露当前模型在“干净平台 + 单次截断”上仍缺什么（见结果说明）。
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as _fm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from deposit_profile import deposit_profile, profile_integral, total_deposit

for _p in ("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
           "/usr/share/fonts/opentype/noto/NotoSansCJK-Bold.ttc",
           r"C:\Windows\Fonts\msyh.ttc", r"C:\Windows\Fonts\msyhbd.ttc",
           r"C:\Windows\Fonts\simhei.ttf"):
    if os.path.exists(_p):
        try:
            _fm.fontManager.addfont(_p)
        except Exception:
            pass
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "PingFang SC",
                                   "Noto Sans CJK SC", "Noto Sans CJK JP",
                                   "Arial Unicode MS", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

# ---------------- 几何与材料（与 v1 一致的粉笔参数） ----------------
L = 0.08
RADIUS = 3.5e-3
RHO = 2300.0
M = RHO * np.pi * RADIUS**2 * L
I_C = M * L**2 / 12.0
GRAV = 9.8

# ---------------- 名义 v2 参数 ----------------
P = dict(
    n0=4.0, theta0_deg=20.0,
    kx=300.0, ky=1.0e6, ktheta=1000.0,   # 软切向、硬法向/转动（候选“平台”区）
    cx=0.0, cy=0.0,
    kn=2.0e4, cn=1.0,
    mu_s=0.62, mu_k=0.51, vs=0.01,
    sigma0=1.0e5, sigma1=30.0, sigma2=5.0,
    v_board=0.5,
)


def rhs_factory(p):
    kx, ky, kth = p["kx"], p["ky"], p["ktheta"]
    cx, cy, cn = p["cx"], p["cy"], p["cn"]
    kn = p["kn"]
    mu_s, mu_k, vs = p["mu_s"], p["mu_k"], p["vs"]
    sig0, sig1, sig2 = p["sigma0"], p["sigma1"], p["sigma2"]
    v_board = p["v_board"]
    th0 = np.deg2rad(p["theta0_deg"])
    n0 = p["n0"]
    # 初始：给定法向压入 δ0=n0/kn；手端 y 由接触几何定；y 弹簧预载 (n0+mg)/ky
    dz0 = n0 / kn
    y0 = L * np.cos(th0) - dz0
    x0 = 0.0
    y_anchor = y0 + (n0 + M * GRAV) / ky
    x_anchor = 0.0
    v0 = v_board
    g0 = n0 * (mu_k + (mu_s - mu_k) * np.exp(-(v0 / vs) ** 2))
    z0 = g0 / sig0

    def rhs(t, s):
        x, y, th, vx, vy, w, z = s
        cth, sth = np.cos(th), np.sin(th)
        delta = L * cth - y                  # 法向压入量 δ ≥ 0
        vx_tip = vx + L * cth * w
        v_rel = vx_tip + v_board             # 尖端相对板面速度
        delta_dot = -L * sth * w - vy
        if delta > 0.0:
            N = max(kn * delta + cn * delta_dot, 0.0)
        else:
            N = 0.0
        if N > 1e-9 * n0:
            g = N * (mu_k + (mu_s - mu_k) * np.exp(-(v_rel / vs) ** 2))
            zdot = v_rel - sig0 * abs(v_rel) * z / g
            F = sig0 * z + sig1 * zdot + sig2 * v_rel   # 摩擦力量值（沿 +v_rel）
            f = -F                                        # 作用在杆上的切向力
        else:
            zdot = 0.0
            f = 0.0
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

    return rhs, [x0, y0, th0, 0.0, 0.0, 0.0, z0], n0, v_board, th0


def simulate(p, t_max=0.01, rtol=1e-8):
    rhs, s0, n0, v_board, th0 = rhs_factory(p)

    def ev_contact(t, s):
        return L * np.cos(s[2]) - s[1]      # δ=0：失去接触

    ev_contact.terminal = True
    ev_contact.direction = -1
    sol = solve_ivp(rhs, (0.0, t_max), s0, method="DOP853", rtol=rtol,
                    atol=1e-10, events=ev_contact, dense_output=True,
                    max_step=1e-5)
    tt = np.linspace(0.0, sol.t[-1], 4001)
    Y = sol.sol(tt)
    x, y, th, vx, vy, w, z = Y
    delta = L * np.cos(th) - y
    N = np.maximum(p["kn"] * delta, 0.0)
    vx_tip = vx + L * np.cos(th) * w
    v_rel = vx_tip + v_board
    g = N * (p["mu_k"] + (p["mu_s"] - p["mu_k"]) *
             np.exp(-(v_rel / p["vs"]) ** 2))
    zdot = np.zeros_like(v_rel)
    ok = N > 1e-9 * n0
    zdot[ok] = (v_rel[ok] - p["sigma0"] * np.abs(v_rel[ok]) * z[ok] / g[ok])
    F_f = np.zeros_like(v_rel)
    F_f[ok] = (p["sigma0"] * z[ok] + p["sigma1"] * zdot[ok]
               + p["sigma2"] * v_rel[ok])
    xi = np.concatenate([[0.0],
                         np.cumsum(0.5 * (v_rel[:-1] + v_rel[1:]) * np.diff(tt))])
    q = np.abs(F_f * v_rel)
    xc, I = deposit_profile(tt, q, xi=xi, nbins=200, sigma_bins=2.0)
    I_int = profile_integral(xc, I) if xc.size else 0.0
    I2 = total_deposit(tt, q)

    info = dict(
        t_end_ms=float(tt[-1] * 1e3),
        touch_end=bool(sol.t_events[0].size > 0),
        N_min=float(N.min()), N_max=float(N.max()),
        N_mean=float(N.mean()),
        N_rel_std=float(N.std() / max(N.mean(), 1e-12)),
        v_rel_min=float(v_rel.min()), v_rel_max=float(v_rel.max()),
        stick_frac=float(np.mean(np.abs(v_rel) < 1e-3)),
        xi_span_mm=float((xi.max() - xi.min()) * 1e3),
        I2=float(I2),
        I_int=float(I_int),
        I_conservation=float(I_int / I2) if I2 > 0 else np.nan,
        n_profile_bins=int(xc.size),
    )
    traj = dict(tt=tt, N=N, v_rel=v_rel, delta=delta, F_f=F_f, xi=xi,
                xc=xc, I=I, th=th, x=x, y=y)
    return info, traj


def plot_case(tag, info, tr, outdir):
    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    t = tr["tt"] * 1e3
    ax[0, 0].plot(t, tr["N"], "g-", lw=1.6)
    ax[0, 0].set_xlabel("t / ms"); ax[0, 0].set_ylabel("N / N")
    ax[0, 0].set_title(f"({tag}a) 法向力 N(t)：{'平台' if info['N_rel_std'] < 0.1 else '明显变化'}")
    ax[0, 0].grid(alpha=0.3)
    ax[0, 1].plot(t, tr["v_rel"], "b-", lw=1.4)
    ax[0, 1].axhline(0.0, color="k", lw=0.8, ls=":")
    ax[0, 1].set_xlabel("t / ms"); ax[0, 1].set_ylabel("v_rel / (m/s)")
    ax[0, 1].set_title(f"({tag}b) 相对滑动速度（stick 占比 {info['stick_frac']*100:.1f}%）")
    ax[0, 1].grid(alpha=0.3)
    ax[1, 0].plot(t, tr["delta"] * 1e3, "m-", lw=1.4)
    ax[1, 0].set_xlabel("t / ms"); ax[1, 0].set_ylabel("δ / mm")
    ax[1, 0].set_title(f"({tag}c) 法向压入量 δ(t)（N = k_n δ）")
    ax[1, 0].grid(alpha=0.3)
    if tr["xc"].size:
        ax[1, 1].plot(tr["xc"] * 1e3, tr["I"], "r-", lw=1.6)
        ax[1, 1].set_xlabel("ξ / mm"); ax[1, 1].set_ylabel("I(ξ) / (J/m)")
        ax[1, 1].set_title(f"({tag}d) I(ξ)（∫I dξ / I_2nd = {info['I_conservation']:.3f}）")
        ax[1, 1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"S4_{tag}_时域与强度.png"), dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="./S4_v2输出")
    ap.add_argument("--tmax", type=float, default=0.01)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    cases = [
        ("基线A", dict(P)),
        ("基线B", dict(P, kx=1000.0, ky=1.0e5, ktheta=100.0)),
        ("基线C", dict(P, v_board=0.05)),
    ]
    rows = []
    for tag, p in cases:
        info, tr = simulate(p, t_max=args.tmax)
        rows.append(dict(case=tag, **{k: p[k] for k in
                                      ("kx", "ky", "ktheta", "v_board")},
                         **info))
        plot_case(tag, info, tr, args.outdir)
        print(f"[{tag}] t_end={info['t_end_ms']:.2f} ms, "
              f"lost={info['touch_end']}, N={info['N_min']:.2f}-{info['N_max']:.2f}, "
              f"rel_std={info['N_rel_std']:.3f}, stick={info['stick_frac']:.3f}, "
              f"ξ={info['xi_span_mm']:.3f} mm, ∫I/I2={info['I_conservation']:.4f}")

    with open(os.path.join(args.outdir, "S4_基线_结果表.csv"), "w",
              newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\n输出目录：{os.path.abspath(args.outdir)}")
    print("注意：基线A/B 使用 v=0.5 m/s 以直接对照 v1 名义工况；实验窗口 0.01–0.2 m/s 需另扫。")


if __name__ == "__main__":
    main()
