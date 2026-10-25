# -*- coding: utf-8 -*-
"""
仿真 S3：LuGre 鬃毛摩擦模型替换 tanh（本构对比）
==================================================
对应《新参考文献与N下降段方向对照笔记》接缝一"可直接执行的动作 ②"：
以 LuGre 重写摩擦项，对比翻转前后的 N(t)、I(ξ) 差异，检验 S2 结论的"绝对性"。

背景
----
S2 已证 tanh 正则化在窗口 ε∈[1e-3,1e-1] 内稳健，但那只是"换种光滑化都差不多"。
LuGre 模型（Åström & Canudas de Wit 2008）给摩擦真实的切向内态记忆：

  ż      = v_rel − σ₀|v_rel|·z / g(v_rel)
  F_fric = σ₀z + σ₁ż + σ₂v_rel
  g(v)   = N·[μ_k + (μ_s − μ_k)·exp(−(v/v_s)²)]   （Coulomb + Stribeck 峰）

物理上它意味着：
  · "黏"相有真实预滑柔度（Heslot 1994 蠕变图像），tanh 只给出刚性切换；
  · Stribeck 静摩擦峰使翻转点附近摩擦先升后降；
  · 力换向滞后于速度换向（鬃毛记忆），滞后期是本仿真的核心观测量。

对比基准：S2 名义工况 k_c=2e4 N/m, ε=0.02, N0=4 N, ω0=0, v=0.5 m/s, μ_k=0.35。
注意：μ_s 暂取 1.15·μ_k 占位——E1（静摩擦时间依赖）测出真值后必须替换。

输出（--outdir，默认 ./S3_LuGre输出/）
  S3_LuGre_结果表.csv     t_flip / t* / I_2nd / 换向滞后期，tanh 基准 + LuGre 扫描
  S3_LuGre_对比.png       N(t)、F_fric(t)、v_rel(t)、I(ξ) 四联对比
  S3_LuGre_内态.png       鬃毛内态 z(t) 与饱和度（预滑记忆）

数值注意事项（本仿真踩过的坑，均已处理）：
  · 内态初值必须继承滑动稳态 z0=g(v0)/σ₀，否则 σ₁ż 项在 t=0 爆发；
  · 分离点 g=0 处鬃毛项无定义，后处理需掩膜；
  · σ₁ 过大（≥300 N·s/m）超出稳定区，双容差报警，网格上限取 100。

运行：
python 仿真S3_LuGre本构对比.py [--outdir 目录]
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp, cumulative_trapezoid
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

# ---------------- 物理参数（与主模型一致） ----------------
L   = 0.08
h   = 0.075
RADIUS = 3.5e-3
RHO = 2300.0
GRAV = 9.8
M   = RHO * np.pi * RADIUS**2 * L
I_O = M * L**2 / 3.0
MU_K = 0.35           # 动摩擦系数（E1/E2 后替换）
V_BOARD = 0.5
EPS_TANH = 2e-2       # tanh 基准（S2 名义值）

# 工况
KC   = 2.0e4          # S1/S2 名义接触刚度
N0   = 4.0
OMEGA0 = 0.0
T_MAX = 0.05

# LuGre 参数网格
SIGMA0_LIST = [1.0e5, 3.0e5]     # 鬃毛刚度 N/m
SIGMA1_LIST = [0.0, 3.0e1, 1.0e2]  # 鬃毛阻尼 N·s/m（300 已证超出稳定区）
SIGMA2 = 5.0                     # 黏性项 N·s/m
V_S = 1.0e-2                     # Stribeck 特征速度 m/s
MUS_RATIO = 1.15                 # μ_s / μ_k 占位比（E1 后替换）

RTOL_PAIR = (1e-9, 1e-12)
_trapz = getattr(np, "trapezoid", None) or np.trapz


# ============================== 模型内核 ==============================
def delta(th):
    return L * np.cos(th) - h


def contact_omega(kc, th):
    return float(np.sqrt(max(kc * (L * np.sin(th)) ** 2 / I_O, 1.0)))


# ---- tanh 基准 ----
def rhs_tanh(t, y, mu):
    th, w = y
    d = delta(th)
    if d <= 0.0:
        return [w, -(M * GRAV * L / 2.0) * np.sin(th) / I_O]
    Nn = KC * d
    v_rel = V_BOARD - L * np.cos(th) * w
    tau = (M * GRAV * L / 2.0) * np.sin(th) \
        - Nn * (L * np.sin(th) + mu * h * np.tanh(v_rel / EPS_TANH))
    return [w, -tau / I_O]


# ---- LuGre ----
def make_rhs_lugre(sigma0, sigma1, mu_s):
    mu_k = MU_K
    def rhs(t, y):
        th, w, z = y
        d = delta(th)
        if d <= 0.0:
            return [w, -(M * GRAV * L / 2.0) * np.sin(th) / I_O, 0.0]
        Nn = KC * d
        v_rel = V_BOARD - L * np.cos(th) * w
        g = Nn * (mu_k + (mu_s - mu_k) * np.exp(-(v_rel / V_S) ** 2))
        zdot = v_rel - sigma0 * abs(v_rel) * z / g
        F_f = sigma0 * z + sigma1 * zdot + SIGMA2 * v_rel
        tau = (M * GRAV * L / 2.0) * np.sin(th) \
            - Nn * L * np.sin(th) - F_f * h
        return [w, -tau / I_O, zdot]
    return rhs


def _sep_event(t, y):
    return delta(y[0])


_sep_event.terminal = True
_sep_event.direction = -1


# ============================== 积分 ==============================
def run(model, rtol=1e-10, t_max=T_MAX, **pars):
    """model ∈ {'tanh','lugre'}；返回 (info, traj)。"""
    th0 = np.arccos(np.clip((h + N0 / KC) / L, -1.0, 1.0))
    if model == "lugre":
        # 继承前半个周期滑动建立的鬃毛内态：z0 = g(v0)/σ₀
        # （z0=0 会让 ż(0)=v(0)≈0.5，σ₁ż 项瞬间爆 ~σ₁v 牛顿，属数值事故）
        Nn0 = KC * delta(th0)
        v0 = V_BOARD - L * np.cos(th0) * OMEGA0
        g0 = Nn0 * (MU_K + (pars["mu_s"] - MU_K) * np.exp(-(v0 / V_S) ** 2))
        y0 = [th0, OMEGA0, g0 / pars["sigma0"]]
    else:
        y0 = [th0, OMEGA0]
    if model == "tanh":
        f = lambda t, y: rhs_tanh(t, y, MU_K)
    else:
        f = make_rhs_lugre(pars["sigma0"], pars["sigma1"], pars["mu_s"])
    sol = solve_ivp(f, (0.0, t_max), y0, method="DOP853", rtol=rtol,
                    atol=1e-12, events=_sep_event, dense_output=True,
                    max_step=min(5e-5, 0.05 * 2 * np.pi / contact_omega(KC, th0)))
    separated = len(sol.t_events[0]) > 0
    t_end = float(sol.t_events[0][0]) if separated else t_max
    tt = np.linspace(0.0, t_end, 4001)
    Y = sol.sol(tt)
    th, w = Y[0], Y[1]
    Nn = KC * np.clip(delta(th), 0.0, None)
    v_rel = V_BOARD - L * np.cos(th) * w
    if model == "tanh":
        F_f = MU_K * Nn * np.tanh(v_rel / EPS_TANH)
    else:
        z = Y[2]
        g = Nn * (MU_K + (pars["mu_s"] - MU_K) * np.exp(-(v_rel / V_S) ** 2))
        F_f = np.zeros_like(v_rel)
        # 分离点附近 Nn→0 ⇒ g→0，鬃毛项无定义；必须按相对接触力掩膜，
        # 只用 g>0 会 retained 极小 g，把 z/g 与 F_f 放大到非物理量级。
        ok = Nn > 1e-6 * N0
        zdot = np.zeros_like(v_rel)
        zdot[ok] = v_rel[ok] - pars["sigma0"] * np.abs(v_rel[ok]) * z[ok] / g[ok]
        F_f[ok] = (pars["sigma0"] * z[ok] + pars["sigma1"] * zdot[ok]
                   + SIGMA2 * v_rel[ok])

    # ---- 翻转与换向滞后期 ----
    def first_cross(arr):
        idx = np.nonzero(arr < 0.0)[0]
        if not idx.size:
            return False, np.nan
        i = int(idx[0])
        tc = float(tt[i - 1] + (tt[i] - tt[i - 1]) * arr[i - 1]
                   / (arr[i - 1] - arr[i])) if i > 0 else float(tt[0])
        return True, tc

    flip_v, t_flip_v = first_cross(v_rel) if v_rel[0] > 0 else (True, 0.0)
    flip_F, t_flip_F = first_cross(F_f) if F_f[0] > 0 else (True, 0.0)
    lag = (t_flip_F - t_flip_v) if (flip_v and flip_F) else np.nan

    I2 = float(_trapz(np.abs(F_f * v_rel), tt))
    ell = float(_trapz(np.abs(v_rel), tt))
    xi = cumulative_trapezoid(v_rel, tt, initial=0.0)
    info = dict(model=model, t_star_ms=t_end * 1e3,
                t_flip_ms=(t_flip_v * 1e3 if flip_v else np.nan),
                flip=flip_v, t_flipF_ms=(t_flip_F * 1e3 if flip_F else np.nan),
                lag_ms=(lag * 1e3 if np.isfinite(lag) else np.nan),
                I_2nd=I2, ell_mm=ell * 1e3,
                omega_star=float(w[-1]),
                xdot_max=float(np.max(L * np.cos(th) * np.abs(w))))
    traj = dict(tt=tt, th=th, w=w, Nn=Nn, v_rel=v_rel, F_f=F_f, xi=xi)
    if model == "lugre":
        traj["z_raw"] = Y[2]
        info.update(sigma0=pars["sigma0"], sigma1=pars["sigma1"],
                    mu_s=pars["mu_s"])
    return info, traj


def run_dual(model, **pars):
    a, traj = run(model, rtol=RTOL_PAIR[0], **pars)
    b, _ = run(model, rtol=RTOL_PAIR[1], **pars)
    a["num_err_t_ms"] = abs(a["t_star_ms"] - b["t_star_ms"])
    return a, traj


# ============================== 输出 ==============================
FIELDNAMES = ["model", "sigma0", "sigma1", "mu_s", "t_star_ms", "t_flip_ms",
              "flip", "t_flipF_ms", "lag_ms", "I_2nd", "ell_mm",
              "omega_star", "xdot_max", "num_err_t_ms"]


def save_csv(rows, path):
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=FIELDNAMES, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if (isinstance(v, float) and not np.isfinite(v)) else v)
                        for k, v in r.items() if k in FIELDNAMES})


def plot_compare(rows_traj, outdir):
    tanh_tr = rows_traj[("tanh",)]
    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    # (a) N(t)
    ax[0, 0].plot(tanh_tr["tt"] * 1e3, tanh_tr["Nn"], "k--", lw=1.8,
                  label="tanh 基准 (ε=0.02)")
    for (key, tr), c in zip(rows_traj.items(),
                            plt.cm.viridis(np.linspace(0, 0.85, len(rows_traj)))):
        if key[0] == "lugre":
            ax[0, 0].plot(tr["tt"] * 1e3, tr["Nn"], color=c, lw=1.4,
                          label=f"LuGre σ0={key[1]:.0e}, σ1={key[2]:.0e}")
    ax[0, 0].set_xlabel("t / ms"); ax[0, 0].set_ylabel("N / N")
    ax[0, 0].set_title("(a) N(t)：本构对比"); ax[0, 0].legend(fontsize=7)
    ax[0, 0].grid(alpha=0.3)
    # (b) F_fric(t)
    ax[0, 1].plot(tanh_tr["tt"] * 1e3, tanh_tr["F_f"], "k--", lw=1.8,
                  label="tanh 基准")
    for (key, tr), c in zip(rows_traj.items(),
                            plt.cm.viridis(np.linspace(0, 0.85, len(rows_traj)))):
        if key[0] == "lugre":
            ax[0, 1].plot(tr["tt"] * 1e3, tr["F_f"], color=c, lw=1.4,
                          label=f"LuGre σ0={key[1]:.0e}, σ1={key[2]:.0e}")
    ax[0, 1].axhline(0, color="gray", lw=0.8)
    ax[0, 1].set_xlabel("t / ms"); ax[0, 1].set_ylabel("F_fric / N")
    ax[0, 1].set_title("(b) 摩擦力 F_fric(t)：Stribeck 峰与换向滞后")
    ax[0, 1].legend(fontsize=7); ax[0, 1].grid(alpha=0.3)
    # (c) v_rel(t)：翻转时刻对比
    ax[1, 0].plot(tanh_tr["tt"] * 1e3, tanh_tr["v_rel"], "k--", lw=1.8,
                  label="tanh 基准")
    for (key, tr), c in zip(rows_traj.items(),
                            plt.cm.viridis(np.linspace(0, 0.85, len(rows_traj)))):
        if key[0] == "lugre":
            ax[1, 0].plot(tr["tt"] * 1e3, tr["v_rel"], color=c, lw=1.4,
                          label=f"LuGre σ0={key[1]:.0e}, σ1={key[2]:.0e}")
    ax[1, 0].axhline(0, color="gray", lw=0.8)
    ax[1, 0].set_xlabel("t / ms"); ax[1, 0].set_ylabel("v_rel / (m/s)")
    ax[1, 0].set_title("(c) 相对滑动速度 v_rel(t)：翻转时刻对比")
    ax[1, 0].legend(fontsize=7); ax[1, 0].grid(alpha=0.3)
    # (d) I(ξ) 叠加
    for tr, lab, ls in [(tanh_tr, "tanh 基准", "k--")]:
        _plot_xi(ax[1, 1], tr, lab, ls, "k")
    for (key, tr), c in zip(rows_traj.items(),
                            plt.cm.viridis(np.linspace(0, 0.85, len(rows_traj)))):
        if key[0] == "lugre":
            _plot_xi(ax[1, 1], tr, f"σ0={key[1]:.0e}, σ1={key[2]:.0e}", "-", c)
    ax[1, 1].set_xlabel("ξ / mm"); ax[1, 1].set_ylabel("I(ξ) / (J/m)")
    ax[1, 1].set_title("(d) 空间强度分布 I(ξ) 叠加")
    ax[1, 1].legend(fontsize=7); ax[1, 1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "S3_LuGre_对比.png"), dpi=150)
    plt.close(fig)


def _plot_xi(ax, tr, label, ls, color):
    xi, Nn, tt = tr["xi"], tr["Nn"], tr["tt"]
    F_f, v_rel = tr["F_f"], tr["v_rel"]
    q = np.abs(F_f * v_rel)                       # 摩擦功率 = 沉积率（J/s）
    xc, hist = deposit_profile(tt, q, xi=xi, nbins=200, sigma_bins=2.0)
    if xc.size < 2:
        return
    ax.plot(xc * 1e3, hist, ls, color=color, lw=1.5, label=label)


def plot_z(traj_z, outdir):
    """单独画 LuGre 鬃毛内态 z(t)：预滑与记忆。"""
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for (key, tr), c in zip(traj_z.items(),
                            plt.cm.viridis(np.linspace(0, 0.85, len(traj_z)))):
        ax[0].plot(tr["tt"] * 1e3, tr["z"], color=c, lw=1.5,
                   label=f"σ0={key[0]:.0e}, σ1={key[1]:.0e}")
        g = tr["Nn"] * (MU_K + (tr["mu_s"] - MU_K) *
                        np.exp(-(tr["v_rel"] / V_S) ** 2))
        sat = np.zeros_like(tr["z"])
        ok = g > 0.0
        sat[ok] = tr["z"][ok] * key[0] / g[ok]
        ax[1].plot(tr["tt"] * 1e3, sat, color=c, lw=1.5,
                   label=f"σ0={key[0]:.0e}, σ1={key[1]:.0e}")
    ax[0].set_xlabel("t / ms"); ax[0].set_ylabel("z / m")
    ax[0].set_title("鬃毛偏转内态 z(t)（预滑记忆）")
    ax[1].set_xlabel("t / ms"); ax[1].set_ylabel("σ₀z / g(v)（饱和≈1）")
    ax[1].set_title("内态饱和度：黏滞→滑移转变")
    for a in ax:
        a.grid(alpha=0.3); a.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "S3_LuGre_内态.png"), dpi=150)
    plt.close(fig)


# ============================== 主流程 ==============================
def scan(outdir, want_plots=True):
    os.makedirs(outdir, exist_ok=True)
    rows, trajs, traj_z = [], {}, {}
    r, tr = run_dual("tanh")
    rows.append(r); trajs[("tanh",)] = tr
    print("  [S3] tanh 基准完成")
    mu_s = MUS_RATIO * MU_K
    for s0 in SIGMA0_LIST:
        for s1 in SIGMA1_LIST:
            r, tr = run_dual("lugre", sigma0=s0, sigma1=s1, mu_s=mu_s)
            rows.append(r); trajs[("lugre", s0, s1)] = tr
            traj_z[(s0, s1)] = dict(tt=tr["tt"], z=tr["z_raw"], Nn=tr["Nn"],
                                    v_rel=tr["v_rel"], mu_s=mu_s)
            print(f"  [S3] LuGre σ0={s0:.0e}, σ1={s1:.0e} 完成")
    save_csv(rows, os.path.join(outdir, "S3_LuGre_结果表.csv"))
    if want_plots:
        plot_compare(trajs, outdir)
        plot_z(traj_z, outdir)
    print_summary(rows, outdir)
    return rows


def print_summary(rows, outdir):
    print("=" * 76)
    print("S3：LuGre 本构对比")
    print(f"基准工况：k_c={KC:.3g} N/m, N0={N0} N, ω0={OMEGA0}, v={V_BOARD} m/s, "
          f"μ_k={MU_K}, μ_s={MUS_RATIO}·μ_k（占位）, v_s={V_S} m/s")
    print("=" * 76)
    hdr = (f"{'模型':>22} {'t*/ms':>9} {'tflip/ms':>9} {'t_flipF/ms':>10} "
           f"{'滞后/ms':>9} {'I2nd/J':>10} {'数差/ms':>9}")
    print(hdr); print("-" * 76)
    for r in rows:
        name = ("tanh(ε=0.02)" if r["model"] == "tanh"
                else f"LuGre σ0={r['sigma0']:.0e},σ1={r['sigma1']:.0e}")
        def fmt(v, d=4):
            return f"{v:.{d}g}" if np.isfinite(v) else "—"
        print(f"{name:>22} {fmt(r['t_star_ms'],4):>9} {fmt(r['t_flip_ms'],4):>9} "
              f"{fmt(r['t_flipF_ms'],4):>10} {fmt(r['lag_ms'],3):>9} "
              f"{fmt(r['I_2nd'],3):>10} {fmt(r['num_err_t_ms'],2):>9}")
    print("-" * 76)
    base = rows[0]
    print("判读提示：")
    for r in rows[1:]:
        name = f"σ0={r['sigma0']:.0e},σ1={r['sigma1']:.0e}"
        dt = ((r["t_flip_ms"] - base["t_flip_ms"]) / base["t_flip_ms"] * 100
              if np.isfinite(r["t_flip_ms"]) and np.isfinite(base["t_flip_ms"])
              else np.nan)
        dI = ((r["I_2nd"] - base["I_2nd"]) / base["I_2nd"] * 100
              if r["I_2nd"] and base["I_2nd"] else np.nan)
        print(f"  · LuGre({name}) vs tanh：t_flip 偏差 {dt:+.1f}%，"
              f"I_2nd 偏差 {dI:+.1f}%，换向滞后 {r['lag_ms']:.4f} ms")
    print(f"\n输出目录：{os.path.abspath(outdir)}")


def main():
    ap = argparse.ArgumentParser(description="S3：LuGre 本构对比")
    ap.add_argument("--outdir", default="./S3_LuGre输出")
    ap.add_argument("--no-plots", action="store_true")
    args = ap.parse_args()
    scan(args.outdir, want_plots=not args.no_plots)


if __name__ == "__main__":
    main()
