# -*- coding: utf-8 -*-
"""
仿真 S2：ε 灵敏度扫描（tanh 正则化速度尺度的稳健性）
=====================================================
对应《新参考文献与N下降段方向对照笔记》接缝一"可直接执行的动作 ①"：
评估预言③（摩擦翻转及强度分布折叠尖峰）对 ε 的依赖。

背景
----
摩擦方向用 tanh(v_rel/ε) 正则化代替 sign(v_rel)，ε = 0.02 m/s 是人为选取的。
本扫描在 k_c ∈ {1e3, 1e4, 2e4, 3.16e4} N/m（S1 已证该区间存在摩擦翻转；
k_c ≥ 1e5 翻转消失，不属于本扫描）上变化 ε = 1e-4 … 1 m/s，检验：

  t_flip   摩擦翻转时刻（预言③的核心量）
  t*       分离时间（预期对 ε 不敏感，作对照）
  I_2nd    本段笔痕强度 μ∫N|v_rel|dt
  ℓ        痕迹长
  I(ξ)     空间强度分布（末端回划折叠尖峰的 ε 敏感性，名义 k_c=2e4）

判读：
  · t_flip 在 ε 跨 2 个数量级内变化 <10% → 正则化稳健，预言③可直接对照实验；
  · 强漂移或翻转消失 → tanh 描述不可靠，须升级 LuGre 鬃毛模型（仿真 S3）。

数值要点：ε → 0 时 tanh 切换变陡（局部刚性），双容差复算（rtol=1e-9/1e-12）
的差值随 ε 的变化即"数值可信下限"，防止把积分误差当作物理敏感性。

输出（--outdir，默认 ./S2_ε扫描输出/）
  S2_ε扫描_结果表.csv     全部输出量 + 相对名义 ε 的偏差 + 诊断列
  S2_ε扫描_收敛曲线.png   t_flip、t*、I_2nd、ℓ 对 ε（4 条 k_c 曲线）
  S2_ε扫描_强度分布.png   名义 k_c=2e4 的 I(ξ) 曲线族（折叠尖峰敏感性）
  S2_ε扫描_诊断.png       相对偏差（±10% 带）与数值误差随 ε 的变化

运行：
python 仿真S2_ε灵敏度扫描.py [--kc 1e3 1e4 2e4 3.16e4] [--outdir 目录]
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
MU_DEFAULT = 0.35
V_BOARD = 0.5
EPS_NOM = 2e-2          # 名义正则化速度尺度 m/s

TH_STAR = np.arccos(h / L)

# ---------------- 扫描网格 ----------------
KC_SCAN   = [1.0e3, 1.0e4, 2.0e4, 3.162278e4]   # 翻转存在区间（S1 结论）
EPS_LIST  = [1e-4, 1e-3, 3e-3, 1e-2, 2e-2, 5e-2, 1e-1, 3e-1, 1.0]
N0        = 4.0          # 交接压力 N
OMEGA0    = 0.0          # 交接角速度 rad/s
KC_PROFILE = 2.0e4       # 记录 I(ξ) 的名义 k_c
EPS_PROFILE = [1e-4, 1e-3, 2e-2, 1e-1]   # 记录 I(ξ) 的 ε 子集
T_MAX     = 0.05
RTOL_PAIR = (1e-9, 1e-12)

_trapz = getattr(np, "trapezoid", None) or np.trapz


# ============================== 模型内核 ==============================
def delta(th):
    return L * np.cos(th) - h


def contact_omega(kc, th):
    return float(np.sqrt(max(kc * (L * np.sin(th)) ** 2 / I_O, 1.0)))


def make_rhs(kc, mu, eps):
    def rhs(t, y):
        th, w = y
        d = delta(th)
        if d <= 0.0:
            return [w, -(M * GRAV * L / 2.0) * np.sin(th) / I_O]
        Nn = kc * d
        v_rel = V_BOARD - L * np.cos(th) * w
        tau = (M * GRAV * L / 2.0) * np.sin(th) \
            - Nn * (L * np.sin(th) + mu * h * np.tanh(v_rel / eps))
        return [w, -tau / I_O]
    return rhs


def _sep_event(t, y):
    return delta(y[0])


_sep_event.terminal = True
_sep_event.direction = -1


# ============================== 单段积分 ==============================
def run_segment(kc, N0, omega0=0.0, mu=MU_DEFAULT, eps=EPS_NOM, rtol=1e-10,
                t_max=T_MAX, keep_traj=False):
    """积分单段 [t0, t*]，返回 (info, traj|None)。eps 为扫描变量。"""
    th0 = np.arccos(np.clip((h + N0 / kc) / L, -1.0, 1.0))
    sol = solve_ivp(make_rhs(kc, mu, eps), (0.0, t_max), [th0, omega0],
                    method="DOP853", rtol=rtol, atol=1e-12,
                    events=_sep_event, dense_output=True,
                    max_step=min(5e-5, 0.05 * 2 * np.pi / contact_omega(kc, th0)))
    separated = len(sol.t_events[0]) > 0
    t_end = float(sol.t_events[0][0]) if separated else t_max

    tt = np.linspace(0.0, t_end, 4001)
    th, w = sol.sol(tt)
    Nn = kc * np.clip(delta(th), 0.0, None)
    v_rel = V_BOARD - L * np.cos(th) * w

    # ---- 翻转检测 ----
    if v_rel[0] > 0.0:
        idx = np.nonzero(v_rel < 0.0)[0]
        flip = idx.size > 0
        if flip:
            i = int(idx[0])
            t_flip = float(tt[i - 1] + (tt[i] - tt[i - 1]) * v_rel[i - 1]
                           / (v_rel[i - 1] - v_rel[i])) if i > 0 else float(tt[0])
        else:
            t_flip = np.nan
    else:
        flip, t_flip = True, 0.0

    I2 = mu * float(_trapz(Nn * np.abs(v_rel), tt))
    ell = float(_trapz(np.abs(v_rel), tt))
    if flip:
        m0 = tt <= t_flip
        I_start = mu * float(_trapz(Nn[m0] * np.abs(v_rel[m0]), tt[m0]))
    else:
        I_start = I2
    xdot = L * np.cos(th) * np.abs(w)

    info = dict(kc=kc, eps=eps, N0=N0, omega0=omega0,
                t_star_ms=t_end * 1e3,
                t_flip_ms=(t_flip * 1e3 if flip else np.nan),
                flip=flip, I_2nd=I2, I_start=I_start,
                I_slide=I2 - I_start, ell_mm=ell * 1e3,
                omega_star=float(w[-1]), xdot_max=float(xdot.max()),
                separated=separated)
    traj = None
    if keep_traj:
        xi = cumulative_trapezoid(v_rel, tt, initial=0.0)   # 板面系痕迹坐标
        traj = dict(tt=tt, th=th, w=w, Nn=Nn, v_rel=v_rel, xi=xi)
    return info, traj


def run_dual(kc, N0, omega0, mu, eps):
    a, traj = run_segment(kc, N0, omega0, mu, eps, rtol=RTOL_PAIR[0],
                          keep_traj=True)
    b, _ = run_segment(kc, N0, omega0, mu, eps, rtol=RTOL_PAIR[1])
    for key, col in (("t_star_ms", "num_err_t_ms"), ("t_flip_ms", "num_err_flip_ms"),
                     ("I_2nd", "num_err_I")):
        va, vb = a[key], b[key]
        a[col] = abs(va - vb) if (np.isfinite(va) and np.isfinite(vb)) else np.nan
    return a, traj


# ============================== 诊断 ==============================
def add_diagnostics(rows):
    """对每个 k_c 组（按 ε 排序）补相对名义 ε 的偏差与敏感性列。"""
    groups = {}
    for r in rows:
        groups.setdefault(r["kc"], []).append(r)
    for g in groups.values():
        g.sort(key=lambda r: r["eps"])
        ref = next((r for r in g if abs(r["eps"] - EPS_NOM) < 1e-15), g[len(g) // 2])
        for key in ("t_flip_ms", "t_star_ms", "I_2nd", "ell_mm"):
            rv = ref[key]
            out = "reldev_" + key.replace("_ms", "").replace("_2nd", "")
            for r in g:
                if np.isfinite(rv) and rv != 0.0 and np.isfinite(r[key]):
                    r[out] = (r[key] - rv) / rv
                else:
                    r[out] = np.nan
        # t_flip 对数斜率（相邻 ε 点）
        ys = np.array([r["t_flip_ms"] for r in g], dtype=float)
        xs = np.array([r["eps"] for r in g], dtype=float)
        slopes = [np.nan]
        for i in range(1, len(g)):
            if np.all(np.isfinite([ys[i - 1], ys[i]])) and ys[i - 1] != 0.0:
                s = ((np.log10(abs(ys[i])) - np.log10(abs(ys[i - 1])))
                     / (np.log10(xs[i]) - np.log10(xs[i - 1])))
            else:
                s = np.nan
            slopes.append(s)
        for r, s in zip(g, slopes):
            r["slope_flip"] = s


FIELDNAMES = ["kc", "eps", "N0", "omega0", "t_star_ms", "t_flip_ms", "flip",
              "I_2nd", "I_start", "I_slide", "ell_mm", "omega_star", "xdot_max",
              "separated", "num_err_t_ms", "num_err_flip_ms", "num_err_I",
              "reldev_t_flip", "reldev_t_star", "reldev_I", "reldev_ell",
              "slope_flip"]


def save_csv(rows, path):
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=FIELDNAMES, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if (isinstance(v, float) and not np.isfinite(v)) else v)
                        for k, v in r.items() if k in FIELDNAMES})


# ============================== 作图 ==============================
def plot_convergence(rows, outdir):
    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    specs = [("t_flip_ms", "t_flip / ms", "摩擦翻转时刻", True),
             ("t_star_ms", "t* / ms", "分离时间（对照）", False),
             ("I_2nd", "I_2nd / J", "本段笔痕强度", False),
             ("ell_mm", "ℓ / mm", "痕迹长", False)]
    markers = ["o", "s", "^", "D"]
    for a, (key, ylab, tit, use_err) in zip(ax.flat, specs):
        for kc, mk in zip(KC_SCAN, markers):
            g = sorted([r for r in rows if r["kc"] == kc], key=lambda r: r["eps"])
            xs = np.array([r["eps"] for r in g])
            ys = np.array([r[key] for r in g], dtype=float)
            es = (np.array([r["num_err_flip_ms"] for r in g], dtype=float)
                  if use_err else None)
            a.errorbar(xs, ys, yerr=es, fmt=mk + "-", lw=1.5, ms=5, capsize=3,
                       label=f"k_c={kc:.3g}")
        a.set_xscale("log"); a.set_yscale("log")
        a.set_xlabel("ε / (m/s)"); a.set_ylabel(ylab); a.set_title(tit)
        a.grid(alpha=0.3, which="both"); a.legend(fontsize=8)
    fig.suptitle(f"S2（N0={N0} N, ω0={OMEGA0}）：输出量对正则化尺度 ε 的灵敏度")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "S2_ε扫描_收敛曲线.png"), dpi=150)
    plt.close(fig)


def plot_profiles(profiles, outdir):
    """I(ξ) 空间强度分布：折叠尖峰对 ε 的敏感性（统一后处理，权重 q·dt）。"""
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.6))
    colors = plt.cm.viridis(np.linspace(0, 1, len(profiles)))
    for (eps, tr), c in zip(sorted(profiles.items()), colors):
        xi, Nn, tt, v_rel = tr["xi"], tr["Nn"], tr["tt"], tr["v_rel"]
        q = MU_DEFAULT * Nn * np.abs(v_rel)          # J/s，与主模型式 (11) 一致
        xc, hist = deposit_profile(tt, q, xi=xi, nbins=200, sigma_bins=2.0)
        if xc.size < 2:
            continue
        err = profile_integral(xc, hist) / total_deposit(tt, q) - 1.0
        print(f"  [I(ξ) 守恒] ε={eps:g}: ∫I dξ/I_2nd = {1.0 + err:.4f}")
        ax[0].plot(xc * 1e3, hist, color=c, lw=1.6, label=f"ε={eps:g}")
        if np.nanmax(hist) > 0:
            ax[1].plot(xc * 1e3, hist / np.nanmax(hist), color=c, lw=1.6,
                       label=f"ε={eps:g}")
    ax[0].set_xlabel("ξ / mm（板面系坐标）"); ax[0].set_ylabel("I(ξ) / (J/m)")
    ax[0].set_title(f"空间强度分布（k_c={KC_PROFILE:.3g}，折叠叠加）")
    ax[1].set_xlabel("ξ / mm"); ax[1].set_ylabel("归一化 I(ξ)")
    ax[1].set_title("归一化对比（峰位与峰形）")
    for a in ax:
        a.grid(alpha=0.3); a.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "S2_ε扫描_强度分布.png"), dpi=150)
    plt.close(fig)


def plot_diagnostics(rows, outdir):
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    markers = ["o", "s", "^", "D"]
    for kc, mk in zip(KC_SCAN, markers):
        g = sorted([r for r in rows if r["kc"] == kc], key=lambda r: r["eps"])
        xs = np.array([r["eps"] for r in g])
        yf = np.array([r["reldev_t_flip"] for r in g], dtype=float)
        yt = np.array([r["reldev_t_star"] for r in g], dtype=float)
        ax[0].plot(xs, yf * 100, mk + "-", lw=1.5, ms=5,
                   label=f"t_flip, k_c={kc:.3g}")
        if kc == KC_SCAN[0]:
            ax[0].plot(xs, yt * 100, "k--", lw=1, label="t*（对照, 首组）")
    ax[0].axhline(10, color="gray", ls=":", lw=1)
    ax[0].axhline(-10, color="gray", ls=":", lw=1)
    ax[0].axhline(0, color="k", lw=0.8)
    ax[0].set_xscale("log")
    ax[0].set_xlabel("ε / (m/s)"); ax[0].set_ylabel("相对偏差 / %")
    ax[0].set_title("相对名义 ε=0.02 的偏差（±10% 带）")
    ax[0].legend(fontsize=7); ax[0].grid(alpha=0.3)
    for kc, mk in zip(KC_SCAN, markers):
        g = sorted([r for r in rows if r["kc"] == kc], key=lambda r: r["eps"])
        xs = np.array([r["eps"] for r in g])
        ne = np.array([r["num_err_flip_ms"] for r in g], dtype=float)
        ax[1].plot(xs, ne, mk + "-", lw=1.5, ms=5, label=f"k_c={kc:.3g}")
    ax[1].set_xscale("log"); ax[1].set_yscale("log")
    ax[1].set_xlabel("ε / (m/s)"); ax[1].set_ylabel("双容差数值误差 / ms")
    ax[1].set_title("数值误差随 ε 的变化（可信下限）")
    ax[1].legend(fontsize=7); ax[1].grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "S2_ε扫描_诊断.png"), dpi=150)
    plt.close(fig)


# ============================== 汇总与判读 ==============================
def print_summary(rows, mu, outdir):
    print("=" * 76)
    print("S2：ε 灵敏度扫描")
    print(f"μ = {mu}   k_c ∈ {KC_SCAN} N/m   N0 = {N0} N, ω0 = {OMEGA0} rad/s")
    print(f"名义 ε = {EPS_NOM} m/s   板速 v = {V_BOARD} m/s")
    print("=" * 76)
    for kc in KC_SCAN:
        g = sorted([r for r in rows if r["kc"] == kc], key=lambda r: r["eps"])
        print(f"\n--- k_c = {kc:.3g} N/m ---")
        print(f"{'ε':>8} {'t*/ms':>9} {'tflip/ms':>9} {'flip':>5} {'I2nd/J':>10} "
              f"{'Δt_flip%':>9} {'数差t*/ms':>10} {'数差flip/ms':>11}")
        for r in g:
            def fmt(v, d=4):
                return f"{v:.{d}g}" if np.isfinite(v) else "—"
            rd = r["reldev_t_flip"] * 100 if np.isfinite(r["reldev_t_flip"]) else np.nan
            print(f"{r['eps']:>8.3g} {fmt(r['t_star_ms'],4):>9} "
                  f"{fmt(r['t_flip_ms'],4):>9} {str(r['flip']):>5} "
                  f"{fmt(r['I_2nd'],3):>10} {fmt(rd,3):>9} "
                  f"{fmt(r['num_err_t_ms'],2):>10} {fmt(r['num_err_flip_ms'],2):>11}")

    print("\n判读提示：")
    robust, fragile = [], []
    for kc in KC_SCAN:
        g = [r for r in rows if r["kc"] == kc]
        window = [r for r in g if 1e-3 <= r["eps"] <= 1e-1]
        devs = [abs(r["reldev_t_flip"]) for r in window
                if np.isfinite(r["reldev_t_flip"])]
        flips = sum(1 for r in window if r["flip"])
        if len(window) >= 3 and devs and max(devs) < 0.10 and flips == len(window):
            robust.append(kc)
        else:
            fragile.append(kc)
    if robust:
        print(f"  · k_c={['%.3g' % k for k in robust]}：ε∈[1e-3,1e-1] 内 t_flip 偏差 <10% "
              f"且翻转始终存在 → 该区间正则化稳健。")
    if fragile:
        print(f"  ⚠ k_c={['%.3g' % k for k in fragile]}：t_flip 漂移超 10% 或翻转丢失 "
              f"→ 预言③对这些 k_c 不稳健，S3 须 LuGre 化。")
    ne_all = [r["num_err_flip_ms"] for r in rows if np.isfinite(r["num_err_flip_ms"])]
    ne_eps = [(r["eps"], r["num_err_flip_ms"]) for r in rows
              if np.isfinite(r["num_err_flip_ms"])]
    if ne_all:
        worst = max(ne_eps, key=lambda x: x[1])
        print(f"  · 双容差最大数值误差 {worst[1]:.2e} ms @ ε={worst[0]:g}；"
              f"ε 更小时若误差持续上升，该 ε 以下不作物理判读。")
    print(f"\n输出目录：{os.path.abspath(outdir)}")


# ============================== 主流程 ==============================
def scan(kc_list, mu, outdir, want_plots=True):
    os.makedirs(outdir, exist_ok=True)
    global KC_SCAN
    KC_SCAN = kc_list
    rows, profiles = [], {}
    for kc in kc_list:
        for eps in EPS_LIST:
            keep = (abs(kc - KC_PROFILE) < 1e-9) and (eps in EPS_PROFILE)
            r, traj = run_dual(kc, N0, OMEGA0, mu, eps)
            if keep and traj is not None:
                profiles[eps] = traj
            rows.append(r)
        print(f"  [S2] k_c = {kc:.3g} N/m 完成（{len(EPS_LIST)} 个 ε ×2 容差）")
    add_diagnostics(rows)
    save_csv(rows, os.path.join(outdir, "S2_ε扫描_结果表.csv"))
    if want_plots:
        plot_convergence(rows, outdir)
        plot_profiles(profiles, outdir)
        plot_diagnostics(rows, outdir)
    print_summary(rows, mu, outdir)
    return rows


def main():
    ap = argparse.ArgumentParser(description="S2：ε 灵敏度扫描")
    ap.add_argument("--kc", type=float, nargs="+", default=KC_SCAN,
                    help="k_c 扫描值（N/m），须在翻转存在区间内")
    ap.add_argument("--mu", type=float, default=MU_DEFAULT)
    ap.add_argument("--outdir", default="./S2_ε扫描输出")
    ap.add_argument("--no-plots", action="store_true")
    args = ap.parse_args()
    scan(args.kc, args.mu, args.outdir, want_plots=not args.no_plots)


if __name__ == "__main__":
    main()
