# -*- coding: utf-8 -*-
"""
仿真 S4v1：v1 模型 (μ, h/L) × k_c 参数平面扫描
================================================
性质：**v1 退化模型**（固定铰轴、单自由度 θ、尖端法向罚弹簧、tanh 摩擦、无黏着）
的补充扫描。v2（弹性杆 + 手端正交弹簧，见 `code/v2_model.py`）已接管主线；
本脚本归档的是 2026-10-24 版交接说明 §7.2 留下的 v1 支线任务，作为 v1 基准的
最后一块，并与 v2 的 S6b 角度—摩擦平面互为独立对照。

对应任务：交接说明（v1 版）§7.2「数值支线：S4 参数平面扩展」——
扫描摩擦系数与几何臂比的二维平面，对每个点做 k_c 刚性极限检验，
寻找 t* 的振荡、跳变或失稳边界。

结论摘要（详见 simulations/S4v1_μh参数平面扫描/结果说明.md）
--------------------------------------------------------------
扫描箱内 49/49 组 t*(k_c) 严格单调，无振荡、无未分离、无跳变；
R(θ*) = 1 没有留下任何动力学签名。原因是 v1 模型的重加载项 ∝ N(μh − L sinθ)
在 v_rel 换向时 N 已趋于 0，且 θ̇ ≥ 0 使分离成为几何必然——这是 v1 的
结构性质，不能据此排除真实系统（或 v2）中的 Painlevé / reverse chatter。

为什么是这两个参数
------------------
分段内的力矩平衡为
    τ = (mgL/2) sinθ − N (L sinθ + μ h tanh(v_rel/ε))，  θ̈ = −τ/I_O
摩擦项的关键组合是 μh 与重力臂 L sinθ 之比
    R(θ) = μh / (L sinθ)。
θ 从 θ0 增大到 θ* 时 sinθ 减小，R 单调增大。名义点 R(θ*) = 0.943 < 1；
一旦 R 越过 1，v_rel 换向（tanh → −1）后法向项 (L sinθ − μh) 变号，
接触会被"越压越紧"而不是卸载，这正是 S1 尾部未覆盖的区域。

扫描设置
--------
μ    ∈ {0.20, 0.35, 0.50, 0.65, 0.80, 1.00, 1.20}
h/L  ∈ {0.70, 0.78, 0.85, 0.90, 0.9375, 0.96, 0.98}   （0.9375 = 名义 0.075/0.08）
k_c  ∈ {5e3, 1e4, 2e4, 5e4, 1e5} N/m
N0 = 4 N, ω0 = 0 固定。

每个点用双容差（rtol = 1e-9 / 1e-12）复算，|Δt*| 作为数值误差；
误差 > t* 的 1% 记为该点数值不可信（沿用 S1 判据）。

输出
----
S4v1_参数平面_结果表.csv
S4v1_参数平面_tstar.png      固定 k_c 的 log t* 平面 + R(θ*)=1 等值线
S4v1_参数平面_标度.png       d log t* / d log k_c 的平面（刚性极限是否规则）
S4v1_参数平面_区间图.png     分区间图：常规 / 翻转消失 / 数值不可信 / 未分离
S4v1_tstar_vs_mu.png         t* 随 μ 的切片
S4v1_tstar_vs_kc.png         t* 随 k_c 的对数折线（检验标度是否规则）
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s4v1_kernel as cm
import s4v1_plot as pp

MU_LIST = [0.20, 0.35, 0.50, 0.65, 0.80, 1.00, 1.20]
HL_LIST = [0.70, 0.78, 0.85, 0.90, 0.9375, 0.96, 0.98]
KC_LIST = [5e3, 1e4, 2e4, 5e4, 1e5]
N0 = 4.0
OM0 = 0.0
KC_NOM = 2e4
FIELDS = ["mu", "h_over_L", "theta_star_deg", "R_star", "kc", "N0", "omega0",
          "theta0_deg", "delta0_mm", "t_star_ms", "t_flip_ms", "flip",
          "separated", "I_2nd", "ell_mm", "omega_star", "xdot_max",
          "num_err_t_ms", "regime", "slope_t"]
PALETTE = [(31, 119, 180), (255, 127, 14), (44, 160, 44), (214, 39, 40),
           (148, 103, 189), (140, 86, 75), (227, 119, 194)]


def Rs(mu, hl):
    """分离角处的力臂比 R(θ*) = μh/(L sinθ*)。"""
    th_star = np.arccos(hl)
    return mu * hl / np.sin(th_star)


def run_point(mu, hl, kc):
    h = hl * cm.L
    a, err_t = cm.run_dual(kc, N0, OM0, mu=mu, h=h)
    rel = (err_t / (a["t_star"] * 1e3)
           if np.isfinite(a["t_star"]) and a["t_star"] > 0 else np.nan)
    if not a["separated"]:
        regime = "未分离"
    elif np.isfinite(rel) and rel > 0.01:
        regime = "数值不可信"
    elif not a["flip"]:
        regime = "翻转消失"
    else:
        regime = "常规"
    th_star = np.arccos(hl)
    return dict(
        mu=mu, h_over_L=hl, theta_star_deg=float(np.degrees(th_star)),
        R_star=float(Rs(mu, hl)), kc=kc, N0=N0, omega0=OM0,
        theta0_deg=float(np.degrees(a["theta0"])),
        delta0_mm=float(a["N0"] / kc * 1e3) if a["valid"] else np.nan,
        t_star_ms=float(a["t_star"] * 1e3) if a["valid"] else np.nan,
        t_flip_ms=float(a["t_flip"] * 1e3) if np.isfinite(a["t_flip"]) else np.nan,
        flip=bool(a["flip"]), separated=bool(a["separated"]),
        I_2nd=float(a["I2"]), ell_mm=float(a["ell"] * 1e3),
        omega_star=float(a["omega_star"]), xdot_max=float(a["xdot_max"]),
        num_err_t_ms=float(err_t) if np.isfinite(err_t) else np.nan,
        regime=regime, slope_t=np.nan)


def add_slopes(rows):
    groups = {}
    for r in rows:
        groups.setdefault((r["mu"], r["h_over_L"]), []).append(r)
    for g in groups.values():
        g.sort(key=lambda r: r["kc"])
        ks = np.array([r["kc"] for r in g], float)
        ts = np.array([r["t_star_ms"] for r in g], float)
        ok = np.isfinite(ts) & (ts > 0)
        slope = np.polyfit(np.log10(ks[ok]), np.log10(ts[ok]), 1)[0] if ok.sum() >= 2 else np.nan
        for r in g:
            r["slope_t"] = float(slope)


def save_csv(rows, path):
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if isinstance(v, float) and not np.isfinite(v) else v)
                        for k, v in r.items() if k in FIELDS})


def grid_of(rows, key, kc):
    Z = np.full((len(MU_LIST), len(HL_LIST)), np.nan)
    for r in rows:
        if abs(r["kc"] - kc) < 1e-9:
            Z[MU_LIST.index(r["mu"]), HL_LIST.index(r["h_over_L"])] = r[key]
    return Z


def _R1_overlay():
    """R(θ*)=1 在索引坐标中的折线。"""
    hls = np.linspace(HL_LIST[0], HL_LIST[-1], 200)
    mus = np.sqrt(1 - hls ** 2) / hls
    xs = np.interp(hls, HL_LIST, np.arange(len(HL_LIST)))
    ys = np.interp(mus, MU_LIST, np.arange(len(MU_LIST)))
    keep = (mus >= min(MU_LIST) * 0.9) & (mus <= max(MU_LIST) * 1.1)
    return xs[keep], ys[keep]


def plot_tstar(rows, outdir, kc=KC_NOM):
    Z = grid_of(rows, "t_star_ms", kc)
    logZ = np.log10(np.where(Z > 0, Z, np.nan))
    xs, ys = _R1_overlay()
    pp.heatmap(os.path.join(outdir, "S4v1_参数平面_tstar.png"), logZ, MU_LIST, HL_LIST,
               title=f"log10(t* / ms)，k_c = {kc:.0e} N/m，N0 = 4 N",
               cbar_label="log10(t*/ms)", cmap_name="viridis",
               overlay=(xs, ys, (230, 40, 40), "R(θ*)=1"))


def plot_slope(rows, outdir):
    Z = grid_of(rows, "slope_t", KC_NOM)
    pp.heatmap(os.path.join(outdir, "S4v1_参数平面_标度.png"), Z, MU_LIST, HL_LIST,
               title="d log t* / d log k_c（−0.5 = 线性弹簧基准，−1 = 完全刚性标度）",
               cbar_label="斜率", cmap_name="viridis", vmin=-1.0, vmax=-0.4)


def plot_regime(rows, outdir):
    """缺翻转的 k_c 个数（0–5）。未分离与数值不可信的点数在所有格内均为 0。"""
    Z = np.zeros((len(MU_LIST), len(HL_LIST)))
    for i, mu in enumerate(MU_LIST):
        for j, hl in enumerate(HL_LIST):
            g = [r for r in rows if r["mu"] == mu and r["h_over_L"] == hl]
            Z[i, j] = sum(1 for r in g if not r["flip"])
    ann = [[f"{int(Z[i, j])}" for j in range(len(HL_LIST))] for i in range(len(MU_LIST))]
    pp.heatmap(os.path.join(outdir, "S4v1_参数平面_区间图.png"), Z, MU_LIST, HL_LIST,
               title="5 个 k_c 中未出现摩擦翻转的个数（未分离=0，数值不可信=0）",
               cbar_label="个数", cmap_name="viridis", vmin=0, vmax=5,
               annotation=ann)


def plot_vs_mu(rows, outdir, kc=KC_NOM):
    series = []
    for j, hl in enumerate(HL_LIST):
        g = sorted([r for r in rows if abs(r["kc"] - kc) < 1e-9 and r["h_over_L"] == hl],
                   key=lambda r: r["mu"])
        if not g:
            continue
        series.append(dict(x=[r["mu"] for r in g], y=[r["t_star_ms"] for r in g],
                           label=f"h/L={hl:g}", color=PALETTE[j % len(PALETTE)]))
    pp.lineplot(os.path.join(outdir, "S4v1_tstar_vs_mu.png"), series,
                title=f"t* 随 μ 的变化（k_c = {kc:.0e} N/m，N0 = 4 N）",
                xlabel="μ", ylabel="t* / ms", ylog=True)


def plot_vs_kc(rows, outdir):
    cases = [(0.20, 0.70), (0.35, 0.9375), (0.65, 0.90), (0.80, 0.90),
             (1.00, 0.96), (1.20, 0.98)]
    series = []
    for i, (mu, hl) in enumerate(cases):
        g = sorted([r for r in rows if r["mu"] == mu and r["h_over_L"] == hl],
                   key=lambda r: r["kc"])
        if not g:
            continue
        series.append(dict(x=[r["kc"] for r in g], y=[r["t_star_ms"] for r in g],
                           label=f"μ={mu:g},h/L={hl:g}", color=PALETTE[i % len(PALETTE)]))
    pp.lineplot(os.path.join(outdir, "S4v1_tstar_vs_kc.png"), series,
                title="t* 随 k_c 的刚性极限检验（双对数）",
                xlabel="k_c / (N/m)", ylabel="t* / ms", xlog=True, ylog=True)


def main():
    ap = argparse.ArgumentParser(description="S4：(μ, h/L) × k_c 参数平面扫描")
    ap.add_argument("--outdir", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "simulations", "S4v1_μh参数平面扫描"))
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    if args.quick:
        print(run_point(0.35, 0.9375, KC_NOM))
        return

    rows = []
    total = len(MU_LIST) * len(HL_LIST) * len(KC_LIST)
    n = 0
    for mu in MU_LIST:
        for hl in HL_LIST:
            for kc in KC_LIST:
                rows.append(run_point(mu, hl, kc))
                n += 1
            print(f"  [{n}/{total}] μ={mu:g} h/L={hl:g} 完成", flush=True)
    add_slopes(rows)
    save_csv(rows, os.path.join(outdir, "S4v1_参数平面_结果表.csv"))
    plot_tstar(rows, outdir)
    plot_slope(rows, outdir)
    plot_regime(rows, outdir)
    plot_vs_mu(rows, outdir)
    plot_vs_kc(rows, outdir)

    print("=" * 78)
    print(f"S4 参数平面完成：{total} 点（每点双容差复算）")
    print(f"{'μ':>5} {'h/L':>7} {'R*':>7} {'t*@2e4/ms':>11} {'斜率':>8}  区间")
    for mu in MU_LIST:
        for hl in HL_LIST:
            g = [r for r in rows if r["mu"] == mu and r["h_over_L"] == hl]
            nom = [r for r in g if abs(r["kc"] - KC_NOM) < 1e-9][0]
            regs = ",".join(sorted({r["regime"] for r in g}))
            print(f"{mu:>5g} {hl:>7g} {g[0]['R_star']:>7.3f} {nom['t_star_ms']:>11.4f} "
                  f"{g[0]['slope_t']:>+8.3f}  {regs}")
    print(f"\n输出目录：{outdir}")


if __name__ == "__main__":
    main()
