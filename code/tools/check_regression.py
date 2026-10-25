# -*- coding: utf-8 -*-
"""回归测试：确认 v1 基准、v2 基线与 I(ξ) 后处理没有被静默改坏。

用法（在本机 Anaconda Python 下）：
    python code/tools/check_regression.py
通过时退出码 0；任一断言失败退出码 1。

基准值来源：
  · v1 名义工况：t* = 1.292 ms, I_2nd = 1.75e-4 J（模型报告/本次复算）
  · v2 基线 A（2026-10-25 A2 预载修正后）：N 涨落 0.25%、I_c/I_plateau ≈ 1.002、
    ∫I dξ/I_2nd ≈ 1、ξ_end ≈ 0.442 mm
  · v2 刚性接触（LCP 极限）：t_stop ≈ 1.550 ms、约束残差 <1e-9、N(0)=N0
  · v2 S6b（2026-10-25）：单事件停止边界 μ_s≈0.205（θ0=20°, k_n=2e4, v=0.5）；
    近刚性罚接触（k_n=2e6）在 θ0=30°、μ_s=0.50 处停止前脱离，而 θ0=20°、μ_s=0.62
    仍为“N变化+截断”；刚性接触 θ0=40° 立即脱离、θ0=30° 停止后强振铃。
"""
from __future__ import annotations

import os
import sys

import numpy as np
from scipy.integrate import solve_ivp

CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, CODE)
import v2_model as vm                      # noqa: E402
from deposit_profile import (deposit_profile, profile_integral,  # noqa: E402
                             total_deposit, trace_coordinate)

FAIL = []
MUS_RATIO = vm.DEFAULT["mu_s"] / vm.DEFAULT["mu_k"]


def check(name, cond, detail=""):
    ok = bool(cond)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAIL.append(name)


# ---------------- 1. I(ξ) 后处理自检 ----------------
def test_deposit_profile():
    t = np.linspace(0.0, 1e-3, 4001)
    v = 0.5
    xi = trace_coordinate(t, np.full_like(t, v))
    q = np.full_like(t, 0.7)
    xc, I = deposit_profile(t, q, xi=xi, nbins=200, sigma_bins=0.0)
    check("deposit_profile 常速强度 I=q/v", abs(np.median(I) - q[0] / v) < 1e-6,
          f"I_med={np.median(I):.6f}")
    err = abs(profile_integral(xc, I) - total_deposit(t, q)) / total_deposit(t, q)
    check("deposit_profile 守恒 ∫I dξ=I2nd", err < 1e-6, f"rel_err={err:.2e}")


# ---------------- 2. v1 名义工况回归 ----------------
def test_v1_nominal():
    L, h = 0.08, 0.075
    rr, rho, g = 3.5e-3, 2300.0, 9.8
    m = rho * np.pi * rr ** 2 * L
    I_O = m * L ** 2 / 3.0
    mu, kc, v, eps = 0.35, 2.0e4, 0.5, 0.02

    def rhs(t, y):
        th, w = y
        d = L * np.cos(th) - h
        if d <= 0.0:
            return [w, -(m * g * L / 2.0) * np.sin(th) / I_O]
        N = kc * d
        v_rel = v - L * np.cos(th) * w
        tau = (m * g * L / 2.0) * np.sin(th) - N * (
            L * np.sin(th) + mu * h * np.tanh(v_rel / eps))
        return [w, -tau / I_O]

    def ev(t, y):
        return L * np.cos(y[0]) - h
    ev.terminal, ev.direction = True, -1
    th0 = np.arccos((h + 4.0 / kc) / L)
    sol = solve_ivp(rhs, (0, 0.05), [th0, 0.0], method="DOP853",
                    rtol=1e-10, atol=1e-12, events=ev, dense_output=True)
    t_star = sol.t_events[0][0]
    tt = np.linspace(0, t_star, 4001)
    th, w = sol.sol(tt)
    N = kc * np.clip(L * np.cos(th) - h, 0, None)
    I2 = mu * np.trapezoid(N * np.abs(v - L * np.cos(th) * w), tt)
    check("v1 t* 回归 1.292 ms", abs(t_star * 1e3 - 1.292) < 0.01,
          f"{t_star*1e3:.4f} ms")
    check("v1 I2nd 回归 1.75e-4 J", abs(I2 / 1.7497e-4 - 1.0) < 0.01,
          f"{I2:.4e} J")


# ---------------- 3. v2 基线回归（2026-10-25 A2 预载修正） ----------------
def test_v2_baseline():
    p = dict(vm.DEFAULT, friction="lugre_springslider", sigma1=0.0, sigma2=0.0)
    info, tr = vm.simulate(p, t_max=0.004)
    xc, I, met = vm.dash_profile(tr, nbins=120, sigma_bins=1.0)
    check("v2 基线：N 平台（相对涨落 <2%）", info["N_rel_std"] < 0.02,
          f"Nrel={info['N_rel_std']*100:.2f}%")
    check("v2 基线：有限强度截断（I_c/I_p > 0.95）",
          met["cutoff_ratio"] > 0.95, f"Ic/Ip={met['cutoff_ratio']:.3f}")
    check("v2 基线：强度守恒", abs(met["I_conservation"] - 1.0) < 1e-3,
          f"∫/I2={met['I_conservation']:.5f}")
    check("v2 基线：划痕长度 0.40–0.50 mm",
          0.40 < info["xi_end_mm"] < 0.50, f"ξ={info['xi_end_mm']:.3f} mm")
    check("v2 基线：t_stop ≈ 1.816 ms",
          abs(info["t_stop_ms"] - 1.816) < 0.005,
          f"t_stop={info['t_stop_ms']:.4f} ms")


# ---------------- 4. v2 刚性接触（LCP 极限）回归 ----------------
def test_v2_rigid():
    p = dict(vm.DEFAULT, contact="rigid")
    info, tr = vm.simulate(p, t_max=0.004)
    xc, I, met = vm.dash_profile(tr, nbins=120, sigma_bins=1.0)
    check("刚性接触：初始 N(0)=N0", abs(info["n0_check"] - p["n0"]) < 1e-6,
          f"N(0)={info['n0_check']:.6f} N")
    check("刚性接触：约束加速度残差 <1e-9",
          info["constraint_residual_max"] < 1e-9,
          f"max|Cdd|={info['constraint_residual_max']:.1e}")
    check("刚性接触：无释放事件（基线）", info["release_events"] == 0,
          f"events={info['release_events']}")
    check("刚性接触：N ≥ 0", info["N_min"] >= 0.0,
          f"Nmin={info['N_min']:.3f} N")
    check("刚性接触：t_stop ≈ 1.550 ms",
          abs(info["t_stop_ms"] - 1.550) < 0.005,
          f"t_stop={info['t_stop_ms']:.4f} ms")
    check("刚性接触：平台+截断（I_c/I_p > 0.95）",
          met["cutoff_ratio"] > 0.95, f"Ic/Ip={met['cutoff_ratio']:.3f}")
    check("刚性接触：强度守恒", abs(met["I_conservation"] - 1.0) < 1e-3,
          f"∫/I2={met['I_conservation']:.5f}")


# ---------------- 5. v2 S6b 角度—摩擦平面回归 ----------------
def test_v2_s6b():
    # 5.1 单事件停止边界：μ_s=0.20 连续滑移，μ_s=0.22 在 25 ms 内停止
    p = dict(vm.DEFAULT, theta0_deg=20.0, mu_s=0.20,
             mu_k=0.20 / MUS_RATIO)
    info, tr = vm.simulate(p, t_max=0.025)
    stopped = info["i_stop"] < len(tr["tt"]) - 1
    check("S6b 连续滑移边界：μ_s=0.20 在 25 ms 内不停止",
          (not stopped) and info["v_rel_min"] > 1e-3,
          f"t_stop={info['t_stop_ms']:.3f} ms, vmin={info['v_rel_min']:.4f}")
    p = dict(vm.DEFAULT, theta0_deg=20.0, mu_s=0.22,
             mu_k=0.22 / MUS_RATIO)
    info, tr = vm.simulate(p, t_max=0.025)
    stopped = info["i_stop"] < len(tr["tt"]) - 1
    check("S6b 停止边界：μ_s=0.22 在 25 ms 内停止（t≈6 ms）",
          stopped and 5.5 < info["t_stop_ms"] < 7.5,
          f"t_stop={info['t_stop_ms']:.3f} ms")

    # 5.2 近刚性罚接触（k_n=2e6）：θ0=30°, μ_s=0.50 停止前脱离
    p = dict(vm.DEFAULT, theta0_deg=30.0, mu_s=0.50,
             mu_k=0.50 / MUS_RATIO, kn=2e6)
    info, tr = vm.simulate(p, t_max=0.006)
    check("S6b 近刚性接触脱离：θ0=30°, μ_s=0.50",
          info["touch_end"] and info["t_end_ms"] < 2.0,
          f"t_end={info['t_end_ms']:.3f} ms, touch_end={info['touch_end']}")

    # 5.3 近刚性罚接触（k_n=2e6）：θ0=20°, μ_s=0.62 仍为“N变化+截断”
    p = dict(vm.DEFAULT, theta0_deg=20.0, mu_s=0.62,
             mu_k=0.62 / MUS_RATIO, kn=2e6)
    info, tr = vm.simulate(p, t_max=0.006)
    xc, I, met = vm.dash_profile(tr, nbins=120, sigma_bins=1.0)
    check("S6b 近刚性接触未脱离：θ0=20°, μ_s=0.62",
          (not info["touch_end"]) and 1.8 < info["t_stop_ms"] < 2.2
          and 0.12 < info["N_rel_std"] < 0.22,
          f"t_stop={info['t_stop_ms']:.3f} ms, Nrel={info['N_rel_std']*100:.1f}%")

    # 5.4 刚性接触：θ0=40° 在 0.1 ms 内接触释放
    p = dict(vm.DEFAULT, theta0_deg=40.0, contact="rigid")
    info, tr = vm.simulate(p, t_max=0.004)
    check("S6b 刚性接触大角度脱离：θ0=40°",
          info["touch_end"] and info["t_end_ms"] < 0.1,
          f"t_end={info['t_end_ms']:.4f} ms")

    # 5.5 刚性接触：θ0=30° 停止后强振铃（N 最大 >10 N、t_stop<1 ms）
    p = dict(vm.DEFAULT, theta0_deg=30.0, contact="rigid")
    info, tr = vm.simulate(p, t_max=0.004)
    vr = tr["v_rel"]
    i = info["i_stop"]
    dt = tr["tt"][1] - tr["tt"][0]
    j1 = min(len(vr), i + max(5, int(0.3e-3 / dt)))
    post = float(np.max(np.abs(vr[i:j1]))) if i < len(vr) - 1 else np.nan
    check("S6b 刚性接触停止后颤振：θ0=30°",
          info["N_max"] > 10.0 and info["t_stop_ms"] < 1.0
          and post > 0.1,
          f"Nmax={info['N_max']:.1f} N, t_stop={info['t_stop_ms']:.3f} ms, "
          f"post|v|={post:.3f} m/s")


if __name__ == "__main__":
    print("=" * 68)
    print("回归测试：v1 基准 / v2 基线 / v2 刚性接触 / v2 S6b 角度—摩擦 / I(ξ) 后处理")
    print("=" * 68)
    test_deposit_profile()
    test_v1_nominal()
    test_v2_baseline()
    test_v2_rigid()
    test_v2_s6b()
    print("-" * 68)
    if FAIL:
        print(f"结果：{len(FAIL)} 项失败 -> {FAIL}")
        sys.exit(1)
    print("结果：全部通过")
