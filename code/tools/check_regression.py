# -*- coding: utf-8 -*-
"""回归测试：确认 v1 基准、v2 基线与 I(ξ) 后处理没有被静默改坏。

用法（在本机 Anaconda Python 下）：
    python code/tools/check_regression.py
通过时退出码 0；任一断言失败退出码 1。

基准值来源：
  · v1 名义工况：t* = 1.292 ms, I_2nd = 1.75e-4 J（模型报告/本次复算）
  · v2 基线 A：N 平台、I_c/I_plateau > 0.8、∫I dξ/I_2nd ≈ 1（S4 结果说明）
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


# ---------------- 3. v2 基线回归 ----------------
def test_v2_baseline():
    p = dict(vm.DEFAULT, friction="lugre_springslider", sigma1=0.0, sigma2=0.0)
    info, tr = vm.simulate(p, t_max=0.004)
    xc, I, met = vm.dash_profile(tr, nbins=120, sigma_bins=1.0)
    check("v2 基线：N 平台（相对涨落 <10%）", info["N_rel_std"] < 0.10,
          f"Nrel={info['N_rel_std']*100:.1f}%")
    check("v2 基线：有限强度截断（I_c/I_p > 0.8）",
          met["cutoff_ratio"] > 0.8, f"Ic/Ip={met['cutoff_ratio']:.3f}")
    check("v2 基线：强度守恒", abs(met["I_conservation"] - 1.0) < 0.01,
          f"∫/I2={met['I_conservation']:.4f}")
    check("v2 基线：划痕长度量级 0.2–0.8 mm",
          0.2 < info["xi_end_mm"] < 0.8, f"ξ={info['xi_end_mm']:.3f} mm")


if __name__ == "__main__":
    print("=" * 68)
    print("回归测试：v1 基准 / v2 基线 / I(ξ) 后处理")
    print("=" * 68)
    test_deposit_profile()
    test_v1_nominal()
    test_v2_baseline()
    print("-" * 68)
    if FAIL:
        print(f"结果：{len(FAIL)} 项失败 -> {FAIL}")
        sys.exit(1)
    print("结果：全部通过")
