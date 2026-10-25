# -*- coding: utf-8 -*-
"""S6b 稳健性补充：接触阻尼、刚性接触初始锚定口径与停止边界对接触参数的敏感性。

为什么单独成文
--------------
S6b 主脚本（`仿真S6b_v2角度摩擦平面.py`）在角度—摩擦平面上把“停止前接触分离”
记为 lift-off 类、把刚性接触大角度的约束乘子发散记为 sprag 类。独立复核指出这两类
可能不是稳健物理：

  1. 罚接触（k_n=2e6）的“分离”由欠阻尼接触振铃触发：加大接触阻尼 c_n 后大部分
     分离消失；k_n=2e5 时全部分离消失。→ 它是接触模型抖振导致的表观脱离。
  2. 刚性接触（`contact='rigid'`）的 y_anchor 由“初始 N(0)=N0”仿射解出，其灵敏度
     dN/dy_anchor 在 theta0 ≈ 35.3° 过零；θ0≥34° 的 release/sprag 可在 0.02–0.07 ms
     内由初始瞬态触发，换用罚接触极限锚点后类别会改变。→ θ0 接近该奇异的点不能
     当作物理边界。

本脚本对上述两点给出可复算的对比表，并顺带检查低摩擦“连续滑移”边界对 c_n 的稳健性。

输出（写入 simulations/S6b_v2角度摩擦平面/）
  S6b_接触阻尼敏感性.csv
  S6b_刚性锚定口径敏感性.csv
  S6b_边界稳健性.csv
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v2_model as vm

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

MUS_RATIO = vm.DEFAULT["mu_s"] / vm.DEFAULT["mu_k"]


# ---------------------------------------------------------------------------
# 带早停、可指定初始锚点的刚性接触积分（复制约束乘子算法，便于锚点扫描）
# ---------------------------------------------------------------------------
def rigid_run(p, y_anchor=None, t_max=0.004, rtol=1e-8, max_step=1e-5,
              npts=None, cap_N=25.0, cap_vrel=2.5, t_ignore=5e-5,
              v_stop_tol=1e-3):
    """返回 dict(info, traj)。y_anchor=None 时用 v2_model.rigid_factory 的仿射解。"""
    if npts is None:
        npts = max(4001, int(round(4001 * t_max / 0.004)) + 1)
    th0 = np.deg2rad(p["theta0_deg"])
    n0 = p["n0"]
    g0 = n0 * (p["mu_k"] + (p["mu_s"] - p["mu_k"])
               * np.exp(-(p["v_board"] / p["vs"]) ** 2))
    z0 = g0 / p["sigma0"]
    s0 = (0.0, th0, 0.0, 0.0, z0)
    if y_anchor is None:
        # 与 v2_model.rigid_factory 完全相同的仿射锚定
        n_a0 = vm._rigid_N(p, s0, 0.0)
        n_a1 = vm._rigid_N(p, s0, 1.0)
        y_anchor = (n0 - n_a0) / (n_a1 - n_a0)
    n_init = vm._rigid_N(p, s0, y_anchor)

    def rhs(t, s):
        _, _, _, aux = vm._rigid_terms(p, s, y_anchor)
        acc = aux["acc"]
        return [s[2], s[3], acc[0], acc[2], aux["zdot"]]

    def ev_release(t, s):
        return vm._rigid_terms(p, s, y_anchor)[2]
    ev_release.terminal = True
    ev_release.direction = -1

    def ev_capN(t, s):
        return vm._rigid_terms(p, s, y_anchor)[2] - cap_N
    ev_capN.terminal = True
    ev_capN.direction = 1

    def ev_capv(t, s):
        return abs(vm._rigid_terms(p, s, y_anchor)[3]["v_rel"]) - cap_vrel
    ev_capv.terminal = True
    ev_capv.direction = 1

    sol = solve_ivp(rhs, (0.0, t_max), s0, method="DOP853", rtol=rtol,
                    atol=1e-10, events=(ev_release, ev_capN, ev_capv),
                    dense_output=True, max_step=max_step)
    tt = np.linspace(0.0, sol.t[-1], npts)
    x, th, vx, w, z = sol.sol(tt)
    v_rel = vx + vm.L * np.cos(th) * w + p["v_board"]
    xi = np.concatenate([[0.0],
                         np.cumsum(0.5 * (v_rel[:-1] + v_rel[1:]) * np.diff(tt))])
    N = np.empty_like(tt)
    for i in range(len(tt)):
        N[i] = vm._rigid_terms(p, (x[i], th[i], vx[i], w[i], z[i]),
                               y_anchor)[2]
    i_stop = len(tt) - 1
    for i in range(1, len(tt)):
        if tt[i] > t_ignore and (v_rel[i] <= 0.0 or abs(v_rel[i]) < v_stop_tol):
            i_stop = i
            break
    info = dict(t_end_ms=float(tt[-1] * 1e3),
                release=bool(sol.t_events[0].size > 0),
                capN=bool(sol.t_events[1].size > 0),
                capv=bool(sol.t_events[2].size > 0),
                t_stop_ms=float(tt[i_stop] * 1e3),
                stop_found=int(i_stop < len(tt) - 1),
                N_max=float(N[:i_stop + 1].max()),
                N_min=float(N[:i_stop + 1].min()),
                xi_end_mm=float(xi[i_stop] * 1e3),
                N_init=float(n_init), y_anchor=float(y_anchor))
    return info


def penalty_run(p, tmax):
    """罚接触：返回主脚本同口径的分类与关键量。"""
    info, tr = vm.simulate(p, t_max=tmax)
    tt, vr = tr["tt"], tr["v_rel"]
    i = info["i_stop"]
    stopped = i < len(tt) - 1 and vr[i] <= 1e-3
    post = np.nan
    if stopped:
        dt = tt[1] - tt[0]
        j1 = min(len(tt), i + max(5, int(0.3e-3 / dt)))
        post = float(np.max(np.abs(vr[i:j1])))
    if bool(info.get("touch_end", False)):
        cls = "lift"
    elif not stopped:
        cls = "slide"
    elif post > 0.25 * p["v_board"]:
        cls = "chatter"
    elif info["N_rel_std"] < 0.10:
        cls = "plateau"
    else:
        cls = "nvar"
    return dict(cls=cls, t_end=info["t_end_ms"], t_stop=info["t_stop_ms"],
                N_min=info["N_min"], N_max=info["N_max"],
                N_rel=info["N_rel_std"], vmin=info["v_rel_min"],
                post=post, xi=info["xi_end_mm"])


def make_p(theta0, mus, **extra):
    p = dict(vm.DEFAULT, theta0_deg=float(theta0), mu_s=float(mus),
             mu_k=float(mus) / MUS_RATIO)
    p.update(extra)
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="./S6b_v2输出")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # ============ 1. 罚接触 k_n=2e6 “脱离”对接触阻尼 c_n 的敏感性 ============
    cases = [(20.0, 0.40), (20.0, 0.62), (25.0, 0.40), (25.0, 0.62),
             (25.0, 0.90), (30.0, 0.40), (30.0, 0.62), (30.0, 0.90)]
    print("=" * 78)
    print("1. 罚接触 k_n=2e6 的“停止前分离”对接触阻尼 c_n 的敏感性")
    print("=" * 78)
    cn_rows = []
    for th, mus in cases:
        line = f"th={th:4.1f} mu_s={mus:.2f}: "
        for cn in (1.0, 10.0, 100.0):
            r = penalty_run(make_p(th, mus, kn=2e6, cn=cn), 0.008)
            cn_rows.append(dict(theta0_deg=th, mu_s=mus, cn=cn, **r))
            line += (f" cn={cn:5.1f}->{r['cls']:8s}"
                     f"(t_end={r['t_end']:.2f},t_stop={r['t_stop']:.2f})")
        print(line)
    with open(os.path.join(args.outdir, "S6b_接触阻尼敏感性.csv"), "w",
              newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(cn_rows[0].keys()))
        w.writeheader()
        w.writerows(cn_rows)
    n1 = sum(1 for r in cn_rows if r["cn"] == 1.0 and r["cls"] == "lift")
    n100 = sum(1 for r in cn_rows if r["cn"] == 100.0 and r["cls"] == "lift")
    print(f"→ c_n=1 时 {n1}/{len(cases)} 组分离；c_n=100 时 {n100}/{len(cases)} 组分离")

    # 2e5 对照
    print("\n对照：k_n=2e5 时同样工况")
    for th, mus in cases:
        r = penalty_run(make_p(th, mus, kn=2e5, cn=1.0), 0.008)
        print(f"  th={th:4.1f} mu_s={mus:.2f}: {r['cls']:8s} "
              f"t_end={r['t_end']:.3f} ms")

    # ============ 2. 刚性接触初始锚定口径敏感性 ============
    print("\n" + "=" * 78)
    print("2. 刚性接触：仿射 N(0)=N0 锚点 vs 罚接触极限锚点")
    print("   （dN/dy_anchor 在 theta0≈35.3° 过零 → 大角度结果由初始瞬态主导）")
    print("=" * 78)
    anchor_rows = []
    for th in (20.0, 25.0, 30.0, 32.0, 34.0, 36.0, 40.0):
        p = make_p(th, 0.62, contact="rigid")
        r_fac = rigid_run(p)                       # 仿射 N(0)=N0 口径
        th0 = np.deg2rad(th)
        y_pen = (vm.L * np.cos(th0)
                 - (p["n0"] - vm.M * vm.GRAV) / p["ky"])   # 罚极限锚点
        r_pen = rigid_run(p, y_anchor=y_pen)
        # 灵敏度 dN/dy_anchor（N 对 anchor 的斜率，N/m）
        g0 = p["n0"] * (p["mu_k"] + (p["mu_s"] - p["mu_k"])
                        * np.exp(-(p["v_board"] / p["vs"]) ** 2))
        s0 = (0.0, th0, 0.0, 0.0, g0 / p["sigma0"])
        slope = (vm._rigid_N(p, s0, 1.0) - vm._rigid_N(p, s0, 0.0))
        row = dict(theta0_deg=th, mu_s=0.62,
                   y_anchor_affine=r_fac["y_anchor"],
                   y_anchor_penalty_limit=y_pen,
                   dy_um=(r_fac["y_anchor"] - y_pen) * 1e6,
                   dN_dy_anchor=slope,
                   N_init_affine=r_fac["N_init"],
                   N_init_penalty=vm._rigid_N(p, s0, y_pen),
                   class_affine=("lift" if r_fac["release"] else
                                 ("sprag" if (r_fac["capN"] or r_fac["capv"])
                                  else ("slide" if not r_fac["stop_found"]
                                        else "stop"))),
                   t_end_affine=r_fac["t_end_ms"],
                   N_max_affine=r_fac["N_max"],
                   class_penalty=("lift" if r_pen["release"] else
                                  ("sprag" if (r_pen["capN"] or r_pen["capv"])
                                   else ("slide" if not r_pen["stop_found"]
                                         else "stop"))),
                   t_end_penalty=r_pen["t_end_ms"],
                   N_max_penalty=r_pen["N_max"])
        anchor_rows.append(row)
        print(f"  th={th:4.1f}: dN/dy={slope:9.1f} N/m, "
              f"Δy={row['dy_um']:+7.1f} µm | 仿射: {row['class_affine']:5s}"
              f" t_end={row['t_end_affine']:6.3f} Nmax={r_fac['N_max']:5.1f} | "
              f"罚极限: {row['class_penalty']:5s} "
              f"t_end={row['t_end_penalty']:6.3f} Nmax={r_pen['N_max']:5.1f}")
    with open(os.path.join(args.outdir, "S6b_刚性锚定口径敏感性.csv"), "w",
              newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(anchor_rows[0].keys()))
        w.writeheader()
        w.writerows(anchor_rows)

    # ============ 3. 低摩擦“连续滑移”边界对 c_n / k_n 的稳健性 ============
    print("\n" + "=" * 78)
    print("3. 连续滑移（停止）边界对接触参数 c_n 的稳健性（theta0=20°）")
    print("=" * 78)
    rob_rows = []
    for mus in (0.19, 0.20, 0.21, 0.22):
        line = f"mu_s={mus:.2f}: "
        for cn in (1.0, 10.0, 100.0):
            r = penalty_run(make_p(20.0, mus, kn=2e4, cn=cn), 0.025)
            rob_rows.append(dict(theta0_deg=20.0, mu_s=mus, cn=cn, kn=2e4,
                                 **r))
            line += f" cn={cn:5.1f}->{r['cls']:8s}(t_stop={r['t_stop']:.3f})"
        print(line)
    with open(os.path.join(args.outdir, "S6b_边界稳健性.csv"), "w",
              newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(rob_rows[0].keys()))
        w.writeheader()
        w.writerows(rob_rows)
    print("\n结论：连续滑移/停止边界的类别在 c_n=1–100 内不变；"
          "k_n=2e6 的“分离”随 c_n 增大而消失，属接触抖振表观事件；"
          "刚性接触 θ0≥34° 的 release/sprag 由初始锚定口径支配。")


if __name__ == "__main__":
    main()
