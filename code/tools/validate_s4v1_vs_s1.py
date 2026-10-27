# -*- coding: utf-8 -*-
"""
S4v1 回归校验：用 S1 已发布 CSV 校验 s4v1_kernel 的积分结果
============================================================
本机 DSH 运行时没有 SciPy，S4v1 的积分器是纯 NumPy 重写（DOPRI5 + 步内二分事件），
与 S1 的 solve_ivp(DOP853) 是不同的实现。本脚本把两者在 S1 的 105 个可分离工况上
逐点比较，作为"口径一致"的证据，而不是一句口头声明。

用法：
    python code/tools/validate_s4v1_vs_s1.py
    python code/tools/validate_s4v1_vs_s1.py --csv simulations/S1_kc刚性极限收敛扫描/fixed_N0/S1_kc扫描_结果表.csv

判据：t* 相对偏差 < 1e-3 视为通过；当前实测最大值约 2e-6。
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import s4v1_kernel as cm

DEFAULT_CSV = os.path.join(HERE, "..", "..", "simulations",
                           "S1_kc刚性极限收敛扫描", "fixed_N0", "S1_kc扫描_结果表.csv")


def main():
    ap = argparse.ArgumentParser(description="S4v1 积分器对 S1 的回归校验")
    ap.add_argument("--csv", default=DEFAULT_CSV)
    ap.add_argument("--tol", type=float, default=1e-3, help="t* 相对偏差通过阈值")
    args = ap.parse_args()

    rows = []
    with open(args.csv, encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            if r["mode"] == "fixed_N0":
                rows.append(r)
    rows = [r for r in rows if r["t_star_ms"].strip() != ""]
    print(f"S1 参考行数（fixed_N0，可分离）= {len(rows)}\n")

    worst = []
    for r in rows:
        kc, N0, om = float(r["kc"]), float(r["N0"]), float(r["omega0"])
        ref_t = float(r["t_star_ms"])
        got = cm.integrate_segment(kc, N0, om, mu=0.35, rtol=1e-9)
        d_t = got["t_star"] * 1e3 - ref_t
        worst.append((abs(d_t) / max(ref_t, 1e-12), kc, N0, om, d_t,
                      got["I2"] - float(r["I_2nd"]),
                      got["ell"] * 1e3 - float(r["ell_mm"])))

    worst.sort(reverse=True)
    print(f"{'kc':>9} {'N0':>4} {'w0':>4} {'Δt*/ms':>11} {'ΔI2':>11} {'Δℓ/mm':>10}")
    for rel, kc, N0, om, d_t, d_I, d_e in worst[:8]:
        print(f"{kc:>9.3g} {N0:>4.0f} {om:>4.0f} {d_t:>11.2e} {d_I:>11.2e} {d_e:>10.2e}")
    print(f"\n最大相对 t* 偏差 = {worst[0][0]:.2e}")
    n_bad = sum(1 for w in worst if w[0] > args.tol)
    print(f"超过阈值 {args.tol:g} 的点数 = {n_bad} / {len(rows)} → {'PASS' if n_bad == 0 else 'FAIL'}")
    return 0 if n_bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
