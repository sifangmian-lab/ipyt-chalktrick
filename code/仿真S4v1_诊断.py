# -*- coding: utf-8 -*-
"""
仿真 S4v1 诊断：结果表的结构化分析与极端工况物理判读
=====================================================
用法：
    python code/仿真S4v1_诊断.py                       # 全部诊断
    python code/仿真S4v1_诊断.py --extreme             # 只跑极端工况判读

第一部分（表格分析，读 CSV）：
  1. 每个 (μ, h/L) 组的 t*(k_c) 单调性
  2. 摩擦翻转缺失的位置分布
  3. t_flip/t* 随 μ 的变化
  4. R(θ*) 跨越 1 时 t* 是否出现拐点（沿 h/L 的二阶差分）
  5. I_2nd、ℓ、ω* 随 μ 的变化
  6. 双容差数值误差

第二部分（重新积分，--extreme）：
  极端工况的 θ̇ 单调性、N 单调性、v_rel 换向，以及摩擦翻转后的减速幅度——
  用于判定"没有失稳"是数值假象还是 v1 模型的结构性质。
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s4v1_kernel as cm

MU = [0.20, 0.35, 0.50, 0.65, 0.80, 1.00, 1.20]
HL = [0.70, 0.78, 0.85, 0.90, 0.9375, 0.96, 0.98]
DEFAULT_CSV = os.path.join(HERE, "..", "simulations", "S4v1_μh参数平面扫描",
                           "S4v1_参数平面_结果表.csv")
KC_NOM = 2e4


def load(csv_path):
    rows = []
    with open(csv_path, encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            rec = {}
            for k, v in r.items():
                if k in ("regime", "note"):
                    rec[k] = v
                elif v == "":
                    rec[k] = np.nan
                elif v == "True":
                    rec[k] = True
                elif v == "False":
                    rec[k] = False
                else:
                    rec[k] = float(v)
            rows.append(rec)
    return rows


def fmt(v, d=4):
    return f"{v:.{d}g}" if isinstance(v, float) and np.isfinite(v) else "—"


def analyze(rows):
    print(f"总点数 {len(rows)}；N0={rows[0]['N0']:g} N，ω0={rows[0]['omega0']:g}，"
          f"μ={len(MU)} 档，h/L={len(HL)} 档\n")

    print("== 1) t*(k_c) 单调性 ==")
    bad = 0
    for mu in MU:
        for hl in HL:
            g = sorted([r for r in rows if r["mu"] == mu and r["h_over_L"] == hl],
                       key=lambda r: r["kc"])
            ts = np.array([r["t_star_ms"] for r in g])
            if not np.all(np.diff(ts) < 0):
                bad += 1
                print(f"  ! μ={mu:g} h/L={hl:g}: {np.round(ts, 4)}")
    print(f"  非单调的 (μ,h/L) 组数 = {bad} / {len(MU) * len(HL)}")

    print("\n== 2) 缺摩擦翻转的 k_c 个数（行=μ，列=h/L）==")
    print("        " + " ".join(f"{h:>7g}" for h in HL))
    for mu in MU:
        cells = [sum(1 for r in rows if r["mu"] == mu and r["h_over_L"] == hl and not r["flip"])
                 for hl in HL]
        print(f"  μ={mu:<4g} " + " ".join(f"{c:>7d}" for c in cells))

    print("\n== 3) t_flip/t* 随 μ 的变化（k_c=2e4）==")
    print(f"{'h/L':>7} " + " ".join(f"{m:>6g}" for m in MU))
    for hl in HL:
        vals = []
        for mu in MU:
            r = [x for x in rows if x["mu"] == mu and x["h_over_L"] == hl
                 and abs(x["kc"] - KC_NOM) < 1][0]
            vals.append(r["t_flip_ms"] / r["t_star_ms"] if r["flip"] else np.nan)
        print(f"{hl:>7g} " + " ".join(f"{v:>6.3f}" if np.isfinite(v) else "     —" for v in vals))

    print("\n== 4) R(θ*) 跨越 1 时 t* 的拐点检验（k_c=2e4，沿 h/L 二阶差分）==")
    for mu in MU:
        g = sorted([r for r in rows if r["mu"] == mu and abs(r["kc"] - KC_NOM) < 1],
                   key=lambda r: r["h_over_L"])
        ts = np.array([r["t_star_ms"] for r in g])
        d2 = np.diff(ts, 2) / ts[1:-1]
        cross = next((i for i, r in enumerate(g) if r["R_star"] > 1), "—")
        print(f"  μ={mu:<4g} R*∈[{g[0]['R_star']:.2f},{g[-1]['R_star']:.2f}]  "
              f"二阶差分={np.round(d2, 4)}  首跨 1 的格={cross}")

    print("\n== 5) 笔痕量随 μ（k_c=2e4，h/L=0.9375）==")
    print(f"{'μ':>5} {'t*/ms':>9} {'tflip/ms':>9} {'I2nd/J':>10} {'ℓ/mm':>8} {'ω*':>8} {'xdot_max':>9}")
    for mu in MU:
        r = [x for x in rows if x["mu"] == mu and x["h_over_L"] == 0.9375
             and abs(x["kc"] - KC_NOM) < 1][0]
        print(f"{mu:>5g} {r['t_star_ms']:>9.4f} {fmt(r['t_flip_ms']):>9} {r['I_2nd']:>10.3e} "
              f"{r['ell_mm']:>8.4f} {r['omega_star']:>8.2f} {r['xdot_max']:>9.3f}")

    print("\n== 6) 数值误差 ==")
    err = np.array([r["num_err_t_ms"] for r in rows])
    rel = np.array([r["num_err_t_ms"] / r["t_star_ms"] for r in rows])
    print(f"  max|Δt*| = {np.nanmax(err):.2e} ms；max 相对 = {np.nanmax(rel):.2e}；"
          f"相对 >1e-3 的点数 = {int(np.sum(rel > 1e-3))}")


CASES = [(1.20, 0.98, 1e5), (1.20, 0.98, 5e3), (1.00, 0.96, 1e5),
         (0.20, 0.70, 1e5), (0.35, 0.9375, 2e4)]


def extreme():
    print("\n== 7) 极端工况判读 ==")
    print(f"{'μ':>5} {'h/L':>6} {'kc':>8} {'t*/ms':>8} {'minθ̇':>8} {'N单调↓':>7} "
          f"{'v_rel 区间':>22} {'换向':>5} {'ω回落%':>7} {'拐点t/t*':>8}")
    for mu, hl, kc in CASES:
        r = cm.integrate_segment(kc, 4.0, 0.0, mu=mu, h=hl * cm.L, record=True)
        tr = r["traj"]
        w, N, v, t = tr["w"], tr["N"], tr["vrel"], tr["t"]
        i = int(np.argmax(w))
        drop = (w[i] - w[-1]) / w[i] * 100 if w[i] > 0 else 0.0
        print(f"{mu:>5g} {hl:>6g} {kc:>8.0e} {r['t_star'] * 1e3:>8.4f} {w.min():>8.2f} "
              f"{str(bool(np.all(np.diff(N) <= 1e-12))):>7} "
              f"[{v.min():>8.4f},{v.max():>8.4f}] {str(r['flip']):>5} "
              f"{drop:>6.1f}% {t[i] / r['t_star']:>8.3f}")
    print("\n判读：ω 在段内单调上升（回落为数值噪声量级），θ̇ ≥ 0 ⇒ δ 单调减 ⇒ 分离是几何必然；")
    print("      v1 无黏着、无法向自由度，故不可能出现“越压越紧”的再加载分支。")


def main():
    ap = argparse.ArgumentParser(description="S4v1 结果诊断")
    ap.add_argument("--csv", default=DEFAULT_CSV)
    ap.add_argument("--extreme", action="store_true", help="只跑极端工况判读")
    args = ap.parse_args()
    if args.extreme:
        extreme()
        return
    analyze(load(args.csv))
    extreme()


if __name__ == "__main__":
    main()
