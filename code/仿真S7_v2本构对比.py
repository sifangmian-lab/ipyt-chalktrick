# -*- coding: utf-8 -*-
"""仿真 S7：v2 本构对比（LuGre-springslider / LuGre-viscous / tanh）。

目的：在相同几何与手端柔度下，比较不同切向摩擦本构对
  N(t)、v_rel(t)、F_f(t)、I(ξ) 与“平台 + 有限截断”的影响。
注意：tanh 无静摩擦，不能真正黏着；该对比用于展示本构敏感性，而非选型定论。
输出：simulations/S7_v2本构对比/ 下的 CSV、图和结果说明。
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as _fm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v2_model as vm

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
                                   "Noto Sans CJK SC", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="./S7_v2输出")
    ap.add_argument("--tmax", type=float, default=0.004)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    base = dict(vm.DEFAULT)
    cases = {
        "LuGre-springslider": dict(base, friction="lugre_springslider"),
        "LuGre-viscous": dict(base, friction="lugre", sigma1=30.0, sigma2=5.0),
        "tanh": dict(base, friction="tanh", eps_tanh=0.02),
    }
    rows, trajs = [], {}
    for tag, p in cases.items():
        info, tr = vm.simulate(p, t_max=args.tmax)
        xc, I, met = vm.dash_profile(tr, nbins=120, sigma_bins=1.0)
        trajs[tag] = (info, tr, xc, I, met)
        rows.append(dict(case=tag, **{k: p[k] for k in
                                      ("friction", "sigma0", "sigma1", "sigma2",
                                       "mu_s", "mu_k", "eps_tanh")},
                         **info, **met))
        print(f"[{tag}] t_stop={info['t_stop_ms']:.3f} ms, "
              f"N={info['N_min']:.2f}-{info['N_max']:.2f}, "
              f"ξ={info['xi_end_mm']:.3f} mm, Ip={met['I_plateau']:.3f}, "
              f"Ic={met['I_c']:.3f}, Ic/Ip={met['cutoff_ratio']:.3f}, "
              f"∫/I2={met['I_conservation']:.4f}")

    with open(os.path.join(args.outdir, "S7_v2本构_结果表.csv"), "w",
              newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    colors = {"LuGre-springslider": "C0", "LuGre-viscous": "C3", "tanh": "k"}
    for tag, (info, tr, xc, I, met) in trajs.items():
        c = colors[tag]
        ax[0, 0].plot(tr["tt"] * 1e3, tr["N"], c, lw=1.5, label=tag)
        ax[0, 1].plot(tr["tt"] * 1e3, tr["v_rel"], c, lw=1.5, label=tag)
        ax[1, 0].plot(tr["tt"] * 1e3, tr["F_f"], c, lw=1.5, label=tag)
        if xc.size:
            ax[1, 1].plot(xc * 1e3, I, c, lw=1.6, label=tag)
    ax[0, 0].set_xlabel("t / ms"); ax[0, 0].set_ylabel("N / N")
    ax[0, 0].set_title("(a) 法向力 N(t)"); ax[0, 0].grid(alpha=0.3); ax[0, 0].legend()
    ax[0, 1].set_xlabel("t / ms"); ax[0, 1].set_ylabel("v_rel / (m/s)")
    ax[0, 1].set_title("(b) 相对滑动速度"); ax[0, 1].grid(alpha=0.3); ax[0, 1].legend()
    ax[1, 0].set_xlabel("t / ms"); ax[1, 0].set_ylabel("F_f / N")
    ax[1, 0].set_title("(c) 摩擦力"); ax[1, 0].grid(alpha=0.3); ax[1, 0].legend()
    ax[1, 1].set_xlabel("ξ / mm"); ax[1, 1].set_ylabel("I(ξ) / (J/m)")
    ax[1, 1].set_title("(d) 单次划痕空间强度"); ax[1, 1].grid(alpha=0.3)
    ax[1, 1].legend()
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "S7_v2本构对比.png"), dpi=150)
    plt.close(fig)
    print(f"\n输出目录：{os.path.abspath(args.outdir)}")


if __name__ == "__main__":
    main()
