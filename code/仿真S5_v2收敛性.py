# -*- coding: utf-8 -*-
"""仿真 S5：v2 数值收敛性与正则化/参数敏感性（k_n、σ0、容差、接触阻尼）。

目的：检查 S4 基线 A 的“N 平台 + 有限强度截断”结论是否对数值参数稳健。
输出：simulations/S5_v2收敛性/ 下的 CSV、图与结果说明。
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


def run_case(tag, p, rtol=1e-8, tmax=0.004):
    info, tr = vm.simulate(p, t_max=tmax, rtol=rtol)
    xc, I, met = vm.dash_profile(tr, nbins=120, sigma_bins=1.0)
    return dict(case=tag, rtol=rtol,
                kn=p["kn"], sigma0=p["sigma0"], cn=p["cn"],
                **info, **met)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="./S5_v2输出")
    ap.add_argument("--tmax", type=float, default=0.004)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    base = dict(vm.DEFAULT)
    rows = []

    # 1) 接触刚度收敛
    for kn in [2e3, 2e4, 2e5, 2e6]:
        rows.append(run_case(f"kn={kn:.0e}", dict(base, kn=kn), tmax=args.tmax))
    # 2) 鬃毛刚度/正则化
    for s0 in [1e4, 1e5, 1e6]:
        rows.append(run_case(f"sigma0={s0:.0e}", dict(base, sigma0=s0),
                             tmax=args.tmax))
    # 3) 积分容差
    for rt in [1e-6, 1e-8, 1e-10]:
        rows.append(run_case(f"rtol={rt:.0e}", dict(base), rtol=rt,
                             tmax=args.tmax))
    # 4) 接触阻尼
    for cn in [0.1, 1.0, 10.0]:
        rows.append(run_case(f"cn={cn:g}", dict(base, cn=cn), tmax=args.tmax))

    with open(os.path.join(args.outdir, "S5_收敛性_结果表.csv"), "w",
              newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    for r in rows:
        print(f"{r['case']:>16s}: t_stop={r['t_stop_ms']:.3f} ms, "
              f"Nrel={r['N_rel_std']*100:.1f}%, ξ={r['xi_end_mm']:.3f} mm, "
              f"Ip={r['I_plateau']:.3f}, Ic={r['I_c']:.3f}, "
              f"Ic/Ip={r['cutoff_ratio']:.3f}, ∫/I2={r['I_conservation']:.4f}")

    # 参考值：基线 A（kn=2e4, sigma0=1e5, rtol=1e-8, cn=1）
    ref = [r for r in rows if r["case"] == "kn=2e+04"][0]

    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    kn_rows = [r for r in rows if r["case"].startswith("kn=")]
    ax[0, 0].plot([r["kn"] for r in kn_rows], [r["t_stop_ms"] for r in kn_rows],
                  "o-", label="t_stop / ms")
    ax[0, 0].set_xscale("log"); ax[0, 0].set_xlabel("k_n / (N/m)")
    ax[0, 0].set_ylabel("t_stop / ms")
    ax[0, 0].set_title("(a) 停止时间 vs 接触刚度")
    ax[0, 0].grid(alpha=0.3, which="both"); ax[0, 0].legend()
    ax[0, 1].plot([r["kn"] for r in kn_rows],
                  [r["cutoff_ratio"] for r in kn_rows], "s-", label="I_c/I_plateau")
    ax[0, 1].plot([r["kn"] for r in kn_rows],
                  [r["N_rel_std"] for r in kn_rows], "^-", label="N 相对涨落")
    ax[0, 1].set_xscale("log"); ax[0, 1].set_xlabel("k_n / (N/m)")
    ax[0, 1].set_title("(b) 截断平台比与 N 涨落 vs k_n")
    ax[0, 1].grid(alpha=0.3, which="both"); ax[0, 1].legend(fontsize=8)
    s_rows = [r for r in rows if r["case"].startswith("sigma0=")]
    ax[1, 0].plot([r["sigma0"] for r in s_rows],
                  [r["I_plateau"] for r in s_rows], "o-", label="I_plateau")
    ax[1, 0].plot([r["sigma0"] for r in s_rows],
                  [r["I_c"] for r in s_rows], "s-", label="I_c")
    ax[1, 0].set_xscale("log"); ax[1, 0].set_xlabel("σ0 / (N/m)")
    ax[1, 0].set_ylabel("I / (J/m)")
    ax[1, 0].set_title("(c) 强度 vs 鬃毛刚度 σ0")
    ax[1, 0].grid(alpha=0.3, which="both"); ax[1, 0].legend(fontsize=8)
    rt_rows = [r for r in rows if r["case"].startswith("rtol=")]
    dev = [max(abs(r["t_stop_ms"] - ref["t_stop_ms"]) / ref["t_stop_ms"],
               abs(r["xi_end_mm"] - ref["xi_end_mm"]) / ref["xi_end_mm"],
               abs(r["I_plateau"] - ref["I_plateau"]) / ref["I_plateau"])
           for r in rt_rows]
    ax[1, 1].semilogx([r["rtol"] for r in rt_rows], dev, "o-")
    ax[1, 1].set_xlabel("rtol"); ax[1, 1].set_ylabel("相对偏差（vs rtol=1e-8）")
    ax[1, 1].set_title("(d) 积分容差收敛")
    ax[1, 1].grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "S5_收敛性.png"), dpi=150)
    plt.close(fig)
    print(f"\n输出目录：{os.path.abspath(args.outdir)}")


if __name__ == "__main__":
    main()
