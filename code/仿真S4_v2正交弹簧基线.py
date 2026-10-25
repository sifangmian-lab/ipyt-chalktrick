# -*- coding: utf-8 -*-
"""仿真 S4：v2 最小模型（弹性杆 + 手端正交弹簧）基线

模型见 `code/v2_model.py` 与 `docs/研究记录/模型v2规格_2026-10-25.md`。
本脚本输出：
  · 名义基线的 N(t)、v_rel(t)、δ(t)、单次划痕 I(ξ)；
  · 截断指标：平台强度 I_plateau、截断强度 I_c、比值 I_c/I_plateau、截断位置 ξ_c；
  · 手端柔度/速度对比（基线 B/C）。
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
from deposit_profile import deposit_profile, profile_integral, total_deposit  # noqa: F401

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


def plot_case(tag, info, tr, xc, I, met, outdir):
    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    t = tr["tt"] * 1e3
    ax[0, 0].plot(t, tr["N"], "g-", lw=1.6)
    ax[0, 0].axvline(info["t_stop_ms"], color="r", ls=":", lw=1.2,
                     label=f"停止点 {info['t_stop_ms']:.3f} ms")
    ax[0, 0].set_xlabel("t / ms"); ax[0, 0].set_ylabel("N / N")
    ax[0, 0].set_title(f"({tag}a) N(t)：{info['N_min']:.2f}–{info['N_max']:.2f} N，"
                       f"相对涨落 {info['N_rel_std']*100:.1f}%")
    ax[0, 0].legend(); ax[0, 0].grid(alpha=0.3)
    ax[0, 1].plot(t, tr["v_rel"], "b-", lw=1.4)
    ax[0, 1].axhline(0, color="k", lw=0.8, ls=":")
    ax[0, 1].axvline(info["t_stop_ms"], color="r", ls=":", lw=1.2)
    ax[0, 1].set_xlabel("t / ms"); ax[0, 1].set_ylabel("v_rel / (m/s)")
    ax[0, 1].set_title(f"({tag}b) 相对滑动速度：停止点后进入黏着")
    ax[0, 1].grid(alpha=0.3)
    ax[1, 0].plot(t, tr["delta"] * 1e3, "m-", lw=1.4)
    ax[1, 0].axvline(info["t_stop_ms"], color="r", ls=":", lw=1.2)
    ax[1, 0].set_xlabel("t / ms"); ax[1, 0].set_ylabel("δ / mm")
    ax[1, 0].set_title(f"({tag}c) 法向压入量 δ(t)（N = k_n δ + 阻尼）")
    ax[1, 0].grid(alpha=0.3)
    if xc.size:
        ax[1, 1].plot(xc * 1e3, I, "r-", lw=1.7, label="I(ξ)")
        ax[1, 1].axhline(met["I_plateau"], color="gray", ls="--", lw=1,
                         label=f"平台 {met['I_plateau']:.2f} J/m")
        ax[1, 1].axhline(met["I_c"], color="k", ls=":", lw=1,
                         label=f"截断 {met['I_c']:.2f} J/m")
        ax[1, 1].set_xlabel("ξ / mm"); ax[1, 1].set_ylabel("I(ξ) / (J/m)")
        ax[1, 1].set_title(f"({tag}d) 单次划痕 I(ξ)：I_c/I_plateau = "
                           f"{met['cutoff_ratio']:.3f}")
        ax[1, 1].legend(); ax[1, 1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"S4_{tag}_N与截断.png"), dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="./S4_v2输出")
    ap.add_argument("--tmax", type=float, default=0.004)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    cases = [
        ("基线A", dict(vm.DEFAULT)),
        ("基线B", dict(vm.DEFAULT, kx=1000.0, ky=1.0e5, ktheta=100.0)),
        ("基线C", dict(vm.DEFAULT, v_board=0.1)),
    ]
    rows = []
    profiles = {}
    for tag, p in cases:
        info, tr = vm.simulate(p, t_max=args.tmax)
        xc, I, met = vm.dash_profile(tr, nbins=120, sigma_bins=1.0)
        profiles[tag] = (xc, I)
        row = dict(case=tag, kx=p["kx"], ky=p["ky"], ktheta=p["ktheta"],
                   v_board=p["v_board"], friction=p["friction"],
                   sigma0=p["sigma0"], sigma1=p["sigma1"], sigma2=p["sigma2"],
                   **info, **met)
        rows.append(row)
        plot_case(tag, info, tr, xc, I, met, args.outdir)
        print(f"[{tag}] t_stop={info['t_stop_ms']:.3f} ms, "
              f"N={info['N_min']:.2f}-{info['N_max']:.2f} (rel {info['N_rel_std']*100:.1f}%), "
              f"ξ_end={info['xi_end_mm']:.3f} mm, I_plateau={met['I_plateau']:.3f}, "
              f"I_c={met['I_c']:.3f}, I_c/I_plateau={met['cutoff_ratio']:.3f}, "
              f"∫I/I2={met['I_conservation']:.4f}")
        if tag == "基线A":
            with open(os.path.join(args.outdir, "S4_基线A_强度剖面.csv"),
                      "w", newline="", encoding="utf-8-sig") as f:
                writer = csv.writer(f)
                writer.writerow(["xi_mm", "I_J_per_m"])
                for a, b in zip(xc, I):
                    writer.writerow([f"{a*1e3:.6e}", f"{b:.6e}"])

    with open(os.path.join(args.outdir, "S4_基线_结果表.csv"), "w",
              newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\n输出目录：{os.path.abspath(args.outdir)}")
    print("说明：基线 A/B 用 v=0.5 m/s 直接对照 v1 名义；C 用 0.05 m/s 看低速。")


if __name__ == "__main__":
    main()
