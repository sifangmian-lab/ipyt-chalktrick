# -*- coding: utf-8 -*-
"""仿真 S6：v2 参数平面（手端柔度、弯曲刚度、板速）与“平台/截断”分类。

输出：
  · (k_x, k_y) 平面：N 相对涨落、I_c/I_plateau、t_stop 的分布；
  · k_θ 与 v 的一维扫描；
  · CSV 全表 + 结果说明（simulations/S6_v2参数平面/）。
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


def run_case(tag, p, tmax=0.004):
    info, tr = vm.simulate(p, t_max=tmax)
    xc, I, met = vm.dash_profile(tr, nbins=120, sigma_bins=1.0)
    return dict(case=tag, kx=p["kx"], ky=p["ky"], ktheta=p["ktheta"],
                v_board=p["v_board"], **info, **met)


def classify(r):
    platform = r["N_rel_std"] < 0.10
    cutoff = r["cutoff_ratio"] > 0.80
    if platform and cutoff:
        return "平台+截断"
    if platform:
        return "平台/弱截断"
    if cutoff:
        return "N变化/有截断"
    return "N变化/无截断"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="./S6_v2输出")
    ap.add_argument("--tmax", type=float, default=0.004)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    base = dict(vm.DEFAULT)
    rows = []

    kxs = [100.0, 300.0, 1000.0]
    kys = [3e5, 1e6, 3e6]
    for kx in kxs:
        for ky in kys:
            p = dict(base, kx=kx, ky=ky)
            rows.append(run_case(f"kx={kx:g},ky={ky:.0e}", p, args.tmax))
    for kth in [100.0, 300.0, 1000.0, 3000.0]:
        rows.append(run_case(f"ktheta={kth:g}", dict(base, ktheta=kth), args.tmax))
    for v in [0.05, 0.1, 0.2, 0.5, 1.0]:
        rows.append(run_case(f"v={v:g}", dict(base, v_board=v), args.tmax))

    for r in rows:
        r["class"] = classify(r)
    with open(os.path.join(args.outdir, "S6_v2参数平面_结果表.csv"), "w",
              newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    for r in rows:
        print(f"{r['case']:>18s}: Nrel={r['N_rel_std']*100:5.1f}%, "
              f"Ic/Ip={r['cutoff_ratio']:.3f}, t_stop={r['t_stop_ms']:.3f} ms, "
              f"ξ={r['xi_end_mm']:.3f} mm -> {r['class']}")

    # 平面图
    grid = {k: np.full((len(kxs), len(kys)), np.nan) for k in
            ("N_rel_std", "cutoff_ratio", "t_stop_ms")}
    for r in rows:
        if r["case"].startswith("kx="):
            i = kxs.index(r["kx"]); j = kys.index(r["ky"])
            for k in grid:
                grid[k][i, j] = r[k]
    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    for a, key, title in zip(
            [ax[0, 0], ax[0, 1], ax[1, 0]],
            ["N_rel_std", "cutoff_ratio", "t_stop_ms"],
            ["(a) N 相对涨落", "(b) I_c/I_plateau", "(c) t_stop / ms"]):
        im = a.imshow(grid[key], origin="lower", aspect="auto", cmap="viridis")
        a.set_xticks(range(len(kys)), [f"{k:.0e}" for k in kys])
        a.set_yticks(range(len(kxs)), [f"{k:g}" for k in kxs])
        a.set_xlabel("k_y / (N/m)"); a.set_ylabel("k_x / (N/m)")
        a.set_title(title)
        fig.colorbar(im, ax=a)
    # k_theta / v 扫描
    kth_rows = [r for r in rows if r["case"].startswith("ktheta=")]
    v_rows = [r for r in rows if r["case"].startswith("v=")]
    ax[1, 1].plot([r["ktheta"] for r in kth_rows],
                  [r["cutoff_ratio"] for r in kth_rows], "o-", label="I_c/I_plateau")
    ax[1, 1].plot([r["ktheta"] for r in kth_rows],
                  [1 - r["N_rel_std"] for r in kth_rows], "s-", label="1 - N 相对涨落")
    ax[1, 1].set_xscale("log"); ax[1, 1].set_xlabel("k_θ / (N·m/rad)")
    ax[1, 1].set_title("(d) 弯曲刚度扫描（k_x=300, k_y=1e6）")
    ax[1, 1].set_ylim(0, 1.05); ax[1, 1].grid(alpha=0.3); ax[1, 1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "S6_v2参数平面.png"), dpi=150)
    plt.close(fig)

    # v 扫描单独图
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    ax[0].semilogx([r["v_board"] for r in v_rows],
                   [r["cutoff_ratio"] for r in v_rows], "o-")
    ax[0].set_xlabel("v / (m/s)"); ax[0].set_ylabel("I_c/I_plateau")
    ax[0].set_title("截断平台比 vs 板速"); ax[0].grid(alpha=0.3, which="both")
    ax[1].semilogx([r["v_board"] for r in v_rows],
                   [r["xi_end_mm"] for r in v_rows], "s-")
    ax[1].set_xlabel("v / (m/s)"); ax[1].set_ylabel("ξ_end / mm")
    ax[1].set_title("划痕长度 vs 板速"); ax[1].grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "S6_板速扫描.png"), dpi=150)
    plt.close(fig)
    print(f"\n输出目录：{os.path.abspath(args.outdir)}")


if __name__ == "__main__":
    main()
