# -*- coding: utf-8 -*-
"""仿真 S5b：v2 法向接触模型升级与刚性极限（LCP / 事件驱动）。

背景
----
S5 报告"k_n 刚性极限不收敛"（k_n 2e3 -> 2e6 时 N 涨落 0.3% -> 47.8%）。
A2 核查发现两件事：
  1. v2_model.rhs_factory 的手端 y 弹簧预载符号写反（y_anchor = y0 + (N0+Mg)/ky），
     使初态净法向力为 +2N0（向上 8 N）。修正为 y0 - (N0-Mg)/ky 后，
     基线 A 的 N 相对涨落从 3.1% 降到 0.25%。
  2. 罚接触在 k_n -> ∞ 时接触力振荡频率 ~ sqrt(k_n/M) 发散、振幅不消失，
     宏观量没有唯一极限；因此"继续增大 k_n"不是求刚性极限的正确做法。
     本脚本用指标 1 的单边约束（contact='rigid'）给出 k_n -> ∞ 的极限参考。

输出
----
  simulations/S5b_v2接触极限/
    S5b_接触模型_结果表.csv
    S5b_接触极限对比.png
    S5b_刚性接触验证.png
    结果说明.md（人工撰写，见仓库）
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

FIELDS = ["case", "model", "kn", "cn", "ky", "rtol", "max_step",
          "t_stop_ms", "N_min", "N_max", "N_mean", "N_rel_std", "v_rel_min",
          "xi_end_mm", "I_plateau", "I_c", "cutoff_ratio", "I_conservation",
          "xi_c_mm", "release_events", "constraint_residual_max", "n0_check"]


def run_case(tag, model, p, rtol=1e-8, max_step=1e-5, tmax=0.004):
    info, tr = vm.simulate(p, t_max=tmax, rtol=rtol, max_step=max_step)
    xc, I, met = vm.dash_profile(tr, nbins=120, sigma_bins=1.0)
    row = dict(
        case=tag, model=model, kn=p.get("kn"), cn=p.get("cn"), ky=p.get("ky"),
        rtol=rtol, max_step=max_step,
        t_stop_ms=info["t_stop_ms"], N_min=info["N_min"], N_max=info["N_max"],
        N_mean=info["N_mean"], N_rel_std=info["N_rel_std"],
        v_rel_min=info["v_rel_min"], xi_end_mm=info["xi_end_mm"],
        I_plateau=met["I_plateau"], I_c=met["I_c"],
        cutoff_ratio=met["cutoff_ratio"], I_conservation=met["I_conservation"],
        xi_c_mm=met["xi_c_mm"],
        release_events=info.get("release_events", ""),
        constraint_residual_max=info.get("constraint_residual_max", ""),
        n0_check=info.get("n0_check", ""),
    )
    return row, xc, I, info, tr, met


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="./S5b_v2输出")
    ap.add_argument("--tmax", type=float, default=0.004)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    base = dict(vm.DEFAULT)
    rows = []
    profiles = {}

    # ---------------- 1. 罚接触：k_n 扫描 × 两档数值分辨率 ----------------
    for kn in [2e3, 2e4, 2e5, 2e6]:
        for rtol, mstep in [(1e-8, 1e-5), (1e-10, 1e-6)]:
            tag = f"penalty_kn={kn:.0e}_rtol={rtol:.0e}_ms={mstep:.0e}"
            row, xc, I, info, tr, met = run_case(
                tag, "penalty", dict(base, kn=kn), rtol=rtol,
                max_step=mstep, tmax=args.tmax)
            rows.append(row)
            if kn == 2e4 and rtol == 1e-8:
                profiles["penalty"] = (xc, I, info, tr, met)

    # ---------------- 2. 罚接触：接触阻尼扫描 ----------------
    for cn in [0.1, 1.0, 10.0]:
        tag = f"penalty_cn={cn:g}"
        row, xc, I, info, tr, met = run_case(
            tag, "penalty", dict(base, cn=cn), tmax=args.tmax)
        rows.append(row)

    # ---------------- 3. 罚接触：k_y 扫描（手端法向柔度） ----------------
    for ky in [1e5, 3e5, 1e6, 3e6]:
        tag = f"penalty_ky={ky:.0e}"
        row, xc, I, info, tr, met = run_case(
            tag, "penalty", dict(base, ky=ky), tmax=args.tmax)
        rows.append(row)

    # ---------------- 4. 刚性接触：双容差 / max_step ----------------
    rigid_std = dict(base, contact="rigid")
    for rtol, mstep in [(1e-8, 1e-5), (1e-10, 1e-5),
                        (1e-10, 1e-6), (1e-12, 1e-6)]:
        tag = f"rigid_rtol={rtol:.0e}_ms={mstep:.0e}"
        row, xc, I, info, tr, met = run_case(
            tag, "rigid", rigid_std, rtol=rtol, max_step=mstep,
            tmax=args.tmax)
        rows.append(row)
        if rtol == 1e-8 and mstep == 1e-5:
            profiles["rigid"] = (xc, I, info, tr, met)

    # ---------------- 5. 刚性接触：k_y 扫描 ----------------
    for ky in [1e5, 3e5, 1e6, 3e6]:
        tag = f"rigid_ky={ky:.0e}"
        row, xc, I, info, tr, met = run_case(
            tag, "rigid", dict(rigid_std, ky=ky), tmax=args.tmax)
        rows.append(row)

    # ---------------- 写 CSV ----------------
    with open(os.path.join(args.outdir, "S5b_接触模型_结果表.csv"), "w",
              newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: r.get(k, "") for k in FIELDS})

    for r in rows:
        rr = r["constraint_residual_max"]
        rr = f"{rr:.1e}" if isinstance(rr, float) else rr
        print(f"{r['case']:>42s}: t_stop={r['t_stop_ms']:8.4f} ms  "
              f"Nrel={r['N_rel_std']*100:6.2f}%  Nmin={r['N_min']:6.3f}  "
              f"xi={r['xi_end_mm']:7.4f} mm  Ic/Ip={r['cutoff_ratio']:6.3f}  "
              f"∫/I2={r['I_conservation']:8.5f}  rel={r['release_events']}  "
              f"res={rr}")

    # ---------------- 图 1：罚 vs 刚性 ----------------
    pen = {r["kn"]: r for r in rows
           if r["model"] == "penalty" and r["rtol"] == 1e-8
           and str(r["case"]).startswith("penalty_kn=")}
    kns = sorted(pen)
    rigid_rows = [r for r in rows if r["model"] == "rigid"
                  and str(r["case"]).startswith("rigid_rtol=")]
    rigid_ref = rigid_rows[0]
    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    ax[0, 0].semilogx(kns, [pen[k]["t_stop_ms"] for k in kns], "o-",
                      label="罚接触（修正预载后）")
    ax[0, 0].axhline(rigid_ref["t_stop_ms"], color="r", ls="--",
                     label=f"刚性接触 t_stop={rigid_ref['t_stop_ms']:.3f} ms")
    ax[0, 0].set_xlabel("k_n / (N/m)"); ax[0, 0].set_ylabel("t_stop / ms")
    ax[0, 0].set_title("(a) 停止时间 vs 接触刚度")
    ax[0, 0].grid(alpha=0.3, which="both"); ax[0, 0].legend(fontsize=8)
    ax[0, 1].semilogx(kns, [pen[k]["xi_end_mm"] for k in kns], "s-",
                      label="罚接触")
    ax[0, 1].axhline(rigid_ref["xi_end_mm"], color="r", ls="--",
                     label=f"刚性接触 ξ={rigid_ref['xi_end_mm']:.3f} mm")
    ax[0, 1].set_xlabel("k_n / (N/m)"); ax[0, 1].set_ylabel("ξ_end / mm")
    ax[0, 1].set_title("(b) 划痕长度 vs 接触刚度")
    ax[0, 1].grid(alpha=0.3, which="both"); ax[0, 1].legend(fontsize=8)
    ax[1, 0].semilogx(kns, [pen[k]["N_rel_std"] * 100 for k in kns], "^-",
                      label="罚接触")
    ax[1, 0].axhline(rigid_ref["N_rel_std"] * 100, color="r", ls="--",
                     label=f"刚性接触 {rigid_ref['N_rel_std']*100:.1f}%")
    ax[1, 0].set_xlabel("k_n / (N/m)"); ax[1, 0].set_ylabel("N 相对涨落 / %")
    ax[1, 0].set_title("(c) N 平台涨落 vs 接触刚度")
    ax[1, 0].grid(alpha=0.3, which="both"); ax[1, 0].legend(fontsize=8)
    pky = [r for r in rows if r["model"] == "penalty"
           and str(r["case"]).startswith("penalty_ky=")]
    rky = [r for r in rows if r["model"] == "rigid"
           and str(r["case"]).startswith("rigid_ky=")]
    ax[1, 1].semilogx([r["ky"] for r in pky], [r["t_stop_ms"] for r in pky],
                      "o-", label="罚接触 k_n=2e4")
    ax[1, 1].semilogx([r["ky"] for r in rky], [r["t_stop_ms"] for r in rky],
                      "s--", label="刚性接触")
    ax[1, 1].set_xlabel("k_y / (N/m)"); ax[1, 1].set_ylabel("t_stop / ms")
    ax[1, 1].set_title("(d) 手端法向柔度扫描")
    ax[1, 1].grid(alpha=0.3, which="both"); ax[1, 1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "S5b_接触极限对比.png"), dpi=150)
    plt.close(fig)

    # ---------------- 图 2：刚性接触验证 ----------------
    xc_p, I_p, info_p, tr_p, met_p = profiles["penalty"]
    xc_r, I_r, info_r, tr_r, met_r = profiles["rigid"]
    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    ax[0, 0].plot(tr_p["tt"] * 1e3, tr_p["N"], label="罚接触 k_n=2e4")
    ax[0, 0].plot(tr_r["tt"] * 1e3, tr_r["N"], "--", label="刚性接触")
    ax[0, 0].set_xlabel("t / ms"); ax[0, 0].set_ylabel("N / N")
    ax[0, 0].set_title("(a) 法向力 N(t)：平台形态")
    ax[0, 0].grid(alpha=0.3); ax[0, 0].legend(fontsize=8)
    ax[0, 1].plot(tr_p["tt"] * 1e3, tr_p["v_rel"], label="罚接触 k_n=2e4")
    ax[0, 1].plot(tr_r["tt"] * 1e3, tr_r["v_rel"], "--", label="刚性接触")
    ax[0, 1].axhline(0, color="k", lw=0.8, ls=":")
    ax[0, 1].set_xlabel("t / ms"); ax[0, 1].set_ylabel("v_rel / (m/s)")
    ax[0, 1].set_title("(b) 相对速度 v_rel(t)")
    ax[0, 1].grid(alpha=0.3); ax[0, 1].legend(fontsize=8)
    if xc_p.size:
        ax[1, 0].plot(xc_p * 1e3, I_p, label="罚接触 k_n=2e4")
    if xc_r.size:
        ax[1, 0].plot(xc_r * 1e3, I_r, "--", label="刚性接触")
    ax[1, 0].set_xlabel("ξ / mm"); ax[1, 0].set_ylabel("I(ξ) / (J/m)")
    ax[1, 0].set_title("(c) 单次划痕强度 I(ξ)")
    ax[1, 0].grid(alpha=0.3); ax[1, 0].legend(fontsize=8)
    # 刚性接触的约束加速度残差随时间（应为数值零）
    vx_r = np.gradient(tr_r["x"], tr_r["tt"])
    w_r = np.gradient(tr_r["th"], tr_r["tt"])
    res = []
    for i in range(len(tr_r["tt"])):
        s = (tr_r["x"][i], tr_r["th"][i], vx_r[i], w_r[i], tr_r["z"][i])
        ri = vm._rigid_terms(rigid_std, s,
                             tr_r["meta"]["y_anchor"])[3]["residual"]
        res.append(abs(ri))
    ax[1, 1].semilogy(tr_r["tt"] * 1e3, np.maximum(res, 1e-18), "-")
    ax[1, 1].set_xlabel("t / ms"); ax[1, 1].set_ylabel("|residual| / (m/s²)")
    ax[1, 1].set_title("(d) 刚性接触约束加速度残差")
    ax[1, 1].grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "S5b_刚性接触验证.png"), dpi=150)
    plt.close(fig)

    # ---------------- 验收判据汇总 ----------------
    print("\n================ 验收判据 ================")
    print(f"[刚性] 约束加速度残差 max = {rigid_ref['constraint_residual_max']:.2e}"
          f"  (<1e-9: {'PASS' if rigid_ref['constraint_residual_max'] < 1e-9 else 'FAIL'})")
    print(f"[刚性] 初始法向力 n0_check = {rigid_ref['n0_check']:.6f} N")
    print(f"[刚性] 接触释放事件数 = {rigid_ref['release_events']}（基线应为 0）")
    print(f"[刚性] N_min = {rigid_ref['N_min']:.4f} N（应 ≥ 0）")
    t_r = [r["t_stop_ms"] for r in rigid_rows]
    print(f"[刚性] 双容差/步长 t_stop 极差 = {max(t_r)-min(t_r):.2e} ms")
    # 罚接触：同一 k_n 两档分辨率的差异
    print("[罚接触] 同 k_n、两档数值分辨率的 t_stop 相对差：")
    for kn in kns:
        a = [r for r in rows if r["model"] == "penalty"
             and r["kn"] == kn and r["rtol"] == 1e-8][0]
        b = [r for r in rows if r["model"] == "penalty"
             and r["kn"] == kn and r["rtol"] == 1e-10][0]
        print(f"   k_n={kn:.0e}: Δt/t = {abs(a['t_stop_ms']-b['t_stop_ms'])/a['t_stop_ms']:.1e}, "
              f"Δξ/ξ = {abs(a['xi_end_mm']-b['xi_end_mm'])/a['xi_end_mm']:.1e}")
    print(f"\n输出目录：{os.path.abspath(args.outdir)}")


if __name__ == "__main__":
    main()
