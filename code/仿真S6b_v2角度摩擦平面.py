# -*- coding: utf-8 -*-
"""仿真 S6b：v2 角度—摩擦系数平面与“连续滑移 / 黏着截断 / 脱离 / 自锁颤振”判别。

背景
----
S6 在固定 theta0=20° 下给出 (k_x,k_y) 平台的宽可达区，并把“角度/摩擦系数平面”
留给下一轮 S6b；S5b 建议用刚性接触（contact='rigid'，LCP/事件驱动）作 k_n→∞
的极限参考。本脚本据此扫描 (theta0, mu_s) 平面，并用三档罚接触刚度
k_n = 2e4 / 2e5 / 2e6 与刚性接触对比，回答：

  1. “平台+有限强度截断”在角度—摩擦平面上的可达区；
  2. 低摩擦端何时不再出现停止点（连续滑移 → 实线，而不是虚线点痕）；
  3. 高倾角/低摩擦端何时在停止前发生接触脱离（lift-off，sprag 型出口）；
  4. 理想刚性接触（k_n→∞）在同一平面上的适用范围与病态区（约束乘子增长/颤振）；
  5. 与 sprag-slip / 静摩擦自锁判据 mu_s > tan(theta0) 及动摩擦判据
     mu_k > tan(theta0) 的定性/定量对照。

输出（默认 --outdir ./S6b_v2输出）
---------------------------------
  S6b_角度摩擦平面_结果表.csv    全平面（4 种接触口径 × theta0 × mu_s）
  S6b_边界与判据.csv             连续滑移边界、脱离边界、刚性接触分界与判据值
  S6b_判据表.csv                 tan(theta0)、静/动摩擦判据对应的 mu_s
  S6b_摩擦比扫描.csv             mu_s/mu_k 比值的影响
  S6b_板速扫描.csv               板速对单事件停止边界的影响（theta0=20°, kn=2e4）
  S6b_数值校验.csv               双容差复算与 tmax 敏感性
  S6b_典型轨迹.csv / .png        四类工况的 N(t)、v_rel(t)、I(ξ)
  S6b_角度摩擦相图.png           四联分类相图（含两条判据曲线）
  S6b_罚模型指标.png             罚模型 N 涨落 / 截断比 / t_stop / 角度切面
  S6b_边界与判据.png             边界—判据对比、刚性相图、板速扫描
  结果说明.md（人工撰写，见仓库）

数值约定
--------
  · 罚接触：kn=2e4/2e5/2e6 N/m，cn=1，t_max=12 ms；
  · 刚性接触：约束乘子 + 释放事件；N>25 N 或 |v_rel|>2.5 m/s 早停，避免
    理想刚性/无冲击律区间的无界振铃拖死积分（该早停本身即“自锁/模型未闭合”信号）；
  · 分类：
      plateau 平台+截断（停止、N 相对涨落<10%、I_c/I_p>0.8）
      nvar    N变化+截断（停止，但平台判据不满足）
      chatter 停止后颤振（停止后 0.3 ms 内 max|v_rel|>0.1 m/s）
      slide   连续滑移（t_max 内无停止点；低摩擦无法把笔尖带到板速）
      lift    接触脱离（停止前发生接触分离）
      sprag   自锁/约束发散（刚性接触早停：N 或 |v_rel| 超限）
  · 摩擦比固定 mu_s/mu_k = 0.62/0.51 = 1.2157（与 S4–S7 名义一致）；另做比值扫描。
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
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from scipy.integrate import solve_ivp

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

MUS_RATIO = vm.DEFAULT["mu_s"] / vm.DEFAULT["mu_k"]     # 1.215686…

CLASS_ORDER = ["plateau", "nvar", "chatter", "slide", "lift", "sprag"]
CLASS_LABEL = {
    "plateau": "平台+截断",
    "nvar": "N变化+截断",
    "chatter": "停止后颤振",
    "slide": "连续滑移（无停止）",
    "lift": "接触脱离（lift-off）",
    "sprag": "自锁/约束发散",
}
CLASS_COLOR = {
    "plateau": "#2e7d32",
    "nvar": "#f9a825",
    "chatter": "#ef6c00",
    "slide": "#1565c0",
    "lift": "#c62828",
    "sprag": "#6a1b9a",
}

MODELS = [
    ("penalty_2e4", "罚接触 k_n=2e4", dict(kn=2e4)),
    ("penalty_2e5", "罚接触 k_n=2e5", dict(kn=2e5)),
    ("penalty_2e6", "罚接触 k_n=2e6", dict(kn=2e6)),
    ("rigid", "刚性接触 LCP", dict(contact="rigid")),
]
THETAS = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 55.0]
MUS = [0.15, 0.19, 0.22, 0.25, 0.30, 0.40, 0.50, 0.62, 0.80, 1.00]


# ---------------------------------------------------------------------------
# 刚性接触：带早停的积分（复制 simulate_rigid 的约束乘子算法，仅加事件）
# ---------------------------------------------------------------------------
def simulate_rigid_capped(p, t_max=0.004, rtol=1e-8, max_step=1e-5, npts=4001,
                          v_stop_tol=1e-3, t_ignore=5e-5,
                          cap_N=25.0, cap_vrel=2.5):
    """返回 (info, traj)，字段与 v2_model.simulate_rigid 兼容，另加早停标志。"""
    rhs, y0, meta = vm.rigid_factory(p)
    y_anchor = meta["y_anchor"]

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

    sol = solve_ivp(rhs, (0.0, t_max), y0, method="DOP853", rtol=rtol,
                    atol=1e-10, events=(ev_release, ev_capN, ev_capv),
                    dense_output=True, max_step=max_step)
    tt = np.linspace(0.0, sol.t[-1], npts)
    x, th, vx, w, z = sol.sol(tt)
    cth, sth = np.cos(th), np.sin(th)
    y = vm.L * cth
    v_rel = vx + vm.L * cth * w + p["v_board"]
    N = np.empty_like(tt)
    residual = np.empty_like(tt)
    for i in range(len(tt)):
        _, _, Ni, aux = vm._rigid_terms(p, (x[i], th[i], vx[i], w[i], z[i]),
                                        y_anchor)
        N[i] = Ni
        residual[i] = aux["residual"]
    F = p["sigma0"] * z
    q = np.abs(F * v_rel)
    xi = np.concatenate([[0.0],
                         np.cumsum(0.5 * (v_rel[:-1] + v_rel[1:]) * np.diff(tt))])
    i_stop = len(tt) - 1
    for i in range(1, len(tt)):
        if tt[i] > t_ignore and (v_rel[i] <= 0.0 or abs(v_rel[i]) < v_stop_tol):
            i_stop = i
            break
    info = dict(
        t_end_ms=float(tt[-1] * 1e3),
        touch_end=bool(sol.t_events[0].size > 0),
        t_stop_ms=float(tt[i_stop] * 1e3), i_stop=int(i_stop),
        N_min=float(N[:i_stop + 1].min()),
        N_max=float(N[:i_stop + 1].max()),
        N_mean=float(N[:i_stop + 1].mean()),
        N_rel_std=float(N[:i_stop + 1].std() / max(N[:i_stop + 1].mean(), 1e-12)),
        v_rel_min=float(v_rel[:i_stop + 1].min()),
        v_rel_max=float(v_rel[:i_stop + 1].max()),
        stick_frac=float(np.mean(np.abs(v_rel[:i_stop + 1]) < v_stop_tol)),
        xi_span_mm=float((xi[:i_stop + 1].max() - xi[:i_stop + 1].min()) * 1e3),
        xi_end_mm=float(xi[i_stop] * 1e3), contact="rigid",
        release_events=int(sol.t_events[0].size),
        constraint_residual_max=(float(np.max(np.abs(residual)))
                                 if residual.size else 0.0),
        n0_check=float(meta["n0_check"]),
        cap_N_hit=bool(sol.t_events[1].size > 0),
        cap_v_hit=bool(sol.t_events[2].size > 0),
    )
    traj = dict(tt=tt, N=N, v_rel=v_rel, delta=np.zeros_like(tt), F_f=F,
                xi=xi, q=q, th=th, x=x, y=y, z=z, i_stop=i_stop, meta=meta)
    return info, traj


# ---------------------------------------------------------------------------
# 单工况运行与分类
# ---------------------------------------------------------------------------
def run_case(model_key, p, tmax, cap_N=25.0, cap_vrel=2.5,
             rtol=1e-8, max_step=1e-5):
    """返回 (row, traj, xc, I)；row 含分类与全部指标。

    npts 按 tmax 等比放大，使时间采样步长与 4 ms 基线一致（~1 µs），
    避免长积分时 t_stop 被粗网格量化。
    """
    npts = max(4001, int(round(4001 * tmax / 0.004)) + 1)
    if model_key == "rigid":
        info, tr = simulate_rigid_capped(p, t_max=tmax, rtol=rtol,
                                         max_step=max_step, cap_N=cap_N,
                                         cap_vrel=cap_vrel, npts=npts)
    else:
        info, tr = vm.simulate(p, t_max=tmax, rtol=rtol, max_step=max_step,
                               npts=npts)
    xc, I, met = vm.dash_profile(tr, nbins=120, sigma_bins=1.0)
    tt, vr = tr["tt"], tr["v_rel"]
    i = info["i_stop"]
    lost = bool(info.get("touch_end", False))
    capN = bool(info.get("cap_N_hit", False))
    capv = bool(info.get("cap_v_hit", False))
    stopped = i < len(tt) - 1 and vr[i] <= 1e-3
    post_osc = np.nan
    if stopped:
        dt = tt[1] - tt[0]
        j1 = min(len(tt), i + max(5, int(0.3e-3 / dt)))
        post_osc = float(np.max(np.abs(vr[i:j1])))
    if lost:
        cls = "lift"
    elif capN or capv:
        cls = "sprag"
    elif not stopped:
        cls = "slide"
    elif post_osc > 0.25 * p["v_board"]:
        cls = "chatter"
    elif info["N_rel_std"] < 0.10 and met["cutoff_ratio"] > 0.80:
        cls = "plateau"
    else:
        cls = "nvar"
    row = dict(
        case=f"{model_key}_th={p['theta0_deg']:g}_mus={p['mu_s']:.3f}",
        model=model_key, theta0_deg=float(p["theta0_deg"]),
        mu_s=float(p["mu_s"]), mu_k=float(p["mu_k"]),
        tan_theta=float(np.tan(np.deg2rad(p["theta0_deg"]))),
        kn=p.get("kn", ""), class_code=cls, class_label=CLASS_LABEL[cls],
        t_end_ms=info["t_end_ms"], t_stop_ms=info["t_stop_ms"],
        stop_found=int(i < len(tt) - 1),
        contact_lost=int(lost), cap_N_hit=int(capN), cap_v_hit=int(capv),
        N_min=info["N_min"], N_max=info["N_max"], N_mean=info["N_mean"],
        N_rel_std=info["N_rel_std"], v_rel_min=info["v_rel_min"],
        v_rel_max=info["v_rel_max"], post_stop_osc=post_osc,
        xi_end_mm=info["xi_end_mm"], xi_span_mm=info["xi_span_mm"],
        I_plateau=met["I_plateau"], I_c=met["I_c"],
        cutoff_ratio=met["cutoff_ratio"],
        I_conservation=met["I_conservation"], xi_c_mm=met["xi_c_mm"],
        release_events=info.get("release_events", ""),
        constraint_residual_max=info.get("constraint_residual_max", ""),
        n0_check=info.get("n0_check", ""),
    )
    return row, tr, xc, I


def make_p(theta0, mu_s, model):
    p = dict(vm.DEFAULT, theta0_deg=float(theta0), mu_s=float(mu_s),
             mu_k=float(mu_s) / MUS_RATIO)
    p.update(model)
    return p


def stop_predicate(row):
    return row["class_code"] not in ("slide", "lift", "sprag", "chatter")


def stopped_or_chatter(row):
    return row["class_code"] in ("plateau", "nvar", "chatter")


def lift_predicate(row):
    return row["class_code"] == "lift"


def sprag_predicate(row):
    return row["class_code"] == "sprag"


def bisect(make_row, lo, hi, pred, niter=9):
    """在 [lo,hi] 上二分单调谓词；返回边界值与两端谓词。"""
    r_lo, r_hi = make_row(lo), make_row(hi)
    f_lo, f_hi = pred(r_lo), pred(r_hi)
    if f_lo == f_hi:
        return None, f_lo, f_hi
    for _ in range(niter):
        mid = 0.5 * (lo + hi)
        if pred(make_row(mid)) == f_lo:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi), f_lo, f_hi


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="./S6b_v2输出")
    ap.add_argument("--tmax-penalty", type=float, default=0.012)
    ap.add_argument("--tmax-rigid", type=float, default=0.004)
    ap.add_argument("--no-bisect", action="store_true",
                    help="跳过边界二分（快速出图用）")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    base = dict(vm.DEFAULT)

    def tmax_of(model_key):
        return args.tmax_rigid if model_key == "rigid" else args.tmax_penalty

    # ============ 1. 主平面：4 种接触口径 × theta0 × mu_s ============
    rows = []
    cache = {}
    print("=" * 76)
    print("S6b：v2 角度—摩擦系数平面（theta0 × mu_s）")
    print(f"theta0 = {THETAS}")
    print(f"mu_s   = {MUS}；mu_s/mu_k = {MUS_RATIO:.4f}")
    print("=" * 76)
    for mkey, mlabel, mpar in MODELS:
        print(f"\n--- {mlabel}（t_max={tmax_of(mkey)*1e3:.0f} ms）---")
        for th in THETAS:
            line = f"th={th:4.1f}: "
            for mus in MUS:
                p = make_p(th, mus, mpar)
                row, tr, xc, I = run_case(mkey, p, tmax_of(mkey))
                row["model_label"] = mlabel
                rows.append(row)
                cache[(mkey, th, mus)] = (row, tr, xc, I)
                line += f" {row['class_code'][:4]:>4s}"
            print(line)
    with open(os.path.join(args.outdir, "S6b_角度摩擦平面_结果表.csv"), "w",
              newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\n平面结果表：{len(rows)} 行 -> S6b_角度摩擦平面_结果表.csv")

    # 分类统计
    print("\n分类统计：")
    for mkey, mlabel, _ in MODELS:
        sub = [r for r in rows if r["model"] == mkey]
        cnt = {c: sum(r["class_code"] == c for r in sub) for c in CLASS_ORDER}
        print(f"  {mlabel:16s}: " + "  ".join(
            f"{CLASS_LABEL[c]}={cnt[c]}" for c in CLASS_ORDER if cnt[c]))

    # 基线与 S6/S4 连续性
    r0 = cache[("penalty_2e4", 20.0, 0.62)][0]
    print(f"\n[连续性] kn=2e4, th=20, mu_s=0.62: t_stop={r0['t_stop_ms']:.3f} ms, "
          f"Nrel={r0['N_rel_std']*100:.2f}%, Ic/Ip={r0['cutoff_ratio']:.3f}, "
          f"∫I/I2={r0['I_conservation']:.5f}（应命中 S6/S4 基线 1.816 ms / 0.25% / 1.002）")
    rr = cache[("rigid", 20.0, 0.62)][0]
    print(f"[连续性] rigid, th=20, mu_s=0.62: t_stop={rr['t_stop_ms']:.3f} ms, "
          f"Nrel={rr['N_rel_std']*100:.2f}%, Ic/Ip={rr['cutoff_ratio']:.3f}, "
          f"res={rr['constraint_residual_max']:.1e}（应命中 S5b 基线 1.550 ms / 5.19%）")

    # ============ 2. mu_s/mu_k 比值扫描 ============
    ratio_rows = []
    for mkey, mlabel, mpar in MODELS:
        for th in (20.0, 35.0):
            for ratio in (1.0, 1.2, 1.5, 2.0):
                p = make_p(th, 0.62, mpar)
                p["mu_k"] = 0.62 / ratio
                row, tr, xc, I = run_case(mkey, p, tmax_of(mkey))
                row["ratio"] = ratio
                ratio_rows.append(row)
                print(f"[比值] {mlabel:14s} th={th:4.1f} mu_s/mu_k={ratio:.1f}: "
                      f"{row['class_code']:8s} tstop={row['t_stop_ms']:7.3f} ms "
                      f"Nrel={row['N_rel_std']*100:6.2f}% Ic/Ip={row['cutoff_ratio']:.3f}")
    with open(os.path.join(args.outdir, "S6b_摩擦比扫描.csv"), "w",
              newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(ratio_rows[0].keys()))
        w.writeheader()
        w.writerows(ratio_rows)

    # ============ 2b. 板速对单事件停止边界的影响 ============
    # 单事件“停止”是虚线形成的必要条件；完整虚/实线相图还需离面—再接触周期。
    print("\n--- 板速对停止边界的影响（theta0=20°, kn=2e4）---")
    v_rows = []
    for v in (0.1, 0.2, 0.5):
        tmax_v = 0.030 if v <= 0.1 else (0.020 if v <= 0.2
                                         else args.tmax_penalty)
        for mus in (0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10, 0.15,
                    0.20, 0.22, 0.25):
            p = make_p(20.0, mus, dict(kn=2e4))
            p["v_board"] = v
            row, tr, xc, I = run_case("penalty_2e4", p, tmax_v)
            row["v_board"] = v
            row["tmax_ms"] = tmax_v * 1e3
            v_rows.append(row)
            print(f"  v={v:.2f} mu_s={mus:.2f}: {row['class_code']:8s} "
                  f"tstop={row['t_stop_ms']:7.3f} ms "
                  f"vmin={row['v_rel_min']:.4f}")
    with open(os.path.join(args.outdir, "S6b_板速扫描.csv"), "w",
              newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(v_rows[0].keys()))
        w.writeheader()
        w.writerows(v_rows)
    # 停止边界 v 标度估计
    print("  停止边界（首个非 slide/lift/sprag 的 mu_s）：")
    for v in (0.1, 0.2, 0.5):
        sub = [r for r in v_rows if r["v_board"] == v]
        stops = [r for r in sub if r["class_code"] in ("plateau", "nvar",
                                                       "chatter")]
        thr = min((r["mu_s"] for r in stops), default=None)
        print(f"    v={v:.2f} m/s -> mu_s*≈{thr}")

    # ============ 3. 边界：连续滑移（罚 2e4）与脱离（罚 2e5/2e6、刚性） ============
    boundary_rows = []

    def add_boundary(kind, model_key, boundary_var, boundary_value,
                     theta0_deg, mu_s, lo, hi, f_lo, f_hi, note=""):
        boundary_rows.append(dict(
            kind=kind, model=model_key, boundary_var=boundary_var,
            boundary_value=("" if boundary_value is None else boundary_value),
            theta0_deg=theta0_deg, mu_s=mu_s,
            bracket_lo=lo, bracket_hi=hi,
            pred_lo=int(f_lo), pred_hi=int(f_hi), note=note))
        bv = "—" if boundary_value is None else f"{boundary_value:.4f}"
        th_txt = (f"{theta0_deg:5.1f}" if isinstance(theta0_deg, (int, float))
                  else str(theta0_deg))
        mu_txt = f"{mu_s:4.2f}" if isinstance(mu_s, (int, float)) else str(mu_s)
        print(f"[边界] {kind:14s} {model_key:12s} "
              f"(theta0={th_txt}, mu_s={mu_txt}): "
              f"{boundary_var}*={bv}  括号({lo:.2f},{hi:.2f}) "
              f"pred {int(f_lo)}->{int(f_hi)} {note}")

    if not args.no_bisect:
        print("\n--- 连续滑移（无停止）边界：罚接触 kn=2e4（对 mu_s 二分）---")
        t_ent = max(args.tmax_penalty, 0.025)   # 临界减速：留足时间判“永不停”
        for th in [10.0, 20.0, 30.0, 40.0, 50.0]:
            def mk(mus, th=th):
                row, _, _, _ = run_case("penalty_2e4", make_p(th, mus, dict(kn=2e4)),
                                        t_ent)
                return row
            mu_b, f_lo, f_hi = bisect(mk, 0.10, 0.32, stop_predicate)
            add_boundary("entrainment", "penalty_2e4", "mu_s", mu_b, th, "",
                         0.10, 0.32, f_lo, f_hi,
                         f"pred=stop: slide(False)->stop(True); tmax={t_ent*1e3:.0f}ms")

        print("\n--- 接触脱离边界：罚接触 kn=2e5 / 2e6（对 mu_s 二分）---")
        for mkey, _, mpar in [MODELS[1], MODELS[2]]:
            for th in [20.0, 25.0, 30.0, 35.0, 40.0]:
                def mk(mus, th=th, mpar=mpar, mkey=mkey):
                    row, _, _, _ = run_case(mkey, make_p(th, mus, mpar),
                                            tmax_of(mkey))
                    return row
                # 先在 mu_s∈[0.15,1.0] 粗扫，确定 lift 区间及上边界括号
                coarse = np.linspace(0.15, 1.00, 10)
                cls = [lift_predicate(mk(m)) for m in coarse]
                if not any(cls):
                    add_boundary("liftoff", mkey, "mu_s", None, th, "",
                                 0.15, 1.0, 0, 0, "mu_s<=1.0 内无脱离")
                    continue
                i_first = cls.index(True)
                i_last = max(i for i, c in enumerate(cls) if c)
                if i_last == len(coarse) - 1:
                    add_boundary("liftoff", mkey, "mu_s", None, th, "",
                                 coarse[i_first], coarse[i_last], 1, 1,
                                 "mu_s<=1.0 仍有脱离（未闭合）")
                    continue
                lo, hi = coarse[i_first], coarse[i_last + 1]
                mu_b, f_lo, f_hi = bisect(mk, lo, hi, lift_predicate)
                add_boundary("liftoff", mkey, "mu_s", mu_b, th, "",
                             lo, hi, f_lo, f_hi,
                             "pred=lift: lift(True)->no-lift(False)")

        print("\n--- 刚性接触：脱离起始角与自锁起始角（对 theta0 二分）---")
        for mus in (0.2, 0.62, 1.0):
            def mk_th(th, mus=mus):
                row, _, _, _ = run_case("rigid", make_p(th, mus,
                                                        dict(contact="rigid")),
                                        args.tmax_rigid)
                return row
            th_b, f_lo, f_hi = bisect(mk_th, 20.0, 50.0, lift_predicate)
            add_boundary("rigid_liftoff", "rigid", "theta0_deg", th_b, "", mus,
                         20.0, 50.0, f_lo, f_hi,
                         "mu_s 固定；pred=lift。注意刚性初值的 N0 锚定口径")
            # 自锁/颤振只在 25–35° 的窄带出现，故括号取 24–35°
            th_s, f_lo, f_hi = bisect(mk_th, 24.0, 35.0, sprag_predicate)
            add_boundary("rigid_sprag_onset", "rigid", "theta0_deg", th_s, "", mus,
                         24.0, 35.0, f_lo, f_hi,
                         "mu_s 固定；pred=sprag(早停)")

    with open(os.path.join(args.outdir, "S6b_边界与判据.csv"), "w",
              newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["kind", "model", "boundary_var",
                                          "boundary_value", "theta0_deg",
                                          "mu_s", "bracket_lo", "bracket_hi",
                                          "pred_lo", "pred_hi", "note"])
        w.writeheader()
        w.writerows(boundary_rows)

    # 判据表
    crit_rows = []
    print("\n--- sprag/自锁判据（mu_s > tan(theta0) 与 mu_k > tan(theta0)）---")
    for th in [15.0, 20.0, 25.0, 30.0, 35.0, 40.0]:
        tt = np.tan(np.deg2rad(th))
        crit_rows.append(dict(theta0_deg=th, tan_theta=tt,
                              mu_s_static_eq=tt,
                              mu_s_kinetic_eq=tt * MUS_RATIO))
        print(f"  th={th:4.1f}: tan={tt:.3f} | mu_s=tan -> {tt:.3f} | "
              f"mu_k=tan -> mu_s={tt*MUS_RATIO:.3f}")
    with open(os.path.join(args.outdir, "S6b_判据表.csv"), "w",
              newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["theta0_deg", "tan_theta",
                                          "mu_s_static_eq", "mu_s_kinetic_eq"])
        w.writeheader()
        w.writerows(crit_rows)

    # ============ 4. 数值校验：双容差 + tmax 敏感性 ============
    print("\n--- 数值校验：双容差复算 ---")
    checks = [
        ("penalty_2e4", dict(kn=2e4), 20.0, 0.62),
        ("penalty_2e4", dict(kn=2e4), 20.0, 0.22),
        ("penalty_2e6", dict(kn=2e6), 30.0, 0.50),
        ("penalty_2e6", dict(kn=2e6), 20.0, 0.62),
        ("rigid", dict(contact="rigid"), 20.0, 0.62),
        ("rigid", dict(contact="rigid"), 40.0, 0.62),
    ]
    check_rows = []
    for mkey, mpar, th, mus in checks:
        p = make_p(th, mus, mpar)
        r1, _, _, _ = run_case(mkey, p, tmax_of(mkey), rtol=1e-8,
                               max_step=1e-5)
        r2, _, _, _ = run_case(mkey, p, tmax_of(mkey), rtol=1e-10,
                               max_step=1e-6)
        dt = abs(r1["t_stop_ms"] - r2["t_stop_ms"]) / max(r1["t_stop_ms"], 1e-9)
        dxi = abs(r1["xi_end_mm"] - r2["xi_end_mm"]) / max(abs(r1["xi_end_mm"]), 1e-9)
        dn = abs(r1["N_rel_std"] - r2["N_rel_std"]) / max(abs(r1["N_rel_std"]), 1e-9)
        same = r1["class_code"] == r2["class_code"]
        row = dict(model=mkey, theta0_deg=th, mu_s=mus,
                   class_rtol8=r1["class_code"], class_rtol10=r2["class_code"],
                   class_same=int(same), t_stop_rtol8=r1["t_stop_ms"],
                   t_stop_rtol10=r2["t_stop_ms"], rel_dt=dt,
                   xi_rtol8=r1["xi_end_mm"], xi_rtol10=r2["xi_end_mm"],
                   rel_dxi=dxi, Nrel_rtol8=r1["N_rel_std"],
                   Nrel_rtol10=r2["N_rel_std"], rel_dNrel=dn)
        check_rows.append(row)
        print(f"  {mkey:12s} th={th:4.1f} mu={mus:.2f}: {r1['class_code']:8s}/"
              f"{r2['class_code']:8s} Δt/t={dt:.1e} Δξ/ξ={dxi:.1e} "
              f"ΔNrel={dn:.1e} {'OK' if same else 'CLASS-DIFF'}")

    print("\n--- 数值校验：低摩擦无停止工况的 tmax 敏感性 ---")
    for mus in (0.19, 0.20, 0.21, 0.22):
        p = make_p(20.0, mus, dict(kn=2e4))
        r12, _, _, _ = run_case("penalty_2e4", p, 0.012)
        r30, _, _, _ = run_case("penalty_2e4", p, 0.030)
        note = ("12ms==" if r12["class_code"] == r30["class_code"] else "12ms!=") + \
               f"{r12['class_code']}/{r30['class_code']}"
        row = dict(model="penalty_2e4", theta0_deg=20.0, mu_s=mus,
                   tmax12_class=r12["class_code"],
                   tmax30_class=r30["class_code"],
                   tmax12_tstop=r12["t_stop_ms"],
                   tmax30_tstop=r30["t_stop_ms"],
                   tmax12_vmin=r12["v_rel_min"],
                   tmax30_vmin=r30["v_rel_min"], note=note)
        check_rows.append(row)
        print(f"  mu_s={mus:.2f}: tmax12={r12['class_code']:8s} "
              f"tstop={r12['t_stop_ms']:7.3f}  tmax30={r30['class_code']:8s} "
              f"tstop={r30['t_stop_ms']:7.3f}  vmin={r30['v_rel_min']:.4f}")
    with open(os.path.join(args.outdir, "S6b_数值校验.csv"), "w",
              newline="", encoding="utf-8-sig") as f:
        field = ["model", "theta0_deg", "mu_s", "class_rtol8", "class_rtol10",
                 "class_same", "t_stop_rtol8", "t_stop_rtol10", "rel_dt",
                 "xi_rtol8", "xi_rtol10", "rel_dxi", "Nrel_rtol8",
                 "Nrel_rtol10", "rel_dNrel", "tmax12_class", "tmax30_class",
                 "tmax12_tstop", "tmax30_tstop", "tmax12_vmin", "tmax30_vmin",
                 "note"]
        w = csv.DictWriter(f, fieldnames=field, extrasaction="ignore")
        w.writeheader()
        w.writerows(check_rows)

    # ============ 5. 典型轨迹 ============
    typical = [
        ("A_平台+截断(罚2e4,th20,mu0.62)", "penalty_2e4",
         dict(kn=2e4), 20.0, 0.62),
        ("B_连续滑移(罚2e4,th20,mu0.19)", "penalty_2e4", dict(kn=2e4), 20.0, 0.19),
        ("C_脱离(罚2e6,th30,mu0.50)", "penalty_2e6", dict(kn=2e6), 30.0, 0.50),
        ("D_自锁(刚性,th34,mu0.62)", "rigid", dict(contact="rigid"), 34.0, 0.62),
    ]
    traj_rows, prof_rows = [], []
    fig, ax = plt.subplots(3, 1, figsize=(10, 9.5), sharex=False)
    for label, mkey, mpar, th, mus in typical:
        row, tr, xc, I = run_case(mkey, make_p(th, mus, mpar), tmax_of(mkey))
        ax[0].plot(tr["tt"] * 1e3, tr["N"], label=f"{label}  [{row['class_code']}]")
        ax[1].plot(tr["tt"] * 1e3, tr["v_rel"], label=label)
        if xc.size:
            ax[2].plot(xc * 1e3, I, label=label)
        for k in range(0, len(tr["tt"]), max(1, len(tr["tt"]) // 400)):
            traj_rows.append(dict(case=label, t_ms=tr["tt"][k] * 1e3,
                                  N=tr["N"][k], v_rel=tr["v_rel"][k],
                                  xi_mm=tr["xi"][k] * 1e3,
                                  class_code=row["class_code"]))
        for k in range(len(xc)):
            prof_rows.append(dict(case=label, xi_mm=xc[k] * 1e3, I=I[k],
                                  class_code=row["class_code"]))
        print(f"[典型] {label}: {row['class_code']:8s} "
              f"tstop={row['t_stop_ms']:.3f} ms, tend={row['t_end_ms']:.3f} ms, "
              f"Nrel={row['N_rel_std']*100:.1f}%, Ic/Ip={row['cutoff_ratio']:.3f}, "
              f"ξ={row['xi_end_mm']:.3f} mm")
    ax[0].set_ylabel("N / N")
    ax[1].set_ylabel("v_rel / (m/s)")
    ax[1].axhline(0.0, color="k", lw=0.8, ls=":")
    ax[2].set_xlabel("ξ / mm")
    ax[2].set_ylabel("I(ξ) / (J/m)")
    ax[0].set_ylim(0.0, 8.0)      # 自锁/脱离工况的瞬态尖峰被截断显示
    ax[1].set_ylim(-0.12, 0.62)
    ax[2].set_ylim(0.0, 2.6)
    for a, t in zip(ax, ["(a) 法向力 N(t)", "(b) 相对速度 v_rel(t)",
                         "(c) 单次划痕强度 I(ξ)"]):
        a.set_title(t)
        a.grid(alpha=0.3)
        a.legend(fontsize=7)
    ax[0].set_xlabel("t / ms")
    ax[1].set_xlabel("t / ms")
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "S6b_典型轨迹.png"), dpi=150)
    plt.close(fig)
    with open(os.path.join(args.outdir, "S6b_典型轨迹.csv"), "w",
              newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["case", "t_ms", "N", "v_rel",
                                          "xi_mm", "class_code"])
        w.writeheader()
        w.writerows(traj_rows)
    with open(os.path.join(args.outdir, "S6b_典型剖面.csv"), "w",
              newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["case", "xi_mm", "I", "class_code"])
        w.writeheader()
        w.writerows(prof_rows)

    # ============ 6. 相图（4 联） ============
    cmap = ListedColormap([CLASS_COLOR[c] for c in CLASS_ORDER])
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    th_cr = np.linspace(3.0, 60.0, 300)
    for ax, (mkey, mlabel, _) in zip(axes.ravel(), MODELS):
        sub = [r for r in rows if r["model"] == mkey]
        X = np.array([r["theta0_deg"] for r in sub])
        Y = np.array([r["mu_s"] for r in sub])
        idx = np.array([CLASS_ORDER.index(r["class_code"]) for r in sub])
        ax.scatter(X, Y, c=idx, cmap=cmap, vmin=-0.5, vmax=len(CLASS_ORDER) - 0.5,
                   marker="s", s=110, edgecolors="k", linewidths=0.4)
        ax.plot(th_cr, np.tan(np.deg2rad(th_cr)), "k-.", lw=1.2,
                label=r"$\mu_s=\tan\theta_0$（静摩擦自锁）")
        ax.plot(th_cr, MUS_RATIO * np.tan(np.deg2rad(th_cr)), "k--", lw=1.2,
                label=r"$\mu_k=\tan\theta_0$（动摩擦判据 $\mu_s=1.216\tan\theta_0$）")
        ax.set_xlim(3, 60)
        ax.set_ylim(0.10, 1.06)
        ax.set_xlabel(r"$\theta_0$ / (°)")
        ax.set_ylabel(r"$\mu_s$")
        ax.set_title(mlabel)
        ax.grid(alpha=0.25)
    handles = [Patch(facecolor=CLASS_COLOR[c], edgecolor="k",
                     label=CLASS_LABEL[c]) for c in CLASS_ORDER]
    handles += [plt.Line2D([], [], color="k", ls="-.", label=r"$\mu_s=\tan\theta_0$"),
                plt.Line2D([], [], color="k", ls="--",
                           label=r"$\mu_k=\tan\theta_0$（$\mu_s=1.216\tan\theta_0$）")]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=8.5,
               frameon=False)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(os.path.join(args.outdir, "S6b_角度摩擦相图.png"), dpi=150)
    plt.close(fig)

    # ============ 7. 罚模型指标图 ============
    ths = np.array(THETAS)
    mus_arr = np.array(MUS)
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))
    # (a) kn=2e4 Nrel, (b) cutoff, (c) t_stop
    for a, key, title in zip(
            [ax[0, 0], ax[0, 1], ax[1, 0]],
            ["N_rel_std", "cutoff_ratio", "t_stop_ms"],
            ["(a) 罚 k_n=2e4：N 相对涨落",
             "(b) 罚 k_n=2e4：I_c/I_plateau", "(c) 罚 k_n=2e4：t_stop / ms"]):
        Z = np.full((len(mus_arr), len(ths)), np.nan)
        for r in rows:
            if r["model"] != "penalty_2e4" or np.isnan(r[key]):
                continue
            if r["class_code"] in ("slide", "lift", "sprag"):
                continue          # 非停止工况没有 t_stop/截断比，留白
            i = int(np.where(ths == r["theta0_deg"])[0][0])
            j = int(np.where(mus_arr == r["mu_s"])[0][0])
            Z[j, i] = r[key]
        if key == "N_rel_std":
            Z = Z * 100
        pc = a.pcolormesh(ths, mus_arr, Z, shading="nearest", cmap="viridis")
        fig.colorbar(pc, ax=a)
        a.set_xlabel(r"$\theta_0$ / (°)")
        a.set_ylabel(r"$\mu_s$")
        a.set_title(title)
    # (d) theta 切面（mu_s=0.62）
    for mkey, mlabel, mpar in MODELS:
        th_cut, nrel_cut = [], []
        for th in np.arange(5.0, 46.0, 2.5):
            row = run_case(mkey, make_p(th, 0.62, mpar), tmax_of(mkey))[0]
            if row["class_code"] not in ("plateau", "nvar", "chatter"):
                continue
            th_cut.append(th)
            nrel_cut.append(row["N_rel_std"] * 100)
        ax[1, 1].plot(th_cut, nrel_cut, "o-", ms=3, label=mlabel)
    ax[1, 1].axhline(10.0, color="k", ls=":", lw=1)
    ax[1, 1].set_xlabel(r"$\theta_0$ / (°)")
    ax[1, 1].set_ylabel("N 相对涨落 / %")
    ax[1, 1].set_title(r"(d) 角度切面（$\mu_s=0.62$）")
    ax[1, 1].grid(alpha=0.3)
    ax[1, 1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "S6b_罚模型指标.png"), dpi=150)
    plt.close(fig)

    # ============ 8. 边界与判据图 ============
    fig, ax = plt.subplots(1, 3, figsize=(19, 5.2))
    th_cr = np.linspace(3.0, 55.0, 300)
    ax[0].plot(th_cr, np.tan(np.deg2rad(th_cr)), "k-.", lw=1.2,
               label=r"$\mu_s=\tan\theta_0$")
    ax[0].plot(th_cr, MUS_RATIO * np.tan(np.deg2rad(th_cr)), "k--", lw=1.2,
               label=r"$\mu_k=\tan\theta_0$")
    styles = {"entrainment": ("s", "#1565c0", "连续滑移边界（罚 2e4）"),
              "liftoff": ("o", "#c62828", "脱离边界（罚 2e5/2e6）")}
    for kind, (mk, col, lab) in styles.items():
        for mkey, mlabel, _ in MODELS:
            sub = [b for b in boundary_rows if b["kind"] == kind
                   and b["model"] == mkey and b["boundary_value"] != ""]
            if sub:
                ax[0].plot([float(b["theta0_deg"]) for b in sub],
                           [float(b["boundary_value"]) for b in sub], mk + "-",
                           color=col, alpha=0.55 if kind == "liftoff" else 1.0,
                           label=f"{lab} {mkey}")
    ax[0].set_xlabel(r"$\theta_0$ / (°)")
    ax[0].set_ylabel(r"$\mu_s$")
    ax[0].set_xlim(3, 55)
    ax[0].set_ylim(0.1, 1.1)
    ax[0].set_title("(a) 边界 vs 判据曲线")
    ax[0].grid(alpha=0.3)
    ax[0].legend(fontsize=7)
    # 刚性边界：横轴 theta0，纵轴 mu_s
    sub = [r for r in rows if r["model"] == "rigid"]
    ax[1].scatter([r["theta0_deg"] for r in sub], [r["mu_s"] for r in sub],
                  c=[CLASS_COLOR[r["class_code"]] for r in sub],
                  marker="s", s=100, edgecolors="k", linewidths=0.4)
    for b in boundary_rows:
        if b["kind"] == "rigid_liftoff" and b["boundary_value"] != "":
            ax[1].axvline(float(b["boundary_value"]), color="#c62828", ls="--",
                          lw=1)
        if b["kind"] == "rigid_sprag_onset" and b["boundary_value"] != "":
            ax[1].axvline(float(b["boundary_value"]), color="#6a1b9a", ls=":",
                          lw=1)
    ax[1].set_xlim(3, 60)
    ax[1].set_ylim(0.10, 1.06)
    ax[1].set_xlabel(r"$\theta_0$ / (°)")
    ax[1].set_ylabel(r"$\mu_s$")
    ax[1].set_title("(b) 刚性接触相图\n红虚线=脱离起始角 紫点线=自锁起始角")
    ax[1].grid(alpha=0.25)
    handles = [Patch(facecolor=CLASS_COLOR[c], edgecolor="k",
                     label=CLASS_LABEL[c]) for c in CLASS_ORDER]
    ax[1].legend(handles=handles, fontsize=7, loc="upper left")
    # (c) 板速—摩擦停止边界（单事件）
    for v in (0.1, 0.2, 0.5):
        sub = [r for r in v_rows if r["v_board"] == v]
        ax[2].scatter([r["mu_s"] for r in sub], [v] * len(sub),
                      c=[CLASS_COLOR[r["class_code"]] for r in sub],
                      marker="s", s=90, edgecolors="k", linewidths=0.4)
    for v in (0.1, 0.2, 0.5):
        sub = [r for r in v_rows if r["v_board"] == v]
        stops = [r for r in sub if r["class_code"] in ("plateau", "nvar",
                                                       "chatter")]
        if stops:
            thr = min(r["mu_s"] for r in stops)
            ax[2].plot([thr, thr], [v * 0.9, v * 1.1], "k-", lw=1.5)
    ax[2].set_xlabel(r"$\mu_s$")
    ax[2].set_ylabel(r"$v$ / (m/s)")
    ax[2].set_xlim(0, 0.22)
    ax[2].set_ylim(0.05, 0.6)
    ax[2].set_title("(c) 板速对停止边界的影响\n（单事件；竖线=首个停止）")
    ax[2].grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "S6b_边界与判据.png"), dpi=150)
    plt.close(fig)

    # ============ 9. 验收小结 ============
    print("\n" + "=" * 76)
    print("S6b 验收小结")
    print("=" * 76)
    print(f"① 平面覆盖：{len(THETAS)} theta0 × {len(MUS)} mu_s × {len(MODELS)} 接触口径")
    print("② 罚接触 2e4 连续滑移边界（mu_s*；括号内为 mu_k* = mu_s*/1.2157）：")
    for b in boundary_rows:
        if b["kind"] == "entrainment" and b["boundary_value"] != "":
            mb = float(b["boundary_value"])
            print(f"   theta0={float(b['theta0_deg']):4.1f}° -> mu_s*={mb:.3f} "
                  f"(mu_k*={mb/MUS_RATIO:.3f})")
    print("③ 罚接触脱离边界（停止前先脱离的最大 mu_s；空=未闭合/无脱离）：")
    for b in boundary_rows:
        if b["kind"] == "liftoff":
            print(f"   {b['model']:12s} theta0={float(b['theta0_deg']):4.1f}° -> "
                  f"mu_s*={b['boundary_value'] if b['boundary_value'] == '' else round(float(b['boundary_value']), 3)}"
                  f" {b['note']}")
    print("③b 板速—停止边界（theta0=20°, kn=2e4）：")
    for v in (0.1, 0.2, 0.5):
        sub = [r for r in v_rows if r["v_board"] == v]
        stops = [r for r in sub if r["class_code"] in ("plateau", "nvar",
                                                       "chatter")]
        thr = min((r["mu_s"] for r in stops), default=None)
        vmin_slide = [r["v_rel_min"] for r in sub if r["class_code"] == "slide"]
        print(f"   v={v:.2f} m/s -> mu_s*≈{thr}"
              + (f"（其下方 mu_s 连续滑移，v_rel_min≈{min(vmin_slide):.3f} m/s）"
                 if vmin_slide else ""))
    print("④ 刚性接触分界：")
    for b in boundary_rows:
        if b["kind"] in ("rigid_liftoff", "rigid_sprag_onset"):
            bv = b["boundary_value"]
            print(f"   {b['kind']:18s} mu_s={float(b['mu_s']):.2f} -> "
                  f"theta*={bv if bv == '' else round(float(bv), 2)}° {b['note']}")
    cls_counts = {}
    for r in rows:
        cls_counts.setdefault(r["model"], {}).setdefault(r["class_code"], 0)
        cls_counts[r["model"]][r["class_code"]] += 1
    for mkey, mlabel, _ in MODELS:
        c = cls_counts.get(mkey, {})
        print(f"⑤ {mlabel:16s}: " + "  ".join(
            f"{CLASS_LABEL[k]}={c.get(k, 0)}" for k in CLASS_ORDER))
    print(f"\n输出目录：{os.path.abspath(args.outdir)}")


if __name__ == "__main__":
    main()
