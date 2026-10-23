# -*- coding: utf-8 -*-
"""
仿真 S1：k_c 刚性极限收敛扫描
==============================================
对应研究记录《新参考文献与N下降段方向对照笔记》第 2 节"刚性极限收敛研究"（优先级最高）。

背景
----
第一阶段模型用罚弹簧 N = k_c·δ 代替刚性接触（δ = L·cosθ − h 为压入量），
k_c = 2×10^4 N/m 是人为选定的标定参数。本扫描把 k_c 在 1e3 ~ 1e7 N/m 区间
系统变化，重复计算单周期 N 下降段，检验输出量随 k_c 的行为：

  t*      分离时间（δ 归零事件时刻）
  t_flip  摩擦翻转时刻（尖端水平速度赶上板速、v_rel 换向，两阶段图像的分界）
  I_2nd   本段笔痕强度  μ∫N|v_rel|dt
  ℓ       痕迹长 ∫|v_rel|dt

判读（对照笔记）：
  · 规则收敛/标度  → 罚模型自洽，k_c 只须"足够大"，标定路线可行；
  · 振荡、跳跃或反常发散 → 刚性极限下 Painlevé 悖论在起作用，
    分离须改用法向-切向 LCP / 碰撞律描述（第二阶段建模入口）。

两种模式
--------
MODE='fixed_N0'     固定交接压力 N0（默认）。交接量是压力而非压入量。
                    注意 δ0 = N0/k_c 随 k_c 增大而缩小，θ0 → θ*，
                    段的几何区间退化，t*、I_2nd 预期按幂律趋于 0——
                    本模式判据是"标度是否规则"而非"是否趋于常数"。
MODE='fixed_delta0' 固定初始压入量 δ0 = 2×10^-4 m（名义点 N0=4 N @ k_c=2e4），
                    N0 = k_c·δ0 随 k_c 增大，段几何保持不变，
                    线性化弹簧-振子给出基准 t* ∝ k_c^(-1/2)，可对照斜率。

数值可信度
----------
每个 (k_c, N0, ω0) 点用两套容差（rtol=1e-9 / 1e-12）各积分一次，
二者之差记为数值误差，防止把"积分不准"误判为"物理发散"。
k_c 越大接触频率越高，max_step 随特征频率自适应收紧。

输出（--outdir，默认 ./S1_kc扫描输出/）
  S1_kc扫描_结果表.csv    全输出量 + 收敛诊断列
  S1_kc扫描_收敛曲线.png  t*、t_flip、I_2nd、ℓ 对 k_c（名义工况，双容差误差棒）
  S1_kc扫描_曲线族.png    N(t)、θ(t) 曲线族（名义工况，按 k_c 着色）
  S1_kc扫描_诊断.png      对数斜率与相对变化（收敛/振荡判读）

运行
----
python 仿真S1_kc刚性极限收敛扫描.py                 # 默认 fixed_N0 模式
python 仿真S1_kc刚性极限收敛扫描.py --mode fixed_delta0
python 仿真S1_kc刚性极限收敛扫描.py --outdir ../输出 --mu 0.40

依赖：numpy、scipy、matplotlib。
注意：本脚本参数与 code/粉笔_N下降段_完整模型.py 保持一致，主模型更新时需同步修改。
"""
from __future__ import annotations

import argparse
import csv
import os

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- 中文字体：优先注册系统 Noto CJK，再按常见本机字体回退 ----
from matplotlib import font_manager as _fm

for _p in ("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
           "/usr/share/fonts/opentype/noto/NotoSansCJK-Bold.ttc"):
    if os.path.exists(_p):
        try:
            _fm.fontManager.addfont(_p)
        except Exception:
            pass
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "PingFang SC",
                                   "Noto Sans CJK SC", "Noto Sans CJK JP",
                                   "Arial Unicode MS", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

# ---------------- 物理参数（与 code/粉笔_N下降段_完整模型.py 一致） ----------------
L   = 0.08            # 粉笔长 m
h   = 0.075           # 轴 O 到板面距离 m（h < L）
RADIUS = 3.5e-3       # 粉笔半径 m
RHO = 2300.0          # 密度 kg/m^3（石膏）
GRAV = 9.8
M   = RHO * np.pi * RADIUS**2 * L
I_O = M * L**2 / 3.0
MU_DEFAULT = 0.35     # 动摩擦系数（E1/E2 标定后替换）
V_BOARD = 0.5         # 板速 m/s
EPS = 0.02            # 摩擦方向正则化速度尺度 m/s

TH_STAR = np.arccos(h / L)
ARM_RATIO_STAR = MU_DEFAULT * h / (L * np.sin(TH_STAR))   # 分离点力臂比 ≈ 0.94（近临界）

# ---------------- 扫描网格（按需修改） ----------------
KC_LIST    = np.logspace(3.0, 7.0, 9)     # 1e3 … 1e7
N0_LIST    = [1.0, 2.0, 4.0, 8.0]         # 交接压力网格 N（fixed_N0 模式）
OM_LIST    = [0.0, 3.0, 6.0]              # 交接角速度网格 rad/s
NOMINAL    = (4.0, 0.0)                   # 名义工况（作图与判读用）
KC_NOM     = 2.0e4
DELTA0_REF = NOMINAL[0] / KC_NOM          # fixed_delta0 模式参考压入量 = 2e-4 m
T_MAX      = 0.05                         # 段积分上限 s
RTOL_PAIR  = (1e-9, 1e-12)                # 双容差复算

# numpy<2 兼容
_trapz = getattr(np, "trapezoid", None) or np.trapz


# ============================== 模型内核（与主模型一致） ==============================
def delta(th):
    """压入量 δ(θ) = L·cosθ − h（>0 表示接触）。"""
    return L * np.cos(th) - h


def contact_omega(kc, th):
    """接触弹簧-转动振子的特征角频率，用于自适应 max_step。"""
    return float(np.sqrt(max(kc * (L * np.sin(th)) ** 2 / I_O, 1.0)))


def make_rhs(kc, mu):
    """右端项：θ̈ = −τ/I_O（θ 增大为顺时针）。"""
    def rhs(t, y):
        th, w = y
        d = delta(th)
        if d <= 0.0:
            return [w, -(M * GRAV * L / 2.0) * np.sin(th) / I_O]
        Nn = kc * d
        v_rel = V_BOARD - L * np.cos(th) * w
        tau = (M * GRAV * L / 2.0) * np.sin(th) \
            - Nn * (L * np.sin(th) + mu * h * np.tanh(v_rel / EPS))
        return [w, -tau / I_O]
    return rhs


def _sep_event(t, y):
    return delta(y[0])


_sep_event.terminal = True
_sep_event.direction = -1


def initial_state(kc, N0):
    """由 N0 = k_c·δ0 反解 θ0；几何不可达时返回 None。"""
    d0 = N0 / kc
    cos0 = (h + d0) / L
    if cos0 >= 1.0:
        return None, d0, (f"几何不可达：δ0={d0*1e3:.2f} mm "
                          f"≥ L−h={(L-h)*1e3:.1f} mm")
    return np.arccos(np.clip(cos0, -1.0, 1.0)), d0, ""


# ============================== 单段积分 ==============================
def run_segment(kc, N0, omega0=0.0, mu=MU_DEFAULT, rtol=1e-10,
                t_max=T_MAX, keep_traj=False):
    """积分单段 [t0, t*]，返回 (info字典, traj字典或None)。"""
    th0, d0, note = initial_state(kc, N0)
    base = dict(kc=kc, N0=N0, omega0=omega0, delta0_mm=d0 * 1e3,
                valid=th0 is not None, note=note)
    if th0 is None:
        base.update(theta0_deg=np.nan, arm_ratio0=np.nan, t_star_ms=np.nan,
                    t_flip_ms=np.nan, flip=False, I_2nd=np.nan, I_start=np.nan,
                    I_slide=np.nan, ell_mm=np.nan, omega_star=np.nan,
                    xdot_max=np.nan, separated=False)
        return base, None

    sol = solve_ivp(make_rhs(kc, mu), (0.0, t_max), [th0, omega0],
                    method="DOP853", rtol=rtol, atol=1e-12,
                    events=_sep_event, dense_output=True,
                    max_step=min(5e-5, 0.05 * 2 * np.pi / contact_omega(kc, th0)))
    separated = len(sol.t_events[0]) > 0
    t_end = float(sol.t_events[0][0]) if separated else t_max

    tt = np.linspace(0.0, t_end, 4001)
    th, w = sol.sol(tt)
    Nn = kc * np.clip(delta(th), 0.0, None)
    v_rel = V_BOARD - L * np.cos(th) * w

    # ---- 摩擦翻转检测：v_rel 首次由正变负（线性插值细化时刻） ----
    if v_rel[0] > 0.0:
        idx = np.nonzero(v_rel < 0.0)[0]
        flip = idx.size > 0
        if flip:
            i = int(idx[0])
            t_flip = float(tt[i - 1] + (tt[i] - tt[i - 1]) * v_rel[i - 1]
                           / (v_rel[i - 1] - v_rel[i])) if i > 0 else float(tt[0])
        else:
            t_flip = np.nan
    else:
        flip, t_flip = True, 0.0   # 交接瞬间相对速度已为负（交接ω0过大）

    # ---- 积分量 ----
    I2 = mu * float(_trapz(Nn * np.abs(v_rel), tt))
    ell = float(_trapz(np.abs(v_rel), tt))
    if flip:
        m0 = tt <= t_flip
        I_start = mu * float(_trapz(Nn[m0] * np.abs(v_rel[m0]), tt[m0]))
    else:
        I_start = I2
    I_slide = I2 - I_start
    xdot = L * np.cos(th) * np.abs(w)

    base.update(theta0_deg=float(np.degrees(th0)),
                arm_ratio0=float(mu * h / (L * np.sin(th0))),
                t_star_ms=t_end * 1e3,
                t_flip_ms=(t_flip * 1e3 if flip else np.nan),
                flip=flip, I_2nd=I2, I_start=I_start, I_slide=I_slide,
                ell_mm=ell * 1e3, omega_star=float(w[-1]),
                xdot_max=float(xdot.max()), separated=separated)
    if not separated:
        base["note"] = (note + "；t_max 内未分离").strip("；")
    traj = dict(tt=tt, th=th, w=w, Nn=Nn, v_rel=v_rel) if keep_traj else None
    return base, traj


def run_dual(kc, N0, omega0, mu):
    """双容差复算：返回 (宽松容差结果含数值误差列, 轨迹)。"""
    a, traj = run_segment(kc, N0, omega0, mu, rtol=RTOL_PAIR[0], keep_traj=True)
    b, _ = run_segment(kc, N0, omega0, mu, rtol=RTOL_PAIR[1])
    if np.isfinite(a["t_star_ms"]) and np.isfinite(b["t_star_ms"]):
        a["num_err_t_ms"] = abs(a["t_star_ms"] - b["t_star_ms"])
    else:
        a["num_err_t_ms"] = np.nan
    if np.isfinite(a["I_2nd"]) and np.isfinite(b["I_2nd"]):
        a["num_err_I"] = abs(a["I_2nd"] - b["I_2nd"])
    else:
        a["num_err_I"] = np.nan
    return a, traj


# ============================== 收敛诊断 ==============================
def add_diagnostics(rows, mode):
    """组内（按 kc 排序）补收敛诊断列。
    fixed_N0 模式按 (N0, omega0) 分组；fixed_delta0 模式 N0=k_c·δ0 随 k_c 变化，
    改按 omega0 分组（同组即同一条 k_c 扫描线）。"""
    groups = {}
    for r in rows:
        gkey = (r["N0"], r["omega0"]) if mode == "fixed_N0" else (0.0, r["omega0"])
        groups.setdefault(gkey, []).append(r)
    for g in groups.values():
        g.sort(key=lambda r: r["kc"])
        for key, out in (("t_star_ms", "slope_t"), ("I_2nd", "slope_I")):
            ys = np.array([r[key] for r in g], dtype=float)
            xs = np.array([r["kc"] for r in g], dtype=float)
            slopes = [np.nan]
            for i in range(1, len(g)):
                if np.all(np.isfinite([ys[i - 1], ys[i]])) and ys[i - 1] != 0.0:
                    s = ((np.log10(abs(ys[i])) - np.log10(abs(ys[i - 1])))
                         / (np.log10(xs[i]) - np.log10(xs[i - 1])))
                else:
                    s = np.nan
                slopes.append(s)
            for r, s in zip(g, slopes):
                r[out] = s
        # t* 的相邻相对变化与振荡计数（变向次数）
        ts = np.array([r["t_star_ms"] for r in g], dtype=float)
        rel, sgn, osc = [np.nan], 0, 0
        for i in range(1, len(g)):
            if np.isfinite(ts[i - 1]) and np.isfinite(ts[i]) and ts[i - 1] != 0.0:
                rel.append(abs((ts[i] - ts[i - 1]) / ts[i - 1]))
                s = np.sign(ts[i] - ts[i - 1])
                if sgn != 0 and s != 0 and s != sgn:
                    osc += 1
                if s != 0:
                    sgn = s
            else:
                rel.append(np.nan)
        for r, v in zip(g, rel):
            r["rel_chg_t"] = v
        g[0]["osc_count"] = osc   # 组的振荡计数记在组首行


# ============================== 输出 ==============================
FIELDNAMES = ["mode", "kc", "N0", "omega0", "delta0_mm", "theta0_deg",
              "arm_ratio0", "t_star_ms", "t_flip_ms", "flip", "I_2nd",
              "I_start", "I_slide", "ell_mm", "omega_star", "xdot_max",
              "separated", "valid", "num_err_t_ms", "num_err_I",
              "slope_t", "slope_I", "rel_chg_t", "osc_count", "note"]


def save_csv(rows, path):
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=FIELDNAMES, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if (isinstance(v, float) and not np.isfinite(v)) else v)
                        for k, v in r.items() if k in FIELDNAMES})


def plot_convergence(nom, outdir):
    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    ks = np.array([r["kc"] for r in nom], dtype=float)
    specs = [("t_star_ms", "t* / ms", "分离时间", True),
             ("t_flip_ms", "t_flip / ms", "摩擦翻转时刻", False),
             ("I_2nd", "I_2nd / J", "本段笔痕强度", False),
             ("ell_mm", "ℓ / mm", "痕迹长", False)]
    for a, (key, ylab, tit, use_err) in zip(ax.flat, specs):
        ys = np.array([r[key] for r in nom], dtype=float)
        es = (np.array([r["num_err_t_ms"] for r in nom], dtype=float)
              if use_err else None)
        a.errorbar(ks, ys, yerr=es, fmt="o-", lw=1.6, ms=5, capsize=3)
        a.set_xscale("log"); a.set_yscale("log")
        a.set_xlabel("k_c / (N/m)"); a.set_ylabel(ylab); a.set_title(tit)
        a.grid(alpha=0.3, which="both")
    fig.suptitle(f"S1 名义工况 (N0={NOMINAL[0]} N, ω0={NOMINAL[1]} rad/s)：输出量对 k_c 的收敛性")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "S1_kc扫描_收敛曲线.png"), dpi=150)
    plt.close(fig)


def plot_family(trajs, outdir):
    ks_sorted = sorted(trajs)
    if not ks_sorted:
        return
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.6))
    norm = plt.Normalize(np.log10(min(ks_sorted)), np.log10(max(ks_sorted)))
    cmap, sm = plt.cm.viridis, plt.cm.ScalarMappable(
        norm=plt.Normalize(np.log10(min(ks_sorted)), np.log10(max(ks_sorted))),
        cmap=plt.cm.viridis)
    sm.set_array([])
    for kc in ks_sorted:
        tr = trajs[kc]
        c = cmap(norm(np.log10(kc)))
        ax[0].plot(tr["tt"] * 1e3, tr["Nn"], color=c, lw=1.5)
        ax[1].plot(tr["tt"] * 1e3, np.degrees(tr["th"]), color=c, lw=1.5)
    ax[0].set_xlabel("t / ms"); ax[0].set_ylabel("N / N")
    ax[0].set_title("N(t) 曲线族（按 k_c 着色）")
    ax[1].set_xlabel("t / ms"); ax[1].set_ylabel("θ / °")
    ax[1].set_title("θ(t) 曲线族")
    for a in ax:
        a.grid(alpha=0.3)
    fig.colorbar(sm, ax=ax, label="log10(k_c / (N/m))")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "S1_kc扫描_曲线族.png"), dpi=150)
    plt.close(fig)


def plot_diagnostics(nom, outdir):
    ks = np.array([r["kc"] for r in nom][1:], dtype=float)
    s_t = np.array([r["slope_t"] for r in nom][1:], dtype=float)
    s_I = np.array([r["slope_I"] for r in nom][1:], dtype=float)
    rel = np.array([r["rel_chg_t"] for r in nom][1:], dtype=float)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].plot(ks, s_t, "o-", label="t* 对数斜率")
    ax[0].plot(ks, s_I, "s-", label="I_2nd 对数斜率")
    ax[0].axhline(-0.5, color="gray", ls="--", lw=1,
                  label="基准 −1/2（线性弹簧-振子）")
    ax[0].set_xscale("log")
    ax[0].set_xlabel("k_c / (N/m)"); ax[0].set_ylabel("d log y / d log k_c")
    ax[0].set_title("对数斜率"); ax[0].legend(); ax[0].grid(alpha=0.3)
    ax[1].plot(ks, rel, "o-")
    ax[1].axhline(0.01, color="gray", ls="--", lw=1, label="1% 判线")
    ax[1].set_xscale("log"); ax[1].set_yscale("log")
    ax[1].set_xlabel("k_c / (N/m)"); ax[1].set_ylabel("|Δt*| / t*")
    ax[1].set_title("相邻点相对变化"); ax[1].legend(); ax[1].grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "S1_kc扫描_诊断.png"), dpi=150)
    plt.close(fig)


# ============================== 汇总与判读提示 ==============================
def print_summary(rows, mode, mu, outdir):
    print("=" * 74)
    print("S1：k_c 刚性极限收敛扫描")
    print(f"模式 = {mode}   μ = {mu}   板速 v = {V_BOARD} m/s   ε = {EPS} m/s")
    print(f"棒质量 m = {M*1e3:.2f} g,  I_O = {I_O:.3e} kg·m²")
    print(f"分离角 θ* = {np.degrees(TH_STAR):.2f}°   "
          f"分离点力臂比 μh/(L·sinθ*) = {mu*h/(L*np.sin(TH_STAR)):.3f}  "
          f"(=1 为临界，越接近 1 越要盯紧)")
    print("=" * 74)
    nom = [r for r in rows if r.get("is_nominal")]
    nom.sort(key=lambda r: r["kc"])
    hdr = (f"{'kc':>9} {'δ0/mm':>7} {'θ0/°':>7} {'t*/ms':>9} {'tflip/ms':>9} "
           f"{'I2nd/J':>10} {'ℓ/mm':>8} {'数差t/ms':>10} {'斜率t':>7}")
    print(hdr); print("-" * 74)
    for r in nom:
        def fmt(v, d=4):
            return f"{v:.{d}g}" if np.isfinite(v) else "—"
        print(f"{r['kc']:>9.3g} {fmt(r['delta0_mm'],3):>7} {fmt(r['theta0_deg'],3):>7} "
              f"{fmt(r['t_star_ms'],4):>9} {fmt(r['t_flip_ms'],4):>9} "
              f"{fmt(r['I_2nd'],3):>10} {fmt(r['ell_mm'],3):>8} "
              f"{fmt(r['num_err_t_ms'],2):>10} {fmt(r.get('slope_t',np.nan),3):>7}")
    print("-" * 74)

    # ---------- 判读提示（不下物理结论，只给判据） ----------
    osc = nom[0].get("osc_count", 0) if nom else 0
    s_t = np.array([r.get("slope_t", np.nan) for r in nom][1:], dtype=float)
    finite_s = s_t[np.isfinite(s_t)]
    rel = np.array([r.get("rel_chg_t", np.nan) for r in nom][1:], dtype=float)
    finite_r = rel[np.isfinite(rel)]
    print("判读提示：")
    if not nom or not np.any(np.isfinite(s_t)):
        print("  · 名义工况数据不足，无法判读。")
    else:
        if osc >= 2:
            print(f"  ⚠ t*(k_c) 出现 {osc} 次变向：提示刚性极限存在振荡/多解，")
            print("    应沿 Painlevé 方向排查（Leine 2002 / Champneys 2016），")
            print("    第二阶段分离表述考虑法向-切向 LCP 或碰撞律。")
        else:
            print("  · t*(k_c) 未见明显振荡。")
        if finite_r.size and np.nanmean(finite_r[-2:]) < 0.01:
            print("  · 后段相邻变化 < 1%：输出对 k_c 不再敏感，罚模型倾向自洽，")
            print("    k_c 取该区间的值即可，无需精确标定。")
        elif finite_r.size:
            print(f"  · 后段相邻变化仍达 {np.nanmean(finite_r[-2:]):.1%}：")
            print("    输出对 k_c 敏感，需记录标度指数或重审分离表述。")
        if finite_s.size:
            print(f"  · t* 对数斜率范围 [{finite_s.min():.2f}, {finite_s.max():.2f}]"
                  + ("，可与 −1/2 基准对照" if mode == "fixed_delta0" else ""))
    ne = np.array([r["num_err_t_ms"] for r in nom], dtype=float)
    ne = ne[np.isfinite(ne)]
    ts = np.array([r["t_star_ms"] for r in nom], dtype=float)
    ts = ts[np.isfinite(ts)]
    if ne.size and ts.size and np.max(ne) > 0.01 * np.max(ts):
        print("  ⚠ 双容差数值误差偏大：先收紧容差/加密步数，再谈物理判读。")
    print(f"\n输出目录：{os.path.abspath(outdir)}")
    print("  S1_kc扫描_结果表.csv / S1_kc扫描_收敛曲线.png / "
          "S1_kc扫描_曲线族.png / S1_kc扫描_诊断.png")


# ============================== 主流程 ==============================
def scan(mode, mu, outdir, want_plots=True):
    os.makedirs(outdir, exist_ok=True)
    rows, trajs = [], {}
    n0_nom = NOMINAL[0] if mode == "fixed_N0" else None
    for kc in KC_LIST:
        n0s = N0_LIST if mode == "fixed_N0" else [kc * DELTA0_REF]
        for N0 in n0s:
            for om in OM_LIST:
                r, traj = run_dual(kc, N0, om, mu)
                r["mode"] = mode
                r["is_nominal"] = bool(
                    (mode == "fixed_delta0" and abs(om - NOMINAL[1]) < 1e-12)
                    or (mode == "fixed_N0" and abs(N0 - n0_nom) < 1e-12
                        and abs(om - NOMINAL[1]) < 1e-12))
                if r["is_nominal"] and traj is not None:
                    trajs[kc] = traj
                rows.append(r)
        print(f"  [S1] k_c = {kc:.3g} N/m 完成 "
              f"（{len(n0s) * len(OM_LIST)} 组 ×2 容差）")
    add_diagnostics(rows, mode)
    save_csv(rows, os.path.join(outdir, "S1_kc扫描_结果表.csv"))
    if want_plots:
        nom = [r for r in rows if r["is_nominal"]]
        plot_convergence(nom, outdir)
        plot_family(trajs, outdir)
        plot_diagnostics(nom, outdir)
    print_summary(rows, mode, mu, outdir)
    return rows


def main():
    ap = argparse.ArgumentParser(description="S1：k_c 刚性极限收敛扫描")
    ap.add_argument("--mode", choices=["fixed_N0", "fixed_delta0"],
                    default="fixed_N0", help="扫描模式（默认 fixed_N0）")
    ap.add_argument("--mu", type=float, default=MU_DEFAULT, help="动摩擦系数")
    ap.add_argument("--outdir", default="./S1_kc扫描输出", help="输出目录")
    ap.add_argument("--kc-lo", type=float, default=3.0, help="log10(k_c) 下限")
    ap.add_argument("--kc-hi", type=float, default=7.0, help="log10(k_c) 上限")
    ap.add_argument("--kc-n", type=int, default=9, help="k_c 点数")
    ap.add_argument("--no-plots", action="store_true", help="只出 CSV 不出图")
    args = ap.parse_args()

    global KC_LIST
    KC_LIST = np.logspace(args.kc_lo, args.kc_hi, args.kc_n)
    scan(args.mode, args.mu, args.outdir, want_plots=not args.no_plots)


if __name__ == "__main__":
    main()
