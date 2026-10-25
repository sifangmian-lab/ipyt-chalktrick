# -*- coding: utf-8 -*-
"""空间沉积强度 I(ξ) 的统一后处理。

背景
----
主模型、S2、S3 原先各自实现 I(ξ) 直方图，权重互不一致：
  · 旧主模型：weights = q*|dξ|   （量纲/守恒错误，多乘 |v_rel|）
  · 旧 S2    ：weights = μN*dt   （量纲/守恒错误，少乘 |v_rel|，转折点发散）
  · 旧 S3    ：weights = |F_f*v_rel|*dt（正确）
本模块统一为：

    q(t) = α |F_f(t) v_rel(t)|                # 单位时间沉积率（功率型，J/s）
    ξ(t) = ∫ v_rel dt                          # 板面系痕迹坐标
    I(ξ) = Σ_{t: ξ(t)=ξ} q(t) / |v_rel(t)| ≈ Σ q dt / Δξ

即直方图权重必须是 q*dt，而不是 q*|dξ| 或 μN*dt。
守恒检验：∫ I(ξ) dξ = ∫ q dt = I_2nd（α=1 时）。

数值细节：
  · 用中点梯形权重，避免端点一阶误差；
  · 回划折叠自然由同一 ξ 分箱累加实现；
  · 轻度高斯平滑默认 σ=2 个 bin，与旧脚本一致（用 numpy 卷积，避免额外依赖）。
"""
from __future__ import annotations

import numpy as np


def _gaussian_smooth(y: np.ndarray, sigma_bins: float) -> np.ndarray:
    if sigma_bins is None or sigma_bins <= 0 or y.size == 0:
        return y
    rad = max(1, int(4.0 * sigma_bins))
    x = np.arange(-rad, rad + 1, dtype=float)
    ker = np.exp(-0.5 * (x / float(sigma_bins)) ** 2)
    ker /= ker.sum()
    ypad = np.pad(y, rad, mode="edge")
    return np.convolve(ypad, ker, mode="same")[rad:-rad]


def trace_coordinate(t: np.ndarray, v_rel: np.ndarray) -> np.ndarray:
    """ξ(t) = ∫ v_rel dt（trapezoid，与旧脚本 cumsum 中点公式一致）。"""
    t = np.asarray(t, dtype=float)
    v_rel = np.asarray(v_rel, dtype=float)
    return np.concatenate([[0.0],
                           np.cumsum(0.5 * (v_rel[:-1] + v_rel[1:]) * np.diff(t))])


def deposit_profile(t, q, xi=None, nbins: int = 200, sigma_bins: float = 2.0):
    """由沉积率 q(t) 与痕迹坐标 ξ(t) 计算空间强度 I(ξ)。

    Parameters
    ----------
    t : (N,) 时间序列（s，单调）
    q : (N,) 单位时间沉积率（J/s；例如 |F_f*v_rel| 或 μN|v_rel|）
    xi : (N,) 痕迹坐标（m）；None 时由 v_rel 无法恢复，需显式传入
    nbins : 分箱数
    sigma_bins : 高斯平滑宽度（bin）；<=0 关闭

    Returns
    -------
    xc : (nbins,) 箱中心（m）
    I  : (nbins,) I(ξ)（J/m）
    """
    t = np.asarray(t, dtype=float)
    q = np.asarray(q, dtype=float)
    if xi is None:
        raise ValueError("必须显式传入 xi=trace_coordinate(t,v_rel)；"
                         "本函数不再从 q 猜 v_rel。")
    xi = np.asarray(xi, dtype=float)
    if t.shape != q.shape or t.shape != xi.shape:
        raise ValueError(f"t/q/xi 形状不一致: {t.shape}, {q.shape}, {xi.shape}")
    if t.size < 2:
        raise ValueError("至少需要 2 个采样点")

    dt = np.diff(t)
    qm = 0.5 * (q[:-1] + q[1:])
    xim = 0.5 * (xi[:-1] + xi[1:])

    x_lo, x_hi = float(np.min(xi)), float(np.max(xi))
    if not np.isfinite([x_lo, x_hi]).all() or x_hi <= x_lo:
        return np.array([]), np.array([])

    edges = np.linspace(x_lo, x_hi, int(nbins) + 1)
    widths = np.diff(edges)
    hist, _ = np.histogram(xim, bins=edges, weights=qm * dt)
    I = hist / widths
    I = _gaussian_smooth(I, sigma_bins)
    xc = 0.5 * (edges[:-1] + edges[1:])
    return xc, I


def profile_integral(xc: np.ndarray, I: np.ndarray) -> float:
    """∫ I(ξ) dξ（用于守恒自检）。

    注意：本函数的分箱均匀，故直接 Σ I·Δξ；不要用梯形法对箱中心积分
    （会漏掉两端各半个箱宽，产生约 1/nbins 的系统偏差）。
    """
    if xc.size < 2 or I.size < 2:
        return 0.0
    dx = float(np.mean(np.diff(xc)))
    return float(np.sum(I) * dx)


def total_deposit(t, q) -> float:
    """∫ q dt（应对应 ∫I dξ，α=1）。"""
    return float(np.trapezoid(q, np.asarray(t, dtype=float)))


if __name__ == "__main__":
    # 自检 1：单调常速、常沉积率 → I=q/v_rel，且守恒
    t = np.linspace(0.0, 1e-3, 5001)
    v0 = 0.5
    xi = trace_coordinate(t, np.full_like(t, v0))
    q = np.full_like(t, 0.7)          # J/s
    xc, I = deposit_profile(t, q, xi=xi, nbins=200, sigma_bins=0.0)
    print("自检1 I 中位数 =", np.median(I), "期望 =", q[0] / v0)
    print("自检1 守恒相对误差 =",
          abs(profile_integral(xc, I) - total_deposit(t, q)) / total_deposit(t, q))

    # 自检 2：一个来回（回划折叠）→ 折叠区强度应约加倍
    tt = np.linspace(0.0, 2e-3, 8001)
    v = 0.5 * np.sign(np.cos(np.pi * (tt / 2e-3)))  # +0.5 前半、-0.5 后半
    xi2 = trace_coordinate(tt, v)
    q2 = np.full_like(tt, 0.7)
    xc2, I2 = deposit_profile(tt, q2, xi=xi2, nbins=200, sigma_bins=1.0)
    print("自检2 ∫I dξ/I_2nd =",
          profile_integral(xc2, I2) / total_deposit(tt, q2))
