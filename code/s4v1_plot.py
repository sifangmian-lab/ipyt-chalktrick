# -*- coding: utf-8 -*-
"""极简 PNG 绘图（Pillow backend）：热图与折线图。

本环境无 matplotlib（DSH 运行时只带 numpy/pandas/Pillow），因此用 PIL 直接绘制。
只实现 S4 需要的两类图：二维参数平面热图、一维折线/对数折线。
"""
from __future__ import annotations

import math
import os

import numpy as np
from PIL import Image, ImageDraw, ImageFont

_FONT = r"C:\Windows\Fonts\msyh.ttc"
_FONT_BOLD = r"C:\Windows\Fonts\msyhbd.ttc"

_CMAPS = {
    "viridis": [(68, 1, 84), (72, 40, 120), (62, 74, 137), (49, 104, 142),
                (38, 130, 142), (31, 158, 137), (53, 183, 121), (109, 205, 89),
                (180, 222, 44), (253, 231, 37)],
    "coolwarm": [(59, 76, 192), (144, 178, 254), (220, 220, 220), (245, 156, 125), (180, 4, 38)],
    "categorical": [(46, 125, 50), (249, 168, 37), (239, 108, 0), (198, 40, 40)],
}


def _font(size, bold=False):
    path = _FONT_BOLD if bold and os.path.exists(_FONT_BOLD) else _FONT
    try:
        return ImageFont.truetype(path, size)
    except Exception:
        return ImageFont.load_default()


def _log_ticks(lo, hi):
    out = []
    if lo <= 0 or hi <= 0:
        return out
    for e in range(int(math.floor(math.log10(lo))), int(math.ceil(math.log10(hi))) + 1):
        for m in (1, 2, 5):
            v = m * 10.0 ** e
            if lo - 1e-12 <= v <= hi + 1e-12:
                out.append(v)
    return out


def _vtext(img, text, font, fill, center):
    """竖排（旋转 90°）文字。"""
    box = font.getbbox(text)
    w, h = box[2] - box[0], box[3] - box[1]
    tmp = Image.new("RGBA", (max(w, 1) + 10, max(h, 1) + 10), (0, 0, 0, 0))
    d = ImageDraw.Draw(tmp)
    d.text((5 - box[0], 5 - box[1]), text, font=font, fill=fill)
    tmp = tmp.rotate(90, expand=True)
    img.paste(tmp, (int(center[0] - tmp.width / 2), int(center[1] - tmp.height / 2)), tmp)


def _halo_text(d, xy, text, font, fill, anchor="rs", halo=(255, 255, 255)):
    for dx, dy in ((-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (1, 1), (-1, 1), (1, -1)):
        d.text((xy[0] + dx, xy[1] + dy), text, font=font, fill=halo, anchor=anchor)
    d.text(xy, text, font=font, fill=fill, anchor=anchor)


def cmap(name, v):
    """连续色标取色，v ∈ [0,1]。"""
    stops = _CMAPS[name]
    v = min(max(float(v), 0.0), 1.0)
    pos = v * (len(stops) - 1)
    i = int(math.floor(pos))
    if i >= len(stops) - 1:
        return stops[-1]
    f = pos - i
    a, b = stops[i], stops[i + 1]
    return tuple(int(round(a[k] + (b[k] - a[k]) * f)) for k in range(3))


def _nice_ticks(lo, hi, n=5):
    if not np.isfinite([lo, hi]).all() or hi <= lo:
        return [lo]
    span = hi - lo
    raw = span / max(n - 1, 1)
    mag = 10 ** math.floor(math.log10(raw))
    for m in (1, 2, 2.5, 5, 10):
        step = m * mag
        if raw <= step:
            break
    first = math.ceil(lo / step) * step
    ticks = []
    t = first
    while t <= hi + 1e-12:
        ticks.append(round(t, 12))
        t += step
    return ticks


def _tick_label(v):
    if v == 0:
        return "0"
    if abs(v) >= 1e4 or abs(v) < 1e-3:
        return f"{v:.3g}"
    s = f"{v:.4g}"
    return s


def heatmap(path, Z, row_labels, col_labels, title="", cbar_label="",
            cmap_name="viridis", vmin=None, vmax=None, overlay=None,
            categorical=False, legend=None, invert_rows=True, figsize=(760, 560),
            annotation=None):
    """Z: (nrow, ncol) 矩阵；Z[i,j] 对应 row_labels[i], col_labels[j]。"""
    W, H = figsize
    ml, mr = 92, 132
    mt = 104 if legend else 64
    mb = 78
    img = Image.new("RGB", (W, H), "white")
    d = ImageDraw.Draw(img)
    f_title = _font(19, bold=True); f_lab = _font(15); f_tick = _font(13)
    f_small = _font(12)

    nr, nc = Z.shape
    pw, ph = W - ml - mr, H - mt - mb
    cw, ch = pw / nc, ph / nr

    Zf = np.array(Z, dtype=float)
    finite = Zf[np.isfinite(Zf)]
    if vmin is None:
        vmin = float(finite.min()) if finite.size else 0.0
    if vmax is None:
        vmax = float(finite.max()) if finite.size else 1.0
    if vmax == vmin:
        vmax = vmin + 1.0

    for i in range(nr):
        for j in range(nc):
            x0 = ml + j * cw
            y0 = mt + (nr - 1 - i if invert_rows else i) * ch
            v = Zf[i, j]
            if categorical:
                if np.isfinite(v):
                    c = _CMAPS[cmap_name][int(v) % len(_CMAPS[cmap_name])]
                else:
                    c = (235, 235, 235)
            else:
                c = (235, 235, 235) if not np.isfinite(v) else cmap(cmap_name, (v - vmin) / (vmax - vmin))
            d.rectangle([x0, y0, x0 + cw, y0 + ch], fill=c, outline=(255, 255, 255), width=1)
            if annotation is not None:
                txt = annotation[i][j]
                if txt != "":
                    lum = 0.299 * c[0] + 0.587 * c[1] + 0.114 * c[2]
                    d.text((x0 + cw / 2, y0 + ch / 2), str(txt), font=f_small,
                           fill=(255, 255, 255) if lum < 140 else (25, 25, 25), anchor="mm")

    # 坐标轴
    d.rectangle([ml, mt, ml + pw, mt + ph], outline=(40, 40, 40), width=1)
    for j, lab in enumerate(col_labels):
        cx = ml + (j + 0.5) * cw
        d.text((cx, mt + ph + 8), str(lab), font=f_tick, fill=(30, 30, 30), anchor="ma")
    for i, lab in enumerate(row_labels):
        cy = mt + (nr - 1 - i if invert_rows else i) * ch + ch / 2
        d.text((ml - 10, cy), str(lab), font=f_tick, fill=(30, 30, 30), anchor="rm")
    d.text((ml + pw / 2, mt + ph + 40), "h / L", font=f_lab, fill=(20, 20, 20), anchor="ma")
    d.text((26, mt + ph / 2), "μ", font=f_lab, fill=(20, 20, 20), anchor="mm")
    d.text((ml + pw / 2, 24), title, font=f_title, fill=(10, 10, 10), anchor="ma")

    # 叠加曲线（数据坐标：列索引, 行索引）
    if overlay:
        xs, ys, col, lab = overlay
        pts = []
        for xv, yv in zip(xs, ys):
            if not (np.isfinite(xv) and np.isfinite(yv)):
                continue
            pts.append((ml + xv * cw, mt + (nr - 1 - yv) * ch))
        if len(pts) >= 2:
            d.line(pts, fill=col, width=3)
        if lab:
            lx, ly = pts[len(pts) // 2]
            tw = int(d.textlength(lab, font=f_small))
            d.rectangle([lx - tw / 2 - 4, ly - 18, lx + tw / 2 + 4, ly - 2],
                        fill=(255, 255, 255), outline=col)
            d.text((lx, ly - 10), lab, font=f_small, fill=col, anchor="mm")

    if legend:
        x = ml
        y = mt - 42
        for col, text in legend:
            d.rectangle([x, y, x + 22, y + 14], fill=col)
            d.text((x + 28, y + 7), text, font=f_small, fill=(20, 20, 20), anchor="lm")
            x += 34 + int(d.textlength(text, font=f_small))

    # 色条
    if not categorical:
        bx0, bx1 = W - mr + 26, W - mr + 54
        d.rectangle([bx0, mt, bx1, mt + ph], fill=(250, 250, 250), outline=(120, 120, 120))
        for k in range(ph):
            v = 1.0 - k / max(ph - 1, 1)
            d.line([(bx0 + 1, mt + k), (bx1 - 1, mt + k)], fill=cmap(cmap_name, v))
        d.rectangle([bx0, mt, bx1, mt + ph], outline=(80, 80, 80))
        for tv in _nice_ticks(vmin, vmax, 5):
            py = mt + ph - (tv - vmin) / (vmax - vmin) * ph
            d.line([(bx1, py), (bx1 + 5, py)], fill=(60, 60, 60), width=1)
            d.text((bx1 + 8, py), _tick_label(tv), font=f_tick, fill=(30, 30, 30), anchor="lm")
        d.text(((bx0 + bx1) / 2, mt + ph + 16), cbar_label, font=f_small, fill=(30, 30, 30), anchor="ma")
    img.save(path)


def lineplot(path, series, title="", xlabel="", ylabel="", xlog=False, ylog=False,
             hlines=(), figsize=(820, 560)):
    """series: [{'x':..., 'y':..., 'label':..., 'color':(r,g,b), 'marker':bool}]"""
    W, H = figsize
    ml, mr, mt, mb = 104, 178, 66, 76
    img = Image.new("RGB", (W, H), "white")
    d = ImageDraw.Draw(img)
    f_title = _font(19, bold=True); f_lab = _font(15); f_tick = _font(13); f_small = _font(12)

    allx = np.concatenate([np.asarray(s["x"], float) for s in series])
    ally = np.concatenate([np.asarray(s["y"], float) for s in series])
    ok = np.isfinite(allx) & np.isfinite(ally)
    allx, ally = allx[ok], ally[ok]
    if xlog:
        allx = allx[allx > 0]; xf = np.log10(allx)
    else:
        xf = allx
    if ylog:
        pos = ally > 0
        ally = ally[pos]
        yf = np.log10(ally)
    else:
        yf = ally
    x0, x1 = float(np.min(xf)), float(np.max(xf))
    y0, y1 = float(np.min(yf)), float(np.max(yf))
    px = 0.06 * (x1 - x0 or 1); py = 0.1 * (y1 - y0 or 1)
    x0, x1 = x0 - px, x1 + px
    y0, y1 = y0 - py, y1 + py

    pw, ph = W - ml - mr, H - mt - mb

    def TX(v):
        v = math.log10(v) if xlog else v
        return ml + (v - x0) / (x1 - x0) * pw

    def TY(v):
        v = math.log10(v) if ylog else v
        return mt + ph - (v - y0) / (y1 - y0) * ph

    d.rectangle([ml, mt, ml + pw, mt + ph], outline=(40, 40, 40), width=1)
    for g in range(1, 5):
        d.line([(ml, mt + ph * g / 5), (ml + pw, mt + ph * g / 5)], fill=(232, 232, 232))
        d.line([(ml + pw * g / 5, mt), (ml + pw * g / 5, mt + ph)], fill=(232, 232, 232))

    if xlog:
        xt = _log_ticks(10 ** x0, 10 ** x1)
    else:
        xt = _nice_ticks(x0, x1, 6)
    for tv in xt:
        if xlog:
            if not (10 ** x0 - 1e-9 <= tv <= 10 ** x1 + 1e-9):
                continue
        elif not (x0 - 1e-9 <= tv <= x1 + 1e-9):
            continue
        cx = TX(tv)
        d.line([(cx, mt + ph), (cx, mt + ph + 5)], fill=(60, 60, 60), width=1)
        d.text((cx, mt + ph + 10), _tick_label(tv), font=f_tick, fill=(30, 30, 30), anchor="ma")
    if ylog:
        yt = _log_ticks(10 ** y0, 10 ** y1)
    else:
        yt = _nice_ticks(y0, y1, 6)
    for tv in yt:
        cy = TY(tv)
        if cy < mt - 1 or cy > mt + ph + 1:
            continue
        d.line([(ml - 5, cy), (ml, cy)], fill=(60, 60, 60), width=1)
        d.text((ml - 9, cy), _tick_label(tv), font=f_tick, fill=(30, 30, 30), anchor="rm")

    for yv, col, lab in hlines:
        cy = TY(yv)
        if mt <= cy <= mt + ph:
            for k in range(ml, ml + pw, 12):
                d.line([(k, cy), (min(k + 6, ml + pw), cy)], fill=col, width=2)
            if lab:
                d.text((ml + pw - 6, cy - 6), lab, font=f_small, fill=col, anchor="rs")

    for s in series:
        xs = np.asarray(s["x"], float); ys = np.asarray(s["y"], float)
        col = s.get("color", (30, 90, 180))
        pts = [(TX(xv), TY(yv)) for xv, yv in zip(xs, ys)
               if np.isfinite(xv) and np.isfinite(yv)
               and (not xlog or xv > 0) and (not ylog or yv > 0)]
        if len(pts) >= 2:
            d.line(pts, fill=col, width=2)
        if s.get("marker", True):
            for (a, b) in pts:
                r = 3
                d.ellipse([a - r, b - r, a + r, b + r], fill=col)

    d.text((ml + pw / 2, mt + ph + 42), xlabel, font=f_lab, fill=(20, 20, 20), anchor="ma")
    d.text((ml + pw / 2, 24), title, font=f_title, fill=(10, 10, 10), anchor="ma")
    _vtext(img, ylabel, f_lab, (20, 20, 20), (30, mt + ph / 2))

    y = mt + 6
    for s in series:
        col = s.get("color", (30, 90, 180))
        d.line([(ml + pw + 16, y + 7), (ml + pw + 40, y + 7)], fill=col, width=3)
        d.text((ml + pw + 46, y + 7), s.get("label", ""), font=f_small, fill=(20, 20, 20), anchor="lm")
        y += 20

    img.save(path)
