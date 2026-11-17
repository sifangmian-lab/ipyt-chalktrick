"""Recompute the simulation studies using the v3 eccentric two-spring model."""
from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict, replace
from itertools import product
import json
from pathlib import Path

import numpy as np
import scipy

from deposit_profile import deposit_profile, profile_integral, total_deposit, trace_coordinate
from observation_model import camera_response, footprint_kernel, profile_to_pixel_mass, render_mass
from v3_single_mark import Benchmark, InfeasiblePreload, Parameters, integrate, require, write_table
from v3_speed_design import assess


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "simulations" / "V3_全量复算_2026-11-17"
BASE = Parameters()
RATIO = BASE.mu_s / BASE.mu_k
CONTACTS = (("soft", 2e4), ("medium", 2e5), ("stiff", 2e6), ("rigid", None))
ANGLES = (5, 10, 15, 20, 25, 30, 35, 40, 45, 55)
FRICTION = (.15, .19, .22, .25, .30, .40, .50, .62, .80, 1.00)
E1_LOWER, E1_UPPER = (h / np.sqrt(60.0**2 - h**2) for h in (33.9, 34.0))
E1 = {"lower": E1_LOWER, "center": (E1_LOWER + E1_UPPER)/2, "upper": E1_UPPER}


def cases():
    yield "S0", "nominal", {}
    for kn, load, omega in product(np.logspace(3, 7, 9), (1, 2, 4, 8), (0, 3, 6)):
        yield "S1", "fixed_preload", dict(kn_n_m=kn, preload_n=load, initial_omega=omega)
    for kn, omega in product(np.logspace(3, 7, 9), (0, 3, 6)):
        yield "S1", "fixed_penetration", dict(kn_n_m=kn, preload_n=kn * .0002, initial_omega=omega)
    for kn, eps in product((1e3, 1e4, 2e4, 31622.78),
                           (1e-4, 1e-3, 3e-3, .01, .02, .05, .1, .3, 1.0)):
        yield "S2", "tanh", dict(kn_n_m=kn, tanh_speed_m_s=eps, mu_k=.35)
    yield "S3", "tanh_reference", dict(tanh_speed_m_s=.02, mu_k=.35)
    for sigma0, sigma1 in product((1e5, 3e5), (0, 30, 100)):
        yield "S3", "lugre", dict(sigma0_n_m=sigma0, sigma1_ns_m=sigma1,
                                  sigma2_ns_m=5, mu_k=.35, mu_s=.35 * 1.15)
    for mu, hl, kn in product((.20, .35, .50, .65, .80, 1.00, 1.20),
                              (.70, .78, .85, .90, .9375, .96, .98),
                              (5e3, 1e4, 2e4, 5e4, 1e5)):
        # The old height ratio fixes the starting angle through cos(theta)=h/L.
        yield "S4v1", "initial_geometry", dict(angle_deg=float(np.degrees(np.arccos(hl))),
                                                mu_k=mu, mu_s=mu * RATIO, kn_n_m=kn)
    yield "S4", "A", {}
    yield "S4", "B", dict(kx_n_m=1000, ky_n_m=1e5)
    yield "S4", "C", dict(board_speed_m_s=.1)
    for key, values in (
        ("kn_n_m", (2e3, 2e4, 2e5, 2e6)), ("sigma0_n_m", (1e4, 1e5, 1e6)),
        ("cn_ns_m", (.1, 1, 10)), ("rtol", (1e-6, 1e-8, 1e-10)),
        ("max_step", (20e-6, 10e-6, 5e-6)),
    ):
        for val in values:
            yield "S5", key, {key: val}
    for kn, tol in product((2e3, 2e4, 2e5, 2e6), (1e-8, 1e-10)):
        yield "S5b", "penalty", dict(kn_n_m=kn, rtol=tol, max_step=1e-6 if tol == 1e-10 else 1e-5)
    for cn in (.1, 1, 10):
        yield "S5b", "damping", dict(kn_n_m=2e6, cn_ns_m=cn)
    for ky, rigid in product((1e5, 3e5, 1e6, 3e6), (False, True)):
        yield "S5b", "vertical_stiffness", dict(ky_n_m=ky, kn_n_m=2e6, rigid=rigid)
    for tol in (1e-8, 1e-10):
        yield "S5b", "rigid", dict(rigid=True, rtol=tol)
    for kx, ky in product((100, 300, 1000), (3e5, 1e6, 3e6)):
        yield "S6", "spring_plane", dict(kx_n_m=kx, ky_n_m=ky)
    for arm in (0.0, .125, .25, .375, .5):
        yield "S6", "support_arm", dict(support_arm_ratio=arm)
    for speed in (.05, .1, .2, .5, 1.0):
        yield "S6", "board_speed", dict(board_speed_m_s=speed)
    for (name, kn), angle, mu_s in product(CONTACTS, ANGLES, FRICTION):
        yield "S6b", name, dict(angle_deg=angle, mu_s=mu_s, mu_k=mu_s/RATIO,
                                 rigid=kn is None, kn_n_m=kn or BASE.kn_n_m)
    for (name, kn), angle, ratio in product(CONTACTS, (20, 35), (1, 1.2, 1.5, 2)):
        yield "S6b", f"{name}_ratio", dict(angle_deg=angle, mu_s=.62, mu_k=.62/ratio,
                                          rigid=kn is None, kn_n_m=kn or BASE.kn_n_m)
    for speed, mu_s in product((.1, .2, .5),
                                (.01, .02, .03, .04, .05, .06, .08, .10, .15, .20, .22, .25,
                                 .35, .37, .40, .50, .62)):
        yield "S6b", "speed_friction", dict(board_speed_m_s=speed, mu_s=mu_s, mu_k=mu_s/RATIO,
                                           horizon={.1: .03, .2: .02, .5: .012}[speed])
    for angle, mu_s, cn in product((20, 25, 30, 35, 40), (.4, .62, .8), (1, 10, 100)):
        yield "S6b", "damping", dict(angle_deg=angle, mu_s=mu_s, mu_k=mu_s/RATIO,
                                      kn_n_m=2e6, cn_ns_m=cn)
    yield "S7", "dry", {}
    yield "S7", "viscous", dict(sigma1_ns_m=30, sigma2_ns_m=5)
    yield "S7", "tanh", dict(tanh_speed_m_s=.02)
    e1_cases = {"reference": {}, **{tag: dict(mu_s=mu) for tag, mu in E1.items()},
                "mu_k_040": dict(mu_s=E1["center"], mu_k=.40),
                "mu_k_060": dict(mu_s=E1["center"], mu_k=.60),
                "vs_005": dict(mu_s=E1["center"], stribeck_speed_m_s=.005),
                "vs_020": dict(mu_s=E1["center"], stribeck_speed_m_s=.020)}
    for speed, (name, changes) in product((.005, .010, .020, .050, .100, .200, .500), e1_cases.items()):
        yield "E1", name, dict(changes, board_speed_m_s=speed)
    for speed in (.1, .5, 1.0):
        yield "S9", "continuous_contact", dict(board_speed_m_s=speed, horizon=.05,
                                               stop_at_velocity_zero=False)


def calculate(study, case, changes, *, keep_series=False):
    changes = changes.copy()
    controls = dict(rtol=1e-8, max_step=20e-6, horizon=.012, initial_omega=0.0,
                    stop_at_velocity_zero=True)
    for key in controls:
        if key in changes:
            controls[key] = changes.pop(key)
    p = replace(BASE, **changes)
    row = {"study": study, "case": case, **asdict(p), **controls}
    row["static_friction_required_n"] = abs(
        (p.preload_n - p.density_kg_m3 * np.pi * p.radius_m**2 * p.length_m * p.gravity_m_s2
         * p.support_arm_ratio / (.5 + p.support_arm_ratio)) * np.tan(np.radians(p.angle_deg)))
    try:
        model = Benchmark(p)
    except InfeasiblePreload:
        return row | {"end_reason": "infeasible_preload"}, None

    if p.tanh_speed_m_s is None:
        f, torque, _ = model.loads(model.initial, 0.0)
        require(np.linalg.norm(f) < 1e-8 and abs(torque) < 1e-9, "Static balance failed.")
    force, torque, z_dot = model.loads(model.initial, p.board_speed_m_s)
    expected = np.r_[model.initial[3:6], force/model.mass, torque/model.inertia, z_dot]
    require(np.allclose(model.rhs(0, model.initial), expected, atol=1e-7, rtol=1e-10),
            "Scalar and sampled loads differ.")

    atol = np.full(7, controls["rtol"] * .01)
    # The bristle state can be sub-nanometre as a stiff contact unloads.
    atol[6] = min(atol[6], 1e-11)
    sol, end, events = integrate(model, method="LSODA", atol=atol, **controls)
    if not controls["stop_at_velocity_zero"]:
        end = next((name for name in ("normal_force_zero", "geometric_gap_opening")
                    if events.get(name) is not None), "time_limit")
    t = np.linspace(0.0, sol.t[-1], max(4001, int(sol.t[-1] / 2e-6) + 1))
    t = np.unique(np.r_[t, sol.t, [v for v in events.values() if v is not None]])
    state = sol.sol(t)
    tip, _, _, _ = model.geometry(state)
    normal, resistance, z_dot, _, speed = model.sample_contact(state, p.board_speed_m_s)
    xi = tip[0] - tip[0, 0] + p.board_speed_m_s * t
    m1_rate = normal * np.abs(speed)
    work_rate = resistance * speed
    m2_rate = work_rate if p.tanh_speed_m_s is not None else (
        work_rate - p.sigma0_n_m * state[6] * z_dot)
    m1 = total_deposit(t, m1_rate)
    energy = (np.zeros_like(t) if p.tanh_speed_m_s is not None else
              .5 * p.sigma0_n_m * state[6]**2)
    storage = energy[-1] - energy[0]
    work_intervals = .5 * (work_rate[:-1] + work_rate[1:]) * np.diff(t)
    m2_intervals = work_intervals - np.diff(energy)
    work = float(work_intervals.sum())
    m2 = float(m2_intervals.sum())
    centres, m1_density = deposit_profile(t, m1_rate, xi, sigma_bins=0)
    # Storage differences resolve the narrow release transient more reliably
    # than sampling z_dot. With bristle damping the residual can be negative.
    edges = np.linspace(xi.min(), xi.max(), len(centres) + 1)
    if np.all(np.diff(xi) > 0):
        m1_intervals = .5 * (m1_rate[:-1] + m1_rate[1:]) * np.diff(t)
        m1_density, m2_density = [
            np.diff(np.interp(edges, xi, np.r_[0, weights.cumsum()])) / np.diff(edges)
            for weights in (m1_intervals, m2_intervals)
        ]
    else:
        midpoints = .5 * (xi[:-1] + xi[1:])
        m2_density = np.histogram(midpoints, edges, weights=m2_intervals)[0] / np.diff(edges)
    coordinate_error = float(np.max(np.abs(xi - trace_coordinate(t, speed))))
    mapping_error = max(abs(profile_integral(centres, m1_density) - m1) / max(m1, 1e-12),
                        abs(profile_integral(centres, m2_density) - m2) /
                        max(np.abs(m2_intervals).sum(), 1e-12))
    energy_error = abs(work - storage - m2) / max(np.abs(work_intervals).sum(), 1e-12)
    rate_error = abs(total_deposit(t, m2_rate) - m2) / max(np.abs(work_intervals).sum(), 1e-12)
    require(mapping_error < 1e-9, "Spatial integral failed.")
    require(coordinate_error < 2e-6, "Geometric position disagrees with integrated speed.")
    if p.rigid:
        require(np.max(abs(tip[1])) < 2e-6, "Rigid-contact gap grew off the constraint.")
    first = events["first_velocity_zero"]
    first_xi = None
    if first is not None:
        first_xi = float(xi[np.searchsorted(t, first)])
    row.update(end_reason=end, end_time_s=t[-1], first_stop_s=first,
               force_zero_s=events.get("normal_force_zero"),
               gap_opening_s=events.get("geometric_gap_opening"),
               first_stop_xi_m=first_xi, span_m=np.ptp(xi), travel_m=total_deposit(t, abs(speed)),
               normal_min_n=normal.min(), normal_max_n=normal.max(), normal_end_n=normal[-1],
               normal_initial_n=normal[0], relative_speed_min_m_s=speed.min(),
               relative_speed_final_m_s=speed[-1], tip_height_max_m=tip[1].max(),
               velocity_downcrossings=len(sol.t_events[0]), m1_integral_n_m=m1, m2_residual_j=m2,
               m2_negative_part_j=total_deposit(t, np.maximum(-m2_rate, 0)),
               friction_work_j=work, bristle_storage_change_j=storage,
               energy_balance_relative_error=energy_error, rate_quadrature_relative_error=rate_error,
               spatial_balance_relative_error=mapping_error,
               coordinate_error_m=coordinate_error)
    if end == "first_velocity_zero":
        for key, density in (("m1", m1_density), ("m2", m2_density)):
            middle = density[40:160].mean()
            row[f"{key}_last_5pct_to_middle"] = density[-10:].mean()/middle if middle > 0 else None
    series = (t, xi, normal, speed, tip[1], centres, m1_density, m2_density) if keep_series else None
    return row, series


def save_rows(path, rows):
    fields = list(dict.fromkeys(key for row in rows for key in row))
    def format_cell(value):
        return f"{value:.10g}" if isinstance(value, (float, np.floating)) else value
    write_table(path, {key: [format_cell(row.get(key)) for row in rows] for key in fields})


def e2_span_um(speed_mm_s, mu_s, fitted):
    mu_k, log_vs = fitted
    row, _ = calculate("E2", "uncertainty", dict(board_speed_m_s=speed_mm_s/1000,
                       mu_s=mu_s, mu_k=mu_k, stribeck_speed_m_s=np.exp(log_vs)/1000,
                       rtol=1e-10))
    require(row["end_reason"] == "first_velocity_zero", "E2 span requires a first velocity zero.")
    return 1e6 * row["first_stop_xi_m"]


def observations(outdir, centres, m1_density, m2_density):
    import matplotlib.pyplot as plt

    from demo_observation_model import register_projection

    prior = ROOT / "simulations" / "S8_有限足迹观测模型"
    observations = {mark: np.loadtxt(prior / f"mark{mark}_comparison.csv",
                                    encoding="utf-8-sig", delimiter=",", skiprows=1, usecols=1)
                    for mark in (44, 51)}
    footprints = {"point": footprint_kernel(0, 0), "upright": footprint_kernel(tilt_dx_dy=0),
                  "tilted": footprint_kernel()}
    rows, scores = [], []
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.6), layout="constrained")
    for (source, density), ax in zip((("M1", m1_density), ("M2", m2_density)), axes):
        for span, (name, kernel) in product((2, 4, 8, 16), footprints.items()):
            grid, mass = profile_to_pixel_mass(centres, density, span)
            total = mass.sum()
            x, _, field = render_mass(mass / total, kernel, grid[0], .7)
            error = float(abs(field.sum() - 1))
            require(error < 1e-10, "Footprint or optical blur lost deposit.")
            for mark, obs in observations.items():
                fit = register_projection(x, field, obs)
                scores.append(dict(deposition=source, span_px=span, footprint=name, mark=mark,
                                   shift_px=fit["x_shift_px"], gain=fit["gain"], nrmse=fit["nrmse"]))
            for mode in ("linear", "saturating", "clipped"):
                gain = (3 if mode == "clipped" else 1) / field.max()
                projection = camera_response(field, gain=gain, mode=mode).sum(axis=0)
                rows.extend(dict(deposition=source, span_px=span, footprint=name, camera=mode,
                                 source_integral=total, mass_balance_relative_error=error,
                                 x_px=xx, projected_intensity=yy) for xx, yy in zip(x, projection))
                if span == 4 and mode == "linear":
                    ax.plot(x, projection/projection.max(), label=name)
        ax.set(title=f"{source}, 4 px source span", xlabel="Position (px)", ylabel="Intensity / peak")
        ax.legend(frameon=False)
    fig.savefig(outdir / "footprints.png", dpi=155)
    plt.close(fig)
    save_rows(outdir / "footprints.csv", rows)
    save_rows(outdir / "photo_comparison.csv", scores)
    return {"scenarios": 72, "photo_comparisons": len(scores), "max_mass_balance_relative_error":
            max(row["mass_balance_relative_error"] for row in rows)}


def plot_results(outdir, rows, post_stop):
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch

    categories = {"infeasible_preload": "#d5d9db", "first_velocity_zero": "#25786f",
                  "time_limit": "#e1ab59", "normal_force_zero": "#b14848",
                  "geometric_gap_opening": "#74428a"}
    fig, axes = plt.subplots(2, 2, figsize=(9.4, 6.2), layout="constrained")
    for (name, _), ax in zip(CONTACTS, axes.flat):
        outcome = {(row["angle_deg"], row["mu_s"]): row["end_reason"]
                   for row in rows if row["study"] == "S6b" and row["case"] == name}
        codes = list(categories)
        z = [[codes.index(outcome[(angle, mu)]) for mu in FRICTION] for angle in ANGLES]
        ax.imshow(z, aspect="auto", origin="lower", interpolation="nearest",
                  cmap=ListedColormap(list(categories.values())), vmin=-.5, vmax=len(codes)-.5)
        ax.set(title=name.capitalize(), xlabel="Static friction coefficient", ylabel="Initial angle (deg)",
               xticks=range(len(FRICTION)), xticklabels=[f"{x:g}" for x in FRICTION],
               yticks=range(0, len(ANGLES), 2), yticklabels=ANGLES[::2])
        ax.tick_params(axis="x", labelsize=7)
    labels = ["Infeasible preload", "First velocity zero", "Time limit", "Normal force zero", "Gap opening"]
    fig.legend([Patch(color=c) for c in categories.values()], labels, loc="outside lower center",
               ncol=3, frameon=False)
    fig.savefig(outdir / "angle_friction.png", dpi=155)
    plt.close(fig)
    fig, axes = plt.subplots(2, 1, figsize=(8.0, 5.0), layout="constrained", sharex=True)
    for speed, trace in post_stop:
        t, _, normal, velocity, _, *_ = trace
        axes[0].plot(t*1e3, velocity, label=f"{speed:g} m/s")
        axes[1].plot(t*1e3, normal)
    axes[0].axhline(0, color=".55", lw=.8)
    axes[0].set(ylabel="Relative velocity (m/s)")
    axes[0].legend(frameon=False)
    axes[1].set(xlabel="Time (ms)", ylabel="Normal force (N)")
    fig.savefig(outdir / "post_stop.png", dpi=155)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=OUTPUT)
    args = parser.parse_args()
    rows, post_stop, traces, profiles = [], [], [], []
    studies = list(cases())
    last = None
    for study, case, changes in studies:
        if study != last:
            print(f"Computing {study}", flush=True)
            last = study
        try:
            row, trace = calculate(study, case, changes,
                                   keep_series=study in ("S0", "S4", "S7", "E1", "S9"))
        except RuntimeError as exc:
            raise RuntimeError(f"{study}/{case} {changes}: {exc}") from exc
        rows.append(row)
        if study == "S0":
            baseline = trace
        if study in ("S4", "S7", "E1"):
            speed = row["board_speed_m_s"]
            profiles.extend(dict(study=study, case=case, board_speed_m_s=speed,
                                 xi_m=x, m1_density_n=m1, m2_density_j_m=m2)
                            for x, m1, m2 in zip(*trace[-3:]))
        if study == "S9":
            speed = row["board_speed_m_s"]
            post_stop.append((speed, trace))
            t, xi, normal, velocity, height, *_ = trace
            samples = np.unique(np.r_[np.linspace(0, len(t)-1, 1201, dtype=int),
                                      np.searchsorted(t, row["first_stop_s"])
                                      if row["first_stop_s"] is not None else 0])
            traces.extend(dict(board_speed_m_s=speed, time_s=a, xi_m=b, normal_n=c,
                               relative_speed_m_s=d, tip_height_m=e)
                          for a, b, c, d, e in zip(t[samples], xi[samples], normal[samples],
                                                  velocity[samples], height[samples]))

    # Recompute the nominal case independently at tighter settings.
    reference = rows[0]
    tight, _ = calculate("S0", "tight", dict(rtol=1e-10, max_step=5e-6))
    for field in ("end_time_s", "span_m", "normal_end_n", "m1_integral_n_m", "m2_residual_j"):
        require(abs(reference[field]/tight[field] - 1) < 1e-4, f"Nominal convergence failed: {field}")
    rows.append(tight)

    import matplotlib
    matplotlib.use("Agg")
    args.outdir.mkdir(parents=True, exist_ok=True)
    save_rows(args.outdir / "cases.csv", rows)
    save_rows(args.outdir / "post_stop.csv", traces)
    save_rows(args.outdir / "profiles.csv", profiles)
    checks = observations(args.outdir, *baseline[-3:])
    print("Computing E2", flush=True)
    assess(args.outdir, E1, BASE.mu_k, e2_span_um, save_rows)
    plot_results(args.outdir, rows, post_stop)
    valid = [row for row in rows if row["end_reason"] != "infeasible_preload"]
    summary = {
        "cases": len(rows),
        "by_study": {study: dict(Counter(row["end_reason"] for row in rows if row["study"] == study))
                     for study in dict.fromkeys(row["study"] for row in rows)},
        "by_outcome": dict(Counter(row["end_reason"] for row in rows)),
        "max_coordinate_error_m": max(row["coordinate_error_m"] for row in valid),
        "max_spatial_balance_relative_error": max(row["spatial_balance_relative_error"] for row in valid),
        "max_energy_balance_relative_error": max(row["energy_balance_relative_error"] for row in valid),
        "max_rate_quadrature_relative_error": max(row["rate_quadrature_relative_error"] for row in valid),
        "max_pre_stop_rate_quadrature_relative_error": max(
            row["rate_quadrature_relative_error"] for row in valid if row["end_reason"] == "first_velocity_zero"),
        "footprints": checks,
        "environment": {"numpy": np.__version__, "scipy": scipy.__version__,
                        "matplotlib": matplotlib.__version__},
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n",
                                              encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
