"""Run the v3 first-slip benchmark for a common, eccentric spring attachment."""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass
import json
from math import cos, exp, sin, tanh
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp

from deposit_profile import total_deposit, trace_coordinate


@dataclass(frozen=True)
class Parameters:
    length_m: float = 0.08
    radius_m: float = 0.0035
    density_kg_m3: float = 2300.0
    gravity_m_s2: float = 9.8
    angle_deg: float = 20.0
    preload_n: float = 4.0
    board_speed_m_s: float = 0.5
    kx_n_m: float = 300.0
    ky_n_m: float = 1.0e6
    kn_n_m: float = 2.0e4
    cn_ns_m: float = 1.0
    mu_s: float = 0.62
    mu_k: float = 0.51
    stribeck_speed_m_s: float = 0.01
    sigma0_n_m: float = 1.0e5
    support_arm_ratio: float = 0.5
    sigma1_ns_m: float = 0.0
    sigma2_ns_m: float = 0.0
    tanh_speed_m_s: float | None = None
    rigid: bool = False


def lever_direction(angle):
    return np.array([np.sin(angle), -np.cos(angle)])


class InfeasiblePreload(ValueError):
    """The specified static configuration cannot support the preload."""


class Benchmark:
    """Centre-of-mass dynamics with both springs attached at the same point."""

    def __init__(self, p):
        self.p = p
        self.mass = p.density_kg_m3 * np.pi * p.radius_m**2 * p.length_m
        self.inertia = self.mass * p.length_m**2 / 12.0
        self.arm = p.length_m / 2.0
        self.support_arm = p.support_arm_ratio * p.length_m
        self.weight = self.mass * p.gravity_m_s2
        self.stiffness = np.array([p.kx_n_m, p.ky_n_m])
        if p.rigid and (p.tanh_speed_m_s is not None or p.sigma1_ns_m or p.sigma2_ns_m):
            raise ValueError("The rigid comparison uses dry LuGre.")

        angle = np.radians(p.angle_deg)
        penetration = 0.0 if p.rigid else p.preload_n / p.kn_n_m
        centre_y = self.arm * np.cos(angle) - penetration
        # Balance torque about the support H with the board at rest.
        tip_fx = -(p.preload_n - self.weight * self.support_arm /
                   (self.arm + self.support_arm)) * np.tan(angle)
        self.initial_friction_n = -tip_fx
        if abs(tip_fx) > p.mu_s * p.preload_n:
            raise InfeasiblePreload(
                f"Static preload requires {abs(tip_fx):.6g} N of friction, "
                f"with capacity {p.mu_s * p.preload_n:.6g} N.")
        self.initial = np.array([0.0, centre_y, angle, 0.0, 0.0, 0.0,
                                 -tip_fx / p.sigma0_n_m])
        self.spring_preload = np.array([-tip_fx, self.weight - p.preload_n])
        self.support0 = self.initial[:2] - self.support_arm * lever_direction(angle)

    def geometry(self, state):
        centre = state[:2]
        angle = state[2]
        tangent = np.array([np.cos(angle), np.sin(angle)])
        lever = self.arm * lever_direction(angle)
        tip = centre + lever
        support = centre - self.support_arm * lever_direction(angle)
        tip_velocity = state[3:5] + self.arm * tangent * state[5]
        return tip, support, tip_velocity, tangent

    def trial_normal(self, tip, tip_velocity):
        return -self.p.kn_n_m * tip[1] - self.p.cn_ns_m * tip_velocity[1]

    def spring_force(self, support):
        offset = support - self.support0.reshape((2,) + (1,) * (support.ndim - 1))
        shape = (2,) + (1,) * (support.ndim - 1)
        return self.spring_preload.reshape(shape) - self.stiffness.reshape(shape) * offset

    def normal_trial(self, state, tip, support, tip_velocity, tangent, spring):
        if not self.p.rigid:
            return self.trial_normal(tip, tip_velocity)
        jy = self.arm * tangent[1]
        free_torque = (-self.arm * tangent[0] * self.p.sigma0_n_m * state[6]
                       - self.support_arm * np.sum(tangent * spring, axis=0))
        free_accel = ((spring[1] - self.weight) / self.mass + jy * free_torque / self.inertia
                      + self.arm * tangent[0] * state[5]**2)
        return -free_accel / (1.0 / self.mass + jy**2 / self.inertia)

    def contact(self, state, tip, tip_velocity, trial, board_speed):
        p = self.p
        normal = np.maximum(trial, 0.0) if p.rigid else np.where(
            tip[1] < 0.0, np.maximum(trial, 0.0), 0.0)
        v_rel = tip_velocity[0] + board_speed
        if p.tanh_speed_m_s is not None:
            resistance = normal * p.mu_k * np.tanh(v_rel / p.tanh_speed_m_s)
            return normal, resistance, np.zeros_like(v_rel), np.zeros_like(v_rel), v_rel
        mu = p.mu_k + (p.mu_s - p.mu_k) * np.exp(-(v_rel / p.stribeck_speed_m_s)**2)
        force_scale = normal * mu
        active = normal > 1.0e-9 * p.preload_n
        bristle_force = p.sigma0_n_m * state[6]
        z_dot = np.where(active, v_rel - np.abs(v_rel) * bristle_force /
                         np.maximum(force_scale, 1.0e-12), 0.0)
        resistance = np.where(active, bristle_force + p.sigma1_ns_m * z_dot +
                              p.sigma2_ns_m * v_rel, 0.0)
        return normal, resistance, z_dot, force_scale, v_rel

    def sample_contact(self, state, board_speed):
        tip, support, velocity, tangent = self.geometry(state)
        trial = self.normal_trial(state, tip, support, velocity, tangent, self.spring_force(support))
        return self.contact(state, tip, velocity, trial, board_speed)

    def loads(self, state, board_speed):
        tip, support, tip_velocity, tangent = self.geometry(state)
        spring = self.spring_force(support)
        trial = self.normal_trial(state, tip, support, tip_velocity, tangent, spring)
        normal, resistance, z_dot, _, _ = self.contact(state, tip, tip_velocity, trial, board_speed)
        tip_force = np.array([-resistance, normal])
        net_force = spring + tip_force - np.array([0.0, self.weight])
        torque = np.dot(tangent, self.arm * tip_force - self.support_arm * spring)
        return net_force, torque, z_dot

    def rhs(self, _, state):
        # The scalar path avoids small array allocations at every solver stage.
        x, y, angle, vx, vy, omega, z = state
        p = self.p
        s, c = sin(angle), cos(angle)
        sx = self.spring_preload[0] - p.kx_n_m * (x - self.support_arm * s - self.support0[0])
        sy = self.spring_preload[1] - p.ky_n_m * (y + self.support_arm * c - self.support0[1])
        spring_torque = -self.support_arm * (c * sx + s * sy)
        bristle = p.sigma0_n_m * z
        if p.rigid:
            jy = self.arm * s
            free_torque = spring_torque - self.arm * c * bristle
            free_accel = ((sy - self.weight) / self.mass + jy * free_torque / self.inertia
                          + self.arm * c * omega**2)
            normal = max(0.0, -free_accel / (1.0 / self.mass + jy**2 / self.inertia))
        else:
            tip_y = y - self.arm * c
            normal = (max(0.0, -p.kn_n_m * tip_y - p.cn_ns_m * (vy + self.arm * s * omega))
                      if tip_y < 0.0 else 0.0)
        velocity = vx + self.arm * c * omega + p.board_speed_m_s
        z_dot = 0.0
        resistance = 0.0
        if p.tanh_speed_m_s is not None:
            resistance = normal * p.mu_k * tanh(velocity / p.tanh_speed_m_s)
        elif normal > 1.0e-9 * p.preload_n:
            mu = p.mu_k + (p.mu_s - p.mu_k) * exp(-(velocity / p.stribeck_speed_m_s)**2)
            z_dot = velocity - abs(velocity) * bristle / max(normal * mu, 1.0e-12)
            resistance = bristle + p.sigma1_ns_m * z_dot + p.sigma2_ns_m * velocity
        torque = spring_torque + self.arm * (-c * resistance + s * normal)
        return (vx, vy, omega, (sx - resistance) / self.mass,
                (sy + normal - self.weight) / self.mass, torque / self.inertia, z_dot)


def integrate(model, max_step, *, rtol=1e-10, atol=1e-12, horizon=0.012,
              stop_at_velocity_zero=True, method="DOP853", initial_omega=0.0):
    def first_stop(_, state):
        _, _, velocity, _ = model.geometry(state)
        return velocity[0] + model.p.board_speed_m_s

    def force_loss(_, state):
        tip, support, velocity, tangent = model.geometry(state)
        return model.normal_trial(state, tip, support, velocity, tangent, model.spring_force(support))

    def gap_opening(_, state):
        tip, _, _, _ = model.geometry(state)
        return -tip[1]

    events = (first_stop, force_loss) if model.p.rigid else (first_stop, force_loss, gap_opening)
    names = ("first_velocity_zero", "normal_force_zero", "geometric_gap_opening")[:len(events)]
    for event in events:
        event.direction = -1
        event.terminal = True
    first_stop.terminal = stop_at_velocity_zero

    initial = model.initial.copy()
    initial[5] = initial_omega
    sol = solve_ivp(model.rhs, (0.0, horizon), initial, method=method,
                    events=events, dense_output=True, max_step=max_step,
                    rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    event_times = {name: float(times[0]) if times.size else None
                   for name, times in zip(names, sol.t_events)}
    end = next((name for name in names if event_times[name] is not None), "time_limit")
    return sol, end, event_times


def analyze(model, sol, end, event_times):
    p = model.p
    t = np.linspace(0.0, sol.t[-1], 4001)
    state = sol.sol(t)
    tip, _, tip_velocity, _ = model.geometry(state)
    normal, resistance, _, g, velocity = model.sample_contact(state, p.board_speed_m_s)
    xi = tip[0] - tip[0, 0] + p.board_speed_m_s * t
    if np.any(np.diff(xi) <= 0):
        raise ValueError("The first-slip profile requires positive displacement in every interval.")

    # M1 uses load times sliding distance. M2 uses dry LuGre dissipation
    # after subtracting bristle energy storage from contact power.
    m1_rate = normal * np.abs(velocity)
    m2_rate = resistance**2 * np.abs(velocity) / np.maximum(g, 1e-12)
    rates = np.stack([m1_rate, m2_rate])
    cumulative = cumulative_trapezoid(rates, t, axis=1, initial=0)
    edges = np.linspace(0.0, xi[-1], 201)
    widths = np.diff(edges)
    densities = np.array([np.diff(np.interp(edges, xi, item)) / widths for item in cumulative])
    midpoint = 0.5 * (edges[:-1] + edges[1:])
    space_totals = densities @ widths
    integrals = np.array([total_deposit(t, rate) for rate in rates])
    work = total_deposit(t, resistance * velocity)
    bristle_energy_change = 0.5 * p.sigma0_n_m * (state[6, -1]**2 - state[6, 0]**2)

    trace = {
        "time_s": t, "xi_m": xi, "x_c_m": state[0], "y_c_m": state[1],
        "theta_deg": np.degrees(state[2]), "x_tip_m": tip[0],
        "y_tip_m": tip[1], "relative_speed_m_s": velocity, "normal_n": normal,
        "tangential_force_n": -resistance, "lugre_z_m": state[6],
        "m1_rate_n_m_s": m1_rate, "m2_rate_w": m2_rate,
    }
    profile = {
        "xi_left_m": edges[:-1], "xi_right_m": edges[1:], "xi_m": midpoint,
        "m1_density_n": densities[0], "m2_density_j_m": densities[1],
    }
    summary = {
        "end_reason": end, "events_s": event_times, "end_time_s": t[-1], "span_m": xi[-1],
        "normal_min_n": normal.min(), "normal_max_n": normal.max(),
        "normal_at_end_n": normal[-1], "tip_y_at_end_m": tip[1, -1],
        "tip_y_max_m": tip[1].max(), "relative_speed_at_end_m_s": velocity[-1],
        "coordinate_error_m": np.max(np.abs(xi - trace_coordinate(t, velocity))),
        "m1_integral_n_m": integrals[0], "m2_dissipation_j": integrals[1],
        "friction_work_j": work, "bristle_energy_change_j": bristle_energy_change,
        "energy_balance_relative_error": abs(work - bristle_energy_change - integrals[1]) / integrals[1],
        "spatial_mapping_relative_error": np.max(np.abs(space_totals / integrals - 1.0)),
    }
    return summary, trace, profile


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def run():
    p = Parameters()
    model = Benchmark(p)
    initial_force, initial_torque, z_dot = model.loads(model.initial, 0.0)
    require(np.max(np.abs(initial_force)) < 1e-10 and abs(initial_torque) < 1e-11 and z_dot == 0.0,
            "Static initial state failed force, torque, or LuGre balance.")

    convergence = []
    for max_step in (20e-6, 10e-6, 5e-6):
        sol, end, event_times = integrate(model, max_step)
        summary, trace, profile = analyze(model, sol, end, event_times)
        require(end == "first_velocity_zero" and summary["normal_min_n"] > 0
                and summary["tip_y_max_m"] < 0 and abs(summary["relative_speed_at_end_m_s"]) < 1e-8,
                f"Unexpected first event or contact state at max_step={max_step:g}.")
        convergence.append({"method": "DOP853", "max_step_s": max_step, "rtol": 1e-10,
                            "atol": 1e-12, "end_time_s": summary["end_time_s"],
                            "span_m": summary["span_m"], "normal_at_end_n": summary["normal_at_end_n"],
                            "m2_dissipation_j": summary["m2_dissipation_j"]})

    for key in ("end_time_s", "span_m", "normal_at_end_n", "m2_dissipation_j"):
        value = convergence[-1][key]
        worst_rel = max(abs(item[key] / value - 1.0) for item in convergence)
        require(worst_rel < 1.0e-4, f"Step-size convergence failed for {key}: {worst_rel:g}.")
    require(summary["coordinate_error_m"] < 1e-9, "Position and integrated speed disagree.")
    require(summary["energy_balance_relative_error"] < 1e-5, "LuGre energy balance failed.")
    require(summary["spatial_mapping_relative_error"] < 1e-10, "Spatial mapping lost the time integral.")

    report = {
        "parameters": asdict(p), "mass_kg": model.mass, "inertia_kg_m2": model.inertia,
        "initial_state": model.initial.tolist(),
        "spring_anchor_m": (model.support0 + model.spring_preload / model.stiffness).tolist(),
        "initial_spring_force_n": model.spring_preload.tolist(),
        "initial_net_force_n": initial_force.tolist(), "initial_net_torque_n_m": float(initial_torque),
        "initial_z_dot_m_s": float(z_dot), "result": summary, "convergence": convergence,
    }
    return report, trace, profile


def write_table(path, columns):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        writer.writerows(zip(*columns.values()))


def plot(path, trace, profile, summary):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    time_ms = trace["time_s"] * 1e3
    xi_mm = profile["xi_m"] * 1e3
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 6.2), layout="constrained")
    fig.suptitle("v3 first-slip benchmark", fontsize=16, fontweight="bold")
    axes[0, 0].plot(time_ms, trace["normal_n"], color="#176d68", label="Normal")
    axes[0, 0].plot(time_ms, -trace["tangential_force_n"], color="#ad5329", label="Friction magnitude")
    axes[0, 0].set(xlabel="Time (ms)", ylabel="Force (N)")
    axes[0, 0].legend(frameon=False)

    tag = f"First stop at {summary['end_time_s'] * 1e3:.3f} ms"
    axes[0, 1].axhline(0, lw=0.7, color="#8a939a")
    axes[0, 1].plot(time_ms, trace["relative_speed_m_s"], color="#315f8a")
    axes[0, 1].scatter([time_ms[-1]], [0], color="#315f8a", s=18, label=tag, zorder=3)
    axes[0, 1].set(xlabel="Time (ms)", ylabel="Tip speed relative to board (m/s)")
    axes[0, 1].legend(frameon=False)

    axes[1, 0].axhline(0, lw=0.8, ls="--", color="#8a939a", label="Board plane")
    axes[1, 0].plot(time_ms, trace["y_tip_m"] * 1e3, color="#315f8a")
    axes[1, 0].set(xlabel="Time (ms)", ylabel="Tip height in penalty model (mm)")
    axes[1, 0].legend(frameon=False)

    for name, label, color in (("m1_density_n", "M1 load–slip", "#176d68"),
                               ("m2_density_j_m", "M2 dissipated energy", "#ad5329")):
        density = profile[name]
        middle = density[40:160].mean()
        axes[1, 1].plot(xi_mm, density / middle, label=label, color=color)
    axes[1, 1].set(xlabel="Distance along board (mm)", ylabel="Proxy / middle mean")
    axes[1, 1].legend(frameon=False)
    for ax in axes.flat:
        ax.grid(alpha=0.22)
        ax.spines[["right", "top"]].set_visible(False)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=Path(__file__).resolve().parents[1] /
                        "simulations" / "V3_单痕基准_2026-11-17")
    args = parser.parse_args()
    report, trace, profile = run()
    args.outdir.mkdir(parents=True, exist_ok=True)
    write_table(args.outdir / "trajectory.csv", trace)
    write_table(args.outdir / "profile.csv", profile)
    (args.outdir / "summary.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n",
                                             encoding="utf-8")
    plot(args.outdir / "benchmark.png", trace, profile, report["result"])
    result = report["result"]
    print(f"{result['end_reason']}: {result['end_time_s'] * 1e3:.6f} ms, "
          f"span={result['span_m'] * 1e3:.6f} mm, N_end={result['normal_at_end_n']:.6f} N")


if __name__ == "__main__":
    main()
